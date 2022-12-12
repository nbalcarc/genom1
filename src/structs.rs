use std::{thread, fs::{File, self}, sync::{Arc, Mutex}};

use rand::Rng;

use crate::{errors::PhyloError, algorithms};

/// Establishes the structure of our phylogenetic tree
#[derive(Debug, Clone)]
pub struct TreeNode {
    id: u8,             // unique identifier used for finding genome paths
    vertex: TreeVertex, // decides the structure of this node
    count: u32,         // the total count of genomes under this node 
}
impl TreeNode {

    /// Initializes a new TreeNode with a TreeVertex::Split
    pub fn new_with_split(id: u8, count: u32) -> Self {
        TreeNode { id: id, vertex: TreeVertex::Split(Vec::new()), count: count }
    }

    /// Initializes a new TreeNode with a TreeVertex::Floor
    pub fn new_with_floor(id: u8, count: u32) -> Self {
        TreeNode { id: id, vertex: TreeVertex::Floor(Vec::new()), count: count }
    }

    /// If we have a floor, we'll switch to a split where one of the children is our current floor
    pub fn split(&mut self, id: u8) {
        let mut floor_node = TreeNode::new_with_floor(id, self.count); //the new node that will point to our current floor
        let cur_floor = std::mem::replace(&mut self.vertex, TreeVertex::new_split()); //retrieve the current floor
        floor_node.vertex = cur_floor; //place the floor into the node
        self.vertex.push_node(floor_node); //push the newly created node containing our old floor into our split
    }
}


/// Decides whether a node turns into a floor or a split
#[derive(Debug, Clone)]
pub enum TreeVertex {
    Split(Vec<TreeNode>),   // used to split into multiple branches
    Floor(Vec<Genome>),     // used to contain a list of genomes
}
impl TreeVertex {

    /// Initializes a new TreeVertex split
    pub fn new_split() -> Self {
        TreeVertex::Split(Vec::new())
    }

    /// Initializes a new TreeVertex floor
    pub fn new_floor() -> Self {
        TreeVertex::Floor(Vec::new())
    }

    /// Push a new node to a Split
    pub fn push_node(&mut self, node: TreeNode) {
        if let TreeVertex::Split(v) = self {
            v.push(node);
        }
    }

    /// Push a new genome to a Floor 
    pub fn push_genome(&mut self, genome: Genome) {
        if let TreeVertex::Floor(v) = self {
            v.push(genome);
        }
    }
}


/// Represents a single genome
#[derive(Debug, Clone)]
pub struct Genome {
    pub path: Vec<u8>,              // the path to reach this genome
    pub dir: String,                // the directory of the genome
    pub kmers: Vec<String>,         // the list of kmers for this genome
    pub closest_relative: Vec<i32>, // the path of the closest relative
    pub closest_distance: u32,      // Levenshtein distance between this genome and its closest relative
}


/// Manages the phylogenetic tree
pub struct PhyloTree {
    root: TreeNode,
    next_index: u8,     // used to decide the next TreeNode id
}
impl PhyloTree {

    /// Create a new phylogenetic tree
    pub fn new() -> Self {
        PhyloTree { root: TreeNode::new_with_floor(0, 0), next_index: 1 }
    }

    /// Push a new genome onto the tree
    pub fn push(&mut self, mut genome: Genome) -> Result<(), PhyloError>{
        //let refs: Vec<&mut TreeNode> = Vec::new(); //store all references to each layer here

        // if we have an empty tree, just push it
        if let TreeVertex::Floor(s) = &mut self.root.vertex {
            genome.path = vec![0];
            if s.len() == 0 {
                println!("exiting early, root");
                s.push(genome);
                self.root.count = 1;
                return Ok(());
            }
        }

        /* First we want to find 8 genomes to compare to, if available */

        let mut heads: Vec<(&TreeNode, u32, Vec<u8>)> = Vec::new(); //keep track of all heads (ref, heads, path)
        let mut genomes: Vec<(usize, &Genome)> = Vec::with_capacity(8);
        heads.push((&mut self.root, 8, vec![0])); //push the root as the first head

        // Find all the genomes to run the kmer check on
        while heads.len() > 0 {
            println!("running with heads: {}", heads[0].1);

            // Repeat once per tuple in the current heads
            for i in 0..heads.len() {
                println!("running on head: {} out of {}", i, heads.len());

                let mut tup = heads[i].clone(); //get the current tuple of information
                let mut new_heads: Vec<(&TreeNode, u32, Vec<u8>)> = Vec::new(); //create new vector to replace current one

                // if the TreeNode has fewer genomes than we have heads
                if &tup.0.count < &tup.1 {
                    tup.1 = tup.0.count; //reduce the number of heads
                }

                // if more splits, then update the list of heads
                // if reach floor, then add to list of genomes and clip heads
                match &tup.0.vertex {
                    TreeVertex::Split(nodes) => { //if we have more splits

                        println!("made it into split, heads: {}", tup.1);

                        // for each node, allocate a certain number of heads to it
                        let count = nodes.len().clone(); //number of nodes on this floor
                        let mut weights: Vec<u32> = Vec::new(); //keep track of how many genomes each node has
                        let mut indices: Vec<u32> = Vec::new(); //indices

                        // prepare the weights, iterate for each node in the floor
                        for i in 0..count {
                            weights.push(nodes[i].count);
                            indices.push(i.try_into().map_err(|_| PhyloError::ConversionError)?);
                        }

                        // get all the branches our heads will go to
                        let branches = algorithms::vec_to_dict(algorithms::random_weighted(weights, indices, tup.1, false));
                        //let mut node_refs = nodes.clone();
                       
                        // iterate through all the branches that will receive heads
                        for branch_index in branches.keys() {

                            let mut this_path = tup.2.clone(); //the whole path up to this current node
                            //let x = nodes[*branch_index as usize].id;

                            // update the local path
                            let this_id = nodes[*branch_index as usize].id.clone(); //retrieve the id of the next node
                            this_path.push(this_id); //push the next id to the path

                            // push the new head to the list of heads
                            let node_ref: &TreeNode = &nodes[*branch_index as usize]; //retrieve a reference to the next node
                            new_heads.push((node_ref, branches[branch_index], this_path));
                        }

                        heads = new_heads;
                    },
                    TreeVertex::Floor(v) => { //if we have a floor of genomes
                        // assign each head its own genome
                        println!("made it to a floor, heads = {}", tup.1);
                        heads.remove(i);

                        // establish vector with an element per genome
                        let mut options: Vec<usize> = Vec::new();
                        for i in 0..v.len() {
                            options.push(i);
                        }

                        // choose randomly which genomes to look at, one iteration per head
                        let mut rand = rand::thread_rng();
                        for i in 0..tup.1 {
                            let chosen = rand.gen_range(0..options.len());
                            println!("option chosen: {}", options[chosen]);
                            genomes.push((options[chosen], &v[chosen]));
                            options.remove(chosen);
                        }

                        dbg!("here");
                        
                        break;
                    }
                }
            }
        }

        //dbg!(&genomes);
        //println!("{}", genomes.len());
        dbg!("here, down");

        // if we've successfully reached the final 8 genomes
        if genomes.len() < 9 {
            // launch the insertion protocol
            // -do real comparisons on all genomes
            // -find the closest relative
            // -resort the tree if need be, and insert the genome
            let distances: Arc<Mutex<Vec<(Vec<u8>, usize)>>> = Arc::new(Mutex::new(Vec::new()));

            let genome_str = fs::read_to_string(&genome.dir).map_err(|_| PhyloError::FileOpenError(String::from(&genome.dir)))?;
            let mut threads: Vec<thread::JoinHandle<()>> = Vec::new();

            // for each genome
            for cur_tup in genomes {
                let cur_genome = cur_tup.1;

                // copy variables that we'll need in the closure
                let dist_arc = distances.clone();
                let genome_str0 = genome_str.clone();
                let cur_genome0 = cur_genome.clone();
                
                let mut new_path = cur_genome0.path;
                println!("old path: {:?}", new_path);
                new_path.push(cur_tup.0.try_into().unwrap());
                println!("new path: {:?}", new_path);

                // launch a new thread for levenshtein distance
                let cur_thread = thread::spawn( move || {
                    let genome_str1 = fs::read_to_string(&cur_genome0.dir).map_err(|_| PhyloError::FileOpenError(String::from(&cur_genome0.dir))).unwrap();
                    dbg!("starting levenshtein");
                    dist_arc.lock().unwrap().push((new_path, algorithms::levenshtein(&genome_str0, &genome_str1)));
                    
                });
                dbg!("idk here?");
                threads.push(cur_thread);
                //let genome_str1 = fs::read_to_string(&cur_genome.dir).map_err(|_| PhyloError::FileOpenError(String::from(&cur_genome.dir)))?;
                //distances.push(algorithms::levenshtein(&genome_str, &genome_str1));
            }

            dbg!("now we waiting");

            // rejoin all threads back together
            let x = threads.len();
            for thr in threads {
                println!("size of threads: {}", x);
                thr.join();
            }
            //dbg!(distances);

            let distances = distances.lock().unwrap().clone();
            let dummy = (vec![0 as u8], std::usize::MAX);

            let best = distances.iter().fold(&dummy, |accum, elem| if elem.1 < accum.1 {elem} else {accum});

        } else {
            // launch the filter protocol
            // -find the closest relative
            // -start from the bottom, looking for the first node to have at most 1/8 of all nodes
            // -if that node is the node we just searched, then just choose the next node down
            // -"recurse"
        }

        dbg!("here, down down");

        Ok(())

    }

}

