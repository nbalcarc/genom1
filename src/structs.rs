use std::{thread, fs, sync::{Arc, Mutex}};

use crate::{errors::PhyloError, algorithms::{self, retrieve_genome}};


/// Establishes the structure of our phylogenetic tree
#[derive(Debug, Clone)]
pub struct TreeNode {
    pub id: u8,             // unique identifier used for finding genome paths
    pub count: u32,         // the total count of genomes under this node 
    pub vertex: TreeVertex, // decides the structure of this node
}
impl TreeNode {

    /// Initializes a new TreeNode with a TreeVertex::Split
    //pub fn new_with_split(id: u8, count: u32) -> Self {
    //    TreeNode { id: id, vertex: TreeVertex::Split(Vec::new()), count: count }
    //}

    /// Initializes a new TreeNode with a TreeVertex::Floor
    pub fn new_with_floor(id: u8, count: u32) -> Self {
        TreeNode { id: id, vertex: TreeVertex::Floor(Vec::new()), count: count }
    }

    /// If we have a floor, we'll switch to a split where one of the children is our current floor
    pub fn split(&mut self, id: u8) {
        if let TreeVertex::Floor(_) = &mut self.vertex {
            let mut floor_node = TreeNode::new_with_floor(id, self.count); //the new node that will point to our current floor
            let cur_floor = std::mem::replace(&mut self.vertex, TreeVertex::new_split()); //retrieve the current floor
            floor_node.vertex = cur_floor; //place the floor into the node
            self.vertex.push_node(floor_node); //push the newly created node containing our old floor into our split
        }
    }

    pub fn find<'s>(&'s self, number_heads: u32) -> Vec<&'s Genome> {
        /* First we want to find 8 genomes to compare to, if available */

        let mut heads: Vec<(&TreeNode, u32)> = Vec::new(); //keep track of all heads (ref, heads)
        let mut genomes: Vec<&Genome> = Vec::with_capacity(number_heads.try_into().unwrap()); //result
        heads.push((self, number_heads)); //push the root as the first head

        let mut thread_rng = rand::thread_rng();

        // Find all the genomes to run the kmer check on
        while heads.len() > 0 {

            // Repeat once per tuple in the current heads
            for i in 0..heads.len() {

                let mut tup = heads[i]; //get the current tuple of information
                let mut new_heads: Vec<(&TreeNode, u32)> = Vec::new(); //create new vector to replace current one

                // if the TreeNode has fewer genomes than we have heads
                if &tup.0.count < &tup.1 {
                    tup.1 = tup.0.count; //reduce the number of heads
                }

                // if more splits, then update the list of heads
                // if reach floor, then add to list of genomes and clip heads
                match &tup.0.vertex {
                    TreeVertex::Split(nodes) => { //if we have more splits

                        // for each node, allocate a certain number of heads to it
                        let count = nodes.len().clone(); //number of nodes on this floor
                        let mut weights: Vec<u32> = Vec::new(); //keep track of how many genomes each node has
                        let mut indices: Vec<u32> = Vec::new(); //indices

                        // prepare the weights, iterate for each node in the floor
                        // this sets the weight of each node as the number of genomes it holds
                        for i in 0u32..(count as u32) {
                            weights.push(nodes[i as usize].count);
                            indices.push(i);
                        }

                        // get all the branches our heads will go to
                        let weight_results = algorithms::random_weighted(indices, weights, tup.1, false);
                        let branches = algorithms::vec_to_dict(weight_results);
                       
                        // iterate through all the branches that will receive heads
                        for branch_index in branches.keys() {

                            // push the new head to the list of heads
                            let node_ref: &TreeNode = &nodes[*branch_index as usize]; //retrieve a reference to the next node
                            new_heads.push((node_ref, branches[branch_index]));
                        }

                        heads = new_heads;
                    },
                    TreeVertex::Floor(v) => { //if we have a floor of genomes
                        // assign each head its own genome
                        heads.remove(i);

                        use rand::seq::SliceRandom;
                        genomes.extend(v.choose_multiple(&mut thread_rng, tup.1 as _));
                        
                        break;
                    }
                }
            }
        }
        genomes        
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
    //pub fn new_floor() -> Self {
    //    TreeVertex::Floor(Vec::new())
    //}

    /// Push a new node to a Split
    pub fn push_node(&mut self, node: TreeNode) {
        if let TreeVertex::Split(v) = self {
            v.push(node);
        }
    }

    // Push a new genome to a Floor 
    //pub fn push_genome(&mut self, genome: Genome) {
    //    if let TreeVertex::Floor(v) = self {
    //        v.push(genome);
    //    }
    //}
}


/// Represents a single genome
#[derive(Debug, Clone)]
pub struct Genome {
    pub path: Vec<u8>,              // the path to reach this genome
    pub dir: String,                // the directory of the genome
    pub kmers: Vec<String>,         // the list of kmers for this genome
    pub closest_distance: usize,      // Levenshtein distance between this genome and its closest relative
}


/// Manages the phylogenetic tree
#[derive(Debug)]
pub struct PhyloTree {
    pub root: TreeNode,
    pub next_index: u8,     // used to decide the next TreeNode id
}
impl PhyloTree {

    /// Create a new phylogenetic tree
    pub fn new() -> Self {
        PhyloTree { root: TreeNode::new_with_floor(0, 0), next_index: 1 }
    }

    /// Push a new genome onto the tree
    pub fn push(&mut self, mut genome: Genome) -> Result<(), PhyloError>{
        println!("======== PUSHING A NEW GENOME NOW ========");

        let root_count_increment: bool;

        // if we have an empty tree, just push it
        if let TreeVertex::Floor(s) = &mut self.root.vertex {
            genome.path = vec![0, 0];
            genome.closest_distance = std::usize::MAX;
            if s.len() == 0 {
                println!("---> CASE 0: First genome, push to top floor");
                s.push(genome);
                self.root.count = 1;
                return Ok(());
            }
            root_count_increment = true;
        } else {
            root_count_increment = false;
        }

        // prepare variables that will be updated each iteration
        let mut checked: Vec<u8> = Vec::new(); //the paths we've checked so far
        let mut num_checked: u32; //the number of nodes we're checking this iteration
        let mut cur = &self.root; //the node we're checking next
        checked.push(cur.id);
        let mut genomes: Vec<&Genome>;

        // find the next set of 8 nodes in this loop
        'main_loop: loop {
            // retrieve the genomes for this node
            genomes = cur.find(8); //retrieve a random set of 8 genomes
            num_checked = cur.count; //update the number of genomes we've looked over

            // decide if we exit or do another iteration ======= THIS IS WHERE WE DECIDE WHETHER TO START THE INSERTION STEP =======
            if num_checked < 9 { //we have enough genomes to start the insertion step
                break 'main_loop;                    

            } else { //compare similarities before starting next iteration or exiting

                // for each genome, calculate the kmer similarity
                let mut distances = Vec::new();
                for cur_genome in &genomes {
                    distances.push((algorithms::kmer_similarity(*cur_genome, &genome), *cur_genome));
                }
                let best_genome = *distances.iter().max_by_key(|a|a.0).unwrap(); //(similarity, ref), the best genome
                let node_path = algorithms::get_full_path(&self.root, &best_genome.1.path)?; //get the full list of nodes leading to the genome's parent

                // check each node to see if we've checked it or not
                for i in 0..node_path.len() {
                    // skip all checked nodes (last one should always be unchecked naturally)
                    // if this node is larger than threshold and next one is smaller than threshold, then choose this one
                    // if this node is larger than threshold and the next one is also, then continue
                    // if this node is smaller than threshold or the last one, then choose this one

                    if checked.contains(&node_path[i].id) { //we've already checked this node
                        continue;
                    }
                    let threshold = num_checked / 8;
                    if node_path[i].count >= threshold && i+1 < node_path.len() && node_path[i+1].count >= threshold {
                        // if our node is large enough to hold the threshold, and if the next node is also large enough then continue
                        continue
                    }

                    // if we hit here, then the node could the the last node before going under the threshold or last node entirely
                    // alternatively, this node must be smaller than the threshold and previous ones were blocked, so run again on this node
                    cur = node_path[i];
                    checked.push(cur.id);
                }

            }

        }

        // if we've successfully reached the final 8 genomes
        // launch the insertion protocol
        // -do real comparisons on all genomes
        // -find the closest relative
        // -resort the tree if need be, and insert the genome
        let distances: Arc<Mutex<Vec<(usize, Vec<u8>)>>> = Arc::new(Mutex::new(Vec::new())); // (distance, &Genome)

        let genome_str = fs::read_to_string(&genome.dir).map_err(|_| PhyloError::FileOpenError(String::from(&genome.dir)))?;
        let mut threads: Vec<thread::JoinHandle<()>> = Vec::new();

        // for each genome, generate a thread that runs the levenshtein algorithm
        for cur_genome in genomes {

            // copy variables that we'll need in the closure
            let dist_arc = distances.clone();
            let genome_str0 = genome_str.clone();
            let genome_path = cur_genome.path.clone();
            let genome_dir = cur_genome.dir.clone();

            // launch a new thread for levenshtein distance
            let cur_thread = thread::spawn( move || {
                let genome_str1 = fs::read_to_string(&genome_dir).map_err(|_| PhyloError::FileOpenError(String::from(genome_dir))).unwrap();
                let distance = algorithms::levenshtein(&genome_str0, &genome_str1);
                dist_arc.lock().unwrap().push((distance,genome_path));
            });
            threads.push(cur_thread);
        }

        // rejoin all threads back together
        for thr in threads {
            thr.join().map_err(|_| PhyloError::GenomeInsertError(String::from("Error joining threads after running Levenshtein")))?;
        }

        let distances = distances.lock().unwrap().clone();
        if let Some((best_dist,best_genome_path)) = distances.iter().min_by_key(|a|a.0) { //(path, distance), the best genome
            let best_genome_mut = retrieve_genome(&mut self.root, &best_genome_path)?;

            // update our new genome
            genome.closest_distance = *best_dist;
            
            // first retrieve the parent node of the CR
            let mut relative_distance = genome.closest_distance as f64 / best_genome_mut.closest_distance as f64;

            // retrieve the parent node of the CR
            let parent_node = algorithms::get_mut_node_and_increment(&mut self.root, &best_genome_path)?;

            // consider an alternative case for Case 1, where the existing genome is the only genome in its floor
            if relative_distance <= 0.85 {
                if let TreeVertex::Floor(ref mut f) = parent_node.vertex {

                    // if the length of the existing floor has only one member
                    if f.len() == 1 {
                        println!("------> CASE 1.1: Closest relative is the only genome in its branch, don't split");
                        relative_distance = 1.0;
                    }
                }
            }

            // CASE 1 TODO consider the case where the closest relative is the only genome in that branch, simply place new genome there
            if relative_distance <= 0.85 { //create a new branch, bring the new genome and its closest relative into it, update genome paths
                println!("---> CASE 1: New branch, put genome and closest relative there");

                parent_node.split(self.next_index);
                
                if let TreeVertex::Split(ref mut s) = parent_node.vertex {

                    // update counts
                    s[0].count -= 1;
                    s[1].count = 1;

                    let mut genome_path = best_genome_path.clone(); //update the path of the newly inserted genome
                    let distance = genome.closest_distance; //retrieve the distance

                    genome_path.remove(genome_path.len()-1); //remove the index
                    genome_path.push(self.next_index+1); //push the new split
                    let mut cr_path = genome_path.clone(); //update the path of the existing node (since the tree will be reorganized)
                    genome_path.push(0); //the new genome will be added to index 0
                    cr_path.push(1); //the cr genome will be added to index 1
                    genome.path = genome_path.clone();
                    let mut closest_relative;

                    // open the first branch
                    if let TreeVertex::Floor(ref mut f_original) = s[0].vertex {
                        closest_relative = f_original.remove(best_genome_path[best_genome_path.len()-1] as usize); //grab the closest relative Genome so we can move it

                        // update all the paths because the vector was just shifted
                        for i in 0..f_original.len() {
                            let mut cur_genome = &mut f_original[i];
                            let mut new_path = cur_genome.path.clone();
                            new_path.remove(new_path.len()-1); //remove the index
                            new_path.push(self.next_index); //push the new split
                            new_path.push(i.try_into().unwrap()); //push the new index
                            cur_genome.path = new_path;
                        } 
                    } else { //this is here to make the compiler happy
                        closest_relative = Genome {
                            path: Vec::new(),
                            dir: String::from(""),
                            kmers: Vec::new(),
                            closest_distance: 0,
                        }
                    }

                    // create the second branch, where the new genome and its CR will reside
                    s.push(TreeNode::new_with_floor(self.next_index + 1, parent_node.count));
                    if let TreeVertex::Floor(ref mut f_new) = s[1].vertex {

                        // update the genome's path to its closest relative (who is now in the same floor as it)
                        closest_relative.path = cr_path;
                        closest_relative.closest_distance = distance; //always want to update the distance because the ratio is 0.85 or less
                        f_new.push(genome);
                        f_new.push(closest_relative);
                    }

                }
                self.next_index = self.next_index + 2;

            // CASE 2
            } else if relative_distance >= 1.17 { //create a new branch, place the new genome there
                println!("---> CASE 2: New branch, put genome there");
                parent_node.split(self.next_index);

                // retrieve the parent node above where the floor node is
                if let TreeVertex::Split(ref mut s) = parent_node.vertex {

                    let mut genome_path = best_genome_path.clone(); //update the path of the newly inserted genome
                    let mut cr_path = best_genome_path.clone(); //update the path of the existing node (since the tree will be reorganized)

                    genome_path.remove(genome_path.len()-1); //remove the index
                    genome_path.push(self.next_index+1); //push the new split
                    genome_path.push(0); //the new genome will be added to index 0

                    cr_path.insert(cr_path.len()-1, self.next_index); //push the new split
                    
                    genome.path = genome_path.clone();

                    // create the second branch, where the new genome and its CR will reside
                    s.push(TreeNode::new_with_floor(self.next_index + 1, 1));

                    // update the paths of the genomes in the original floor
                    if let TreeVertex::Floor(ref mut f) = s[0].vertex {
                        for cur_genome in f {
                            let mut new_path = cur_genome.path.clone(); //fix the location of this genome
                            let length = new_path.len();
                            new_path.insert(length-1, self.next_index);
                            cur_genome.path = new_path;
                            // NOTE: since the ratio is 1.17 or greater, it means the new genome doesn't beat out the old closest distance
                            //       as such, don't update the closest distance on our cur_genome
                        }
                    }

                    // push the new genome onto its branch
                    if let TreeVertex::Floor(ref mut f) = s[1].vertex {
                        f.push(genome);
                    }
                }
                self.next_index = self.next_index + 2;

            // CASE 3
            } else { //place the new genome in the same branch as its closest relative
                println!("---> CASE 3: Put genome in same branch as closest relative");
                if let TreeVertex::Floor(ref mut f) = parent_node.vertex {
                    let mut new_path = best_genome_path.clone(); //update the path of the newly inserted genome
                    new_path.remove(new_path.len()-1);
                    new_path.push(f.len().try_into().unwrap());
                    
                    let closest_relative = &mut f[best_genome_path[best_genome_path.len()-1] as usize];
                    if closest_relative.closest_distance > genome.closest_distance { //if we need to update cr's closest distance, do it here
                        closest_relative.closest_distance = genome.closest_distance;
                    }
                    genome.path = new_path;
                    f.push(genome);
                }
            }

            // the function get_mut_node_and_increment fights our root's count when we have a floor, so we account for that here
            if root_count_increment {
                self.root.count += 1;
            }
            return Ok(());
        }

        Err(PhyloError::GenomeInsertError(String::from("Distances vector was empty, could find no nodes to compare to")))
    }
}

