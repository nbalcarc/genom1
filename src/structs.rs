use std::thread;

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
    pub fn push(&mut self, genome: Genome) -> Result<(), PhyloError>{
        //let refs: Vec<&mut TreeNode> = Vec::new(); //store all references to each layer here

        // if we have an empty tree, just push it
        if let TreeVertex::Floor(s) = &mut self.root.vertex {
            s.push(genome);
            return Ok(());
        }

        /* First we want to find 8 genomes to compare to, if available */

        let mut heads: Vec<(&TreeNode, u32, Vec<u8>)> = Vec::new(); //keep track of all heads (ref, heads, path)
        let mut genomes: Vec<&Genome>;
        heads.push((&mut self.root, 8, vec![0])); //push the root as the first head

        // Find all the genomes to run the kmer check on
        while heads.len() > 0 {
            println!("running with heads: {}", heads[0].1);

            // Repeat once per tuple in the current heads
            for i in 0..heads.len() {
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

                            //node_refs[*branch_index as usize] = node_ref;

                            // push the new head to the list of heads
                            //new_heads.push((&nodes[*branch_index as usize], branches[branch_index], this_path))

                            //this_path.push(nodes[*branch_index as usize].id.clone());
                            //new_heads.push((&mut nodes[*branch_index as usize], branches[branch_index], this_path));
                        }


                        //for branch_index in branches {
                        //    new_heads.push(&mut nodes[branch_index], )
                        //}

                        heads = new_heads;

                        

                    },
                    TreeVertex::Floor(v) => { //if we have a floor of genomes
                        // assign each head its own genome
                        println!("made it here, heads = {}", tup.1);
                        heads.remove(i);
                        
                        break;
                    }
                }



            }
        }
        Ok(())

    }

}







/*


    /// Push a new genome onto the tree
    pub fn push(&mut self, genome: Genome) -> Result<(), PhyloError>{
        //let refs: Vec<&mut TreeNode> = Vec::new(); //store all references to each layer here

        /* First we want to find 8 genomes to compare to, if available */

        let mut heads: Vec<(&mut TreeNode, u32, Vec<u8>)> = Vec::new(); //keep track of all heads (ref, heads, path)
        let mut genomes: Vec<&Genome>;
        heads.push((&mut self.root, 8, vec![0])); //push the root as the first head

        // Find all the genomes to run the kmer check on
        loop {

            // Repeat once per tuple in the current heads
            for tup in &mut heads {
                let mut new_heads: Vec<(&mut TreeNode, u32, Vec<u8>)> = Vec::new(); //create new vector to replace current one

                // if the TreeNode has fewer genomes than we have heads
                if &tup.0.count < &tup.1 {
                    tup.1 = tup.0.count; //reduce the number of heads
                }

                match &mut tup.0.vertex {
                    TreeVertex::Split(nodes) => { //if we have more splits
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
                       
                        // iterate through all the branches that will receive heads
                        for branch_index in branches.keys() {

                            let mut this_path = tup.2.clone(); //the whole path up to this current node
                            //let x = nodes[*branch_index as usize].id;

                            // update the local path
                            let this_id = nodes[*branch_index as usize].id.clone(); //retrieve the id of the next node
                            this_path.push(this_id); //push the next id to the path

                            let mut node_ref = nodes[*branch_index as usize]; //retrieve a reference to the next node
                            new_heads.push((&mut node_ref, branches[branch_index], this_path));
                            nodes[*branch_index as usize] = node_ref;

                            // push the new head to the list of heads
                            //new_heads.push((&nodes[*branch_index as usize], branches[branch_index], this_path))

                            //this_path.push(nodes[*branch_index as usize].id.clone());
                            //new_heads.push((&mut nodes[*branch_index as usize], branches[branch_index], this_path));
                        }


                        //for branch_index in branches {
                        //    new_heads.push(&mut nodes[branch_index], )
                        //}

                        

                        todo!();
                    },
                    TreeVertex::Floor(v) => { //if we have a floor of genomes
                        // assign each head its own genome
                        
                        break;
                    }
                }



            }


*/


