use std::thread;

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
        PhyloTree { root: TreeNode::new_with_floor(0, 0), next_index: 0 }
    }

    /// Push a new genome onto the tree
    pub fn push(&mut self, genome: Genome) {
        //let refs: Vec<&mut TreeNode> = Vec::new(); //store all references to each layer here
        let mut paths: Vec<(u8, u32)> = Vec::new();
        paths.push((self.root.id, 8)); //push the root as the first head

    }

}


