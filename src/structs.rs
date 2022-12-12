/// Establishes the structure of our phylogenetic tree
pub struct TreeNode {
    id: u8,             // unique identifier used for finding genome paths
    vertex: TreeVertex, // decides the structure of this node
    count: u32,         // the total count of genomes under this node 
}


/// Decides whether a node turns into a floor or a split
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
}


/// Represents a single genome
#[derive(Debug)]
pub struct Genome {
    pub path: Vec<u8>,              // the path to reach this genome
    pub dir: String,                // the directory of the genome
    pub kmers: Vec<String>,         // the list of kmers for this genome
    pub closest_relative: Vec<i32>, // the path of the closest relative
    pub closest_distance: u32,      // Levenshtein distance between this genome and its closest relative
}


/// Manages the phylogenetic tree
pub struct PhyloTree {
    root: TreeVertex,
    index: u8,          // used to decide the next TreeNode id
}
impl PhyloTree {

    /// Create a new phylogenetic tree
    pub fn new() -> Self {
        PhyloTree { root: TreeVertex::new_floor(), index: 0 }
    }

    /// Push a new genome onto the tree
    pub fn push(&mut self, genome: Genome) {
        let refs: Vec<&mut TreeNode> = Vec::new();

    }

}


