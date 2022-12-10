/// Establishes the structure of our phylogenetic tree
pub struct TreeNode {
    id: i32,            // unique identifier used for finding genome paths
    vertex: TreeVertex, // decides the structure of this node
}


/// Decides whether a node turns into a floor or a split
pub enum TreeVertex {
    Split(i32, Vec<TreeNode>),  // used to split into multiple branches
    Floor(i32, Vec<Genome>),    // used to contain a list of genomes
}


/// Represents a single genome
pub struct Genome {
    path: Vec<i32>,             // the path to reach this genome
    dir: String,                // the directory of the genome
    ngrams: Vec<String>,        // the list of ngrams for this genome
    closest_relative: Vec<i32>, // the path of the closest relative
    closest_distance: u32,      // Levenshtein distance between this genome and its closest relative
}


/// Manages the phylogenetic tree
pub struct PhyloTree {
    root: TreeVertex,
}



