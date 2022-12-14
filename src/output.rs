use std::{path::Path, fs::{self, File}, io::Write};

use crate::{structs::{TreeNode, TreeVertex}, errors::PhyloError};


/// Produce an output file from a TreeNode
pub fn output_tree(root: &TreeNode) -> Result<(), PhyloError> {
    let path = Path::new("phylo_tree.txt");
    if path.exists() {
        fs::remove_file("phylo_tree.txt").map_err(|_| PhyloError::FileDeleteError)?; //if the output file exists already, override it
    }
    let mut file = File::create(path).map_err(|_| PhyloError::FileOpenError(String::from("Error opening the output file")))?;
    output_tree_recursive(root, &mut file, 0)
}

/// Internal recursive function that handles tree construction without worrying about initial conditions
fn output_tree_recursive(root: &TreeNode, file: &mut File, tabs: usize) -> Result<(), PhyloError> {
    match &root.vertex { //first find the type of node we're dealing with
        TreeVertex::Split(s) => {
            file.write_all((vec![' '; tabs*4].iter().collect::<String>() + "split:\n").as_bytes()) //write the node type
                .map_err(|_| PhyloError::FileWriteError)?; //write the node type

            // for each node in this split, run the function again
            for node in s {
                output_tree_recursive(node, file, tabs + 1)?;
            }
        },
        TreeVertex::Floor(f) => {
            file.write_all((vec![' '; tabs*4].iter().collect::<String>() + "floor:\n").as_bytes()) //write the node type
                .map_err(|_| PhyloError::FileWriteError)?; //write the node type
            
            // for each genome in this floor, print to file
            for genome in f {
                let mut slash_loc = genome.dir.rfind('/').ok_or(PhyloError::PathError(String::from(&genome.dir)))?; //filter out the file name
                let mut dir: String = genome.dir[..slash_loc].into();
                slash_loc = dir.rfind('/').ok_or(PhyloError::PathError(String::from(&genome.dir)))?; //filter out the upper directories
                dir = dir[slash_loc+1..].into();

                // write to file
                file.write_all((vec![' '; tabs*4+4].iter().collect::<String>() + &dir + "\n").as_bytes()).map_err(|_| PhyloError::FileWriteError)?;
            }
        }
    }
    Ok(())
}
