use crate::structs::PhyloTree;

use std::{path::Path, env, fs};

mod algorithms;
mod errors;
mod structs;

fn main() {
    let result = algorithms::levenshtein("kit", "glimmen");
    println!("{}", result);
    println!("{}", algorithms::levenshtein("kitten", "alderkitten"));
    println!("{}", algorithms::levenshtein("alderkitten", "kitten"));
    dbg!(algorithms::file_size("/home/terrior/Programming/genome-tree/src/test.txt"));

    //algorithms::generate_kmers(file_dir, k, num);
    let x = env::current_dir().unwrap().to_str().unwrap().to_owned();
    let dir = env::current_dir().unwrap().as_path().display().to_string();
    println!("{}", x);
    println!("{}", dir);
    let paths = fs::read_dir(x).unwrap();
    for path in paths {
        println!("Name: {}", path.unwrap().path().display())
    }


    let tree = PhyloTree::new();
}
