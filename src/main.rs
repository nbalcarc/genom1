use crate::structs::{PhyloTree, Genome};
use std::{path::Path, env, fs::{self, File}, os::unix::prelude::FileExt};

mod algorithms;
mod errors;
mod structs;


/// A function dedicated to testing functionality
fn testing() {
    let mut ve = vec![(0, 1), (2, 3), (4, 5)];
    for t in &mut ve {
        t.0 = 5;
        t.1 = 6;
    }
    //dbg!(ve);

    let items = vec![1,  2, 3,  4,   5,   6, 7,  8];
    let probs = vec![40, 3, 18, 100, 200, 1, 35, 16];
    //dbg!(algorithms::random_weighted(items, probs, 10, true));

    let mut x = Box::new(3);
    let y = &mut x;
    let z = x;

    let mut buff = vec![1; 32];
    let file = File::open("/home/terrior/Programming/genome-tree/src/test1.txt").unwrap();
    file.read_exact_at(&mut buff, 0);
    //dbg!(buff);

    let kmers = algorithms::generate_kmers("/home/terrior/Programming/genome-tree/genomes/Abaca_bunchy_top_virus/GCF_000872625.1_ViralMultiSegProj28697_genomic.fna", 12, 64).unwrap();
    let kmers1 = algorithms::generate_kmers("/home/terrior/Programming/genome-tree/genomes/Abaca_bunchy_top_virus/GCF_000872625.1_ViralMultiSegProj28697_genomic.fna", 12, 64).unwrap();
    let kmers2 = algorithms::generate_kmers("/home/terrior/Programming/genome-tree/genomes/Akhmeta_virus/GCF_006452035.1_ASM645203v1_genomic.fna", 12, 64).unwrap();

    let mut genome = Genome {
        path: Vec::new(),
        dir: String::from("/home/terrior/Programming/genome-tree/genomes/Abaca_bunchy_top_virus/GCF_000872625.1_ViralMultiSegProj28697_genomic.fna"),
        kmers: kmers,
        closest_relative: Vec::new(),
        closest_distance: 0,
    };
    
    let mut genome1 = Genome {
        path: Vec::new(),
        dir: String::from("/home/terrior/Programming/genome-tree/genomes/Abaca_bunchy_top_virus/GCF_000872625.1_ViralMultiSegProj28697_genomic.fna"),
        kmers: kmers1,
        closest_relative: Vec::new(),
        closest_distance: 0,
    };

    let mut genome2 = Genome {
        path: Vec::new(),
        dir: String::from("/home/terrior/Programming/genome-tree/genomes/Akhmeta_virus/GCF_006452035.1_ASM645203v1_genomic.fna"),
        kmers: kmers2,
        closest_relative: Vec::new(),
        closest_distance: 0,
    };

    dbg!(algorithms::kmer_similarity(&genome, &genome));
    dbg!(algorithms::kmer_similarity(&genome, &genome1));
    dbg!(algorithms::kmer_similarity(&genome, &genome2));


    //for i in 0..30000 {
    //    println!("{} of {}", i, 30000);
    //    for j in 0..30000 {
    //
    //    }
    //}

}


/// Handle all the tree generation
fn tree_generation() {
    let mut tree = structs::PhyloTree::new();
    
    let dir = env::current_dir().unwrap().to_str().unwrap().to_owned(); //get the current working directory
    let paths = fs::read_dir(dir + "/genomes").unwrap(); //get all paths in the genomes directory

    // iterate through all genome folders
    for path in paths {
        //println!("Name: {}", path.unwrap().path().display());
        // iterate through all files in the genome folder
        for file in fs::read_dir(path.unwrap().path().to_str().unwrap()).unwrap() {
            let file_path = String::from(file.unwrap().path().to_str().unwrap());
            if !&file_path.ends_with(".fna") {
                continue;
            }
            //println!("{}", file_path);
            let kmers = algorithms::generate_kmers(&file_path, 12, 64).unwrap();
            let mut genome = Genome {
                path: Vec::new(),
                dir: file_path,
                kmers: kmers,
                closest_relative: Vec::new(),
                closest_distance: 0,
            };

            tree.push(genome);

        }
    }
}


/// Entry point
fn main() {
    //testing();
    tree_generation();
}
