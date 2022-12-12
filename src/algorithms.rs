use std::{fs::{self, File}, io::Read, os::unix::prelude::FileExt};
use rand::Rng;

use crate::errors::PhyloError;

/// Calculate the Levenshtein distance between two strings
pub fn levenshtein(first: &str, second: &str) -> usize {
    let long: &[u8];
    let short: &[u8];

    // first figure out the longest and shortest strings
    if first.len() > second.len() {
        long = first.as_bytes();
        short = second.as_bytes();
    } else {
        short = first.as_bytes();
        long = second.as_bytes();
    }

    // get the sizes of both strings
    let longd = long.len();
    let shortd = short.len();

    // initialize the two vectors we need, sized according to short to reduce mem usage
    let mut prev = vec![0; short.len() + 1];
    let mut cur = vec![0; short.len() + 1];

    // initialize all elements in the starter vector
    for i in 0..shortd {
        prev[i] = i;
    }

    // iterate once for every letter in the long word
    for y in 0..longd {
        cur[0] = y+1; //set the character index

        // iterate once for every letter in the short word (size of the arrays)
        for x in 0..shortd {
            let del_cost = prev[x+1] + 1;   // generate cost of deletion
            let ins_cost = cur[x] + 1;      // generate cost of insertion
            let sub_cost: usize;

            // generate cost of substitution
            if long[y] == short[x] {
                sub_cost = prev[x];
            } else {
                sub_cost = prev[x] + 1;
            }

            // insert the minimum cost into the array
            cur[x+1] = *[del_cost, ins_cost, sub_cost].iter().min().unwrap();
        }

        // move the current vector to the prev location so that it can be looked at next iteration
        std::mem::swap(&mut cur, &mut prev);
    }

    return prev[shortd];
    
}


/// Retrieves the size of a file in bytes
pub fn file_size(dir: &str) -> Result<u64, PhyloError> {
    Ok(fs::metadata(dir).map_err(|_| PhyloError::FileReadError(String::from(dir)))?.len())
}


/// Generates k-mers (n-grams) for the given file
pub fn generate_kmers(file_dir: &str, k: u32, num: u32) -> Result<Vec<String>, PhyloError> {
    let size = file_size(file_dir)? as usize;
    let mut ret = Vec::with_capacity(size);     // create a vector with the right size
    let mut loc: Vec<usize> = Vec::with_capacity(num as usize);     // decide all locations to generate kmers at
    let mut rng = rand::thread_rng();

    // Ensure the file is sized satisfactorily
    if size < k as usize {
        return Err(PhyloError::FileTooSmall(String::from(file_dir)));
    }

    // Find all locations we'll take kmers at
    for _ in 0..num {
        loc.push(rng.gen_range(0..(size - k as usize - 1)));
    }
    loc.sort();

    // Prepare the file and buffer
    let file = File::open(file_dir).map_err(|_| PhyloError::FileOpenError(String::from(file_dir)))?; // Open the genome file
    let mut buffer: Vec<u8> = vec![0; k.try_into().map_err(|_| PhyloError::KTooBig(k))?];

    // For each kmer
    for i in 0..num {
        file.read_exact_at(&mut buffer, loc[i as usize] as u64).map_err(|_| PhyloError::FileReadError(String::from(file_dir)))?; // Read from the specified location
        ret.push(String::from_utf8(buffer.clone()).map_err(|_| PhyloError::FileReadError(String::from(file_dir)))?);    // Push to return vector
    }

    Ok(ret)
}

