use std::{fs::{self, File}, io::Read, os::unix::prelude::FileExt, collections::HashMap};
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


/// Given a list of numbers and weights, choose a random element; if limitless is false then the probability acts like a limit
pub fn random_weighted(elems: Vec<u32>, probabilities: Vec<u32>, rounds: u32, limitless: bool) -> Vec<u32> {

    // return if input is invalid
    if elems.len() != probabilities.len() {
        return Vec::new();
    }

    // preparing
    let mut total: u32 = 0;
    let mut new_probabilities: Vec<u32> = Vec::new();
    let mut amounts: Vec<u32> = vec![0; elems.len()];


    // first total up the weights
    for prob in &probabilities {
        total += prob;
        new_probabilities.push(total); //tally up all the probabilities up until now
    }

    // return if input is invalid
    if total < rounds {
        return Vec::new();
    }
    let mut rng = rand::thread_rng();
    let mut ret: Vec<u32> = Vec::new();

    // for the number of items we want to generate
    for _ in 0..rounds {
        'limitless_iter: loop {
            let gen = rng.gen_range(0..total); //generate the new random number

            // for every item
            for i in 0..new_probabilities.len() {
                if gen < new_probabilities[i] { //we found the index

                    // if we're limiting the amounts and we've reached our limit, generate again
                    if !limitless && amounts[i] >= probabilities[i] {
                        continue 'limitless_iter;
                    }

                    ret.push(elems[i]); //push the selected element to the return vector
                    amounts[1] += 1;
                    break;
                }
            }
            break;
        }           
    }
    
    ret
}


/// Convert a vector of items into a hashmap where every key is tied to the number of its occurrences in the vector
pub fn vec_to_dict(elems: Vec<u32>) -> HashMap<u32, u32> {
    let mut ret = HashMap::new();

    for elem in elems {
        if ret.contains_key(&elem) {
            *ret.get_mut(&elem).unwrap() += 1;
        } else {
            ret.insert(elem, 0);
        }
    }
    
    ret
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

