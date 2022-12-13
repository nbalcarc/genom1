use std::{fs::{self, File}, io::{Read, BufReader, BufRead}, os::unix::prelude::FileExt, collections::HashMap};
use rand::Rng;

use crate::{errors::PhyloError, structs::{Genome, TreeNode, TreeVertex}};

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

    // declare variables outside loop so we don't have to reallocate them
    let mut del_cost;
    let mut ins_cost;
    let mut sub_cost;

    // iterate once for every letter in the long word
    for y in 0..longd {
        //println!("{} of {}!", y, longd);
        cur[0] = y+1; //set the character index

        // iterate once for every letter in the short word (size of the arrays)
        for x in 0..shortd {
            del_cost = prev[x+1] + 1;   // generate cost of deletion
            ins_cost = cur[x] + 1;      // generate cost of insertion

            // generate cost of substitution
            if long[x+1] == short[x] {
                sub_cost = prev[x];
            } else {
                sub_cost = prev[x] + 1;
            }

            // insert the minimum cost into the array
            //cur[x+1] = *[del_cost, ins_cost, sub_cost].iter().min().unwrap();
            cur[x+1] = if del_cost < ins_cost && del_cost < sub_cost {
                del_cost
            } else if ins_cost < del_cost && ins_cost < sub_cost {
                ins_cost
            } else {
                sub_cost
            }
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


/// Check how many kmers apply to the given genome
pub fn kmer_similarity(host: &Genome, guest: &Genome) -> u32 {
    let mut ret = 0;
    let kmers = host.kmers.clone();
   
    // open file
    let mut buffer: Vec<u8> = vec![0; 1024];
    let file = File::open(&guest.dir).unwrap();
    let mut reader = BufReader::new(file);

    // read all data
    let mut buff = Vec::new();
    reader.read_to_end(&mut buff);
    let all = String::from_utf8(buff).unwrap();

    // for every kmer
    for i in 0..kmers.len() {
        if all.contains(&kmers[i]) {
            ret += 1;
        }
    }

    ret
}


/// Retrieve a genome from the tree
pub fn retrieve_genome<'a>(root: &'a mut TreeNode, path: &Vec<u8>) -> Result<&'a mut Genome, PhyloError> {
    if root.id != path[0] {
        return Err(PhyloError::SearchGenomeError);
    }

    let mut paths = path.clone(); //prepare for looping
    let mut cur = root;
    paths.reverse();
    paths.remove(paths.len()-1);

    while paths.len() > 0 {
        match &mut cur.vertex {
            TreeVertex::Floor(f) => { //we hit a floor
                if paths.len() == 1 { //if we hit a floor and we only have one index left
                    return Ok(&mut f[paths[0] as usize]);
                } else { //something didn't match up
                    return Err(PhyloError::SearchGenomeError);
                }
            },
            TreeVertex::Split(s) => { //we hit a split
                for i in 0..s.len() {
                    if s[i].id == paths[paths.len()-1] { //if we found the next node in the path
                        cur = &mut s[i];
                        paths.remove(paths.len()-1);
                        break;
                    }
                    // if we found no node with the given id
                    return Err(PhyloError::SearchGenomeError);
                }
                return Err(PhyloError::SearchGenomeError);
            }
        }
    }

    Err(PhyloError::SearchGenomeError)
}


/// Retrieves all TreeNodes leading up to the path, excludes the final one
pub fn get_full_path<'a>(root: &'a TreeNode, path: &Vec<u8>) -> Result<Vec<&'a TreeNode>, PhyloError> {
    let mut ret = Vec::new();

    if root.id != path[0] { //return early if path already invalid
        return Err(PhyloError::SearchNodeError);
    }

    ret.push(root);

    // for each part of the path, find the node that corresponds to it and push it to the return vector
    let mut cur = root;
    'main_loop: for i in 1..path.len()-1 { //exclude the last index

        // find what type of vertex we're working with
        match &cur.vertex {
            TreeVertex::Split(s) => { //split found
                if i >= path.len()-2 { //too close to the edge
                    return Err(PhyloError::SearchNodeError);
                }
                // we aren't too close to the edge, find the next node with the given id
                // find the node to traverse next
                for node in s {
                    if node.id == path[i] { //found the next node
                        ret.push(node); //since we found a valid node, push it to ret
                        cur = node;
                        continue 'main_loop;
                    }
                    continue;
                }
                // if we reach here, then we couldn't find the node
                return Err(PhyloError::SearchNodeError);

            },
            TreeVertex::Floor(_) => { //floor found
                if i != path.len() - 2 { //if this isn't the second-to-last index
                    return Err(PhyloError::SearchNodeError);
                }
                // the for loop will boot us out into the return after this iteration
            }
        }
    }
    Ok(ret)
}


/// Return a mutable reference to a given node
pub fn get_mut_node_and_increment<'a>(root: &'a mut TreeNode, path: &Vec<u8>) -> Result<&'a mut TreeNode, PhyloError> {
    if root.id != path[0] { //return early if path already invalid
        return Err(PhyloError::SearchNodeError);
    }

    // for each part of the path, find the node that corresponds to it and push it to the return vector
    let mut cur = root;
    cur.count = cur.count + 1;
    'main_loop: for i in 1..path.len()-1 { //exclude the last index
        // find what type of vertex we're working with
        match  cur.vertex { // &mut cur.vertex didn't work, don't know why exactly
            // this just tales it by mutable ref, but I think that helps with your scoping problem
            // You were borrowing cur.vertex for the entire match expression, so you tried to return while it was
            // already dborrowed
            // wait so how did you make it not borrow the whole time, i have had this problem before
            // ref: make the thing I matched in this pattern a reference to the value
            // mut: make the thing I matched in this pattern mutable
            // -> ref mut: make the thing I matched in this pattern a mutable reference to the value
            // so why does it drop the reference to cur up at the top of the match?
            // I'm not 100% on it tbh
            // But basically, you were borrowing on line 305, that borrow didnt end till 349, and that was the problem
            // makes sense
            // I think since this line takes it by reference, it doesnt need to
            // I aslo think this would 
            TreeVertex::Split(ref mut s) => { //split found
                // What does this comment even mean? lol it means a split shouldn't be 1 before an index, it should be a floor
                if i >= path.len()-2 { //too close to the edge, split cannot exist where a floor should be right before an index
                    return Err(PhyloError::SearchNodeError);
                }
                // we aren't too close to the edge, find the next node with the given id
                // find the node to traverse next
                for node in s {
                    if node.id == path[i] { //found the next node
                        //std::mem::replace(cur, node); // Oh this is *not* what you want, I'm almost positve
                        cur = node; // This was a mutable ref to a mutable ref; not what you're looking for
                        cur.count = cur.count + 1; //increment the count
                        continue 'main_loop;
                    }
                    continue;
                }
                // if we reach here, then we couldn't find the node
                return Err(PhyloError::SearchNodeError);

            },
            TreeVertex::Floor(_) => { //floor found
                if i != path.len() - 2 { //if this isn't the second-to-last index
                    return Err(PhyloError::SearchNodeError);
                }
                // the for loop will boot us out into the return after this iteration
                // THIS IS BAD
                // DO YOU NOT REMEMBER THE LAST PROBLEM WE HAD?
                // At the end of this block, curr.vertex will *still be borrowed*
                // That's a nono
                // That's aliasing baby
                //return Err(PhyloError::SearchNodeError); // This fixed it; I don't know what your errors mean, but something has to happen
                break 'main_loop;
            }
        }
    }
    Ok(cur) 
}




