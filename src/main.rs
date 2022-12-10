mod algorithms;
mod errors;
mod structs;

fn main() {
    let result = algorithms::levenshtein("kit", "glimmen");
    println!("{}", result);
    println!("{}", algorithms::levenshtein("kitten", "alderkitten"));
    println!("{}", algorithms::levenshtein("alderkitten", "kitten"));
    dbg!(algorithms::file_size("/home/terrior/Programming/genome-tree/src/test.txt"));
}
