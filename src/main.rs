mod algorithms;

fn main() {
    let result = algorithms::levenshtein("kit", "glimmen");
    println!("{}", result);
    println!("{}", algorithms::levenshtein("kitten", "alderkitten"));
    println!("{}", algorithms::levenshtein("alderkitten", "kitten"));
}
