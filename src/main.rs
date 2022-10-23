mod levenshtein;


fn main() {
    println!("Hello, world!");
    levenshtein::levenshtein_distance("kitten", "sitting");
}
