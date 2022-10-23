mod levenshtein;

fn main() {
    let result = levenshtein::levenshtein_distance("kit", "glimmen");
    println!("{}", result);
}
