pub fn levenshtein_distance(first: &str, second: &str) -> i32 {
    let long: &str;
    let short: &str;

    // first figure out the longest and shortest strings
    if first.len() > second.len() {
        long = first;
        short = second;
    } else {
        short = first;
        long = second;
    }
    
    let matrix: [Vec<i32>; 3] = [Vec::with_capacity(short.len()), Vec::with_capacity(short.len()), Vec::with_capacity(short.len())];
    return 0;
    
}
