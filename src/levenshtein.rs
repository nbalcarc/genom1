pub fn levenshtein_distance(first: &str, second: &str) -> usize {
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

    // create all the vectors we'll need, all sized to the smaller string
    // let mut store = Vec::with_capacity(short.len());
    let mut last = vec![0; short.len()];
    let mut cur = vec![0; short.len()];

    // initialize all elements
    for i in 0..shortd {
        cur[i] = i;
    }

    // iterate once for every letter in the long word
    for y in 1..longd+1 {
        last.clone_from(&cur); //move cur to last
        cur[0] = y; //set the character index

        for x in 1..shortd+1 {

            let subcost: usize;
            if long[y-1] == short[x-1] {
                subcost = 0;
            } else {
                subcost = 1;
            }

            cur[x] = *[cur[x-1] + 1, last[x] + 1, last[x-1] + subcost].iter().min().unwrap();
        }
    }


    println!("{:?}", last);    
    println!("{:?}", cur);    

    return cur[shortd - 1];
    
}
