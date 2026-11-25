mod bench;
mod best;
mod bundle;
mod cand;
mod collapse;
mod data;
mod follow;
mod large;
mod scan;
mod sw;
mod truth;
use std::path::Path;

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let dir = Path::new(&a[2]);
    if a[1] == "bench" {
        return bench::run(a[2].parse().unwrap()); // bench ROWS
    }
    if a[1] == "best" {
        return best::run(Path::new(&a[2]), 6); // best DIR
    }
    if a[1] == "collapse" {
        return collapse::run(Path::new(&a[2]), 4); // collapse DIR
    }
    let threads: usize = std::env::var("THREADS").map_or(6, |x| x.parse().unwrap());
    let arms: Vec<&str> = if dir.join("members.rows").exists() { vec!["reps", "members"] } else { vec!["reps"] };
    match a[1].as_str() {
        // truth DIR [FIRST_N_TARGETS]
        "truth" => truth::run(dir, &arms, threads, a.get(3).map_or(0, |x| x.parse().unwrap())),
        // cand DIR ARM SUB OUT [OPTIONS]: ARM reps, reps_ga or members; SUB keeps 1 in 2^SUB units
        "cand" => cand::run(
            dir, &a[3], a[4].parse().unwrap(), Path::new(&a[5]), threads, &arms, a.get(6).map_or("", |x| x),
        ),
        // follow DIR OUT OPTION...: e.g. 1100101/h1, bits14, top500
        "follow" => {
            let o: Vec<&str> = a[4..].iter().map(|x| &x[..]).collect();
            follow::run(dir, &o, Path::new(&a[3]), threads, &arms)
        }
        // scan DIR OUT.tsv K:SECOND[,K:SECOND...] [SUB,SUB...]: the method end to end. THREADS (1 for timing), LIMIT, WEAK=1, FLANK, MEMBERS=1
        "scan" => {
            let pair = |x: &str| x.split_once(':').map(|(a, b)| (a.parse().unwrap(), b.parse().unwrap())).unwrap();
            let configs: Vec<(usize, usize)> = a[4].split(',').map(pair).collect();
            let subs: Vec<u32> = a.get(5).map_or(vec![0], |x| x.split(',').map(|y| y.parse().unwrap()).collect());
            let limit = std::env::var("LIMIT").map_or(0, |x| x.parse().unwrap());
            scan::run(dir, Path::new(&a[3]), &configs, &subs, threads, limit, std::env::var("WEAK").is_ok(), std::env::var("FLANK").ok().map(|x| x.parse().unwrap()), std::env::var("MEMBERS").is_ok_and(|x| x == "1"))
        }
        // bundle-build OUT_DIR T < rows
        "bundle-build" => bundle::build(dir, a[3].parse().unwrap()),
        // bundle-query DIR QUERIES OUT K SECOND
        "bundle-query" => bundle::query(dir, Path::new(&a[3]), Path::new(&a[4]), a[5].parse().unwrap(), a[6].parse().unwrap(), threads),
        "truthsub" => large::truthsub(dir, threads), // truthsub DIR
        // pick IDS.u64 OUT_PREFIX CHUNKS < FASTA
        "pick" => large::pick(dir, &a[3], a[4].parse().unwrap()),
        x => panic!("unknown step {x}"),
    }
}
