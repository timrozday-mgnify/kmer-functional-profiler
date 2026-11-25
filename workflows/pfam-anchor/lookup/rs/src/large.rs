//! For one very large family (hpc/kfp-ablations/pfam-anchor-large).
//!   pick      FASTA on stdin: the records whose MGYP id is in a file of u64 ids -> N chunk FASTAs named by integer id
//!   truthsub  every query's best Smith-Waterman score over the targets, and over 1 in 4, 16 and 64 of them
use crate::cand::hash;
use crate::data::*;
use crate::sw;
use std::collections::HashSet;
use std::fmt::Write as _;
use std::fs::File;
use std::io::{BufRead, BufWriter, Write};
use std::path::Path;
use std::sync::Mutex;

pub const SUBS: [u32; 4] = [0, 2, 4, 6];
pub fn kept(unit: u32, sub: u32) -> bool {
    hash(unit as u64) & ((1 << sub) - 1) == 0
}

pub fn pick(ids: &Path, out: &str, chunks: usize) {
    let ids: HashSet<u64> = u64s(ids).into_iter().collect();
    let mut w: Vec<_> = (0..chunks)
        .map(|i| BufWriter::with_capacity(1 << 20, File::create(format!("{out}{i:03}.faa")).unwrap()))
        .collect();
    let (mut n, mut seen, mut on) = (0usize, 0u64, false);
    let mut r = std::io::BufReader::with_capacity(1 << 22, std::io::stdin().lock());
    let mut line = vec![];
    while r.read_until(b'\n', &mut line).expect("reading stdin") > 0 {
        if line[0] == b'>' {
            seen += 1;
            let name = line[1..].split(|c| c.is_ascii_whitespace()).next().unwrap();
            let digits = name.strip_prefix(b"MGYP").unwrap_or(name);
            let id: u64 = std::str::from_utf8(digits).unwrap().parse().unwrap_or_else(|_| panic!("not an MGYP name: {}", String::from_utf8_lossy(name)));
            on = ids.contains(&id);
            if on {
                n += 1;
                writeln!(w[n % chunks], ">{id}").unwrap();
            }
        } else if on {
            w[n % chunks].write_all(&line).unwrap();
        }
        line.clear();
    }
    w.iter_mut().for_each(|x| x.flush().unwrap());
    // fewer records than the release holds means a cut stream: run.sh checks `seen`
    println!("{seen}\t{n}\t{}", ids.len());
}

pub fn truthsub(dir: &Path, threads: usize) {
    let (q, t) = (Seqs::load(dir, "q"), Seqs::load(dir, "reps"));
    let mut order: Vec<u32> = (0..t.len() as u32).collect();
    order.sort_unstable_by_key(|&i| t.off[i as usize + 1] - t.off[i as usize]);
    let batches: Vec<&[u32]> = order.chunks(sw::LANES).collect();
    let qs: Vec<&[u8]> = (0..q.len()).map(|i| q.get(i)).collect();
    let all = Mutex::new(vec![[(0i16, 0u32); 4]; qs.len()]);
    let t0 = std::time::Instant::now();
    let next = std::sync::atomic::AtomicUsize::new(0);
    parallel(threads, threads, |_| {
        let mut mine = vec![[(0i16, 0u32); 4]; qs.len()];
        let mut out = vec![[0i16; sw::LANES]; qs.len()];
        loop {
            let b = next.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
            if b >= batches.len() {
                break;
            }
            let ts: Vec<&[u8]> = batches[b].iter().map(|&i| t.get(i as usize)).collect();
            for (qb, ob) in qs.chunks(48).zip(out.chunks_mut(48)) {
                sw::batch(&ts, qb, ob);
            }
            for (k, &u) in batches[b].iter().enumerate() {
                for (s, &sub) in SUBS.iter().enumerate() {
                    if !kept(u, sub) {
                        continue;
                    }
                    for (m, o) in mine.iter_mut().zip(&out) {
                        if o[k] > m[s].0 {
                            m[s] = (o[k], 1);
                        } else if o[k] == m[s].0 {
                            m[s].1 += 1;
                        }
                    }
                }
            }
            if b % 50_000 == 0 {
                eprintln!("batch {b}/{} {:.0} s", batches.len(), t0.elapsed().as_secs_f64());
            }
        }
        for (a, m) in all.lock().unwrap().iter_mut().zip(mine) {
            for s in 0..4 {
                if m[s].0 > a[s].0 {
                    a[s] = m[s];
                } else if m[s].0 == a[s].0 {
                    a[s].1 += m[s].1;
                }
            }
        }
    });
    let mut o = String::from("sub\tqid\tbest\tties\n");
    for (i, a) in all.into_inner().unwrap().iter().enumerate() {
        for (s, sub) in SUBS.iter().enumerate() {
            writeln!(o, "{sub}\t{i}\t{}\t{}", a[s].0, a[s].1).unwrap();
        }
    }
    std::fs::write(dir.join("truthsub.tsv"), o).unwrap();
    let cells = qs.iter().map(|x| x.len() as f64).sum::<f64>() * t.seq.len() as f64;
    let s = t0.elapsed().as_secs_f64();
    eprintln!("{} queries, {} targets, {s:.0} s, {:.2} GCUPS a thread", qs.len(), t.len(), cells / s / 1e9 / threads as f64);
}
