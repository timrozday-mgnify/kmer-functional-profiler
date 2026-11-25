//! Every query's best target outside its own group (genome): best.tsv of qid, score, ties, targets (up to 50).
use crate::data::*;
use crate::sw;
use std::fmt::Write as _;
use std::path::Path;
use std::sync::Mutex;

pub fn run(dir: &Path, threads: usize) {
    let (q, t) = (Seqs::load(dir, "q"), Seqs::load(dir, "reps"));
    let (qg, tg) = (u32s(&dir.join("q.group")), u32s(&dir.join("reps.group")));
    let mut order: Vec<u32> = (0..t.len() as u32).collect();
    order.sort_unstable_by_key(|&i| t.off[i as usize + 1] - t.off[i as usize]);
    let batches: Vec<&[u32]> = order.chunks(sw::LANES).collect();
    let qs: Vec<&[u8]> = (0..q.len()).map(|i| q.get(i)).collect();
    let best: Vec<Mutex<(i16, usize, Vec<u32>)>> = qs.iter().map(|_| Mutex::new((1, 0, vec![]))).collect();
    parallel(batches.len(), threads, |b| {
        let ts: Vec<&[u8]> = batches[b].iter().map(|&i| t.get(i as usize)).collect();
        let mut out = vec![[0i16; sw::LANES]; qs.len()];
        for (qb, ob) in qs.chunks(48).zip(out.chunks_mut(48)) {
            sw::batch(&ts, qb, ob);
        }
        for (i, o) in out.iter().enumerate() {
            for (k, &u) in batches[b].iter().enumerate() {
                if tg[u as usize] == qg[i] {
                    continue;
                }
                let mut m = best[i].lock().unwrap();
                if o[k] > m.0 {
                    *m = (o[k], 1, vec![u]);
                } else if o[k] == m.0 {
                    m.1 += 1;
                    if m.2.len() < 50 {
                        m.2.push(u);
                    }
                }
            }
        }
    });
    let mut s = String::from("qid\tscore\tties\ttargets\n");
    for (i, m) in best.iter().enumerate() {
        let m = m.lock().unwrap();
        let l: Vec<String> = m.2.iter().map(|u| u.to_string()).collect();
        writeln!(s, "{i}\t{}\t{}\t{}", m.0, m.1, l.join(",")).unwrap();
    }
    std::fs::write(dir.join("best.tsv"), s).unwrap();
}
