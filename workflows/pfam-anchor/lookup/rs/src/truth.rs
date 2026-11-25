//! Exhaustive truth: every query against every MGnify90 representative, Smith-Waterman.
//!
//! truth.best   per query: best score over all targets (i16), and how many reach it (u32)
//! truth.fam    per query: the score of each of its family's units (i16), in `family_units` order
use crate::data::*;
use crate::sw;
use std::path::Path;
use std::sync::atomic::{AtomicI16, AtomicUsize, Ordering::Relaxed};
use std::sync::Mutex;

pub fn fam_offsets(q: &Queries, units: &[Vec<u32>]) -> Vec<usize> {
    let mut o = vec![0];
    for &f in &q.fam {
        o.push(o.last().unwrap() + units[f as usize].len());
    }
    o
}

pub fn run(dir: &Path, arms: &[&str], threads: usize, limit: usize) {
    let fams = families(dir);
    let q = Queries::load(dir, &fams);
    let rows: Vec<Rows> = arms.iter().map(|a| Rows::load(dir, a, &fams)).collect();
    let units = family_units(&fams, &rows.iter().collect::<Vec<_>>());
    drop(rows);
    let t = Seqs::load(dir, "reps");
    let off = fam_offsets(&q, &units);
    let nt = if limit > 0 { limit } else { t.len() };
    // slot[f][unit]: its place in units[f]
    let slot: Vec<Vec<u32>> = units
        .iter()
        .map(|u| {
            let mut s = vec![u32::MAX; t.len()];
            u.iter().enumerate().for_each(|(i, &x)| s[x as usize] = i as u32);
            s
        })
        .collect();
    let fam_scores: Vec<AtomicI16> = (0..*off.last().unwrap()).map(|_| AtomicI16::new(0)).collect();
    // target order: by length, so that a batch pads little
    let mut order: Vec<u32> = (0..nt as u32).collect();
    order.sort_unstable_by_key(|&i| t.off[i as usize + 1] - t.off[i as usize]);
    let batches: Vec<&[u32]> = order.chunks(sw::LANES).collect();
    let nq = q.len();
    let qs: Vec<&[u8]> = (0..nq).map(|i| q.seqs.get(i)).collect();
    let best = Mutex::new(vec![(0i16, 0u32); nq]);
    let done = AtomicUsize::new(0);
    let cells: f64 = qs.iter().map(|x| x.len()).sum::<usize>() as f64
        * batches.iter().map(|b| b.iter().map(|&i| t.get(i as usize).len()).sum::<usize>()).sum::<usize>() as f64;
    eprintln!("{nq} queries x {nt} targets, {:.2e} cells, {threads} threads", cells);
    let t0 = std::time::Instant::now();
    parallel(threads, threads, |_| {
        let mut mine = vec![(0i16, 0u32); nq];
        let mut out = vec![[0i16; sw::LANES]; nq];
        loop {
            let b = done.fetch_add(1, Relaxed);
            if b >= batches.len() {
                break;
            }
            let ids = batches[b];
            let ts: Vec<&[u8]> = ids.iter().map(|&i| t.get(i as usize)).collect();
            for (qb, ob) in qs.chunks(48).zip(out.chunks_mut(48)) {
                sw::batch(&ts, qb, ob);
            }
            for (i, o) in out.iter().enumerate() {
                let (f, at) = (q.fam[i] as usize, off[i]);
                for (k, &u) in ids.iter().enumerate() {
                    let s = o[k];
                    if s > mine[i].0 {
                        mine[i] = (s, 1);
                    } else if s == mine[i].0 {
                        mine[i].1 += 1;
                    }
                    let sl = slot[f][u as usize];
                    if sl != u32::MAX {
                        fam_scores[at + sl as usize].store(s, Relaxed);
                    }
                }
            }
            if b % 2000 == 0 {
                eprintln!("batch {b}/{} {:.0} s", batches.len(), t0.elapsed().as_secs_f64());
            }
        }
        let mut all = best.lock().unwrap();
        for (a, m) in all.iter_mut().zip(mine) {
            if m.0 > a.0 {
                *a = m;
            } else if m.0 == a.0 {
                a.1 += m.1;
            }
        }
    });
    let s = t0.elapsed().as_secs_f64();
    eprintln!("{s:.0} s, {:.1} GCUPS in all, {:.2} a thread", cells / s / 1e9, cells / s / 1e9 / threads as f64);
    let mut o = vec![];
    for (b, n) in best.into_inner().unwrap() {
        o.extend_from_slice(&b.to_le_bytes());
        o.extend_from_slice(&n.to_le_bytes());
    }
    std::fs::write(dir.join("truth.best"), o).unwrap();
    let mut o = Vec::with_capacity(2 * fam_scores.len());
    for x in &fam_scores {
        o.extend_from_slice(&x.load(Relaxed).to_le_bytes());
    }
    std::fs::write(dir.join("truth.fam"), o).unwrap();
}
