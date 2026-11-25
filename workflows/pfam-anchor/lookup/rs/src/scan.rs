//! The residue scan end to end, one thread: a family's GA representative rows by column, the read's residues scored in
//! quarter bits, the top K rows by histogram, then columns (BLOSUM62) -> top 50 -> banded gapped -> the call.
use crate::cand::fam;
use crate::data::*;
use crate::sw;
use std::fmt::Write as _;
use std::path::Path;
use std::time::Instant;

const TILE: usize = 16384;

pub fn kadane(it: impl Iterator<Item = i16>) -> i16 {
    let (mut cur, mut best) = (0i16, 0i16);
    for s in it {
        cur = (cur + s).max(0);
        best = best.max(cur);
    }
    best
}

/// Scores of every row: quarter bits of the residues that match. `cols`: (the column's bytes, residue, weight).
pub fn score_rows(cols: &[(&[u8], u8, u16)], out: &mut [u16]) {
    out.fill(0);
    for lo in (0..out.len()).step_by(TILE) {
        let hi = (lo + TILE).min(out.len());
        let acc = &mut out[lo..hi];
        for &(c, q, w) in cols {
            // q and w are copies: read through `cols` they stop the loop vectorising (1.14 ns a row and column)
            for (x, &v) in acc.iter_mut().zip(&c[lo..hi]) {
                *x = x.wrapping_add(if v == q { w } else { 0 });
            }
        }
    }
}

/// The score of the k-th best row (1 at least, so that a row without a match is never kept).
pub fn bar(scores: &[u16], k: usize) -> u16 {
    let top = scores.iter().copied().max().unwrap_or(0) as usize;
    let mut h = vec![0u32; top + 1];
    scores.iter().for_each(|&s| h[s as usize] += 1);
    let mut n = 0;
    for s in (1..=top).rev() {
        n += h[s] as usize;
        if n >= k {
            return s as u16;
        }
    }
    1
}

/// `flank`: residues either side of the domain for the banded score, or the whole protein.
/// `configs`: K:SECOND pairs; `subs`: keep 1 in 2^sub units; `limit`: the first queries only (0: all).
/// scan.tsv: sub, k, second, qid, score, ties, units. scan_times.tsv: core-seconds by stage, summed over threads.
pub fn run(dir: &Path, out: &Path, configs: &[(usize, usize)], subs: &[u32], threads: usize, limit: usize, weak: bool, flank: Option<usize>) {
    let fams = families(dir);
    let q = Queries::load(dir, &fams);
    let rows = Rows::load(dir, "reps", &fams);
    let units = family_units(&fams, &[&rows]);
    let ga = if weak { vec![1u8; rows.unit.len()] } else { bytes(&dir.join("reps.ga")) };
    let t = Seqs::load(dir, "reps");
    let bl = blosum();
    let nq = if limit > 0 { limit.min(q.len()) } else { q.len() };
    let mut o = String::from("sub\tk\tsecond\tqid\tscore\tties\tunits\n");
    let mut times = String::from("sub\tk\tsecond\tthreads\treads\trows_read\tcells\trows_kept\tscan_s\ttop_s\tcolumns_s\tgapped_s\twall_s\n");
    for &sub in subs {
        let mut sums = vec![[0f64; 8]; configs.len()]; // rows, cells, kept, scan, top, columns, gapped, wall
        for fi in 0..fams.len() {
            let mine: Vec<usize> = (0..nq).filter(|&i| q.fam[i] as usize == fi).collect();
            if mine.is_empty() {
                continue;
            }
            let f = fam(&fams[fi], fi, &rows, Some(&ga), &units[fi], sub);
            let n = f.rows.len();
            if n == 0 {
                continue;
            }
            let mut m = vec![NONE; n * f.len];
            for (r, row) in f.rows.iter().enumerate() {
                for (c, &a) in row.iter().enumerate() {
                    m[c * n + r] = a;
                }
            }
            for (ci, &(k, second)) in configs.iter().enumerate() {
                let done = std::sync::Mutex::new((vec![], [0f64; 8]));
                let next = std::sync::atomic::AtomicUsize::new(0);
                let wall = Instant::now();
                parallel(threads, threads, |_| {
                    let mut s = vec![0u16; n];
                    let (mut buf, mut lines, mut a) = (vec![], vec![], [0f64; 8]);
                    loop {
                        let x = next.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                        let Some(&i) = mine.get(x) else { break };
                        let (row, pos, qs) = (q.row(i), q.rowpos(i), q.seqs.get(i));
                        let cols: Vec<usize> = (0..f.len).filter(|&c| row[c] < 20).collect();
                        let t0 = Instant::now();
                        let cw: Vec<(&[u8], u8, u16)> = cols.iter()
                            .map(|&c| (&m[c * n..(c + 1) * n], row[c], (4.0 * f.bits[c][row[c] as usize]).round() as u16))
                            .collect();
                        score_rows(&cw, &mut s);
                        let t1 = Instant::now();
                        let b = bar(&s, k.min(n));
                        let cand: Vec<u32> = (0..n as u32).filter(|&r| s[r as usize] >= b).collect();
                        let t2 = Instant::now();
                        let mut v: Vec<(i16, u32)> = cand.iter().map(|&r| {
                            let x = f.rows[r as usize];
                            let sc = |c: &usize| if x[*c] < 21 { bl[row[*c] as usize][x[*c] as usize] } else { -4 };
                            (kadane(cols.iter().map(sc)), r)
                        }).collect();
                        v.sort_unstable_by(|a, b| b.cmp(a));
                        let cut = v.get(second.min(v.len()).saturating_sub(1)).map_or(0, |x| x.0);
                        let t3 = Instant::now();
                        let mut best: (i16, Vec<u32>) = (0, vec![]);
                        for &(_, r) in v.iter().take_while(|x| x.0 >= cut) {
                            let (x, tp) = (f.rows[r as usize], f.pos[r as usize]);
                            let mut d: Vec<i32> = cols.iter().filter(|&&c| x[c] == row[c]).map(|&c| tp[c] as i32 - pos[c] as i32).collect();
                            d.sort_unstable();
                            let Some(&d) = d.get(d.len() / 2) else { continue };
                            let u = f.units[f.slot[r as usize] as usize];
                            let (tt, d) = match flank {
                                None => (t.get(u as usize), d),
                                // only the domain's own residues, inserts included, as an index without proteins holds
                                Some(w) => {
                                    let on = || tp.iter().copied().filter(|&p| p != u16::MAX);
                                    let (a, b) = (on().min().unwrap() as usize, on().max().unwrap() as usize + 1);
                                    let p = t.get(u as usize);
                                    let a = a.saturating_sub(w);
                                    (&p[a..(b + w).min(p.len())], d - a as i32)
                                }
                            };
                            let lo = (d - 8).clamp(0, tt.len() as i32);
                            let hi = (qs.len() as i32 + d + 8).clamp(lo, tt.len() as i32);
                            let sc = sw::score(qs, &tt[lo as usize..hi as usize], &bl, Some((d - 8 - lo, d + 8 - lo)), &mut buf);
                            if sc > best.0 {
                                best = (sc, vec![u]);
                            } else if sc == best.0 && !best.1.contains(&u) {
                                best.1.push(u);
                            }
                        }
                        let t4 = Instant::now();
                        let d = |a: Instant, b: Instant| (b - a).as_secs_f64();
                        for (x, y) in a.iter_mut().zip([n as f64, (n * cols.len()) as f64, cand.len() as f64, d(t0, t1), d(t1, t2), d(t2, t3), d(t3, t4), 0.0]) {
                            *x += y;
                        }
                        let l: Vec<String> = best.1.iter().take(200).map(|u| u.to_string()).collect();
                        lines.push((i, format!("{sub}\t{k}\t{second}\t{i}\t{}\t{}\t{}\n", best.0, best.1.len(), l.join(","))));
                    }
                    let mut g = done.lock().unwrap();
                    g.0.extend(lines);
                    g.1.iter_mut().zip(a).for_each(|(x, y)| *x += y);
                });
                let wall = wall.elapsed().as_secs_f64();
                let (mut lines, a) = done.into_inner().unwrap();
                lines.sort();
                lines.iter().for_each(|x| o.push_str(&x.1));
                sums[ci].iter_mut().zip(a).for_each(|(x, y)| *x += y);
                sums[ci][7] += wall;
            }
        }
        for (&(k, second), a) in configs.iter().zip(&sums) {
            writeln!(times, "{sub}\t{k}\t{second}\t{threads}\t{nq}\t{:.0}\t{:.0}\t{:.0}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}", a[0], a[1], a[2], a[3], a[4], a[5], a[6], a[7]).unwrap();
            eprintln!("sub {sub} K {k} second {second} threads {threads}: scan {:.3} + top {:.3} ns a row, {:.1} GB/s a thread; columns {:.0} + gapped {:.0} us a read",
                a[3] * 1e9 / a[0], a[4] * 1e9 / a[0], a[1] / a[3] / 1e9, a[5] * 1e6 / nq as f64, a[6] * 1e6 / nq as f64);
        }
    }
    std::fs::write(out, o).unwrap();
    std::fs::write(out.with_extension("times.tsv"), times).unwrap();
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn scores_and_bar() {
        let n = 2 * TILE + 5;
        let cols: Vec<Vec<u8>> = (0..7u64).map(|c| (0..n as u64).map(|r| (crate::cand::hash(r * 7 + c) % 21) as u8).collect()).collect();
        let w: Vec<(&[u8], u8, u16)> = cols.iter().enumerate().map(|(i, c)| (&c[..], i as u8, 3 + i as u16)).collect();
        let mut s = vec![0; n];
        score_rows(&w, &mut s);
        for r in 0..n {
            assert_eq!(s[r], w.iter().filter(|x| x.0[r] == x.1).map(|x| x.2).sum::<u16>());
        }
        let mut v = s.clone();
        v.sort_unstable_by(|a, b| b.cmp(a));
        assert_eq!(bar(&s, 500), v[499].max(1));
    }
}
