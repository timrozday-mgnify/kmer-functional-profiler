//! Second stage: what to do with a first stage's candidates before the alignment that scores them.
//!
//!   window    best local BLOSUM62 run over the columns within 8 of a seed (seeds only)
//!   columns   the same over every column of the read: the anchored full diagonal, HMM gaps allowed
//!   diagonal  ungapped, in sequence: the seed's diagonal over the whole ORF and protein
//!   gapped    banded Smith-Waterman, 8 either side of that diagonal
//!   sw        Smith-Waterman, the score of the call
//!
//! follow.tsv: option, stage, rule, qid, rows scored, units kept, best and near as in cand.rs.
//! calls.tsv: option, qid, best sw score among candidates, units tied at it (up to 200).
use crate::cand::*;
use crate::data::*;
use crate::sw;
use crate::truth::fam_offsets;
use std::collections::HashMap;
use std::fmt::Write as _;
use std::path::Path;
use std::sync::Mutex;

const STAGES: [&str; 5] = ["window", "columns", "diagonal", "gapped", "sw"];
const GAP: i16 = -4; // a column where the target has no residue

fn kadane(it: impl Iterator<Item = i16>) -> i16 {
    let (mut cur, mut best) = (0i16, 0i16);
    for s in it {
        cur = (cur + s).max(0);
        best = best.max(cur);
    }
    best
}

/// The first stage to run: a seed (with the hits wanted), or the top units by bits over every column.
enum First {
    Seed(Seed, u16),
    Top(usize),
}
fn first(o: &str) -> First {
    let (s, h) = o.split_once("/h").map_or((o, 1), |x| (x.0, x.1.parse().unwrap()));
    if let Some(k) = s.strip_prefix("top") {
        First::Top(k.parse().unwrap())
    } else if let Some(t) = s.strip_prefix("bits") {
        First::Seed(Seed::Bits(t.parse().unwrap()), h)
    } else if let Some(t) = s.strip_prefix("ruler") {
        First::Seed(Seed::Ruler(t.parse().unwrap()), h)
    } else {
        First::Seed(Seed::Shape(shape(s)), h) // 1100101, say
    }
}

pub fn run(dir: &Path, options: &[&str], out: &Path, threads: usize, truth_arms: &[&str]) {
    let fams = families(dir);
    let q = Queries::load(dir, &fams);
    let units = {
        let r: Vec<Rows> = truth_arms.iter().map(|a| Rows::load(dir, a, &fams)).collect();
        family_units(&fams, &r.iter().collect::<Vec<_>>())
    };
    let rows = Rows::load(dir, "reps", &fams);
    let t = Seqs::load(dir, "reps");
    let off = fam_offsets(&q, &units);
    let truth: Vec<i16> = u16s(&dir.join("truth.fam")).into_iter().map(|x| x as i16).collect();
    let tsv = Mutex::new(String::from("option\tstage\trule\tqid\trows\tkept\tbest\tnear\n"));
    let calls = Mutex::new(String::from("option\tqid\tscore\tdiag_equal\tgapped_equal\tcands\tties\tident\tcols\tunits\n"));
    let secs = Mutex::new(HashMap::<String, (f64, usize)>::new());
    let bl = blosum();
    parallel(fams.len(), threads, |fi| {
        let f = fam(&fams[fi], fi, &rows, None, &units[fi], 0);
        let mine: Vec<usize> = (0..q.len()).filter(|&i| q.fam[i] as usize == fi).collect();
        if f.units.is_empty() || mine.is_empty() {
            return;
        }
        let (mut o, mut c_o) = (String::new(), String::new());
        let mut time: HashMap<String, (f64, usize)> = HashMap::new();
        for opt in options {
            let fs = first(opt);
            // row-level index: (key, row, column)
            let mut ix: Vec<(u64, u32, u16)> = vec![];
            if let First::Seed(seed, _) = &fs {
                let mut ks = vec![];
                for (r, row) in f.rows.iter().enumerate() {
                    ks.clear();
                    keys(seed, row, &f.bits, None, &mut ks);
                    ix.extend(ks.iter().map(|k| (k.0, r as u32, k.1)));
                }
                ix.sort_unstable();
            }
            let mut seeds: Vec<Vec<u16>> = vec![vec![]; f.rows.len()];
            for &i in &mine {
                let (row, pos, qs) = (q.row(i), q.rowpos(i), q.seqs.get(i));
                let best = f.uslot.iter().map(|&s| truth[off[i] + s as usize]).max().unwrap();
                // candidate rows
                let mut cand: Vec<u32> = vec![];
                match &fs {
                    First::Seed(seed, h) => {
                        let mut ks = vec![];
                        keys(seed, row, &f.bits, None, &mut ks);
                        for (k, _) in &ks {
                            let a = ix.partition_point(|p| p.0 < *k);
                            for p in ix[a..].iter().take_while(|p| p.0 == *k) {
                                if seeds[p.1 as usize].is_empty() {
                                    cand.push(p.1);
                                }
                                seeds[p.1 as usize].push(p.2);
                            }
                        }
                        cand.retain(|&r| {
                            let ok = seeds[r as usize].len() >= *h as usize;
                            if !ok {
                                seeds[r as usize].clear();
                            }
                            ok
                        });
                    }
                    First::Top(k) => {
                        let cols: Vec<usize> = (0..f.len).filter(|&c| row[c] < 20).collect();
                        let s: Vec<f32> = f.rows.iter()
                            .map(|r| cols.iter().filter(|&&c| r[c] == row[c]).map(|&c| f.bits[c][row[c] as usize]).sum())
                            .collect();
                        let mut u = vec![0f32; f.units.len()];
                        s.iter().zip(&f.slot).for_each(|(&x, &sl)| u[sl as usize] = u[sl as usize].max(x));
                        u.sort_by(|a, b| b.total_cmp(a));
                        let bar = u[(*k).min(u.len()) - 1];
                        cand.extend((0..s.len() as u32).filter(|&r| s[r as usize] > 0.0 && s[r as usize] >= bar));
                    }
                }
                let cols: Vec<usize> = (0..f.len).filter(|&c| row[c] < 21).collect();
                let mut st: [Vec<i16>; 5] = Default::default();
                let mut clock = |name: &str, n: usize, g: &mut dyn FnMut()| {
                    let t0 = std::time::Instant::now();
                    g();
                    let e = time.entry(name.to_string()).or_default();
                    e.0 += t0.elapsed().as_secs_f64();
                    e.1 += n;
                };
                let score_col = |r: &[u8], c: usize| if r[c] < 21 { bl[row[c] as usize][r[c] as usize] } else { GAP };
                if matches!(fs, First::Seed(..)) {
                    clock("window", cand.len(), &mut || {
                        st[0] = cand.iter().map(|&r| {
                            let t = f.rows[r as usize];
                            seeds[r as usize].iter().map(|&c| {
                                let (a, b) = ((c as usize).saturating_sub(8), (c as usize + 16).min(f.len));
                                kadane((a..b).filter(|&x| row[x] < 21).map(|x| score_col(t, x)))
                            }).max().unwrap()
                        }).collect();
                    });
                } else {
                    st[0] = vec![0; cand.len()];
                }
                clock("columns", cand.len(), &mut || {
                    st[1] = cand.iter().map(|&r| kadane(cols.iter().map(|&c| score_col(f.rows[r as usize], c)))).collect();
                });
                // each row's diagonal: the seed's, or the one most shared columns lie on
                let diag: Vec<i32> = cand.iter().map(|&r| {
                    let (t, tp) = (f.rows[r as usize], f.pos[r as usize]);
                    let d = |c: usize| tp[c] as i32 - pos[c] as i32;
                    match seeds[r as usize].get(seeds[r as usize].len() / 2) {
                        Some(&c) => d(c as usize),
                        None => {
                            let mut v: Vec<i32> = cols.iter().filter(|&&c| t[c] == row[c]).map(|&c| d(c)).collect();
                            v.sort_unstable();
                            v[v.len() / 2]
                        }
                    }
                }).collect();
                let target = |r: u32| t.get(f.units[f.slot[r as usize] as usize] as usize);
                clock("diagonal", cand.len(), &mut || {
                    st[2] = cand.iter().zip(&diag).map(|(&r, &d)| {
                        let tt = target(r);
                        let lo = (-d).max(0);
                        let hi = (qs.len() as i32).min(tt.len() as i32 - d);
                        kadane((lo..hi).map(|x| bl[qs[x as usize] as usize][tt[(x + d) as usize] as usize]))
                    }).collect();
                });
                let mut buf = vec![];
                clock("gapped", cand.len(), &mut || {
                    st[3] = cand.iter().zip(&diag).map(|(&r, &d)| {
                        // only the stretch of target the band can reach
                        let tt = target(r);
                        let a = (d - 8).clamp(0, tt.len() as i32);
                        let b = (qs.len() as i32 + d + 8).clamp(a, tt.len() as i32);
                        sw::score(qs, &tt[a as usize..b as usize], &bl, Some((d - 8 - a, d + 8 - a)), &mut buf)
                    }).collect();
                });
                clock("sw", cand.len(), &mut || {
                    st[4] = cand.iter().map(|&r| sw::score(qs, target(r), &bl, None, &mut buf)).collect();
                });
                // per unit, the best of its rows
                let mut per: HashMap<u32, [i16; 5]> = HashMap::new();
                for (k, &r) in cand.iter().enumerate() {
                    let e = per.entry(f.slot[r as usize]).or_insert([0; 5]);
                    for s in 0..5 {
                        e[s] = e[s].max(st[s][k]);
                    }
                }
                cand.iter().for_each(|&r| seeds[r as usize].clear());
                let per: Vec<[i16; 5]> = {
                    let mut v: Vec<_> = per.into_iter().collect();
                    v.sort_unstable();
                    let top = v.iter().map(|x| x.1[4]).max().unwrap_or(0);
                    let ties: Vec<u32> = v.iter().filter(|x| x.1[4] == top).map(|x| f.units[x.0 as usize]).collect();
                    let list: Vec<String> = ties.iter().take(200).map(|u| u.to_string()).collect();
                    let eq = |s: usize| v.iter().filter(|x| x.1[4] == top && x.1[s] == top).count();
                    let a = ties.first().map_or((0, 0, 0, 0), |&u| sw::align(qs, t.get(u as usize), &bl, None));
                    writeln!(c_o, "{opt}\t{i}\t{top}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", eq(2), eq(3), v.len(), ties.len(), a.1, a.2, list.join(",")).unwrap();
                    v.into_iter().map(|x| x.1).collect()
                };
                for (s, stage) in STAGES.iter().enumerate() {
                    if s == 0 && matches!(fs, First::Top(_)) {
                        continue;
                    }
                    let top = per.iter().map(|x| x[s]).max().unwrap_or(0);
                    let mut v: Vec<i16> = per.iter().map(|x| x[s]).collect();
                    v.sort_unstable_by(|a, b| b.cmp(a));
                    let mut emit = |rule: String, bar: i16| {
                        let kept: Vec<&[i16; 5]> = per.iter().filter(|x| x[s] >= bar).collect();
                        let b = kept.iter().any(|x| x[4] == best);
                        let n = kept.iter().any(|x| x[4] as f32 >= 0.95 * best as f32);
                        writeln!(o, "{opt}\t{stage}\t{rule}\t{i}\t{}\t{}\t{}\t{}", cand.len(), kept.len(), b as u8, n as u8).unwrap();
                    };
                    emit("all".into(), i16::MIN);
                    for d in [0, 5, 10, 20, 30] {
                        emit(format!("best-{d}"), top - d);
                    }
                    for k in [10usize, 50, 200] {
                        emit(format!("top{k}"), if v.is_empty() { 0 } else { v[k.min(v.len()) - 1] });
                    }
                    for k in [20, 30, 40, 50] {
                        emit(format!("at{k}"), k);
                    }
                }
            }
        }
        tsv.lock().unwrap().push_str(&o);
        calls.lock().unwrap().push_str(&c_o);
        let mut all = secs.lock().unwrap();
        for (k, v) in time {
            let e = all.entry(k).or_default();
            e.0 += v.0;
            e.1 += v.1;
        }
    });
    std::fs::create_dir_all(out).unwrap();
    std::fs::write(out.join("follow.tsv"), &*tsv.lock().unwrap()).unwrap();
    std::fs::write(out.join("calls.tsv"), &*calls.lock().unwrap()).unwrap();
    let mut s = String::from("stage\tseconds\trows\n");
    for (k, v) in secs.lock().unwrap().iter() {
        writeln!(s, "{k}\t{:.4}\t{}", v.0, v.1).unwrap();
    }
    std::fs::write(out.join("stage_seconds.tsv"), s).unwrap();
}
