//! First stage: which units of the read's family to look at. Every option, scored against truth.
//!
//! Output (TSV, one row per option and query): option, qid, hits (postings touched), cands (units),
//! best (a unit with the family's best exhaustive score is among them), near (one within 5% of it),
//! tied (share of the units at the best score that are among them), src (the source cluster is; simulated
//! reads only). Index sizes go to index.tsv.
use crate::data::*;
use crate::truth::fam_offsets;
use std::collections::HashMap;
use std::fmt::Write as _;
use std::path::Path;
use std::sync::Mutex;

/// One family's rows of one arm, for the units kept.
pub struct Fam<'a> {
    pub len: usize,
    pub rows: Vec<&'a [u8]>,
    pub pos: Vec<&'a [u16]>,
    pub slot: Vec<u32>,      // unit slot of each row
    pub units: Vec<u32>,     // unit of each slot
    pub uslot: Vec<u32>,     // its place in the truth's list
    pub bits: Vec<[f32; AA]>, // -log2 of the residue's share among the rows, by column
    pub entropy: Vec<f32>,
}

pub fn hash(mut x: u64) -> u64 {
    x = (x ^ (x >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
    x = (x ^ (x >> 27)).wrapping_mul(0x94d049bb133111eb);
    x ^ (x >> 31)
}

pub fn fam<'a>(f: &Family, i: usize, r: &'a Rows, ga: Option<&[u8]>, truth_units: &[u32], sub: u32) -> Fam<'a> {
    let keep = |u: u32| crate::large::kept(u, sub);
    let (mut rows, mut pos, mut unit) = (vec![], vec![], vec![]);
    for x in r.row0[i]..r.row0[i + 1] {
        if keep(r.unit[x]) && ga.is_none_or(|g| g[x] == 1) {
            let b = r.byte0[i] + (x - r.row0[i]) * f.len;
            rows.push(&r.codes[b..b + f.len]);
            pos.push(&r.pos[b..b + f.len]);
            unit.push(r.unit[x]);
        }
    }
    let mut units = unit.clone();
    units.sort_unstable();
    units.dedup();
    let slot = unit.iter().map(|u| units.binary_search(u).unwrap() as u32).collect();
    let uslot = units.iter().map(|u| truth_units.binary_search(u).unwrap() as u32).collect();
    let mut count = vec![[0u32; AA]; f.len]; // not f32, which stops counting at 16.8 M
    for x in &rows {
        for (c, &a) in x.iter().enumerate() {
            count[c][a as usize] += 1;
        }
    }
    let mut bits = vec![[0f32; AA]; f.len];
    let mut entropy = vec![0f32; f.len];
    for c in 0..f.len {
        let n = count[c][..20].iter().sum::<u32>() as f32;
        for a in 0..20 {
            let p = (count[c][a] as f32 + 0.05) / (n + 1.0);
            bits[c][a] = -p.log2();
            entropy[c] -= p * p.log2();
        }
    }
    Fam { len: f.len, rows, pos, slot, units, uslot, bits, entropy }
}

// ---------------------------------------------------------------- seeds

#[derive(Clone)]
pub enum Seed {
    /// Residues at fixed offsets.
    Shape(Vec<usize>),
    /// Consecutive offsets until the residues' bits reach T.
    Bits(f32),
    /// Offsets of a ruler, as many as reach T.
    Ruler(f32),
    /// `Bits`, at a T per row: the highest of LEVELS with an expected -ln P(miss) of theta.
    Target(f32),
}
pub const RULER: [usize; 8] = [0, 1, 4, 9, 15, 22, 32, 34];
pub const LEVELS: [f32; 8] = [6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 20.0, 24.0];
const CAP: usize = 14; // residues in a seed, at most
const THETA_ID: f32 = 0.8; // identity at which a target's recall is estimated

fn seed_len(row: &[u8], bits: &[[f32; AA]], c: usize, t: f32) -> Option<usize> {
    let mut b = 0.0;
    for k in 0..CAP {
        let a = *row.get(c + k)?;
        if a >= 20 {
            return None;
        }
        b += bits[c + k][a as usize];
        if b >= t {
            return Some(k + 1);
        }
    }
    Some(CAP) // a conserved stretch: the cap, short of T
}

fn key_of(row: &[u8], c: usize, offs: impl Iterator<Item = usize>, salt: u64) -> Option<u64> {
    let mut x = hash(c as u64 * 64 + salt);
    for o in offs {
        let a = *row.get(c + o)?;
        if a >= 20 {
            return None;
        }
        x = hash(x ^ (a as u64 + 1));
    }
    Some(x)
}

/// Every key of a row (key, start column). `level`: for `Target`, the row's level or all of them.
pub fn keys(seed: &Seed, row: &[u8], bits: &[[f32; AA]], level: Option<usize>, out: &mut Vec<(u64, u16)>) {
    for c in 0..row.len() {
        match seed {
            Seed::Shape(o) => {
                if let Some(k) = key_of(row, c, o.iter().copied(), 0) {
                    out.push((k, c as u16));
                }
            }
            Seed::Bits(t) => {
                if let Some(k) = seed_len(row, bits, c, *t).and_then(|n| key_of(row, c, 0..n, 0)) {
                    out.push((k, c as u16));
                }
            }
            Seed::Ruler(t) => {
                let mut b = 0.0;
                for k in 0..RULER.len() {
                    match row.get(c + RULER[k]) {
                        Some(&a) if a < 20 => b += bits[c + RULER[k]][a as usize],
                        _ => break,
                    }
                    if b >= *t || k + 1 == RULER.len() {
                        let x = key_of(row, c, RULER[..=k].iter().copied(), 0).unwrap();
                        out.push((x, c as u16));
                        break;
                    }
                }
            }
            Seed::Target(_) => {
                let range = level.map_or(0..LEVELS.len(), |l| l..l + 1);
                for l in range {
                    if let Some(k) = seed_len(row, bits, c, LEVELS[l]).and_then(|n| key_of(row, c, 0..n, 1 + l as u64)) {
                        out.push((k, c as u16));
                    }
                }
            }
        }
    }
}

/// The highest level at which a 40-column read of the target, at THETA_ID, is expected to keep
/// `theta` seeds' worth of evidence; the lowest level if none.
pub fn level_of(row: &[u8], bits: &[[f32; AA]], theta: f32) -> usize {
    let n = row.iter().filter(|&&a| a < 20).count().max(1) as f32;
    (0..LEVELS.len())
        .rev()
        .find(|&l| {
            let q: f32 = (0..row.len())
                .filter_map(|c| seed_len(row, bits, c, LEVELS[l]))
                .map(|k| THETA_ID.powi(k as i32))
                .sum();
            q * (40.0f32.min(n) / n) >= theta
        })
        .unwrap_or(0)
}

pub struct Index {
    pub post: Vec<(u64, u32)>, // key, unit slot; sorted, no repeats
    pub keys: usize,
}
pub fn index(seed: &Seed, f: &Fam, step: usize) -> Index {
    let mut post = vec![];
    let mut ks = vec![];
    for (row, &slot) in f.rows.iter().zip(&f.slot) {
        ks.clear();
        let level = if let Seed::Target(th) = seed { Some(level_of(row, &f.bits, *th)) } else { None };
        keys(seed, row, &f.bits, level, &mut ks);
        post.extend(ks.iter().filter(|k| k.1 as usize % step == 0).map(|k| (k.0, slot)));
    }
    post.sort_unstable();
    post.dedup();
    let n = post.chunk_by(|a, b| a.0 == b.0).count();
    Index { post, keys: n }
}
/// Units with at least `h` distinct seeds shared with the query -> (slot, seeds).
pub fn lookup(ix: &Index, qkeys: &[(u64, u16)], count: &mut [u16], touched: &mut Vec<u32>) -> usize {
    touched.clear();
    let mut hits = 0;
    for (k, _) in qkeys {
        let a = ix.post.partition_point(|p| p.0 < *k);
        for p in ix.post[a..].iter().take_while(|p| p.0 == *k) {
            hits += 1;
            if count[p.1 as usize] == 0 {
                touched.push(p.1);
            }
            count[p.1 as usize] += 1;
        }
    }
    hits
}

// ---------------------------------------------------------------- options

pub fn shape(s: &str) -> Vec<usize> {
    s.bytes().enumerate().filter(|x| x.1 == b'1').map(|x| x.0).collect()
}
fn seed_options() -> Vec<(String, Seed, usize)> {
    let mut v = vec![];
    for (name, s) in [
        ("k3", "111"), ("k4", "1111"), ("k5", "11111"), ("k6", "111111"), ("k8", "11111111"),
        ("golomb4a", "10110001"), ("golomb4b", "1100101"), ("golomb5", "110010000101"),
        ("golomb6", "110010000010100001"), ("sp5", "1101011"), ("sp6", "11010010011"),
    ] {
        let steps: &[usize] = if matches!(name, "k5" | "golomb5" | "golomb4b") { &[1, 2, 4] } else { &[1] };
        for &st in steps {
            v.push((format!("{name}/s{st}"), Seed::Shape(shape(s)), st));
        }
    }
    for t in [6.0, 8.0, 9.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 24.0] {
        v.push((format!("bits{t}/s1"), Seed::Bits(t), 1));
        v.push((format!("ruler{t}/s1"), Seed::Ruler(t), 1));
    }
    for t in [8.0, 12.0, 16.0] {
        for st in [2, 4] {
            v.push((format!("bits{t}/s{st}"), Seed::Bits(t), st));
        }
    }
    for th in [2.0, 4.0, 6.0, 8.0, 12.0, 16.0] {
        v.push((format!("target{th}/s1"), Seed::Target(th), 1));
    }
    v
}

struct Score<'a> {
    truth: Vec<i16>, // by slot
    best: i16,
    tied: usize,
    src: u32, // slot of the source unit
    qid: usize,
    out: &'a mut String,
}
impl Score<'_> {
    fn emit(&mut self, option: &str, hits: usize, cands: impl Iterator<Item = u32>) {
        let (mut n, mut tied, mut near, mut src) = (0, 0, false, false);
        for s in cands {
            n += 1;
            let t = self.truth[s as usize];
            tied += (t == self.best) as usize;
            near |= t as f32 >= 0.95 * self.best as f32;
            src |= s == self.src;
        }
        writeln!(
            self.out, "{option}\t{}\t{hits}\t{n}\t{}\t{}\t{:.4}\t{}",
            self.qid, (tied > 0) as u8, near as u8, tied as f32 / self.tied as f32, src as u8
        )
        .unwrap();
    }
}

#[allow(clippy::too_many_arguments)]
pub fn run(dir: &Path, arm: &str, sub: u32, out: &Path, threads: usize, truth_arms: &[&str], only: &str) {
    let fams = families(dir);
    let q = Queries::load(dir, &fams);
    let units = {
        let r: Vec<Rows> = truth_arms.iter().map(|a| Rows::load(dir, a, &fams)).collect();
        family_units(&fams, &r.iter().collect::<Vec<_>>())
    };
    let (base, ga_only) = arm.strip_suffix("_ga").map_or((arm, false), |x| (x, true));
    let rows = Rows::load(dir, base, &fams);
    let ga = ga_only.then(|| bytes(&dir.join(format!("{base}.ga"))));
    let src = u32s(&dir.join("q.src"));
    let off = fam_offsets(&q, &units);
    let truth: Vec<i16> = u16s(&dir.join("truth.fam")).into_iter().map(|x| x as i16).collect();
    assert_eq!(truth.len(), *off.last().unwrap());
    let tsv = Mutex::new(String::from("option\tqid\thits\tcands\tbest\tnear\ttied\tsrc\n"));
    let idx = Mutex::new(String::from("option\tfamily\tunits\trows\tresidues\tpostings\tkeys\n"));
    let info = Mutex::new(String::from("qid\tfamily_units\tcolumns\tself_bits\tbest_score\ttied_units\n"));
    let secs = Mutex::new(HashMap::<String, (f64, usize)>::new());
    let wanted = |o: &str| only.is_empty() || only.split(',').any(|w| o.starts_with(w));
    parallel(fams.len(), threads, |fi| {
        let f = fam(&fams[fi], fi, &rows, ga.as_deref(), &units[fi], sub);
        if f.units.is_empty() {
            return;
        }
        let n = f.units.len();
        let residues: usize = f.rows.iter().map(|r| r.iter().filter(|&&a| a < 20).count()).sum();
        let mut o = String::new();
        let mut ixs = String::new();
        let mut inf = String::new();
        let mine: Vec<usize> = (0..q.len()).filter(|&i| q.fam[i] as usize == fi).collect();
        let score = |i: usize| -> Option<(Vec<i16>, i16, usize)> {
            let t: Vec<i16> = f.uslot.iter().map(|&s| truth[off[i] + s as usize]).collect();
            let best = *t.iter().max().unwrap();
            let tied = t.iter().filter(|&&x| x == best).count();
            Some((t, best, tied))
        };
        let src_slot = |i: usize| f.units.binary_search(&src[i]).map_or(u32::MAX, |x| x as u32);
        let scores: Vec<_> = mine.iter().map(|&i| score(i).unwrap()).collect();
        for (&i, s) in mine.iter().zip(&scores) {
            let row = q.row(i);
            let cols = row.iter().filter(|&&a| a < 20).count();
            let sb: f32 = (0..f.len).filter(|&c| row[c] < 20).map(|c| f.bits[c][row[c] as usize]).sum();
            writeln!(inf, "{i}\t{n}\t{cols}\t{sb:.1}\t{}\t{}", s.1, s.2).unwrap();
        }
        writeln!(ixs, "rows\t{}\t{n}\t{}\t{residues}\t{residues}\t0", fams[fi].acc, f.rows.len()).unwrap();

        // seeds
        let mut count = vec![0u16; n];
        let mut touched = vec![];
        let mut ks = vec![];
        for (name, seed, step) in seed_options() {
            if !wanted(&name) {
                continue;
            }
            let ix = index(&seed, &f, step);
            writeln!(ixs, "{name}\t{}\t{n}\t{}\t{residues}\t{}\t{}", fams[fi].acc, f.rows.len(), ix.post.len(), ix.keys).unwrap();
            let t0 = std::time::Instant::now();
            let mut spent = 0.0;
            for (&i, s) in mine.iter().zip(&scores) {
                let t1 = std::time::Instant::now();
                ks.clear();
                keys(&seed, q.row(i), &f.bits, None, &mut ks);
                let hits = lookup(&ix, &ks, &mut count, &mut touched);
                spent += t1.elapsed().as_secs_f64();
                let mut sc = Score { truth: s.0.clone(), best: s.1, tied: s.2, src: src_slot(i), qid: i, out: &mut o };
                sc.emit(&format!("{name}/h1"), hits, touched.iter().copied());
                sc.emit(&format!("{name}/h2"), hits, touched.iter().copied().filter(|&u| count[u as usize] >= 2));
                touched.iter().for_each(|&u| count[u as usize] = 0);
            }
            let _ = t0;
            let mut m = secs.lock().unwrap();
            let e = m.entry(name).or_default();
            e.0 += spent;
            e.1 += mine.len();
        }
        if !wanted("scan") {
            tsv.lock().unwrap().push_str(&o);
            idx.lock().unwrap().push_str(&ixs);
            info.lock().unwrap().push_str(&inf);
            return;
        }

        // residues: every row scored on the read's columns. A column's residue counts its bits.
        let mut by_bits: Vec<f32> = vec![];
        let col_count: Vec<[u32; AA]> = (0..f.len)
            .map(|c| {
                let mut k = [0u32; AA];
                f.rows.iter().for_each(|r| k[r[c] as usize] += 1);
                k
            })
            .collect();
        for min_bits in [1.0f32, 2.0, 3.0, 4.0, 5.0] {
            let p: usize = f.rows.iter().map(|r| (0..f.len).filter(|&c| r[c] < 20 && f.bits[c][r[c] as usize] >= min_bits).count()).sum();
            writeln!(ixs, "rare{min_bits}\t{}\t{n}\t{}\t{residues}\t{p}\t0", fams[fi].acc, f.rows.len()).unwrap();
        }
        let bl = blosum();
        let theta = 0.8f32;
        for (&i, s) in mine.iter().zip(&scores) {
            let row = q.row(i);
            let cols: Vec<usize> = (0..f.len).filter(|&c| row[c] < 20).collect();
            let bit = |c: usize| f.bits[c][row[c] as usize];
            let mut by_col = cols.clone();
            by_col.sort_by(|&a, &b| bit(b).total_cmp(&bit(a)));
            let mut by_ent = cols.clone();
            by_ent.sort_by(|&a, &b| f.entropy[b].total_cmp(&f.entropy[a]));
            let self_q: f32 = cols.iter().map(|&c| bit(c)).sum();
            let mut sc = Score { truth: s.0.clone(), best: s.1, tied: s.2, src: src_slot(i), qid: i, out: &mut o };
            // one score per unit, the best of its rows
            let per_unit = |g: &dyn Fn(&[u8]) -> f32| -> Vec<f32> {
                let mut u = vec![f32::NEG_INFINITY; n];
                for (r, &sl) in f.rows.iter().zip(&f.slot) {
                    let x = g(r);
                    if x > u[sl as usize] {
                        u[sl as usize] = x;
                    }
                }
                u
            };
            let timed = |name: &str, g: &dyn Fn(&[u8]) -> f32| -> Vec<f32> {
                let t = std::time::Instant::now();
                let u = per_unit(g);
                let mut m = secs.lock().unwrap();
                let e = m.entry(name.to_string()).or_default();
                e.0 += t.elapsed().as_secs_f64();
                e.1 += f.rows.len();
                u
            };
            let all = timed("scan_bits", &|r| cols.iter().filter(|&&c| r[c] == row[c]).map(|&c| bit(c)).sum());
            let top = |k: usize| by_col[..k.min(by_col.len())].to_vec();
            let top16 = top(16);
            let top8 = top(8);
            let b16 = timed("scan_bits16", &|r| top16.iter().filter(|&&c| r[c] == row[c]).map(|&c| bit(c)).sum());
            let b8 = timed("scan_bits8", &|r| top8.iter().filter(|&&c| r[c] == row[c]).map(|&c| bit(c)).sum());
            let e16: Vec<usize> = by_ent[..16.min(by_ent.len())].to_vec();
            let ent = timed("scan_ent16", &|r| e16.iter().filter(|&&c| r[c] == row[c]).count() as f32);
            let ident = timed("scan_ident", &|r| cols.iter().filter(|&&c| r[c] == row[c]).count() as f32);
            let by_blosum = timed("scan_blosum", &|r| cols.iter().filter(|&&c| r[c] < 21).map(|&c| bl[row[c] as usize][r[c] as usize] as f32).sum());
            let miss = (1.0 - theta).log2();
            let llr = timed("scan_llr", &|r| {
                cols.iter()
                    .map(|&c| {
                        if r[c] == row[c] {
                            let p = (-bit(c)).exp2();
                            ((theta + (1.0 - theta) * p) / p).log2()
                        } else {
                            miss
                        }
                    })
                    .sum()
            });
            // a target's own bits over the read's columns: what a perfect match would score
            let self_t = per_unit(&|r| cols.iter().filter(|&&c| r[c] < 20).map(|&c| f.bits[c][r[c] as usize]).sum());
            let lists: [(&str, &Vec<f32>, usize); 7] = [
                ("bits", &all, n * cols.len()), ("bits16", &b16, n * 16), ("bits8", &b8, n * 8), ("ent16", &ent, n * 16),
                ("ident", &ident, n * cols.len()), ("blosum", &by_blosum, n * cols.len()), ("llr", &llr, n * cols.len()),
            ];
            for (name, u, hits) in lists {
                by_bits.clear();
                by_bits.extend_from_slice(u);
                by_bits.sort_by(|a, b| b.total_cmp(a));
                let mut rank = |label: String, k: usize, sc: &mut Score| {
                    let bar = by_bits[k.clamp(1, n) - 1].max(0.0);
                    sc.emit(&label, hits, (0..n as u32).filter(|&x| u[x as usize] >= bar && u[x as usize] > 0.0));
                };
                for k in [1usize, 20, 50, 100, 200, 500, 1000] {
                    rank(format!("scan_{name}/top{k}"), k, &mut sc);
                }
                for p in [0.005, 0.01, 0.02] {
                    rank(format!("scan_{name}/share{p}"), (p * n as f64).ceil() as usize, &mut sc);
                }
            }
            for (name, u) in [("bits", &all), ("bits16", &b16), ("llr", &llr)] {
                for t in [10.0f32, 15.0, 20.0, 25.0, 30.0, 40.0] {
                    sc.emit(&format!("scan_{name}/T{t}"), n, (0..n as u32).filter(|&x| u[x as usize] >= t));
                }
                for d in [0.0f32, 4.0, 8.0, 12.0] {
                    let t = (n as f32).log2() + d;
                    sc.emit(&format!("scan_{name}/N+{d}"), n, (0..n as u32).filter(|&x| u[x as usize] >= t));
                }
            }
            for p in [0.2f32, 0.3, 0.4, 0.5, 0.6] {
                sc.emit(&format!("scan_bits/selfq{p}"), n, (0..n as u32).filter(|&x| all[x as usize] >= p * self_q));
                sc.emit(&format!("scan_bits/selft{p}"), n, (0..n as u32).filter(|&x| {
                    let x = x as usize;
                    all[x] > 0.0 && all[x] >= p * self_t[x]
                }));
            }
            // within a margin of the read's best, in bits
            for (name, u) in [("bits", &all), ("llr", &llr)] {
                let top = u.iter().copied().fold(0.0, f32::max);
                for d in [5.0f32, 10.0, 20.0, 30.0] {
                    sc.emit(&format!("scan_{name}/best-{d}"), n, (0..n as u32).filter(|&x| u[x as usize] > 0.0 && u[x as usize] >= top - d));
                }
            }
            // posting lists of rare residues only
            for min_bits in [2.0f32, 3.0, 4.0] {
                let rare: Vec<usize> = cols.iter().copied().filter(|&c| bit(c) >= min_bits).collect();
                let hits: usize = rare.iter().map(|&c| col_count[c][row[c] as usize] as usize).sum();
                let u = timed(&format!("rare{min_bits}"), &|r| rare.iter().filter(|&&c| r[c] == row[c]).map(|&c| bit(c)).sum());
                for t in [8.0f32, 12.0, 16.0, 20.0, 25.0] {
                    sc.emit(&format!("rare{min_bits}/T{t}"), hits, (0..n as u32).filter(|&x| u[x as usize] >= t));
                }
                for k in [100usize, 500] {
                    by_bits.clear();
                    by_bits.extend_from_slice(&u);
                    by_bits.sort_by(|a, b| b.total_cmp(a));
                    let bar = by_bits[k.clamp(1, n) - 1];
                    sc.emit(&format!("rare{min_bits}/top{k}"), hits, (0..n as u32).filter(|&x| u[x as usize] > 0.0 && u[x as usize] >= bar));
                }
            }
        }
        eprintln!("PF{:05} {n} units {} queries", fams[fi].acc, mine.len());
        tsv.lock().unwrap().push_str(&o);
        idx.lock().unwrap().push_str(&ixs);
        info.lock().unwrap().push_str(&inf);
    });
    std::fs::create_dir_all(out).unwrap();
    std::fs::write(out.join("cand.tsv"), &*tsv.lock().unwrap()).unwrap();
    std::fs::write(out.join("index.tsv"), &*idx.lock().unwrap()).unwrap();
    std::fs::write(out.join("queries.tsv"), &*info.lock().unwrap()).unwrap();
    let mut s = String::from("option\tseconds\tunits_of_work\n");
    let m = secs.lock().unwrap();
    let mut names: Vec<_> = m.keys().collect();
    names.sort();
    for k in names {
        writeln!(s, "{k}\t{:.6}\t{}", m[k].0, m[k].1).unwrap();
    }
    std::fs::write(out.join("seconds.tsv"), s).unwrap();
}

#[cfg(test)]
mod tests {
    use super::*;
    /// A row is found by its own seeds, whatever the seed; a deleted column breaks the seeds over it only.
    #[test]
    fn own_seeds() {
        let bits = vec![[3.0f32; AA]; 40];
        let row: Vec<u8> = (0..40u8).map(|i| (i * 7) % 20).collect();
        let mut gap = row.clone();
        gap[20] = NONE;
        for seed in [Seed::Shape(shape("1100101")), Seed::Bits(10.0), Seed::Ruler(10.0), Seed::Target(4.0)] {
            let f = Fam {
                len: 40, rows: vec![&row[..], &gap[..]], pos: vec![], slot: vec![0, 1], units: vec![5, 9], uslot: vec![0, 1],
                bits: bits.clone(), entropy: vec![0.0; 40],
            };
            let ix = index(&seed, &f, 1);
            let (mut ks, mut count, mut touched) = (vec![], vec![0u16; 2], vec![]);
            keys(&seed, &row, &bits, None, &mut ks);
            lookup(&ix, &ks, &mut count, &mut touched);
            assert!(touched.contains(&0), "a row must find itself");
            assert!(count[1] < count[0] && count[1] > 0, "the gap breaks some seeds, not all");
        }
    }
}
