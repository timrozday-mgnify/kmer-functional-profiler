//! Pseudo-alignment: k-mer scores standing in for an ungapped alignment of a read pair's
//! translated frames against a leaf's representative (docs/tree-index-plan.md, Progress log:
//! pseudo-alignment trial). The representative is reduced to a sketch, as an index would store
//! it (sampled k-mers, positions coarsened or dropped); the read keeps all its k-mers.

use crate::threads;

/// How a (read, representative) pair is scored.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Method {
    /// Most matched k-mers on one diagonal band (positions stored, coarsened by `quantum`).
    Diag,
    /// Longest chain of matched k-mers in the same order in both (ranks stored, no positions).
    Order,
    /// Distinct shared k-mer pairs at most `gap` apart, the gap binned (no positions stored).
    Pairs,
    /// [`Method::Diag`]'s `top` bands scanned exactly (ungapped, `matrix`): needs the sequence.
    Anchored,
    /// As [`Method::Anchored`], but banded Smith-Waterman with affine gaps (open 11, extend 1)
    /// within `gap` diagonals of each band's centre: indels.
    Gapped,
}

/// A [`Method`] and its sketch.
#[derive(Clone, Debug)]
pub struct Params {
    pub method: Method,
    pub k: usize,
    /// Residue byte to class; 255 breaks k-mers (stops, unknowns).
    pub classes: [u8; 256],
    /// Representative k-mers kept: the least hash per window of this many (1: all) ...
    pub minimizer: usize,
    /// ... and of those, hashes divisible by this (1: all).
    pub modulo: u32,
    /// Representative positions stored as position / quantum.
    pub quantum: u32,
    /// Diagonal band width ([`Method::Diag`], [`Method::Anchored`]).
    pub band: u32,
    /// [`Method::Pairs`]: k-mers paired up to this many positions apart, the gap stored as
    /// gap / `gap_bin`; [`Method::Gapped`]: the band's half-width in diagonals.
    pub gap: usize,
    pub gap_bin: usize,
    /// The pair filter (`pair_route`): k-mer pairs at least this many positions apart.
    pub min_gap: usize,
    /// The pair filter: k-mers per feature, 2 (pairs) or 3 (triples; consecutive ones
    /// `min_gap` or more apart, the first to the last at most `gap`).
    pub combo: usize,
    /// [`Method::Anchored`]: bands scanned.
    pub top: usize,
}

fn mix(mut x: u64) -> u64 {
    x = (x ^ (x >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    x = (x ^ (x >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    x ^ (x >> 31)
}

/// (hash, position) of every k-mer of `seq` whose residues all have a class.
fn kmers(seq: &[u8], p: &Params) -> Vec<(u32, u32)> {
    let (mut code, mut run) = (0u64, 0usize);
    let mut out = Vec::with_capacity(seq.len());
    for (i, &c) in seq.iter().enumerate() {
        let x = p.classes[usize::from(c)];
        if x == 255 {
            run = 0;
            continue;
        }
        code = (code << 5 | u64::from(x)) & ((1 << (5 * p.k)) - 1);
        run += 1;
        if run >= p.k {
            out.push((mix(code ^ (p.k as u64) << 60) as u32, (i + 1 - p.k) as u32));
        }
    }
    out
}

/// The representative's stored k-mers, in position order (positions coarsened).
fn sampled(seq: &[u8], p: &Params) -> Vec<(u32, u32)> {
    let all = kmers(seq, p);
    let mut out: Vec<(u32, u32)> = Vec::new();
    if p.minimizer <= 1 || all.is_empty() {
        out = all;
    } else {
        for w in all.windows(p.minimizer.min(all.len())) {
            let m = *w.iter().min_by_key(|x| x.0).unwrap();
            if out.last() != Some(&m) {
                out.push(m);
            }
        }
    }
    out.retain(|x| x.0 % p.modulo == 0);
    out.iter_mut().for_each(|x| x.1 /= p.quantum);
    out
}

fn pair_hash(a: u32, b: u32, d: usize) -> u32 {
    mix(u64::from(a) << 32 ^ u64::from(b) ^ (d as u64) << 58) as u32
}

/// [`Method::Pairs`]' features of `seq`, sorted and distinct (sampled by `modulo` only).
fn pairs(seq: &[u8], p: &Params) -> Vec<u32> {
    let all = kmers(seq, p);
    let mut out = Vec::new();
    for (i, &(a, pa)) in all.iter().enumerate() {
        for &(b, pb) in all[i + 1..]
            .iter()
            .take_while(|x| (x.1 - pa) as usize <= p.gap)
        {
            let h = pair_hash(a, b, (pb - pa) as usize / p.gap_bin);
            if h % p.modulo == 0 {
                out.push(h);
            }
        }
    }
    out.sort_unstable();
    out.dedup();
    out
}

/// A representative as stored: by hash for lookups, and its features' count.
struct Sketch {
    by_hash: Vec<(u32, u32)>, // (hash, position / quantum or rank)
    pairs: Vec<u32>,
}

impl Sketch {
    fn new(seq: &[u8], p: &Params) -> Self {
        if p.method == Method::Pairs {
            return Self {
                by_hash: Vec::new(),
                pairs: pairs(seq, p),
            };
        }
        let mut by_hash = sampled(seq, p);
        if p.method == Method::Order {
            by_hash
                .iter_mut()
                .enumerate()
                .for_each(|(r, x)| x.1 = r as u32);
        }
        by_hash.sort_unstable();
        Self {
            by_hash,
            pairs: Vec::new(),
        }
    }

    fn len(&self) -> usize {
        self.by_hash.len() + self.pairs.len()
    }

    fn hits(&self, h: u32) -> &[(u32, u32)] {
        let lo = self.by_hash.partition_point(|x| x.0 < h);
        let hi = lo + self.by_hash[lo..].partition_point(|x| x.0 == h);
        &self.by_hash[lo..hi]
    }
}

/// Read positions with a match per diagonal band (diagonal = rep position - read position +
/// read length): a position counts once per band, so repeats matching repeats do not pile up.
fn bands(read: &[u8], rep_len: usize, s: &Sketch, p: &Params) -> Vec<u32> {
    let n = (rep_len + read.len()) / p.band as usize + 2;
    let (mut c, mut last) = (vec![0u32; n], vec![u32::MAX; n]);
    for (h, i) in kmers(read, p) {
        for &(_, j) in s.hits(h) {
            let b = ((j * p.quantum) as usize + read.len() - i as usize) / p.band as usize;
            if last[b] != i {
                (c[b], last[b]) = (c[b] + 1, i);
            }
        }
    }
    c
}

fn best_band(c: &[u32]) -> u32 {
    c.windows(2).map(|w| w[0] + w[1]).max().unwrap_or(0)
}

/// Longest chain of matches increasing in both read position and representative rank.
fn chain(read: &[u8], s: &Sketch, p: &Params) -> u32 {
    let mut tails: Vec<u32> = Vec::new(); // least rank ending a chain of each length
    for (h, _) in kmers(read, p) {
        // one read position: its ranks in decreasing order, so none chains to another
        for &(_, r) in s.hits(h).iter().rev() {
            let at = tails.partition_point(|&t| t < r);
            if at == tails.len() {
                tails.push(r);
            } else {
                tails[at] = r;
            }
        }
    }
    tails.len() as u32
}

fn shared_pairs(read: &[u8], s: &Sketch, p: &Params) -> u32 {
    pairs(read, p)
        .iter()
        .filter(|h| s.pairs.binary_search(h).is_ok())
        .count() as u32
}

/// Exact ungapped best segment over the diagonals of the `top` best bands (pairs of bands).
fn anchored(read: &[u8], rep: &[u8], s: &Sketch, p: &Params, matrix: &[i8]) -> i32 {
    let c = bands(read, rep.len(), s, p);
    let mut order: Vec<usize> = (0..c.len().saturating_sub(1))
        .filter(|&b| c[b] + c[b + 1] > 0)
        .collect();
    order.sort_unstable_by_key(|&b| std::cmp::Reverse(c[b] + c[b + 1]));
    let mut best = 0;
    for &b in order.iter().take(p.top) {
        let lo = b * p.band as usize;
        for d in lo..lo + 2 * p.band as usize {
            // rep j = read i + d - read.len()
            let mut run = 0i32;
            for i in 0..read.len() {
                let Some(j) = (i + d).checked_sub(read.len()) else {
                    continue;
                };
                if j >= rep.len() {
                    break;
                }
                run = (run + i32::from(matrix[usize::from(read[i]) << 8 | usize::from(rep[j])]))
                    .max(0);
                best = best.max(run);
            }
        }
    }
    best
}

const GAP_OPEN: i32 = 11;
const GAP_EXTEND: i32 = 1;

/// The best band-pairs of [`bands`], most matches first.
fn top_bands(c: &[u32], top: usize) -> Vec<usize> {
    let mut order: Vec<usize> = (0..c.len().saturating_sub(1))
        .filter(|&b| c[b] + c[b + 1] > 0)
        .collect();
    order.sort_unstable_by_key(|&b| std::cmp::Reverse(c[b] + c[b + 1]));
    order.truncate(top);
    order
}

/// Local alignment score (affine gaps) of `a` against `b` within `w` diagonals of diagonal
/// `d0` (b's j = a's i + d0), Gotoh's recurrences on the band only.
fn banded(a: &[u8], b: &[u8], d0: i64, w: i64, matrix: &[i8]) -> i32 {
    let width = (2 * w + 1) as usize;
    let neg = i32::MIN / 4;
    // row i: column k is j = i + d0 - w + k
    let (mut h, mut f) = (vec![0i32; width + 2], vec![neg; width + 2]);
    let mut best = 0;
    for (i, &x) in a.iter().enumerate() {
        let (mut nh, mut nf) = (vec![0i32; width + 2], vec![neg; width + 2]);
        let mut e = neg; // gap along b (left neighbour in this row)
        for k in 0..width {
            let j = i as i64 + d0 - w + k as i64;
            if j < 0 || j >= b.len() as i64 {
                e = neg;
                continue;
            }
            // previous row: (i-1, j-1) is column k, (i-1, j) is column k+1
            let diag = h[k] + i32::from(matrix[usize::from(x) << 8 | usize::from(b[j as usize])]);
            nf[k] = (h[k + 1] - GAP_OPEN).max(f[k + 1]) - GAP_EXTEND;
            let cell = diag.max(e).max(nf[k]).max(0);
            nh[k] = cell;
            best = best.max(cell);
            e = (cell - GAP_OPEN).max(e) - GAP_EXTEND;
        }
        (h, f) = (nh, nf);
    }
    best
}

/// [`Method::Gapped`]: banded affine-gap scans around the `top` best bands' centres.
fn gapped(read: &[u8], rep: &[u8], s: &Sketch, p: &Params, matrix: &[i8]) -> i32 {
    let c = bands(read, rep.len(), s, p);
    top_bands(&c, p.top)
        .into_iter()
        .map(|b| {
            let centre = (b * p.band as usize + p.band as usize) as i64 - read.len() as i64;
            banded(read, rep, centre, p.gap as i64 + p.band as i64, matrix)
        })
        .max()
        .unwrap_or(0)
}

/// Each pair's score: read `reads[pair_read[i]]` against representative `reps[pair_rep[i]]`.
/// Pairs sorted by representative build each sketch once per run of it.
pub fn pseudo_scores<S: AsRef<[u8]> + Sync>(
    reads: &[S],
    reps: &[S],
    pair_read: &[u32],
    pair_rep: &[u32],
    p: &Params,
    matrix: &[i8],
) -> Vec<i32> {
    let n = pair_read.len().min(pair_rep.len());
    let mut out = vec![0i32; n];
    let per_thread = n.div_ceil(threads()).max(1);
    std::thread::scope(|scope| {
        for (t, chunk) in out.chunks_mut(per_thread).enumerate() {
            scope.spawn(move || {
                let mut cur: Option<(u32, Sketch)> = None;
                for (k, o) in chunk.iter_mut().enumerate() {
                    let i = t * per_thread + k;
                    let (r, q) = (pair_read[i], pair_rep[i]);
                    let rep = reps[q as usize].as_ref();
                    if cur.as_ref().is_none_or(|c| c.0 != q) {
                        cur = Some((q, Sketch::new(rep, p)));
                    }
                    let s = &cur.as_ref().unwrap().1;
                    let read = reads[r as usize].as_ref();
                    *o = match p.method {
                        Method::Diag => best_band(&bands(read, rep.len(), s, p)) as i32,
                        Method::Order => chain(read, s, p) as i32,
                        Method::Pairs => shared_pairs(read, s, p) as i32,
                        Method::Anchored => anchored(read, rep, s, p, matrix),
                        Method::Gapped => gapped(read, rep, s, p, matrix),
                    };
                }
            });
        }
    });
    out
}

/// `seq` with tandem repeats masked as `x`: at each position the shortest period (1-3) with at
/// least four copies there, all its copies (leftmost first, as the regex
/// `(.{1,3}?)\1{3,}` substitutes; tier 1's `verify.mask`).
pub fn mask_repeats(seq: &[u8]) -> Vec<u8> {
    let mut out = seq.to_vec();
    let mut i = 0;
    while i < seq.len() {
        let mut run = 0;
        for u in 1..=3 {
            if i + 4 * u > seq.len() {
                break;
            }
            let mut copies = 1;
            while i + (copies + 1) * u <= seq.len()
                && seq[i + copies * u..i + (copies + 1) * u] == seq[i..i + u]
            {
                copies += 1;
            }
            if copies >= 4 {
                run = copies * u;
                break;
            }
        }
        if run > 0 {
            out[i..i + run].fill(b'x');
            i += run;
        } else {
            i += 1;
        }
    }
    out
}

/// A read pair's twelve frames (each mate's six, forward then reverse), masked and joined by
/// stops, for each pair (`r1[i]`, `r2[i]`), in parallel: what tier 2 aligns and routes.
pub fn pair_frames<S: AsRef<[u8]> + Sync>(
    r1: &[S],
    r2: &[S],
    code: crate::GeneticCode,
) -> Vec<Vec<u8>> {
    let n = r1.len().min(r2.len());
    let mut out = vec![Vec::new(); n];
    let per_thread = n.div_ceil(threads()).max(1);
    std::thread::scope(|scope| {
        for (t, chunk) in out.chunks_mut(per_thread).enumerate() {
            scope.spawn(move || {
                for (k, o) in chunk.iter_mut().enumerate() {
                    let i = t * per_thread + k;
                    let mut joined = Vec::new();
                    for mate in [r1[i].as_ref(), r2[i].as_ref()] {
                        for f in crate::six_frames(mate, code) {
                            if !joined.is_empty() {
                                joined.push(b'*');
                            }
                            joined.extend_from_slice(&f);
                        }
                    }
                    *o = mask_repeats(&joined);
                }
            });
        }
    });
    out
}

/// Dense routing for tier 2 (Progress log: tier 2, dense routing): every read's candidate
/// leaves by seed votes. `leaves`' k-mers (`p`'s k and classes; those whose hash is
/// divisible by `p.modulo`, all at 1: a sparse index for the pair filter after tier 1) are
/// indexed, those held more than `max_occ` times dropped; each read votes, each of its positions once
/// per (leaf, diagonal band of width `p.band`), and a leaf's votes are its best two adjacent
/// bands'. Returns (read, leaf, votes) for each read's `top` leaves with at least `min_votes`
/// (ties at the cut kept), reads in order.
pub fn route_votes<S: AsRef<[u8]> + Sync>(
    reads: &[S],
    leaves: &[S],
    p: &Params,
    min_votes: u32,
    top: usize,
    max_occ: usize,
) -> (Vec<u32>, Vec<u32>, Vec<u32>) {
    let mut index: Vec<(u32, u32, u32)> = Vec::new(); // (hash, leaf, position)
    for (l, seq) in leaves.iter().enumerate() {
        let modulo = p.modulo.max(1);
        index.extend(
            kmers(seq.as_ref(), p)
                .into_iter()
                .filter(|&(h, _)| h % modulo == 0)
                .map(|(h, i)| (h, l as u32, i)),
        );
    }
    index.sort_unstable();
    // drop k-mers held more than max_occ times (repeats, low complexity left after masking)
    let mut kept = Vec::with_capacity(index.len());
    let mut i = 0;
    while i < index.len() {
        let j = i + index[i..].partition_point(|x| x.0 == index[i].0);
        if j - i <= max_occ {
            kept.extend_from_slice(&index[i..j]);
        }
        i = j;
    }
    let index = kept;
    let band = i64::from(p.band.max(1));
    let n = reads.len();
    let per_thread = n.div_ceil(threads()).max(1);
    let mut out: Vec<Vec<(u32, u32, u32)>> = vec![Vec::new(); n.div_ceil(per_thread).max(1)];
    std::thread::scope(|scope| {
        for (t, slot) in out.iter_mut().enumerate() {
            let index = &index;
            scope.spawn(move || {
                let mut votes: std::collections::HashMap<(u32, i64), (u32, u32)> =
                    std::collections::HashMap::new();
                let first = t * per_thread;
                let mine = &reads[first..((t + 1) * per_thread).min(n)];
                for (r, read) in mine
                    .iter()
                    .enumerate()
                    .map(|(i, x)| (first + i, x.as_ref()))
                {
                    votes.clear();
                    for (h, i) in kmers(read, p) {
                        let lo = index.partition_point(|x| x.0 < h);
                        for &(_, l, j) in index[lo..].iter().take_while(|x| x.0 == h) {
                            let b = (i64::from(j) - i64::from(i)).div_euclid(band);
                            let e = votes.entry((l, b)).or_insert((0, u32::MAX));
                            if e.1 != i {
                                *e = (e.0 + 1, i); // a read position once per band
                            }
                        }
                    }
                    let mut best: std::collections::HashMap<u32, u32> =
                        std::collections::HashMap::new();
                    for (&(l, b), &(c, _)) in &votes {
                        let next = votes.get(&(l, b + 1)).map_or(0, |x| x.0);
                        let v = best.entry(l).or_insert(0);
                        *v = (*v).max(c + next);
                    }
                    let mut got: Vec<(u32, u32)> =
                        best.into_iter().filter(|&(_, v)| v >= min_votes).collect();
                    got.sort_unstable_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(&b.0)));
                    if got.len() > top {
                        let cut = got[top - 1].1;
                        got.retain(|x| x.1 >= cut);
                    }
                    slot.extend(got.into_iter().map(|(l, v)| (r as u32, l, v)));
                }
            });
        }
    });
    let (mut a, mut b, mut c) = (Vec::new(), Vec::new(), Vec::new());
    for (r, l, v) in out.into_iter().flatten() {
        a.push(r);
        b.push(l);
        c.push(v);
    }
    (a, b, c)
}

/// Pair features of `seq` (k-mer pairs `p.min_gap` to `p.gap` positions apart, in order, the
/// gap not hashed: an indel between them keeps the pair; with `p.combo` 3, triples; with
/// `keep`, only k-mers starting where it is non-zero), in position order; with `sample`, window minimizers over them
/// (`p.minimizer` > 1) and then those divisible by `p.modulo`. Distinct, sorted.
fn pair_features(seq: &[u8], p: &Params, sample: bool, keep: Option<&[u8]>) -> Vec<u32> {
    let all: Vec<(u32, u32)> = kmers(seq, p)
        .into_iter()
        .filter(|x| keep.is_none_or(|m| m.get(x.1 as usize).is_some_and(|&b| b != 0)))
        .collect();
    let mut feats = Vec::new();
    for (i, &(a, pa)) in all.iter().enumerate() {
        for &(b, pb) in all[i + 1..]
            .iter()
            .take_while(|x| (x.1 - pa) as usize <= p.gap)
            .filter(|x| (x.1 - pa) as usize >= p.min_gap)
        {
            if p.combo < 3 {
                feats.push(pair_hash(a, b, 0));
                continue;
            }
            let ab = pair_hash(a, b, 0);
            for &(c, _) in all[i + 1..]
                .iter()
                .take_while(|x| (x.1 - pa) as usize <= p.gap)
                .filter(|x| x.1 >= pb && (x.1 - pb) as usize >= p.min_gap)
            {
                feats.push(pair_hash(ab, c, 1));
            }
        }
    }
    if sample && p.minimizer > 1 && !feats.is_empty() {
        let mut mins = Vec::new();
        for w in feats.windows(p.minimizer.min(feats.len())) {
            let m = *w.iter().min().unwrap();
            if mins.last() != Some(&m) {
                mins.push(m);
            }
        }
        feats = mins;
    }
    let modulo = p.modulo.max(1);
    feats.retain(|h| h % modulo == 0);
    feats.sort_unstable();
    feats.dedup();
    feats
}

/// The pair filter after tier 1: every read's candidate representatives by shared pair
/// features (two k-mers in order within a window, hashed as one, so one shared feature is
/// already two seeds close together). The representatives' features are sampled
/// (`pair_features`: minimizers and/or modulo, a sparse index) and indexed, those held more
/// than `max_occ` times dropped; a read's features are all looked up (those divisible by
/// `p.modulo`, which any indexed one is). Returns (read, rep, score) for each read's `top`
/// representatives scoring at least `min_hits` (ties at the cut kept). The score is the
/// distinct shared features, or with `min_bits` each feature's information in quarter-bits,
/// log2(representatives / representatives holding it) (promiscuous features carry little),
/// features below `min_bits` dropped.
pub fn pair_route<S: AsRef<[u8]> + Sync>(
    reads: &[S],
    reps: &[S],
    p: &Params,
    min_hits: u32,
    top: usize,
    max_occ: usize,
    min_bits: Option<f32>,
) -> (Vec<u32>, Vec<u32>, Vec<u32>) {
    let mut index: Vec<(u32, u32)> = Vec::new(); // (feature, rep)
    for (r, seq) in reps.iter().enumerate() {
        index.extend(
            pair_features(seq.as_ref(), p, true, None)
                .into_iter()
                .map(|h| (h, r as u32)),
        );
    }
    pair_route_index(reads, index, p, min_hits, top, max_occ, min_bits)
}

/// Each sequence's sampled pair features, flattened: (sequence, feature).
pub fn pair_features_many<S: AsRef<[u8]> + Sync>(seqs: &[S], p: &Params) -> (Vec<u32>, Vec<u32>) {
    let (mut a, mut b) = (Vec::new(), Vec::new());
    for (i, s) in seqs.iter().enumerate() {
        for h in pair_features(s.as_ref(), p, true, None) {
            a.push(i as u32);
            b.push(h);
        }
    }
    (a, b)
}

/// Each sequence's pair features from the k-mers starting where its mask is non-zero only
/// (e.g. the 3-mers inside a group's stored keys), unsampled: (sequence, feature).
pub fn pair_features_masked_many<S: AsRef<[u8]> + Sync, M: AsRef<[u8]> + Sync>(
    seqs: &[S],
    masks: &[M],
    p: &Params,
) -> (Vec<u32>, Vec<u32>) {
    let (mut a, mut b) = (Vec::new(), Vec::new());
    for (i, (s, m)) in seqs.iter().zip(masks).enumerate() {
        for h in pair_features(s.as_ref(), p, false, Some(m.as_ref())) {
            a.push(i as u32);
            b.push(h);
        }
    }
    (a, b)
}

/// [`pair_route`] against a given index of (feature, rep) entries (e.g. a group's features
/// held by a share of its leaves, a soft intersection built by the caller).
pub fn pair_route_index<S: AsRef<[u8]> + Sync>(
    reads: &[S],
    mut index: Vec<(u32, u32)>,
    p: &Params,
    min_hits: u32,
    top: usize,
    max_occ: usize,
    min_bits: Option<f32>,
) -> (Vec<u32>, Vec<u32>, Vec<u32>) {
    index.sort_unstable();
    index.dedup();
    let mut reps: Vec<u32> = index.iter().map(|x| x.1).collect();
    reps.sort_unstable();
    reps.dedup();
    let n_reps = reps.len().max(1) as f32;
    let mut kept: Vec<(u32, u32, u32)> = Vec::with_capacity(index.len()); // (feature, rep, weight)
    let mut i = 0;
    while i < index.len() {
        let j = i + index[i..].partition_point(|x| x.0 == index[i].0);
        let w = match min_bits {
            None => Some(1),
            Some(b) => {
                let bits = (n_reps / (j - i) as f32).log2();
                (bits >= b).then(|| ((bits * 4.0).round() as u32).max(1))
            }
        };
        if let (true, Some(w)) = (j - i <= max_occ, w) {
            kept.extend(index[i..j].iter().map(|&(h, q)| (h, q, w)));
        }
        i = j;
    }
    let index = kept;
    let n = reads.len();
    let per_thread = n.div_ceil(threads()).max(1);
    let mut out: Vec<Vec<(u32, u32, u32)>> = vec![Vec::new(); n.div_ceil(per_thread).max(1)];
    std::thread::scope(|scope| {
        for (t, slot) in out.iter_mut().enumerate() {
            let index = &index;
            scope.spawn(move || {
                let mut hits: std::collections::HashMap<u32, u32> =
                    std::collections::HashMap::new();
                let first = t * per_thread;
                for (r, read) in reads[first..((t + 1) * per_thread).min(n)]
                    .iter()
                    .enumerate()
                    .map(|(i, x)| (first + i, x.as_ref()))
                {
                    hits.clear();
                    for h in pair_features(read, p, false, None) {
                        let lo = index.partition_point(|x| x.0 < h);
                        for &(_, q, w) in index[lo..].iter().take_while(|x| x.0 == h) {
                            *hits.entry(q).or_insert(0) += w;
                        }
                    }
                    let mut got: Vec<(u32, u32)> = hits
                        .iter()
                        .filter(|x| *x.1 >= min_hits)
                        .map(|(&q, &c)| (q, c))
                        .collect();
                    got.sort_unstable_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(&b.0)));
                    if got.len() > top {
                        let cut = got[top - 1].1;
                        got.retain(|x| x.1 >= cut);
                    }
                    slot.extend(got.into_iter().map(|(q, c)| (r as u32, q, c)));
                }
            });
        }
    });
    let (mut a, mut b, mut c) = (Vec::new(), Vec::new(), Vec::new());
    for (r, q, h) in out.into_iter().flatten() {
        a.push(r);
        b.push(q);
        c.push(h);
    }
    (a, b, c)
}

/// Each representative's sampled pair features (`pair_route`'s index entries).
pub fn pair_feature_counts<S: AsRef<[u8]> + Sync>(reps: &[S], p: &Params) -> Vec<u32> {
    reps.iter()
        .map(|s| pair_features(s.as_ref(), p, true, None).len() as u32)
        .collect()
}

/// Features each representative's sketch stores.
pub fn sketch_sizes<S: AsRef<[u8]> + Sync>(reps: &[S], p: &Params) -> Vec<u32> {
    reps.iter()
        .map(|r| Sketch::new(r.as_ref(), p).len() as u32)
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn params(method: Method) -> Params {
        let mut classes = [255u8; 256];
        for (i, c) in b"ACDEFGHIKLMNPQRSTVWY".iter().enumerate() {
            classes[usize::from(*c)] = i as u8;
        }
        Params {
            method,
            k: 3,
            classes,
            minimizer: 1,
            modulo: 1,
            quantum: 1,
            band: 1,
            gap: 4,
            gap_bin: 1,
            min_gap: 1,
            combo: 2,
            top: 1,
        }
    }

    #[test]
    fn scores_follow_the_shared_diagonal() {
        let rep: &[u8] = b"MKVLATGHEWRSPQ";
        let read: &[u8] = b"QQKVLATGHE"; // rep[1..9] after two others
        let other: &[u8] = b"HEATGVLKQQ";
        let mut m = vec![-1i8; 1 << 16];
        for c in b'A'..=b'Z' {
            m[usize::from(c) << 8 | usize::from(c)] = 2;
        }
        let score = |method, read: &[u8]| {
            pseudo_scores(&[read], &[rep], &[0], &[0], &params(method), &m)[0]
        };
        assert_eq!(score(Method::Diag, read), 6); // KVL VLA LAT ATG TGH GHE
        assert_eq!(score(Method::Order, read), 6);
        assert_eq!(score(Method::Anchored, read), 16);
        assert_eq!(score(Method::Gapped, read), 16); // no gap needed: the ungapped score
        assert!(score(Method::Pairs, read) > score(Method::Pairs, other));
        assert_eq!(score(Method::Order, other), 1); // ATG
        // a deletion in the read: the gapped scan bridges it (10 + 9 matches, a gap of 12)
        let long: &[u8] = b"MKVLATWQPRGHEWRSPQDN";
        let del: &[u8] = b"MKVLATWQPRHEWRSPQDN"; // without its G
        let one = |method| pseudo_scores(&[del], &[long], &[0], &[0], &params(method), &m)[0];
        assert_eq!(one(Method::Anchored), 20); // the longer side only
        assert_eq!(one(Method::Gapped), 2 * 19 - 12);
        // repeats: the shortest period with four copies, as verify.mask
        assert_eq!(mask_repeats(b"MKVLTHTHTHTHAQ"), b"MKVLxxxxxxxxAQ");
        assert_eq!(mask_repeats(b"MKVAAAAGCVCVCVCVQ"), b"MKVxxxxGxxxxxxxxQ");
        assert_eq!(mask_repeats(b"GCVCVCVCQ"), b"GCVCVCVCQ");
        // dense routing: the read votes for the leaf sharing its diagonal, not the other
        let (r, l, v) = route_votes(&[read], &[other, rep], &params(Method::Diag), 2, 4, 100);
        assert_eq!((r, l, v), (vec![0], vec![1], vec![6]));
        // a representative with no k-mer (masked, or shorter than k) scores 0
        let mut p = params(Method::Anchored);
        p.minimizer = 4;
        assert_eq!(
            pseudo_scores(&[read], &[b"xxxxMK" as &[u8]], &[0], &[0], &p, &m)[0],
            0
        );
    }
}
