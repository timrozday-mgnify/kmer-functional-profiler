//! Seed alignment: a pseudo-alignment made into a full alignment (docs/pfam-anchor-plan.md,
//! "Tier 2 as a pseudo-alignment"). Exact protein k-mer seeds between a query stretch and a
//! target, the ones within `band` diagonals of the best-supported diagonal, chained increasing
//! in both sequences; then the gaps between chained pairs filled (along the diagonal where both
//! sides skip alike, else an affine global DP of the two stretches) and the chain extended
//! ungapped from both ends with X-drop. The Python oracle is `anchor.py`'s `seed_chain` and
//! `build_alignment`.

use std::collections::HashMap;

use crate::threads;

/// A target's k-mers to their positions.
type KmerIndex<'a> = HashMap<&'a [u8], Vec<u32>>;

/// Affine gap costs, DIAMOND's (a gap of g residues costs `GAP_OPEN + g * GAP_EXTEND`).
pub const GAP_OPEN: i32 = 11;
pub const GAP_EXTEND: i32 = 1;

/// Settings of [`seed_align`].
#[derive(Clone, Copy, Debug)]
pub struct SeedParams {
    pub k: usize,
    /// Diagonals kept on either side of the best-supported one.
    pub band: i64,
    /// Score below the best at which an ungapped extension stops.
    pub x_drop: i32,
}

/// One alignment: residue ranges (0-based, end exclusive) in the query stretch and target, its
/// raw score (substitutions and gaps) and its BLAST traceback (BTOP: runs of identities as
/// numbers, a substitution as query then target residue, `X-` a query-only residue, `-X` a
/// target-only one).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SeedAlignment {
    pub query_start: u32,
    pub query_end: u32,
    pub target_start: u32,
    pub target_end: u32,
    pub raw: i32,
    pub btop: String,
}

fn score(matrix: &[i8], a: u8, b: u8) -> i32 {
    i32::from(matrix[usize::from(a) << 8 | usize::from(b)])
}

/// Target k-mers to their positions.
fn index(target: &[u8], k: usize) -> KmerIndex<'_> {
    let mut out: KmerIndex<'_> = HashMap::new();
    if target.len() >= k {
        for y in 0..=target.len() - k {
            out.entry(&target[y..y + k]).or_default().push(y as u32);
        }
    }
    out
}

/// The chained residue pairs (query, target), increasing in both.
fn chain(query: &[u8], idx: &KmerIndex<'_>, p: &SeedParams) -> Vec<(i64, i64)> {
    if query.len() < p.k {
        return Vec::new();
    }
    let mut hits: Vec<(i64, i64)> = Vec::new();
    for x in 0..=query.len() - p.k {
        if let Some(ys) = idx.get(&query[x..x + p.k]) {
            hits.extend(ys.iter().map(|&y| (x as i64, i64::from(y))));
        }
    }
    if hits.is_empty() {
        return Vec::new();
    }
    // the best-supported diagonal: most seeds, the smallest diagonal on a tie (as np.unique)
    let mut counts: HashMap<i64, u32> = HashMap::new();
    for &(x, y) in &hits {
        *counts.entry(y - x).or_default() += 1;
    }
    let best = counts
        .iter()
        .max_by(|a, b| a.1.cmp(b.1).then(b.0.cmp(a.0)))
        .map(|(&d, _)| d)
        .unwrap();
    let mut pairs: Vec<(i64, i64)> = hits
        .iter()
        .filter(|&&(x, y)| (y - x - best).abs() <= p.band)
        .flat_map(|&(x, y)| (0..p.k as i64).map(move |k| (x + k, y + k)))
        .collect();
    pairs.sort_unstable();
    pairs.dedup();
    // longest chain strictly increasing in both: by x, then y descending, so a strictly
    // increasing run of y can't hold two pairs with the same x (patience sort)
    pairs.sort_unstable_by(|a, b| a.0.cmp(&b.0).then(b.1.cmp(&a.1)));
    let (mut tails, mut tail_y, mut prev) = (Vec::new(), Vec::new(), vec![usize::MAX; pairs.len()]);
    for (n, &(_, y)) in pairs.iter().enumerate() {
        let lo = tail_y.partition_point(|&t| t < y);
        prev[n] = if lo > 0 { tails[lo - 1] } else { usize::MAX };
        if lo == tails.len() {
            tails.push(n);
            tail_y.push(y);
        } else {
            tails[lo] = n;
            tail_y[lo] = y;
        }
    }
    let mut out = Vec::with_capacity(tails.len());
    let mut n = tails.last().copied().unwrap_or(usize::MAX);
    while n != usize::MAX {
        out.push(pairs[n]);
        n = prev[n];
    }
    out.reverse();
    out
}

/// Global alignment of `a` and `b` with affine gaps (Gotoh), as residue pairs relative to their
/// starts. Ties go to a substitution, then a query-only gap.
fn global(a: &[u8], b: &[u8], matrix: &[i8]) -> Vec<(i64, i64)> {
    const NEG: i32 = i32::MIN / 4;
    let (n, m) = (a.len(), b.len());
    let w = m + 1;
    // h: best ending in a pair; e: ending with a target-only gap; f: with a query-only gap
    let (mut h, mut e, mut f) = (
        vec![NEG; (n + 1) * w],
        vec![NEG; (n + 1) * w],
        vec![NEG; (n + 1) * w],
    );
    h[0] = 0;
    for (j, cell) in e.iter_mut().enumerate().take(m + 1).skip(1) {
        *cell = -(GAP_OPEN + j as i32 * GAP_EXTEND);
    }
    for i in 1..=n {
        f[i * w] = -(GAP_OPEN + i as i32 * GAP_EXTEND);
        for j in 1..=m {
            let diag = h[(i - 1) * w + j - 1]
                .max(e[(i - 1) * w + j - 1])
                .max(f[(i - 1) * w + j - 1]);
            h[i * w + j] = diag + score(matrix, a[i - 1], b[j - 1]);
            e[i * w + j] = (h[i * w + j - 1] - GAP_OPEN - GAP_EXTEND)
                .max(f[i * w + j - 1] - GAP_OPEN - GAP_EXTEND)
                .max(e[i * w + j - 1] - GAP_EXTEND);
            f[i * w + j] = (h[(i - 1) * w + j] - GAP_OPEN - GAP_EXTEND)
                .max(e[(i - 1) * w + j] - GAP_OPEN - GAP_EXTEND)
                .max(f[(i - 1) * w + j] - GAP_EXTEND);
        }
    }
    let (mut i, mut j) = (n, m);
    let at = |s: &Vec<i32>, i: usize, j: usize| s[i * w + j];
    let mut state = {
        let (x, y, z) = (at(&h, n, m), at(&f, n, m), at(&e, n, m));
        if x >= y && x >= z {
            0
        } else if y >= z {
            1
        } else {
            2
        }
    };
    let mut out = Vec::new();
    while i > 0 || j > 0 {
        match state {
            0 if i > 0 && j > 0 => {
                out.push((i as i64 - 1, j as i64 - 1));
                let (x, y, z) = (
                    at(&h, i - 1, j - 1),
                    at(&f, i - 1, j - 1),
                    at(&e, i - 1, j - 1),
                );
                state = if x >= y && x >= z {
                    0
                } else if y >= z {
                    1
                } else {
                    2
                };
                i -= 1;
                j -= 1;
            }
            1 if i > 0 => {
                let cur = at(&f, i, j);
                state = if cur == at(&f, i - 1, j) - GAP_EXTEND && i > 1 {
                    1
                } else if cur == at(&h, i - 1, j) - GAP_OPEN - GAP_EXTEND {
                    0
                } else {
                    2
                }; // fmt: skip
                i -= 1;
            }
            _ if j > 0 => {
                let cur = at(&e, i, j);
                state = if cur == at(&e, i, j - 1) - GAP_EXTEND && j > 1 {
                    2
                } else if cur == at(&h, i, j - 1) - GAP_OPEN - GAP_EXTEND {
                    0
                } else {
                    1
                }; // fmt: skip
                j -= 1;
            }
            _ => {
                state = if i > 0 { 1 } else { 2 };
            }
        }
    }
    out.reverse();
    out
}

/// The full alignment of `query` (one stop-free stretch) against `target`, from its seed chain;
/// None with fewer than two chained pairs.
fn align_one(
    query: &[u8],
    target: &[u8],
    idx: &KmerIndex<'_>,
    p: &SeedParams,
    matrix: &[i8],
) -> Option<SeedAlignment> {
    let seeds = chain(query, idx, p);
    if seeds.len() < 2 {
        return None;
    }
    let mut pairs: Vec<(i64, i64)> = vec![seeds[0]];
    for w in seeds.windows(2) {
        let ((x1, y1), (x2, y2)) = (w[0], w[1]);
        let (dx, dy) = (x2 - x1, y2 - y1);
        if dx == dy {
            pairs.extend((1..dx).map(|k| (x1 + k, y1 + k)));
        } else if dx > 1 && dy > 1 {
            let a = &query[(x1 + 1) as usize..x2 as usize];
            let b = &target[(y1 + 1) as usize..y2 as usize];
            pairs.extend(
                global(a, b, matrix)
                    .into_iter()
                    .map(|(i, j)| (x1 + 1 + i, y1 + 1 + j)),
            );
        }
        pairs.push((x2, y2));
    }
    // ungapped X-drop extension from both ends
    let extend = |x0: i64, y0: i64, step: i64| -> Vec<(i64, i64)> {
        let (mut s, mut best, mut keep, mut k) = (0, 0, 0, 1);
        while (0..query.len() as i64).contains(&(x0 + step * k))
            && (0..target.len() as i64).contains(&(y0 + step * k))
        {
            s += score(
                matrix,
                query[(x0 + step * k) as usize],
                target[(y0 + step * k) as usize],
            );
            if s > best {
                (best, keep) = (s, k);
            } else if best - s > p.x_drop {
                break;
            }
            k += 1;
        }
        (1..=keep).map(|m| (x0 + step * m, y0 + step * m)).collect()
    };
    let (first, last) = (pairs[0], *pairs.last().unwrap());
    let mut left = extend(first.0, first.1, -1);
    left.reverse();
    let right = extend(last.0, last.1, 1);
    let pairs: Vec<(i64, i64)> = left.into_iter().chain(pairs).chain(right).collect();
    // BTOP and raw score
    let (mut btop, mut run, mut raw) = (String::new(), 0u32, 0i32);
    let flush = |btop: &mut String, run: &mut u32| {
        if *run > 0 {
            btop.push_str(&run.to_string());
            *run = 0;
        }
    };
    for (n, &(x, y)) in pairs.iter().enumerate() {
        if n > 0 {
            let (px, py) = pairs[n - 1];
            let (gq, gt) = (x - px - 1, y - py - 1);
            for (g, query_side) in [(gq, true), (gt, false)] {
                if g > 0 {
                    flush(&mut btop, &mut run);
                    raw -= GAP_OPEN + g as i32 * GAP_EXTEND;
                    for r in 1..=g {
                        if query_side {
                            btop.push(query[(px + r) as usize] as char);
                            btop.push('-');
                        } else {
                            btop.push('-');
                            btop.push(target[(py + r) as usize] as char);
                        }
                    }
                }
            }
        }
        let (a, b) = (query[x as usize], target[y as usize]);
        raw += score(matrix, a, b);
        if a == b {
            run += 1;
        } else {
            flush(&mut btop, &mut run);
            btop.push(a as char);
            btop.push(b as char);
        }
    }
    flush(&mut btop, &mut run);
    let (first, last) = (pairs[0], *pairs.last().unwrap());
    Some(SeedAlignment {
        query_start: first.0 as u32,
        query_end: last.0 as u32 + 1,
        target_start: first.1 as u32,
        target_end: last.1 as u32 + 1,
        raw,
        btop,
    })
}

/// Each pair's best alignment: query `queries[pair_query[i]]` (stop-free stretches joined by
/// `*`, as a read's frames are) against target `targets[pair_target[i]]`, the stretch with the
/// longest seed chain; coordinates on the joined query. Pairs sorted by target index each
/// target once per run of it.
pub fn seed_align<S: AsRef<[u8]> + Sync>(
    queries: &[S],
    targets: &[S],
    pair_query: &[u32],
    pair_target: &[u32],
    p: &SeedParams,
    matrix: &[i8],
) -> Vec<Option<SeedAlignment>> {
    let n = pair_query.len().min(pair_target.len());
    let mut out: Vec<Option<SeedAlignment>> = vec![None; n];
    let per_thread = n.div_ceil(threads()).max(1);
    std::thread::scope(|scope| {
        for (t, chunk) in out.chunks_mut(per_thread).enumerate() {
            scope.spawn(move || {
                let mut cur: Option<(u32, KmerIndex<'_>)> = None;
                for (k, o) in chunk.iter_mut().enumerate() {
                    let i = t * per_thread + k;
                    let q = pair_target[i];
                    let target = targets[q as usize].as_ref();
                    if cur.as_ref().is_none_or(|c| c.0 != q) {
                        cur = Some((q, index(target, p.k)));
                    }
                    let idx = &cur.as_ref().unwrap().1;
                    let query = queries[pair_query[i] as usize].as_ref();
                    let mut best: Option<(usize, SeedAlignment)> = None;
                    let mut offset = 0usize;
                    for stretch in query.split(|&c| c == b'*') {
                        let chained = chain(stretch, idx, p).len();
                        if chained >= 2 && best.as_ref().is_none_or(|b| chained > b.0) {
                            if let Some(mut a) = align_one(stretch, target, idx, p, matrix) {
                                a.query_start += offset as u32;
                                a.query_end += offset as u32;
                                best = Some((chained, a));
                            }
                        }
                        offset += stretch.len() + 1;
                    }
                    *o = best.map(|b| b.1);
                }
            });
        }
    });
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn blosum_like() -> Vec<i8> {
        // identity +5, substitution -1 (enough to test the mechanics)
        let mut m = vec![-1i8; 1 << 16];
        for c in 0..256usize {
            m[c << 8 | c] = 5;
        }
        m
    }

    fn p() -> SeedParams {
        SeedParams {
            k: 4,
            band: 8,
            x_drop: 15,
        }
    }

    #[test]
    fn exact_query_aligns_end_to_end() {
        let m = blosum_like();
        let got = seed_align(
            &[b"ACDEFGHIK".as_slice()],
            &[b"WWACDEFGHIKLM".as_slice()],
            &[0],
            &[0],
            &p(),
            &m,
        );
        let a = got[0].clone().unwrap();
        assert_eq!(
            (a.query_start, a.query_end, a.target_start, a.target_end),
            (0, 9, 2, 11)
        );
        assert_eq!(a.btop, "9");
        assert_eq!(a.raw, 45);
    }

    #[test]
    fn substitution_and_indel_are_in_the_btop() {
        let m = blosum_like();
        // query has one extra residue (Y) between two seeded stretches, and a substitution
        let query = b"ACDEFGHYIKLMNPQRSTVW";
        let target = b"ACDEFGHIKLMNPQRSAVW";
        let a = seed_align(
            &[query.as_slice()],
            &[target.as_slice()],
            &[0],
            &[0],
            &p(),
            &m,
        )[0]
        .clone()
        .unwrap();
        assert_eq!(a.btop, "7Y-9TA2");
        assert_eq!(a.raw, 18 * 5 - 1 - (GAP_OPEN + GAP_EXTEND));
    }

    #[test]
    fn the_stretch_with_the_longest_chain_wins_and_keeps_its_offset() {
        let m = blosum_like();
        let a = seed_align(
            &[b"WWWW*ACDEFGH".as_slice()],
            &[b"ACDEFGH".as_slice()],
            &[0],
            &[0],
            &p(),
            &m,
        )[0]
        .clone()
        .unwrap();
        assert_eq!((a.query_start, a.query_end), (5, 12));
    }

    #[test]
    fn global_places_a_gap_between_anchors() {
        let m = blosum_like();
        assert_eq!(global(b"AXC", b"AC", &m), vec![(0, 0), (2, 1)]);
    }
}
