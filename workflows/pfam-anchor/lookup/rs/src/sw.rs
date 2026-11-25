//! Smith-Waterman, BLOSUM62, gaps 11 + 1 a residue (DIAMOND's and blastp's defaults).
use crate::data::{blosum, AA, NONE};

pub const OPEN: i16 = 12; // first residue of a gap
pub const EXT: i16 = 1;
pub const LANES: usize = 16;
type V = [i16; LANES];

#[inline(always)]
fn max(a: V, b: V) -> V {
    std::array::from_fn(|k| if a[k] > b[k] { a[k] } else { b[k] })
}
#[inline(always)]
fn sub(a: V, b: i16) -> V {
    std::array::from_fn(|k| a[k] - b)
}

/// Best local score of every query against each of up to 16 targets: `out[query][lane]`.
pub fn batch(targets: &[&[u8]], queries: &[&[u8]], out: &mut [V]) {
    let bl = blosum();
    let n = targets.iter().map(|t| t.len()).max().unwrap_or(0);
    let m: usize = queries.iter().map(|q| q.len()).sum();
    let mut h = vec![[0i16; LANES]; m];
    let mut e = vec![[-1000i16; LANES]; m];
    out.iter_mut().for_each(|o| *o = [0; LANES]);
    for j in 0..n {
        let col: [u8; LANES] =
            std::array::from_fn(|k| targets.get(k).and_then(|t| t.get(j)).copied().unwrap_or(NONE));
        let sp: [V; AA] = std::array::from_fn(|a| std::array::from_fn(|k| bl[a][col[k] as usize]));
        let mut i = 0;
        for (q, o) in queries.iter().zip(out.iter_mut()) {
            let (mut diag, mut f, mut best) = ([0i16; LANES], [-1000i16; LANES], *o);
            for &c in q.iter() {
                let s = &sp[c as usize];
                let old = h[i];
                let ei = max(sub(e[i], EXT), sub(old, OPEN));
                let d: V = std::array::from_fn(|k| diag[k] + s[k]);
                let v = max(max(max(d, [0; LANES]), ei), f);
                best = max(best, v);
                f = max(sub(f, EXT), sub(v, OPEN));
                diag = old;
                h[i] = v;
                e[i] = ei;
                i += 1;
            }
            *o = best;
        }
    }
}

/// Score, identities, alignment columns and target start of the best local alignment, with
/// traceback. With a band: only cells with lo <= j - i <= hi (target minus query position).
pub fn align(q: &[u8], t: &[u8], bl: &[[i16; AA]; AA], band: Option<(i32, i32)>) -> (i16, u32, u32, usize) {
    let (m, n) = (q.len(), t.len());
    let w = n + 1;
    const BIG: i16 = -20000;
    let mut h = vec![0i16; (m + 1) * w];
    let mut e = vec![BIG; (m + 1) * w]; // gap in query (moves along target)
    let mut f = vec![BIG; (m + 1) * w];
    let (mut best, mut bi, mut bj) = (0, 0, 0);
    for i in 1..=m {
        let (lo, hi) = match band {
            Some((a, b)) => ((i as i32 + a).max(1) as usize, (i as i32 + b).min(n as i32)),
            None => (1, n as i32),
        };
        for j in lo..=(hi.max(0) as usize) {
            let x = i * w + j;
            e[x] = (e[x - 1] - EXT).max(h[x - 1] - OPEN);
            f[x] = (f[x - w] - EXT).max(h[x - w] - OPEN);
            let v = (h[x - w - 1] + bl[q[i - 1] as usize][t[j - 1] as usize]).max(0).max(e[x]).max(f[x]);
            h[x] = v;
            if v > best {
                (best, bi, bj) = (v, i, j);
            }
        }
    }
    let (mut ident, mut cols) = (0, 0);
    let (mut i, mut j) = (bi, bj);
    let mut state = 0;
    while i > 0 && j > 0 {
        let x = i * w + j;
        match state {
            0 => {
                if h[x] == 0 {
                    break;
                }
                if h[x] == e[x] {
                    state = 1;
                } else if h[x] == f[x] {
                    state = 2;
                } else {
                    ident += (q[i - 1] == t[j - 1] && q[i - 1] < 20) as u32;
                    cols += 1;
                    i -= 1;
                    j -= 1;
                }
            }
            1 => {
                cols += 1;
                if e[x] == h[x - 1] - OPEN {
                    state = 0;
                }
                j -= 1;
            }
            _ => {
                cols += 1;
                if f[x] == h[x - w] - OPEN {
                    state = 0;
                }
                i -= 1;
            }
        }
    }
    (best, ident, cols, j)
}

/// Score only, linear memory; with a band as in `align`.
pub fn score(q: &[u8], t: &[u8], bl: &[[i16; AA]; AA], band: Option<(i32, i32)>, buf: &mut Vec<i16>) -> i16 {
    let (m, n) = (q.len(), t.len());
    buf.clear();
    buf.resize(2 * (n + 1), 0);
    let (h, f) = buf.split_at_mut(n + 1);
    f.fill(-20000);
    let mut best = 0;
    for i in 1..=m {
        let (lo, hi) = match band {
            Some((a, b)) => ((i as i32 + a).max(1), (i as i32 + b).min(n as i32)),
            None => (1, n as i32),
        };
        if lo > hi {
            continue;
        }
        let (lo, hi) = (lo as usize, hi as usize);
        let s = &bl[q[i - 1] as usize];
        let mut diag = h[lo - 1];
        let mut e = -20000i16;
        let mut left = 0i16;
        if band.is_some() {
            h[lo - 1] = 0; // outside the band from here on
        }
        for j in lo..=hi {
            let old = h[j];
            e = (e - EXT).max(left - OPEN);
            f[j] = (f[j] - EXT).max(old - OPEN);
            let v = (diag + s[t[j - 1] as usize]).max(0).max(e).max(f[j]);
            h[j] = v;
            left = v;
            diag = old;
            best = best.max(v);
        }
    }
    best
}

#[cfg(test)]
mod tests {
    use super::*;
    fn rnd(state: &mut u64, n: usize) -> Vec<u8> {
        (0..n)
            .map(|_| {
                *state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
                ((*state >> 33) % 20) as u8
            })
            .collect()
    }
    #[test]
    fn batch_equals_align_and_score() {
        let bl = blosum();
        let mut s = 7u64;
        let mut t: Vec<Vec<u8>> = (0..16).map(|k| rnd(&mut s, 30 + 11 * k)).collect();
        let q: Vec<Vec<u8>> = (0..9).map(|k| rnd(&mut s, 20 + 5 * k)).collect();
        // homologous with gaps
        let mut x = q[8].clone();
        x.drain(20..23);
        t[3] = [rnd(&mut s, 40), x, rnd(&mut s, 12)].concat();
        let (tr, qr): (Vec<&[u8]>, Vec<&[u8]>) = (t.iter().map(|x| &x[..]).collect(), q.iter().map(|x| &x[..]).collect());
        let mut out = vec![[0i16; LANES]; 9];
        batch(&tr, &qr, &mut out);
        let mut buf = vec![];
        for (i, q) in qr.iter().enumerate() {
            for (k, t) in tr.iter().enumerate() {
                let a = align(q, t, &bl, None);
                assert_eq!(out[i][k], a.0);
                assert_eq!(out[i][k], score(q, t, &bl, None, &mut buf));
                assert!(a.1 <= a.2);
                assert!(score(q, t, &bl, Some((-8, 8)), &mut buf) <= a.0);
                let wide = Some((-(q.len() as i32), t.len() as i32));
                assert_eq!(score(q, t, &bl, wide, &mut buf), a.0);
                assert_eq!(align(q, t, &bl, Some((-8, 8))).0, score(q, t, &bl, Some((-8, 8)), &mut buf));
            }
        }
        assert!(out[8][3] > 200);
    }
}
