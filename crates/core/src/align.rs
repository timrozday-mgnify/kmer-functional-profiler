//! Pairwise protein identity by overlap alignment (free end gaps), for benchmarks that need
//! a distance independent of k-mers (docs/tree-index-plan.md, phase 0: edge cost by identity).

use crate::threads;

const MATCH: i32 = 2;
const MISMATCH: i32 = -1;
const GAP: i32 = -2;

/// Identity of `a` and `b` aligned with free end gaps: matches over the columns between the
/// first and last aligned residue pair, and that column count (0 when nothing aligns).
pub fn identity(a: &[u8], b: &[u8]) -> (f32, u32) {
    let (n, m) = (a.len(), b.len());
    let width = m + 1;
    // score and traceback (0 diagonal, 1 up = gap in b, 2 left = gap in a)
    let mut score = vec![0i32; (n + 1) * width];
    let mut from = vec![0u8; (n + 1) * width];
    for i in 1..=n {
        for j in 1..=m {
            let s = if a[i - 1].eq_ignore_ascii_case(&b[j - 1]) {
                MATCH
            } else {
                MISMATCH
            };
            let diag = score[(i - 1) * width + j - 1] + s;
            let up = score[(i - 1) * width + j] + GAP;
            let left = score[i * width + j - 1] + GAP;
            let (best, dir) = if diag >= up && diag >= left {
                (diag, 0)
            } else if up >= left {
                (up, 1)
            } else {
                (left, 2)
            };
            score[i * width + j] = best;
            from[i * width + j] = dir;
        }
    }
    // free end gaps: the best cell on the last row or column
    let (mut i, mut j) = (n, m);
    let mut best = i32::MIN;
    for jj in 0..=m {
        if score[n * width + jj] > best {
            (best, i, j) = (score[n * width + jj], n, jj);
        }
    }
    for ii in 0..=n {
        if score[ii * width + m] > best {
            (best, i, j) = (score[ii * width + m], ii, m);
        }
    }
    let (mut matches, mut columns, mut tail) = (0u32, 0u32, 0u32);
    let mut started = false; // columns after the last aligned pair are not counted
    while i > 0 && j > 0 {
        match from[i * width + j] {
            0 => {
                started = true;
                columns += 1 + tail;
                tail = 0;
                matches += u32::from(a[i - 1].eq_ignore_ascii_case(&b[j - 1]));
                (i, j) = (i - 1, j - 1);
            }
            1 => {
                tail += u32::from(started);
                i -= 1;
            }
            _ => {
                tail += u32::from(started);
                j -= 1;
            }
        }
    }
    if columns == 0 {
        return (0.0, 0);
    }
    (matches as f32 / columns as f32, columns)
}

/// [`identity`] of each pair (`a[i]`, `b[i]`), in parallel.
pub fn identities<S: AsRef<[u8]> + Sync>(a: &[S], b: &[S]) -> (Vec<f32>, Vec<u32>) {
    let n = a.len().min(b.len());
    let mut ident = vec![0f32; n];
    let mut cols = vec![0u32; n];
    let per_thread = n.div_ceil(threads()).max(1);
    std::thread::scope(|scope| {
        for (t, (id, co)) in ident
            .chunks_mut(per_thread)
            .zip(cols.chunks_mut(per_thread))
            .enumerate()
        {
            scope.spawn(move || {
                for k in 0..id.len() {
                    let i = t * per_thread + k;
                    (id[k], co[k]) = identity(a[i].as_ref(), b[i].as_ref());
                }
            });
        }
    });
    (ident, cols)
}

/// Best ungapped local alignment score of `a` against `b` over every diagonal (the maximum
/// segment sum along each), residue pairs scored by `matrix` (256 x 256 by byte, row `a`'s).
pub fn ungapped(a: &[u8], b: &[u8], matrix: &[i8]) -> i32 {
    let mut best = 0;
    for start in 0..a.len() + b.len() {
        // diagonals from the top row (a from 0, b from start - a.len()) then down the first
        // column (a from start, b from 0)
        let (mut i, mut j) = if start < a.len() {
            (a.len() - 1 - start, 0)
        } else {
            (0, start - a.len())
        };
        let mut run = 0i32;
        while i < a.len() && j < b.len() {
            run = (run + i32::from(matrix[usize::from(a[i]) << 8 | usize::from(b[j])])).max(0);
            best = best.max(run);
            i += 1;
            j += 1;
        }
    }
    best
}

/// [`ungapped`] for each pair `(a[i], b[i])`, over the configured threads.
pub fn ungapped_many<S: AsRef<[u8]> + Sync>(a: &[S], b: &[S], matrix: &[i8]) -> Vec<i32> {
    let n = a.len().min(b.len());
    let mut out = vec![0i32; n];
    let per_thread = n.div_ceil(threads()).max(1);
    std::thread::scope(|scope| {
        for (t, chunk) in out.chunks_mut(per_thread).enumerate() {
            scope.spawn(move || {
                for (k, o) in chunk.iter_mut().enumerate() {
                    let i = t * per_thread + k;
                    *o = ungapped(a[i].as_ref(), b[i].as_ref(), matrix);
                }
            });
        }
    });
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ungapped_finds_the_best_diagonal() {
        let mut m = vec![-1i8; 1 << 16];
        for c in b'A'..=b'Z' {
            m[usize::from(c) << 8 | usize::from(c)] = 2;
        }
        assert_eq!(ungapped(b"MKVLAT", b"MKVLAT", &m), 12);
        assert_eq!(ungapped(b"QQMKVLQQ", b"MKVL", &m), 8); // the shifted diagonal
        assert_eq!(ungapped(b"MKXVL", b"MKAVL", &m), 7); // through a mismatch: 2+2-1+2+2
        assert_eq!(ungapped(b"", b"MK", &m), 0);
    }

    #[test]
    fn identity_of_simple_pairs() {
        assert_eq!(identity(b"MKVLAT", b"MKVLAT"), (1.0, 6));
        assert_eq!(identity(b"MKVLAT", b"MKALAT"), (5.0 / 6.0, 6));
        // a fragment aligns with free end gaps
        assert_eq!(identity(b"VLA", b"MKVLATQ"), (1.0, 3));
        assert_eq!(identity(b"", b"MK"), (0.0, 0));
    }
}
