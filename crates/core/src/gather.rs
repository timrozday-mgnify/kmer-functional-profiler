//! Greedy gather: assign each hit k-mer to one unit, as `sourmash gather` does.

use std::cmp::Ordering;
use std::collections::BinaryHeap;

use crate::Error;

/// Units gather takes, in the order taken (rank = position).
#[derive(Debug, Default, PartialEq, Eq)]
pub struct Gathered {
    pub unit: Vec<u32>,
    /// Hit k-mers assigned to the unit: its k-mers no earlier unit took.
    pub kmers_unique: Vec<u32>,
}

/// Heap entry, popped highest score first, then lowest unit.
struct Entry {
    score: f64,
    unit: u32,
}

impl Ord for Entry {
    fn cmp(&self, other: &Self) -> Ordering {
        self.score
            .total_cmp(&other.score)
            .then(other.unit.cmp(&self.unit))
    }
}

impl PartialOrd for Entry {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl PartialEq for Entry {
    fn eq(&self, other: &Self) -> bool {
        self.cmp(other) == Ordering::Equal
    }
}

impl Eq for Entry {}

/// Gather over (unit, hash) pairs: repeatedly take the unit with the most untaken hit
/// k-mers divided by its sampling rate `t_g[unit]`, and give it those k-mers. Units left
/// with none are explained away. Scores only fall as k-mers are taken, so a lazy heap
/// rechecks only a stale top. Repeated pairs count once.
pub fn gather(units: &[u32], hashes: &[u64], t_g: &[f64]) -> Result<Gathered, Error> {
    if units.len() != hashes.len() {
        return Err(Error::LengthMismatch);
    }
    let n = t_g.len();
    // Dense k-mer ids: rank of each pair's hash among the distinct hashes.
    let mut by_hash: Vec<(u64, usize)> = hashes.iter().copied().zip(0..).collect();
    by_hash.sort_unstable();
    let mut kmer = vec![0u32; hashes.len()];
    let mut n_kmers = 0u32;
    for (i, &(hash, pair)) in by_hash.iter().enumerate() {
        if i > 0 && hash != by_hash[i - 1].0 {
            n_kmers += 1;
        }
        kmer[pair] = n_kmers;
    }
    drop(by_hash);
    // Each unit's k-mers, CSR: members[start[u]..end[u]] are its untaken ones.
    let mut start = vec![0usize; n + 1];
    for &u in units {
        let u = u as usize;
        if u >= n {
            return Err(Error::UnitOutOfRange(u as u64));
        }
        start[u + 1] += 1;
    }
    for u in 0..n {
        start[u + 1] += start[u];
    }
    let mut end = start[..n].to_vec();
    let mut members = vec![0u32; units.len()];
    for (&u, &k) in units.iter().zip(&kmer) {
        members[end[u as usize]] = k;
        end[u as usize] += 1;
    }
    drop(kmer);
    let mut heap = BinaryHeap::new();
    for u in 0..n {
        let mine = &mut members[start[u]..end[u]];
        mine.sort_unstable();
        let mut len = 0;
        for j in 0..mine.len() {
            if j == 0 || mine[j] != mine[j - 1] {
                mine[len] = mine[j];
                len += 1;
            }
        }
        end[u] = start[u] + len;
        if len > 0 {
            #[allow(clippy::cast_possible_truncation)] // n <= t_g.len() units, ids fit u32
            heap.push(Entry {
                score: len as f64 / t_g[u],
                unit: u as u32,
            });
        }
    }
    let mut taken = vec![false; n_kmers as usize + 1];
    let mut out = Gathered::default();
    while let Some(Entry { unit, .. }) = heap.pop() {
        let u = unit as usize;
        let mine = &mut members[start[u]..end[u]];
        let mut len = 0;
        for j in 0..mine.len() {
            if !taken[mine[j] as usize] {
                mine[len] = mine[j];
                len += 1;
            }
        }
        end[u] = start[u] + len;
        if len == 0 {
            continue;
        }
        let score = len as f64 / t_g[u];
        if heap.peek().is_some_and(|top| score < top.score) {
            heap.push(Entry { score, unit });
            continue;
        }
        for &k in &members[start[u]..end[u]] {
            taken[k as usize] = true;
        }
        out.unit.push(unit);
        #[allow(clippy::cast_possible_truncation)] // a unit has far fewer than 2^32 k-mers
        out.kmers_unique.push(len as u32);
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    const UNITS: [u32; 10] = [0, 0, 0, 0, 1, 1, 2, 2, 3, 3];
    const HASHES: [u64; 10] = [1, 2, 3, 4, 3, 4, 4, 5, 6, 7];

    #[test]
    fn explains_away_shared_kmers() {
        // 1 is a subset of 0: explained away. 2 keeps only hash 5; 3 shares nothing.
        let got = gather(&UNITS, &HASHES, &[0.01; 4]).unwrap();
        assert_eq!(got.unit, [0, 3, 2]);
        assert_eq!(got.kmers_unique, [4, 2, 1]);
        // Sampled 10x more sparsely, unit 1's 2 k-mers stand for more than unit 0's 4.
        let got = gather(&UNITS, &HASHES, &[0.01, 0.001, 0.01, 0.01]).unwrap();
        assert_eq!(got.unit, [1, 0, 3, 2]);
        assert_eq!(got.kmers_unique, [2, 2, 2, 1]);
    }

    #[test]
    fn ties_go_to_the_lowest_unit_and_repeats_count_once() {
        let got = gather(&[2, 1, 1, 2], &[9, 9, 9, 9], &[1.0; 3]).unwrap();
        assert_eq!((got.unit, got.kmers_unique), (vec![1], vec![1]));
    }

    #[test]
    fn rejects_bad_input() {
        assert!(gather(&[0], &[], &[1.0]).is_err());
        assert!(gather(&[1], &[5], &[1.0]).is_err());
    }
}
