//! Connected components of units linked by shared k-mers (the query's `components` stage).

use crate::Error;

/// Hash-range passes over the pairs: each sorts ~1/`PASSES` of them, so memory is the
/// union-find (4 bytes a unit) plus 1/`PASSES` of the pairs (12 bytes each).
const PASSES: u32 = 16;

fn find(parent: &mut [u32], mut u: u32) -> u32 {
    while parent[u as usize] != u {
        let up = parent[parent[u as usize] as usize];
        parent[u as usize] = up; // path halving
        u = up;
    }
    u
}

/// Component label of each of units `0..n`, linked when they hold a hash in common among
/// the (`units`, `hashes`) pairs. Components are numbered in order of their smallest unit;
/// a unit in no pair is a component of its own.
pub fn components(units: &[u32], hashes: &[u64], n: usize) -> Result<Vec<u32>, Error> {
    if units.len() != hashes.len() {
        return Err(Error::LengthMismatch);
    }
    if let Some(&u) = units.iter().find(|&&u| u as usize >= n) {
        return Err(Error::UnitOutOfRange(u64::from(u)));
    }
    let size = u32::try_from(n).map_err(|_| Error::UnitOutOfRange(n as u64))?;
    let mut parent: Vec<u32> = (0..size).collect();
    let mut part: Vec<(u64, u32)> = Vec::new();
    for pass in 0..u64::from(PASSES) {
        part.clear();
        part.extend(
            hashes
                .iter()
                .zip(units)
                .filter(|&(&h, _)| h % u64::from(PASSES) == pass)
                .map(|(&h, &u)| (h, u)),
        );
        part.sort_unstable();
        for pair in part.windows(2) {
            if pair[0].0 == pair[1].0 {
                let (a, b) = (find(&mut parent, pair[0].1), find(&mut parent, pair[1].1));
                if a != b {
                    // the smaller root wins, so a root is its component's smallest unit
                    parent[a.max(b) as usize] = a.min(b);
                }
            }
        }
    }
    let mut label = vec![u32::MAX; n];
    let mut next = 0;
    for u in 0..size {
        let root = find(&mut parent, u) as usize;
        if label[root] == u32::MAX {
            label[root] = next;
            next += 1;
        }
        label[u as usize] = label[root];
    }
    Ok(label)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn links_units_sharing_a_hash() {
        // 0-2 share 7, 2-3 share 9; 1 alone; 4 in no pair.
        let got = components(&[0, 2, 2, 3, 1], &[7, 7, 9, 9, 5], 5).unwrap();
        assert_eq!(got, [0, 1, 0, 0, 2]);
    }

    #[test]
    fn rejects_bad_input() {
        assert!(components(&[0], &[], 1).is_err());
        assert!(components(&[3], &[1], 2).is_err());
    }
}
