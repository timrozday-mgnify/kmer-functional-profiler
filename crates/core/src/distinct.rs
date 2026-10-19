//! Distinct `u64` hashes from batches, kept sorted in hash ranges (the query's sketch of a
//! sample's sampled hashes, plan phase 6, step 40).
//!
//! Each range holds a sorted, deduplicated `total` and a `pending` run of batch hashes,
//! merged into the total once pending outgrows a quarter of it: memory is ~1.25x the
//! distinct hashes, and a merge's transient is one range's.

use crate::threads;

/// Range bits: 256 ranges, so a merge's transient is ~1/256 of the total.
const RANGE_BITS: u32 = 8;
/// Pending hashes below which a range never merges (a merge's fixed cost).
const MIN_PENDING: usize = 1 << 16;

#[derive(Debug, Default)]
struct Range {
    total: Vec<u64>,
    pending: Vec<u64>,
}

impl Range {
    fn merge(&mut self) {
        if self.pending.is_empty() {
            return;
        }
        self.pending.sort_unstable();
        self.pending.dedup();
        let mut out = Vec::with_capacity(self.total.len() + self.pending.len());
        let (mut a, mut b) = (self.total.iter().peekable(), self.pending.iter().peekable());
        loop {
            let next = match (a.peek(), b.peek()) {
                (Some(&&x), Some(&&y)) => {
                    if x <= y {
                        a.next();
                        if x == y {
                            b.next();
                        }
                        x
                    } else {
                        b.next();
                        y
                    }
                }
                (Some(&&x), None) => {
                    a.next();
                    x
                }
                (None, Some(&&y)) => {
                    b.next();
                    y
                }
                (None, None) => break,
            };
            out.push(next);
        }
        self.total = out;
        self.pending = Vec::new();
    }
}

/// Distinct hashes at most `max_hash`, in `2^RANGE_BITS` consecutive hash ranges.
#[derive(Debug)]
pub struct Distinct {
    ranges: Vec<Range>,
    shift: u32,
}

impl Distinct {
    #[must_use]
    pub fn new(max_hash: u64) -> Self {
        Self {
            ranges: (0..1 << RANGE_BITS).map(|_| Range::default()).collect(),
            shift: (64 - max_hash.leading_zeros()).saturating_sub(RANGE_BITS),
        }
    }

    /// Number of ranges.
    #[must_use]
    pub fn len(&self) -> usize {
        self.ranges.len()
    }

    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.ranges.is_empty()
    }

    fn range(&self, hash: u64) -> usize {
        usize::try_from(hash >> self.shift)
            .unwrap_or(usize::MAX)
            .min(self.ranges.len() - 1)
    }

    /// Add a batch of hashes (any order, repeats allowed); ranges merge in parallel.
    pub fn add(&mut self, hashes: &[u64]) {
        for &h in hashes {
            let r = self.range(h);
            self.ranges[r].pending.push(h);
        }
        let per_thread = self.ranges.len().div_ceil(threads());
        std::thread::scope(|scope| {
            for chunk in self.ranges.chunks_mut(per_thread) {
                scope.spawn(move || {
                    for range in chunk {
                        if range.pending.len() > (range.total.len() / 4).max(MIN_PENDING) {
                            range.merge();
                        }
                    }
                });
            }
        });
    }

    /// Range `i`'s distinct hashes, sorted; the range is emptied (so each is taken once).
    pub fn take(&mut self, i: usize) -> Vec<u64> {
        let range = &mut self.ranges[i];
        range.merge();
        std::mem::take(&mut range.total)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;

    proptest! {
        #[test]
        fn takes_the_sorted_distinct_hashes_range_by_range(
            batches in proptest::collection::vec(proptest::collection::vec(0..1u64 << 20, 0..3000), 0..8),
        ) {
            let mut d = Distinct::new((1 << 20) - 1);
            for b in &batches {
                d.add(b);
            }
            let got: Vec<u64> = (0..d.len()).flat_map(|i| d.take(i)).collect();
            let mut want: Vec<u64> = batches.concat();
            want.sort_unstable();
            want.dedup();
            prop_assert_eq!(got, want);
        }
    }

    #[test]
    fn merges_keep_the_total_distinct() {
        let mut d = Distinct::new(u64::MAX);
        let batch: Vec<u64> = (0..200_000u64)
            .map(|i| i.wrapping_mul(0x9E37_79B9_7F4A_7C15))
            .collect();
        d.add(&batch);
        d.add(&batch); // repeats merge away
        let got: Vec<u64> = (0..d.len()).flat_map(|i| d.take(i)).collect();
        let mut want = batch;
        want.sort_unstable();
        assert_eq!(got, want);
        assert!(d.take(0).is_empty()); // taken once
    }
}
