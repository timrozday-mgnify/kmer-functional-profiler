//! Read routing (`kfp-pileup route`): the distinct k-mers each read pair shares with each
//! block of a reference, through a [`PackedTable`] whose set values are block ids.

use crate::{Error, PackedTable, threads};

/// One row per (pair, block it shares a k-mer with): the distinct k-mers shared.
#[derive(Debug, Default, PartialEq, Eq)]
pub struct RouteHits {
    pub read: Vec<u64>,
    pub block: Vec<u32>,
    pub kmers: Vec<u32>,
}

/// Hashes per thread below which `route_hits` does not split further.
const MIN_PER_THREAD: usize = 1 << 16;

/// Count, per pair in `reads` and block, the distinct `hashes` the table gives that block.
/// `reads` must be non-decreasing (a pair's hashes contiguous, as [`crate::FastxHits`]
/// yields them); rows come in read order, blocks ascending within a read.
///
/// Runs of whole pairs go to separate threads, joined in order.
pub fn route_hits(
    table: &PackedTable<'_>,
    reads: &[u64],
    hashes: &[u64],
) -> Result<RouteHits, Error> {
    if hashes.len() != reads.len() {
        return Err(Error::LengthMismatch);
    }
    if reads.windows(2).any(|w| w[1] < w[0]) {
        return Err(Error::UnsortedReads);
    }
    // chunk ends moved forward to a pair boundary, so no pair spans two threads
    let per_thread = reads.len().div_ceil(threads()).max(MIN_PER_THREAD);
    let mut bounds = vec![0];
    while let Some(&lo) = bounds.last().filter(|&&lo| lo < reads.len()) {
        let mut hi = (lo + per_thread).min(reads.len());
        while hi < reads.len() && reads[hi] == reads[hi - 1] {
            hi += 1;
        }
        bounds.push(hi);
    }
    let parts: Vec<RouteHits> = std::thread::scope(|scope| {
        let handles: Vec<_> = bounds
            .windows(2)
            .map(|w| {
                let (r, h) = (&reads[w[0]..w[1]], &hashes[w[0]..w[1]]);
                scope.spawn(move || route_serial(table, r, h))
            })
            .collect();
        handles
            .into_iter()
            .map(|h| h.join().expect("route thread panicked"))
            .collect()
    });
    let mut out = RouteHits::default();
    for part in parts {
        out.read.extend(part.read);
        out.block.extend(part.block);
        out.kmers.extend(part.kmers);
    }
    Ok(out)
}

fn route_serial(table: &PackedTable<'_>, reads: &[u64], hashes: &[u64]) -> RouteHits {
    let mut out = RouteHits::default();
    let mut seen: Vec<(u32, u64)> = Vec::new(); // (block, hash) of the current pair
    let mut flush = |read: u64, seen: &mut Vec<(u32, u64)>| {
        seen.sort_unstable();
        seen.dedup();
        for run in seen.chunk_by(|a, b| a.0 == b.0) {
            out.read.push(read);
            out.block.push(run[0].0);
            #[allow(clippy::cast_possible_truncation)] // a pair has far fewer than 2^32 k-mers
            out.kmers.push(run.len() as u32);
        }
        seen.clear();
    };
    for (i, (&read, &hash)) in reads.iter().zip(hashes).enumerate() {
        if i > 0 && read != reads[i - 1] {
            flush(reads[i - 1], &mut seen);
        }
        if let Some(set) = table.lookup(hash) {
            #[allow(clippy::cast_possible_truncation)] // block ids are far below 2^32
            seen.extend(table.values(set).map(|block| (block as u32, hash)));
        }
    }
    if let Some(&last) = reads.last() {
        flush(last, &mut seen);
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Column;

    // As packed.rs's tests: keys are the top 8 bits, 4-bit bucket and fingerprint.
    // (0, 3) -> set 0 = {block 1}, (0, 9) -> set 1 = {blocks 0, 2}, (2, 5) -> set 0.
    fn table() -> PackedTable<'static> {
        PackedTable {
            max_hash: u64::MAX >> 1,
            shift: 56,
            fp_bits: 4,
            offsets: Column::U8(&[0, 2, 2, 3, 3]),
            fingerprints: Column::U16(&[3, 9, 5]),
            set_ids: Column::U32(&[0, 1, 0]),
            set_offsets: Column::U8(&[0, 1, 3]),
            set_values: Column::U64(&[1, 0, 2]),
        }
    }

    fn hash(bucket: u64, fp: u64) -> u64 {
        ((bucket << 4) | fp) << 56
    }

    #[test]
    fn distinct_kmers_per_pair_and_block() {
        let t = table();
        let (h03, h09, h25, miss) = (hash(0, 3), hash(0, 9), hash(2, 5), hash(0, 4));
        // pair 4: h03 twice (one distinct k-mer), h25 (block 1 again), h09 (blocks 0, 2);
        // pair 5: misses only; pair 6: h09
        let reads = [4, 4, 4, 4, 5, 6];
        let hashes = [h03, h03, h25, h09, miss, h09];
        let got = route_hits(&t, &reads, &hashes).unwrap();
        assert_eq!(got.read, [4, 4, 4, 6, 6]);
        assert_eq!(got.block, [0, 1, 2, 0, 2]);
        assert_eq!(got.kmers, [1, 2, 1, 1, 1]);
    }

    #[test]
    fn threaded_equals_serial_and_keeps_pairs_whole() {
        let t = table();
        let n = 5 * MIN_PER_THREAD + 7;
        // pairs of 3 hashes, so chunk ends fall inside pairs
        let reads: Vec<u64> = (0..n as u64).map(|i| i / 3).collect();
        let hashes: Vec<u64> = (0..n as u64).map(|i| hash(i % 4, (i / 4) % 16)).collect();
        let serial = route_serial(&t, &reads, &hashes);
        assert!(!serial.read.is_empty());
        assert_eq!(route_hits(&t, &reads, &hashes).unwrap(), serial);
    }

    #[test]
    fn rejects_unsorted_reads_and_length_mismatch() {
        let t = table();
        assert!(route_hits(&t, &[2, 1], &[0, 0]).is_err());
        assert!(route_hits(&t, &[0], &[]).is_err());
    }
}
