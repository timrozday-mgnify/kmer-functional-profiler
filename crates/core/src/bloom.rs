//! Blocked Bloom filter over k-mer hashes, for the partitioned index build.
//!
//! The filter is a byte array of 64-byte blocks. Each key sets `PROBES` bits inside one
//! block chosen from its hash, so an insert or a lookup touches one cache line; at 10 bits
//! per key the false-positive rate is about 1% (0.8% for an unblocked filter).

use crate::Error;
use crate::kmers::hash_kmer;
use crate::threads;

/// Bytes per block (one cache line).
pub const BLOCK_BYTES: usize = 64;
/// Bits set per key.
pub const PROBES: u32 = 7;
const PROBE_BITS: u32 = 9; // log2(BLOCK_BYTES * 8)

/// Block index and the in-block bit positions of `hash`.
#[inline]
fn locate(hash: u64, n_blocks: u64) -> (u64, u64) {
    // Candidate hashes all lie below t_max, so remix before taking bits.
    let h = hash_kmer(hash);
    let block = ((u128::from(h) * u128::from(n_blocks)) >> 64) as u64;
    (block, hash_kmer(h))
}

#[inline]
fn bit(probes: u64, i: u32) -> (usize, u8) {
    let b = (probes >> (i * PROBE_BITS)) & ((1 << PROBE_BITS) - 1);
    ((b >> 3) as usize, 1 << (b & 7))
}

fn n_blocks(bits: &[u8]) -> Result<u64, Error> {
    if bits.is_empty() || bits.len() % BLOCK_BYTES != 0 {
        return Err(Error::BloomSize(bits.len()));
    }
    Ok((bits.len() / BLOCK_BYTES) as u64)
}

/// Adds `hashes` to the filter `bits` (a non-empty multiple of `BLOCK_BYTES` long).
///
/// Each thread owns a disjoint run of blocks and scans every hash, setting only bits in its
/// own blocks: no atomics, and the random writes are spread over the threads.
pub fn bloom_insert(bits: &mut [u8], hashes: &[u64]) -> Result<(), Error> {
    let n = n_blocks(bits)?;
    let per_thread = (n as usize).div_ceil(threads());
    std::thread::scope(|scope| {
        for (t, part) in bits.chunks_mut(per_thread * BLOCK_BYTES).enumerate() {
            let lo = (t * per_thread) as u64;
            let hi = lo + (part.len() / BLOCK_BYTES) as u64;
            scope.spawn(move || {
                for &hash in hashes {
                    let (block, probes) = locate(hash, n);
                    if (lo..hi).contains(&block) {
                        let base = (block - lo) as usize * BLOCK_BYTES;
                        for i in 0..PROBES {
                            let (byte, mask) = bit(probes, i);
                            part[base + byte] |= mask;
                        }
                    }
                }
            });
        }
    });
    Ok(())
}

/// Whether each hash may be in the filter (never false for an inserted hash).
pub fn bloom_contains(bits: &[u8], hashes: &[u64]) -> Result<Vec<bool>, Error> {
    let n = n_blocks(bits)?;
    let mut out = vec![false; hashes.len()];
    let per_thread = hashes.len().div_ceil(threads()).max(1);
    std::thread::scope(|scope| {
        for (keys, found) in hashes.chunks(per_thread).zip(out.chunks_mut(per_thread)) {
            scope.spawn(move || {
                for (&hash, f) in keys.iter().zip(found) {
                    let (block, probes) = locate(hash, n);
                    let base = block as usize * BLOCK_BYTES;
                    *f = (0..PROBES).all(|i| {
                        let (byte, mask) = bit(probes, i);
                        bits[base + byte] & mask != 0
                    });
                }
            });
        }
    });
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Distinct small values (like hashes below t_max) for disjoint `range`s.
    fn keys(range: std::ops::Range<u64>) -> Vec<u64> {
        range.map(|i| hash_kmer(i) >> 20).collect()
    }

    #[test]
    fn no_false_negatives_and_about_one_percent_false_hits() {
        let inserted = keys(0..50_000);
        let mut bits = vec![0u8; (50_000 * 10 / 8usize).div_ceil(BLOCK_BYTES) * BLOCK_BYTES];
        bloom_insert(&mut bits, &inserted).unwrap();
        assert!(bloom_contains(&bits, &inserted).unwrap().iter().all(|&f| f));
        let others = keys(50_000..250_000);
        let rate = bloom_contains(&bits, &others)
            .unwrap()
            .iter()
            .filter(|&&f| f)
            .count() as f64
            / 200_000.0;
        assert!((0.005..0.02).contains(&rate), "false-positive rate {rate}");
    }

    #[test]
    fn threaded_insert_matches_sequential() {
        let hashes = keys(0..10_000);
        let mut bits = vec![0u8; 97 * BLOCK_BYTES]; // not a multiple of the thread count
        bloom_insert(&mut bits, &hashes).unwrap();
        let mut expected = vec![0u8; bits.len()];
        for &h in &hashes {
            let (block, probes) = locate(h, 97);
            for i in 0..PROBES {
                let (byte, mask) = bit(probes, i);
                expected[block as usize * BLOCK_BYTES + byte] |= mask;
            }
        }
        assert_eq!(bits, expected);
    }

    #[test]
    fn rejects_partial_blocks() {
        assert!(bloom_insert(&mut [0u8; 65], &[1]).is_err());
        assert!(bloom_contains(&[], &[1]).is_err());
    }
}
