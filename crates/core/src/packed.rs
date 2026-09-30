//! Lookup in a packed tier-2 table (`PackedTable` in the Python index), read in place.
//!
//! A hash's key is `hash >> shift`: its top `bucket_bits` address `offsets` directly and
//! its low `fp_bits` are the stored fingerprint. Bucket `b` holds the sorted fingerprints
//! `fingerprints[offsets[b]..offsets[b + 1]]` (~4 on average) and their value-set ids; set
//! `s` is `set_values[set_offsets[s]..set_offsets[s + 1]]`, each value `unit << PIN_BITS |
//! pin_q`. So a lookup touches one bucket and one set, and no keys are materialised: the
//! arrays can be memory-mapped `.npy` files.

use crate::{Error, threads};

/// Low bits of a set value holding the k-mer's quantised `p_in` in that unit.
pub const PIN_BITS: u32 = 4;

/// An unsigned array of any width, read as `u64`. The index stores each array in the
/// smallest dtype that holds its values.
#[derive(Clone, Copy, Debug)]
pub enum Column<'a> {
    U8(&'a [u8]),
    U16(&'a [u16]),
    U32(&'a [u32]),
    U64(&'a [u64]),
}

impl Column<'_> {
    #[inline]
    fn get(&self, i: usize) -> u64 {
        match self {
            Self::U8(v) => u64::from(v[i]),
            Self::U16(v) => u64::from(v[i]),
            Self::U32(v) => u64::from(v[i]),
            Self::U64(v) => v[i],
        }
    }

    #[inline]
    fn len(&self) -> usize {
        match self {
            Self::U8(v) => v.len(),
            Self::U16(v) => v.len(),
            Self::U32(v) => v.len(),
            Self::U64(v) => v.len(),
        }
    }
}

/// A packed hash -> value-set table over borrowed arrays.
#[derive(Clone, Copy, Debug)]
pub struct PackedTable<'a> {
    pub max_hash: u64,
    /// `64 - lead_bits - bucket_bits - fp_bits`.
    pub shift: u32,
    pub fp_bits: u32,
    pub offsets: Column<'a>,
    pub fingerprints: Column<'a>,
    pub set_ids: Column<'a>,
    pub set_offsets: Column<'a>,
    pub set_values: Column<'a>,
}

impl PackedTable<'_> {
    /// Value-set id of `hash`, or `None`. False hits (another key with the same bucket and
    /// fingerprint) occur at ~keys per bucket x 2^-fp_bits per lookup.
    #[inline]
    pub fn lookup(&self, hash: u64) -> Option<usize> {
        if hash > self.max_hash {
            return None;
        }
        let key = hash.checked_shr(self.shift).unwrap_or(0);
        let bucket = usize::try_from(key.checked_shr(self.fp_bits).unwrap_or(0)).ok()?;
        if bucket + 1 >= self.offsets.len() {
            return None;
        }
        let fp = key & (1u64.checked_shl(self.fp_bits).unwrap_or(0).wrapping_sub(1));
        let (start, end) = (self.offsets.get(bucket), self.offsets.get(bucket + 1));
        // Fingerprints are sorted within a bucket, but a bucket holds only a few.
        (start as usize..end as usize)
            .find(|&i| self.fingerprints.get(i) == fp)
            .map(|i| self.set_ids.get(i) as usize)
    }

    /// Values of set `s`.
    #[inline]
    fn values(&self, s: usize) -> impl Iterator<Item = u64> + '_ {
        (self.set_offsets.get(s) as usize..self.set_offsets.get(s + 1) as usize)
            .map(|i| self.set_values.get(i))
    }
}

/// Hits of sampled query hashes, one row per (hash occurrence, unit it counts for).
#[derive(Debug, Default, PartialEq, Eq)]
pub struct UnitHits {
    pub unit: Vec<u32>,
    pub hash: Vec<u64>,
    pub read: Vec<u64>,
    /// The k-mer's quantised `p_in` in that unit.
    pub pin_q: Vec<u8>,
    /// Units the hash counts for (the rows it has).
    pub holders: Vec<u32>,
}

/// Hashes per thread below which `unit_hits` does not split further: a thread's start-up
/// is worth ~10^4 lookups.
const MIN_PER_THREAD: usize = 1 << 14;

/// Expand `hashes` (with their `reads`) to unit hits through `table`, in input order and
/// set order. A hash counts for a unit only if it passes that unit's `max_hash_g`, which
/// also discards most fingerprint false hits.
///
/// Lookups are random reads into the table, bound by memory latency, so contiguous runs of
/// hashes go to separate threads and their rows are joined in order.
pub fn unit_hits(
    table: &PackedTable<'_>,
    max_hash_g: &[u64],
    hashes: &[u64],
    reads: &[u64],
) -> Result<UnitHits, Error> {
    if hashes.len() != reads.len() {
        return Err(Error::LengthMismatch);
    }
    let per_thread = hashes.len().div_ceil(threads()).max(MIN_PER_THREAD);
    let parts = std::thread::scope(|scope| {
        let handles: Vec<_> = hashes
            .chunks(per_thread)
            .zip(reads.chunks(per_thread))
            .map(|(h, r)| scope.spawn(move || unit_hits_serial(table, max_hash_g, h, r)))
            .collect();
        handles
            .into_iter()
            .map(|h| h.join().expect("lookup thread panicked"))
            .collect::<Result<Vec<_>, _>>()
    })?;
    let mut out = UnitHits::default();
    for part in parts {
        out.unit.extend(part.unit);
        out.hash.extend(part.hash);
        out.read.extend(part.read);
        out.pin_q.extend(part.pin_q);
        out.holders.extend(part.holders);
    }
    Ok(out)
}

fn unit_hits_serial(
    table: &PackedTable<'_>,
    max_hash_g: &[u64],
    hashes: &[u64],
    reads: &[u64],
) -> Result<UnitHits, Error> {
    let mut out = UnitHits::default();
    let mut kept: Vec<(u32, u8)> = Vec::new();
    for (&hash, &read) in hashes.iter().zip(reads) {
        let Some(set) = table.lookup(hash) else {
            continue;
        };
        kept.clear();
        for value in table.values(set) {
            let unit = value >> PIN_BITS;
            let limit = usize::try_from(unit)
                .ok()
                .and_then(|u| max_hash_g.get(u))
                .ok_or(Error::UnitOutOfRange(unit))?;
            if hash <= *limit {
                #[allow(clippy::cast_possible_truncation)] // unit ids fit u32; p_in is 4 bits
                kept.push((unit as u32, (value & ((1 << PIN_BITS) - 1)) as u8));
            }
        }
        #[allow(clippy::cast_possible_truncation)] // a set holds far fewer than 2^32 units
        let holders = kept.len() as u32;
        for &(unit, pin_q) in &kept {
            out.unit.push(unit);
            out.hash.push(hash);
            out.read.push(read);
            out.pin_q.push(pin_q);
            out.holders.push(holders);
        }
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    // Keys are the top 8 bits (shift 56): 4-bit bucket, 4-bit fingerprint; buckets 0-3.
    // (bucket, fp) -> set: (0, 3) -> sets[0], (0, 9) -> sets[1], (2, 5) -> sets[2].
    fn table<'a>(sets: &'a [u32], values: &'a [u64]) -> PackedTable<'a> {
        PackedTable {
            max_hash: u64::MAX >> 1,
            shift: 56,
            fp_bits: 4,
            offsets: Column::U8(&[0, 2, 2, 3, 3]),
            fingerprints: Column::U16(&[3, 9, 5]),
            set_ids: Column::U32(sets),
            set_offsets: Column::U8(&[0, 2, 3]),
            set_values: Column::U64(values),
        }
    }

    fn hash(bucket: u64, fp: u64) -> u64 {
        ((bucket << 4) | fp) << 56
    }

    #[test]
    fn lookup_finds_bucket_and_fingerprint() {
        let t = table(&[0, 1, 0], &[0, 0, 0]);
        assert_eq!(t.lookup(hash(0, 3)), Some(0));
        assert_eq!(t.lookup(hash(0, 9) | 0xff), Some(1)); // low bits below the key ignored
        assert_eq!(t.lookup(hash(2, 5)), Some(0));
        assert_eq!(t.lookup(hash(0, 4)), None);
        assert_eq!(t.lookup(hash(1, 3)), None); // empty bucket
        assert_eq!(t.lookup(u64::MAX), None); // above max_hash
    }

    #[test]
    fn unit_hits_expands_sets_and_checks_max_hash_g() {
        // Set 0: unit 1 (pin 2), unit 2 (pin 15); set 1: unit 0 (pin 0).
        let values = [(1 << PIN_BITS) | 2, (2 << PIN_BITS) | 15, 0];
        let t = table(&[0, 1, 0], &values);
        let (h03, h09, h25) = (hash(0, 3), hash(0, 9), hash(2, 5));
        // Unit 2 keeps only hashes up to h03, so h25 counts for unit 1 alone.
        let max_hash_g = [u64::MAX, u64::MAX, h03];
        let got = unit_hits(
            &t,
            &max_hash_g,
            &[h03, hash(0, 4), h25, h09],
            &[7, 8, 9, 10],
        )
        .unwrap();
        assert_eq!(got.unit, [1, 2, 1, 0]);
        assert_eq!(got.hash, [h03, h03, h25, h09]);
        assert_eq!(got.read, [7, 7, 9, 10]);
        assert_eq!(got.pin_q, [2, 15, 2, 0]);
        assert_eq!(got.holders, [2, 2, 1, 1]);
    }

    #[test]
    fn unit_hits_threaded_equals_serial() {
        let values = [(1 << PIN_BITS) | 2, (2 << PIN_BITS) | 15, 0];
        let t = table(&[0, 1, 0], &values);
        let max_hash_g = [u64::MAX, u64::MAX, hash(0, 3)];
        // Several chunks of every bucket and fingerprint, hits and misses mixed.
        let n = 5 * MIN_PER_THREAD + 7;
        let hashes: Vec<u64> = (0..n as u64).map(|i| hash(i % 4, (i / 4) % 16)).collect();
        let reads: Vec<u64> = (0..n as u64).collect();
        let serial = unit_hits_serial(&t, &max_hash_g, &hashes, &reads).unwrap();
        assert!(!serial.unit.is_empty());
        assert_eq!(unit_hits(&t, &max_hash_g, &hashes, &reads).unwrap(), serial);
    }

    #[test]
    fn unit_hits_rejects_units_beyond_max_hash_g() {
        let t = table(&[0, 1, 0], &[(5 << PIN_BITS), 0, 0]);
        assert!(unit_hits(&t, &[u64::MAX], &[hash(0, 3)], &[0]).is_err());
        assert!(unit_hits(&t, &[u64::MAX], &[hash(0, 3)], &[]).is_err());
    }
}
