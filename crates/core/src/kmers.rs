use crate::translate::six_frames_into;
use crate::{Alphabet, Error, FrameMode, GeneticCode};

/// Hashes a packed k-mer with the splitmix64 finalizer.
///
/// A bijection on `u64`, so distinct k-mers of one (alphabet, k) never collide.
// ponytail: splitmix64 mixes well enough for threshold sampling; swap for xxh3 if
// sampled fractions ever show bias against the binomial expectation.
#[inline]
pub fn hash_kmer(kmer: u64) -> u64 {
    let mut z = kmer.wrapping_add(0x9E37_79B9_7F4A_7C15);
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

/// FracMinHash threshold keeping a `fraction` of hashes: a hash is kept iff `hash <= max_hash`.
///
/// `fraction >= 1` keeps everything; the smallest threshold (0) still keeps hash 0.
pub fn max_hash(fraction: f64) -> u64 {
    if fraction >= 1.0 {
        u64::MAX
    } else {
        // `as` saturates: negative and NaN fractions give 0.
        ((fraction * 2f64.powi(64)) as u64).saturating_sub(1)
    }
}

/// Validated k-mer extraction settings shared by the protein and DNA paths.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct KmerParams {
    k: usize,
    alphabet: Alphabet,
    code: GeneticCode,
    frames: FrameMode,
    max_hash: u64,
}

impl KmerParams {
    /// Settings for k-mers of length `k`; genetic code 11, stop-free frames, no sampling.
    ///
    /// # Errors
    /// [`Error::InvalidK`] if `k` is 0 or the packed k-mer would not fit in a `u64`.
    pub fn new(k: usize, alphabet: Alphabet) -> Result<Self, Error> {
        let max = alphabet.max_k();
        if !(1..=max).contains(&k) {
            return Err(Error::InvalidK { k, max });
        }
        Ok(Self {
            k,
            alphabet,
            code: GeneticCode::default(),
            frames: FrameMode::default(),
            max_hash: u64::MAX,
        })
    }

    #[must_use]
    pub fn with_genetic_code(mut self, code: GeneticCode) -> Self {
        self.code = code;
        self
    }

    #[must_use]
    pub fn with_frames(mut self, frames: FrameMode) -> Self {
        self.frames = frames;
        self
    }

    #[must_use]
    pub fn with_max_hash(mut self, max_hash: u64) -> Self {
        self.max_hash = max_hash;
        self
    }

    pub fn k(&self) -> usize {
        self.k
    }

    pub fn alphabet(&self) -> Alphabet {
        self.alphabet
    }

    pub fn genetic_code(&self) -> GeneticCode {
        self.code
    }

    pub fn frames(&self) -> FrameMode {
        self.frames
    }

    pub fn max_hash(&self) -> u64 {
        self.max_hash
    }
}

/// Calls `emit` with the hash of every k-mer of `protein` that passes the threshold, in order.
///
/// K-mers never span a byte outside the alphabet (stop `*`, `X`, ambiguity codes).
#[inline]
pub fn protein_kmers(protein: &[u8], params: &KmerParams, mut emit: impl FnMut(u64)) {
    let bits = params.alphabet.bits();
    let width = params.k as u32 * bits;
    let mask = if width == u64::BITS {
        u64::MAX
    } else {
        (1 << width) - 1
    };
    let (mut packed, mut run) = (0u64, 0usize);
    for &aa in protein {
        let Some(code) = params.alphabet.encode(aa) else {
            run = 0;
            continue;
        };
        packed = ((packed << bits) | u64::from(code)) & mask;
        run += 1;
        if run >= params.k {
            let hash = hash_kmer(packed);
            if hash <= params.max_hash {
                emit(hash);
            }
        }
    }
}

/// Sampled k-mer hits as parallel columns, one row per kept k-mer.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Hits {
    /// Read (pair) index, counted from 0 across the whole input.
    pub read: Vec<u64>,
    /// 0 for the first mate or a single-end read, 1 for the second mate.
    pub mate: Vec<u8>,
    /// Reading frame, 0–2 forward and 3–5 reverse complement.
    pub frame: Vec<u8>,
    pub hash: Vec<u64>,
}

impl Hits {
    pub fn len(&self) -> usize {
        self.hash.len()
    }

    pub fn is_empty(&self) -> bool {
        self.hash.is_empty()
    }
}

/// Translates reads and collects their sampled k-mer hashes, reusing buffers across reads.
#[derive(Debug, Clone)]
pub struct DnaScanner {
    params: KmerParams,
    rc: Vec<u8>,
    frames: [Vec<u8>; 6],
}

impl DnaScanner {
    pub fn new(params: KmerParams) -> Self {
        Self {
            params,
            rc: Vec::new(),
            frames: Default::default(),
        }
    }

    pub fn params(&self) -> &KmerParams {
        &self.params
    }

    /// Appends the hits of one read's kept frames to `hits`.
    pub fn scan(&mut self, dna: &[u8], read: u64, mate: u8, hits: &mut Hits) {
        six_frames_into(dna, self.params.code, &mut self.rc, &mut self.frames);
        let stop_free = self.params.frames == FrameMode::StopFree;
        for (frame, aa) in (0u8..).zip(&self.frames) {
            if stop_free && aa.contains(&b'*') {
                continue;
            }
            protein_kmers(aa, &self.params, |hash| {
                hits.read.push(read);
                hits.mate.push(mate);
                hits.frame.push(frame);
                hits.hash.push(hash);
            });
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;

    fn scan(dna: &[u8], params: KmerParams) -> Hits {
        let mut hits = Hits::default();
        DnaScanner::new(params).scan(dna, 0, 0, &mut hits);
        hits
    }

    #[test]
    fn protein_yields_one_hash_per_window_and_splits_at_stops() {
        let params = KmerParams::new(5, Alphabet::Protein).unwrap();
        let mut n = 0;
        protein_kmers(b"MKVLAAGIVL*MKVLAX", &params, |_| n += 1);
        assert_eq!(n, 6 + 1);
    }

    #[test]
    fn rejects_k_that_does_not_fit() {
        assert!(KmerParams::new(12, Alphabet::Protein).is_ok());
        assert!(KmerParams::new(13, Alphabet::Protein).is_err());
        assert!(KmerParams::new(16, Alphabet::Murphy10).is_ok());
        assert!(KmerParams::new(0, Alphabet::Dayhoff).is_err());
    }

    #[test]
    fn max_hash_bounds() {
        assert_eq!(max_hash(1.0), u64::MAX);
        assert_eq!(max_hash(0.5), (1 << 63) - 1);
        assert_eq!(max_hash(0.0), 0);
        assert_eq!(max_hash(f64::NAN), 0);
    }

    #[test]
    fn stop_filter_drops_frames_with_stops() {
        let params = KmerParams::new(3, Alphabet::Protein).unwrap();
        // Frame 0: M A W * (dropped by the filter); other frames are stop-free.
        let dna = b"ATGGCCTGGTAA";
        let all = scan(dna, params.with_frames(FrameMode::All));
        let filtered = scan(dna, params);
        assert!(all.frame.contains(&0));
        assert!(!filtered.frame.contains(&0));
    }

    /// Codons for each amino acid under genetic code 11.
    fn synonymous_codons(aa: u8) -> Vec<[u8; 3]> {
        let bases = *b"TCAG";
        let mut out = Vec::new();
        for a in bases {
            for b in bases {
                for c in bases {
                    let mut p = Vec::new();
                    crate::translate(&[a, b, c], GeneticCode::Bacterial, &mut p);
                    if p[0] == aa {
                        out.push([a, b, c]);
                    }
                }
            }
        }
        out
    }

    fn back_translate(protein: &str, picks: &[u8]) -> Vec<u8> {
        protein
            .bytes()
            .zip(picks)
            .flat_map(|(aa, &pick)| {
                let codons = synonymous_codons(aa);
                codons[pick as usize % codons.len()]
            })
            .collect()
    }

    proptest! {
        #[test]
        fn synonymous_changes_leave_frame_zero_hashes_unchanged(
            protein in "[ACDEFGHIKLMNPQRSTVWY]{0,80}",
            picks1 in prop::collection::vec(any::<u8>(), 80),
            picks2 in prop::collection::vec(any::<u8>(), 80),
            k in 1usize..=12,
        ) {
            let params = KmerParams::new(k, Alphabet::Protein).unwrap().with_frames(FrameMode::All);
            let frame0 = |hits: Hits| -> Vec<u64> {
                hits.frame.iter().zip(hits.hash).filter(|(f, _)| **f == 0).map(|(_, h)| h).collect()
            };
            let a = frame0(scan(&back_translate(&protein, &picks1), params));
            let b = frame0(scan(&back_translate(&protein, &picks2), params));
            prop_assert_eq!(a.len(), protein.len().saturating_sub(k - 1));
            prop_assert_eq!(a, b);
        }

        #[test]
        fn lower_thresholds_keep_a_subset(
            dna in "[ACGTN]{0,300}",
            t1 in any::<u64>(),
            t2 in any::<u64>(),
            k in 1usize..=16,
        ) {
            let (lo, hi) = (t1.min(t2), t1.max(t2));
            let params = KmerParams::new(k, Alphabet::Murphy10).unwrap();
            let sparse = scan(dna.as_bytes(), params.with_max_hash(lo));
            let dense = scan(dna.as_bytes(), params.with_max_hash(hi));
            let expected: Vec<_> = dense.frame.iter().zip(&dense.hash).filter(|(_, h)| **h <= lo).collect();
            let got: Vec<_> = sparse.frame.iter().zip(&sparse.hash).collect();
            prop_assert_eq!(got, expected);
        }
    }
}
