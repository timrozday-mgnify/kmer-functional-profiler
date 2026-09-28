use std::str::FromStr;

use crate::Error;

/// Amino acids for codons in NCBI order (T, C, A, G at each position), genetic code 11.
const CODE_11: &[u8; 64] = b"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
/// Genetic code 4: as code 11 except TGA codes for Trp.
const CODE_4: &[u8; 64] = b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";

/// Marks non-ACGTU bases; any codon containing one translates to `X`.
const INVALID: u8 = u8::MAX;
const BASE: [u8; 256] = {
    let mut table = [INVALID; 256];
    let bases = b"TCAGUtcagu";
    let codes = [0, 1, 2, 3, 0, 0, 1, 2, 3, 0];
    let mut i = 0;
    while i < bases.len() {
        table[bases[i] as usize] = codes[i];
        i += 1;
    }
    table
};

/// NCBI translation table. Start codons are irrelevant here: reads are never ORF-called.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
#[non_exhaustive]
pub enum GeneticCode {
    /// Table 11: bacteria, archaea, plastids.
    #[default]
    Bacterial,
    /// Table 4: Mycoplasma and Spiroplasma (TGA = Trp).
    Mycoplasma,
}

impl GeneticCode {
    fn table(self) -> &'static [u8; 64] {
        match self {
            Self::Bacterial => CODE_11,
            Self::Mycoplasma => CODE_4,
        }
    }
}

impl TryFrom<u8> for GeneticCode {
    type Error = Error;

    fn try_from(table: u8) -> Result<Self, Error> {
        match table {
            11 => Ok(Self::Bacterial),
            4 => Ok(Self::Mycoplasma),
            _ => Err(Error::UnknownGeneticCode(table)),
        }
    }
}

/// Which of the six reading frames of a read are kept.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
#[non_exhaustive]
pub enum FrameMode {
    /// Only frames without a stop codon: ~1.5 of 6 frames for a coding read at moderate GC.
    #[default]
    StopFree,
    /// All six frames; k-mers are still split at stops.
    All,
}

impl FromStr for FrameMode {
    type Err = Error;

    fn from_str(s: &str) -> Result<Self, Error> {
        match s.to_ascii_lowercase().as_str() {
            "stopfree" => Ok(Self::StopFree),
            "all" => Ok(Self::All),
            _ => Err(Error::UnknownFrameMode(s.to_owned())),
        }
    }
}

/// Appends the translation of every full codon of `dna` to `out`; stops are `*`, unknown `X`.
pub fn translate(dna: &[u8], code: GeneticCode, out: &mut Vec<u8>) {
    let table = code.table();
    out.extend(dna.chunks_exact(3).map(|c| {
        let (a, b, d) = (
            BASE[c[0] as usize],
            BASE[c[1] as usize],
            BASE[c[2] as usize],
        );
        if a == INVALID || b == INVALID || d == INVALID {
            b'X'
        } else {
            table[(a as usize) << 4 | (b as usize) << 2 | d as usize]
        }
    }));
}

/// Reverse complement, uppercased; non-ACGT bases become `N`.
pub fn reverse_complement(dna: &[u8], out: &mut Vec<u8>) {
    out.clear();
    out.extend(dna.iter().rev().map(|b| match b.to_ascii_uppercase() {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' | b'U' => b'A',
        _ => b'N',
    }));
}

/// Translates all six frames into `frames`, reusing its buffers and `rc` as scratch.
///
/// Frames 0–2 start at offsets 0–2 of `dna`; frames 3–5 at offsets 0–2 of its reverse
/// complement, so frame `f` of `dna` is frame `(f + 3) % 6` of its reverse complement.
pub(crate) fn six_frames_into(
    dna: &[u8],
    code: GeneticCode,
    rc: &mut Vec<u8>,
    frames: &mut [Vec<u8>; 6],
) {
    reverse_complement(dna, rc);
    for (f, frame) in frames.iter_mut().enumerate() {
        frame.clear();
        let strand = if f < 3 { dna } else { rc.as_slice() };
        if let Some(s) = strand.get(f % 3..) {
            translate(s, code, frame);
        }
    }
}

/// Translates all six frames of `dna`: frames 0–2 start at offsets 0–2 of `dna`, frames 3–5
/// at offsets 0–2 of its reverse complement.
pub fn six_frames(dna: &[u8], code: GeneticCode) -> [Vec<u8>; 6] {
    let mut frames = Default::default();
    six_frames_into(dna, code, &mut Vec::new(), &mut frames);
    frames
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;

    fn tr(dna: &[u8], code: GeneticCode) -> Vec<u8> {
        let mut out = Vec::new();
        translate(dna, code, &mut out);
        out
    }

    #[test]
    fn translates_known_codons() {
        assert_eq!(tr(b"ATGGCCTGGTAAtga", GeneticCode::Bacterial), b"MAW**");
        assert_eq!(tr(b"ATGTGA", GeneticCode::Mycoplasma), b"MW");
        assert_eq!(tr(b"ATGNCCAUGGC", GeneticCode::Bacterial), b"MXM");
    }

    proptest! {
        #[test]
        fn frames_are_symmetric_under_reverse_complement(dna in "[ACGTN]{0,200}") {
            let fwd = six_frames(dna.as_bytes(), GeneticCode::Bacterial);
            let mut rc = Vec::new();
            reverse_complement(dna.as_bytes(), &mut rc);
            let rev = six_frames(&rc, GeneticCode::Bacterial);
            for f in 0..6 {
                prop_assert_eq!(&fwd[f], &rev[(f + 3) % 6]);
            }
        }
    }
}
