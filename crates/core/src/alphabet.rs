use std::str::FromStr;

use crate::Error;

/// Marks bytes that are not amino acids of the alphabet (stops, `X`, ambiguity codes).
const INVALID: u8 = u8::MAX;

/// Alphabet used to encode k-mers: amino acids (reduced alphabets merge residues into
/// classes), or nucleotides for indexes of DNA targets.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
#[non_exhaustive]
pub enum Alphabet {
    /// The 20 standard amino acids, 5 bits per residue.
    #[default]
    Protein,
    /// Murphy et al. (2000) 10-letter alphabet, 4 bits per residue.
    Murphy10,
    /// Dayhoff 6-letter alphabet, 3 bits per residue.
    Dayhoff,
    /// Nucleotides A, C, G, T (U read as T), 2 bits per base. Sequences are hashed as
    /// canonical k-mers without translation (see [`crate::dna_kmers`]).
    Dna,
}

const PROTEIN: [u8; 256] = classes(&[
    "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W",
    "Y",
]);
const MURPHY10: [u8; 256] = classes(&["LVIM", "C", "A", "G", "ST", "P", "FYW", "EDNQ", "KR", "H"]);
const DAYHOFF: [u8; 256] = classes(&["AGPST", "C", "DENQ", "HKR", "ILMV", "FWY"]);
// Complement is 3 - code (A <-> T, C <-> G), which canonical k-mers rely on.
const DNA: [u8; 256] = classes(&["A", "C", "G", "TU"]);

/// Builds an ASCII lookup table mapping each residue (either case) to its class index.
const fn classes(groups: &[&str]) -> [u8; 256] {
    let mut table = [INVALID; 256];
    let mut g = 0;
    while g < groups.len() {
        let residues = groups[g].as_bytes();
        let mut i = 0;
        while i < residues.len() {
            table[residues[i] as usize] = g as u8;
            table[residues[i].to_ascii_lowercase() as usize] = g as u8;
            i += 1;
        }
        g += 1;
    }
    table
}

impl Alphabet {
    /// Bits used per residue when packing k-mers into a `u64`.
    pub fn bits(self) -> u32 {
        match self {
            Self::Protein => 5,
            Self::Murphy10 => 4,
            Self::Dayhoff => 3,
            Self::Dna => 2,
        }
    }

    /// Largest k whose packed k-mer fits in a `u64`.
    pub fn max_k(self) -> usize {
        (u64::BITS / self.bits()) as usize
    }

    /// Whether sequences are nucleotides, hashed untranslated, rather than amino acids.
    pub fn is_nucleotide(self) -> bool {
        self == Self::Dna
    }

    /// Class index of a residue byte, or `None` for stops, `X`, `N` and ambiguity codes.
    #[inline]
    pub fn encode(self, aa: u8) -> Option<u8> {
        let table = match self {
            Self::Protein => &PROTEIN,
            Self::Murphy10 => &MURPHY10,
            Self::Dayhoff => &DAYHOFF,
            Self::Dna => &DNA,
        };
        let code = table[aa as usize];
        (code != INVALID).then_some(code)
    }
}

impl FromStr for Alphabet {
    type Err = Error;

    fn from_str(s: &str) -> Result<Self, Error> {
        match s.to_ascii_lowercase().as_str() {
            "protein" => Ok(Self::Protein),
            "murphy10" => Ok(Self::Murphy10),
            "dayhoff" => Ok(Self::Dayhoff),
            "dna" => Ok(Self::Dna),
            _ => Err(Error::UnknownAlphabet(s.to_owned())),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn encodes_every_standard_residue_in_range() {
        for alphabet in [Alphabet::Protein, Alphabet::Murphy10, Alphabet::Dayhoff] {
            let n = 1u8 << alphabet.bits();
            for &aa in b"ACDEFGHIKLMNPQRSTVWYacdefghiklmnpqrstvwy" {
                assert!(
                    alphabet.encode(aa).is_some_and(|c| c < n),
                    "{alphabet:?} {aa}"
                );
            }
            for &aa in b"*XBZJUO-" {
                assert_eq!(alphabet.encode(aa), None);
            }
        }
    }

    #[test]
    fn dna_codes_complement_to_three() {
        let code = |b| Alphabet::Dna.encode(b).unwrap();
        assert_eq!(code(b'A') + code(b'T'), 3);
        assert_eq!(code(b'c') + code(b'G'), 3);
        assert_eq!(code(b'U'), code(b'T'));
        assert_eq!(Alphabet::Dna.encode(b'N'), None);
        assert_eq!(Alphabet::Dna.max_k(), 32);
        assert!("dna".parse::<Alphabet>().unwrap().is_nucleotide());
    }

    #[test]
    fn reduced_alphabets_merge_classes() {
        assert_eq!(
            Alphabet::Murphy10.encode(b'L'),
            Alphabet::Murphy10.encode(b'I')
        );
        assert_eq!(
            Alphabet::Dayhoff.encode(b'A'),
            Alphabet::Dayhoff.encode(b'T')
        );
        assert_ne!(
            Alphabet::Protein.encode(b'L'),
            Alphabet::Protein.encode(b'I')
        );
    }
}
