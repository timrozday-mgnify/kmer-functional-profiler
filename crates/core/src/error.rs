/// Errors from parameter validation and FASTX input.
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum Error {
    #[error("k = {k} is outside 1..={max} for this alphabet")]
    InvalidK { k: usize, max: usize },
    #[error("unknown alphabet {0:?} (expected protein, murphy10 or dayhoff)")]
    UnknownAlphabet(String),
    #[error("unsupported genetic code {0} (expected 11 or 4)")]
    UnknownGeneticCode(u8),
    #[error("unknown frame mode {0:?} (expected stopfree, all, edges or edges:M)")]
    UnknownFrameMode(String),
    #[error("group ids must be non-decreasing and one per sequence")]
    UnsortedGroups,
    #[error("Bloom filter of {0} bytes is not a non-empty multiple of 64")]
    BloomSize(usize),
    #[error("mate files have different read counts")]
    MateCountMismatch,
    #[error("hashes and reads differ in length")]
    LengthMismatch,
    #[error("set value names unit {0}, beyond the unit table")]
    UnitOutOfRange(u64),
    #[error(transparent)]
    Fastx(#[from] needletail::errors::ParseError),
}
