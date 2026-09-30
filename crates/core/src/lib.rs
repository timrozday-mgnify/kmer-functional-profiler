//! Pure-Rust kernels for the k-mer functional profiler: FASTX streaming,
//! translation, amino-acid k-mer packing, hashing and FracMinHash sampling.
//! No PyO3 here, so the final CLI reuses them unchanged.

mod alphabet;
mod bloom;
mod error;
mod fastx;
mod gather;
mod kmers;
mod packed;
mod translate;

pub use alphabet::Alphabet;
pub use bloom::{BLOCK_BYTES, bloom_contains, bloom_insert};
pub use error::Error;
pub use fastx::FastxHits;
pub use gather::{Gathered, gather};
pub use kmers::{DnaScanner, Hits, KmerParams, distinct_kmers, hash_kmer, max_hash, protein_kmers};
pub use packed::{Column, PIN_BITS, PackedTable, UnitHits, unit_hits};
pub use translate::{FrameMode, GeneticCode, reverse_complement, six_frames, translate};

/// Worker threads for the parallel kernels: the CPUs this process may use.
pub(crate) fn threads() -> usize {
    std::thread::available_parallelism().map_or(1, |n| n.get())
}

/// Crate version, re-exported to Python as `_core.__version__`.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");
