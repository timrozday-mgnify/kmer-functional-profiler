//! Pure-Rust kernels for the k-mer functional profiler: FASTX streaming,
//! translation, amino-acid and canonical nucleotide k-mer packing, hashing and FracMinHash
//! sampling.
//! No PyO3 here, so the final CLI reuses them unchanged.

mod align;
mod alphabet;
mod bloom;
mod components;
mod diamond;
mod distinct;
mod error;
mod fastx;
mod gather;
mod kmers;
mod packed;
mod route;
mod translate;
mod tree;

pub use align::{identities, identity};
pub use alphabet::Alphabet;
pub use bloom::{BLOCK_BYTES, bloom_contains, bloom_insert};
pub use components::components;
pub use diamond::{DiamondBatch, DiamondRefs, Dropped, HitFilter};
pub use distinct::Distinct;
pub use error::Error;
pub use fastx::{BaseStats, FastxHits};
pub use gather::{Gathered, gather};
pub use kmers::{
    DnaScanner, Hits, KmerParams, distinct_kmers, dna_kmers, hash_kmer, max_hash, protein_kmers,
    sequence_kmers,
};
pub use packed::{Column, PIN_BITS, PackedTable, UnitHits, unit_hits};
pub use route::{RouteHits, route_hits};
pub use translate::{
    EDGES_MIN_LEN, FrameMode, GeneticCode, reverse_complement, six_frames, translate,
};
pub use tree::{
    Kruskal, Layout, Pairs, Schemes, kruskal, layout, overlap_pairs, pair_weights, schemes,
};

/// Worker threads for the parallel kernels: the CPUs this process may use.
pub(crate) fn threads() -> usize {
    std::thread::available_parallelism().map_or(1, |n| n.get())
}

/// Crate version, re-exported to Python as `_core.__version__`.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");
