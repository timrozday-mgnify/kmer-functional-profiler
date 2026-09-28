//! Pure-Rust kernels for the k-mer functional profiler. No PyO3 here, so the
//! final CLI reuses them unchanged.

/// Crate version, re-exported to Python as `_core.__version__`.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");
