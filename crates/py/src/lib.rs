//! Python bindings: batch kernels returning numpy columns, with the GIL released.

use std::path::PathBuf;
use std::sync::{Mutex, PoisonError};

use kmer_functional_profiler_core::{self as kfp, DnaScanner, Error, Hits, KmerParams};
use numpy::IntoPyArray;
use pyo3::exceptions::{PyOSError, PyValueError};
use pyo3::prelude::*;
use pyo3::pybacked::PyBackedBytes;
use pyo3::types::{PyBytes, PyDict};

fn to_py_err(e: Error) -> PyErr {
    match e {
        Error::Fastx(_) => PyOSError::new_err(e.to_string()),
        _ => PyValueError::new_err(e.to_string()),
    }
}

fn params(
    k: usize,
    alphabet: &str,
    genetic_code: u8,
    frames: &str,
    max_hash: u64,
) -> PyResult<KmerParams> {
    let build = || -> Result<KmerParams, Error> {
        Ok(KmerParams::new(k, alphabet.parse()?)?
            .with_genetic_code(genetic_code.try_into()?)
            .with_frames(frames.parse()?)
            .with_max_hash(max_hash))
    };
    build().map_err(to_py_err)
}

fn hits_dict(py: Python<'_>, hits: Hits) -> PyResult<Bound<'_, PyDict>> {
    let dict = PyDict::new(py);
    dict.set_item("read", hits.read.into_pyarray(py))?;
    dict.set_item("mate", hits.mate.into_pyarray(py))?;
    dict.set_item("frame", hits.frame.into_pyarray(py))?;
    dict.set_item("hash", hits.hash.into_pyarray(py))?;
    Ok(dict)
}

/// Six-frame translation of one DNA sequence (frames 0-2 forward, 3-5 reverse complement).
#[pyfunction]
#[pyo3(signature = (seq, genetic_code = 11))]
fn translate_frames<'py>(
    py: Python<'py>,
    seq: &[u8],
    genetic_code: u8,
) -> PyResult<Vec<Bound<'py, PyBytes>>> {
    let code = genetic_code.try_into().map_err(to_py_err)?;
    Ok(kfp::six_frames(seq, code)
        .iter()
        .map(|f| PyBytes::new(py, f))
        .collect())
}

/// FracMinHash threshold keeping `fraction` of hashes (keep iff hash <= max_hash).
#[pyfunction]
#[pyo3(name = "max_hash")]
fn max_hash_py(fraction: f64) -> u64 {
    kfp::max_hash(fraction)
}

/// Sampled k-mer hits of DNA reads; `read` is the index into `seqs`.
#[pyfunction]
#[pyo3(signature = (seqs, k, *, alphabet = "protein", genetic_code = 11, frames = "stopfree", max_hash = u64::MAX))]
fn hash_dna<'py>(
    py: Python<'py>,
    seqs: Vec<PyBackedBytes>,
    k: usize,
    alphabet: &str,
    genetic_code: u8,
    frames: &str,
    max_hash: u64,
) -> PyResult<Bound<'py, PyDict>> {
    let params = params(k, alphabet, genetic_code, frames, max_hash)?;
    let hits = py.detach(|| {
        let mut scanner = DnaScanner::new(params);
        let mut hits = Hits::default();
        for (read, seq) in (0u64..).zip(&seqs) {
            scanner.scan(seq, read, 0, &mut hits);
        }
        hits
    });
    hits_dict(py, hits)
}

/// Sampled k-mer hashes of protein sequences: columns `seq` (index into `seqs`) and `hash`.
#[pyfunction]
#[pyo3(signature = (seqs, k, *, alphabet = "protein", max_hash = u64::MAX))]
fn hash_proteins<'py>(
    py: Python<'py>,
    seqs: Vec<PyBackedBytes>,
    k: usize,
    alphabet: &str,
    max_hash: u64,
) -> PyResult<Bound<'py, PyDict>> {
    let params = params(k, alphabet, 11, "all", max_hash)?;
    let (seq_ids, hashes) = py.detach(|| {
        let (mut seq_ids, mut hashes) = (Vec::new(), Vec::new());
        for (i, seq) in (0u64..).zip(&seqs) {
            kfp::protein_kmers(seq, &params, |h| {
                seq_ids.push(i);
                hashes.push(h);
            });
        }
        (seq_ids, hashes)
    });
    let dict = PyDict::new(py);
    dict.set_item("seq", seq_ids.into_pyarray(py))?;
    dict.set_item("hash", hashes.into_pyarray(py))?;
    Ok(dict)
}

/// Iterator over FASTA/FASTQ (optionally paired, gzip/zstd) yielding dicts of hit columns.
#[pyclass(name = "FastxHits")]
struct PyFastxHits {
    inner: Mutex<kfp::FastxHits>,
}

#[pymethods]
impl PyFastxHits {
    #[new]
    #[pyo3(signature = (
        r1, r2 = None, *, k, alphabet = "protein", genetic_code = 11, frames = "stopfree",
        max_hash = u64::MAX, batch_reads = 100_000,
    ))]
    #[allow(clippy::too_many_arguments)]
    fn new(
        r1: PathBuf,
        r2: Option<PathBuf>,
        k: usize,
        alphabet: &str,
        genetic_code: u8,
        frames: &str,
        max_hash: u64,
        batch_reads: usize,
    ) -> PyResult<Self> {
        let params = params(k, alphabet, genetic_code, frames, max_hash)?;
        let inner =
            kfp::FastxHits::open(&r1, r2.as_deref(), params, batch_reads).map_err(to_py_err)?;
        Ok(Self {
            inner: Mutex::new(inner),
        })
    }

    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__<'py>(&self, py: Python<'py>) -> PyResult<Option<Bound<'py, PyDict>>> {
        let next = py.detach(|| {
            self.inner
                .lock()
                .unwrap_or_else(PoisonError::into_inner)
                .next()
        });
        next.map(|hits| hits_dict(py, hits.map_err(to_py_err)?))
            .transpose()
    }

    /// Reads (pairs) consumed so far.
    #[getter]
    fn n_reads(&self) -> u64 {
        self.inner
            .lock()
            .unwrap_or_else(PoisonError::into_inner)
            .n_reads()
    }
}

#[pymodule]
fn _core(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add("__version__", kfp::VERSION)?;
    m.add_function(wrap_pyfunction!(translate_frames, m)?)?;
    m.add_function(wrap_pyfunction!(max_hash_py, m)?)?;
    m.add_function(wrap_pyfunction!(hash_dna, m)?)?;
    m.add_function(wrap_pyfunction!(hash_proteins, m)?)?;
    m.add_class::<PyFastxHits>()?;
    Ok(())
}
