//! Python bindings: batch kernels returning numpy columns, with the GIL released.

use std::path::PathBuf;
use std::sync::{Mutex, PoisonError};

use kmer_functional_profiler_core::{
    self as kfp, Column, DnaScanner, Error, Hits, KmerParams, PackedTable,
};
use numpy::{IntoPyArray, PyArray1, PyReadonlyArray1, PyReadwriteArray1};
use pyo3::exceptions::{PyOSError, PyValueError};
use pyo3::prelude::*;
use pyo3::pybacked::PyBackedBytes;
use pyo3::types::{PyAny, PyBytes, PyDict};

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

/// Distinct k-mers per run of proteins sharing a (non-decreasing) group id: columns
/// `group` and `n_kmers`, one row per run.
#[pyfunction]
#[pyo3(signature = (seqs, groups, k, *, alphabet = "protein", max_hash = u64::MAX))]
fn distinct_kmers<'py>(
    py: Python<'py>,
    seqs: Vec<PyBackedBytes>,
    groups: PyReadonlyArray1<'py, u32>,
    k: usize,
    alphabet: &str,
    max_hash: u64,
) -> PyResult<Bound<'py, PyDict>> {
    let params = params(k, alphabet, 11, "all", max_hash)?;
    let groups = groups.as_slice()?;
    let (ids, counts) = py
        .detach(|| kfp::distinct_kmers(&seqs, groups, &params))
        .map_err(to_py_err)?;
    let dict = PyDict::new(py);
    dict.set_item("group", ids.into_pyarray(py))?;
    dict.set_item("n_kmers", counts.into_pyarray(py))?;
    Ok(dict)
}

/// Adds `hashes` to the blocked Bloom filter `bits` in place (length a multiple of 64).
#[pyfunction]
fn bloom_insert(
    py: Python<'_>,
    mut bits: PyReadwriteArray1<'_, u8>,
    hashes: PyReadonlyArray1<'_, u64>,
) -> PyResult<()> {
    let (bits, hashes) = (bits.as_slice_mut()?, hashes.as_slice()?);
    py.detach(|| kfp::bloom_insert(bits, hashes))
        .map_err(to_py_err)
}

/// Whether each hash may be in the blocked Bloom filter `bits`.
#[pyfunction]
fn bloom_contains<'py>(
    py: Python<'py>,
    bits: PyReadonlyArray1<'py, u8>,
    hashes: PyReadonlyArray1<'py, u64>,
) -> PyResult<Bound<'py, PyArray1<bool>>> {
    let (bits, hashes) = (bits.as_slice()?, hashes.as_slice()?);
    let found = py
        .detach(|| kfp::bloom_contains(bits, hashes))
        .map_err(to_py_err)?;
    Ok(found.into_pyarray(py))
}

/// A borrowed 1-d unsigned numpy array of any width (memory-mapped `.npy` files included).
enum Borrowed<'py> {
    U8(PyReadonlyArray1<'py, u8>),
    U16(PyReadonlyArray1<'py, u16>),
    U32(PyReadonlyArray1<'py, u32>),
    U64(PyReadonlyArray1<'py, u64>),
}

impl<'py> Borrowed<'py> {
    fn new(array: &Bound<'py, PyAny>) -> PyResult<Self> {
        if let Ok(a) = array.extract() {
            return Ok(Self::U8(a));
        }
        if let Ok(a) = array.extract() {
            return Ok(Self::U16(a));
        }
        if let Ok(a) = array.extract() {
            return Ok(Self::U32(a));
        }
        Ok(Self::U64(array.extract().map_err(|_| {
            PyValueError::new_err("expected a 1-d unsigned integer array")
        })?))
    }

    fn column(&self) -> PyResult<Column<'_>> {
        Ok(match self {
            Self::U8(a) => Column::U8(a.as_slice()?),
            Self::U16(a) => Column::U16(a.as_slice()?),
            Self::U32(a) => Column::U32(a.as_slice()?),
            Self::U64(a) => Column::U64(a.as_slice()?),
        })
    }
}

/// The arrays of a Python `PackedTable`, borrowed in place, and its layout.
struct Packed<'py> {
    arrays: [Borrowed<'py>; 5],
    max_hash: u64,
    shift: u32,
    fp_bits: u32,
}

impl<'py> Packed<'py> {
    fn new(table: &Bound<'py, PyAny>) -> PyResult<Self> {
        let array = |name: &str| Borrowed::new(&table.getattr(name)?);
        Ok(Self {
            arrays: [
                array("offsets")?,
                array("fingerprints")?,
                array("set_ids")?,
                array("set_offsets")?,
                array("set_values")?,
            ],
            max_hash: table.getattr("max_hash")?.extract()?,
            shift: table.getattr("shift")?.extract()?,
            fp_bits: table.getattr("fp_bits")?.extract()?,
        })
    }

    fn table(&self) -> PyResult<PackedTable<'_>> {
        let [offsets, fingerprints, set_ids, set_offsets, set_values] = &self.arrays;
        Ok(PackedTable {
            max_hash: self.max_hash,
            shift: self.shift,
            fp_bits: self.fp_bits,
            offsets: offsets.column()?,
            fingerprints: fingerprints.column()?,
            set_ids: set_ids.column()?,
            set_offsets: set_offsets.column()?,
            set_values: set_values.column()?,
        })
    }
}

/// Value-set id per hash in a `PackedTable`, or -1.
#[pyfunction]
fn packed_lookup<'py>(
    py: Python<'py>,
    table: &Bound<'py, PyAny>,
    hashes: PyReadonlyArray1<'py, u64>,
) -> PyResult<Bound<'py, PyArray1<i64>>> {
    let packed = Packed::new(table)?;
    let (table, hashes) = (packed.table()?, hashes.as_slice()?);
    #[allow(clippy::cast_possible_wrap)] // set ids are far below 2^63
    let ids: Vec<i64> = py.detach(|| {
        hashes
            .iter()
            .map(|&h| table.lookup(h).map_or(-1, |s| s as i64))
            .collect()
    });
    Ok(ids.into_pyarray(py))
}

/// Unit hits of sampled hashes through a `PackedTable`: columns `unit`, `hash`, `read`,
/// `pin_q` and `holders`, one row per (hash occurrence, unit passing its `max_hash_g`).
#[pyfunction]
fn unit_hits<'py>(
    py: Python<'py>,
    table: &Bound<'py, PyAny>,
    max_hash_g: PyReadonlyArray1<'py, u64>,
    hashes: PyReadonlyArray1<'py, u64>,
    reads: PyReadonlyArray1<'py, u64>,
) -> PyResult<Bound<'py, PyDict>> {
    let packed = Packed::new(table)?;
    let table = packed.table()?;
    let (max_hash_g, hashes, reads) = (
        max_hash_g.as_slice()?,
        hashes.as_slice()?,
        reads.as_slice()?,
    );
    let hits = py
        .detach(|| kfp::unit_hits(&table, max_hash_g, hashes, reads))
        .map_err(to_py_err)?;
    let dict = PyDict::new(py);
    dict.set_item("unit", hits.unit.into_pyarray(py))?;
    dict.set_item("hash", hits.hash.into_pyarray(py))?;
    dict.set_item("read", hits.read.into_pyarray(py))?;
    dict.set_item("pin_q", hits.pin_q.into_pyarray(py))?;
    dict.set_item("holders", hits.holders.into_pyarray(py))?;
    Ok(dict)
}

/// Greedy gather over (unit, hash) pairs: the units taken, in order, and their k-mers.
#[pyfunction]
fn gather<'py>(
    py: Python<'py>,
    units: PyReadonlyArray1<'py, u32>,
    hashes: PyReadonlyArray1<'py, u64>,
    t_g: PyReadonlyArray1<'py, f64>,
) -> PyResult<Bound<'py, PyDict>> {
    let (units, hashes, t_g) = (units.as_slice()?, hashes.as_slice()?, t_g.as_slice()?);
    let got = py
        .detach(|| kfp::gather(units, hashes, t_g))
        .map_err(to_py_err)?;
    let dict = PyDict::new(py);
    dict.set_item("unit", got.unit.into_pyarray(py))?;
    dict.set_item("kmers_unique", got.kmers_unique.into_pyarray(py))?;
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
        max_hash = u64::MAX, batch_reads = 100_000, min_qual = 0,
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
        min_qual: u8,
    ) -> PyResult<Self> {
        let params = params(k, alphabet, genetic_code, frames, max_hash)?;
        let inner = kfp::FastxHits::open(&r1, r2.as_deref(), params, batch_reads, min_qual)
            .map_err(to_py_err)?;
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
    m.add_function(wrap_pyfunction!(distinct_kmers, m)?)?;
    m.add_function(wrap_pyfunction!(bloom_insert, m)?)?;
    m.add_function(wrap_pyfunction!(bloom_contains, m)?)?;
    m.add("BLOOM_BLOCK_BYTES", kfp::BLOCK_BYTES)?;
    m.add_function(wrap_pyfunction!(packed_lookup, m)?)?;
    m.add_function(wrap_pyfunction!(unit_hits, m)?)?;
    m.add_function(wrap_pyfunction!(gather, m)?)?;
    m.add_class::<PyFastxHits>()?;
    Ok(())
}
