use pyo3::prelude::*;

#[pymodule]
fn _core(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add("__version__", kmer_functional_profiler_core::VERSION)?;
    Ok(())
}
