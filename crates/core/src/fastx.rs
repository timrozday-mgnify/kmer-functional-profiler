use std::fmt;
use std::path::Path;

use needletail::{FastxReader, parse_fastx_file};

use crate::{DnaScanner, Error, Hits, KmerParams};

/// Streams FASTA/FASTQ (plain, gzip or zstd; single or paired) as batches of k-mer hits.
///
/// Each item holds the hits of up to `batch_reads` reads (pairs). A batch can be empty
/// when none of its reads had a sampled k-mer; [`FastxHits::n_reads`] counts all reads.
pub struct FastxHits {
    r1: Box<dyn FastxReader>,
    r2: Option<Box<dyn FastxReader>>,
    scanner: DnaScanner,
    batch_reads: usize,
    n_reads: u64,
}

impl fmt::Debug for FastxHits {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("FastxHits")
            .field("paired", &self.r2.is_some())
            .field("params", self.scanner.params())
            .field("n_reads", &self.n_reads)
            .finish_non_exhaustive()
    }
}

impl FastxHits {
    /// Opens `r1` and, for paired input, its mate file `r2`.
    ///
    /// # Errors
    /// [`Error::Fastx`] if a file cannot be opened or its format is not recognised.
    pub fn open(
        r1: &Path,
        r2: Option<&Path>,
        params: KmerParams,
        batch_reads: usize,
    ) -> Result<Self, Error> {
        Ok(Self {
            r1: parse_fastx_file(r1)?,
            r2: r2.map(parse_fastx_file).transpose()?,
            scanner: DnaScanner::new(params),
            batch_reads: batch_reads.max(1),
            n_reads: 0,
        })
    }

    /// Reads (pairs) consumed so far.
    pub fn n_reads(&self) -> u64 {
        self.n_reads
    }

    /// Scans the next read (pair) into `hits`; `Ok(false)` at end of input.
    // ponytail: mates are paired by position only; read IDs are not compared.
    fn scan_next(&mut self, hits: &mut Hits) -> Result<bool, Error> {
        let Some(rec1) = self.r1.next().transpose()? else {
            if self.r2.as_mut().is_some_and(|r2| r2.next().is_some()) {
                return Err(Error::MateCountMismatch);
            }
            return Ok(false);
        };
        self.scanner.scan(&rec1.seq(), self.n_reads, 0, hits);
        if let Some(r2) = &mut self.r2 {
            let rec2 = r2.next().transpose()?.ok_or(Error::MateCountMismatch)?;
            self.scanner.scan(&rec2.seq(), self.n_reads, 1, hits);
        }
        self.n_reads += 1;
        Ok(true)
    }
}

impl Iterator for FastxHits {
    type Item = Result<Hits, Error>;

    fn next(&mut self) -> Option<Self::Item> {
        let start = self.n_reads;
        let mut hits = Hits::default();
        for _ in 0..self.batch_reads {
            match self.scan_next(&mut hits) {
                Ok(true) => {}
                Ok(false) => break,
                Err(e) => return Some(Err(e)),
            }
        }
        (self.n_reads > start).then_some(Ok(hits))
    }
}
