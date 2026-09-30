use std::fmt;
use std::path::Path;

use needletail::{FastxReader, parse_fastx_file};

use crate::{DnaScanner, Error, Hits, KmerParams, threads};

/// Streams FASTA/FASTQ (plain, gzip or zstd; single or paired) as batches of k-mer hits.
///
/// Each item holds the hits of up to `batch_reads` reads (pairs). A batch can be empty
/// when none of its reads had a sampled k-mer; [`FastxHits::n_reads`] counts all reads.
/// A batch's records are read on one thread and scanned on all, in contiguous runs of
/// reads joined in order, so hits come in input order at any thread count.
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

    /// Appends the next read (pair)'s sequences to `seqs`, ending each at `ends`;
    /// `Ok(false)` at end of input.
    // ponytail: mates are paired by position only; read IDs are not compared.
    // ponytail: reading and decompression stay on one thread (~10% of the scan's work on
    // gzip input); a reader thread overlapping the scan if that becomes the ceiling.
    fn read_next(&mut self, seqs: &mut Vec<u8>, ends: &mut Vec<usize>) -> Result<bool, Error> {
        let Some(rec1) = self.r1.next().transpose()? else {
            if self.r2.as_mut().is_some_and(|r2| r2.next().is_some()) {
                return Err(Error::MateCountMismatch);
            }
            return Ok(false);
        };
        seqs.extend_from_slice(&rec1.seq());
        ends.push(seqs.len());
        if let Some(r2) = &mut self.r2 {
            let rec2 = r2.next().transpose()?.ok_or(Error::MateCountMismatch)?;
            seqs.extend_from_slice(&rec2.seq());
            ends.push(seqs.len());
        }
        Ok(true)
    }

    /// Reads and scans the next batch.
    fn next_batch(&mut self) -> Result<Option<Hits>, Error> {
        let (mut seqs, mut ends) = (Vec::new(), Vec::new());
        let mut n = 0;
        while n < self.batch_reads && self.read_next(&mut seqs, &mut ends)? {
            n += 1;
        }
        if n == 0 {
            return Ok(None);
        }
        let mates = if self.r2.is_some() { 2 } else { 1 };
        let per_thread = n.div_ceil(threads()).max(MIN_READS_PER_THREAD);
        let (first, seqs, ends) = (self.n_reads, &seqs, &ends);
        let parts: Vec<Hits> = std::thread::scope(|scope| {
            let handles: Vec<_> = (0..n)
                .step_by(per_thread)
                .map(|lo| {
                    let mut scanner = self.scanner.clone();
                    scope.spawn(move || {
                        let mut hits = Hits::default();
                        for read in lo..(lo + per_thread).min(n) {
                            for mate in 0..mates {
                                let i = read * mates + mate;
                                let start = if i == 0 { 0 } else { ends[i - 1] };
                                #[allow(clippy::cast_possible_truncation)] // mate is 0 or 1
                                scanner.scan(
                                    &seqs[start..ends[i]],
                                    first + read as u64,
                                    mate as u8,
                                    &mut hits,
                                );
                            }
                        }
                        hits
                    })
                })
                .collect();
            handles
                .into_iter()
                .map(|h| h.join().expect("scan thread panicked"))
                .collect()
        });
        self.n_reads += n as u64;
        let mut hits = Hits::default();
        for part in parts {
            hits.append(part);
        }
        Ok(Some(hits))
    }
}

/// Reads per thread below which a batch is not split further.
const MIN_READS_PER_THREAD: usize = 1 << 12;

impl Iterator for FastxHits {
    type Item = Result<Hits, Error>;

    fn next(&mut self) -> Option<Self::Item> {
        self.next_batch().transpose()
    }
}
