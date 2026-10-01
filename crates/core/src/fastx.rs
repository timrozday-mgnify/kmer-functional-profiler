use std::fmt;
use std::path::Path;
use std::sync::mpsc::{self, Receiver};
use std::thread::{self, JoinHandle};

use needletail::parse_fastx_file;

use crate::{DnaScanner, Error, Hits, KmerParams, threads};

/// Streams FASTA/FASTQ (plain, gzip or zstd; single or paired) as batches of k-mer hits.
///
/// Each item holds the hits of up to `batch_reads` reads (pairs). A batch can be empty
/// when none of its reads had a sampled k-mer; [`FastxHits::n_reads`] counts all reads.
/// A batch's records are read on one thread and scanned on all, in contiguous runs of
/// reads joined in order, so hits come in input order at any thread count.
pub type BatchResult = Result<(Vec<u8>, Vec<usize>, usize), Error>;

pub struct FastxHits {
    rx: Receiver<BatchResult>,
    _reader_thread: Option<JoinHandle<()>>,
    scanner: DnaScanner,
    mates: u8,
    n_reads: u64,
}

impl fmt::Debug for FastxHits {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("FastxHits")
            .field("paired", &(self.mates == 2))
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
        let mut r1_reader = parse_fastx_file(r1)?;
        let mut r2_reader = r2.map(parse_fastx_file).transpose()?;
        let batch_reads = batch_reads.max(1);
        let mates = if r2_reader.is_some() { 2 } else { 1 };

        let (tx, rx) = mpsc::sync_channel(2);

        let reader_thread = thread::spawn(move || {
            loop {
                let mut seqs = Vec::new();
                let mut ends = Vec::new();
                let mut n = 0;
                let mut error = None;

                while n < batch_reads {
                    let rec1 = match r1_reader.next() {
                        Some(Ok(rec)) => rec,
                        Some(Err(e)) => {
                            error = Some(Error::from(e));
                            break;
                        }
                        None => {
                            if r2_reader.as_mut().is_some_and(|r2| r2.next().is_some()) {
                                error = Some(Error::MateCountMismatch);
                            }
                            break;
                        }
                    };
                    seqs.extend_from_slice(&rec1.seq());
                    ends.push(seqs.len());

                    if let Some(r2) = &mut r2_reader {
                        let rec2 = match r2.next() {
                            Some(Ok(rec)) => rec,
                            Some(Err(e)) => {
                                error = Some(Error::from(e));
                                break;
                            }
                            None => {
                                error = Some(Error::MateCountMismatch);
                                break;
                            }
                        };
                        seqs.extend_from_slice(&rec2.seq());
                        ends.push(seqs.len());
                    }
                    n += 1;
                }

                if let Some(err) = error {
                    let _ = tx.send(Err(err));
                    break;
                }
                if n == 0 {
                    break;
                }
                if tx.send(Ok((seqs, ends, n))).is_err() {
                    break;
                }
            }
        });

        Ok(Self {
            rx,
            _reader_thread: Some(reader_thread),
            scanner: DnaScanner::new(params),
            mates,
            n_reads: 0,
        })
    }

    /// Reads (pairs) consumed so far.
    pub fn n_reads(&self) -> u64 {
        self.n_reads
    }

    /// Reads and scans the next batch.
    fn next_batch(&mut self) -> Result<Option<Hits>, Error> {
        let (seqs, ends, n) = match self.rx.recv() {
            Ok(Ok(batch)) => batch,
            Ok(Err(e)) => return Err(e),
            Err(_) => return Ok(None),
        };

        let mates = self.mates as usize;
        let per_thread = n.div_ceil(threads()).max(MIN_READS_PER_THREAD);
        let (first, seqs_ref, ends_ref) = (self.n_reads, &seqs, &ends);
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
                                let start = if i == 0 { 0 } else { ends_ref[i - 1] };
                                #[allow(clippy::cast_possible_truncation)] // mate is 0 or 1
                                scanner.scan(
                                    &seqs_ref[start..ends_ref[i]],
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
