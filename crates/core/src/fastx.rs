use std::fmt;
use std::path::Path;
use std::sync::LazyLock;
use std::sync::mpsc::{self, Receiver};
use std::thread::{self, JoinHandle};

use needletail::parse_fastx_file;

use crate::{DnaScanner, Error, Hits, KmerParams, threads};

/// A batch's concatenated sequences, the end of each, its reads (pairs) and base counts.
type BatchResult = Result<(Vec<u8>, Vec<usize>, usize, BaseStats), Error>;

/// Bases read, for converting k-mer depth into read bases and for the error thinning of
/// k-mer counts.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct BaseStats {
    /// All bases, both mates.
    pub bases: u64,
    /// Bases that are not A, C, G or T after the quality mask: no k-mer spans them.
    pub lost: u64,
    /// Sum of Phred error probabilities over the other bases (0 for FASTA).
    pub expected_errors: f64,
}

impl BaseStats {
    fn add(&mut self, other: Self) {
        self.bases += other.bases;
        self.lost += other.lost;
        self.expected_errors += other.expected_errors;
    }

    /// Counts the bases of `seq` (already masked) with their Phred+33 qualities.
    fn count(&mut self, seq: &[u8], qual: Option<&[u8]>) {
        self.bases += seq.len() as u64;
        let called = |b: &u8| matches!(b, b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't');
        match qual {
            Some(qual) => {
                for (b, &q) in seq.iter().zip(qual) {
                    if called(b) {
                        self.expected_errors += PHRED_ERROR[q as usize];
                    } else {
                        self.lost += 1;
                    }
                }
            }
            None => self.lost += seq.iter().filter(|b| !called(b)).count() as u64,
        }
    }
}

/// Error probability of each Phred+33 quality byte.
static PHRED_ERROR: LazyLock<[f64; 256]> = LazyLock::new(|| {
    std::array::from_fn(|q| 10f64.powf(-(q.saturating_sub(33) as f64) / 10.0).min(1.0))
});

/// Streams FASTA/FASTQ (plain, gzip or zstd; single or paired) as batches of k-mer hits.
///
/// Each item holds the hits of up to `batch_reads` reads (pairs). A batch can be empty
/// when none of its reads had a sampled k-mer; [`FastxHits::n_reads`] counts all reads.
/// A batch's records are read on one thread and scanned on all, in contiguous runs of
/// reads joined in order, so hits come in input order at any thread count.
pub struct FastxHits {
    rx: Receiver<BatchResult>,
    reader: Option<JoinHandle<()>>,
    scanner: DnaScanner,
    mates: u8,
    n_reads: u64,
    base_stats: BaseStats,
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
    /// Bases with Phred quality below `min_qual` become `N` (0 = off; FASTA is unaffected),
    /// so they remove only the k-mers spanning their codon and never create a stop.
    ///
    /// # Errors
    /// [`Error::Fastx`] if a file cannot be opened or its format is not recognised.
    pub fn open(
        r1: &Path,
        r2: Option<&Path>,
        params: KmerParams,
        batch_reads: usize,
        min_qual: u8,
    ) -> Result<Self, Error> {
        let mut r1_reader = parse_fastx_file(r1)?;
        let mut r2_reader = r2.map(parse_fastx_file).transpose()?;
        let batch_reads = batch_reads.max(1);
        let mates = if r2_reader.is_some() { 2 } else { 1 };

        let (tx, rx) = mpsc::sync_channel(2);

        // ponytail: mates are paired by position only; read IDs are not compared.
        // Reading and decompression run one batch ahead of the scan, on their own thread.
        let reader = thread::spawn(move || {
            loop {
                let mut seqs = Vec::new();
                let mut ends = Vec::new();
                let mut n = 0;
                let mut stats = BaseStats::default();
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
                    let start = seqs.len();
                    push_masked(&mut seqs, &rec1.seq(), rec1.qual(), min_qual);
                    stats.count(&seqs[start..], rec1.qual());
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
                        let start = seqs.len();
                        push_masked(&mut seqs, &rec2.seq(), rec2.qual(), min_qual);
                        stats.count(&seqs[start..], rec2.qual());
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
                if tx.send(Ok((seqs, ends, n, stats))).is_err() {
                    break;
                }
            }
        });

        Ok(Self {
            rx,
            reader: Some(reader),
            scanner: DnaScanner::new(params),
            mates,
            n_reads: 0,
            base_stats: BaseStats::default(),
        })
    }

    /// Reads (pairs) consumed so far.
    pub fn n_reads(&self) -> u64 {
        self.n_reads
    }

    /// Bases of the reads consumed so far.
    pub fn base_stats(&self) -> BaseStats {
        self.base_stats
    }

    /// Reads and scans the next batch.
    fn next_batch(&mut self) -> Result<Option<Hits>, Error> {
        let (seqs, ends, n, stats) = match self.rx.recv() {
            Ok(Ok(batch)) => batch,
            Ok(Err(e)) => return Err(e),
            Err(_) => {
                // Disconnected: the reader finished, or panicked, which must not pass as
                // end of input.
                if let Some(Err(panic)) = self.reader.take().map(JoinHandle::join) {
                    std::panic::resume_unwind(panic);
                }
                return Ok(None);
            }
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
        self.base_stats.add(stats);
        let mut hits = Hits::default();
        for part in parts {
            hits.append(part);
        }
        Ok(Some(hits))
    }
}

/// Appends `seq` to `out`, with bases of quality below `min_qual` (Phred+33) as `N`.
fn push_masked(out: &mut Vec<u8>, seq: &[u8], qual: Option<&[u8]>, min_qual: u8) {
    match qual {
        Some(qual) if min_qual > 0 => {
            let cut = min_qual.saturating_add(33);
            out.extend(
                seq.iter()
                    .zip(qual)
                    .map(|(&b, &q)| if q < cut { b'N' } else { b }),
            );
        }
        _ => out.extend_from_slice(seq),
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn base_stats_count_lost_bases_and_phred_errors() {
        let mut stats = BaseStats::default();
        // Q10 (+), Q20 (5), Q30 (?), and an N whose quality is ignored.
        stats.count(b"ACGN", Some(b"+5?I"));
        assert_eq!((stats.bases, stats.lost), (4, 1));
        assert!((stats.expected_errors - (0.1 + 0.01 + 0.001)).abs() < 1e-12);
        stats.count(b"acNNt", None); // FASTA: no errors known
        assert_eq!((stats.bases, stats.lost), (9, 3));
        assert!((stats.expected_errors - 0.111).abs() < 1e-12);
    }
}
