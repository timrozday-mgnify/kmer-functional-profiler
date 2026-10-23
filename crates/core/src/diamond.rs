//! DIAMOND tabular alignments (`kfp-pileup call`'s crude pass) into per-residue columns:
//! one row per alignment and one per aligned reference column, read from the BTOP string as
//! `kfp_pileup.pileup` did in Python. Fields: `qseqid sseqid bitscore qstart qend sstart send
//! qframe btop`.
//!
//! A read's alignments are kept or dropped together, before any is expanded ([`HitFilter`]):
//! routed DIAMOND runs apply `--top` within each block of the reference, so a read routed to
//! two blocks brings the runners-up of both, and only the read's lines together say which are
//! within `--top` of its best.

use std::collections::HashMap;

use crate::Error;

/// The reference: sequence names to ids, their residues and units.
#[derive(Debug, Default)]
pub struct DiamondRefs {
    ids: HashMap<Box<[u8]>, u32>,
    residues: Vec<u8>,
    offsets: Vec<usize>,
    units: Vec<u32>,
}

impl DiamondRefs {
    /// Sequence `i` is `names[i]` with residues `seqs[i]`, of unit `units[i]`.
    pub fn new<N: AsRef<[u8]>, S: AsRef<[u8]>>(
        names: &[N],
        seqs: &[S],
        units: &[u32],
    ) -> Result<Self, Error> {
        if names.len() != seqs.len() || names.len() != units.len() {
            return Err(Error::LengthMismatch);
        }
        let mut refs = Self {
            ids: HashMap::with_capacity(names.len()),
            residues: Vec::new(),
            offsets: vec![0],
            units: units.to_vec(),
        };
        for (i, (name, seq)) in names.iter().zip(seqs).enumerate() {
            let id = u32::try_from(i).map_err(|_| Error::LengthMismatch)?;
            refs.ids.insert(name.as_ref().into(), id);
            refs.residues.extend_from_slice(seq.as_ref());
            refs.offsets.push(refs.residues.len());
        }
        Ok(refs)
    }

    fn seq(&self, id: u32) -> &[u8] {
        let i = id as usize;
        &self.residues[self.offsets[i]..self.offsets[i + 1]]
    }
}

/// Which of a read's alignments are kept, against its best (highest bit score). A line is
/// dropped by the first rule it fails, in this order; the defaults keep everything.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HitFilter {
    /// Percent below the read's best score, over all its lines (DIAMOND's `--top`, on bit
    /// scores).
    pub top: f64,
    /// Bits below the read's best score.
    pub max_bits_below: f64,
    /// Bits below the read's best alignment to the same unit (other members of a unit).
    pub unit_bits_below: f64,
    /// Least identity of the read's best alignment (matched residues over substitutions,
    /// insertions and deletions); a read below it is dropped whole.
    pub min_identity: f64,
}

impl Default for HitFilter {
    fn default() -> Self {
        Self {
            top: f64::INFINITY,
            max_bits_below: f64::INFINITY,
            unit_bits_below: f64::INFINITY,
            min_identity: 0.0,
        }
    }
}

/// Lines dropped by each rule of [`HitFilter`], and reads (and their lines) by the floor.
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub struct Dropped {
    pub top: u64,
    pub margin: u64,
    pub unit: u64,
    pub floor_reads: u64,
    pub floor_lines: u64,
}

/// A batch of alignments, as columns. `read` indexes `names`; `hit` indexes the alignment
/// columns. Events have no reference column: `/` a base missing from the read, `\` an extra
/// base, `I` an inserted residue; `ev_col` the column they precede, `ev_key` the read
/// position key (-1 for `/`).
#[derive(Debug, Default, PartialEq)]
pub struct DiamondBatch {
    pub names: Vec<String>,
    pub read: Vec<u32>,
    pub seq: Vec<u32>,
    pub score: Vec<f64>,
    pub start: Vec<i32>,
    pub end: Vec<i32>,
    pub reverse: Vec<bool>,
    pub ins: Vec<i32>,
    pub shifts: Vec<i32>,
    pub hit: Vec<u32>,
    pub col: Vec<i32>,
    pub qpos: Vec<i64>,
    pub rres: Vec<u8>,
    pub sres: Vec<u8>,
    pub ev_hit: Vec<u32>,
    pub ev_col: Vec<i32>,
    pub ev_kind: Vec<u8>,
    pub ev_key: Vec<i64>,
    pub dropped: Dropped,
}

fn bad(line: &[u8], what: &str) -> Error {
    Error::DiamondLine(format!("{what}: {}", String::from_utf8_lossy(line)))
}

fn fields(line: &[u8]) -> Result<Vec<&[u8]>, Error> {
    let f: Vec<&[u8]> = line.split(|&b| b == b'\t').collect();
    if f.len() < 9 {
        return Err(bad(line, "fewer than 9 fields"));
    }
    Ok(f)
}

/// Identity of an alignment from its BTOP: matched residues over matches, substitutions,
/// insertions and deletions (frameshifts not counted).
fn identity(btop: &[u8]) -> f64 {
    let (mut matched, mut other, mut i) = (0u64, 0u64, 0);
    while i < btop.len() {
        if btop[i].is_ascii_digit() {
            let j = i + btop[i..].iter().take_while(|b| b.is_ascii_digit()).count();
            matched += std::str::from_utf8(&btop[i..j]).map_or(0, |n| n.parse().unwrap_or(0));
            i = j;
        } else {
            other += u64::from(btop[i] != b'/' && btop[i] != b'\\');
            i += 2;
        }
    }
    if matched + other == 0 {
        return 0.0;
    }
    matched as f64 / (matched + other) as f64
}

fn field<T: std::str::FromStr>(f: &[u8], line: &[u8]) -> Result<T, Error> {
    std::str::from_utf8(f)
        .ok()
        .and_then(|s| s.parse().ok())
        .ok_or_else(|| bad(line, "a malformed number"))
}

impl DiamondBatch {
    /// Parse whole reads from the start of `buf` until `batch_reads` reads are taken or the
    /// buffer runs out, keeping each read's lines that pass `filter`. A read is whole once a
    /// line of another read follows it, or at the end of `buf` when `eof`. Returns the bytes
    /// consumed (the first line of the next read, or of the last read if it may go on).
    pub fn parse(
        refs: &DiamondRefs,
        filter: &HitFilter,
        buf: &[u8],
        batch_reads: usize,
        eof: bool,
    ) -> Result<(Self, usize), Error> {
        let mut batch = Self::default();
        let mut lines: Vec<&[u8]> = Vec::new();
        let (mut pos, mut read_start, mut qseqid) = (0, 0, &b""[..]);
        while pos < buf.len() {
            let (line, next) = match buf[pos..].iter().position(|&b| b == b'\n') {
                Some(len) => (&buf[pos..pos + len], pos + len + 1),
                None if eof => (&buf[pos..], buf.len()),
                None => break, // a partial line: the caller reads more
            };
            if line.is_empty() {
                pos = next;
                continue;
            }
            let q = line.split(|&b| b == b'\t').next().unwrap_or_default();
            if q != qseqid {
                if !lines.is_empty() {
                    batch.read(refs, filter, &lines)?;
                    lines.clear();
                    if batch.names.len() >= batch_reads {
                        return Ok((batch, pos)); // the previous read is whole, the batch full
                    }
                }
                (qseqid, read_start) = (q, pos);
            }
            lines.push(line);
            pos = next;
        }
        if eof {
            if !lines.is_empty() {
                batch.read(refs, filter, &lines)?;
            }
            return Ok((batch, buf.len()));
        }
        Ok((batch, read_start)) // the last read may continue past the buffer: give it back
    }

    /// One read's lines: those `filter` keeps, expanded.
    fn read(
        &mut self,
        refs: &DiamondRefs,
        filter: &HitFilter,
        lines: &[&[u8]],
    ) -> Result<(), Error> {
        let mut parsed = Vec::with_capacity(lines.len());
        for &line in lines {
            let f = fields(line)?;
            let seq = *refs
                .ids
                .get(f[1])
                .ok_or_else(|| bad(line, "a subject not in the reference"))?;
            let score: f64 = field(f[2], line)?;
            parsed.push((f, seq, refs.units[seq as usize], score));
        }
        let best = parsed.iter().map(|p| p.3).fold(f64::NEG_INFINITY, f64::max);
        let best_identity = parsed
            .iter()
            .filter(|p| p.3 == best)
            .map(|p| identity(p.0[8]))
            .fold(0.0, f64::max);
        if best_identity < filter.min_identity {
            self.dropped.floor_reads += 1;
            self.dropped.floor_lines += lines.len() as u64;
            return Ok(());
        }
        let mut unit_best: HashMap<u32, f64> = HashMap::new();
        for p in &parsed {
            let b = unit_best.entry(p.2).or_insert(p.3);
            *b = b.max(p.3);
        }
        let top_floor = if filter.top.is_finite() {
            best * (1.0 - filter.top / 100.0)
        } else {
            f64::NEG_INFINITY
        };
        self.names
            .push(String::from_utf8_lossy(parsed[0].0[0]).into_owned());
        let read = (self.names.len() - 1) as u32;
        for (line, (f, seq, unit, score)) in lines.iter().zip(&parsed) {
            if *score < top_floor {
                self.dropped.top += 1;
            } else if *score < best - filter.max_bits_below {
                self.dropped.margin += 1;
            } else if *score < unit_best[unit] - filter.unit_bits_below {
                self.dropped.unit += 1;
            } else {
                self.alignment(refs, line, f, read, *seq, *score)?;
            }
        }
        Ok(())
    }

    /// One alignment's row, columns and events.
    fn alignment(
        &mut self,
        refs: &DiamondRefs,
        line: &[u8],
        f: &[&[u8]],
        read: u32,
        seq: u32,
        score: f64,
    ) -> Result<(), Error> {
        let (qstart, sstart, frame): (i64, i32, i64) =
            (field(f[3], line)?, field(f[5], line)?, field(f[7], line)?);
        let hit = u32::try_from(self.read.len()).map_err(|_| bad(line, "too many alignments"))?;
        self.read.push(read);
        self.seq.push(seq);
        self.score.push(score);
        self.start.push(sstart - 1);
        self.reverse.push(frame < 0);
        // Read position key of a base: twice its position on the original read, plus 1 on
        // the reverse strand (qstart is then the higher coordinate).
        let (rev, sign) = (i64::from(frame < 0), if frame < 0 { -1 } else { 1 });
        let key = |nt: i64| 2 * (qstart - 1 + sign * nt) + rev;
        let reference = refs.seq(seq);
        let (mut col, mut nt, mut ins, mut shifts) = (sstart - 1, 0i64, 0, 0);
        let btop = f[8];
        let mut i = 0;
        while i < btop.len() {
            if btop[i].is_ascii_digit() {
                let j = i + btop[i..].iter().take_while(|b| b.is_ascii_digit()).count();
                let n: i32 = field(&btop[i..j], line)?;
                let (a, b) = (col as usize, (col + n) as usize);
                let run = reference
                    .get(a..b)
                    .ok_or_else(|| bad(line, "an alignment past the reference's end"))?;
                for m in 0..n {
                    self.hit.push(hit);
                    self.col.push(col + m);
                    self.qpos.push(key(nt) + 6 * sign * i64::from(m));
                }
                self.rres.extend_from_slice(run);
                self.sres.extend_from_slice(run);
                (col, nt) = (col + n, nt + 3 * i64::from(n));
                i = j;
                continue;
            }
            let (q, s) = (
                btop[i],
                *btop.get(i + 1).ok_or_else(|| bad(line, "odd BTOP"))?,
            );
            i += 2;
            if q == b'/' {
                self.event(hit, col, b'/', -1);
                shifts += 1;
                nt -= 1;
            } else if q == b'\\' {
                self.event(hit, col, b'\\', key(nt));
                shifts += 1;
                nt += 1;
            } else if s == b'-' {
                self.event(hit, col, b'I', key(nt));
                ins += 1;
                nt += 3;
            } else if q == b'-' {
                self.column(hit, col, -1, b'-', s);
                col += 1;
            } else {
                self.column(hit, col, key(nt), q, s);
                col += 1;
                nt += 3;
            }
        }
        self.end.push(col);
        self.ins.push(ins);
        self.shifts.push(shifts);
        Ok(())
    }

    fn column(&mut self, hit: u32, col: i32, qpos: i64, r: u8, s: u8) {
        self.hit.push(hit);
        self.col.push(col);
        self.qpos.push(qpos);
        self.rres.push(r);
        self.sres.push(s);
    }

    fn event(&mut self, hit: u32, col: i32, kind: u8, key: i64) {
        self.ev_hit.push(hit);
        self.ev_col.push(col);
        self.ev_kind.push(kind);
        self.ev_key.push(key);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn refs() -> DiamondRefs {
        DiamondRefs::new(&["s0", "s1"], &["MKVLAAGIVLW", "MKVLSAGIILW"], &[0, 1]).unwrap()
    }

    #[test]
    fn btop_runs_substitutions_and_indels() {
        // 4 matches, A->S, an inserted L, a deleted residue (A), 4 matches; frame +1
        let text = b"r1\ts0\t30\t1\t33\t1\t11\t1\t4SAL--A4\n";
        let (b, used) =
            DiamondBatch::parse(&refs(), &HitFilter::default(), text, 10, true).unwrap();
        assert_eq!(used, text.len());
        assert_eq!(b.names, ["r1"]);
        assert_eq!((b.start[0], b.end[0], b.ins[0]), (0, 10, 1));
        assert_eq!(b.col, [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]);
        assert_eq!(b.rres, b"MKVLS-GIVL");
        assert_eq!(b.sres, b"MKVLAAGIVL");
        assert_eq!(b.qpos, [0, 6, 12, 18, 24, -1, 36, 42, 48, 54]);
        assert_eq!(
            (b.ev_col, b.ev_kind, b.ev_key),
            (vec![5], b"I".to_vec(), vec![30])
        );
    }

    #[test]
    fn whole_reads_only_and_batches_by_reads() {
        let text = b"r1\ts0\t30\t1\t33\t1\t11\t1\t11\nr1\ts1\t25\t1\t33\t1\t11\t1\t4AS6\n\
                     r2\ts0\t28\t1\t33\t1\t11\t1\t11\nr2\ts1";
        // not at eof: r2's second line is cut, so r2 is given back
        let (b, used) =
            DiamondBatch::parse(&refs(), &HitFilter::default(), text, 10, false).unwrap();
        assert_eq!(
            (b.names, b.read.len(), b.hit.len()),
            (vec!["r1".to_owned()], 2, 22)
        );
        assert_eq!(&text[used..used + 3], b"r2\t");
        // a full batch stops before the next read
        let (b, used) = DiamondBatch::parse(&refs(), &HitFilter::default(), text, 1, true).unwrap();
        assert_eq!((b.names.len(), &text[used..used + 2]), (1, &b"r2"[..]));
        // at eof the last read is whole
        let (b, used) = DiamondBatch::parse(
            &refs(),
            &HitFilter::default(),
            &text[..text.len() - 5],
            10,
            true,
        )
        .unwrap();
        assert_eq!((b.names.len(), used), (2, text.len() - 5));
    }

    fn kept(filter: &HitFilter, text: &[u8]) -> (Vec<u32>, Vec<f64>, Dropped) {
        // s0, s1 of unit 0 (members), s2 of unit 1
        let refs = DiamondRefs::new(
            &["s0", "s1", "s2"],
            &["MKVLAAGIVLW", "MKVLSAGIILW", "MKVLSAGIVLW"],
            &[0, 0, 1],
        )
        .unwrap();
        let (b, _) = DiamondBatch::parse(&refs, filter, text, 10, true).unwrap();
        (b.seq, b.score, b.dropped)
    }

    #[test]
    fn top_reapplied_over_a_reads_lines_from_two_blocks() {
        // block A's best 100, block B's best 80 with its runner-up 75: B's lines are each
        // within 10% of B's best, but not of the read's
        let text = b"r1\ts0\t100\t1\t33\t1\t11\t1\t11\nr1\ts1\t92\t1\t33\t1\t11\t1\t11\n\
                     r1\ts2\t80\t1\t33\t1\t11\t1\t11\nr1\ts2\t75\t1\t33\t1\t11\t1\t11\n";
        let top = HitFilter {
            top: 10.0,
            ..HitFilter::default()
        };
        let (seq, score, dropped) = kept(&top, text);
        assert_eq!((seq, score), (vec![0, 1], vec![100.0, 92.0]));
        assert_eq!(dropped.top, 2);
        // the same read unrouted (A's lines only) keeps the same
        let (seq, _, _) = kept(&top, &text[..text.len() / 2]);
        assert_eq!(seq, [0, 1]);
    }

    #[test]
    fn margin_unit_and_identity_floor() {
        let text = b"r1\ts0\t50\t1\t33\t1\t11\t1\t11\nr1\ts1\t44\t1\t33\t1\t11\t1\t11\n\
                     r1\ts2\t41\t1\t33\t1\t11\t1\t11\n";
        let margin = HitFilter {
            max_bits_below: 8.0,
            ..HitFilter::default()
        };
        assert_eq!(kept(&margin, text).0, [0, 1]); // 41 is 9 bits below
        // s1 is unit 0's runner-up, 6 bits below s0; s2 is another unit's best
        let unit = HitFilter {
            unit_bits_below: 5.0,
            ..HitFilter::default()
        };
        let (seq, _, dropped) = kept(&unit, text);
        assert_eq!((seq, dropped.unit), (vec![0, 2], 1));
        // best alignment at 8 of 11 identical (3 substitutions): 0.727
        let low = b"r1\ts0\t50\t1\t33\t1\t11\t1\t2SA2SA2SA2\nr2\ts0\t50\t1\t33\t1\t11\t1\t11\n";
        let floor = HitFilter {
            min_identity: 0.75,
            ..HitFilter::default()
        };
        let (seq, _, dropped) = kept(&floor, low);
        assert_eq!(
            (seq.len(), dropped.floor_reads, dropped.floor_lines),
            (1, 1, 1)
        );
        let floor = HitFilter {
            min_identity: 0.7,
            ..HitFilter::default()
        };
        assert_eq!(kept(&floor, low).0.len(), 2);
    }

    #[test]
    fn unknown_subject_is_an_error() {
        let text = b"r1\tx\t30\t1\t33\t1\t11\t1\t11\n";
        assert!(DiamondBatch::parse(&refs(), &HitFilter::default(), text, 10, true).is_err());
    }
}
