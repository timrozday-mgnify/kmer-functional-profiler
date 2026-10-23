//! DIAMOND tabular alignments (`kfp-pileup call`'s crude pass) into per-residue columns:
//! one row per alignment and one per aligned reference column, read from the BTOP string as
//! `kfp_pileup.pileup` did in Python. Fields: `qseqid sseqid bitscore qstart qend sstart send
//! qframe btop`.

use std::collections::HashMap;

use crate::Error;

/// The reference: sequence names to ids, and their residues.
#[derive(Debug, Default)]
pub struct DiamondRefs {
    ids: HashMap<Box<[u8]>, u32>,
    residues: Vec<u8>,
    offsets: Vec<usize>,
}

impl DiamondRefs {
    /// Sequence `i` is `names[i]` with residues `seqs[i]`.
    pub fn new<N: AsRef<[u8]>, S: AsRef<[u8]>>(names: &[N], seqs: &[S]) -> Result<Self, Error> {
        if names.len() != seqs.len() {
            return Err(Error::LengthMismatch);
        }
        let mut refs = Self {
            ids: HashMap::with_capacity(names.len()),
            residues: Vec::new(),
            offsets: vec![0],
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
}

fn bad(line: &[u8], what: &str) -> Error {
    Error::DiamondLine(format!("{what}: {}", String::from_utf8_lossy(line)))
}

fn field<T: std::str::FromStr>(f: &[u8], line: &[u8]) -> Result<T, Error> {
    std::str::from_utf8(f)
        .ok()
        .and_then(|s| s.parse().ok())
        .ok_or_else(|| bad(line, "a malformed number"))
}

impl DiamondBatch {
    /// Parse whole reads from the start of `buf` until `batch_reads` reads are taken or the
    /// buffer runs out. A read is whole once a line of another read follows it, or at the end
    /// of `buf` when `eof`. Returns the bytes consumed (the next read's first line).
    pub fn parse(
        refs: &DiamondRefs,
        buf: &[u8],
        batch_reads: usize,
        eof: bool,
    ) -> Result<(Self, usize), Error> {
        let mut batch = Self::default();
        let mut last: &[u8] = &[];
        let (mut pos, mut read_start) = (0, 0);
        // the batch's columns and events as they stood when the current read began
        let mut mark = batch.marks();
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
            let qseqid = line.split(|&b| b == b'\t').next().unwrap_or_default();
            if qseqid != last && batch.names.len() >= batch_reads {
                return Ok((batch, pos)); // the previous read is whole, and the batch full
            }
            batch.line(refs, line, &mut last, &mut read_start, pos, &mut mark)?;
            pos = next;
        }
        if eof && pos >= buf.len() {
            return Ok((batch, buf.len()));
        }
        // the last read may continue past the buffer: give it back
        batch.truncate(mark);
        Ok((batch, read_start))
    }

    fn marks(&self) -> [usize; 3] {
        [self.read.len(), self.hit.len(), self.ev_hit.len()]
    }

    fn truncate(&mut self, [hits, cols, events]: [usize; 3]) {
        self.names.pop();
        for v in [
            &mut self.start,
            &mut self.end,
            &mut self.ins,
            &mut self.shifts,
        ] {
            v.truncate(hits);
        }
        self.read.truncate(hits);
        self.seq.truncate(hits);
        self.score.truncate(hits);
        self.reverse.truncate(hits);
        self.hit.truncate(cols);
        self.col.truncate(cols);
        self.qpos.truncate(cols);
        self.rres.truncate(cols);
        self.sres.truncate(cols);
        self.ev_hit.truncate(events);
        self.ev_col.truncate(events);
        self.ev_kind.truncate(events);
        self.ev_key.truncate(events);
    }

    /// One line. A new read records where it began (`read_start`, the
    /// batch's lengths in `mark`) so a read cut by the buffer's end can be given back.
    fn line<'a>(
        &mut self,
        refs: &DiamondRefs,
        line: &'a [u8],
        last: &mut &'a [u8],
        read_start: &mut usize,
        pos: usize,
        mark: &mut [usize; 3],
    ) -> Result<(), Error> {
        let f: Vec<&[u8]> = line.split(|&b| b == b'\t').collect();
        if f.len() < 9 {
            return Err(bad(line, "fewer than 9 fields"));
        }
        if f[0] != *last {
            *last = f[0];
            *read_start = pos;
            *mark = self.marks();
            self.names.push(String::from_utf8_lossy(f[0]).into_owned());
        }
        let seq = *refs
            .ids
            .get(f[1])
            .ok_or_else(|| bad(line, "a subject not in the reference"))?;
        let (qstart, sstart, frame): (i64, i32, i64) =
            (field(f[3], line)?, field(f[5], line)?, field(f[7], line)?);
        let hit = u32::try_from(self.read.len()).map_err(|_| bad(line, "too many alignments"))?;
        self.read.push((self.names.len() - 1) as u32);
        self.seq.push(seq);
        self.score.push(field(f[2], line)?);
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
        DiamondRefs::new(&["s0", "s1"], &["MKVLAAGIVLW", "MKVLSAGIILW"]).unwrap()
    }

    #[test]
    fn btop_runs_substitutions_and_indels() {
        // 4 matches, A->S, an inserted L, a deleted residue (A), 4 matches; frame +1
        let text = b"r1\ts0\t30\t1\t33\t1\t11\t1\t4SAL--A4\n";
        let (b, used) = DiamondBatch::parse(&refs(), text, 10, true).unwrap();
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
        let (b, used) = DiamondBatch::parse(&refs(), text, 10, false).unwrap();
        assert_eq!(
            (b.names, b.read.len(), b.hit.len()),
            (vec!["r1".to_owned()], 2, 22)
        );
        assert_eq!(&text[used..used + 3], b"r2\t");
        // a full batch stops before the next read
        let (b, used) = DiamondBatch::parse(&refs(), text, 1, true).unwrap();
        assert_eq!((b.names.len(), &text[used..used + 2]), (1, &b"r2"[..]));
        // at eof the last read is whole
        let (b, used) = DiamondBatch::parse(&refs(), &text[..text.len() - 5], 10, true).unwrap();
        assert_eq!((b.names.len(), used), (2, text.len() - 5));
    }

    #[test]
    fn unknown_subject_is_an_error() {
        let text = b"r1\tx\t30\t1\t33\t1\t11\t1\t11\n";
        assert!(DiamondBatch::parse(&refs(), text, 10, true).is_err());
    }
}
