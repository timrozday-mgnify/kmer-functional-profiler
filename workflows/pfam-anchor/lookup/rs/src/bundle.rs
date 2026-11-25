//! The self-contained index: `rows.bin`, a block per family, read one family at a time.
//!
//! Block (format 2): "KFAB", acc u32, columns u32, rows u32, schemes u32, members u32, then six lengths u64, then
//!   bits      columns x 20 bytes: quarter bits of a residue at a column, from the representatives' rows before collapsing
//!   heads     rows x 10: accession 5 bytes, hmm_from u16, sequence length u16, flanks (left << 4 | right)
//!   labels    schemes x rows u16: 0 is unlabelled
//!   ops       the rows' CIGARs, one after another: op << 6 | length - 1, with lengths over 62 as 63 and a u16; 0xFF ends a row
//!   seq       the rows' residues, flanks included, 5 bits each, each row from a whole byte
//!   clusters  members x 5 bytes: the last `members` rows are members of other proteins' clusters, sorted by cluster;
//!             this is each one's representative. The first stage scans representatives only; the members of the
//!             clusters it keeps join the column and banded stages (scan.rs, MEMBERS).
//! All little-endian. Accessions fit 5 bytes: 1.1e12.
use crate::cand::hash;
use crate::data::*;
use crate::scan::{bar, kadane, score_rows};
use crate::sw;
use std::fmt::Write as _;
use std::fs::File;
use std::io::{Read, Seek, SeekFrom, Write};
use std::path::Path;
use std::sync::Mutex;

const MAGIC: &[u8; 4] = b"KFAB";
const HEAD: usize = 10;
const TOP: usize = 24 + 6 * 8; // the block's fixed header

#[derive(Clone, Default, Debug, PartialEq)]
pub struct Row {
    pub acc: u64,
    pub hmm_from: u16, // 1-based
    pub left: u8,
    pub right: u8,
    pub ops: Vec<(u8, u16)>, // 0 M, 1 I, 2 D
    pub seq: Vec<u8>,
    pub labels: Vec<u16>,
    pub cluster: u64, // a member's representative; 0 for a representative
}
impl Row {
    /// The cluster the row stands for.
    pub fn unit(&self) -> u64 {
        if self.cluster == 0 { self.acc } else { self.cluster }
    }
    /// Residues on the model's columns, and each one's position in `seq` (u16::MAX: none).
    pub fn columns(&self, len: usize) -> Result<(Vec<u8>, Vec<u16>), String> {
        let (mut codes, mut pos) = (vec![NONE; len], vec![u16::MAX; len]);
        let (mut c, mut p) = (self.hmm_from as usize - 1, self.left as usize);
        for &(op, n) in &self.ops {
            let n = n as usize;
            if op != 1 && c + n > len {
                return Err(format!("{}: column {} of {len}", self.acc, c + n));
            }
            if op != 2 && p + n > self.seq.len() {
                return Err(format!("{}: residue {} of {}", self.acc, p + n, self.seq.len()));
            }
            match op {
                0 => {
                    codes[c..c + n].copy_from_slice(&self.seq[p..p + n]);
                    (0..n).for_each(|k| pos[c + k] = (p + k) as u16);
                    c += n;
                    p += n;
                }
                1 => p += n,
                _ => c += n,
            }
        }
        if p + self.right as usize != self.seq.len() {
            return Err(format!("{}: the CIGAR uses {p} residues and {} flank, of {}", self.acc, self.right, self.seq.len()));
        }
        Ok((codes, pos))
    }
}

pub fn pack5(seq: &[u8], out: &mut Vec<u8>) {
    let (mut acc, mut n) = (0u32, 0);
    for &a in seq {
        acc |= (a as u32) << n;
        n += 5;
        while n >= 8 {
            out.push(acc as u8);
            acc >>= 8;
            n -= 8;
        }
    }
    if n > 0 {
        out.push(acc as u8);
    }
}
pub fn unpack5(b: &[u8], len: usize) -> Vec<u8> {
    (0..len)
        .map(|i| {
            let (at, sh) = (i * 5 / 8, i * 5 % 8);
            let w = b[at] as u16 | (b.get(at + 1).copied().unwrap_or(0) as u16) << 8;
            ((w >> sh) & 31) as u8
        })
        .collect()
}
fn pack_ops(ops: &[(u8, u16)], out: &mut Vec<u8>) {
    for &(op, n) in ops {
        assert!(op < 3 && n > 0);
        if n <= 63 {
            out.push(op << 6 | (n - 1) as u8);
        } else {
            out.push(op << 6 | 63);
            out.extend_from_slice(&n.to_le_bytes());
        }
    }
    out.push(0xFF);
}
fn unpack_ops(b: &[u8], at: &mut usize) -> Option<Vec<(u8, u16)>> {
    let mut v = vec![];
    loop {
        let x = *b.get(*at)?;
        *at += 1;
        if x == 0xFF {
            return Some(v);
        }
        let n = x & 63;
        if n < 63 {
            v.push((x >> 6, n as u16 + 1));
        } else {
            v.push((x >> 6, u16::from_le_bytes(b.get(*at..*at + 2)?.try_into().unwrap())));
            *at += 2;
        }
    }
}

// ------------------------------------------------------------- a stream of rows, as bundle.py writes it

struct In<R: Read>(R);
impl<R: Read> In<R> {
    fn n<const N: usize>(&mut self) -> Option<[u8; N]> {
        let mut b = [0; N];
        self.0.read_exact(&mut b).ok().map(|_| b)
    }
    fn u32(&mut self) -> Option<u32> {
        self.n().map(u32::from_le_bytes)
    }
    fn u16(&mut self) -> Option<u16> {
        self.n().map(u16::from_le_bytes)
    }
    fn row(&mut self, schemes: usize) -> Option<Row> {
        let mut r = Row { acc: u64::from_le_bytes(self.n()?), hmm_from: self.u16()?, ..Row::default() };
        [r.left, r.right] = self.n()?;
        for _ in 0..self.u16()? {
            let [op] = self.n()?;
            r.ops.push((op, self.u16()?));
        }
        let mut seq = vec![0; self.u16()? as usize];
        self.0.read_exact(&mut seq).ok()?;
        r.seq = seq;
        for _ in 0..schemes {
            r.labels.push(self.u16()?);
        }
        r.cluster = u64::from_le_bytes(self.n()?);
        Some(r)
    }
}

/// Greedy collapse, longest first: a row joins a kept one with the same labels that matches at least `t` of its residues.
// ponytail: all against all within a label, so a family of millions of rows needs a clusterer in front (ANNOTATION.md)
fn collapse(cols: &[Vec<u8>], rows: &[Row], t: f32) -> Vec<usize> {
    let n: Vec<usize> = cols.iter().map(|c| c.iter().filter(|&&a| a < 20).count()).collect();
    let mut order: Vec<usize> = (0..rows.len()).collect();
    order.sort_by_key(|&i| (std::cmp::Reverse(n[i]), rows[i].acc, i));
    if t <= 0.0 {
        return (0..rows.len()).collect();
    }
    let mut kept: std::collections::HashMap<&[u16], Vec<usize>> = Default::default();
    for i in order {
        let need = (t * n[i] as f32).ceil() as usize;
        let k = kept.entry(&rows[i].labels).or_default();
        if !k.iter().any(|&j| cols[j].iter().zip(&cols[i]).filter(|(a, b)| a == b && **b < 20).count() >= need) {
            k.push(i);
        }
    }
    let mut v: Vec<usize> = kept.into_values().flatten().collect();
    v.sort_unstable();
    v
}

/// `rows`: representatives, then members sorted by cluster.
fn block(acc: u32, len: usize, schemes: usize, rows: &[Row], bits: &[u8]) -> Vec<u8> {
    let (mut heads, mut labels, mut ops, mut seq, mut clusters) = (vec![], vec![], vec![], vec![], vec![]);
    let members = rows.iter().filter(|r| r.cluster > 0).count();
    assert!(rows[..rows.len() - members].iter().all(|r| r.cluster == 0), "members must come last");
    assert!(rows[rows.len() - members..].is_sorted_by_key(|r| r.cluster), "members must be sorted by cluster");
    for r in rows {
        assert!(r.left < 16 && r.right < 16 && r.acc < 1 << 40 && r.seq.len() < 65536);
        heads.extend_from_slice(&r.acc.to_le_bytes()[..5]);
        heads.extend_from_slice(&r.hmm_from.to_le_bytes());
        heads.extend_from_slice(&(r.seq.len() as u16).to_le_bytes());
        heads.push(r.left << 4 | r.right);
        pack_ops(&r.ops, &mut ops);
        pack5(&r.seq, &mut seq);
        if r.cluster > 0 {
            assert!(r.cluster < 1 << 40);
            clusters.extend_from_slice(&r.cluster.to_le_bytes()[..5]);
        }
    }
    for s in 0..schemes {
        rows.iter().for_each(|r| labels.extend_from_slice(&r.labels[s].to_le_bytes()));
    }
    let mut o = MAGIC.to_vec();
    for x in [acc, len as u32, rows.len() as u32, schemes as u32, members as u32] {
        o.extend_from_slice(&x.to_le_bytes());
    }
    let parts: [&[u8]; 6] = [bits, &heads, &labels, &ops, &seq, &clusters];
    parts.iter().for_each(|p| o.extend_from_slice(&(p.len() as u64).to_le_bytes()));
    parts.iter().for_each(|p| o.extend_from_slice(p));
    o
}

pub struct Block {
    pub acc: u32,
    pub len: usize,
    pub reps: usize, // rows[..reps] are representatives, rows[reps..] members
    pub bits: Vec<u8>,
    pub rows: Vec<Row>,
}
pub fn read_block(b: &[u8]) -> Result<Block, String> {
    let u = |at: usize| -> u32 { u32::from_le_bytes(b[at..at + 4].try_into().unwrap()) };
    if b.len() < TOP || &b[..4] != MAGIC {
        return Err("not a block".into());
    }
    let (acc, len, n, schemes, members) = (u(4), u(8) as usize, u(12) as usize, u(16) as usize, u(20) as usize);
    let sizes: Vec<usize> = (0..6).map(|i| u64::from_le_bytes(b[24 + 8 * i..32 + 8 * i].try_into().unwrap()) as usize).collect();
    if TOP + sizes.iter().sum::<usize>() != b.len() || sizes[0] != 20 * len || sizes[1] != HEAD * n || sizes[2] != 2 * n * schemes
        || members > n || sizes[5] != 5 * members {
        return Err(format!("PF{acc:05}: the block's sizes do not add up"));
    }
    let mut at = TOP;
    let part: Vec<&[u8]> = sizes.iter().map(|&s| { at += s; &b[at - s..at] }).collect();
    let (mut rows, mut o, mut q) = (Vec::with_capacity(n), 0, 0);
    for i in 0..n {
        let h = &part[1][HEAD * i..HEAD * (i + 1)];
        let sl = u16::from_le_bytes([h[7], h[8]]) as usize;
        let bytes = part[4].get(q..q + (5 * sl).div_ceil(8)).ok_or("short sequences")?;
        q += bytes.len();
        rows.push(Row {
            acc: h[..5].iter().rev().fold(0, |a, &x| a << 8 | x as u64),
            hmm_from: u16::from_le_bytes([h[5], h[6]]),
            left: h[9] >> 4,
            right: h[9] & 15,
            ops: unpack_ops(part[3], &mut o).ok_or("short CIGARs")?,
            seq: unpack5(bytes, sl),
            labels: (0..schemes).map(|s| { let x = 2 * (s * n + i); u16::from_le_bytes([part[2][x], part[2][x + 1]]) }).collect(),
            cluster: match i.checked_sub(n - members) {
                Some(j) => part[5][5 * j..5 * j + 5].iter().rev().fold(0, |a, &x| a << 8 | x as u64),
                None => 0,
            },
        });
    }
    if o != part[3].len() || q != part[4].len() {
        return Err(format!("PF{acc:05}: bytes left over"));
    }
    Ok(Block { acc, len, reps: n - members, bits: part[0].to_vec(), rows })
}

/// bundle-build OUT_DIR T < stream: (acc u32, columns u32, rows u32, schemes u32, rows...)*. Writes rows.bin and families.tsv
/// (accession, columns, representative rows, before collapsing, offset, bytes, member rows, members left out).
pub fn build(out: &Path, t: f32) {
    let mut f = std::io::BufWriter::with_capacity(1 << 20, File::create(out.join("rows.bin")).unwrap());
    let mut i = In(std::io::BufReader::with_capacity(1 << 20, std::io::stdin().lock()));
    let (mut tsv, mut at, mut last) = (String::new(), 0u64, 0u32);
    while let Some(acc) = i.u32() {
        let (len, n, schemes) = (i.u32().unwrap() as usize, i.u32().unwrap() as usize, i.u32().unwrap() as usize);
        assert!(acc > last, "families out of order");
        last = acc;
        let all: Vec<Row> = (0..n).map(|_| i.row(schemes).expect("the stream ends inside a family")).collect();
        let (rows, mut members): (Vec<Row>, Vec<Row>) = all.into_iter().partition(|r| r.cluster == 0);
        assert!(members.is_empty() || t <= 0.0, "members are only kept uncollapsed (collapse 0)");
        let cols: Vec<Vec<u8>> = rows.iter().map(|r| r.columns(len).unwrap_or_else(|e| panic!("PF{acc:05} {e}")).0).collect();
        let mut count = vec![[0u32; AA]; len];
        cols.iter().for_each(|x| x.iter().enumerate().for_each(|(c, &a)| count[c][a as usize] += 1));
        let mut bits = vec![0u8; 20 * len];
        for c in 0..len {
            let s = count[c][..20].iter().sum::<u32>() as f32;
            for a in 0..20 {
                let p = (count[c][a] as f32 + 0.05) / (s + 1.0);
                bits[20 * c + a] = (-4.0 * p.log2()).round().clamp(0.0, 255.0) as u8;
            }
        }
        let mut kept: Vec<Row> = collapse(&cols, &rows, t).into_iter().map(|k| rows[k].clone()).collect();
        // a member whose representative has no row in the family is never reached: left out
        let mut reps: Vec<u64> = kept.iter().map(|r| r.acc).collect();
        reps.sort_unstable();
        let before = members.len();
        members.retain(|r| reps.binary_search(&r.cluster).is_ok());
        members.iter().for_each(|r| { r.columns(len).unwrap_or_else(|e| panic!("PF{acc:05} {e}")); });
        members.sort_by_key(|r| r.cluster);
        let (nr, nm) = (kept.len(), members.len());
        kept.append(&mut members);
        let b = block(acc, len, schemes, &kept, &bits);
        f.write_all(&b).unwrap();
        writeln!(tsv, "{acc}\t{len}\t{nr}\t{}\t{at}\t{}\t{nm}\t{}", rows.len(), b.len(), before - nm).unwrap();
        at += b.len() as u64;
    }
    f.flush().unwrap();
    std::fs::write(out.join("families.tsv"), tsv).unwrap();
}

/// bundle-query DIR QUERIES OUT K SECOND: QUERIES is the same stream with no schemes, a row's accession being the query's id.
/// OUT: qid, score, self score, tied clusters, accession of the first, then each scheme's commonest label among the tied rows.
/// A cluster is named by its representative's accession, also when a member's row scored best.
pub fn query(dir: &Path, queries: &Path, out: &Path, k: usize, second: usize, threads: usize) {
    let index: Vec<Vec<u64>> = std::fs::read_to_string(dir.join("families.tsv")).unwrap().lines()
        .map(|l| l.split('\t').map(|x| x.parse().unwrap()).collect()).collect();
    let mut f = File::open(dir.join("rows.bin")).unwrap();
    let mut i = In(std::io::BufReader::new(File::open(queries).unwrap()));
    let bl = blosum();
    let mut o = String::new();
    let (mut t_read, mut t_decode, mut t_search) = (0.0, 0.0, 0.0);
    while let Some(acc) = i.u32() {
        let (len, n, _) = (i.u32().unwrap() as usize, i.u32().unwrap() as usize, i.u32().unwrap());
        let qs: Vec<Row> = (0..n).map(|_| i.row(0).expect("the queries end inside a family")).collect();
        let Some(x) = index.iter().find(|x| x[0] == acc as u64) else { continue }; // a family with no row in the index
        assert_eq!(x[1] as usize, len, "PF{acc:05}: the reads' model has {len} columns, the index's {}", x[1]);
        let t0 = std::time::Instant::now();
        let mut b = vec![0; x[5] as usize];
        f.seek(SeekFrom::Start(x[4])).unwrap();
        f.read_exact(&mut b).expect("rows.bin is shorter than families.tsv says");
        let t1 = std::time::Instant::now();
        let bk = read_block(&b).unwrap_or_else(|e| panic!("{e}"));
        assert_eq!((bk.acc, bk.len), (acc, len));
        drop(b);
        let nr = bk.reps; // the first stage sees representatives only
        let rc: Vec<(Vec<u8>, Vec<u16>)> = bk.rows.iter().map(|r| r.columns(len).unwrap()).collect();
        let mclusters: Vec<u64> = bk.rows[nr..].iter().map(|r| r.cluster).collect();
        let mut m = vec![NONE; nr * len];
        rc[..nr].iter().enumerate().for_each(|(r, x)| x.0.iter().enumerate().for_each(|(c, &a)| m[c * nr + r] = a));
        let t2 = std::time::Instant::now();
        let lines = Mutex::new(vec![]);
        let next = std::sync::atomic::AtomicUsize::new(0);
        parallel(threads, threads, |_| {
            let (mut s, mut buf, mut mine) = (vec![0u16; nr], vec![], vec![]);
            loop {
                let Some(q) = qs.get(next.fetch_add(1, std::sync::atomic::Ordering::Relaxed)) else { break };
                let (row, pos) = q.columns(len).unwrap_or_else(|e| panic!("PF{acc:05}, read {e}: another Pfam release?"));
                let cols: Vec<usize> = (0..len).filter(|&c| row[c] < 20).collect();
                let cw: Vec<(&[u8], u8, u16)> = cols.iter()
                    .map(|&c| (&m[c * nr..(c + 1) * nr], row[c], bk.bits[20 * c + row[c] as usize] as u16)).collect();
                score_rows(&cw, &mut s);
                let lo = bar(&s, k.min(nr));
                let mut keep: Vec<usize> = (0..nr).filter(|&r| s[r] >= lo).collect();
                if !mclusters.is_empty() {
                    let mut us: Vec<u64> = keep.iter().map(|&r| bk.rows[r].acc).collect();
                    us.sort_unstable();
                    us.dedup();
                    for u in us {
                        let (a, e) = (mclusters.partition_point(|&x| x < u), mclusters.partition_point(|&x| x <= u));
                        keep.extend(nr + a..nr + e);
                    }
                }
                let mut v: Vec<(i16, usize)> = keep.into_iter().map(|r| {
                    let x = &rc[r].0;
                    (kadane(cols.iter().map(|&c| if x[c] < 21 { bl[row[c] as usize][x[c] as usize] } else { -4 })), r)
                }).collect();
                v.sort_unstable_by(|a, b| b.cmp(a));
                let cut = v.get(second.min(v.len()).saturating_sub(1)).map_or(0, |x| x.0);
                let mut best: (i16, Vec<usize>) = (0, vec![]);
                for &(_, r) in v.iter().take_while(|x| x.0 >= cut) {
                    let (x, tp) = &rc[r];
                    let mut d: Vec<i32> = cols.iter().filter(|&&c| x[c] == row[c]).map(|&c| tp[c] as i32 - pos[c] as i32).collect();
                    d.sort_unstable();
                    let Some(&d) = d.get(d.len() / 2) else { continue };
                    let tt = &bk.rows[r].seq;
                    let a = (d - 8).clamp(0, tt.len() as i32);
                    let e = (q.seq.len() as i32 + d + 8).clamp(a, tt.len() as i32);
                    let sc = sw::score(&q.seq, &tt[a as usize..e as usize], &bl, Some((d - 8 - a, d + 8 - a)), &mut buf);
                    if sc > best.0 {
                        best = (sc, vec![r]);
                    } else if sc == best.0 && sc > 0 {
                        best.1.push(r);
                    }
                }
                let own: i32 = q.seq.iter().filter(|&&a| a < 20).map(|&a| bl[a as usize][a as usize] as i32).sum();
                let mut tied: Vec<u64> = best.1.iter().map(|&r| bk.rows[r].unit()).collect();
                tied.sort_unstable_by_key(|&u| (hash(u), u)); // the cluster named is any of the tied, always the same
                tied.dedup();
                let mut l = format!("{}\t{}\t{own}\t{}\t{}", q.acc, best.0, tied.len(), tied.first().map_or(0, |&u| u));
                for sch in 0..bk.rows.first().map_or(0, |r| r.labels.len()) {
                    let mut c: std::collections::BTreeMap<u16, usize> = Default::default();
                    best.1.iter().for_each(|&r| *c.entry(bk.rows[r].labels[sch]).or_default() += 1);
                    // commonest; a label before none, then the lower number, on a tie
                    let top = c.iter().max_by_key(|(&x, &n)| (n, x > 0, std::cmp::Reverse(x))).map_or(0, |x| *x.0);
                    write!(l, "\t{top}").unwrap();
                }
                mine.push((q.acc, l));
            }
            lines.lock().unwrap().extend(mine);
        });
        let mut l = lines.into_inner().unwrap();
        l.sort();
        l.iter().for_each(|x| { o.push_str(&x.1); o.push('\n') });
        t_read += (t1 - t0).as_secs_f64();
        t_decode += (t2 - t1).as_secs_f64();
        t_search += t2.elapsed().as_secs_f64();
    }
    std::fs::write(out, o).unwrap();
    eprintln!("read {t_read:.2} s, decode {t_decode:.2} s, search {t_search:.2} s on {threads} threads");
}

#[cfg(test)]
mod tests {
    use super::*;
    fn rows() -> Vec<Row> {
        (0..50u64).map(|i| {
            let long = 60 + (i * 7 % 90) as u16;
            let ops = vec![(0, long), (1, 3), (0, 5), (2, 64), (0, 63), (0, 62)];
            let (left, right) = ((i % 11) as u8, (i % 4) as u8);
            let n = left as usize + long as usize + 3 + 5 + 63 + 62 + right as usize;
            // the last 10 are members of rows 1-10's clusters
            let cluster = if i >= 40 { (i - 39) * 21_000_000_003 } else { 0 };
            Row { acc: i * 21_000_000_003, hmm_from: 1 + (i % 5) as u16, left, right, ops, labels: vec![i as u16, 65535 - i as u16], cluster,
                  seq: (0..n as u64).map(|k| (hash(i * 1000 + k) % 21) as u8).collect() }
        }).collect()
    }
    #[test]
    fn round_trip() {
        let (r, bits) = (rows(), vec![7u8; 20 * 400]);
        let mut b = block(44, 400, 2, &r, &bits);
        let x = read_block(&b).unwrap();
        assert_eq!((x.acc, x.len, &x.bits, &x.rows), (44, 400, &bits, &r));
        assert!(read_block(&b[..b.len() - 1]).is_err(), "a cut block must not pass");
        b[16] = 1;
        assert!(read_block(&b).is_err(), "nor a wrong number of schemes");
        assert_eq!(x.reps, 40);
        assert_eq!(r[45].unit(), r[6].acc);
        assert_eq!(r[5].unit(), r[5].acc);
    }
    #[test]
    #[should_panic(expected = "members must come last")]
    fn members_last() {
        let mut r = rows();
        r.swap(0, 45);
        block(44, 400, 2, &r, &[7u8; 20 * 400]);
    }
    #[test]
    fn columns_are_checked() {
        let r = rows().remove(3);
        assert!(r.columns(400).is_ok());
        assert!(r.columns(200).is_err(), "past the model");
        assert!(Row { seq: r.seq[1..].to_vec(), ..r.clone() }.columns(400).is_err(), "a residue short");
    }
    #[test]
    fn collapse_keeps_labels_apart() {
        let r = rows().remove(1);
        let other = Row { labels: vec![9, 9], ..r.clone() };
        let v = vec![r.clone(), r, other];
        let c: Vec<Vec<u8>> = v.iter().map(|x| x.columns(400).unwrap().0).collect();
        assert_eq!(collapse(&c, &v, 0.7).len(), 2);
        assert_eq!(collapse(&c, &v, 0.0).len(), 3);
    }
}
