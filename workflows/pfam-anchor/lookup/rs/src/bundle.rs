//! The self-contained index: `rows.bin`, a block per family, read one family at a time.
//!
//! Block (format 4): "KFAB", acc u32, columns u32, rows u32, schemes u32, members u32, then seven lengths u64, then
//!   bits      columns x 20 bytes: quarter bits of a residue at a column, from the representatives' rows before collapsing
//!   heads     representatives x 10: accession 5 bytes, hmm_from u16, sequence length u16, flanks (left << 4 | right);
//!             then members x 5, the same without the accession (a member is known by its cluster)
//!   labels    schemes x rows u16: 0 is unlabelled
//!   ops       the rows' CIGARs, one after another: op << 6 | length - 1, with lengths over 62 as 63 and a u16; 0xFF ends a row.
//!             A member whose CIGAR is its parent's cut to the member's columns has 0xFD and a LEB128 column count
//!   seq       representatives' residues, flanks included, 5 bits each, each row from a whole byte; then each member's
//!             literal residues (see runs), the same way
//!   parents   per member, zigzag LEB128 of its parent's row less the previous member's: the last `members` rows are
//!             members of other proteins' clusters, sorted by cluster, and the parent is the row of the cluster's
//!             representative each is stored against, which names the cluster. The first stage
//!             scans representatives only; the members of the clusters it keeps join the column and banded stages.
//!   runs      per member, pairs (copied, literal) covering its residues in order: a copied residue is the parent's at
//!             the same place (`mapping`), a literal one is read from seq. A pair of under 63 copied and at most 3
//!             literal is one byte, literal << 6 | copied; otherwise 63, then both as LEB128. Members differ from their
//!             representative at a few residues, so most are copied.
//! All little-endian. Accessions fit 5 bytes: 1.1e12.
use crate::cand::hash;
use crate::data::*;
use crate::scan::{kadane, scan_batch};
use crate::sw;
use std::fmt::Write as _;
use std::fs::File;
use std::io::{Read, Seek, SeekFrom, Write};
use std::path::Path;
use std::sync::Mutex;

const MAGIC: &[u8; 4] = b"KFAB";
const HEAD: usize = 10;
const TOP: usize = 24 + 7 * 8; // the block's fixed header
const CLIPPED: u8 = 0xFD; // a member's CIGAR is its parent's, cut to the member's columns
const ESCAPE: u8 = 63; // a pair of runs that does not fit a byte

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
    /// Calls `f(column, position in seq)` for each residue on a match column, checking the row against a model of
    /// `len` columns.
    pub fn walk(&self, len: usize, mut f: impl FnMut(usize, usize)) -> Result<(), String> {
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
                    (0..n).for_each(|k| f(c + k, p + k));
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
        Ok(())
    }
    /// Residues on the model's columns, and each one's position in `seq` (u16::MAX: none).
    pub fn columns(&self, len: usize) -> Result<(Vec<u8>, Vec<u16>), String> {
        let (mut codes, mut pos) = (vec![NONE; len], vec![u16::MAX; len]);
        self.columns_into(&mut codes, &mut pos)?;
        Ok((codes, pos))
    }
    /// `columns` into buffers of the model's length, which are overwritten.
    pub fn columns_into(&self, codes: &mut [u8], pos: &mut [u16]) -> Result<(), String> {
        codes.fill(NONE);
        pos.fill(u16::MAX);
        self.walk(codes.len(), |c, p| {
            codes[c] = self.seq[p];
            pos[c] = p as u16;
        })
    }
}

/// Where each of a row's residues falls in another row (its parent), given the parent's position at each column
/// (`Row::columns`): through the columns both cover, and for flanks, inserts and columns the parent lacks, at the same
/// offset from the nearest such column before (after, for the left flank). None when the rows share no column.
pub fn mapping(r: &Row, parent_pos: &[u16]) -> Vec<Option<usize>> {
    let (mut c, mut p, mut anchors) = (r.hmm_from as usize - 1, r.left as usize, vec![]);
    for &(op, n) in &r.ops {
        for _ in 0..n {
            match op {
                0 => {
                    if let Some(&q) = parent_pos.get(c).filter(|&&q| q != u16::MAX) {
                        anchors.push((p, q as usize));
                    }
                    c += 1;
                    p += 1;
                }
                1 => p += 1,
                _ => c += 1,
            }
        }
    }
    let mut k = 0;
    (0..r.seq.len())
        .map(|p| {
            while k + 1 < anchors.len() && anchors[k + 1].0 <= p {
                k += 1;
            }
            let &(a, q) = anchors.get(k)?;
            (q + p).checked_sub(a)
        })
        .collect()
}

fn zigzag(x: i64) -> usize {
    ((x << 1) ^ (x >> 63)) as usize
}
fn unzigzag(x: usize) -> i64 {
    (x >> 1) as i64 ^ -((x & 1) as i64)
}

fn leb(mut x: usize, out: &mut Vec<u8>) {
    while x >= 128 {
        out.push(x as u8 | 128);
        x >>= 7;
    }
    out.push(x as u8);
}
fn unleb(b: &[u8], at: &mut usize) -> Option<usize> {
    let (mut x, mut sh) = (0usize, 0);
    loop {
        let y = *b.get(*at)?;
        *at += 1;
        x |= ((y & 127) as usize) << sh;
        if y < 128 {
            return Some(x);
        }
        sh += 7;
        if sh > 28 {
            return None;
        }
    }
}

/// A CIGAR's ops over columns c0..c1 only (0-based, `from` its first column): inserts kept only between two kept columns.
pub fn clip(ops: &[(u8, u16)], from: usize, c0: usize, c1: usize) -> Vec<(u8, u16)> {
    let (mut c, mut out): (usize, Vec<(u8, u16)>) = (from, vec![]);
    for &(op, n) in ops {
        let n = n as usize;
        let (a, b) = if op == 1 { (c, c) } else { (c.max(c0), (c + n).min(c1)) };
        let keep = if op == 1 { c > c0 && c < c1 } else { a < b };
        if keep {
            let k = if op == 1 { n } else { b - a };
            match out.last_mut() {
                Some(x) if x.0 == op => x.1 += k as u16,
                _ => out.push((op, k as u16)),
            }
        }
        if op != 1 {
            c += n;
        }
    }
    out
}
/// Columns a CIGAR covers (M and D).
fn span(ops: &[(u8, u16)]) -> usize {
    ops.iter().filter(|x| x.0 != 1).map(|x| x.1 as usize).sum()
}

fn put_pair(a: usize, b: usize, out: &mut Vec<u8>) {
    if a < ESCAPE as usize && b <= 3 {
        out.push((b as u8) << 6 | a as u8);
    } else {
        out.push(ESCAPE);
        leb(a, out);
        leb(b, out);
    }
}
fn get_pair(s: &[u8], at: &mut usize) -> Option<(usize, usize)> {
    let x = *s.get(*at)?;
    *at += 1;
    if x == ESCAPE {
        Some((unleb(s, at)?, unleb(s, at)?))
    } else {
        Some(((x & 63) as usize, (x >> 6) as usize))
    }
}

/// A member against its parent: the runs, and the literal residues.
fn diff(r: &Row, parent: &Row, parent_pos: &[u16]) -> (Vec<u8>, Vec<u8>) {
    let m = mapping(r, parent_pos);
    let same: Vec<bool> = (0..r.seq.len()).map(|p| m[p].is_some_and(|q| parent.seq.get(q) == Some(&r.seq[p]))).collect();
    let (mut runs, mut lit, mut p) = (vec![], vec![], 0);
    while p < same.len() {
        let a = same[p..].iter().take_while(|&&x| x).count();
        let b = same[p + a..].iter().take_while(|&&x| !x).count();
        put_pair(a, b, &mut runs);
        lit.extend_from_slice(&r.seq[p + a..p + a + b]);
        p += a + b;
    }
    (runs, lit)
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

/// `rows`: representatives, then members sorted by cluster, each with a representative row of its cluster. A member is
/// stored against whichever of its representative's rows leaves the fewest literal residues.
fn block(acc: u32, len: usize, schemes: usize, rows: &[Row], bits: &[u8]) -> Vec<u8> {
    let (mut heads, mut labels, mut ops, mut seq, mut parents, mut runs) = (vec![], vec![], vec![], vec![], vec![], vec![]);
    let members = rows.iter().filter(|r| r.cluster > 0).count();
    let nr = rows.len() - members;
    assert!(rows[..nr].iter().all(|r| r.cluster == 0), "members must come last");
    assert!(rows[nr..].is_sorted_by_key(|r| r.cluster), "members must be sorted by cluster");
    let mut of: std::collections::HashMap<u64, Vec<usize>> = Default::default();
    rows[..nr].iter().enumerate().for_each(|(i, r)| of.entry(r.acc).or_default().push(i));
    let mut ppos: Vec<Option<Vec<u16>>> = vec![None; nr];
    let mut last = 0i64;
    for r in rows {
        assert!(r.left < 16 && r.right < 16 && r.acc < 1 << 40 && r.seq.len() < 65536);
        if r.cluster == 0 {
            heads.extend_from_slice(&r.acc.to_le_bytes()[..5]);
        }
        heads.extend_from_slice(&r.hmm_from.to_le_bytes());
        heads.extend_from_slice(&(r.seq.len() as u16).to_le_bytes());
        heads.push(r.left << 4 | r.right);
        if r.cluster == 0 {
            pack_ops(&r.ops, &mut ops);
            pack5(&r.seq, &mut seq);
            continue;
        }
        let cands = of.get(&r.cluster).unwrap_or_else(|| panic!("PF{acc:05}: {}'s representative {} has no row", r.acc, r.cluster));
        let (i, (ru, lit)) = cands.iter().map(|&i| {
            let pp = ppos[i].get_or_insert_with(|| rows[i].columns(len).expect("a representative's row").1);
            (i, diff(r, &rows[i], pp))
        }).min_by_key(|(i, (_, lit))| (lit.len(), *i)).unwrap();
        let (c0, w) = (r.hmm_from as usize - 1, span(&r.ops));
        if clip(&rows[i].ops, rows[i].hmm_from as usize - 1, c0, c0 + w) == r.ops {
            ops.push(CLIPPED);
            leb(w, &mut ops);
        } else {
            pack_ops(&r.ops, &mut ops);
        }
        leb(zigzag(i as i64 - last), &mut parents);
        last = i as i64;
        runs.extend_from_slice(&ru);
        pack5(&lit, &mut seq);
    }
    for s in 0..schemes {
        rows.iter().for_each(|r| labels.extend_from_slice(&r.labels[s].to_le_bytes()));
    }
    let mut o = MAGIC.to_vec();
    for x in [acc, len as u32, rows.len() as u32, schemes as u32, members as u32] {
        o.extend_from_slice(&x.to_le_bytes());
    }
    let parts: [&[u8]; 7] = [bits, &heads, &labels, &ops, &seq, &parents, &runs];
    parts.iter().for_each(|p| o.extend_from_slice(&(p.len() as u64).to_le_bytes()));
    parts.iter().for_each(|p| o.extend_from_slice(p));
    o
}

#[cfg_attr(not(test), allow(dead_code))] // tests: the whole block, against `layout` and its parts
pub struct Block {
    pub acc: u32,
    pub len: usize,
    pub reps: usize, // rows[..reps] are representatives, rows[reps..] members
    pub bits: Vec<u8>,
    pub rows: Vec<Row>,
}
/// The whole block decoded: `layout`, then every representative and every cluster's members.
#[cfg_attr(not(test), allow(dead_code))]
pub fn read_block(b: &[u8]) -> Result<Block, String> {
    let l = layout(b)?;
    let mut rows = (0..l.reps).map(|i| l.rep(i)).collect::<Result<Vec<_>, _>>()?;
    for k in 0..l.clusters.len() {
        let m = l.members(k, &rows)?;
        rows.extend(m);
    }
    Ok(Block { acc: l.acc, len: l.len, reps: l.reps, bits: l.part[0].to_vec(), rows })
}

/// A cluster with members: its representative's accession, its members (indices among the block's members), and
/// where in the parts its first member's bytes start (CIGARs, residues, parents, runs) with the parent row before it.
#[derive(Clone, Copy, Debug)]
pub struct Cluster {
    pub acc: u64,
    pub first: usize,
    pub n: usize,
    at: [usize; 4],
    last: i64,
}

/// Where a block's rows are, from one pass that decodes nothing: representatives can then be decoded in parallel
/// (`rep`), and a cluster's members only when a read keeps the cluster (`members`).
pub struct Layout<'a> {
    pub acc: u32,
    pub len: usize,
    pub n: usize,
    pub reps: usize,
    schemes: usize,
    part: Vec<&'a [u8]>,
    rep_at: Vec<(usize, usize)>, // a representative's CIGAR and residues
    pub clusters: Vec<Cluster>,  // in the block's order, which is by accession
}

fn skip_ops(b: &[u8], at: &mut usize) -> Option<()> {
    loop {
        let x = *b.get(*at)?;
        *at += 1;
        if x == 0xFF {
            return Some(());
        }
        if x & 63 == 63 {
            *at += 2;
        }
    }
}

pub fn layout(b: &[u8]) -> Result<Layout<'_>, String> {
    let u = |at: usize| -> u32 { u32::from_le_bytes(b[at..at + 4].try_into().unwrap()) };
    if b.len() < TOP || &b[..4] != MAGIC {
        return Err("not a block".into());
    }
    let (acc, len, n, schemes, members) = (u(4), u(8) as usize, u(12) as usize, u(16) as usize, u(20) as usize);
    let sizes: Vec<usize> = (0..7).map(|i| u64::from_le_bytes(b[24 + 8 * i..32 + 8 * i].try_into().unwrap()) as usize).collect();
    if TOP + sizes.iter().sum::<usize>() != b.len() || sizes[0] != 20 * len || sizes[1] != HEAD * n - 5 * members || sizes[2] != 2 * n * schemes
        || members > n {
        return Err(format!("PF{acc:05}: the block's sizes do not add up"));
    }
    let mut at = TOP;
    let part: Vec<&[u8]> = sizes.iter().map(|&s| { at += s; &b[at - s..at] }).collect();
    let nr = n - members;
    let mut l = Layout { acc, len, n, reps: nr, schemes, part, rep_at: Vec::with_capacity(nr), clusters: vec![] };
    let (mut o, mut q, mut ra, mut pa, mut last) = (0, 0, 0, 0, 0i64);
    for i in 0..n {
        let sl = l.head(i).2;
        if i < nr {
            l.rep_at.push((o, q));
            skip_ops(l.part[3], &mut o).ok_or("short CIGARs")?;
            q += (5 * sl).div_ceil(8);
            continue;
        }
        let state = ([o, q, pa, ra], last);
        if l.part[3].get(o) == Some(&CLIPPED) {
            o += 1;
            unleb(l.part[3], &mut o).ok_or("short CIGARs")?;
        } else {
            skip_ops(l.part[3], &mut o).ok_or("short CIGARs")?;
        }
        let pi = last + unzigzag(unleb(l.part[5], &mut pa).ok_or("short parents")?);
        last = pi;
        let pi = usize::try_from(pi).ok().filter(|&x| x < nr).ok_or_else(|| format!("PF{acc:05}: a member's parent is not a representative"))?;
        let cl = l.head(pi).0;
        match l.clusters.last_mut() {
            Some(c) if c.acc == cl => c.n += 1,
            _ => l.clusters.push(Cluster { acc: cl, first: i - nr, n: 1, at: state.0, last: state.1 }),
        }
        let (mut p, mut lit) = (0, 0);
        while p < sl {
            let (a, e) = get_pair(l.part[6], &mut ra).ok_or("short runs")?;
            if p + a + e > sl {
                return Err(format!("PF{acc:05}: a member's runs pass its end"));
            }
            lit += e;
            p += a + e;
        }
        q += (5 * lit).div_ceil(8);
    }
    if o > l.part[3].len() || q > l.part[4].len() {
        return Err(format!("PF{acc:05}: short {}", if q > l.part[4].len() { "sequences" } else { "CIGARs" }));
    }
    if o != l.part[3].len() || q != l.part[4].len() || pa != l.part[5].len() || ra != l.part[6].len() {
        return Err(format!("PF{acc:05}: bytes left over"));
    }
    if !l.clusters.is_sorted_by_key(|c| c.acc) {
        return Err(format!("PF{acc:05}: members are not sorted by cluster"));
    }
    Ok(l)
}

impl Layout<'_> {
    /// Row i's accession (0 for a member: a member's head has none, only its last 5 bytes), hmm_from, residues and
    /// flanks.
    fn head(&self, i: usize) -> (u64, u16, usize, u8) {
        let nr = self.reps;
        let (acc, h) = match i.checked_sub(nr) {
            None => (self.part[1][HEAD * i..HEAD * i + 5].iter().rev().fold(0, |a, &x| a << 8 | x as u64), &self.part[1][HEAD * i + 5..HEAD * (i + 1)]),
            Some(j) => (0, &self.part[1][HEAD * nr + 5 * j..HEAD * nr + 5 * (j + 1)]),
        };
        (acc, u16::from_le_bytes([h[0], h[1]]), u16::from_le_bytes([h[2], h[3]]) as usize, h[4])
    }
    fn row_head(&self, i: usize) -> Row {
        let (acc, hmm_from, sl, fl) = self.head(i);
        let n = self.n;
        let part = &self.part[2];
        Row {
            acc,
            hmm_from,
            left: fl >> 4,
            right: fl & 15,
            ops: vec![],
            seq: vec![0; sl],
            labels: (0..self.schemes).map(|s| { let x = 2 * (s * n + i); u16::from_le_bytes([part[x], part[x + 1]]) }).collect(),
            cluster: 0,
        }
    }
    /// Representative i.
    pub fn rep(&self, i: usize) -> Result<Row, String> {
        let mut r = self.row_head(i);
        let (mut o, q) = self.rep_at[i];
        r.ops = unpack_ops(self.part[3], &mut o).ok_or("short CIGARs")?;
        let bytes = self.part[4].get(q..q + (5 * r.seq.len()).div_ceil(8)).ok_or("short sequences")?;
        r.seq = unpack5(bytes, r.seq.len());
        Ok(r)
    }
    /// The members of cluster k, given the representatives (`rep`, all of them).
    pub fn members(&self, k: usize, reps: &[Row]) -> Result<Vec<Row>, String> {
        let (acc, len, nr, part) = (self.acc, self.len, self.reps, &self.part);
        let c = self.clusters[k];
        let ([mut o, mut q, mut pa, mut ra], mut last) = (c.at, c.last);
        let mut ppos: Vec<(usize, Vec<u16>)> = vec![]; // the parents' positions by column: a cluster has few parents
        let mut out = Vec::with_capacity(c.n);
        for j in c.first..c.first + c.n {
            let mut r = self.row_head(nr + j);
            let sl = r.seq.len();
            let mut clipped = None;
            if part[3].get(o) == Some(&CLIPPED) {
                o += 1; // the parent's, cut below
                clipped = Some(unleb(part[3], &mut o).ok_or("short CIGARs")?);
            } else {
                r.ops = unpack_ops(part[3], &mut o).ok_or("short CIGARs")?;
            }
            let pi = last + unzigzag(unleb(part[5], &mut pa).ok_or("short parents")?);
            last = pi;
            let pi = usize::try_from(pi).ok().filter(|&x| x < nr).ok_or_else(|| format!("PF{acc:05}: a member's parent is not a representative"))?;
            let parent = &reps[pi];
            r.cluster = parent.acc;
            if let Some(w) = clipped {
                let c0 = r.hmm_from as usize - 1;
                r.ops = clip(&parent.ops, parent.hmm_from as usize - 1, c0, c0 + w);
                if span(&r.ops) != w {
                    return Err(format!("PF{acc:05}: a member's columns are not all its parent's"));
                }
            }
            let at = match ppos.iter().position(|x| x.0 == pi) {
                Some(x) => x,
                None => {
                    ppos.push((pi, parent.columns(len)?.1));
                    ppos.len() - 1
                }
            };
            let m = mapping(&r, &ppos[at].1);
            let (mut p, mut literal_at) = (0, vec![]);
            while p < sl {
                let (a, e) = get_pair(part[6], &mut ra).ok_or("short runs")?;
                if p + a + e > sl {
                    return Err(format!("PF{acc:05}: a member's runs pass its end"));
                }
                for k in p..p + a {
                    r.seq[k] = *m[k].and_then(|x| parent.seq.get(x)).ok_or_else(|| format!("PF{acc:05}: a copied residue outside the parent"))?;
                }
                literal_at.extend(p + a..p + a + e);
                p += a + e;
            }
            let bytes = part[4].get(q..q + (5 * literal_at.len()).div_ceil(8)).ok_or("short sequences")?;
            q += bytes.len();
            let got = unpack5(bytes, literal_at.len());
            for (k, x) in literal_at.into_iter().zip(got) {
                r.seq[k] = x;
            }
            out.push(r);
        }
        Ok(out)
    }
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

/// Reads to a batch of the first stage (`scan_batch`).
const BATCH: usize = 16;

/// Asks for the cache line at `p` ahead of its use: the column scores read rows at random places.
#[inline(always)]
fn prefetch(p: *const u8) {
    // SAFETY: a prefetch is a hint and never faults, whatever the address
    #[cfg(target_arch = "x86_64")]
    unsafe {
        std::arch::x86_64::_mm_prefetch(p as *const i8, std::arch::x86_64::_MM_HINT_T0)
    }
    #[cfg(target_arch = "aarch64")]
    unsafe {
        std::arch::asm!("prfm pldl1keep, [{0}]", in(reg) p, options(nostack, readonly, preserves_flags))
    }
}
const AHEAD: usize = 16; // rows prefetched ahead

/// A pointer that threads may write through at places no other thread writes.
struct Shared(*mut u8);
unsafe impl Sync for Shared {}

/// bundle-query DIR QUERIES OUT K SECOND [MARGIN]: QUERIES is the same stream with no schemes, a row's accession being
/// the query's id. OUT: qid, score, self score, tied clusters, accession of the first, then each scheme's commonest
/// label among the tied rows. A cluster is named by its representative's accession, also when a member's row scored best.
/// The banded score sees the best SECOND rows by column score, and with MARGIN > 0 only those within MARGIN of the
/// read's best column score.
pub fn query(dir: &Path, queries: &Path, out: &Path, k: usize, second: usize, margin: i16, threads: usize) {
    use std::sync::atomic::{AtomicU64, AtomicUsize, Ordering::Relaxed};
    use std::sync::OnceLock;
    use std::time::Instant;
    let index: Vec<Vec<u64>> = std::fs::read_to_string(dir.join("families.tsv")).unwrap().lines()
        .map(|l| l.split('\t').map(|x| x.parse().unwrap()).collect()).collect();
    let mut f = File::open(dir.join("rows.bin")).unwrap();
    let mut i = In(std::io::BufReader::new(File::open(queries).unwrap()));
    let bl = blosum();
    let mut o = String::new();
    let (mut t_read, mut t_decode, mut t_search) = (0.0, 0.0, 0.0);
    // CPU time by stage, summed over threads, in ns: first stage, members decoded, column scores, banded
    let cpu: [AtomicU64; 4] = Default::default();
    // reads searched, and rows they took to the column and banded scores
    let rows: [AtomicU64; 3] = Default::default();
    let ns = |t: Instant| t.elapsed().as_nanos() as u64;
    while let Some(acc) = i.u32() {
        let (len, n, _) = (i.u32().unwrap() as usize, i.u32().unwrap() as usize, i.u32().unwrap());
        let qs: Vec<Row> = (0..n).map(|_| i.row(0).expect("the queries end inside a family")).collect();
        let Some(x) = index.iter().find(|x| x[0] == acc as u64) else { continue }; // a family with no row in the index
        assert_eq!(x[1] as usize, len, "PF{acc:05}: the reads' model has {len} columns, the index's {}", x[1]);
        let t0 = Instant::now();
        let mut b = vec![0; x[5] as usize];
        f.seek(SeekFrom::Start(x[4])).unwrap();
        f.read_exact(&mut b).expect("rows.bin is shorter than families.tsv says");
        let t1 = Instant::now();
        let lay = layout(&b).unwrap_or_else(|e| panic!("{e}"));
        assert_eq!((lay.acc, lay.len), (acc, len));
        let nr = lay.reps; // the first stage sees representatives only
        let reps: Vec<Row> = par_map(nr, if nr < 4096 { 1 } else { threads }, |r| {
            let x = lay.rep(r).and_then(|x| x.walk(len, |_, _| {}).map(|_| x));
            x.unwrap_or_else(|e| panic!("PF{acc:05} {e}"))
        });
        let bits = &lay.part[0];
        // the representatives' residues by column, each row written by one thread
        // the representatives' residues by column (the first stage) and by row (the column scores), each row written
        // by one thread
        let (mut m, mut rm) = (vec![NONE; nr * len], vec![NONE; nr * len]);
        let (mp, rp) = (Shared(m.as_mut_ptr()), Shared(rm.as_mut_ptr()));
        parallel(nr.div_ceil(4096), threads.min(nr.div_ceil(4096)), |t| {
            let (mp, rp) = (&mp, &rp);
            for (r, x) in reps.iter().enumerate().skip(t * 4096).take(4096) {
                // SAFETY: walk checked the columns (c < len), and row r is this thread's only
                x.walk(len, |c, p| unsafe {
                    *mp.0.add(c * nr + r) = x.seq[p];
                    *rp.0.add(r * len + c) = x.seq[p];
                }).unwrap();
            }
        });
        // a cluster's members and their residues by column, decoded when a read first keeps the cluster
        let members: Vec<OnceLock<(Vec<Row>, Vec<u8>)>> = (0..lay.clusters.len()).map(|_| OnceLock::new()).collect();
        let cluster = |k: usize| -> &(Vec<Row>, Vec<u8>) {
            members[k].get_or_init(|| {
                let t = Instant::now();
                let x = lay.members(k, &reps).unwrap_or_else(|e| panic!("{e}"));
                let mut codes = vec![NONE; x.len() * len];
                for (r, ch) in x.iter().zip(codes.chunks_mut(len)) {
                    r.walk(len, |c, p| ch[c] = r.seq[p]).unwrap_or_else(|e| panic!("PF{acc:05} {e}"));
                }
                cpu[1].fetch_add(ns(t), Relaxed);
                (x, codes)
            })
        };
        let member = |j: usize| -> &Row {
            let k = lay.clusters.partition_point(|c| c.first + c.n <= j);
            &cluster(k).0[j - lay.clusters[k].first]
        };
        let row_of = |r: usize| if r < nr { &reps[r] } else { member(r - nr) };
        let t2 = Instant::now();
        let lines = Mutex::new(vec![]);
        let next = AtomicUsize::new(0);
        parallel(threads, threads, |_| {
            let (mut buf, mut mine, mut x, mut tp) = (vec![], vec![], vec![NONE; len], vec![0u16; len]);
            loop {
                let lo = next.fetch_add(BATCH, Relaxed);
                let Some(batch) = qs.get(lo..(lo + BATCH).min(qs.len())) else { break };
                if batch.is_empty() {
                    break;
                }
                let t = Instant::now();
                let rp: Vec<(Vec<u8>, Vec<u16>)> = batch.iter()
                    .map(|q| q.columns(len).unwrap_or_else(|e| panic!("PF{acc:05}, read {e}: another Pfam release?"))).collect();
                let cws: Vec<Vec<(usize, u8, u16)>> = rp.iter()
                    .map(|(row, _)| (0..len).filter(|&c| row[c] < 20).map(|c| (c, row[c], bits[20 * c + row[c] as usize] as u16)).collect())
                    .collect();
                let kept = scan_batch(&m, nr, &cws, k.clamp(1, nr.max(1)));
                cpu[0].fetch_add(ns(t), Relaxed);
                for ((q, (row, pos)), keep) in batch.iter().zip(&rp).zip(kept) {
                    let t = Instant::now();
                    let cols: Vec<usize> = (0..len).filter(|&c| row[c] < 20).collect();
                    let col_score = |x: &[u8]| kadane(cols.iter().map(|&c| if x[c] < 21 { bl[row[c] as usize][x[c] as usize] } else { -4 }));
                    // (column score, row): representatives, then the members of their clusters, a cluster at a time
                    let (c0, c1) = (cols.first().copied().unwrap_or(0), cols.last().copied().unwrap_or(0));
                    let mut v: Vec<(i16, usize)> = keep.iter().enumerate().map(|(i, &r)| {
                        if let Some(&n) = keep.get(i + AHEAD) {
                            prefetch(rm[n as usize * len + c0..].as_ptr());
                            prefetch(rm[n as usize * len + c1..].as_ptr());
                        }
                        (col_score(&rm[r as usize * len..(r as usize + 1) * len]), r as usize)
                    }).collect();
                    let mut members_t = 0;
                    if !lay.clusters.is_empty() {
                        let mut us: Vec<u64> = keep.iter().map(|&r| reps[r as usize].acc).collect();
                        us.sort_unstable();
                        us.dedup();
                        for u in us {
                            if let Ok(k) = lay.clusters.binary_search_by_key(&u, |c| c.acc) {
                                let tm = Instant::now();
                                let codes = &cluster(k).1;
                                members_t += ns(tm);
                                let first = nr + lay.clusters[k].first;
                                v.extend(codes.chunks(len).enumerate().map(|(t, x)| (col_score(x), first + t)));
                            }
                        }
                    }
                    rows[0].fetch_add(1, Relaxed);
                    rows[1].fetch_add(v.len() as u64, Relaxed);
                    v.sort_unstable_by(|a, b| b.cmp(a));
                    let mut cut = v.get(second.min(v.len()).saturating_sub(1)).map_or(0, |x| x.0);
                    if margin > 0 {
                        cut = cut.max(v.first().map_or(0, |x| x.0).saturating_sub(margin));
                    }
                    rows[2].fetch_add(v.iter().take_while(|x| x.0 >= cut).count() as u64, Relaxed);
                    cpu[2].fetch_add(ns(t).saturating_sub(members_t), Relaxed);
                    let t = Instant::now();
                    let mut best: (i16, Vec<usize>) = (0, vec![]);
                    for &(_, r) in v.iter().take_while(|x| x.0 >= cut) {
                        let rw = row_of(r);
                        rw.columns_into(&mut x, &mut tp).unwrap();
                        let mut d: Vec<i32> = cols.iter().filter(|&&c| x[c] == row[c]).map(|&c| tp[c] as i32 - pos[c] as i32).collect();
                        d.sort_unstable();
                        let Some(&d) = d.get(d.len() / 2) else { continue };
                        let tt = &rw.seq;
                        let a = (d - 8).clamp(0, tt.len() as i32);
                        let e = (q.seq.len() as i32 + d + 8).clamp(a, tt.len() as i32);
                        let sc = sw::score(&q.seq, &tt[a as usize..e as usize], &bl, Some((d - 8 - a, d + 8 - a)), &mut buf);
                        if sc > best.0 {
                            best = (sc, vec![r]);
                        } else if sc == best.0 && sc > 0 {
                            best.1.push(r);
                        }
                    }
                    cpu[3].fetch_add(ns(t), Relaxed);
                    let own: i32 = q.seq.iter().filter(|&&a| a < 20).map(|&a| bl[a as usize][a as usize] as i32).sum();
                    let mut tied: Vec<u64> = best.1.iter().map(|&r| row_of(r).unit()).collect();
                    tied.sort_unstable_by_key(|&u| (hash(u), u)); // the cluster named is any of the tied, always the same
                    tied.dedup();
                    let mut l = format!("{}\t{}\t{own}\t{}\t{}", q.acc, best.0, tied.len(), tied.first().map_or(0, |&u| u));
                    for sch in 0..if lay.n > 0 { lay.schemes } else { 0 } {
                        let mut c: std::collections::BTreeMap<u16, usize> = Default::default();
                        best.1.iter().for_each(|&r| *c.entry(row_of(r).labels[sch]).or_default() += 1);
                        // commonest; a label before none, then the lower number, on a tie
                        let top = c.iter().max_by_key(|(&x, &n)| (n, x > 0, std::cmp::Reverse(x))).map_or(0, |x| *x.0);
                        write!(l, "\t{top}").unwrap();
                    }
                    mine.push((q.acc, l));
                }
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
    let c: Vec<f64> = cpu.iter().map(|x| x.load(Relaxed) as f64 / 1e9).collect();
    let r: Vec<f64> = rows.iter().map(|x| x.load(Relaxed) as f64).collect();
    eprintln!("read {t_read:.2} s, decode {t_decode:.2} s, search {t_search:.2} s on {threads} threads; \
               search CPU: first stage {:.1} s, members {:.1} s, columns {:.1} s, banded {:.1} s; \
               rows a read: {:.0} to the column score, {:.0} to the banded", c[0], c[1], c[2], c[3], r[1] / r[0].max(1.0), r[2] / r[0].max(1.0));
}

#[cfg(test)]
mod tests {
    use super::*;
    /// What reading gives back: members without their accessions.
    fn known(rows: &[Row]) -> Vec<Row> {
        rows.iter().map(|r| Row { acc: if r.cluster > 0 { 0 } else { r.acc }, ..r.clone() }).collect()
    }
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
        assert_eq!((x.acc, x.len, &x.bits, &x.rows), (44, 400, &bits, &known(&r)));
        assert!(read_block(&b[..b.len() - 1]).is_err(), "a cut block must not pass");
        b[16] = 1;
        assert!(read_block(&b).is_err(), "nor a wrong number of schemes");
        assert_eq!(x.reps, 40);
        assert_eq!(r[45].unit(), r[6].acc);
        assert_eq!(r[5].unit(), r[5].acc);
    }
    #[test]
    fn members_are_stored_as_differences() {
        // a representative of 3 domains' worth of residues, and a member: 4 substitutions, an insert of 2 where the
        // representative has none, a start 5 columns later, and other flanks
        let rnd = |n: u64, s: u64| -> Vec<u8> { (0..n).map(|k| (hash(s * 7919 + k) % 20) as u8).collect() };
        let core = rnd(120, 1);
        let rep = Row { acc: 7, hmm_from: 1, left: 10, right: 10, ops: vec![(0, 120)], labels: vec![],
                        seq: [rnd(10, 2), core.clone(), rnd(10, 3)].concat(), cluster: 0 };
        let mut mc = core[5..].to_vec();
        for k in [3, 40, 41, 90] {
            mc[k] = (mc[k] + 1) % 20;
        }
        mc.splice(60..60, [4u8, 9]);
        let mem = Row { acc: 8, hmm_from: 6, left: 7, right: 10, ops: vec![(0, 60), (1, 2), (0, 55)], labels: vec![],
                        seq: [rnd(7, 4), mc, rep.seq[130..140].to_vec()].concat(), cluster: 7 };
        assert!(mem.columns(120).is_ok());
        let v = vec![rep.clone(), mem.clone()];
        let b = block(1, 120, 0, &v, &[0u8; 2400]);
        let x = read_block(&b).unwrap();
        assert_eq!(x.rows, known(&v));
        // literals: 7 left flank + 4 substituted + 2 inserted; the right flank is the representative's, so copied
        let (_, lit) = diff(&mem, &rep, &rep.columns(120).unwrap().1);
        assert_eq!(lit.len(), 13);
        let plain = block(1, 120, 0, &[rep.clone(), Row { cluster: 0, ..mem.clone() }], &[0u8; 2400]);
        assert!(b.len() + 40 < plain.len(), "{} against {} stored plainly", b.len(), plain.len());
        let mut bad = b.clone();
        let at = bad.len() - 1;
        bad[at] = 0x7f; // the last run: past the member's end
        assert!(read_block(&bad).is_err());
    }
    #[test]
    fn cigars_and_runs_as_differences() {
        // columns 0..10 M, insert 3, 10..14 deleted, 14..30 M; cut to columns 5..20
        let ops = vec![(0, 10), (1, 3), (2, 4), (0, 16)];
        assert_eq!(clip(&ops, 0, 5, 20), vec![(0, 5), (1, 3), (2, 4), (0, 6)]);
        assert_eq!(clip(&ops, 0, 10, 20), vec![(2, 4), (0, 6)], "an insert at the cut is dropped");
        assert_eq!(clip(&ops, 0, 0, 30), ops);
        let mut b = vec![];
        for (a, e) in [(0, 0), (62, 3), (63, 0), (5, 4), (1000, 2)] {
            put_pair(a, e, &mut b);
        }
        assert_eq!(b.len(), 2 + 3 + 3 + 4, "two pairs in a byte each");
        let mut at = 0;
        let got: Vec<_> = (0..5).map(|_| get_pair(&b, &mut at).unwrap()).collect();
        assert_eq!(got, vec![(0, 0), (62, 3), (63, 0), (5, 4), (1000, 2)]);
        // a member whose CIGAR is its representative's cut to its columns: 2 bytes of CIGAR, and back exactly
        let s: Vec<u8> = (0..60u64).map(|k| (hash(k) % 20) as u8).collect();
        let rep = Row { acc: 3, hmm_from: 1, left: 0, right: 0, ops: vec![(0, 20), (1, 2), (0, 38)], seq: s.clone(), labels: vec![], cluster: 0 };
        let mem = Row { acc: 4, hmm_from: 6, left: 0, right: 0, ops: vec![(0, 15), (1, 2), (0, 30)], seq: s[5..52].to_vec(), labels: vec![], cluster: 3 };
        let v = vec![rep.clone(), mem];
        let blk = block(1, 58, 0, &v, &[0u8; 20 * 58]);
        assert_eq!(read_block(&blk).unwrap().rows, known(&v));
        let ops_part = u64::from_le_bytes(blk[24 + 3 * 8..32 + 3 * 8].try_into().unwrap());
        let rep_ops = { let mut x = vec![]; pack_ops(&rep.ops, &mut x); x.len() as u64 };
        assert_eq!(ops_part, rep_ops + 2);
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
