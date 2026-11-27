//! The flat files `prep.py pack` writes.
use std::fs;
use std::path::Path;

pub const NONE: u8 = 21; // no residue in the column
pub const AA: usize = 22;

/// BLOSUM62 in the order ACDEFGHIKLMNPQRSTVWY, then X (20) and none (21).
pub fn blosum() -> [[i16; AA]; AA] {
    const B: [&str; 20] = [
        "4 0 -2 -1 -2 0 -2 -1 -1 -1 -1 -2 -1 -1 -1 1 0 0 -3 -2",
        "0 9 -3 -4 -2 -3 -3 -1 -3 -1 -1 -3 -3 -3 -3 -1 -1 -1 -2 -2",
        "-2 -3 6 2 -3 -1 -1 -3 -1 -4 -3 1 -1 0 -2 0 -1 -3 -4 -3",
        "-1 -4 2 5 -3 -2 0 -3 1 -3 -2 0 -1 2 0 0 -1 -2 -3 -2",
        "-2 -2 -3 -3 6 -3 -1 0 -3 0 0 -3 -4 -3 -3 -2 -2 -1 1 3",
        "0 -3 -1 -2 -3 6 -2 -4 -2 -4 -3 0 -2 -2 -2 0 -2 -3 -2 -3",
        "-2 -3 -1 0 -1 -2 8 -3 -1 -3 -2 1 -2 0 0 -1 -2 -3 -2 2",
        "-1 -1 -3 -3 0 -4 -3 4 -3 2 1 -3 -3 -3 -3 -2 -1 3 -3 -1",
        "-1 -3 -1 1 -3 -2 -1 -3 5 -2 -1 0 -1 1 2 0 -1 -2 -3 -2",
        "-1 -1 -4 -3 0 -4 -3 2 -2 4 2 -3 -3 -2 -2 -2 -1 1 -2 -1",
        "-1 -1 -3 -2 0 -3 -2 1 -1 2 5 -2 -2 0 -1 -1 -1 1 -1 -1",
        "-2 -3 1 0 -3 0 1 -3 0 -3 -2 6 -2 0 0 1 0 -3 -4 -2",
        "-1 -3 -1 -1 -4 -2 -2 -3 -1 -3 -2 -2 7 -1 -2 -1 -1 -2 -4 -3",
        "-1 -3 0 2 -3 -2 0 -3 1 -2 0 0 -1 5 1 0 -1 -2 -2 -1",
        "-1 -3 -2 0 -3 -2 0 -3 2 -2 -1 0 -2 1 5 -1 -1 -3 -3 -2",
        "1 -1 0 0 -2 0 -1 -2 0 -2 -1 1 -1 0 -1 4 1 -2 -3 -2",
        "0 -1 -1 -1 -2 -2 -2 -1 -1 -1 -1 0 -1 -1 -1 1 5 0 -2 -2",
        "0 -1 -3 -2 -1 -3 -3 3 -2 1 1 -3 -2 -2 -3 -2 0 4 -3 -1",
        "-3 -2 -4 -3 1 -2 -2 -3 -3 -2 -1 -4 -4 -2 -3 -3 -2 -3 11 2",
        "-2 -2 -3 -2 3 -3 2 -1 -2 -1 -1 -2 -3 -1 -2 -2 -2 -1 2 7",
    ];
    let mut m = [[-1i16; AA]; AA];
    for (i, l) in B.iter().enumerate() {
        for (j, v) in l.split(' ').enumerate() {
            m[i][j] = v.parse().unwrap();
        }
    }
    for i in 0..AA {
        m[i][21] = -50;
        m[21][i] = -50;
    }
    m
}

pub fn bytes(p: &Path) -> Vec<u8> {
    fs::read(p).unwrap_or_else(|e| panic!("{}: {e}", p.display()))
}
pub fn u16s(p: &Path) -> Vec<u16> {
    bytes(p).chunks_exact(2).map(|c| u16::from_le_bytes([c[0], c[1]])).collect()
}
pub fn u32s(p: &Path) -> Vec<u32> {
    bytes(p).chunks_exact(4).map(|c| u32::from_le_bytes(c.try_into().unwrap())).collect()
}
pub fn u64s(p: &Path) -> Vec<u64> {
    bytes(p).chunks_exact(8).map(|c| u64::from_le_bytes(c.try_into().unwrap())).collect()
}

pub struct Seqs {
    pub off: Vec<u64>,
    pub seq: Vec<u8>,
}
impl Seqs {
    pub fn load(dir: &Path, name: &str) -> Self {
        Seqs { off: u64s(&dir.join(format!("{name}.off"))), seq: bytes(&dir.join(format!("{name}.seq"))) }
    }
    pub fn len(&self) -> usize {
        self.off.len() - 1
    }
    pub fn get(&self, i: usize) -> &[u8] {
        &self.seq[self.off[i] as usize..self.off[i + 1] as usize]
    }
}

pub struct Family {
    pub acc: u32,
    pub len: usize,
    pub reps: usize,
    pub members: usize,
}
pub fn families(dir: &Path) -> Vec<Family> {
    fs::read_to_string(dir.join("families.tsv"))
        .unwrap()
        .lines()
        .map(|l| {
            let c: Vec<usize> = l.split('\t').map(|x| x.parse().unwrap()).collect();
            Family { acc: c[0] as u32, len: c[1], reps: c[2], members: c[3] }
        })
        .collect()
}

/// One arm's rows ("reps" or "members"), all families, in memory.
pub struct Rows {
    pub codes: Vec<u8>,
    pub pos: Vec<u16>,
    pub unit: Vec<u32>,
    pub row0: Vec<usize>,  // first row of each family
    pub byte0: Vec<usize>, // and its first byte
}
impl Rows {
    pub fn load(dir: &Path, arm: &str, fams: &[Family]) -> Self {
        let (mut row0, mut byte0) = (vec![0], vec![0]);
        for f in fams {
            let n = if arm == "reps" { f.reps } else { f.members };
            row0.push(row0.last().unwrap() + n);
            byte0.push(byte0.last().unwrap() + n * f.len);
        }
        let r = Rows {
            codes: bytes(&dir.join(format!("{arm}.rows"))),
            pos: u16s(&dir.join(format!("{arm}.pos"))),
            unit: u32s(&dir.join(format!("{arm}.unit"))),
            row0,
            byte0,
        };
        assert_eq!(r.codes.len(), *r.byte0.last().unwrap());
        assert_eq!(r.pos.len(), r.codes.len());
        assert_eq!(r.unit.len(), *r.row0.last().unwrap());
        r
    }
}

pub struct Queries {
    pub seqs: Seqs,
    pub fam: Vec<u32>,
    pub codes: Vec<u8>,
    pub pos: Vec<u16>,
    pub byte0: Vec<usize>,
}
impl Queries {
    pub fn load(dir: &Path, fams: &[Family]) -> Self {
        let fam = u32s(&dir.join("q.fam"));
        let mut byte0 = vec![0];
        for &f in &fam {
            byte0.push(byte0.last().unwrap() + fams[f as usize].len);
        }
        let q = Queries {
            seqs: Seqs::load(dir, "q"),
            fam,
            codes: bytes(&dir.join("q.rows")),
            pos: u16s(&dir.join("q.pos")),
            byte0,
        };
        assert_eq!(q.codes.len(), *q.byte0.last().unwrap());
        q
    }
    pub fn len(&self) -> usize {
        self.fam.len()
    }
    pub fn row(&self, i: usize) -> &[u8] {
        &self.codes[self.byte0[i]..self.byte0[i + 1]]
    }
    pub fn rowpos(&self, i: usize) -> &[u16] {
        &self.pos[self.byte0[i]..self.byte0[i + 1]]
    }
}

/// Units of each family, sorted: those with a row in any arm given.
pub fn family_units(fams: &[Family], arms: &[&Rows]) -> Vec<Vec<u32>> {
    (0..fams.len())
        .map(|f| {
            let mut u: Vec<u32> =
                arms.iter().flat_map(|r| r.unit[r.row0[f]..r.row0[f + 1]].iter().copied()).collect();
            u.sort_unstable();
            u.dedup();
            u
        })
        .collect()
}

/// Run `f(i)` for i in 0..n on `threads` threads, in any order.
/// `f(0..n)` in order, computed on `threads` threads.
pub fn par_map<T: Send>(n: usize, threads: usize, f: impl Fn(usize) -> T + Sync) -> Vec<T> {
    let chunk = n.div_ceil(threads.max(1) * 8).max(1);
    let parts: Vec<std::sync::Mutex<Vec<T>>> = (0..n.div_ceil(chunk)).map(|_| Default::default()).collect();
    parallel(parts.len(), threads, |p| *parts[p].lock().unwrap() = (p * chunk..((p + 1) * chunk).min(n)).map(&f).collect());
    parts.into_iter().flat_map(|x| x.into_inner().unwrap()).collect()
}

pub fn parallel(n: usize, threads: usize, f: impl Fn(usize) + Sync) {
    let next = std::sync::atomic::AtomicUsize::new(0);
    std::thread::scope(|s| {
        for _ in 0..threads {
            s.spawn(|| loop {
                let i = next.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                if i >= n {
                    break;
                }
                f(i)
            });
        }
    });
}
