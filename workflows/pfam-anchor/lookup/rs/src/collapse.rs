//! Greedy collapse of the representatives' GA rows within a family: a row joins the first kept row that matches at least
//! T of its residues, column by column. Longest first. collapse.tsv: t, family, rows, kept, residues, residues_kept;
//! collapse.<t>.units: the units with a kept row (u32), by family.
use crate::data::*;
use std::fmt::Write as _;
use std::path::Path;
use std::sync::Mutex;

pub fn run(dir: &Path, threads: usize) {
    let fams = families(dir);
    let rows = Rows::load(dir, "reps", &fams);
    let ga = bytes(&dir.join("reps.ga"));
    let out = Mutex::new(String::from("t\tfamily\trows\tkept\tresidues\tresidues_kept\n"));
    for t in [0.9f32, 0.8, 0.7, 0.6, 0.5] {
        let units: Vec<Mutex<Vec<u32>>> = fams.iter().map(|_| Mutex::new(vec![])).collect();
        parallel(fams.len(), threads, |i| {
            let len = fams[i].len;
            let mut v: Vec<(usize, usize)> = (rows.row0[i]..rows.row0[i + 1])
                .filter(|&r| ga[r] == 1)
                .map(|r| {
                    let b = rows.byte0[i] + (r - rows.row0[i]) * len;
                    (rows.codes[b..b + len].iter().filter(|&&a| a < 20).count(), r)
                })
                .collect();
            v.sort_unstable_by(|a, b| b.cmp(a));
            let row = |r: usize| {
                let b = rows.byte0[i] + (r - rows.row0[i]) * len;
                &rows.codes[b..b + len]
            };
            let mut kept: Vec<usize> = vec![];
            let (mut res, mut res_kept) = (0, 0);
            for &(n, r) in &v {
                let x = row(r);
                let need = (t * n as f32).ceil() as usize;
                let joins = kept.iter().any(|&k| row(k).iter().zip(x).map(|(&a, &b)| (a == b && b < 20) as u16).sum::<u16>() as usize >= need);
                res += n;
                if !joins {
                    kept.push(r);
                    res_kept += n;
                }
            }
            writeln!(out.lock().unwrap(), "{t}\tPF{:05}\t{}\t{}\t{res}\t{res_kept}", fams[i].acc, v.len(), kept.len()).unwrap();
            let mut u: Vec<u32> = kept.iter().map(|&r| rows.unit[r]).collect();
            u.sort_unstable();
            u.dedup();
            *units[i].lock().unwrap() = u;
        });
        let mut o = vec![];
        for u in &units {
            let u = u.lock().unwrap();
            o.extend_from_slice(&(u.len() as u32).to_le_bytes());
            u.iter().for_each(|x| o.extend_from_slice(&x.to_le_bytes()));
        }
        std::fs::write(dir.join(format!("collapse.{t}.units")), o).unwrap();
        eprintln!("{t} done");
    }
    std::fs::write(dir.join("collapse.tsv"), &*out.lock().unwrap()).unwrap();
}
