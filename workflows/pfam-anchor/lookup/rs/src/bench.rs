//! What a scan costs when it is laid out for it: one thread, a family of N rows held by column, K columns
//! of the read compared, weights added (bits in eighths, one byte), in tiles of 32,768 rows. Plain Rust, no intrinsics.
pub fn run(n: usize) {
    let mut seed = 1u64;
    let mut rnd = || {
        seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        (seed >> 59) as u8
    };
    for k in [8usize, 16, 40] {
        let cols: Vec<Vec<u8>> = (0..k).map(|_| (0..n).map(|_| rnd()).collect()).collect();
        let query: Vec<u8> = (0..k).map(|_| rnd()).collect();
        let reads = 10;
        let mut kept = 0usize;
        let t = std::time::Instant::now();
        for r in 0..reads {
            let mut acc = [0u8; 32768];
            for lo in (0..n).step_by(32768) {
                let hi = (lo + 32768).min(n);
                let a = &mut acc[..hi - lo];
                a.fill(0);
                for (c, &q) in cols.iter().zip(&query) {
                    let (q, w) = (q ^ (r as u8 & 1), 24u8);
                    for (x, &v) in a.iter_mut().zip(&c[lo..hi]) {
                        *x = x.saturating_add(if v == q { w } else { 0 });
                    }
                }
                kept += a.iter().filter(|&&x| x >= 120).count();
            }
        }
        let s = t.elapsed().as_secs_f64() / reads as f64;
        println!("{n} rows x {k} columns: {:.3} ms a read, {:.3} ns a row, {:.1} GB/s ({kept})", s * 1e3, s * 1e9 / n as f64, (n * k) as f64 / s / 1e9);
    }
}
