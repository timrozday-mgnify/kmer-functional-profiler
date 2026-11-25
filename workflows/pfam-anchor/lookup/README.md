# Residue scan against MGnify90

A prototype of the anchored lookup: score a read's residues against every MGnify90 representative of its Pfam family, column by
column, weighted by information bits; keep the best rows; score those by columns, then with a banded alignment.

- `rs/`: the prototype (Rust, no dependencies). `cargo test --release` runs its tests; `src/main.rs` lists the commands.
- `pack.py`: domains and reads on HMM columns, as the flat files `rs/` reads.
- `large.py`: one very large family at full release size, driven by `hpc/kfp-ablations/pfam-anchor-large/run.sh`.
