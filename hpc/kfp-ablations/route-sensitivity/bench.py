"""Routing's cost on a laptop-sized input (plan, phase 13, routing): time and peak memory.

``make OUT``: a reference of MGnify proteins (``--proteins``, the first ``--residues``;
components of 9 sequences) and ``--pairs`` simulated 150 bp pairs, half from reference
proteins back-translated at ~5% nucleotide divergence, half random DNA. It measures cost,
not recall: real reads hit more of the reference.

``run OUT``: ``kfp-pileup route``'s stages (table build, hashing and lookup, total with the
FASTQ writing) with their wall time and the process's peak RSS. Laptop-safe at 2 M pairs
against 20 M residues (cap threads with RAYON_NUM_THREADS / POLARS_MAX_THREADS).

    python bench.py make mid --proteins data/mgnify/mgy_clusters_head.faa --pairs 200000
    python bench.py run mid
"""

import argparse
import gzip
import json
import random
import resource
import sys
import tempfile
import time
from pathlib import Path

import polars as pl

from kfp_pileup.route import PART_KMERS, blocks, kmer_table, pair_hits, route
from kmer_functional_profiler.index import read_fasta

CODE = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"  # TCAG order
READ, INSERT = 150, 300


def make(out: Path, proteins: Path, residues: int, pairs: int, seed: int = 1) -> None:
    rng = random.Random(seed)
    seqs, total = [], 0
    for _, seq in read_fasta(proteins):
        seqs.append(seq)
        total += len(seq)
        if total >= residues:
            break
    (out / "ref").mkdir(parents=True, exist_ok=True)
    (out / "idx").mkdir(exist_ok=True)
    with open(out / "ref" / "ref.fa", "w") as fa, open(out / "ref" / "refs.tsv", "w") as tsv:
        tsv.write("seq_id\trole\tunit\tcomponent\tlength\n")
        for i, s in enumerate(seqs):
            fa.write(f">s{i}\n{s}\n")
            tsv.write(f"s{i}\ttarget\t{i // 3}\t{i // 9}\t{len(s)}\n")
    (out / "idx" / "meta.json").write_text(json.dumps({"params": {"k": 11, "alphabet": "protein"}}))
    back: dict[str, list[str]] = {}
    for i, aa in enumerate(CODE):
        back.setdefault(aa, []).append("TCAG"[i // 16] + "TCAG"[i // 4 % 4] + "TCAG"[i % 4])
    long_enough = [s for s in seqs if 3 * len(s) >= INSERT]
    comp, qual = str.maketrans("ACGT", "TGCA"), "I" * READ
    with (
        gzip.open(out / "r1.fq.gz", "wt", compresslevel=1) as f1,
        gzip.open(out / "r2.fq.gz", "wt", compresslevel=1) as f2,
    ):
        for n in range(pairs):
            if n % 2 == 0:
                p = rng.choice(long_enough)
                start = rng.randrange(len(p) - INSERT // 3 + 1)
                cds = "".join(rng.choice(back.get(a, back["A"])) for a in p[start:][: INSERT // 3])
                frag = "".join(rng.choice("ACGT") if rng.random() < 0.05 else c for c in cds)
            else:
                frag = "".join(rng.choice("ACGT") for _ in range(INSERT))
            f1.write(f"@p{n}/1\n{frag[:READ]}\n+\n{qual}\n")
            f2.write(f"@p{n}/2\n{frag[-READ:].translate(comp)[::-1]}\n+\n{qual}\n")
    print(f"{len(seqs)} sequences, {total} residues, {pairs} pairs")


def run(out: Path, k: int, part_kmers: int) -> None:
    start = last = time.perf_counter()

    def lap(name: str) -> None:
        nonlocal last
        now = time.perf_counter()
        peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        peak /= 1e9 if sys.platform == "darwin" else 1e6  # bytes on macOS, KiB on Linux
        print(f"{name:14s} {now - last:7.1f} s   peak RSS {peak:.1f} GB", flush=True)
        last = now

    refs = pl.read_csv(out / "ref" / "refs.tsv", separator="\t")
    sequences = dict(read_fasta(out / "ref" / "ref.fa"))
    seqs = [sequences[s].encode() for s in refs["seq_id"]]
    block = blocks(refs, 2_000_000)
    reads = (out / "r1.fq.gz", out / "r2.fq.gz")
    lap("read reference")
    with tempfile.TemporaryDirectory(dir=out) as scratch:
        table = kmer_table(seqs, block, k, "protein", Path(scratch), part_kmers)
        lap("table")
        rows, _, pairs = pair_hits(table, reads, k, "protein")
        lap("hash + lookup")
        del table
    print(f"  {pairs} pairs, {rows['read'].n_unique()} with a hit, {rows.height} rows")
    counts = route(out / "ref", out / "idx", reads, out / "route", k=k, part_kmers=part_kmers)
    lap("route() total")
    print(f"  {json.dumps(counts)}\n  all stages {time.perf_counter() - start:.1f} s")


def main() -> None:
    parser = argparse.ArgumentParser(description="Routing's cost on a laptop-sized input")
    sub = parser.add_subparsers(dest="command", required=True)
    m = sub.add_parser("make")
    m.add_argument("out", type=Path)
    m.add_argument("--proteins", type=Path, required=True)
    m.add_argument("--residues", type=int, default=20_000_000)
    m.add_argument("--pairs", type=int, default=200_000)
    r = sub.add_parser("run")
    r.add_argument("out", type=Path)
    r.add_argument("--k", type=int, default=11)
    r.add_argument("--part-kmers", type=int, default=PART_KMERS)
    args = parser.parse_args()
    if args.command == "make":
        make(args.out, args.proteins, args.residues, args.pairs)
    else:
        run(args.out, args.k, args.part_kmers)


if __name__ == "__main__":
    main()
