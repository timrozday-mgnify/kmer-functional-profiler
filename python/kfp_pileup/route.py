"""Read routing ahead of DIAMOND (plan, phase 13, step 3b's "Afterwards"; 2026-10-21).

Each read pair is searched against the part of the per-sample reference it shares k-mers
with, not the whole of it:

- the reference's components (``refs.tsv``) are packed, in order, into blocks of about
  ``block_letters`` residues; a component is never split, so a block holds whole components;
- every k-mer of the reference (the index's alphabet; its k unless ``k`` is given;
  unsampled; ``route.fa``'s sequences when ``reference`` wrote more per unit for routing,
  ``--route-members``) gives the blocks that hold it; the reads are hashed in all six frames (a
  frameshifted read keeps its k-mers on either side of the shift), and a pair goes to every
  block it shares at least ``min_kmers`` distinct k-mers with. A smaller ``k`` routes more
  distant relatives (more siblings, more homologues); a larger ``min_kmers`` fewer;
- pairs with no hit are dropped: they could align only as distant homologues, which the EM
  sends to the novel class (64% of aligned pairs on ``pileup-cost``).

DIAMOND then runs once per block on that block's pairs, so a read's search grows with its
block, not with the reference. ``letters.txt`` holds the whole reference's residues, passed
as DIAMOND's ``--dbsize`` so that e-values do not change with the split.

Writes per block ``b<i>.fa`` (its reference sequences) and, where it has pairs,
``b<i>_1.fa.gz`` and ``b<i>_2.fa.gz`` (its pairs' mates as FASTA, names cut at the first
word: DIAMOND reads neither qualities nor comments, which are about half of a FASTQ.gz);
``letters.txt``; ``route.json`` (counts); and ``routes.parquet``, every (pair, block) hit
with its distinct k-mers (``read_name``, ``pair`` in input order, ``block``, ``kmers``),
routed or not, so that another ``min_kmers`` is a filter, not a rerun.
"""

import gzip
import itertools
import json
import resource
import tempfile
from collections.abc import Iterator, Sequence
from contextlib import ExitStack
from itertools import islice
from pathlib import Path
from typing import Final

import numpy as np
import polars as pl

from kfp_pileup.reference import REF_FASTA, REFS, ROUTE_FASTA, ROUTE_SEQS
from kmer_functional_profiler import _core
from kmer_functional_profiler.index import (
    PackedPart,
    PackedTable,
    key_shift,
    packed_layout,
    read_fasta,
)

BLOCK_LETTERS: Final = 2_000_000  # residues per block (hadza10M's reference: ~10 blocks)
BATCH_READS: Final = 100_000
PART_KMERS: Final = 25_000_000  # reference residues per hash range of the table build
FP_BITS: Final = 32  # fingerprint bits: false hits ~keys per bucket x 2^-32 per lookup
HASH_MAX: Final = 2**64 - 1
# route-sensitivity (2026-10-21): the cheapest setting keeping >= 0.99 of the pairs the
# unrouted EM gives a unit, in every sample (spike-protein x 10, hadza10M). min_kmers applies
# per block, whose size is fixed, so chance hits per block do not grow with the reference
ROUTE_K: Final = 8
MIN_KMERS: Final = 3  # distinct k-mers a pair shares with a block to be sent there
LETTERS: Final = "letters.txt"
ROUTE: Final = "route.json"
ROUTES: Final = "routes.parquet"
NAME_CHUNK: Final = 1_000_000


def blocks(refs: pl.DataFrame, block_letters: int) -> np.ndarray:
    """Each reference sequence's block: components in order of appearance, a new block
    started where the residues so far pass a multiple of ``block_letters``."""
    comp = refs.group_by("component", maintain_order=True).agg(letters=pl.col("length").sum())
    start = comp["letters"].cum_sum() - comp["letters"]
    comp = comp.with_columns(
        block=(start // block_letters).rank("dense").cast(pl.Int64) - 1
    )  # dense: a component above block_letters skips no number
    return refs.join(comp, on="component", how="left", maintain_order="left")["block"].to_numpy()


def routing_sequences(
    ref: Path, refs: pl.DataFrame, sequences: dict[str, str], block: np.ndarray
) -> tuple[list[bytes], np.ndarray]:
    """The sequences routing hashes and each one's block: ``route.fa`` (``route_members``
    per unit, from ``reference``) if there is one, else the reference's own, with blocks
    by component, as the reference's."""
    if not (ref / ROUTE_SEQS).exists():
        return [sequences[s].encode() for s in refs["seq_id"]], block
    of = dict(zip(refs["component"].to_list(), block.tolist(), strict=False))
    more = pl.read_csv(ref / ROUTE_SEQS, separator="\t", infer_schema_length=None)
    found = dict(read_fasta(ref / ROUTE_FASTA))
    return (
        [found[s].encode() for s in more["seq_id"]],
        np.array([of[c] for c in more["component"].to_list()], dtype=np.int64),
    )


def kmer_table(
    seqs: Sequence[bytes],
    block: np.ndarray,
    k: int,
    alphabet: str,
    scratch: Path,
    part_kmers: int = PART_KMERS,
) -> PackedTable:
    """Every reference k-mer's blocks, as a ``PackedTable`` (values: block ids; ~8 bytes per
    k-mer). Built one hash range at a time, each range's k-mers hashed on their own and its
    part saved to ``scratch``, then joined memory-mapped: the peak is one range's build
    (~170 bytes per k-mer), not the whole reference's."""
    n = max(sum(len(s) for s in seqs), 1)  # at least the reference's k-mers
    layout = packed_layout(HASH_MAX, FP_BITS, n)
    shift = key_shift(layout)
    parts = -(-n // part_kmers)
    # range edges on key boundaries, so parts hold disjoint, increasing keys
    edges = [(i * (HASH_MAX + 1) // parts) >> shift << shift for i in range(parts)]
    edges.append(HASH_MAX + 1)
    built = []
    for i, (lo, hi) in enumerate(itertools.pairwise(edges)):
        got = _core.hash_proteins(seqs, k, alphabet=alphabet, min_hash=lo, max_hash=hi - 1)
        part = PackedPart.build(got["hash"], block.astype(np.uint64)[got["seq"]], layout)
        del got
        part.save(str(scratch / f"part{i}"))
        built.append(PackedPart.load(str(scratch / f"part{i}")))
    return PackedTable.concat(layout, built)


def pair_hits(
    table: PackedTable,
    reads: Sequence[str | Path],
    k: int,
    alphabet: str,
    genetic_code: int = 11,
    batch_reads: int = BATCH_READS,
) -> tuple[pl.DataFrame, int]:
    """(read, block, kmers) for every pair that hits the table, and the number of pairs.
    Read ids number pairs in input order and never span batches; the lookup and the
    distinct-k-mer count per pair and block run in Rust, on every core (``route_hits``)."""
    hits = _core.FastxHits(
        *map(str, reads),
        k=k,
        alphabet=alphabet,
        genetic_code=genetic_code,
        frames="all",
        batch_reads=batch_reads,
    )
    u64 = np.uint64
    rows = pl.concat(
        [pl.DataFrame(schema={"read": pl.UInt64, "block": pl.UInt32, "kmers": pl.UInt32})]
        + [
            pl.DataFrame(
                _core.route_hits(
                    table, b["read"].astype(u64, copy=False), b["hash"].astype(u64, copy=False)
                )
            )
            for b in hits
        ]
    )
    return rows, hits.n_reads


def _records(path: Path) -> Iterator[tuple[bytes, bytes]]:
    """FASTQ records as (name, sequence): the name's first word, without ``@``."""
    with open(path, "rb") as f:
        gz = f.read(2) == b"\x1f\x8b"
    with gzip.open(path, "rb") if gz else open(path, "rb") as f:
        while record := list(islice(f, 4)):
            yield record[0][1:].split(maxsplit=1)[0], record[1].rstrip(b"\r\n")


def route(
    ref_dir: str | Path,
    index_dir: str | Path,
    reads: Sequence[str | Path],
    out_dir: str | Path,
    block_letters: int = BLOCK_LETTERS,
    k: int = ROUTE_K,
    min_kmers: int = MIN_KMERS,
    genetic_code: int = 11,
    batch_reads: int = BATCH_READS,
    part_kmers: int = PART_KMERS,
) -> dict[str, int]:
    """Split the reference in ``ref_dir`` into blocks and the pairs in ``reads`` (two
    FASTQs, mate order) by the blocks they share at least ``min_kmers`` distinct k-mers
    with (``k``: 0 for the index's); write them to ``out_dir``; return counts."""
    ref, out = Path(ref_dir), Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    params = json.loads((Path(index_dir) / "meta.json").read_text())["params"]
    if params["alphabet"] == "dna":
        raise ValueError("routing is for amino-acid references (DIAMOND); bwa-mem2 is not routed")
    refs = pl.read_csv(ref / REFS, separator="\t", infer_schema_length=None)
    sequences = dict(read_fasta(ref / REF_FASTA))
    block = blocks(refs, block_letters)
    seqs, seq_block = routing_sequences(ref, refs, sequences, block)
    for b in np.unique(block):
        with open(out / f"b{b}.fa", "w") as fasta:
            fasta.writelines(
                f">{s}\n{sequences[s]}\n" for s in refs.filter(pl.Series(block == b))["seq_id"]
            )
    (out / LETTERS).write_text(f"{refs['length'].sum()}\n")

    k = k or params["k"]
    with tempfile.TemporaryDirectory(prefix=".table-", dir=out) as scratch:
        table = kmer_table(seqs, seq_block, k, params["alphabet"], Path(scratch), part_kmers)
        rows, pairs = pair_hits(table, reads, k, params["alphabet"], genetic_code, batch_reads)
        n_kmers = len(table.fingerprints)
        del table  # its parts are memory-mapped from the scratch directory
    rows = rows.sort("read", "block")
    routed = rows.filter(pl.col("kmers") >= min_kmers)

    # one handle per block with pairs, per mate; a large reference has hundreds of blocks
    soft, hard = resource.getrlimit(resource.RLIMIT_NOFILE)
    resource.setrlimit(resource.RLIMIT_NOFILE, (max(soft, min(hard, 4_096)), hard))
    read_ids, block_ids = routed["read"].to_numpy(), routed["block"].to_numpy()
    hit_ids = rows["read"].unique(maintain_order=True).to_numpy()
    names: list[pl.Series] = []  # Arrow string chunks: a Python list of 10^8 names is ~6 GB
    chunk: list[str] = []
    for mate, path in enumerate(reads, 1):
        with ExitStack() as stack:
            files: dict[int, gzip.GzipFile] = {}
            j = h = 0
            for i, (name, seq) in enumerate(_records(Path(path))):
                if mate == 1 and h < len(hit_ids) and hit_ids[h] == i:
                    chunk.append(name.decode())
                    h += 1
                    if len(chunk) == NAME_CHUNK:
                        names.append(pl.Series(chunk, dtype=pl.String))
                        chunk = []
                while j < len(read_ids) and read_ids[j] == i:
                    b = int(block_ids[j])
                    if b not in files:
                        fa = out / f"b{b}_{mate}.fa.gz"
                        files[b] = stack.enter_context(gzip.open(fa, "wb", compresslevel=1))
                    files[b].write(b">" + name + b"\n" + seq + b"\n")
                    j += 1
    # every hit, routed or not, so another min_kmers needs no rerun
    names.append(pl.Series(chunk, dtype=pl.String))
    named = pl.DataFrame({"read": hit_ids, "read_name": pl.concat(names)})
    rows.join(named, on="read").select(
        "read_name", pair="read", block="block", kmers="kmers"
    ).write_parquet(out / ROUTES, compression="zstd")
    counts = {
        "k": k,
        "min_kmers": min_kmers,
        "pairs": pairs,
        "routed_pairs": routed["read"].n_unique(),
        "pair_blocks": routed.height,
        "blocks": int(block.max()) + 1,
        "blocks_with_pairs": routed["block"].n_unique(),
        "letters": int(refs["length"].sum()),
        "kmers": n_kmers,
    }
    counts["dropped_pairs"] = pairs - counts["routed_pairs"]
    (out / ROUTE).write_text(json.dumps(counts, indent=1) + "\n")
    return counts
