"""Read routing ahead of DIAMOND (plan, phase 13, step 3b's "Afterwards"; 2026-10-21).

Each read pair is searched against the part of the per-sample reference it shares k-mers
with, not the whole of it:

- the reference's components (``refs.tsv``) are packed, in order, into blocks of about
  ``block_letters`` residues; a component is never split, so a block holds whole components;
- every k-mer of the reference (the index's alphabet; its k unless ``k`` is given;
  unsampled) gives the blocks that hold it; the reads are hashed in all six frames (a
  frameshifted read keeps its k-mers on either side of the shift), and a pair goes to every
  block it shares at least ``min_kmers`` distinct k-mers with. A smaller ``k`` routes more
  distant relatives (more siblings, more homologues); a larger ``min_kmers`` fewer;
- pairs with no hit are dropped: they could align only as distant homologues, which the EM
  sends to the novel class (64% of aligned pairs on ``pileup-cost``).

DIAMOND then runs once per block on that block's pairs, so a read's search grows with its
block, not with the reference. ``letters.txt`` holds the whole reference's residues, passed
as DIAMOND's ``--dbsize`` so that e-values do not change with the split.

Writes ``b<i>.fa``, ``b<i>_1.fq.gz`` and ``b<i>_2.fq.gz`` per block with pairs, ``letters.txt``
and ``route.json`` (counts).
"""

import gzip
import json
import resource
from collections.abc import Iterator, Sequence
from contextlib import ExitStack
from itertools import islice
from pathlib import Path
from typing import Final

import numpy as np
import polars as pl

from kfp_pileup.reference import REF_FASTA, REFS
from kmer_functional_profiler import _core
from kmer_functional_profiler.index import read_fasta

BLOCK_LETTERS: Final = 2_000_000  # residues per block (hadza10M's reference: ~10 blocks)
BATCH_READS: Final = 100_000
MIN_KMERS: Final = 1  # distinct k-mers a pair shares with a block to be sent there
LETTERS: Final = "letters.txt"
ROUTE: Final = "route.json"


def blocks(refs: pl.DataFrame, block_letters: int) -> np.ndarray:
    """Each reference sequence's block: components in order of appearance, a new block
    started where the residues so far pass a multiple of ``block_letters``."""
    comp = refs.group_by("component", maintain_order=True).agg(letters=pl.col("length").sum())
    start = comp["letters"].cum_sum() - comp["letters"]
    comp = comp.with_columns(
        block=(start // block_letters).rank("dense").cast(pl.Int64) - 1
    )  # dense: a component above block_letters skips no number
    return refs.join(comp, on="component", how="left", maintain_order="left")["block"].to_numpy()


def kmer_table(
    seqs: Sequence[bytes], block: np.ndarray, k: int, alphabet: str
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Every k-mer's blocks, as sorted distinct hashes, offsets and block ids (CSR)."""
    got = _core.hash_proteins(seqs, k, alphabet=alphabet)
    table = (
        pl.DataFrame({"hash": got["hash"], "block": block[got["seq"]]})
        .unique()
        .sort("hash", "block")
    )
    hashes, first = np.unique(table["hash"].to_numpy(), return_index=True)
    return hashes, np.append(first, table.height), table["block"].to_numpy()


def lookup(
    table: tuple[np.ndarray, np.ndarray, np.ndarray], hashes: np.ndarray, reads: np.ndarray
) -> pl.DataFrame:
    """(read, block, kmers) rows for a batch: the distinct k-mers a read shares with each
    block it hits."""
    keys, offsets, values = table
    at = np.minimum(np.searchsorted(keys, hashes), max(len(keys) - 1, 0))
    found = keys[at] == hashes if len(keys) else np.zeros(len(hashes), dtype=bool)
    at, reads, hashes = at[found], reads[found], hashes[found]
    n = offsets[at + 1] - offsets[at]
    within = np.arange(n.sum()) - np.repeat(np.cumsum(n) - n, n)
    return (
        pl.DataFrame(
            {
                "read": np.repeat(reads, n),
                "block": values[np.repeat(offsets[at], n) + within],
                "hash": np.repeat(hashes, n),
            }
        )
        .group_by("read", "block")
        .agg(kmers=pl.col("hash").n_unique().cast(pl.UInt32))
    )


def pair_hits(
    table: tuple[np.ndarray, np.ndarray, np.ndarray],
    reads: Sequence[str | Path],
    k: int,
    alphabet: str,
    genetic_code: int = 11,
    batch_reads: int = BATCH_READS,
) -> tuple[pl.DataFrame, int]:
    """(read, block, kmers) for every pair that hits the table (read ids number pairs in
    input order and never span batches), and the number of pairs."""
    hits = _core.FastxHits(
        *map(str, reads),
        k=k,
        alphabet=alphabet,
        genetic_code=genetic_code,
        frames="all",
        batch_reads=batch_reads,
    )
    rows = pl.concat(
        [pl.DataFrame(schema={"read": pl.UInt64, "block": pl.Int64, "kmers": pl.UInt32})]
        + [lookup(table, batch["hash"], batch["read"]) for batch in hits]
    )
    return rows, hits.n_reads


def _records(path: Path) -> Iterator[bytes]:
    """FASTQ records, four lines each, gzipped or plain."""
    with open(path, "rb") as f:
        gz = f.read(2) == b"\x1f\x8b"
    with gzip.open(path, "rb") if gz else open(path, "rb") as f:
        while record := b"".join(islice(f, 4)):
            yield record


def route(
    ref_dir: str | Path,
    index_dir: str | Path,
    reads: Sequence[str | Path],
    out_dir: str | Path,
    block_letters: int = BLOCK_LETTERS,
    k: int = 0,
    min_kmers: int = MIN_KMERS,
    genetic_code: int = 11,
    batch_reads: int = BATCH_READS,
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
    seqs = [sequences[s].encode() for s in refs["seq_id"]]
    block = blocks(refs, block_letters)
    for b in np.unique(block):
        with open(out / f"b{b}.fa", "w") as fasta:
            fasta.writelines(
                f">{s}\n{sequences[s]}\n" for s in refs.filter(pl.Series(block == b))["seq_id"]
            )
    (out / LETTERS).write_text(f"{refs['length'].sum()}\n")

    k = k or params["k"]
    table = kmer_table(seqs, block, k, params["alphabet"])
    rows, pairs = pair_hits(table, reads, k, params["alphabet"], genetic_code, batch_reads)
    routed = rows.filter(pl.col("kmers") >= min_kmers).sort("read", "block")

    # one handle per block with pairs, per mate; a large reference has hundreds of blocks
    soft, hard = resource.getrlimit(resource.RLIMIT_NOFILE)
    resource.setrlimit(resource.RLIMIT_NOFILE, (max(soft, min(hard, 4_096)), hard))
    read_ids, block_ids = routed["read"].to_numpy(), routed["block"].to_numpy()
    for mate, path in enumerate(reads, 1):
        with ExitStack() as stack:
            files: dict[int, gzip.GzipFile] = {}
            j = 0
            for i, record in enumerate(_records(Path(path))):
                while j < len(read_ids) and read_ids[j] == i:
                    b = int(block_ids[j])
                    if b not in files:
                        name = out / f"b{b}_{mate}.fq.gz"
                        files[b] = stack.enter_context(gzip.open(name, "wb", compresslevel=1))
                    files[b].write(record)
                    j += 1
    counts = {
        "k": k,
        "min_kmers": min_kmers,
        "pairs": pairs,
        "routed_pairs": routed["read"].n_unique(),
        "pair_blocks": routed.height,
        "blocks": int(block.max()) + 1,
        "blocks_with_pairs": routed["block"].n_unique(),
        "letters": int(refs["length"].sum()),
        "kmers": len(table[0]),
    }
    counts["dropped_pairs"] = pairs - counts["routed_pairs"]
    (out / ROUTE).write_text(json.dumps(counts, indent=1) + "\n")
    return counts
