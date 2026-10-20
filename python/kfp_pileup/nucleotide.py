"""Amino-acid alignments brought back to nucleotides (plan, phase 13, step 6b).

DIAMOND aligns in amino-acid space; each read is put back into nucleotide space along its
alignment, unchanged. Reference column c is nucleotide columns 3c, 3c+1, 3c+2 of a virtual
coding sequence (its bases unknown); each aligned residue contributes its read codon, read
from the FASTQ (recruited reads only, one pass) at the read-position key the crude pass
records. A deleted residue is three deleted bases, an inserted residue three inserted bases
before the next codon column, a frameshift a 1-nt indel: ``/`` (a base missing from the read)
deletes the next column's first base, ``\\`` inserts the skipped read base before it. Weights
are steps 8-9's (``em_posterior``).

Outputs:

- ``pileup_nt.parquet``: (seq_id, col, allele, weight), allele a base, ``-`` (deleted) or
  ``+`` and the inserted bases (before ``col``); a diagnostic.
- ``consensus.fna``: per unit, the weighted consensus base per column of its best-supported
  sequence (``N`` where no read covers it), named by its ``seq_id``: the protein arm's
  per-sample nucleotide reference for step 6c. ``consensus.gff`` (the whole sequence as one
  CDS, for ``bcftools csq``), ``consensus.genes.fna`` (the same as prodigal genes, for
  inStrain) and ``consensus.tsv`` (sequence, unit, covered columns, depth).
- ``recruited_1.fq.gz``, ``recruited_2.fq.gz``: both mates of every pair with a target
  alignment, for re-mapping to the consensus (``kfp-pileup handoff``).

ponytail: the consensus takes substitutions only, so its columns stay 3c..3c+2 and its frame
the protein's: indels, frameshifts included, are left to step 6c's callers, which see them
against it. Indels in the diagnostic pileup are placed as DIAMOND put them (to within a codon),
not left-normalised: the callers re-align. Fold majority indels into the consensus if reads
fail to map across them.
"""

import gzip
from collections.abc import Iterator, Sequence
from pathlib import Path
from typing import Final

import numpy as np
import polars as pl

from kfp_pileup.spill import Spill
from kmer_functional_profiler.reference import BASES, CODE_11

MIN_WEIGHT: Final = 1e-3
_COMP: Final = np.frombuffer(bytes.maketrans(b"ACGTacgt", b"TGCAtgca"), np.uint8)
_CODE: Final = np.frombuffer(CODE_11.encode(), np.uint8)
_BASE_INDEX: Final = np.full(256, 4, np.int64)
for _i, _b in enumerate(BASES):
    _BASE_INDEX[ord(_b)] = _BASE_INDEX[ord(_b.lower())] = _i


def _fastq(path: str | Path) -> Iterator[tuple[str, str, str]]:
    """(header line, sequence, quality) of a FASTQ, gzipped or not."""
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as f:
        while header := f.readline().rstrip("\n"):
            seq = f.readline().rstrip("\n")
            f.readline()
            yield header, seq, f.readline().rstrip("\n")


def _recruit(
    reads: Sequence[str | Path], names: set[str], out: Path
) -> tuple[pl.DataFrame, np.ndarray]:
    """Both mates of every pair whose name (either mate's first word) is in ``names``,
    written to ``recruited_<mate>.fq.gz``; returns (mate, read_name, off, length) and the
    reads' bases, concatenated."""
    rows, chunks, off = [], [], 0
    with (
        gzip.open(out / "recruited_1.fq.gz", "wt", compresslevel=1) as f1,
        gzip.open(out / "recruited_2.fq.gz", "wt", compresslevel=1) as f2,
    ):
        files = (f1, f2)
        for pair in zip(*(_fastq(r) for r in reads), strict=True):
            ids = [h[1:].split(maxsplit=1)[0] for h, _, _ in pair]
            if not names.intersection(ids):
                continue
            for mate, ((header, seq, qual), name, f) in enumerate(
                zip(pair, ids, files, strict=True)
            ):
                f.write(f"{header}\n{seq}\n+\n{qual}\n")
                rows.append((mate, name, off, len(seq)))
                chunks.append(seq.upper())
                off += len(seq)
    table = pl.DataFrame(
        rows,
        schema={"mate": pl.UInt8, "read_name": pl.String, "off": pl.Int64, "rlen": pl.Int64},
        orient="row",
    )
    return table, np.frombuffer("".join(chunks).encode(), np.uint8)


def _bases(buf: np.ndarray, off: np.ndarray, key: np.ndarray, k: np.ndarray) -> np.ndarray:
    """Base k (0, 1, 2 along the frame) from read position ``key``; complemented on the
    reverse strand."""
    rev = (key & 1).astype(bool)
    pos = (key >> 1) + np.where(rev, -k, k)
    b = buf[off + pos]
    out: np.ndarray = np.where(rev, _COMP[b], b)
    return out


def nucleotide_pileup(
    hits: pl.DataFrame,
    spill: Spill,
    refs: pl.DataFrame,
    reads: Sequence[str | Path],
    out: Path,
) -> dict[str, int]:
    """Step 6b from the crude pass's target columns and events (``spill``, read one
    partition at a time) and steps 8-9's weights (``hits``: assignments with
    ``em_posterior``). Writes the module's outputs to ``out``; returns counts,
    ``codon_mismatch`` the read codons whose translation is not the residue the alignment
    shows (0 unless coordinates are wrong; codons with N excepted)."""
    w = hits.filter(pl.col("role") == "target", pl.col("em_posterior") > MIN_WEIGHT).select(
        "hit", "read_name", "mate", "seq", weight="em_posterior"
    )
    table, buf = _recruit(reads, set(w["read_name"]), out)
    w = w.join(table, on=["mate", "read_name"]).select("hit", "seq", "weight", "off")
    piles, codon_mismatch, frameshifts = [], 0, 0
    for _, cols, events in spill.parts():  # all of a sequence's rows are in one partition
        pile, mismatched, shifted = _pile_part(cols.drop("seq"), events.drop("seq"), w, buf)
        piles.append(pile)
        codon_mismatch += mismatched
        frameshifts += shifted
    pile = pl.concat(piles)
    ids = refs.select("seq", "seq_id", "unit", "cluster_rep", "protein_id", "length")
    pile.join(ids.select("seq", "seq_id"), on="seq").drop("seq").sort(
        "seq_id", "col", "allele"
    ).write_parquet(out / "pileup_nt.parquet")

    # consensus: per unit, its best-supported sequence; bases only (see the module note)
    bases = pile.filter(pl.col("allele").is_in(list("ACGT")))
    best = (
        bases.group_by("seq")
        .agg(support=pl.col("weight").sum())
        .join(ids, on="seq")
        .sort("support", "seq", descending=[True, False])
        .unique("unit", keep="first")
        .sort("seq")
    )
    top = (
        bases.sort("weight", "allele", descending=[True, False])  # ties: the first allele
        .group_by("seq", "col")
        .agg(allele=pl.col("allele").first(), depth=pl.col("weight").sum())
    )
    top_by_seq = top.partition_by("seq", as_dict=True)
    none = top.clear()
    fna, gff, genes, info = [], ["##gff-version 3\n"], [], []
    for row in best.iter_rows(named=True):
        length = 3 * row["length"]
        cons = np.full(length, ord("N"), np.uint8)
        at = top_by_seq.get((row["seq"],), none)
        cons[at["col"].to_numpy()] = np.frombuffer("".join(at["allele"]).encode(), np.uint8)
        text, sid = cons.tobytes().decode(), row["seq_id"]
        fna.append(f">{sid} unit={row['unit']} cluster_rep={row['cluster_rep']}\n{text}\n")
        genes.append(f">{sid}_1 # 1 # {length} # 1 # ID={sid}_1\n{text}\n")
        gff.extend(
            f"{sid}\tkfp\t{kind}\t1\t{length}\t.\t+\t{phase}\t{attrs}\n"
            for kind, phase, attrs in (
                ("gene", ".", f"ID=gene:{sid};Name={sid};biotype=protein_coding"),
                (
                    "mRNA",
                    ".",
                    f"ID=transcript:{sid};Parent=gene:{sid};Name={sid};biotype=protein_coding",
                ),
                ("exon", ".", f"Parent=transcript:{sid}"),
                ("CDS", "0", f"Parent=transcript:{sid}"),
            )
        )
        info.append(
            (sid, row["unit"], row["cluster_rep"], row["protein_id"], length, at.height,
             float(at["depth"].sum()) / length)
        )  # fmt: skip
    (out / "consensus.fna").write_text("".join(fna))
    (out / "consensus.genes.fna").write_text("".join(genes))
    (out / "consensus.gff").write_text("".join(gff))
    pl.DataFrame(
        info,
        schema=["seq_id", "unit", "cluster_rep", "protein_id", "length", "covered", "depth"],
        orient="row",
    ).write_csv(out / "consensus.tsv", separator="\t")
    return {
        "recruited_mates": table.height,
        "consensus_sequences": best.height,
        "frameshift_alignments": frameshifts,
        "codon_mismatch": codon_mismatch,
    }


def _pile_part(
    cols: pl.DataFrame, events: pl.DataFrame, w: pl.DataFrame, buf: np.ndarray
) -> tuple[pl.DataFrame, int, int]:
    """One partition's nucleotide pileup (seq, col, allele, weight), its codon mismatches
    and its alignments with a frameshift."""
    # codons of aligned residues, and the translation check
    codons = cols.filter(pl.col("qpos") >= 0).join(w, on="hit")
    off, key = codons["off"].to_numpy(), codons["qpos"].to_numpy()
    b = np.stack([_bases(buf, off, key, np.full(len(key), k)) for k in range(3)])
    idx = _BASE_INDEX[b]
    ok = (idx < 4).all(0)
    aa = _CODE[16 * idx[0] + 4 * idx[1] + idx[2] & 63]
    residue = codons["residue"].to_numpy() & 0xDF  # upper case
    codon_mismatch = int((ok & (aa != residue)).sum())

    # one row per (alignment, nucleotide column): read bases, deletions as "-"
    c = codons["col"].to_numpy().astype(np.int64)
    base_rows = pl.DataFrame(
        {
            "hit": np.tile(codons["hit"].to_numpy(), 3),
            "col": np.concatenate([3 * c + k for k in range(3)]),
            "allele": np.frombuffer(b.tobytes(), "S1").astype("U1"),
        },
        schema={"hit": pl.UInt32, "col": pl.Int64, "allele": pl.String},
    )
    deleted = cols.filter(pl.col("qpos") < 0).join(w.select("hit"), on="hit")
    del_rows = pl.concat(
        [
            deleted.select("hit", col=3 * pl.col("col").cast(pl.Int64) + k, allele=pl.lit("-"))
            for k in range(3)
        ]
    )
    ev = events.join(w, on="hit")
    shift_del = ev.filter(pl.col("kind") == ord("/")).select(
        "hit", col=3 * pl.col("col").cast(pl.Int64), gone=pl.lit(True)
    )
    rows = (
        pl.concat([base_rows, del_rows])
        .join(shift_del, on=["hit", "col"], how="left")
        .with_columns(allele=pl.when(pl.col("gone")).then(pl.lit("-")).otherwise("allele"))
        .drop("gone")
    )
    # insertions: a residue's codon, or the base a frameshift skipped
    ins = ev.filter(pl.col("kind") != ord("/"))
    ioff, ikey = ins["off"].to_numpy(), ins["qpos"].to_numpy()
    n_bases = np.where(ins["kind"].to_numpy() == ord("I"), 3, 1)
    ib = [_bases(buf, ioff, ikey, np.full(len(ikey), k)) for k in range(3)]
    inserted = [
        "+" + bytes(int(ib[k][i]) for k in range(n_bases[i])).decode() for i in range(ins.height)
    ]
    ins_rows = ins.select(
        "hit", col=3 * pl.col("col").cast(pl.Int64), allele=pl.Series(inserted, dtype=pl.String)
    )
    pile = (
        pl.concat([rows, ins_rows])
        .join(w.select("hit", "seq", "weight"), on="hit")
        .group_by("seq", "col", "allele")
        .agg(pl.col("weight").sum())
    )
    frameshifts = int(ev.filter(pl.col("kind") != ord("I"))["hit"].n_unique())
    return pile, codon_mismatch, frameshifts
