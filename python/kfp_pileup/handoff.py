"""The hand-off to existing variant and haplotype callers (plan, phase 13, steps 6b-6c).

Reads mapped by bwa-mem2 (to ``consensus.fna`` for amino-acid targets, to ``ref.fa`` for
nucleotide ones) are kept where step 8 assigns their pair to the unit of the sequence they
map to: per mate, its best alignment (``AS``) to a target sequence of the pair's most
responsible unit, if that unit's responsibility is at least ``min_posterior`` (the
novel-relative share excluded, as ``em_posterior`` is). The kept record becomes primary, its
MAPQ the posterior's Phred (capped at 60), and mate fields are set from the kept mate (or
marked unmapped when the mate is not kept). Output: SAM, name-grouped, for ``samtools sort``.
"""

import math
import re
from collections.abc import Iterator
from pathlib import Path
from typing import Final, TextIO

import polars as pl

from kfp_pileup.reference import REFS

MIN_POSTERIOR: Final = 0.5
MAX_MAPQ: Final = 60
_CIGAR: Final = re.compile(r"(\d+)([MIDNSHP=X])")
_MATE_SUFFIX: Final = re.compile(r"/[12]$")
_COMPLEMENT: Final = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def _ref_span(cigar: str) -> int:
    return sum(int(n) for n, op in _CIGAR.findall(cigar) if op in "MDN=X")


def _records(path: Path) -> Iterator[list[str] | str]:
    """Header lines as strings, then alignment records as field lists."""
    with open(path) as f:
        for line in f:
            yield line if line.startswith("@") else line.rstrip("\n").split("\t")


def assigned_units(pileup_dir: Path) -> dict[str, tuple[int, float]]:
    """Per read pair (mate suffixes stripped): its most responsible target unit and that
    responsibility, from ``assignments.parquet``."""
    a = pl.read_parquet(pileup_dir / "assignments.parquet")
    per_mate = (
        a.filter(pl.col("role") == "target")
        .with_columns(pair=pl.col("read_name").str.replace(_MATE_SUFFIX.pattern, ""))
        .group_by("pair", "mate", "unit")
        .agg(p=pl.col("em_posterior").fill_null(0.0).sum())
        .group_by("pair", "unit")
        .agg(pl.col("p").max())
        .sort("p", descending=True)
        .unique("pair", keep="first")
    )
    return {pair: (unit, p) for pair, unit, p in per_mate.select("pair", "unit", "p").iter_rows()}


def handoff(
    pileup_dir: str | Path,
    ref_dir: str | Path,
    sam: str | Path,
    out: str | Path,
    min_posterior: float = MIN_POSTERIOR,
) -> dict[str, int]:
    """Filter ``sam`` (bwa-mem2, name-grouped) to step 8's assignments; write ``out``."""
    units = assigned_units(Path(pileup_dir))
    refs = pl.read_csv(Path(ref_dir) / REFS, separator="\t", infer_schema=False)
    unit_of = {
        s: int(u) for s, u, r in refs.select("seq_id", "unit", "role").iter_rows() if r == "target"
    }
    kept = pairs = 0
    with open(out, "w") as f:
        group: list[list[str]] = []
        for rec in _records(Path(sam)):
            if isinstance(rec, str):
                f.write(rec)
                continue
            if group and rec[0] != group[0][0]:
                kept += _write_pair(f, group, units, unit_of, min_posterior)
                pairs += 1
                group = []
            group.append(rec)
        if group:
            kept += _write_pair(f, group, units, unit_of, min_posterior)
            pairs += 1
    return {"pairs": pairs, "kept_mates": kept}


def _write_pair(
    f: TextIO,
    group: list[list[str]],
    units: dict[str, tuple[int, float]],
    unit_of: dict[str, int],
    min_posterior: float,
) -> int:
    unit, p = units.get(_MATE_SUFFIX.sub("", group[0][0]), (-1, 0.0))
    if p < min_posterior:
        return 0
    # each mate's primary sequence and quality, to fill records that omit them
    primary = {}
    for r in group:
        if not int(r[1]) & 0x900 and r[9] != "*":
            primary[int(r[1]) & 0xC0] = (r[9], r[10], bool(int(r[1]) & 16))
    best: dict[int, tuple[int, list[str]]] = {}
    for r in group:
        flag = int(r[1])
        if flag & 4 or unit_of.get(r[2]) != unit:
            continue
        mate = flag & 0xC0
        score = next((int(t[5:]) for t in r[11:] if t.startswith("AS:i:")), 0)
        if mate not in best or score > best[mate][0]:
            best[mate] = (score, r)
    chosen = {m: list(r) for m, (_, r) in best.items()}
    mapq = str(min(MAX_MAPQ, round(-10 * math.log10(max(1 - p, 1e-6)))))
    for m, r in chosen.items():
        flag = int(r[1]) & ~0x900 & ~0x2 & ~0x8 & ~0x20
        reverse = bool(flag & 16)
        if r[9] == "*" and m in primary:
            seq, qual, strand = primary[m]
            if strand != reverse:
                seq, qual = seq.translate(_COMPLEMENT)[::-1], qual[::-1]
            r[9], r[10] = seq, qual
        other = chosen.get(m ^ 0xC0) if m else None
        if other is not None:
            same = other[2] == r[2]
            flag |= 0x20 if int(other[1]) & 16 else 0
            flag |= 0x2 if same else 0
            r[6], r[7] = ("=" if same else other[2]), other[3]
            if same:
                start = min(int(r[3]), int(other[3]))
                end = max(int(r[3]) + _ref_span(r[5]), int(other[3]) + _ref_span(other[5]))
                r[8] = str((end - start) * (1 if int(r[3]) <= int(other[3]) else -1))
            else:
                r[8] = "0"
        else:
            flag |= 0x8 if flag & 1 else 0
            r[6], r[7], r[8] = "=", r[3], "0"
        r[1], r[4] = str(flag), mapq
        r[11:] = [t for t in r[11:] if not t.startswith(("SA:Z:", "XA:Z:"))]
        f.write("\t".join(r) + "\n")
    return len(chosen)
