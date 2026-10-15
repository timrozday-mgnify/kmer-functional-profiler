"""Per-sample reference from a profile (plan, phase 13, step 2).

Three roles, one FASTA:

- ``target``: units the profile calls present (``present_prob`` >= ``min_present``);
- ``neighbour``: the other hit units in the targets' components (``component``), decoys
  that catch reads from a sibling the profile did not call;
- ``decoy``: the targets' sequences reversed, a null for the score threshold.

Sequence ids are ``s0``, ``s1``, …, safe as DIAMOND subject and SAM reference names;
``refs.tsv`` maps them to role, unit, ``cluster_rep``, member and component.
"""

import json
from collections.abc import Sequence
from pathlib import Path
from typing import Final

import polars as pl

from kmer_functional_profiler.index import molecule
from kmer_functional_profiler.store import MEMBERS, fetch, write_fasta

MIN_PRESENT: Final = 0.5
REF_FASTA: Final = "ref.fa"
REFS: Final = "refs.tsv"
MOLECULE: Final = "molecule.txt"  # nucleotide or protein: which aligner the pipeline runs


def build_reference(
    profile: str | Path,
    index_dir: str | Path,
    out_dir: str | Path,
    min_present: float = MIN_PRESENT,
    members: int = MEMBERS,
    store: Sequence[str | Path] = (),
) -> dict[str, int]:
    """Write ``ref.fa``, ``refs.tsv`` and ``molecule.txt`` to ``out_dir``; return counts."""
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    hit = pl.read_csv(
        profile, separator="\t", columns=["unit", "present_prob", "component"]
    ).with_columns(pl.col("unit").cast(pl.UInt32))
    targets = hit.filter(pl.col("present_prob") >= min_present)
    if targets.height == 0:
        raise ValueError(f"no units with present_prob >= {min_present} in {profile}")
    neighbours = hit.filter(
        pl.col("component").is_in(targets["component"].implode()),
        ~pl.col("unit").is_in(targets["unit"].implode()),
    )
    roles = pl.concat(
        [
            targets.select("unit", "component", role=pl.lit("target")),
            neighbours.select("unit", "component", role=pl.lit("neighbour")),
        ]
    )
    seqs = fetch(index_dir, roles["unit"].to_list(), members, store).join(roles, on="unit")
    decoys = seqs.filter(pl.col("role") == "target").with_columns(
        pl.col("sequence").str.reverse(), role=pl.lit("decoy")
    )
    refs = pl.concat([seqs, decoys]).with_columns(
        seq_id=pl.format("s{}", pl.int_range(pl.len())),
        length=pl.col("sequence").str.len_bytes(),
    )
    write_fasta(refs, out / REF_FASTA, ids="seq_id")
    refs.drop("sequence").write_csv(out / REFS, separator="\t")
    meta = json.loads((Path(index_dir) / "meta.json").read_text())
    (out / MOLECULE).write_text(molecule(meta["params"]["alphabet"]) + "\n")
    counts = refs.group_by("role").len()
    return dict(zip(counts["role"], counts["len"], strict=True))
