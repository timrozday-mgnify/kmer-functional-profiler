"""Per-sample reference from a profile (plan, phase 13, step 2).

Three roles, one FASTA:

- ``target``: units the profile calls present (``present_prob`` >= ``min_present``);
- ``neighbour``: the other hit units in the targets' components (``component``), decoys
  that catch reads from a sibling the profile did not call; at most ``neighbours`` per
  target in a component, those with the most tier-2 hits (a real sample's giant component
  holds most hit units: millions of sequences and a REFERENCE killed at 256 GB);
- ``decoy``: the targets' sequences reversed, a null for the score threshold.

``withhold`` units (``cluster_rep`` names) are left out of every role, as if the profile had
never seen them: the benchmark's withheld sibling, whose reads the decoys must catch.

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
NEIGHBOURS: Final = 5  # per target in a component
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
    withhold: Sequence[str] = (),
    neighbours: int = NEIGHBOURS,
) -> dict[str, int]:
    """Write ``ref.fa``, ``refs.tsv`` and ``molecule.txt`` to ``out_dir``; return counts."""
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    hit = (
        pl.read_csv(
            profile,
            separator="\t",
            columns=["unit", "cluster_rep", "present_prob", "component", "hits"],
            schema_overrides={"cluster_rep": pl.String},
        )
        .filter(~pl.col("cluster_rep").is_in(list(withhold)))
        .drop("cluster_rep")
        .with_columns(pl.col("unit").cast(pl.UInt32))
    )
    targets = hit.filter(pl.col("present_prob") >= min_present)
    if targets.height == 0:
        raise ValueError(f"no units with present_prob >= {min_present} in {profile}")
    # ponytail: ranked by hits, a proxy for the reads a sibling could take; rank by k-mers
    # shared with a target (the hit table's links) if called siblings slip past the cap
    per_component = targets.group_by("component").agg(cap=pl.len() * neighbours)
    others = (
        hit.join(per_component, on="component")
        .filter(~pl.col("unit").is_in(targets["unit"].implode()))
        .sort(["hits", "unit"], descending=[True, False])
    )
    others = others.filter(pl.int_range(pl.len()).over("component") < pl.col("cap"))
    roles = pl.concat(
        [
            targets.select("unit", "component", role=pl.lit("target")),
            others.select("unit", "component", role=pl.lit("neighbour")),
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
