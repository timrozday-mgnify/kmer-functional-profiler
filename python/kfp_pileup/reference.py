"""Per-sample reference from a profile (plan, phase 13, step 2).

Two roles, one FASTA:

- ``target``: units the profile calls present at a false discovery rate: the top units by
  ``present_prob`` whose expected false share is at most ``max_fdr`` (``present_q`` <=
  ``max_fdr``; q computed from ``present_prob`` for profiles without ``present_q``). No
  aai or ani is needed; a profile queried with ``--min-identity`` carries its identity test
  in ``present_prob`` (plan, phase 13, 2026-10-20);
- ``neighbour``: the other hit units in the targets' components (``component``), decoys
  that catch reads from a sibling the profile did not call; at most ``neighbours`` per
  target in a component, those with the most tier-1 hits (a real sample's giant component
  holds most hit units: millions of sequences and a REFERENCE killed at 256 GB).

Reversed decoys were dropped (plan, phase 13, 2026-10-19): on ``spike-protein`` they kept
0-21 mates per sample off targets, and the EM's novel class takes chance alignments.

``withhold`` units (``cluster_rep`` names) are left out of every role, as if the profile had
never seen them: the benchmark's withheld sibling, whose reads neighbours and the EM's
novel class must catch.

Sequence ids are ``s0``, ``s1``, …, safe as DIAMOND subject and SAM reference names;
``refs.tsv`` maps them to role, unit, ``cluster_rep``, member and component.
"""

import json
from collections.abc import Sequence
from pathlib import Path
from typing import Final

import polars as pl

from kmer_functional_profiler.index import molecule
from kmer_functional_profiler.query import fdr_q
from kmer_functional_profiler.store import MEMBERS, fetch, write_fasta

MAX_FDR: Final = 0.05  # targets: expected false share among them (the user's choice)
NEIGHBOURS: Final = 5  # per target in a component
REF_FASTA: Final = "ref.fa"
REFS: Final = "refs.tsv"
MOLECULE: Final = "molecule.txt"  # nucleotide or protein: which aligner the pipeline runs


def build_reference(
    profile: str | Path,
    index_dir: str | Path,
    out_dir: str | Path,
    max_fdr: float = MAX_FDR,
    members: int = MEMBERS,
    store: Sequence[str | Path] = (),
    withhold: Sequence[str] = (),
    neighbours: int = NEIGHBOURS,
) -> dict[str, int]:
    """Write ``ref.fa``, ``refs.tsv`` and ``molecule.txt`` to ``out_dir``; return counts."""
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    with open(profile) as f:
        has_q = "present_q" in f.readline().rstrip("\n").split("\t")
    hit = pl.read_csv(
        profile,
        separator="\t",
        columns=["unit", "cluster_rep", "present_prob", "component", "hits"]
        + (["present_q"] if has_q else []),
        schema_overrides={"cluster_rep": pl.String},
    )
    if not has_q:  # a profile from before present_q: the same q over the same rows
        hit = hit.with_columns(present_q=pl.Series(fdr_q(hit["present_prob"].to_numpy())))
    hit = (
        hit.filter(~pl.col("cluster_rep").is_in(list(withhold)))
        .drop("cluster_rep")
        .with_columns(pl.col("unit").cast(pl.UInt32))
    )
    targets = hit.filter(pl.col("present_q") <= max_fdr)
    if targets.height == 0:
        raise ValueError(f"no units with present_q <= {max_fdr} in {profile}")
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
    refs = seqs.with_columns(
        seq_id=pl.format("s{}", pl.int_range(pl.len())),
        length=pl.col("sequence").str.len_bytes(),
    )
    write_fasta(refs, out / REF_FASTA, ids="seq_id")
    refs.drop("sequence").write_csv(out / REFS, separator="\t")
    meta = json.loads((Path(index_dir) / "meta.json").read_text())
    (out / MOLECULE).write_text(molecule(meta["params"]["alphabet"]) + "\n")
    counts = refs.group_by("role").len()
    return dict(zip(counts["role"], counts["len"], strict=True))
