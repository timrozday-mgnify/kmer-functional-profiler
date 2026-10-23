"""Presence calibration at full-index scale (plan, phase 13 / phase 7: presence at a false
discovery rate): `--background-share` against spikes of the index's own units in real reads.

Subcommands (run by ``main.nf``; each prints what it wrote):

- ``pick``: ``--units`` units drawn uniformly from the index, keeping those whose
  representative is ``--min-len``-``--max-len`` residues; fetched with up to ``--members``
  sequences each. Writes ``targets.fa`` and ``groups.tsv`` (member, unit: ``spike.py
  simulate``'s inputs) and ``units.tsv`` (unit, ``rep_len``, ``n_members``, ``m_g`` and
  ``uniq`` = ``pin_sum`` / ``m_g``, the mean ``p_in`` of the unit's sampled k-mers: the share
  of its members holding each, 1 for a singleton; not whether other units hold them).
- ``reverse``: a null sample, the reads reversed (not complemented): same count, lengths and
  composition, no coding homology, so every unit it calls present is a false call.
- ``score``: one spiked sample's profile against its truth, beside the unspiked and null
  samples' profiles at the same share -> one row per q threshold (``Q``): ``called``,
  ``recall`` (spiked units called; also by depth, identity and ``uniq`` bins),
  ``new_false`` (units called here but not unspiked, not spiked: the spike's relatives, an
  upper bound), ``false_share`` = new_false / (spiked called + new_false), ``new_relative``
  (new_false in a component, of this profile, that holds a spiked unit: relatives the strain
  may be nearer than its own unit), ``false_share_unrelated`` (false_share without them: the
  q cut moving, chance), ``unspiked_called`` and ``null_called``. Before any presence call
  (the same in every q row): ``recall_hit1``, ``recall_hit2`` (spiked units with >= 1, >= 2
  hit k-mers, shared ones included: what a routing threshold can reach), ``recall_gathered``
  (gather gave them an own k-mer: what presence can reach), ``recall_hit1`` by depth and
  identity bins, and ``unspiked_hit1``, ``_hit2``, ``null_hit1``, ``_hit2``: those samples'
  units at the same thresholds (a routing threshold's reference size and chance hits). A
  profile lists every hit unit only with ``query --aai``; without, the hit counts are of
  gathered units.
- ``decide``: every score row -> ``summary.tsv`` (per share and q, the worst seed's
  false_share and the mean recalls) and ``decision.tsv``: the largest share passing the rule
  at q = ``--q`` in every seed: false_share <= ``--max-false`` and null_called <=
  ``--max-null`` x unspiked_called.
"""

import argparse
from pathlib import Path
from typing import Final

import numpy as np
import polars as pl

from kmer_functional_profiler.index import UnitTable
from kmer_functional_profiler.store import fetch

Q: Final = (0.01, 0.05, 0.1)
DEPTH_BINS: Final = [0.5, 1.0, 2.0, 5.0]  # x coverage of the CDS
UNIQ_BINS: Final = [0.25, 0.5, 0.75]


def pick(args: argparse.Namespace) -> None:
    table = UnitTable(args.index)
    rng = np.random.default_rng(args.seed)
    # ponytail: uniform over units, 3x oversampled for the length filter; stratify by
    # uniqueness if the uniform draw leaves a bin of units.tsv too thin to read
    ids = np.unique(rng.choice(table.height, min(3 * args.units, table.height), replace=False))
    seqs = fetch(args.index, ids.tolist(), args.members, args.store or ())
    reps = seqs.filter(pl.col("rep")).with_columns(rep_len=pl.col("sequence").str.len_bytes())
    keep = (
        reps.filter(pl.col("rep_len").is_between(args.min_len, args.max_len))
        .sample(fraction=1.0, shuffle=True, seed=args.seed)
        .head(args.units)
        .select("unit", "rep_len")
    )
    # members under --min-member-frac of the representative's length are dropped: strains of
    # fragments cover few of the unit's sampled k-mers (21% of spikes, plan 2026-10-23)
    seqs = seqs.join(keep, on="unit").filter(
        pl.col("sequence").str.len_bytes() >= args.min_member_frac * pl.col("rep_len")
    )
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    with open(out / "targets.fa", "w") as fa, open(out / "groups.tsv", "w") as groups:
        for i, (unit, seq) in enumerate(seqs.select("unit", "sequence").iter_rows()):
            fa.write(f">t{i}\n{seq}\n")  # simulate's member names; groups give the unit
            groups.write(f"t{i}\t{unit}\n")
    units = keep.sort("unit")["unit"].to_numpy()
    m_g = table["m_g"][units]
    pl.DataFrame(
        {
            "unit": units,
            "n_members": table["n_members"][units],
            "m_g": m_g,
            "uniq": table["pin_sum"][units] / np.maximum(m_g, 1),
        }
    ).join(keep, on="unit").write_csv(out / "units.tsv", separator="\t")
    print(f"{keep.height} units ({seqs.height} sequences) of {ids.size} drawn -> {out}")


def reverse(args: argparse.Namespace) -> None:
    import gzip

    for src, dst in zip(args.reads, args.out, strict=True):
        with gzip.open(src, "rt") as f, gzip.open(dst, "wt", compresslevel=1) as g:
            for n, line in enumerate(f):
                g.write(line[-2::-1] + "\n" if n % 4 in (1, 3) else line)
    print(f"reversed {args.reads} -> {args.out}")


def _called(profile: Path, q: float) -> pl.Series:
    return (
        pl.scan_csv(profile, separator="\t")
        .filter(pl.col("present_q") <= q)
        .select(pl.col("unit").cast(pl.Int64))
        .collect()["unit"]
    )


def _bin(col: str, edges: list[float]) -> pl.Expr:
    """Bins labelled by their lower edge: 0 (below the first edge), then each edge."""
    labels = ["0"] + [f"{e:g}" for e in edges]
    return pl.col(col).cut(edges, labels=labels, left_closed=True).cast(pl.String)


def score(args: argparse.Namespace) -> None:
    truth = (
        pl.read_csv(args.truth, separator="\t", columns=["unit", "identity_bin", "depth"])
        .with_columns(pl.col("unit").cast(pl.Int64))
        .join(pl.read_csv(args.units, separator="\t").select("unit", "uniq"), on="unit")
        .with_columns(
            depth_bin=_bin("depth", DEPTH_BINS),
            uniq_bin=_bin("uniq", UNIQ_BINS),
            identity_bin=pl.col("identity_bin").cast(pl.String),
        )
    )
    spiked = truth["unit"]
    components = (
        pl.scan_csv(args.profile, separator="\t")
        .select(pl.col("unit").cast(pl.Int64), "component", "kmers_hit", "kmers_unique")
        .collect()
    )
    evidence = truth.join(components, on="unit", how="left").with_columns(
        pl.col("kmers_hit", "kmers_unique").fill_null(0)
    )
    reach: dict[str, object] = {
        "recall_hit1": (evidence["kmers_hit"] >= 1).mean(),
        "recall_hit2": (evidence["kmers_hit"] >= 2).mean(),
        "recall_gathered": (evidence["kmers_unique"] >= 1).mean(),
    }
    for name, path in (("unspiked", args.unspiked), ("null", args.null)):
        hits = pl.scan_csv(path, separator="\t").select("kmers_hit").collect()["kmers_hit"]
        reach[f"{name}_hit1"], reach[f"{name}_hit2"] = (
            int((hits >= 1).sum()),
            int((hits >= 2).sum()),
        )
    for col in ("depth_bin", "identity_bin"):
        by = evidence.group_by(col).agg((pl.col("kmers_hit") >= 1).mean()).sort(col)
        for key, frac in by.iter_rows():
            reach[f"recall_hit1_{col.removesuffix('_bin')}_{key}"] = frac
    spiked_components = components.filter(pl.col("unit").is_in(spiked.implode()))["component"]
    relatives = components.filter(pl.col("component").is_in(spiked_components.implode()))["unit"]
    rows = []
    for q in Q:
        called = _called(args.profile, q)
        unspiked = _called(args.unspiked, q)
        hit = truth.with_columns(called=pl.col("unit").is_in(called.implode()))
        new_false = called.filter(
            ~called.is_in(unspiked.implode()) & ~called.is_in(spiked.implode())
        )
        relative = int(new_false.is_in(relatives.implode()).sum())
        unrelated = new_false.len() - relative
        tp = int(hit["called"].sum())
        row: dict[str, object] = {
            "sample": args.sample,
            "share": args.share,
            "q": q,
            "spiked": truth.height,
            "called": called.len(),
            "recall": tp / max(truth.height, 1),
            "new_false": new_false.len(),
            "false_share": new_false.len() / max(tp + new_false.len(), 1),
            "new_relative": relative,
            "false_share_unrelated": unrelated / max(tp + unrelated, 1),
            "unspiked_called": unspiked.len(),
            "null_called": _called(args.null, q).len(),
            **reach,
        }
        for col in ("depth_bin", "identity_bin", "uniq_bin"):
            for key, frac, n in (
                hit.group_by(col).agg(pl.col("called").mean(), pl.len()).sort(col).iter_rows()
            ):
                row[f"recall_{col.removesuffix('_bin')}_{key}"] = frac
                row[f"n_{col.removesuffix('_bin')}_{key}"] = n
        rows.append(row)
    pl.DataFrame(rows).write_csv(args.out, separator="\t")
    print(f"{args.sample} share {args.share} -> {args.out}")


def decide(args: argparse.Namespace) -> None:
    scores = pl.concat([pl.read_csv(p, separator="\t") for p in args.scores], how="diagonal")
    scores = scores.with_columns(
        null_ratio=pl.col("null_called") / pl.col("unspiked_called").clip(lower_bound=1),
        passes=(pl.col("false_share") <= args.max_false)
        & (pl.col("null_called") <= args.max_null * pl.col("unspiked_called")),
    )
    recalls = [c for c in scores.columns if c.startswith("recall")]
    summary = (
        scores.group_by("share", "q")
        .agg(
            pl.len().alias("seeds"),
            pl.col("false_share").max().alias("false_share_worst"),
            pl.col("new_false").mean(),
            pl.col("new_relative").mean(),
            pl.col("false_share_unrelated").max().alias("false_share_unrelated_worst"),
            pl.col("null_called").max().alias("null_called_worst"),
            pl.col("null_ratio").max().alias("null_ratio_worst"),
            pl.col("unspiked_called").first(),
            *[
                pl.col(c).first()
                for c in scores.columns
                if c.endswith(("_hit1", "_hit2")) and c not in recalls
            ],
            pl.col("passes").all(),
            *[pl.col(c).mean() for c in recalls],
        )
        .sort("q", "share")
    )
    summary.write_csv(args.out, separator="\t")
    scores.write_csv(args.out.with_name("scores.tsv"), separator="\t")
    at_q = summary.filter(pl.col("q") == args.q)
    passing = at_q.filter(pl.col("passes")).sort("share")
    chosen = passing.tail(1)
    decision = at_q.select(
        "share", "q", "passes", "false_share_worst", "null_ratio_worst", "unspiked_called", "recall"
    ).with_columns(chosen=pl.col("share").is_in(chosen["share"].implode()))
    decision.write_csv(args.decision, separator="\t")
    with pl.Config(tbl_cols=-1, tbl_rows=-1, tbl_width_chars=250):
        print(decision)
    print(f"chosen share: {chosen['share'].to_list() or 'none passes'}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = parser.add_subparsers(dest="command", required=True)
    p = sub.add_parser("pick")
    p.add_argument("--index", type=Path, required=True)
    p.add_argument("--store", type=Path, nargs="*")
    p.add_argument("--units", type=int, default=2000)
    p.add_argument("--members", type=int, default=3)
    p.add_argument("--min-len", type=int, default=200)
    p.add_argument("--max-len", type=int, default=600)
    p.add_argument("--min-member-frac", type=float, default=0.9)
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--out", type=Path, default=Path("."))
    p = sub.add_parser("reverse")
    p.add_argument("--reads", type=Path, nargs=2, required=True)
    p.add_argument("--out", type=Path, nargs=2, required=True)
    p = sub.add_parser("score")
    p.add_argument("--sample", required=True)
    p.add_argument("--share", required=True)
    p.add_argument("--truth", type=Path, required=True)
    p.add_argument("--units", type=Path, required=True)
    p.add_argument("--profile", type=Path, required=True)
    p.add_argument("--unspiked", type=Path, required=True)
    p.add_argument("--null", type=Path, required=True)
    p.add_argument("--out", type=Path, default=Path("scores.tsv"))
    p = sub.add_parser("decide")
    p.add_argument("scores", type=Path, nargs="+")
    p.add_argument("--q", type=float, default=0.05)
    p.add_argument("--max-false", type=float, default=0.05)
    p.add_argument("--max-null", type=float, default=0.01)
    p.add_argument("--out", type=Path, default=Path("summary.tsv"))
    p.add_argument("--decision", type=Path, default=Path("decision.tsv"))
    args = parser.parse_args()
    {"pick": pick, "reverse": reverse, "score": score, "decide": decide}[args.command](args)


if __name__ == "__main__":
    main()
