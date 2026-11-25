"""Per-model settings of kfpfam's prefilter, from the thresholds: how its top mode keeps every
model's weakest passing domains. Runs in the kfpfam image (numpy).

- ``bars``: thresholds.tsv -> ``ga.tsv`` (model, seq, dom), for ``filter seeds``, and
  ``gapped.txt``, a gapped bar per model for the filter's FILTER_GAPPED_FILE: a short repeat's
  domain scores under the one bar for the rest. With ``--slack``, every model's bar is that many
  bits under its weakest passing domain in a protein, where that is over ``--flat``.
- ``choose``: ``filter seeds``' table -> a seed threshold T per model, for FILTER_T_FILE: the
  highest at which a 50-residue stretch that just passes holds no word of the index with
  probability under exp(-theta).

From choose_seeds.py and recall_score.py of the B11 arms (kfpfam's recall modes, on reads).
"""

import argparse
import csv
import math
from pathlib import Path

import numpy as np

# The words of one stretch overlap, so they fail together: the ratio for this shape, over Pfam 38.2
OVERLAP = {"10110001": 0.936, "1100101": 0.935}
NT, LOWEST, STEP = 25, 4.5, 0.25  # `filter seeds`' grid of T
SLACK = 4.0  # bits: HMMER sums over alignments, and the filter's is one of them
# `filter seeds` puts the whole threshold in 50 residues. That is right for a read, and for Pfam's
# thresholds (about 25 bits) it worked on proteins, but KOfam's hundreds of bits are spread over a
# long model, so a model's threshold counts up to this. Lower only adds seeds.
MOST = 30.0


def threshold(seq: float, dom: float) -> float:
    """The score in bits a model's weakest passing domain is taken to have."""
    finite = [x for x in (dom, seq) if math.isfinite(x)]
    return finite[0] if finite else MOST  # none: it shows only with --weak-evalue


def weakest_raw(length: int, bits: float, orf: int = 50) -> float:
    """Raw bits of a domain of that score: plus what entering the model and the ORF cost it."""
    return (
        bits
        + math.log2(length * (length + 1) / 2)
        + 2 * math.log2((orf + 3) / 2)
        - math.log2(orf + 1)
    )


def rows(path: Path) -> list[dict[str, str]]:
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def bars(a: argparse.Namespace) -> None:
    ga, gapped = [], []
    for m, t in zip(rows(a.models), rows(a.thresholds), strict=True):
        bits = threshold(float(t["seq"]), float(t["dom"]))
        ga.append(f"{m['acc']}\t{min(bits, MOST)}\t{min(bits, MOST)}\n")
        bar = weakest_raw(int(m["length"]), min(bits, MOST)) - SLACK
        if a.slack:
            bar = max(weakest_raw(int(m["length"]), bits, 100) - a.slack, min(a.flat, bar))
        gapped.append(f"{bar:.3f}\n")
    a.ga.write_text("".join(ga))
    a.gapped.write_text("".join(gapped))


def choose(a: argparse.Namespace) -> None:
    quality = np.fromfile(a.table, "<f4").reshape(-1, 2, NT)[:, 1] * OVERLAP[a.shape]
    quality = np.maximum.accumulate(quality[:, ::-1], 1)[:, ::-1]  # it falls as T rises
    grid = LOWEST + STEP * np.arange(NT)
    # a model that misses theta even at the lowest T gets that
    t = np.array([np.interp(-a.theta, -q, grid) for q in quality])
    np.savetxt(a.out, t, fmt="%.3f")
    print(f"theta {a.theta}: T 5/50/95% {np.quantile(t, [0.05, 0.5, 0.95]).round(2)}")


def main() -> None:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    sub = p.add_subparsers(required=True)
    b = sub.add_parser("bars")
    b.set_defaults(run=bars)
    for name in ("models", "thresholds", "ga", "gapped"):
        b.add_argument(f"--{name}", type=Path, required=True)
    b.add_argument("--slack", type=float, default=0, help="bits; 0: the one bar, but for repeats")
    b.add_argument("--flat", type=float, default=21.25, help="--slack's floor")
    c = sub.add_parser("choose")
    c.set_defaults(run=choose)
    c.add_argument("--table", type=Path, required=True)
    c.add_argument("--shape", choices=sorted(OVERLAP), required=True)
    c.add_argument("--theta", type=float, required=True)
    c.add_argument("--out", type=Path, required=True)
    a = p.parse_args()
    a.run(a)


if __name__ == "__main__":
    main()
