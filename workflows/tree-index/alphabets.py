"""Tree index: reduced amino-acid alphabets for distant hits, before any build
(docs/tree-index-plan.md, phase 6b).

Under spike.py's BLOSUM62 strains (a substituted residue's replacement weighted 2^(S/2)), per
alphabet:

- ``within``: the share of substitutions that stay within a class;
- ``bits``: chance-match bits per residue, -log2 sum_c Q_c^2 (Q_c the class's frequency in
  ``--members``), so two random k-mers match with probability 2^-(k bits): what decoy seeds
  and background hits scale with;
- ``p_<id>``: a residue's chance to keep its class at identity id, id + (1 - id) within;
- ``cost_0.7``: -log2 p_0.7 / bits, survival lost at 70% identity per chance bit; the
  ranking at matched bits (lower is better);
- per bit budget B, k = round(B / bits) (capped at the alphabet's max k), the k-mer's survival
  p^k and the predicted miss rate of a 187-aa protein at 1x (~63% of positions covered) with
  one indexed k-mer per window of 8: a protein is missed when no indexed k-mer survives
  whole in a read, ``(1 - c p^k)^n``, n = 2 (L - k + 1) / (w + 1), independent k-mers.

``opt<n>``: an n-class alphabet that maximises log p_0.7 / bits (survival at 70% identity
per chance bit, so at matched bits), from greedy merges then single-residue moves.

Writes ``alphabets.tsv``; ``--print-opt`` prints the opt alphabets' classes.
"""

import argparse
import math
from pathlib import Path
from typing import Final

import numpy as np
import polars as pl

AMINO: Final = "ACDEFGHIKLMNPQRSTVWY"

# crates/core/src/alphabet.rs's Murphy10 and Dayhoff, then its REDUCED in the same order
BUILTIN: Final = {
    "murphy10": "LVIM C A G ST P FYW EDNQ KR H",  # Murphy, Wallqvist and Levy 2000
    "dayhoff": "AGPST C DENQ HKR ILMV FWY",
}
REDUCED: Final = {
    "murphy8": "LVIMC AG ST P FYW EDNQ KR H",
    "murphy12": "LVIM C A G ST P FY W EQ DN KR H",
    "murphy15": "LVIM C A G S T P FY W E Q D N KR H",
    "diamond11": "KREDQN C G H ILV M F Y W P STA",  # Buchfink et al. 2015
    "seb6": "AST CP DEHKNQR FWY G ILMV",  # Edgar 2004, Table 1
    "seb8": "AST C DHN EKQR FWY G ILMV P",
    "seb10": "AST C DN EQ FY G HW ILMV KR P",
    "seb14": "A C D EQ FY G H IV KR LM N P ST W",
    "sev10": "AST C DEN FY G H ILMV KQR P W",
    "lia10": "AC DE FWY G HN IV KQR LM P ST",
    "lib10": "AST C DEQ FWY G HN IV KR LM P",
    # optimise() on phase 0's quarter sample (2026-10-30)
    "opt4": "CFILMVWY DEHKNQR PST AG",
    "opt6": "FHWY ILMV ACT DEQ GNS KPR",
    "opt8": "FHWY DEQ KNR AC GS IV LM PT",
    "opt10": "FHWY AC DN EQ IV KR LM PT G S",
    "opt12": "FWY AC EQ HN IV KR LM D G P S T",
    "opt14": "FWY AC EQ HN LM D G I K P R S T V",
}
CLASSES: Final = {n: tuple(s.split()) for n, s in (BUILTIN | REDUCED).items()}
IDENTITIES: Final = (0.85, 0.8, 0.75, 0.7)
L, COVER, WINDOW = 187, 1 - math.exp(-1), 8


def substitution() -> np.ndarray:
    """P(b | a substituted), 20 x 20, spike.py's BLOSUM62 weights."""
    from Bio.Align import substitution_matrices

    m = substitution_matrices.load("BLOSUM62")
    s = lambda a, b: m[a][b]  # type: ignore[index]  # noqa: E731
    w = np.array([[0.0 if a == b else 2.0 ** (s(a, b) / 2) for b in AMINO] for a in AMINO])
    return w / w.sum(1, keepdims=True)


def composition(members: Path) -> np.ndarray:
    seqs = pl.read_parquet(members, columns=["sequence"])["sequence"]
    counts = np.bincount(np.frombuffer("".join(seqs).encode(), np.uint8), minlength=256)
    q = counts[np.frombuffer(AMINO.encode(), np.uint8)].astype(float)
    return q / q.sum()


def label(classes: list[int]) -> str:
    """Classes as space-separated residue groups, largest first."""
    groups: dict[int, str] = {}
    for aa, c in zip(AMINO, classes, strict=True):
        groups[c] = groups.get(c, "") + aa
    return " ".join(sorted(groups.values(), key=lambda g: (-len(g), g)))


def parse(spec: str) -> list[int]:
    cls = {aa: i for i, g in enumerate(spec.split()) for aa in g}
    assert sorted(cls) == sorted(AMINO), spec
    return [cls[aa] for aa in AMINO]


def stats(classes: list[int], q: np.ndarray, sub: np.ndarray) -> tuple[float, float]:
    """(within, bits) of an alphabet given as each residue's class."""
    c = np.array(classes)
    same = c[:, None] == c[None, :]
    within = float((q * (sub * same).sum(1)).sum())
    qc = np.bincount(c, weights=q)
    return within, -math.log2(float((qc**2).sum()))


def objective(classes: list[int], q: np.ndarray, sub: np.ndarray, ident: float = 0.7) -> float:
    within, bits = stats(classes, q, sub)
    return math.log(ident + (1 - ident) * within) / bits


def optimise(n: int, q: np.ndarray, sub: np.ndarray) -> list[int]:
    """Greedy merges from 20 classes to n (best objective each step), then single-residue
    moves while any improves it. ponytail: local search, an exhaustive or annealed search if
    the opt alphabets look stuck."""
    c = list(range(20))
    while len(set(c)) > n:
        ids = sorted(set(c))
        merged = lambda p: [p[0] if x == p[1] else x for x in c]  # noqa: B023, E731
        pairs = [(i, j) for i in ids for j in ids if i < j]
        c = merged(max(pairs, key=lambda p: objective(merged(p), q, sub)))
    improved = True
    while improved:
        improved = False
        for r in range(20):
            for t in sorted(set(c)):
                trial = c.copy()
                trial[r] = t
                if len(set(trial)) == n and objective(trial, q, sub) > objective(c, q, sub) + 1e-12:
                    c, improved = trial, True
    return c


def main() -> None:
    p = argparse.ArgumentParser(description=(__doc__ or "").split("\n\n")[0])
    p.add_argument("--members", type=Path, required=True, help="sequences for the composition")
    p.add_argument("--opt", default="4,5,6,7,8,9,10,11,12,13,14", help="opt alphabet sizes")
    p.add_argument("--budgets", default="15,20,25", help="chance bits per k-mer")
    p.add_argument("--out", type=Path, default=Path("alphabets.tsv"))
    p.add_argument("--print-opt", action="store_true")
    a = p.parse_args()
    q, sub = composition(a.members), substitution()
    alphabets = {"protein": list(range(20))}
    alphabets |= {n: parse(s) for n, s in (BUILTIN | REDUCED).items() if not n.startswith("opt")}
    for n in map(int, a.opt.split(",")):
        got = alphabets[f"opt{n}"] = optimise(n, q, sub)
        if f"opt{n}" in REDUCED and label(got) != label(parse(REDUCED[f"opt{n}"])):
            print(f"opt{n} differs from REDUCED's (another composition): {label(got)}")
    rows = []
    for name, c in alphabets.items():
        within, bits = stats(c, q, sub)
        n_cls = len(set(c))
        row = {"alphabet": name, "classes": n_cls, "within": within, "bits": bits}
        row |= {f"p_{i}": i + (1 - i) * within for i in IDENTITIES}
        row["cost_0.7"] = -math.log2(row["p_0.7"]) / bits
        max_k = 64 // max(math.ceil(math.log2(n_cls)), 1)
        for b in map(int, a.budgets.split(",")):
            k = min(max(round(b / bits), 1), max_k)
            row[f"k@{b}"], row[f"bits@{b}"] = k, k * bits
            n_idx = 2 * (L - k + 1) / (WINDOW + 1)
            for i in IDENTITIES:
                s = (i + (1 - i) * within) ** k
                row[f"miss_{i}@{b}"] = (1 - COVER * s) ** n_idx
        row["layout"] = label(c)
        rows.append(row)
    df = pl.DataFrame(rows).sort("cost_0.7")
    df.write_csv(a.out, separator="\t", float_precision=4)
    with pl.Config(tbl_rows=-1, tbl_cols=-1, tbl_width_chars=250, float_precision=3):
        shown = [c for c in df.columns if c.startswith(("k@", "miss_0.7@", "miss_0.8@"))]
        print(df.select("alphabet", "classes", "within", "bits", "p_0.7", "cost_0.7", *shown))
    if a.print_opt:
        for name, c in alphabets.items():
            if name.startswith("opt"):
                print(f"{name}\t{label(c)}")


if __name__ == "__main__":
    main()
