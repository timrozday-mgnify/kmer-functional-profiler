"""Which HMM columns tell MGnify90 units apart, and how to score and seed on them.

docs/pfam-anchor-plan.md, phase 1b. No read simulation: reads are windows of real proteins held
out of the reference, on the HMM columns of their full-length alignments.

- ``known``: a window of a unit's non-reference member; truth is its own unit.
- ``novel``: a window of a unit held out whole; truth is the remaining unit nearest over the
  whole domain, binned by that identity.

Each family is a matrix of units x HMM columns (one domain per unit: the representative's best
at GA, else the best member's), 0 where a unit has no residue. Column statistics come from the
reference only. Writes to ``--out``: ``scores.tsv``, ``seeds.tsv``, ``columns.tsv``,
``residues.tsv``, ``pairs.tsv``, ``reads.parquet`` and ``speed.json``.
"""

import argparse
import json
import re
import time
from pathlib import Path

import numpy as np
import polars as pl
from Bio.Align import substitution_matrices
from numpy.typing import NDArray

AMINO = "ACDEFGHIKLMNPQRSTVWY"
CODE = np.zeros(256, dtype=np.uint8)
CODE[np.frombuffer(AMINO.encode(), dtype=np.uint8)] = np.arange(1, 21)
OTHER = 21  # any other letter: never matches
# HMMER's amino acid background, order AMINO.
BACKGROUND = np.array(
    [
        0.0787945, 0.0151600, 0.0535222, 0.0668298, 0.0397062, 0.0695071, 0.0229198, 0.0590092,
        0.0594422, 0.0963728, 0.0237718, 0.0414386, 0.0482904, 0.0395639, 0.0540978, 0.0683364,
        0.0540687, 0.0673417, 0.0114135, 0.0304133,
    ]
)  # fmt: skip
OP = re.compile(r"(\d+)([MID])")
MIN_OVERLAP = 0.6  # of the window's columns, for a unit to be scored at all
THETA = 0.8  # chance a residue is unchanged, in the likelihood ratio


def column_residues(
    hmm_from: int, ali_from: int, cigar: str, sequence: str, length: int
) -> tuple[NDArray[np.uint8], NDArray[np.int32]]:
    """One domain on the HMM's columns (index = column, 1-based): residue codes (0: none) and each
    residue's position in the protein (-1: none). Inserted residues have no column."""
    codes = np.zeros(length + 1, dtype=np.uint8)
    positions = np.full(length + 1, -1, dtype=np.int32)
    residues = CODE[np.frombuffer(sequence.encode(), dtype=np.uint8)]
    residues[residues == 0] = OTHER
    column, position = hmm_from, ali_from - 1
    for n, op in OP.findall(cigar):
        n = int(n)
        if op == "M":
            codes[column : column + n] = residues[position : position + n]
            positions[column : column + n] = np.arange(position, position + n)
            column += n
            position += n
        elif op == "I":
            position += n
        else:
            column += n
    return codes, positions


def read_hmms(path: str) -> dict[int, NDArray[np.float64]]:
    """Per accession, match emissions: columns x 20, row 0 unused."""
    out: dict[int, NDArray[np.float64]] = {}
    acc, rows, inside = None, [], False
    with open(path) as f:
        for line in f:
            if line.startswith("ACC"):
                acc = int(line.split()[1][2:].split(".")[0])
            elif line.startswith("HMM "):
                inside, rows = True, [np.zeros(20)]
            elif line.startswith("//"):
                out[acc] = np.array(rows)
                inside = False
            elif inside:
                c = line.split()
                if len(c) >= 22 and c[0].isdigit():
                    rows.append(np.exp(-np.array(c[1:21], dtype=float)))
    return out


class Family:
    """A family's reference matrix, its statistics, and its test rows."""

    def __init__(self, acc: int, hits: pl.DataFrame, emissions: NDArray, rng, n_test: int) -> None:
        self.acc = acc
        self.emissions = emissions
        self.length = emissions.shape[0] - 1
        best = hits.sort("score", descending=True).unique(
            ["protein_id"], keep="first", maintain_order=True
        )
        chosen = (
            best.with_columns(is_rep=pl.col("protein_id") == pl.col("cluster_rep"))
            .sort(["is_rep", "score"], descending=True)
            .unique(["cluster_rep"], keep="first", maintain_order=True)
            .sort("cluster_rep")
        )
        self.units = chosen["cluster_rep"].to_numpy()
        matrix = np.zeros((chosen.height, self.length + 1), dtype=np.uint8)
        for i, r in enumerate(chosen.select("hmm_from", "ali_from", "cigar", "sequence").rows()):
            matrix[i], _ = column_residues(*r, self.length)
        self.all_rows = matrix
        held = np.sort(
            rng.choice(len(self.units), min(n_test, len(self.units) // 10), replace=False)
        )
        self.held_rows = held
        self.ref = np.delete(matrix, held, axis=0)
        self.ref_units = np.delete(self.units, held)
        self.held_info = chosen[held.tolist()].select(
            "protein_id", "ali_from", "hmm_from", "cigar", "sequence"
        )
        taken = set(chosen["protein_id"].to_list())
        others = best.filter(
            ~pl.col("protein_id").is_in(list(taken))
            & pl.col("cluster_rep").is_in(self.ref_units.tolist())
        )
        self.members = others.sample(min(n_test, others.height), seed=int(rng.integers(1 << 31)))
        self.statistics()

    def statistics(self) -> None:
        """Per column and residue: frequencies among reference units (with a pseudocount), and
        entropy, conservation and retention per column."""
        counts = np.zeros((self.length + 1, 22))
        for a in range(1, 21):
            counts[:, a] = (self.ref == a).sum(0)
        n = counts[:, 1:21].sum(1, keepdims=True)
        self.q = np.zeros((self.length + 1, 22))
        self.q[:, 1:21] = (counts[:, 1:21] + BACKGROUND) / (n + 1)
        p = self.q[:, 1:21]
        self.entropy = -(p * np.log2(p)).sum(1)
        e = np.clip(self.emissions, 1e-9, None)
        self.conservation = (e * np.log2(e / BACKGROUND)).sum(1)
        self.conservation[0] = 0
        self.hmm_q = np.zeros((self.length + 1, 22))
        self.hmm_q[:, 1:21] = e
        self.pairs = self.nearest_pairs()
        same, seen = self.pairs
        self.retention = (same.sum(0) + 1) / (seen.sum(0) + 2)

    def nearest_pairs(self, n: int = 300) -> tuple[NDArray[np.bool_], NDArray[np.bool_]]:
        """Per sampled reference unit and column: same residue as its nearest other unit, and
        whether both have a residue there."""
        sample = np.linspace(0, len(self.ref) - 1, min(n, len(self.ref))).astype(int)
        same, seen = [], []
        for i in sample:
            identity = whole_identity(self.ref, self.ref[i])
            identity[i] = -1
            j = int(identity.argmax())
            both = (self.ref[i] > 0) & (self.ref[j] > 0)
            same.append(both & (self.ref[i] == self.ref[j]))
            seen.append(both)
        return np.array(same), np.array(seen)


def whole_identity(
    matrix: NDArray[np.uint8], row: NDArray[np.uint8], min_share: float = 0.5
) -> NDArray:
    """Identity of each matrix row with ``row`` over their shared columns; -1 when they share under
    ``min_share`` of the row's residues."""
    both = (matrix > 0) & (row > 0)
    shared = both.sum(1)
    identity = ((matrix == row) & both).sum(1) / np.maximum(shared, 1)
    identity[shared < min_share * (row > 0).sum()] = -1
    return identity


def windows(row: NDArray[np.uint8], width: int, rng) -> int | None:
    """A random start whose window of columns is at least 75% occupied, or None."""
    occupied = row > 0
    cumulative = np.concatenate(([0], np.cumsum(occupied)))
    if len(row) <= width:
        return None
    filled = cumulative[width:] - cumulative[:-width]
    good = np.flatnonzero(filled >= 0.75 * width)
    return int(rng.choice(good)) if len(good) else None


def scorers(fam: Family, cols: NDArray, r: NDArray) -> dict[str, NDArray]:
    """Every score of a read (residues ``r`` at ``cols``) against all reference units; -inf where
    the overlap is too short."""
    m = fam.ref[:, cols]
    valid = (r > 0) & (r != OTHER)
    present = (m > 0) & valid
    match = (m == r) & valid
    overlap = present.sum(1)
    enough = overlap >= MIN_OVERLAP * valid.sum()
    out = {"identity": match.sum(1) / np.maximum(overlap, 1)}
    index = np.where(valid, r, 0)
    weights = {
        "conservation": np.maximum(fam.conservation[cols], 0.01),
        "inverse_conservation": 1 / (0.25 + np.maximum(fam.conservation[cols], 0)),
        "entropy": np.maximum(fam.entropy[cols], 0.01),
    }
    for name, w in weights.items():
        out[name] = (match * w).sum(1) / np.maximum((present * w).sum(1), 1e-9)
    for name, q in (("rarity_units", fam.q), ("rarity_hmm", fam.hmm_q)):
        q_r = np.maximum(q[cols, index], 1e-6)
        out[name] = (match * -np.log2(q_r)).sum(1)
        gain = np.log2((THETA + (1 - THETA) * q_r) / q_r)
        out[f"ratio_{name.split('_')[1]}"] = (match * gain).sum(1) + (
            (present & ~match).sum(1) * np.log2(1 - THETA)
        )
    out["blosum62"] = BLOSUM[index[None, :].repeat(len(m), 0), m].sum(1)
    for v in out.values():
        v[~enough] = -np.inf
    return out


m62 = substitution_matrices.load("BLOSUM62")
BLOSUM = np.zeros((22, 22))
for i_, a_ in enumerate(AMINO, 1):
    for j_, b_ in enumerate(AMINO, 1):
        BLOSUM[i_, j_] = m62[a_][b_]


def seed_sets(fam: Family, rng) -> dict[str, list[NDArray]]:
    """Per scheme, the family's seeds as arrays of columns."""
    L = fam.length
    schemes: dict[str, list[NDArray]] = {}
    for k in (3, 4, 5, 6, 8):
        schemes[f"contiguous_k{k}"] = [np.arange(c, c + k) for c in range(1, L - k + 2)]
    for pattern in ("11011", "1101011", "110100101", "11010010011"):
        offsets = np.flatnonzero(np.frombuffer(pattern.encode(), dtype=np.uint8) == 49)
        name = f"spaced_{pattern}"
        schemes[name] = [c + offsets for c in range(1, L - len(pattern) + 2)]
    ranks = {
        "entropy": fam.entropy,
        "conservation": fam.conservation,
        "retention": fam.retention,
        "efficient": fam.entropy / np.maximum(-np.log2(np.clip(fam.retention, 1e-3, 0.999)), 0.02),
    }
    for width in (12, 24):
        for k in (3, 4, 5):
            for name, rank in ranks.items():
                seen = set()
                for c in range(1, L - width + 2):
                    top = c + np.argsort(-rank[c : c + width], kind="stable")[:k]
                    seen.add(tuple(sorted(top.tolist())))
                schemes[f"top_{name}_w{width}_k{k}"] = [np.array(s) for s in sorted(seen)]
            seen = set()
            for c in range(1, L - width + 2):
                seen.add(tuple(sorted((c + rng.choice(width, k, replace=False)).tolist())))
            schemes[f"random_w{width}_k{k}"] = [np.array(s) for s in sorted(seen)]
            schemes[f"greedy_w{width}_k{k}"] = greedy_seeds(fam, width, k)
    return schemes


def greedy_seeds(fam: Family, width: int, k: int) -> list[NDArray]:
    """In windows at half-width steps, add columns one at a time by new entropy (splitting the
    reference's units further) per unit of retention lost among nearest pairs."""
    same, seen = fam.pairs
    ok = same | ~seen  # a pair with no residue there says nothing
    seeds = set()
    for c in range(1, max(2, fam.length - width + 2), max(1, width // 2)):
        labels = np.zeros(len(fam.ref), dtype=np.int64)
        alive = np.ones(len(same), dtype=bool)
        chosen: list[int] = []
        h = 0.0
        for _ in range(k):
            best, best_value, best_labels = None, -np.inf, None
            for col in range(c, min(c + width, fam.length + 1)):
                if col in chosen:
                    continue
                combined = labels * 22 + fam.ref[:, col]
                _, new_labels, counts = np.unique(combined, return_inverse=True, return_counts=True)
                p = counts / counts.sum()
                gain = -(p * np.log2(p)).sum() - h
                kept = (alive & ok[:, col]).sum() + 1
                loss = max(np.log2((alive.sum() + 1) / kept), 0.02)
                value = gain / loss
                if value > best_value:
                    best, best_value, best_labels = col, value, new_labels
            if best is None:
                break
            chosen.append(best)
            labels = best_labels
            counts = np.bincount(labels)
            p = counts[counts > 0] / counts.sum()
            h = -(p * np.log2(p)).sum()
            alive &= ok[:, best]
        if len(chosen) == k:
            seeds.add(tuple(sorted(chosen)))
    return [np.array(s) for s in sorted(seeds)]


def run_family(fam: Family, a, rng, rows, seed_rows, read_rows, speed, frags, frag_rows) -> None:
    schemes = seed_sets(fam, rng)
    reads = []
    for i, held in enumerate(fam.held_rows):
        row = fam.all_rows[held]
        identity = whole_identity(fam.ref, row)
        if identity.max() < 0:
            continue
        reads.append(
            ("novel", row, identity, float(identity.max()), fam.held_info.row(i, named=True))
        )
    for m in fam.members.iter_rows(named=True):
        row, _ = column_residues(
            m["hmm_from"], m["ali_from"], m["cigar"], m["sequence"], fam.length
        )
        unit = np.flatnonzero(fam.ref_units == m["cluster_rep"])
        if len(unit) == 0:
            continue
        truth = np.full(len(fam.ref), -1.0)
        truth[unit] = 1.0
        full = whole_identity(fam.ref, row)
        reads.append(("known", row, truth, float(max(full[unit[0]], 0)), m))
    for scenario, row, truth_score, identity, info in reads:
        start = windows(row, a.width, rng)
        if start is None:
            continue
        cols = np.arange(max(1, start), start + a.width)
        r = row[cols]
        truth = truth_score >= truth_score.max() - 1e-12
        have = (fam.ref[:, cols] > 0).sum(1) >= MIN_OVERLAP * (r > 0).sum()
        if not (truth & have).any():
            continue
        _, positions = column_residues(
            info["hmm_from"], info["ali_from"], info["cigar"], info["sequence"], fam.length
        )
        held_positions = positions[cols]
        filled = held_positions[held_positions >= 0]
        lo_position, hi_position = int(filled.min()), int(filled.max()) + 1
        index = len(read_rows)
        read_rows.append(
            (
                index,
                fam.acc,
                scenario,
                info["protein_id"],
                int(cols[0]),
                identity,
                lo_position,
                hi_position,
                info["sequence"][lo_position:hi_position],
            )
        )
        full_row = row
        if frags is not None:
            hit = frags.get((index, fam.acc))
            if hit is None:
                frag_rows.append((index, scenario, identity, False, 0, 0.0))
                speed["reads"] += 1
                continue
            row, frag_positions = column_residues(*hit, read_rows[-1][-1], fam.length)
            placed = np.flatnonzero(row > 0)
            agree = positions[placed] == frag_positions[placed] + lo_position
            frag_rows.append((index, scenario, identity, True, len(placed), float(agree.mean())))
            cols = np.arange(placed.min(), placed.max() + 1)
            r = row[cols]
        t0 = time.perf_counter()
        scores = scorers(fam, cols, r)
        speed["seconds"] += time.perf_counter() - t0
        speed["bytes"] += fam.ref.shape[0] * len(cols) * len(scores)
        speed["reads"] += 1
        oracle_right = False
        for name, s in scores.items():
            top = s == s.max()
            right = bool((top & truth).any())
            # identity over the whole domain of what was called, against the best possible
            whole = whole_identity(fam.ref, full_row)
            called = whole[top & (whole >= 0)]
            regret = float(whole[truth].max() - called.mean()) if len(called) else np.nan
            rows.append((fam.acc, scenario, identity, name, right, int(top.sum()), len(s), regret))
            if name == a.oracle:
                oracle_right = right
        lo, hi = cols[0], cols[-1]
        for scheme, seeds in schemes.items():
            count = np.zeros(len(fam.ref), dtype=np.int32)
            used = 0
            for seed in seeds:
                if seed[0] < lo or seed[-1] > hi:
                    continue
                residues = row[seed]
                if (residues == 0).any() or (residues == OTHER).any():
                    continue
                count += (fam.ref[:, seed] == residues).all(1)
                used += 1
            for need in (1, 2):
                candidates = count >= need
                seed_rows.append((
                    fam.acc, scenario, identity, scheme, need, bool((candidates & truth).any()),
                    oracle_right, float(candidates.mean()), int(candidates.sum()), len(fam.ref),
                    used, len(seeds) / max(fam.length, 1),
                ))  # fmt: skip


def read_fragment_hits(path: str) -> dict[tuple[int, int], tuple]:
    """Per (read index, family), the read's best hit on that family: (hmm_from, ali_from, cigar).
    Windows are scored within their source family; another family's hit only counts as a hit."""
    return {
        (r[0], r[1]): tuple(r[2:])
        for r in pl.read_parquet(path)
        .sort("score", descending=True)
        .select("protein_id", "pfam_accession", "hmm_from", "ali_from", "cigar")
        .unique(["protein_id", "pfam_accession"], keep="first", maintain_order=True)
        .rows()
    }


def column_tables(families: list[Family], out: Path) -> None:
    columns, residues, pairs = [], [], []
    for fam in families:
        same, seen = fam.pairs
        for j in range(1, fam.length + 1):
            columns.append(
                (
                    fam.acc,
                    j,
                    float(fam.entropy[j]),
                    float(fam.conservation[j]),
                    float(fam.retention[j]),
                    float((fam.ref[:, j] > 0).mean()),
                )
            )
        for a in range(1, 21):
            at_a = fam.ref == a
            sampled = np.linspace(0, len(fam.ref) - 1, same.shape[0]).astype(int)
            mask = at_a[sampled] & seen
            kept = int((mask & same).sum())
            rarity = -np.log2(np.maximum(fam.q[:, a], 1e-9))
            bits = float((rarity[None, :] * mask).sum())
            residues.append((fam.acc, AMINO[a - 1], int(mask.sum()), kept, bits, int(at_a.sum())))
        ref = fam.ref
        use = np.flatnonzero((ref > 0).mean(0) >= 0.5)
        h = {}
        for j in use:
            h[j] = entropy_of(ref[:, j])
        for j in use:
            for d in (1, 2, 3, 4, 5, 7, 8, 11, 15, 20, 30):
                k = j + d
                if k in h:
                    joint = entropy_of(ref[:, j].astype(np.int64) * 22 + ref[:, k])
                    pairs.append((fam.acc, d, h[j] + h[k] - joint, min(h[j], h[k]), joint))
    pl.DataFrame(
        columns,
        schema=["family", "column", "entropy", "conservation", "retention", "occupancy"],
        orient="row",
    ).write_csv(out / "columns.tsv", separator="\t")
    pl.DataFrame(
        residues,
        schema=["family", "residue", "seen_in_pairs", "kept_in_pairs", "bits", "count"],
        orient="row",
    ).write_csv(out / "residues.tsv", separator="\t")
    pl.DataFrame(
        pairs, schema=["family", "distance", "mutual", "smaller_entropy", "joint"], orient="row"
    ).write_csv(out / "pairs.tsv", separator="\t")


def entropy_of(values: NDArray) -> float:
    _, counts = np.unique(values, return_counts=True)
    p = counts / counts.sum()
    return float(-(p * np.log2(p)).sum())


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--hits", required=True, help="pfam-reannotate hits of the members")
    p.add_argument("--members", required=True)
    p.add_argument("--hmm", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--width", type=int, default=40, help="HMM columns per read window")
    p.add_argument("--per-family", type=int, default=100, help="reads per family and scenario")
    p.add_argument(
        "--families", type=int, default=0, help="only the first N families (for a trial)"
    )
    p.add_argument("--oracle", default="identity")
    p.add_argument("--seed", type=int, default=1)
    p.add_argument(
        "--fragments", help="hits of reads.parquet's fragments: use their own alignments"
    )
    a = p.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(a.seed)
    emissions = read_hmms(a.hmm)
    members = pl.read_parquet(a.members, columns=["protein_id", "cluster_rep", "sequence"])
    hits = pl.read_parquet(a.hits).filter("above_ga").join(members, on="protein_id")
    accessions = sorted(hits["pfam_accession"].unique().to_list())
    if a.families:
        accessions = accessions[: a.families]
    rows, seed_rows, read_rows, families = [], [], [], []
    speed = {"seconds": 0.0, "bytes": 0, "reads": 0}
    frags = read_fragment_hits(a.fragments) if a.fragments else None
    frag_rows: list[tuple] = []
    for acc in accessions:
        t0 = time.perf_counter()
        fam = Family(
            acc, hits.filter(pl.col("pfam_accession") == acc), emissions[acc], rng, a.per_family
        )
        families.append(fam)
        run_family(fam, a, rng, rows, seed_rows, read_rows, speed, frags, frag_rows)
        print(
            f"PF{acc:05d} units {len(fam.units)} length {fam.length} "
            f"reads so far {speed['reads']} {time.perf_counter() - t0:.0f} s",
            flush=True,
        )
    pl.DataFrame(
        rows,
        schema=["family", "scenario", "identity", "score", "right", "tied", "units", "regret"],
        orient="row",
    ).write_csv(out / "scores.tsv", separator="\t")
    pl.DataFrame(
        seed_rows,
        schema=[
            "family",
            "scenario",
            "identity",
            "scheme",
            "need",
            "found",
            "oracle_right",
            "share",
            "candidates",
            "units",
            "seeds_used",
            "seeds_per_column",
        ],
        orient="row",
    ).write_csv(out / "seeds.tsv", separator="\t")
    pl.DataFrame(
        read_rows,
        schema=[
            "read",
            "family",
            "scenario",
            "protein_id",
            "start_column",
            "identity",
            "lo",
            "hi",
            "sequence",
        ],
        orient="row",
    ).write_parquet(out / "reads.parquet")
    if frags is not None:
        pl.DataFrame(
            frag_rows,
            schema=["read", "scenario", "identity", "hit", "columns", "same_column"],
            orient="row",
        ).write_csv(out / "fragments.tsv", separator="\t")
    column_tables(families, out)
    (out / "speed.json").write_text(json.dumps(speed))


if __name__ == "__main__":
    main()
