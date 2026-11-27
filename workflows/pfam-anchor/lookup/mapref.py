"""Every read against the proteins the lookup called: the first measurement of the mapper
(docs/pfam-anchor-plan.md, "Mapping every read to the lookup's proteins").

  mapref.py ref     --reads S_pfam_lookup.reads.parquet --store 'members.bucket*.parquet'
                    --out ref.faa [--floor 0.4]
  mapref.py measure --reads S_pfam_lookup.reads.parquet --ref ref.faa
                    --diamond r1.tsv r2.tsv --out PREFIX
  mapref.py check   the self-check

ref: the representative of every cluster a read was called to (its best call at or over --floor,
as bundle.py profile counts them), from the release's members tables (kfp's store).

measure reads DIAMOND lines of each mate, ``qseqid sseqid pident length bitscore btop``, and
writes PREFIX.summary.json (reference size, mates and pairs aligned, a pair's lookup protein
against its mates' best hits) and PREFIX.identity.tsv (each mate's best line, by identity bin,
whether the pair has a lookup call, and which seed shapes have a hit in the alignment's runs of
identities: the real counterpart of the plan's iid table).
"""

import argparse
import json
import re
import sys
from pathlib import Path

import polars as pl

from kmer_functional_profiler.store import read_units

# the plan's shapes; a seed hits where every 1 lies on an identity within one gapless run
SHAPES = {
    "k4": ["1111"],
    "k5": ["11111"],
    "k6": ["111111"],
    "golomb4_8": ["10110001"],
    "s323": ["11100111"],
    "two5_9": ["110100011", "110010101"],
    "golomb5_12": ["110010000101"],
}
BINS = [0, 40, 50, 60, 70, 80, 90, 95, 100.01]
BTOP = re.compile(r"(\d+)|(..)")


def mate_id(name: str) -> str:
    """A read's pair id: its first word less any /1 or /2."""
    name = name.split()[0]
    return name[:-2] if name[-2:] in ("/1", "/2") else name


def runs(btop: str) -> list[str]:
    """An alignment's gapless runs as identity strings ('1' identical, '0' substituted)."""
    out, cur = [], []
    for n, pair in BTOP.findall(btop):
        if n:
            cur.append("1" * int(n))
        elif "-" in pair or "/" in pair or "\\" in pair:  # a gap or a frameshift ends a run
            out.append("".join(cur))
            cur = []
        else:
            cur.append("0")
    out.append("".join(cur))
    return [r for r in out if r]


def seeded(rs: list[str], shapes: list[str]) -> bool:
    for r in rs:
        for s in shapes:
            pos = [i for i, c in enumerate(s) if c == "1"]
            for o in range(len(r) - len(s) + 1):
                if all(r[o + i] == "1" for i in pos):
                    return True
    return False


def called(path: Path, floor: float) -> pl.DataFrame:
    """Each read's best call at or over the floor: read, pair, protein (an int), similarity."""
    x = pl.read_parquet(path, columns=["read", "protein", "score", "similarity"])
    return (
        x.filter(pl.col("protein").is_not_null() & (pl.col("similarity") >= floor))
        .sort("score", descending=True)
        .unique("read", keep="first")
        .with_columns(
            pair=pl.col("read").map_elements(mate_id, return_dtype=pl.String),
            protein=pl.col("protein").str.slice(4).cast(pl.Int64),
        )
    )


def ref(a: argparse.Namespace) -> None:
    c = called(a.reads, a.floor)
    reps = c["protein"].unique().sort().to_list()
    # ponytail: representatives only; members (kfp's fetch, spread by length) once the
    # identity table says reads sit between a representative and its members
    s = read_units(a.store, reps).filter(pl.col("protein_id") == pl.col("cluster_rep"))
    if s.height < len(reps):
        print(
            f"warning: {len(reps) - s.height} of {len(reps)} proteins not in the store",
            file=sys.stderr,
        )
    with open(a.out, "w") as f:
        for pid, seq in s.select("protein_id", "sequence").iter_rows():
            f.write(f">{pid}\n{seq}\n")
    print(
        json.dumps(
            {
                "reads_called": c.height,
                "proteins": len(reps),
                "written": s.height,
                "residues": int(s["sequence"].str.len_bytes().sum()),
            }
        )
    )


def best_lines(path: Path, mate: int) -> pl.DataFrame:
    names = ["qseqid", "sseqid", "pident", "length", "bitscore", "btop"]
    d = pl.read_csv(
        path,
        separator="\t",
        has_header=False,
        new_columns=names,
        schema_overrides={"sseqid": pl.Int64, "btop": pl.String},
        raise_if_empty=False,
    )
    if d.is_empty():
        return pl.DataFrame(schema={**d.schema, "mate": pl.Int8, "pair": pl.String})
    return (
        d.sort("bitscore", descending=True)
        .unique("qseqid", keep="first")
        .with_columns(
            mate=pl.lit(mate, pl.Int8),
            pair=pl.col("qseqid").map_elements(mate_id, return_dtype=pl.String),
        )
    )


def measure(a: argparse.Namespace) -> None:
    c = called(a.reads, a.floor)
    lookup = (
        c.sort("score", descending=True)
        .unique("pair", keep="first")
        .select("pair", lookup="protein")
    )
    b = pl.concat([best_lines(p, m) for m, p in enumerate(a.diamond, 1)])
    b = b.join(lookup, on="pair", how="left")
    shape_hits = {k: [] for k in SHAPES}
    for btop in b["btop"]:
        rs = runs(btop)
        for k, v in SHAPES.items():
            shape_hits[k].append(seeded(rs, v))
    b = b.with_columns([pl.Series(f"seed_{k}", v) for k, v in shape_hits.items()]).with_columns(
        bin=pl.col("pident").cut(BINS[1:-1], left_closed=True),
        pair_called=pl.col("lookup").is_not_null(),
        same_protein=pl.col("lookup") == pl.col("sseqid"),
    )
    seeds = [f"seed_{k}" for k in SHAPES]
    table = (
        b.group_by("pair_called", "bin")
        .agg(
            mates=pl.len(),
            same_protein=pl.col("same_protein").mean(),
            **{s: pl.col(s).mean() for s in seeds},
        )
        .sort("pair_called", "bin")
    )
    table.write_csv(f"{a.out}.identity.tsv", separator="\t")
    residues = proteins = 0
    with open(a.ref) as f:
        for line in f:
            if line.startswith(">"):
                proteins += 1
            else:
                residues += len(line.strip())
    called_pairs = set(lookup["pair"])
    aligned_pairs = set(b["pair"])
    summary = {
        "reference_proteins": proteins,
        "reference_residues": residues,
        "mates_aligned": b.height,
        "pairs_aligned": len(aligned_pairs),
        "pairs_called_by_lookup": len(called_pairs),
        "called_pairs_aligned": len(called_pairs & aligned_pairs),
        "aligned_pairs_not_called": len(aligned_pairs - called_pairs),
        "mates_of_called_pairs_on_lookup_protein": float(
            b.filter("pair_called")["same_protein"].mean() or 0
        ),
        "mates_by_identity": {
            k: n for k, n in b["bin"].cast(pl.String).value_counts().sort("bin").iter_rows()
        },
        "seeded_all_mates": {s: float(b[s].mean() or 0) for s in seeds},
    }
    Path(f"{a.out}.summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    print(json.dumps(summary, indent=1))


def check() -> None:
    assert mate_id("ERR1.5/1 x") == "ERR1.5" and mate_id("ERR1.5 merged") == "ERR1.5"
    assert runs("4AG3-K2") == ["11110111", "11"]
    assert runs("5\\A3") == ["11111", "111"]  # DIAMOND's frameshift ends the run
    assert seeded(["1101111"], SHAPES["k4"]) and not seeded(["1110111"], SHAPES["k4"])
    g = SHAPES["golomb4_8"]
    assert seeded(["10110001"], g) and not seeded(["1011000"], g)
    print("ok")


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("ref")
    r.add_argument("--reads", type=Path, required=True)
    r.add_argument("--store", type=Path, nargs="+", required=True)
    r.add_argument("--out", type=Path, required=True)
    r.add_argument("--floor", type=float, default=0.4)
    m = sub.add_parser("measure")
    m.add_argument("--reads", type=Path, required=True)
    m.add_argument("--ref", type=Path, required=True)
    m.add_argument("--diamond", type=Path, nargs="+", required=True)
    m.add_argument("--out", required=True)
    m.add_argument("--floor", type=float, default=0.4)
    sub.add_parser("check")
    a = p.parse_args()
    {"ref": ref, "measure": measure, "check": lambda _: check()}[a.cmd](a)


if __name__ == "__main__":
    main()
