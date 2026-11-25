"""Inputs of the anchored-lookup comparison: targets and queries on HMM columns, as flat files.

  prep.py targets   MGnify90 slice -> hmmsearch of every member -> domains.parquet (GA and weak)
  prep.py queries   real (B11 Parquets), simulated (windows of members) and Zymo (genome genes)
  prep.py pack      -> bin/ for the Rust prototype

Targets are the 630,402 MGnify90 clusters of results-hpc/tree-index-slice: every cluster holding a
protein that the release annotates with one of 40 Pfam families, so those families are complete.
Everything is Pfam 38.2, as the B11 Parquets used here are.
"""

import argparse
import glob
import json
import sys
from pathlib import Path

import numpy as np
import polars as pl

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(HERE.parent))
import prepare  # noqa: E402
from columns import OP  # noqa: E402

OUT = ROOT / "data/pfam-anchor/lookup"
SLICE = ROOT / "results-hpc/tree-index-slice/results"
NAR = Path.home() / "Documents/nar-review-response/data"
B11 = NAR / "B11/results_main_nail"
ZYMO = {  # run: reference bundle of B9's truth
    "ERR2984773": "std",
    "SRR12324253": "std",
    "SRR12829160": "std",
    "ERR2935805": "std",
    "SRR23926920": "d6331",
    "SRR32629731": "d6331",
}
AA = "ACDEFGHIKLMNPQRSTVWY"  # code 0-19; 20: X and others; 21: no residue in the column
CODE = np.full(256, 20, np.uint8)
for i_, a_ in enumerate(AA):
    CODE[ord(a_)] = i_
NONE = 21
SHARD = 400_000  # proteins per hmmsearch
WINDOW = 50  # residues: a 150 nt read
IDENTITIES = (1.0, 0.95, 0.9, 0.8, 0.7, 0.6)
PER_FAMILY = 40  # windows per family and identity


def hmm_lengths() -> dict[int, int]:
    return hmm_lengths_of(OUT / "hmms.hmm")


def hmm_lengths_of(path: Path) -> dict[int, int]:
    out, acc, length = {}, None, None
    with open(path) as f:
        for line in f:
            if line.startswith("LENG"):
                length = int(line.split()[1])
            elif line.startswith("ACC"):
                acc = int(line.split()[1][2:].split(".")[0])
            elif line.startswith("//"):
                out[acc] = length
    return out


def search(name: str, rows, threads: int) -> pl.DataFrame:
    """hmmsearch of (integer id, sequence) rows against the 40 models, by prepare.py's
    thresholds."""
    prepare.RUNTIME, prepare.THREADS, prepare.MAX_GB = "none", threads, 6.0
    out = OUT / f"{name}.parquet"
    if not out.exists():
        faa = OUT / f"{name}.faa"
        faa.write_text("".join(f">{i}\n{s}\n" for i, s in rows))
        prepare.search(argparse.Namespace(out=OUT, threads=threads), faa, out, "members")
        faa.unlink()
    return pl.read_parquet(out)


def targets(a: argparse.Namespace) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    for f in ("hmms.hmm", "ga.tsv"):
        (OUT / f).write_bytes((ROOT / "data/pfam-anchor" / f).read_bytes())
    members = pl.read_parquet(SLICE / "members.parquet", columns=["protein_id", "sequence"])
    parts = []
    for i in range(0, members.height, SHARD):
        part = members.slice(i, SHARD)
        parts.append(search(f"members{i // SHARD:02d}", part.iter_rows(), a.threads))
        print(i, parts[-1].height, flush=True)
    pl.concat(parts).write_parquet(OUT / "domains.parquet")
    for i in range(len(parts)):
        (OUT / f"members{i:02d}.parquet").unlink()


def blosum_conditional() -> np.ndarray:
    """P(b | a, substituted), from BLOSUM62's implied joint probabilities."""
    from Bio.Align import substitution_matrices

    m = substitution_matrices.load("BLOSUM62")
    from columns import BACKGROUND

    s = np.array([[m[x][y] for y in AA] for x in AA], float)
    joint = np.outer(BACKGROUND, BACKGROUND) * np.exp(0.3176 * s)
    np.fill_diagonal(joint, 0)
    return joint / joint.sum(1, keepdims=True)


def mutate(seq: str, identity: float, cond: np.ndarray, rng) -> str:
    """Substitutions to ``identity``; below 1, each window has an indel with chance 1 - identity."""
    s = np.array([AA.find(c) for c in seq])
    hit = rng.random(len(s)) > identity
    for i in np.flatnonzero(hit & (s >= 0)):
        s[i] = rng.choice(20, p=cond[s[i]])
    out = "".join(AA[x] if x >= 0 else "X" for x in s)
    if rng.random() < 1 - identity:
        at, n = int(rng.integers(5, len(out) - 5)), int(rng.integers(1, 4))
        insert = "".join(rng.choice(list(AA), n))
        out = out[:at] + (insert + out[at:] if rng.random() < 0.5 else out[at + n :])
    return out


def queries(a: argparse.Namespace) -> None:
    wanted = pl.read_csv(OUT / "ga.tsv", separator="\t")["pfam_accession"].to_list()
    names = [f"PF{x:05d}" for x in wanted]
    keep = [
        "query",
        "set",
        "run",
        "identity",
        "source",
        "source_unit",
        "pfam_accession",
        "hmm_from",
        "ali_from",
        "cigar",
        "sequence",
    ]

    # real reads: the aligned residues only, which is what the pipeline stores
    real = []
    for p in sorted(glob.glob(str(B11 / "*/function-summary/pfam/*_alignments.parquet"))):
        run = Path(p).name.split("_")[0]
        d = pl.read_parquet(p).filter(pl.col("pfam").is_in(names))
        real.append(
            d.select(
                query=pl.lit(f"{run}:")
                + pl.col("read")
                + ":"
                + pl.col("pfam")
                + ":"
                + pl.col("frame").cast(str)
                + ":"
                + pl.col("domain").cast(str),
                set=pl.lit("real"),
                run=pl.lit(run),
                identity=pl.lit(None, pl.Float64),
                source=pl.lit(None, pl.Int64),
                source_unit=pl.lit(None, pl.Int64),
                pfam_accession=pl.col("pfam").str.slice(2).cast(pl.Int32),
                hmm_from="hmm_from",
                ali_from=pl.lit(1, pl.Int32),
                cigar="cigar",
                sequence="residues",
            )
        )
    real = (
        pl.concat(real)
        .select(keep)
        .with_row_index("protein_id")
        .with_columns(pl.col("protein_id").cast(pl.Int64))
    )
    # Re-aligned to the models the targets were aligned to: at least PF12715 (379 columns here, a
    # read at
    # 388 there) and PF07580 are other versions in the B11 Parquets, and columns must agree on both
    # sides.
    # ga here is the aligned residues', so about 3% fall under it.
    again = search("real", real.select("protein_id", "sequence").iter_rows(), a.threads).filter(
        "above_ga"
    )
    again = again.sort("score", descending=True).unique(
        ["protein_id", "pfam_accession"], keep="first"
    )
    n = real.height
    real = (
        real.drop("hmm_from", "ali_from", "cigar")
        .join(again, on=["protein_id", "pfam_accession"])
        .select(keep)
    )
    print(f"real: {real.height} of {n} realigned at GA")

    # simulated: 50-residue windows of members (not representatives), over a GA domain
    rng = np.random.default_rng(1)
    cond = blosum_conditional()
    members = pl.read_parquet(SLICE / "members.parquet")
    domains = pl.read_parquet(OUT / "domains.parquet").filter("above_ga")
    pool = (
        domains.filter(pl.col("ali_to") - pl.col("ali_from") >= 24)
        .join(members, on="protein_id")
        .filter(pl.col("sequence").str.len_chars() >= WINDOW + 4)
        .filter(pl.col("protein_id") != pl.col("cluster_rep"))
    )
    rows, meta = [], []
    for (acc,), g in pool.group_by("pfam_accession", maintain_order=True):
        for identity in IDENTITIES:
            pick = g.sample(min(PER_FAMILY, g.height), seed=int(rng.integers(1 << 31)))
            for pid, unit, af, at, seq in pick.select(
                "protein_id", "cluster_rep", "ali_from", "ali_to", "sequence"
            ).iter_rows():
                lo = max(0, af - 1 + 25 - WINDOW)  # at least 25 residues of domain
                hi = min(len(seq) - WINDOW, at - 25)
                if hi < lo:
                    continue
                start = int(rng.integers(lo, hi + 1))
                s = mutate(seq[start : start + WINDOW], identity, cond, rng)
                meta.append((len(meta), acc, identity, pid, unit, s))
                rows.append((len(rows), s))
    sim = hits_of("sim", rows, meta, a.threads)

    # Zymo genes: the genes' own domains (the truth side), searched with the same models
    zymo = []
    for bundle in ("std", "d6331"):
        for f in sorted((NAR / f"B9/truth/refs/{bundle}/genes").glob("*.faa")):
            name = None
            with open(f) as fh:
                for line in fh:
                    if line.startswith(">"):
                        name, seqs = line[1:].split()[0], []
                        zymo.append([bundle, f.stem, name, seqs])
                    else:
                        seqs.append(line.strip().rstrip("*"))
    zymo = [(i, b, g, n, "".join(s)) for i, (b, g, n, s) in enumerate(zymo)]
    d = search("zymo", ((r[0], r[4].replace("*", "X")) for r in zymo if len(r[4]) >= 20), a.threads)
    genes = pl.DataFrame(
        zymo, schema=["protein_id", "bundle", "genome", "gene", "sequence"], orient="row"
    )
    genes.join(d.filter("above_ga"), on="protein_id").write_parquet(OUT / "zymo_genes.parquet")

    q = pl.concat([real, sim.select(keep)]).with_row_index("qid")
    q.write_parquet(OUT / "queries.parquet")
    print(
        q.group_by("set", "identity")
        .agg(pl.len(), pl.col("sequence").str.len_chars().mean())
        .sort("set", "identity")
    )


def hits_of(name: str, rows: list, meta: list, threads: int) -> pl.DataFrame:
    """The windows' own GA alignments to the family they were cut from."""
    d = search(name, rows, threads).filter("above_ga")
    m = pl.DataFrame(
        meta,
        schema=["protein_id", "pfam_accession", "identity", "source", "source_unit", "sequence"],
        orient="row",
        schema_overrides={"pfam_accession": pl.Int32},
    )
    print(f"{name}: {len(meta)} windows")
    best = d.sort("score", descending=True).unique(["protein_id", "pfam_accession"], keep="first")
    return m.join(best, on=["protein_id", "pfam_accession"]).with_columns(
        query=pl.lit(f"{name}:") + pl.col("protein_id").cast(str),
        set=pl.lit(name),
        run=pl.lit(None, pl.String),
    )


def rows_of(d: pl.DataFrame, lengths: dict[int, int], f_rows, f_pos) -> None:
    """Append each domain's row of residue codes, and of protein positions (65535: none)."""
    for acc, hmm_from, ali_from, cigar, seq in d.select(
        "pfam_accession", "hmm_from", "ali_from", "cigar", "sequence"
    ).iter_rows():
        n = lengths[acc]
        row, pos = np.full(n, NONE, np.uint8), np.full(n, 65535, np.uint16)
        codes = CODE[np.frombuffer(seq.encode(), np.uint8)]
        c, p = hmm_from - 1, ali_from - 1
        for k, op in OP.findall(cigar):
            k = int(k)
            if op == "M":
                row[c : c + k] = codes[p : p + k]
                pos[c : c + k] = np.arange(p, p + k)
                c, p = c + k, p + k
            elif op == "I":
                p += k
            else:
                c += k
        f_rows.write(row.tobytes())
        f_pos.write(pos.tobytes())


def pack_queries(q: pl.DataFrame, b: Path, fams: list[int], lengths: dict[int, int]) -> None:
    fam_index = {f: i for i, f in enumerate(fams)}
    (b / "q.fam").write_bytes(
        np.array([fam_index[x] for x in q["pfam_accession"]], np.uint32).tobytes()
    )
    seqs = [CODE[np.frombuffer(s.encode(), np.uint8)] for s in q["sequence"]]
    (b / "q.off").write_bytes(
        np.concatenate([[0], np.cumsum([len(s) for s in seqs])]).astype(np.uint64).tobytes()
    )
    (b / "q.seq").write_bytes(np.concatenate(seqs).tobytes())
    with open(b / "q.rows", "wb") as f1, open(b / "q.pos", "wb") as f2:
        rows_of(q, lengths, f1, f2)


def pack(a: argparse.Namespace) -> None:
    """bin/: reps' sequences; rows per arm, by family; queries. All little-endian, no headers;
    sizes are in meta.json."""
    b = OUT / "bin"
    b.mkdir(exist_ok=True)
    lengths = hmm_lengths()
    fams = sorted(lengths)
    members = pl.read_parquet(
        SLICE / "members.parquet", columns=["protein_id", "cluster_rep", "sequence"]
    )
    reps = members.filter(pl.col("protein_id") == pl.col("cluster_rep")).sort("cluster_rep")
    unit_of = {u: i for i, u in enumerate(reps["cluster_rep"].to_list())}
    offsets = np.concatenate([[0], np.cumsum(reps["sequence"].str.len_bytes().to_numpy())]).astype(
        np.uint64
    )
    (b / "reps.off").write_bytes(offsets.tobytes())
    with open(b / "reps.seq", "wb") as f:
        for s in reps["sequence"]:
            f.write(CODE[np.frombuffer(s.encode(), np.uint8)].tobytes())
    pl.DataFrame({"unit": reps["cluster_rep"]}).write_parquet(b / "units.parquet")

    domains = pl.read_parquet(OUT / "domains.parquet")
    meta = {"families": fams, "lengths": [lengths[f] for f in fams], "units": reps.height}
    # reps: each representative's domains, GA and weak. members: distinct rows of every member, per
    # unit
    d = domains.join(members, on="protein_id").sort(
        "pfam_accession", "cluster_rep", "protein_id", "domain"
    )
    arms = {"reps": d.filter(pl.col("protein_id") == pl.col("cluster_rep"))}
    if a.members:
        arms["members"] = (
            d.filter("above_ga")
            .with_columns(
                key=pl.col("sequence").str.slice(
                    pl.col("ali_from") - 1, pl.col("ali_to") - pl.col("ali_from") + 1
                )
            )
            .unique(
                ["pfam_accession", "cluster_rep", "key", "hmm_from"],
                keep="first",
                maintain_order=True,
            )
            .drop("key")
        )
    for arm, x in arms.items():
        n_rows = []
        with open(b / f"{arm}.rows", "wb") as f1, open(b / f"{arm}.pos", "wb") as f2:
            for acc in fams:
                g = x.filter(pl.col("pfam_accession") == acc)
                n_rows.append(g.height)
                rows_of(g, lengths, f1, f2)
        unit = np.array([unit_of[u] for u in x["cluster_rep"]], np.uint32)
        (b / f"{arm}.unit").write_bytes(unit.tobytes())
        (b / f"{arm}.ga").write_bytes(x["above_ga"].cast(pl.UInt8).to_numpy().tobytes())
        meta[arm] = n_rows
        print(
            arm,
            x.height,
            "rows",
            sum(n * lengths[f] for n, f in zip(n_rows, fams, strict=False)) / 1e6,
            "MB",
        )

    q = pl.read_parquet(OUT / "queries.parquet").sort("qid")
    pack_queries(q, b, fams, lengths)
    src = [unit_of.get(u, 2**32 - 1) for u in q["source_unit"]]
    (b / "q.src").write_bytes(np.array(src, np.uint32).tobytes())
    meta["queries"] = q.height
    (b / "meta.json").write_text(json.dumps(meta))
    # the same for Rust, one line a family: accession, length, rows of each arm
    (b / "families.tsv").write_text(
        "".join(
            f"{f}\t{lengths[f]}\t{meta['reps'][i]}\t{meta.get('members', [0] * len(fams))[i]}\n"
            for i, f in enumerate(fams)
        )
    )

    # the Zymo genes against the same targets, in their own directory
    z = OUT / "binz"
    z.mkdir(exist_ok=True)
    for f in b.iterdir():
        if not f.name.startswith(("q.", "truth.", "meta")):
            (z / f.name).unlink(missing_ok=True)
            (z / f.name).symlink_to(f)
    genes = (
        pl.read_parquet(OUT / "zymo_genes.parquet")
        .sort("pfam_accession", "protein_id", "domain")
        .with_row_index("gid")
    )
    genes.write_parquet(OUT / "zymo_queries.parquet")
    pack_queries(genes, z, fams, lengths)
    print("zymo gene domains", genes.height)


def main() -> None:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("step", choices=["targets", "queries", "pack"])
    p.add_argument("--threads", type=int, default=6)
    p.add_argument("--members", action="store_true", help="pack every member's rows too")
    a = p.parse_args()
    prepare.watchdog(8)
    {"targets": targets, "queries": queries, "pack": pack}[a.step](a)


if __name__ == "__main__":
    main()
