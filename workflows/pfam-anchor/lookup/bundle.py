"""The self-contained annotated index: build it, check it, and profile a sample with it.

  bundle.py build   --hits H.parquet --seqs S.parquet --hmm Pfam-A.hmm --out DIR
                    [--label ko=KO.tsv ...] [--collapse 0.7 | --clusters C.parquet --collapse 0]
                    [--kfpfam-index X.tidx --kfpfam-prep DIR] [--mgnify 2026_07] [--pfam 37.4]
  bundle.py verify  DIR            every file's SHA-256
  bundle.py profile --bundle DIR --parquet SAMPLE_pfam_alignments.parquet... --out PREFIX
                    [--floor 0.4] [--hmm Pfam-A.hmm]

A bundle:
  manifest.json   format, releases, flank, collapse, schemes, and each file's bytes and SHA-256
  rows.bin        residues, alignments and labels, a block per family (rs/src/bundle.rs)
  families.tsv    where each family's block is
  labels.<scheme>.tsv   family, label name: a family's names in number order (1 on; 0 is
                  unlabelled). A protein's several labels are one name, joined by +. A family
                  keeps its 65,535 commonest names; rows of others are unlabelled, and the
                  manifest counts them
  hmm/Pfam-A.hmm, kfpfam/index.tidx, kfpfam/prep/   what kfpfam and HMMER need to align reads
                  to the same models

H: protein_id, pfam_accession, hmm_from, ali_from, ali_to, cigar, above_ga (pfam-reannotate's
columns); only GA domains are taken. S: protein_id, sequence. Labels: protein_id<TAB>label,
several lines a protein if it has several. Unlabelled proteins stay in: without them, reads of
unlabelled proteins go to labelled relatives (ANNOTATION.md).

C: protein_id, cluster_rep. With it, H's proteins that are not their cluster's representative
are members: the first stage scans representatives only, and the members of the clusters it
keeps join the later stages (REPORT.md, "Members second pass"). A call names a cluster by its
representative. Members are only kept uncollapsed, as collapsing at 70% would merge each into
its representative.

profile writes PREFIX.<scheme>.tsv (label, reads, mean similarity), PREFIX.reads.parquet (every
alignment's call, weak ones too) and PREFIX.summary.json. A read counts once a scheme, for its
best-scoring alignment; calls under --floor of the read's own score are dropped.
"""

import argparse
import hashlib
import json
import os
import shutil
import struct
import subprocess
import sys
from pathlib import Path

import numpy as np
import polars as pl
import pyarrow.parquet as pq

HERE = Path(__file__).resolve().parent
LOOKUP = Path(os.environ.get("LOOKUP_BIN", HERE / "rs/target/release/lookup"))
# 2: members; 3: as differences; 4: no member accessions, CIGAR diffs (rs/src/bundle.rs);
# 5: labels numbered within a family
FORMAT = 5
LABELS = 65535  # a row's label is a u16
FLANK = 10
AA = "ACDEFGHIKLMNPQRSTVWY"
CODE = np.full(256, 20, np.uint8)
for _i, _a in enumerate(AA):
    CODE[ord(_a)] = _i
OPS = {"M": 0, "I": 1, "D": 2}


def sha(p: Path) -> str:
    with open(p, "rb") as f:
        return hashlib.file_digest(f, "sha256").hexdigest()


def cigar_ops(cigar: str) -> list[tuple[int, int]]:
    out, n = [], 0
    for c in cigar:
        if c.isdigit():
            n = 10 * n + int(c)
        else:
            out.append((OPS[c], n))
            n = 0
    return out


def record(
    acc: int,
    hmm_from: int,
    left: int,
    right: int,
    cigar: str,
    seq: str,
    labels: tuple[int, ...] = (),
    cluster: int = 0,
) -> bytes:
    ops = cigar_ops(cigar)
    return b"".join(
        (
            struct.pack("<QHBBH", acc, hmm_from, left, right, len(ops)),
            *(struct.pack("<BH", o, n) for o, n in ops),
            struct.pack("<H", len(seq)),
            CODE[np.frombuffer(seq.encode(), np.uint8)].tobytes(),
            struct.pack(f"<{len(labels)}H", *labels),
            struct.pack("<Q", cluster),  # a member's representative; 0 for a representative
        )
    )


def models_of(hmm: Path) -> dict[int, tuple[int, int]]:
    """Accession -> columns, version."""
    out, acc, length = {}, None, None
    with open(hmm) as f:
        for line in f:
            if line.startswith("LENG"):
                length = int(line.split()[1])
            elif line.startswith("ACC"):
                acc, _, version = line.split()[1][2:].partition(".")
                acc, version = int(acc), int(version or 0)
            elif line.startswith("//"):
                out[acc] = length, version
    return out


def place(src: Path, dst: Path) -> None:
    """A copy, not a link: rebuilding the source later must not change the bundle."""
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(src, dst)


# The rows a part is built from, one per domain: residues are the domain and up to FLANK either
# side, left and right say how many of them are flank. cluster_rep is the protein's cluster (null:
# no members). label_<scheme>: a protein's labels in that scheme, sorted and joined by +, or null.
PART = [
    "protein_id",
    "pfam_accession",
    "hmm_from",
    "cigar",
    "left",
    "right",
    "residues",
    "cluster_rep",
]


def rows_of(a: argparse.Namespace, lengths: dict[int, int]) -> tuple[pl.DataFrame, list[str]]:
    """build's inputs (whole hits, sequences, label TSVs, clusters) as part rows."""
    hits = pl.read_parquet(
        a.hits,
        columns=[
            "protein_id",
            "pfam_accession",
            "hmm_from",
            "ali_from",
            "ali_to",
            "cigar",
            "above_ga",
        ],
    ).filter("above_ga")
    unknown = set(hits["pfam_accession"].unique()) - set(lengths)
    if unknown:
        sys.exit(
            f"{len(unknown)} families of {a.hits} are not in {a.hmm}, PF{min(unknown):05d} for "
            "one: another release?"
        )
    seqs = pl.read_parquet(a.seqs, columns=["protein_id", "sequence"])
    d = hits.join(seqs, on="protein_id", how="left")
    if d["sequence"].null_count():
        sys.exit(f"{d['sequence'].null_count()} domains have no sequence in {a.seqs}")
    n = pl.col("sequence").str.len_bytes()
    d = d.with_columns(
        left=pl.min_horizontal(pl.lit(FLANK), pl.col("ali_from") - 1),
        right=pl.min_horizontal(pl.lit(FLANK), n - pl.col("ali_to")),
    ).with_columns(
        residues=pl.col("sequence").str.slice(
            pl.col("ali_from") - 1 - pl.col("left"),
            pl.col("ali_to") - pl.col("ali_from") + 1 + pl.col("left") + pl.col("right"),
        )
    )
    if a.clusters:
        c = pl.read_parquet(a.clusters, columns=["protein_id", "cluster_rep"])
        d = d.join(c, on="protein_id", how="left")
        if d["cluster_rep"].null_count():
            sys.exit(f"{d['cluster_rep'].null_count()} domains have no cluster in {a.clusters}")
    else:
        d = d.with_columns(cluster_rep=pl.lit(None, pl.Int64))
    schemes = []
    for spec in a.label:
        name, path = spec.split("=", 1)
        x = pl.read_csv(
            path,
            separator="\t",
            has_header=False,
            new_columns=["protein_id", "label"],
            schema_overrides={"label": pl.String},
        )
        x = x.group_by("protein_id").agg(
            pl.col("label").unique().sort().str.join("+").alias(f"label_{name}")
        )
        d = d.join(x, on="protein_id", how="left")
        schemes.append(name)
    return d.select(*PART, *(f"label_{s}" for s in schemes)), schemes


def label_names(d: pl.DataFrame, scheme: str) -> pl.DataFrame:
    """pfam_accession, label_<scheme>, number: each family's label names, numbered from 1."""
    c = f"label_{scheme}"
    fam = "pfam_accession"
    return (
        d.lazy()
        .filter(pl.col(c).is_not_null())
        .group_by(fam, c)
        .len()
        .sort(fam, "len", c, descending=[False, True, False])
        .filter(pl.int_range(pl.len()).over(fam) < LABELS)
        .sort(fam, c)
        .with_columns(number=pl.int_range(1, pl.len() + 1).over(fam))
        .drop("len")
        .collect()
    )


def write_part(
    d: pl.DataFrame,
    schemes: list[str],
    lengths: dict[int, int],
    collapse: float,
    out: Path,
) -> None:
    """Part rows -> out/rows.bin, families.tsv and labels.<scheme>.tsv, through bundle-build,
    family by family."""
    unknown = set(d["pfam_accession"].unique()) - set(lengths)
    if unknown:
        sys.exit(f"{len(unknown)} families are not in the models, PF{min(unknown):05d} for one")
    members = d["cluster_rep"].is_not_null() & (d["cluster_rep"] != d["protein_id"])
    if collapse > 0 and members.any():
        sys.exit("members are only kept uncollapsed: collapse 0")
    cols = [f"l{i}" for i in range(len(schemes))]
    out.mkdir(parents=True, exist_ok=True)
    dropped = {}
    for s, c in zip(schemes, cols, strict=True):
        names = label_names(d, s)
        label = f"label_{s}"
        names.select("pfam_accession", label).write_csv(
            out / f"labels.{s}.tsv", separator="\t", include_header=False, quote_style="never"
        )
        d = d.join(names, on=["pfam_accession", label], how="left").with_columns(
            pl.col("number").cast(pl.Int64).alias(c)
        )
        dropped[s] = d.filter(pl.col(label).is_not_null() & pl.col(c).is_null()).height
        if dropped[s]:
            print(f"warning: {s}: {dropped[s]} rows have labels past their family's {LABELS}")
        d = d.with_columns(pl.col(c).fill_null(0)).drop("number", label)  # 0: unlabelled
    (out / "dropped.json").write_text(json.dumps(dropped))
    d = d.with_columns(cluster=pl.when(members).then(pl.col("cluster_rep")).otherwise(0))
    # a member row only where it differs from the cluster's others: residues, start or labels
    domain = pl.col("residues").str.slice(
        pl.col("left"), pl.col("residues").str.len_bytes() - pl.col("left") - pl.col("right")
    )
    m = (
        d.filter(pl.col("cluster") > 0)
        .sort("pfam_accession", "cluster", "protein_id")
        .with_columns(key=domain)
        .unique(
            ["pfam_accession", "cluster", "key", "hmm_from", *cols],
            keep="first",
            maintain_order=True,
        )
        .drop("key")
    )
    d = pl.concat([d.filter(pl.col("cluster") == 0), m]).sort(
        "pfam_accession", pl.col("cluster") > 0, "cluster", "protein_id", "hmm_from"
    )
    p = subprocess.Popen([LOOKUP, "bundle-build", out, str(collapse)], stdin=subprocess.PIPE)
    assert p.stdin
    try:
        for (acc,), g in d.group_by("pfam_accession", maintain_order=True):
            p.stdin.write(struct.pack("<IIII", acc, lengths[acc], g.height, len(schemes)))
            for pid, hf, cg, left, right, res, cl, *lab in g.select(
                "protein_id", "hmm_from", "cigar", "left", "right", "residues", "cluster", *cols
            ).iter_rows():
                p.stdin.write(record(pid, hf, left, right, cg, res, tuple(lab), cl))
        p.stdin.close()
    except BrokenPipeError:
        pass
    if p.wait():
        sys.exit("bundle-build failed")


def finish(a: argparse.Namespace, out: Path, parts: list[Path]) -> None:
    """Parts -> one bundle in out: rows.bin and families.tsv joined (offsets shifted), labels, the
    models, kfpfam's files and the manifest."""
    # a part's schemes are in the order its rows have them
    drops = [json.loads((part / "dropped.json").read_text()) for part in parts]
    schemes = list(drops[0])
    for part, x in zip(parts, drops, strict=True):
        if list(x) != schemes:
            sys.exit(f"{part} has labels {list(x)}, and {parts[0]} {schemes}")
    dropped = {s: sum(x[s] for x in drops) for s in schemes}
    if parts == [out]:
        (out / "dropped.json").unlink()
    else:
        seen: set[int] = set()
        with open(out / "rows.bin", "wb") as rows, open(out / "families.tsv", "w") as tsv:
            at = 0
            for part in parts:
                for line in (part / "families.tsv").read_text().splitlines():
                    c = line.split("\t")
                    if int(c[0]) in seen:
                        sys.exit(f"PF{int(c[0]):05d} is in two parts")
                    seen.add(int(c[0]))
                    c[4] = str(int(c[4]) + at)
                    tsv.write("\t".join(c) + "\n")
                with open(part / "rows.bin", "rb") as f:
                    shutil.copyfileobj(f, rows, 1 << 24)
                at = rows.tell()
        for s in schemes:
            with open(out / f"labels.{s}.tsv", "wb") as f:
                for part in parts:
                    f.write((part / f"labels.{s}.tsv").read_bytes())
    place(a.hmm, out / "hmm/Pfam-A.hmm")
    if a.kfpfam_index:
        place(a.kfpfam_index, out / "kfpfam/index.tidx")
    if a.kfpfam_prep:
        for f in sorted(a.kfpfam_prep.rglob("*")):
            if f.is_file():
                place(f, out / "kfpfam/prep" / f.relative_to(a.kfpfam_prep))
    fams = pl.read_csv(out / "families.tsv", separator="\t", has_header=False)
    files = {
        str(f.relative_to(out)): {"bytes": f.stat().st_size, "sha256": sha(f)}
        for f in sorted(out.rglob("*"))
        if f.is_file()
    }
    prep = out / "kfpfam/prep/hmm.sha256"  # written by kfpfam_setup.sh
    if prep.exists() and prep.read_text().split()[0] != files["hmm/Pfam-A.hmm"]["sha256"]:
        sys.exit(f"{a.kfpfam_prep} was built from other models than {a.hmm}")
    manifest = {
        "format": FORMAT,
        "pfam": a.pfam,
        "mgnify": a.mgnify,
        "flank": FLANK,
        "collapse": a.collapse,
        "schemes": schemes,
        "rows_label_dropped": dropped,
        "families": fams.height,
        "rows": int(fams["column_3"].sum()),
        "rows_before_collapse": int(fams["column_4"].sum()),
        "member_rows": int(fams["column_7"].sum()),
        "member_rows_unreached": int(fams["column_8"].sum()),  # their representative has no row
        "files": files,
    }
    (out / "manifest.json").write_text(json.dumps(manifest, indent=1))
    b = files["rows.bin"]["bytes"]
    n = manifest["rows"] + manifest["member_rows"]
    print(
        f"{manifest['rows']} rows of {manifest['rows_before_collapse']} in {fams.height} families, "
        f"and {manifest['member_rows']} member rows "
        f"({manifest['member_rows_unreached']} unreached), "
        f"{b / max(n, 1):.1f} bytes a row, "
        f"{sum(f['bytes'] for f in files.values()) / 1e9:.2f} GB in all"
    )


def empty(out: Path) -> None:
    if out.exists() and any(out.iterdir()):
        sys.exit(f"{out} is not empty")
    out.mkdir(parents=True, exist_ok=True)


def build(a: argparse.Namespace) -> None:
    if a.clusters and a.collapse > 0:
        sys.exit("--clusters keeps members, which collapsing would merge away: add --collapse 0")
    empty(a.out)
    lengths = {k: v[0] for k, v in models_of(a.hmm).items()}
    d, schemes = rows_of(a, lengths)
    write_part(d, schemes, lengths, a.collapse, a.out)
    finish(a, a.out, [a.out])


def part_schemes(rows: list[Path]) -> list[str]:
    return [
        c.removeprefix("label_") for c in pq.read_schema(rows[0]).names if c.startswith("label_")
    ]


def build_part(a: argparse.Namespace) -> None:
    """Part rows (lookup-index's SPLIT) -> one part, its families whole."""
    empty(a.out)
    lengths = {k: v[0] for k, v in models_of(a.hmm).items()}
    schemes = part_schemes(a.rows)
    d = pl.read_parquet(a.rows, columns=[*PART, *(f"label_{s}" for s in schemes)])
    write_part(d, schemes, lengths, a.collapse, a.out)


def assemble(a: argparse.Namespace) -> None:
    empty(a.out)
    finish(a, a.out, a.parts)


def manifest_of(d: Path) -> dict:
    m = json.loads((d / "manifest.json").read_text())
    if m["format"] != FORMAT:
        sys.exit(f"{d} is format {m['format']}, and this reads {FORMAT}")
    for name, f in m["files"].items():
        p = d / name
        if not p.is_file() or p.stat().st_size != f["bytes"]:
            sys.exit(f"{p} is missing, or not {f['bytes']} bytes")
    return m


def verify(a: argparse.Namespace) -> None:
    m = manifest_of(a.bundle)
    bad = [n for n, f in m["files"].items() if sha(a.bundle / n) != f["sha256"]]
    if bad:
        sys.exit(f"changed since it was built: {', '.join(bad)}")
    print(
        f"{len(m['files'])} files as built: Pfam {m['pfam']}, MGnify {m['mgnify']}, {m['rows']} "
        "rows"
    )


def profile(a: argparse.Namespace) -> None:
    m = manifest_of(a.bundle)
    if a.hmm and sha(a.hmm) != m["files"]["hmm/Pfam-A.hmm"]["sha256"]:
        sys.exit(f"{a.hmm} is not the bundle's Pfam (release {m['pfam']}): columns would not agree")
    models = models_of(a.bundle / "hmm/Pfam-A.hmm")
    lengths = {k: v[0] for k, v in models.items()}
    parts = []
    for p in a.parquet:
        have = pq.read_schema(p).names
        x = pl.read_parquet(
            p,
            columns=[
                "read",
                "pfam",
                "hmm_from",
                "cigar",
                "residues",
                *({"pfam_version", "above_ga"} & set(have)),
            ],
        )
        x = x.with_columns(acc=pl.col("pfam").str.slice(2, 5).cast(pl.UInt32))
        unknown = set(x["acc"].unique()) - set(models)
        if unknown:
            sys.exit(
                f"{len(unknown)} of {p}'s families are not in the bundle's Pfam {m['pfam']}, "
                f"PF{min(unknown):05d} for one"
            )
        if "pfam_version" in have:
            got = x.select("acc", "pfam_version").unique()
            bad = [(f, v) for f, v in got.iter_rows() if models[f][1] != v]
            if bad:
                f, v = min(bad)
                sys.exit(
                    f"{p}: {len(bad)} of {got.height} families are other versions than the "
                    f"bundle's Pfam {m['pfam']}, PF{f:05d}.{v} for one, {models[f][1]} here"
                )
        elif not a.hmm:
            print(
                f"warning: {p} does not say which models it was aligned to; pass --hmm to check "
                "them",
                file=sys.stderr,
            )
        if a.ga_only and "above_ga" in have:
            x = x.filter("above_ga")
        parts.append(x.select("read", "pfam", "hmm_from", "cigar", "residues", "acc"))
    d = pl.concat(parts).sort("acc", maintain_order=True).with_row_index("qid")
    tmp = Path(f"{a.out}.queries.bin")
    with open(tmp, "wb") as f:
        for (acc,), g in d.group_by("acc", maintain_order=True):
            f.write(struct.pack("<IIII", acc, lengths[acc], g.height, 0))
            for qid, hf, cg, seq in g.select("qid", "hmm_from", "cigar", "residues").iter_rows():
                f.write(record(qid, hf, 0, 0, cg, seq))
    calls = Path(f"{a.out}.calls.tsv")
    r = subprocess.run(
        [LOOKUP, "bundle-query", a.bundle, tmp, calls, str(a.k), str(a.second)],
        env={**os.environ, "THREADS": str(a.threads)},
    )
    if r.returncode:
        sys.exit("bundle-query failed")
    names = ["qid", "score", "self", "ties", "acc", *m["schemes"]]
    c = pl.read_csv(
        calls,
        separator="\t",
        has_header=False,
        new_columns=names,
        schema_overrides=[pl.Int64] * len(names),
    )
    tmp.unlink()
    calls.unlink()
    c = c.with_columns(
        qid=pl.col("qid").cast(pl.UInt32), similarity=pl.col("score") / pl.col("self").clip(1)
    )
    x = d.select("qid", "read", "pfam").join(c, on="qid", how="left")
    x = x.with_columns(fam=pl.col("pfam").str.slice(2, 5).cast(pl.Int64))
    for s in m["schemes"]:
        names = pl.read_csv(
            a.bundle / f"labels.{s}.tsv",
            separator="\t",
            has_header=False,
            quote_char=None,
            raise_if_empty=False,
            schema={"fam": pl.Int64, "name": pl.String},
        ).with_columns(pl.int_range(1, pl.len() + 1).over("fam").alias(s))
        x = (
            x.join(names, on=["fam", s], how="left", maintain_order="left")
            .with_columns(pl.col("name").alias(s))
            .drop("name")
        )
    x = x.drop("fam")
    x = x.with_columns(
        protein=pl.when(pl.col("acc") > 0).then(
            pl.format("MGYP{}", pl.col("acc").cast(pl.String).str.zfill(12))
        )
    ).drop("acc")
    x.write_parquet(f"{a.out}.reads.parquet")
    best = (
        x.filter(pl.col("similarity") >= a.floor)
        .sort("score", "qid", descending=[True, False])
        .unique("read", keep="first", maintain_order=True)
    )
    reads = d["read"].n_unique()
    summary = {
        "alignments": d.height,
        "alignments_in_families_indexed": int(x["score"].is_not_null().sum()),
        "reads": reads,
        "reads_called": best.height,
        "reads_under_floor": reads - best.height,
        "floor": a.floor,
        "pfam": m["pfam"],
        "mgnify": m["mgnify"],
        "labelled_reads": {},
    }
    for s in m["schemes"]:
        t = best.filter(pl.col(s).is_not_null())
        summary["labelled_reads"][s] = t.height
        t = (
            t.group_by(pl.col(s).alias("label"))
            .agg(reads=pl.len(), similarity=pl.col("similarity").mean().round(4))
            .sort(["reads", "label"], descending=[True, False])
        )
        t.write_csv(f"{a.out}.{s}.tsv", separator="\t")
    Path(f"{a.out}.summary.json").write_text(json.dumps(summary, indent=1))
    print(json.dumps(summary))


def main() -> None:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    sub = p.add_subparsers(required=True)
    b = sub.add_parser("build")
    b.set_defaults(run=build)
    for x in ("hits", "seqs", "hmm", "out"):
        b.add_argument(f"--{x}", type=Path, required=True)
    b.add_argument("--label", action="append", default=[], metavar="SCHEME=TSV")
    b.add_argument("--collapse", type=float, default=0.7, help="0 keeps every row")
    b.add_argument("--clusters", type=Path, help="protein_id, cluster_rep: index members too")
    b.add_argument("--kfpfam-index", type=Path)
    b.add_argument("--kfpfam-prep", type=Path)
    b.add_argument("--mgnify", default="")
    b.add_argument("--pfam", default="")
    bp = sub.add_parser("build-part")
    bp.set_defaults(run=build_part)
    bp.add_argument("--rows", type=Path, nargs="+", required=True, help="part rows (bundle.PART)")
    bp.add_argument("--hmm", type=Path, required=True)
    bp.add_argument("--out", type=Path, required=True)
    bp.add_argument("--collapse", type=float, default=0.0)
    asm = sub.add_parser("assemble")
    asm.set_defaults(run=assemble)
    asm.add_argument("--parts", type=Path, nargs="+", required=True)
    asm.add_argument("--hmm", type=Path, required=True)
    asm.add_argument("--out", type=Path, required=True)
    asm.add_argument("--collapse", type=float, default=0.0, help="what the parts were built with")
    asm.add_argument("--kfpfam-index", type=Path)
    asm.add_argument("--kfpfam-prep", type=Path)
    asm.add_argument("--mgnify", default="")
    asm.add_argument("--pfam", default="")
    v = sub.add_parser("verify")
    v.set_defaults(run=verify)
    v.add_argument("bundle", type=Path)
    q = sub.add_parser("profile")
    q.set_defaults(run=profile)
    q.add_argument("--bundle", type=Path, required=True)
    q.add_argument("--parquet", type=Path, nargs="+", required=True)
    q.add_argument("--out", required=True)
    q.add_argument("--floor", type=float, default=0.4)
    q.add_argument(
        "--ga-only",
        action="store_true",
        help="leave out the alignments under the gathering thresholds (pfam_weak_evalue)",
    )
    q.add_argument("--hmm", type=Path, help="the models the Parquets were aligned to, to check")
    q.add_argument("--k", type=int, default=2000, help="rows the scan keeps (REPORT.md)")
    q.add_argument("--second", type=int, default=1000, help="of those, rows to the banded score")
    q.add_argument("--threads", type=int, default=4)
    a = p.parse_args()
    a.run(a)


if __name__ == "__main__":
    main()
