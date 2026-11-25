"""Held-out windows through the bundle: only pairs kfpfam (run from the bundle's files) passes, and
alignments under GA too.
Prints recall over all KO windows, aligned or not. Needs heldout.py run, and kfpfam's pairs in
annot_heldout/pairs.tsv."""

import subprocess
import sys
from pathlib import Path

import polars as pl

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import heldout as h  # noqa: E402

d, held = h.data()
w = h.windows(d, held)
hits = h.search("windows", list(w.select("w", "sequence").iter_rows())).join(
    w.select(pl.col("w").alias("protein_id"), "sequence"), on="protein_id"
)
models = (
    pl.read_csv(h.H / "b0.7/kfpfam/prep/models.tsv", separator="\t")["acc"]
    .str.slice(2)
    .cast(pl.Int32)
    .to_list()
)
pairs = pl.read_csv(h.H / "pairs.tsv", separator="\t", has_header=False).select(
    a=pl.col("column_1"), b=pl.col("column_2")
)
print(pairs.head(2).rows(), pairs["a"].max(), pairs["b"].max())
orf, model = ("b", "a") if pairs["a"].max() > 12_000 else ("a", "b")
p = pairs.select(
    protein_id=pl.col(orf).cast(pl.Int64),
    pfam_accession=pl.col(model).replace_strict(dict(enumerate(models)), return_dtype=pl.Int32),
)
ga = hits.filter("above_ga")
kept = ga.join(p.unique(), on=["protein_id", "pfam_accession"], how="semi")
print(
    f"kfpfam passes {kept.height} of {ga.height} GA alignments ({kept.height / ga.height:.4f}), "
    f"{p.height / w.height:.1f} candidates a window"
)
ko = dict(zip(d["protein_id"].to_list(), d["ko"].to_list(), strict=False))
truth = {str(i): set(ko[x]) for i, x in w.select("w", "protein_id").iter_rows()}
total = sum(map(bool, truth.values()))
arms = {"GA, all": ga, "GA, kfpfam's": kept}
for e in (0.01, 1.0, 10.0):  # E at Z = 10^6 models x reads, as searched
    arms[f"E <= {e:g}"] = hits.filter(pl.col("above_ga") | (pl.col("seq_evalue") <= e))
# as the pipeline can give them: the domain's own E too, and only kfpfam's candidates
for e in (1.0, 10.0):
    both = hits.filter(
        pl.col("above_ga") | ((pl.col("seq_evalue") <= e) & (pl.col("i_evalue") <= e))
    )
    arms[f"E, i-E <= {e:g}"] = both
    arms["... kfpfam's"] = both.join(p.unique(), on=["protein_id", "pfam_accession"], how="semi")
    if e == 1:
        arms["... kfpfam's, 1"] = arms.pop("... kfpfam's")
for name, x in arms.items():
    x.select(
        read=pl.col("protein_id").cast(pl.String),
        pfam=pl.format("PF{}", pl.col("pfam_accession").cast(pl.String).str.zfill(5)),
        hmm_from="hmm_from",
        cigar="cigar",
        residues=pl.col("sequence").str.slice(
            pl.col("ali_from") - 1, pl.col("ali_to") - pl.col("ali_from") + 1
        ),
    ).write_parquet(h.H / "arm.parquet")
    r = subprocess.run(
        [
            sys.executable,
            HERE.parent / "bundle.py",
            "profile",
            "--bundle",
            h.H / "b0.7",
            "--parquet",
            h.H / "arm.parquet",
            "--out",
            h.H / "arm",
            "--hmm",
            h.HMM,
        ],
        capture_output=True,
        text=True,
    )
    assert r.returncode == 0, r.stderr
    c = pl.read_parquet(h.H / "arm.reads.parquet")
    for floor in (0.0, 0.4):
        b = (
            c.filter(pl.col("similarity") >= floor)
            .sort("score", "qid", descending=[True, False])
            .unique("read", keep="first", maintain_order=True)
        )
        called = {k: set(v.split("+")) for k, v in b.select("read", "ko").iter_rows() if v}
        right = sum(bool(truth[k] & v) for k, v in called.items())
        print(
            f"{name:14s} floor {floor}: windows aligned {c['read'].n_unique()}, called "
            f"{len(called)}, recall of all KO windows {right / total:.3f}, precision "
            f"{right / len(called):.3f}"
        )
