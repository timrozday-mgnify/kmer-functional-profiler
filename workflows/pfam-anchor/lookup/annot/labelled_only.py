"""KO precision when the index holds only KO-labelled proteins (no unlabelled ones to compete) ->
printed.
labelled_only.py build ; lookup best OUT/annot_ko ; labelled_only.py score"""

import sys
from pathlib import Path

import numpy as np
import polars as pl

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE.parents[1]))
from analyze import blosum  # noqa: E402
from prep import CODE, OUT  # noqa: E402

D, K = OUT / "annot", OUT / "annot_ko"
d = pl.read_parquet(D / "proteins.parquet").with_row_index("i")
lab = d.filter("has_ko")
if sys.argv[1] == "build":
    K.mkdir(exist_ok=True)
    seqs = [CODE[np.frombuffer(s.encode(), np.uint8)] for s in lab["sequence"]]
    (K / "reps.off").write_bytes(
        np.concatenate([[0], np.cumsum([len(s) for s in seqs])]).astype(np.uint64).tobytes()
    )
    (K / "reps.seq").write_bytes(np.concatenate(seqs).tobytes())
    g = np.fromfile(D / "reps.group", "<u4")[lab["i"].to_numpy()]
    (K / "reps.group").write_bytes(g.tobytes())
    for f in ("q.off", "q.seq", "q.group"):
        (K / f).write_bytes((D / f).read_bytes())
else:
    ko, all_ko = lab["ko"].to_list(), d["ko"].to_list()
    w = pl.read_parquet(D / "windows.parquet")
    self_ = np.diag(blosum())
    seq = np.frombuffer((D / "q.seq").read_bytes(), np.uint8).reshape(-1, 50)
    rows = []
    for name, p, labels in (("all proteins", D, all_ko), ("KO proteins only", K, ko)):
        b = pl.read_csv(p / "best.tsv", separator="\t", schema_overrides={"targets": pl.String})
        for (pid, _, inp), (qid, sc, _, t) in zip(w.iter_rows(), b.iter_rows(), strict=False):
            if inp < 25 or not t:
                continue
            got = [tuple(labels[int(x)]) for x in t.split(",")]
            call = set(max(set(got), key=got.count))
            own = set(all_ko[pid])
            f = sc / max(1, sum(int(self_[c]) for c in seq[qid] if c < 20))
            rows.append((name, f, bool(own), bool(call), bool(own & call)))
    o = pl.DataFrame(rows, schema=["index", "sim", "has", "called", "right"], orient="row")
    for floor in (0.0, 0.4, 0.6):
        x = o.filter(pl.col("sim") >= floor)
        print(
            floor,
            x.group_by("index")
            .agg(
                pl.len(),
                (pl.col("right").sum() / pl.col("has").sum()).alias("recall"),
                (pl.col("right").sum() / pl.col("called").sum()).alias("precision"),
            )
            .sort("index")
            .rows(),
        )
