"""Tree index: a family-complete slice of the release (docs/tree-index-plan.md, Progress log:
precision and read classification). A cluster sample (1 in N MGnify90 clusters) almost never
holds a protein's neighbouring units, so it cannot show reads confused between near units;
a slice keeps every unit of chosen Pfam families.

Subcommands:

- ``families``: from ``--pfam`` (protein_id, pfam_accession), ``--families`` accessions at
  random (``--seed``) among those with ``--min-proteins`` to ``--max-proteins`` annotated
  proteins. Writes ``families.tsv`` (accession, proteins), ``chosen.parquet`` (their
  protein ids) and ``pfam.parquet`` (their rows: Pfam precision against these families).
- ``take``: one members bucket (``protein_id``, ``cluster_rep``, ...; clusters whole within a
  bucket) -> every member of each cluster holding a chosen protein, as ``--out``.
- ``concat``: the takes -> one members table, by cluster.
- ``mini``: a laptop-sized slice of a built slice (``--index`` and its ``--members``): whole
  clumps of the index's spanning tree (every edge kept, so each protein's layer group and its
  neighbouring groups come whole), a share (``--proteins`` over all) of the clumps taken
  systematically in size order, so the clump a random protein sits in is distributed as in the
  slice. Proteins with no
  index position (no k-mer) are dropped. Writes ``--out`` (members) and ``mini.json``.
"""

import argparse
import json
from pathlib import Path

import duckdb
import numpy as np
import polars as pl
from scipy.sparse import coo_array
from scipy.sparse.csgraph import connected_components


def families(a: argparse.Namespace) -> None:
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    con = duckdb.connect()
    con.sql(f"SET threads = {a.threads}")
    con.sql(
        f"CREATE TABLE fam AS SELECT pfam_accession, count(DISTINCT protein_id) AS proteins "
        f"FROM read_parquet('{a.pfam}') GROUP BY 1 "
        f"HAVING proteins BETWEEN {a.min_proteins} AND {a.max_proteins}"
    )
    # a seeded hash order: reproducible whatever the threads
    chosen = con.sql(
        f"SELECT * FROM fam ORDER BY hash(pfam_accession || '{a.seed}') LIMIT {a.families}"
    ).pl()
    chosen.sort("pfam_accession").write_csv(out / "families.tsv", separator="\t")
    accs = ", ".join(f"'{x}'" for x in chosen["pfam_accession"])
    rows = con.sql(
        f"SELECT DISTINCT protein_id, pfam_accession FROM read_parquet('{a.pfam}') "
        f"WHERE pfam_accession IN ({accs}) ORDER BY protein_id"
    ).pl()
    rows.write_parquet(out / "pfam.parquet")
    rows.select("protein_id").unique().sort("protein_id").write_parquet(out / "chosen.parquet")
    info = {"families": chosen.height, "proteins": rows["protein_id"].n_unique()}
    (out / "families.json").write_text(json.dumps(info, indent=1) + "\n")
    print(json.dumps(info))


def take(a: argparse.Namespace) -> None:
    con = duckdb.connect()
    con.sql(f"SET threads = {a.threads}")
    con.sql(
        f"COPY (SELECT m.* FROM read_parquet('{a.members}') m WHERE m.cluster_rep IN ("
        f"SELECT cluster_rep FROM read_parquet('{a.members}') "
        f"WHERE protein_id IN (SELECT protein_id FROM read_parquet('{a.chosen}'))) "
        f"ORDER BY m.cluster_rep, m.protein_id) TO '{a.out}' (FORMAT parquet)"
    )
    n = con.sql(f"SELECT count(*), count(DISTINCT cluster_rep) FROM read_parquet('{a.out}')")
    print(json.dumps(dict(zip(("proteins", "units"), n.fetchone() or (0, 0), strict=True))))


def concat(a: argparse.Namespace) -> None:
    parts = ", ".join(f"'{x}'" for x in a.parts)
    duckdb.sql(
        f"COPY (SELECT * FROM read_parquet([{parts}]) ORDER BY cluster_rep, protein_id) "
        f"TO '{a.out}' (FORMAT parquet)"
    )


def mini(a: argparse.Namespace) -> None:
    idx = Path(a.index)
    ea, eb = (np.load(idx / f"{x}.npy").astype(np.int64) for x in ("edge_a", "edge_b"))
    n = len(np.load(idx / "pos_len.npy", mmap_mode="r"))
    _, clump = connected_components(
        coo_array((np.ones(len(ea)), (ea, eb)), shape=(n, n)), directed=False
    )
    prot = pl.read_parquet(idx / "proteins.parquet", columns=["protein_id", "pos"]).filter(
        pl.col("pos") >= 0
    )
    per = np.bincount(clump[prot["pos"].to_numpy()], minlength=int(clump.max()) + 1)
    # systematic: clumps in size order (ties at random), every 1/share-th from a random start,
    # so the few large clumps that dominate a protein's view are kept at the slice's rate
    rng = np.random.default_rng(a.seed)
    share = min(a.proteins / per.sum(), 1.0)
    rank = np.empty(len(per), dtype=np.int64)
    rank[np.lexsort((rng.random(len(per)), -per))] = np.arange(len(per))
    u = rng.random()
    keep = np.floor((rank + 1 + u) * share) > np.floor((rank + u) * share)
    chosen = prot.filter(pl.Series(keep[clump[prot["pos"].to_numpy()]]))["protein_id"]
    members = pl.read_parquet(a.members).filter(pl.col("protein_id").is_in(chosen.implode()))
    members.sort("cluster_rep", "protein_id").write_parquet(a.out)
    sizes = per[keep & (per > 0)]
    info = {
        "proteins": members.height,
        "units": members["cluster_rep"].n_unique(),
        "clumps": int(len(sizes)),
        "largest_clump": int(sizes.max(initial=0)),
        # the clump a random protein sits in, here and in the whole slice
        "clump_seen_by_protein": float((sizes**2).sum() / max(sizes.sum(), 1)),
        "clump_seen_by_protein_slice": float((per**2).sum() / max(per.sum(), 1)),
    }
    (Path(a.out).parent / "mini.json").write_text(json.dumps(info, indent=1) + "\n")
    print(json.dumps(info))


def main() -> None:
    parser = argparse.ArgumentParser(description=(__doc__ or "").split("\n\n")[0])
    sub = parser.add_subparsers(dest="command", required=True)
    p = sub.add_parser("families")
    p.add_argument("--pfam", required=True)
    p.add_argument("--families", type=int, default=40)
    p.add_argument("--min-proteins", type=int, default=20_000)
    p.add_argument("--max-proteins", type=int, default=150_000)
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--threads", type=int, default=2)
    p.add_argument("--out", required=True)
    p = sub.add_parser("take")
    p.add_argument("--members", required=True, help="one members bucket")
    p.add_argument("--chosen", required=True, help="families' chosen.parquet")
    p.add_argument("--threads", type=int, default=2)
    p.add_argument("--out", required=True)
    p = sub.add_parser("concat")
    p.add_argument("--parts", nargs="+", required=True)
    p.add_argument("--out", required=True)
    p = sub.add_parser("mini")
    p.add_argument("--index", required=True, help="the slice's built index (merge.py's)")
    p.add_argument("--members", required=True, help="the slice's members.parquet")
    p.add_argument("--proteins", type=int, default=300_000)
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--out", required=True)
    a = parser.parse_args()
    {"families": families, "take": take, "concat": concat, "mini": mini}[a.command](a)


if __name__ == "__main__":
    main()
