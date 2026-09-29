"""kMermaid on the fmh benchmark's KO units, with kMermaid's own training and classification.

kMermaid ships a model of 32,308 RefSeq protein clusters with no KO labels, so the benchmark
retrains it with KOs as clusters (the route its README gives for new databases, through
``get_kmer_dict_cluster`` as in ``_retrain_kmermaid_.py``). Clusters are capped at
``--max-members`` proteins, sampled at random: the shipped model averages ~55 proteins per
cluster, and uncapped KEGG KOs (~9 M proteins) would not fit in memory as its Python dicts.

- ``train --members members.parquet``: ``kmer_model.pkl`` and ``cluster_names.pkl``.
- ``classify --model DIR --reads reads.fastq``: ``kmermaid.tsv``, one row per classified read
  (``seq_name``, ``cluster_rep``, ``prot_name``, ``score``), as the ``kmermaid`` CLI writes it.
  The CLI resolves model paths inside its package directory, so it is called as a library.

Runs in the kMermaid container (``containers/kmermaid``), not the project's venv.
"""

import argparse
import pickle
from pathlib import Path

import pandas as pd
from kmermaid.command_line import MIN_LEN, SEGMENT_LENGTH, K, basepairs, gencode
from kmermaid.kmermaid import get_kmer_dict_cluster, proc_classify_fastq


def train(args: argparse.Namespace) -> None:
    members = pd.read_parquet(args.members, columns=["protein_id", "cluster_rep", "sequence"])
    members = (
        members.sample(frac=1, random_state=args.seed).groupby("cluster_rep").head(args.max_members)
    )
    bac_d = {p: [s] for p, s in zip(members["protein_id"], members["sequence"], strict=True)}
    clsd = members.groupby("cluster_rep")["protein_id"].apply(list).to_dict()
    dc, dfreq = get_kmer_dict_cluster(bac_d, clsd, K)
    # _retrain_kmermaid_.py: {kmer: {cluster: mean count per member}}
    model = {kmer: dict(zip(dc[kmer], dfreq[kmer], strict=True)) for kmer in dc}
    out = Path(args.out)
    out.mkdir(exist_ok=True)
    with open(out / "kmer_model.pkl", "wb") as f:
        pickle.dump(model, f)
    with open(out / "cluster_names.pkl", "wb") as f:
        pickle.dump({ko: ko for ko in clsd}, f)


def classify(args: argparse.Namespace) -> None:
    model = Path(args.model)
    with open(model / "cluster_names.pkl", "rb") as f:
        names = pickle.load(f)
    with open(model / "kmer_model.pkl", "rb") as f:
        dc = pickle.load(f)
    with open(args.reads) as reads:
        proc_classify_fastq(reads, args.out, names, SEGMENT_LENGTH, MIN_LEN, gencode,
                            basepairs, K, dc, False)  # fmt: skip


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="step", required=True)
    p = sub.add_parser("train")
    p.add_argument("--members", required=True)
    p.add_argument("--max-members", type=int, default=50)
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--out", default="kmermaid_model")
    p = sub.add_parser("classify")
    p.add_argument("--model", required=True)
    p.add_argument("--reads", required=True, help="uncompressed FASTQ")
    p.add_argument("--out", default="kmermaid.tsv")
    args = parser.parse_args()
    {"train": train, "classify": classify}[args.step](args)


if __name__ == "__main__":
    main()
