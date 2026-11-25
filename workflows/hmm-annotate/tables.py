"""hmm-annotate's small tables, in the standard library.

- ``thresholds``: each model's score thresholds, in kfpfam's model order -> TSV of model, name, seq,
  dom. From the HMM file's own lines (``--cutoff gathering|trusted|noise``), or from ``--table``: a
  TSV with a header line and model, threshold, score_type first, as KOfam's ko_list is. ``full``
  sets the sequence score's threshold and ``domain`` the domain's; ``-`` (none) never passes.
- ``summary``: the chunks' stats -> one total row.
- ``check-summary``: the chunks' checks -> a row ``all``, then one per model, worst recall first.
"""

import argparse
import csv
import gzip
import sys
from collections import Counter, defaultdict
from pathlib import Path

INF = float("inf")
LINES = {"gathering": "GA", "trusted": "TC", "noise": "NC"}


def opened(path: Path):
    return gzip.open(path, "rt") if path.suffix == ".gz" else open(path)


def cutoffs(hmm: Path, line: str) -> dict[str, tuple[float, float]]:
    """Name -> (seq, dom)."""
    out, name = {}, ""
    with opened(hmm) as f:
        for x in f:
            if x.startswith("NAME "):
                name = x.split()[1]
            elif x.startswith(f"{line} "):
                seq, dom = x.replace(";", " ").split()[1:3]
                out[name] = float(seq), float(dom)
    return out


def table(path: Path) -> dict[str, tuple[float, float]]:
    """Model -> (seq, dom)."""
    out = {}
    with opened(path) as f:
        for row in list(csv.reader(f, delimiter="\t"))[1:]:
            model, threshold, kind = row[:3]
            if threshold == "-":
                out[model] = INF, INF
            elif kind == "full":
                out[model] = float(threshold), -INF
            elif kind == "domain":
                out[model] = -INF, float(threshold)
            else:
                sys.exit(f"{path}: {model}'s score type is {kind!r}, not full or domain")
    return out


def thresholds(a: argparse.Namespace) -> None:
    with open(a.models) as f:
        models = [line.rstrip("\n").split("\t")[:2] for line in f][1:]
    if len({m for m, _ in models}) < len(models):
        sys.exit(f"{a.models}: two models have one accession")
    if a.table:
        by_model = table(a.table)
        found = [by_model.get(m, by_model.get(n)) for m, n in models]
    else:
        by_name = cutoffs(a.hmm, LINES[a.cutoff])
        found = [by_name.get(n) for _, n in models]
    missing = [m for (m, _), t in zip(models, found, strict=True) if t is None]
    if missing:
        sys.exit(
            f"{len(missing)} of {len(models)} models have no threshold in "
            f"{a.table or f'{a.hmm} ({a.cutoff})'}, {missing[0]} for one"
        )
    with open(a.out, "w") as out:
        out.write("model\tname\tseq\tdom\n")
        for (m, n), (seq, dom) in zip(models, found, strict=True):
            out.write(f"{m}\t{n}\t{seq}\t{dom}\n")
    never = sum(t == (INF, INF) for t in found)
    print(f"{len(models)} models, {never} with no threshold, which never pass")


def rows(path: Path) -> list[dict[str, str]]:
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def summary(a: argparse.Namespace) -> None:
    total: Counter = Counter()
    for p in a.stats:
        for k, v in rows(p)[0].items():
            total[k] += int(v)
    text = "\t".join(total) + "\n" + "\t".join(map(str, total.values())) + "\n"
    a.out.write_text(text)
    print(text, end="")


def check_summary(a: argparse.Namespace) -> None:
    models: dict[str, Counter] = defaultdict(Counter)
    for p in a.checks:
        for r in rows(p):
            models[r.pop("model")].update({k: int(v) for k, v in r.items()})
    every: Counter = Counter()
    for c in models.values():
        every.update(c)

    def recall(c: Counter) -> float:
        return c["domains_found"] / c["domains_truth"] if c["domains_truth"] else 1.0

    columns = [f"{k}_{n}" for k in ("pairs", "domains") for n in ("truth", "found", "extra")]
    ordered = sorted(models.items(), key=lambda x: (recall(x[1]), -x[1]["domains_truth"], x[0]))
    with open(a.out, "w") as out:
        out.write("\t".join(["model", *columns, "domain_recall"]) + "\n")
        for model, c in [("all", every), *ordered]:
            out.write("\t".join([model, *(str(c[k]) for k in columns), f"{recall(c):.4f}"]) + "\n")
    print(f"domain recall {recall(every):.4f} over {every['domains_truth']} domains")


def main() -> None:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    sub = p.add_subparsers(required=True)
    t = sub.add_parser("thresholds")
    t.set_defaults(run=thresholds)
    t.add_argument("--hmm", type=Path, required=True)
    t.add_argument("--models", type=Path, required=True, help="models.tsv of kfpfam's prep")
    t.add_argument("--cutoff", choices=sorted(LINES), default="gathering")
    t.add_argument("--table", type=Path)
    t.add_argument("--out", type=Path, default=Path("thresholds.tsv"))
    s = sub.add_parser("summary")
    s.set_defaults(run=summary)
    s.add_argument("stats", type=Path, nargs="+")
    s.add_argument("--out", type=Path, default=Path("summary.tsv"))
    c = sub.add_parser("check-summary")
    c.set_defaults(run=check_summary)
    c.add_argument("checks", type=Path, nargs="+")
    c.add_argument("--out", type=Path, default=Path("check.tsv"))
    a = p.parse_args()
    a.run(a)


if __name__ == "__main__":
    main()
