"""hmm-annotate runs of the same proteins by different engines -> TSV: cost, and domains missed.

Usage: bench.py RESULTS... (the first is the truth, ``--engine hmmsearch``; all ``--domains true``)

``cpu_h`` is CPU time: kfpfam's or nail's, by bash's ``time``, and HMMER's, by pyhmmer's process.
``core_h`` is what the tasks held, cores by time. ``speedup`` is of ``cpu_h``.
"""

import re
import sys
from pathlib import Path

import polars as pl

KEY = ["protein_id", "model", "env_from", "env_to"]
CPU = re.compile(r"^cpu_s ([\d.]+) ([\d.]+)$", re.M)


def domains(run: Path) -> pl.LazyFrame:
    return pl.scan_parquet(run / "domains/*.parquet").filter("passes").select(KEY)


def row(run: Path, truth: Path) -> dict[str, object]:
    total = pl.read_csv(run / "summary.tsv", separator="\t").row(0, named=True)
    logs = [p.read_text() for p in (run / "logs").glob("*.log")]
    seconds = [sum(map(float, m)) for log in logs for m in CPU.findall(log)]
    kernels = {m[1] for log in logs if (m := re.search(r"^(\S+): \d+ ORFs", log, re.M))}
    trace = pl.read_csv(run / "trace.tsv", separator="\t", infer_schema=False).filter(
        pl.col("name").str.contains(r"^((MEMBER_)?SEARCH|NAIL)"),
        pl.col("status").is_in(["COMPLETED", "CACHED"]),
    )
    cores = (pl.col("realtime").cast(pl.Float64) * pl.col("cpus").cast(pl.Float64)).sum() / 3.6e6
    rss = pl.col("peak_rss").cast(pl.Float64, strict=False).max() / 1e9
    here, there = domains(run), domains(truth)
    return {
        "run": run.name,
        "proteins": total["proteins"],
        "residues": total["residues"],
        "pairs": total["pairs"],
        "filter_cpu_h": sum(seconds) / 3600,
        "hmmer_cpu_h": total["hmmer_core_s"] / 3600,
        "cpu_h": (sum(seconds) + total["hmmer_core_s"]) / 3600,
        "core_h": trace.select(cores).item(),
        "peak_gb": trace.select(rss).item(),
        "domains": total["domains_passing"],
        "missed": there.join(here, on=KEY, how="anti").select(pl.len()).collect().item(),
        "extra": here.join(there, on=KEY, how="anti").select(pl.len()).collect().item(),
        "kernel": ",".join(sorted(kernels)) or None,
        "cpus": "; ".join(sorted({log.split("\n")[0].split(": ")[-1] for log in logs})),
    }


def main() -> None:
    runs = [Path(a) for a in sys.argv[1:]]
    if not runs:
        sys.exit(__doc__)
    t = pl.DataFrame([row(r, runs[0]) for r in runs], infer_schema_length=None)
    if t.select("proteins", "residues").n_unique() > 1:
        sys.exit("the runs are not of the same proteins")
    speedup = (t["cpu_h"][0] / pl.col("cpu_h")).alias("speedup")
    missed = (100 * pl.col("missed") / t["domains"][0]).alias("missed_pct")
    t = t.with_columns(speedup, missed)
    sys.stdout.write(t.write_csv(separator="\t", float_precision=3))


if __name__ == "__main__":
    main()
