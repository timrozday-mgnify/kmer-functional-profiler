#!/usr/bin/env bash
#SBATCH --job-name=kfp-zi-warm-query-cost
#SBATCH --output=%x.%j.log
#SBATCH --time=7-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --export=ALL
source ../common.sh
# The zero-inflated fit started from the plain EM's coverages (plan, phase 7, step 45): the
# three release tiers (release-tiers) on the pooled gut sample at 10 M and 100 M pairs, the
# shipped fits (no --greedy, no --aai), warm (the new default) in warm/ and --zi-cold-start
# (the old start) in cold/, each run from its own directory: -resume does not key a task on
# query_args, and finds cached tasks through its launch directory's .nextflow/ cache (even
# in another work directory), so in one directory the cold arm would reuse the warm cells.
# Read <arm>/results/query_cost.tsv: fit_zi_wall_s, zi_em_iterations, zi_em_unconverged_units and
# zi_em_max_change, warm against cold per cell (cost.tsv puts them side by side). Then
# agreement.tsv: per cell, the units whose coverage_zi / present_zi / abundance_zi differ
# between the arms by > 1e-3 (relative for coverage and abundance, absolute for presence),
# and their share of the cell's abundance_zi. Plain EM is the same in both (coverage_em
# identical is the check that only the start changed).
TIERS=${TIERS:-slim standard large}
for t in $TIERS; do
    [[ -f ../release-tiers/results-$t/index/meta.json ]] ||
        { echo "release-tiers has not finished (no results-$t/index)" >&2; exit 1; }
done
tiers=$(printf '%s,' $TIERS)
tiers=${tiers%,}
[[ $tiers == *,* ]] && tiers="{$tiers}"  # a one-name brace group matches nothing in Nextflow
# Config files, not the command line: Nextflow reads a value starting with '--' as an option,
# and passes --index false as the string 'false', which the workflow reads as true (it then
# extracts the release before querying).
wf=$(cd "$WF" && pwd)
# the pooled sample once, for both arms (pileup-cost's pooled100M: ENA runs in ../reads/)
fetch_pairs pooled100M 100000000 ERR7738575 ERR7746321
reads="$PWD/pooled100M_1.fastq.gz,$PWD/pooled100M_2.fastq.gz"
for arm in ${ARMS:-warm cold}; do
    mkdir -p "$arm"
    printf 'params.index = false\nparams.query_args = "%s"\n' \
        "$([[ $arm == cold ]] && echo --zi-cold-start)" > "$arm/arm.config"
    (cd "$arm" && nextflow run "$wf/mgnify-subset" -c ../../site.config -c arm.config --query true \
        --query_indexes "../../release-tiers/results-$tiers/index" \
        --query_reads "$reads" --query_ladder 10000000,100000000 \
        --query_draws 0 --query_scratch '' --query_preload false --query_memory '128 GB' \
        --outdir results -resume)
done

# Profiles are not published: find each cell's in the work directory from each arm's trace.
"$REPO/.venv/bin/python" - <<'PY'
import csv, glob, re
import polars as pl

def profiles(arm):
    """(pairs, tier) -> a completed cell's profile, from an arm's trace."""
    out = {}
    for row in csv.DictReader(open(f"{arm}/results/trace.tsv"), delimiter="\t"):
        m = re.match(r"QUERY \((\d+) pairs, ([^,]+),", row["name"])
        got = glob.glob(f"{arm}/work/{row['hash']}*/*.profile.tsv") if m and row["exit"] == "0" else []
        if got:
            out[m.groups()] = got[0]
    return out

cols = ["coverage_em", "coverage_zi", "present_zi", "abundance_zi"]
warm, cold = profiles("warm"), profiles("cold")
rows = []
for key in sorted(warm.keys() & cold.keys()):
    w, c = (pl.read_csv(p, separator="\t", infer_schema_length=None)
            .select("unit", pl.col(cols).cast(pl.Float64))
            for p in (warm[key], cold[key]))
    j = w.join(c, on="unit", suffix="_cold")
    row = {"pairs": int(key[0]), "tier": key[1], "units": j.height}
    total = j["abundance_zi_cold"].sum()
    for col in cols:
        a, b = pl.col(col), pl.col(f"{col}_cold")
        diff = (a - b).abs() if col == "present_zi" else (a - b).abs() / pl.max_horizontal(a.abs(), b.abs(), 1e-12)
        off = j.filter(diff > 1e-3)
        row[f"{col}_differ"] = off.height
        row[f"{col}_max"] = j.select(diff.max()).item()
        if col == "abundance_zi":
            row["abundance_share_differ"] = off["abundance_zi_cold"].sum() / total if total else 0.0
    rows.append(row)
pl.DataFrame(rows).write_csv("agreement.tsv", separator="\t", float_precision=6)

keep = ["index", "pairs", "end_exit", "fit_zi_wall_s", "fit_em_wall_s", "total_wall_s",
        "zi_em_iterations", "zi_em_unconverged_units", "zi_em_max_change", "zi_em_block_rounds"]
w, c = (pl.read_csv(f"{a}/results/query_cost.tsv", separator="\t", infer_schema_length=None)
        for a in ("warm", "cold"))
w, c = (t.select([k for k in keep if k in t.columns]) for t in (w, c))
w.join(c, on=["index", "pairs"], suffix="_cold", how="full", coalesce=True).sort("index", "pairs") \
    .write_csv("cost.tsv", separator="\t", float_precision=3)
PY
