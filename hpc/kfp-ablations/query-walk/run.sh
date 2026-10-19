#!/usr/bin/env bash
#SBATCH --job-name=kfp-query-walk
#SBATCH --output=%x.%j.log
#SBATCH --time=7-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --export=ALL
# Query memory after tier 1 is read by hash range and the unit columns in unit order (plan,
# phase 6, step 40 (c)-(d)): query-memory's cells again, at a lower ladder of limits. A new
# run directory, as -resume would reuse query-memory's QUERY tasks (same script).
# Read results/query_cost.tsv as query-memory's: the smallest limit with end_exit 0 per cell
# (gate: <= 32 GB for standard at 100 M pairs), the new stages sketch (the extra read pass)
# and walk (tier 1 read in order), and wall against query-memory's at its passing limit.
# Then identical.tsv: each passing cell's profile against query-memory's (same tier and
# pairs, any limit) on the columns both have, equal to 1e-6 relative (as the golden tests).
# Cells run --greedy (no EM: coverage_greedy from gather), into results-greedy/: the cost of
# the query itself (plan, step 40, 2026-10-19). FIT=em sbatch run.sh: the shipped fits, into
# results/. KEEP=1: the pages read stay in the page cache (--keep-index-cache), the plan's
# (c2) arm, into <out>-keep/. -resume reuses a cell whose script is unchanged, whatever the
# code: after a code change, RESUME= sbatch run.sh (the profile comparison runs again too).
source ../common.sh
TIERS=${TIERS:-slim standard large}
out=results-greedy args=' --greedy'  # the leading space: not read as a flag by Nextflow
[[ ${FIT:-greedy} == em ]] && out=results args=''
[[ -n ${KEEP:-} ]] && out=$out-keep args="$args --keep-index-cache"
[[ -n $args ]] && args=(--query_args "$args") || args=()
for t in $TIERS; do
    [[ -f ../release-tiers/results-$t/index/meta.json ]] ||
        { echo "release-tiers has not finished (no results-$t/index)" >&2; exit 1; }
done
tiers=$(printf '%s,' $TIERS)
tiers=${tiers%,}
[[ $tiers == *,* ]] && tiers="{$tiers}"  # a one-name brace group matches nothing in Nextflow
nextflow run "$WF/mgnify-subset" -c ../site.config --query true --index false \
    --query_indexes "../release-tiers/results-$tiers/index" \
    --query_run ERR7738575,ERR7746321 --query_ladder 10000000,100000000 \
    --query_draws 0 --query_scratch '' --query_preload false ${args[@]+"${args[@]}"} \
    --query_memory '16 GB,24 GB,32 GB,48 GB,64 GB,128 GB' \
    --outdir "$out" ${RESUME--resume}

# Profiles are not published: find each cell's in its work directory from the traces.
"$REPO/.venv/bin/python" - "$out/trace.tsv" ../query-memory/results/trace.tsv \
    > "$out/identical.tsv" <<'EOF'
import csv, glob, re, sys
import polars as pl
from polars.testing import assert_frame_equal

def profiles(trace, work):
    """(pairs, tier) -> a passing cell's profile, from a run's trace and work directory."""
    out = {}
    for row in csv.DictReader(open(trace), delimiter="\t"):
        m = re.match(r"QUERY \((\d+) pairs, ([^,]+),", row["name"])
        got = glob.glob(f"{work}/{row['hash']}*/*.profile.tsv") if m else []
        if got:
            out.setdefault(m.groups(), []).append(got[0])
    return out

new, old = profiles(sys.argv[1], "work"), profiles(sys.argv[2], "../query-memory/work")
print("pairs\ttier\tcells\tidentical")
for key in sorted(new):
    if key not in old:
        continue
    want = pl.read_csv(old[key][0], separator="\t", infer_schema_length=None)
    ok = 0
    for path in new[key]:
        got = pl.read_csv(path, separator="\t", infer_schema_length=None)
        both = [c for c in got.columns if c in want.columns]
        try:
            assert_frame_equal(
                got.select(both), want.select(both), check_exact=False, rel_tol=1e-6, abs_tol=1e-9
            )
            ok += 1
        except AssertionError:
            pass
    print(*key, len(new[key]), ok, sep="\t")
EOF
