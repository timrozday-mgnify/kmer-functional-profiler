#!/usr/bin/env bash
# The whole comparison, one step at a time, on a 24 GB laptop: about 1.5 h, peak 8 GB (each step watched or capped).
#   ./run.sh [STEP...]     steps: targets queries pack truth cand sub follow scan tables annot
set -euo pipefail
cd "$(dirname "$0")/../../.."                      # the repository
H=workflows/pfam-anchor/lookup
OUT=data/pfam-anchor/lookup
PY="uv run --offline --no-sync python"
L=$H/rs/target/release/lookup
steps=("$@")
[ ${#steps[@]} -gt 0 ] || steps=(targets queries pack truth cand sub follow scan tables annot)
(cd $H/rs && cargo test -q --release --offline && cargo build -q --release --offline)
for s in "${steps[@]}"; do
  echo "== $s"
  case $s in
    targets) $PY $H/prep.py targets ;;             # 12 min: every member of the slice against 40 models
    queries) $PY $H/prep.py queries ;;
    pack) $PY $H/prep.py pack --members ;;
    truth)                                          # 8 and 40 min on 6 threads
      THREADS=6 $L truth $OUT/binz
      THREADS=6 $L truth $OUT/bin ;;
    cand)
      THREADS=4 $L cand $OUT/bin reps 0 $OUT/full_reps
      THREADS=4 $L cand $OUT/bin reps_ga 0 $OUT/full_reps_ga golomb4b,bits14,scan
      THREADS=2 $L cand $OUT/bin members 0 $OUT/full_members golomb4b,k5,bits14,scan ;;
    sub)
      for n in 2 4; do THREADS=4 $L cand $OUT/bin reps $n $OUT/sub${n}_reps; done ;;
    follow)
      THREADS=4 $L follow $OUT/bin $OUT/follow 1100101/h1 110010000101/h1 bits8/h1 bits10/h1 top100 top500 ;;   # 23 min, mostly full SW
    scan)                                           # alone: it is timed
      THREADS=1 $L scan $OUT/bin $OUT/scan.tsv 500:200,100:200 ;;
    annot)                                          # KOfam and dbCAN labels of 21 gut genomes
      $PY $H/annot/genomes.py stats; $L best $OUT/annot; $PY $H/annot/genomes.py score
      $L collapse $OUT/bin; $PY $H/annot/collapse.py; $PY $H/annot/sizes.py
      $PY $H/test_bundle.py; $PY $H/bundle_check.py
      $PY $H/annot/heldout.py align; $PY $H/annot/heldout.py run   # 25 min, mostly hmmsearch
      K=$HOME/Documents/kfpfam/filter/target/release/filter; E=$OUT/annot_heldout   # reach.py: kfpfam from the bundle's files; needs orfs.faa, the windows
      FILTER_MIN_GAPPED=21.25 $K join $E/b0.7/kfpfam/prep $E/b0.7/kfpfam/index.tidx 17 16.5 400000 $E/orfs.faa $E/pairs.tsv; $PY $H/annot/reach.py ;;
    tables)
      $PY $H/sizing.py; $PY $H/check_truth.py; $PY $H/analyze.py; $PY $H/stages.py; $PY $H/zymo.py; $PY $H/members.py; $PY $H/extrapolate.py; $PY $H/loo.py; $PY $H/scancheck.py ;;
    *) echo "unknown step $s" >&2; exit 1 ;;
  esac
done
