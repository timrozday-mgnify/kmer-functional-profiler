# presence-calibration

`kfp query --background-share` calibrated at full-index scale (plan, phase 13 / phase 7:
presence at a false discovery rate). The spike benchmark's 2,000-family index is too small for
the share's cap to bind, and real reads have no truth, so the spikes here are strains of the
index's own units in a real background, with two negatives of known truth.

```text
PICK      calibrate.py pick: --units units drawn uniformly from the index (representative
          200-600 aa), --members sequences each from the store; units.tsv (n_members, m_g,
          uniq: the mean share of a unit's sampled k-mers it holds alone)
SIMULATE  spike-benchmark's spike.py simulate: every picked unit carries a strain (a member
          at --identities), lognormal depth (~0.05-20x), appended to --background
REVERSE   calibrate.py reverse: the background reversed, a null sample (no coding homology)
QUERY     every sample (spiked seeds, unspiked background, null) at every --shares value
SCORE     calibrate.py score: per spiked sample, share and q (0.01, 0.05, 0.1): recall (by
          depth, identity, uniqueness), new_false (called here, not unspiked, not spiked),
          false_share, null_called, unspiked_called
DECIDE    calibrate.py decide -> summary.tsv, decision.tsv (the largest share passing the
          rule at --q in every seed), scores.tsv
```

The rule: `false_share` <= `--max_false` (0.05) and `null_called` <= `--max_null` (0.01) x
`unspiked_called`, at `present_q` <= `--q` (0.05), in every seed. `new_false` counts every new
call that is not a spiked unit, so `false_share` is an upper bound (it includes calls that
shift with the extra reads, not only the spikes' relatives).

```bash
nextflow run workflows/presence-calibration -profile test   # fixture: plumbing only
nextflow run workflows/presence-calibration --index /abs/index --store '/abs/members.bucket*.parquet' \
    --background r1.fq.gz,r2.fq.gz
```

On HPC: `hpc/kfp-ablations/presence-calibration/run.sh` (needs `release-tiers` and
`pileup-cost`'s reads).
