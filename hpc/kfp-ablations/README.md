# kfp-ablations: the phase-7 ablation runs on Slurm

One directory per run of `workflows/fmh-benchmark/README.md`, "Ablations". Every path is
relative, for this layout:

```
<base>/kmer-functional-profiler/                 the repo, after workflows/setup.sh
<base>/kmer-functional-profiler_runs/            <runs>
    singularity_cache/
    mgnify-proteins-subset/results/lin100-dense/  members.parquet, pfam.parquet (1 in 100)
    kfp-ablations/                               this directory
        site.config  common.sh                   shared by every run
        fmh-benchmark-data/                      shared storeDir, created by 00-dbs
        <run>/run.sh                             results/, work/ and the head job's log land here
```

Submit each run from its own directory (`run.sh` checks this):

```bash
cd kfp-ablations/00-dbs && sbatch run.sh     # first, alone: Zenodo inputs, Pfam, hostile, host, decoy
# when it has finished, all at once:
for r in floor-novaseq floor-miseq floor-perfect floor-1.65M floor-26.4M index reads-novaseq reads-miseq host; do
    (cd kfp-ablations/$r && sbatch run.sh)
done
# when floor-novaseq has finished:
cd kfp-ablations/floor-query-cost && sbatch run.sh
# release tiers (plan, phase 7, step 21), then their query cost when tiers-novaseq has finished:
cd kfp-ablations/tiers-novaseq && sbatch run.sh
cd kfp-ablations/tiers-query-cost && sbatch run.sh
# aai survival model (plan, phase 7, step 25), independent of the rest; no 00-dbs needed.
# It is on the dev-aai branch, not main: check that branch out in the repo first.
cd kfp-ablations/aai-model && sbatch run.sh
```

The ablation configs are read from the repo (`workflows/fmh-benchmark/ablations/`), so a
`git pull` there changes them. Re-submitting a run resumes it (`-resume`). Head jobs ask for
7 days; lower `--time` in `run.sh` if your partition caps it.
