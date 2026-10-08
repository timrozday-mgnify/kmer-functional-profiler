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
# standard tier at ~120 GB and frame mode on the release tiers (plan, phase 7, step 28), all at once:
for r in standard-novaseq frames-novaseq frames-miseq; do (cd kfp-ablations/$r && sbatch run.sh); done
# then, each when its NovaSeq run has finished:
cd kfp-ablations/standard-query-cost && sbatch run.sh
cd kfp-ablations/frames-query-cost && sbatch run.sh   # four frame modes, one after another
# EM convergence fix (plan, phase 7, step 29), once tiers-novaseq's indexes exist:
cd kfp-ablations/em-query-cost && sbatch run.sh
# EM convergence from 5 starts on two deep soil metagenomes (step 30), as em-query-cost:
cd kfp-ablations/em-starts-soil && sbatch run.sh
# aai survival model (plan, phase 7, step 25), independent of the rest; no 00-dbs needed:
cd kfp-ablations/aai-model && sbatch run.sh
# aai under the step-34 model against nearest-member identity (plan, phase 7, step 34), once
# floor-novaseq and floor-1.65M have finished: their published profiles rescored, arms
# '+step34', '+step34_nounion', '+step34_ends1' (no index rebuilt, no query rerun):
cd kfp-ablations/aai-rescore && sbatch run.sh
# the three release tiers at full scale (slim, standard, large), independent of the rest;
# reads the whole release (RELEASE=<local mirror> to use one), ~1 TB of buckets in results-*/:
cd kfp-ablations/release-tiers && sbatch run.sh
# abundance under divergence (plan, phase 7, step 35): finished runs' indexes queried again
# with the shipped coverage_zi and its interval (no index rebuilt); after floor-novaseq,
# floor-1.65M and tiers-novaseq respectively:
for r in abundance-novaseq abundance-1.65M abundance-tiers; do (cd kfp-ablations/$r && sbatch run.sh); done
```

The ablation configs are read from the repo (`workflows/fmh-benchmark/ablations/`), so a
`git pull` there changes them. Re-submitting a run resumes it (`-resume`). Head jobs ask for
7 days; lower `--time` in `run.sh` if your partition caps it.
