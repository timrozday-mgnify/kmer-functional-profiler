# aai-alpha

Fits the shape `alpha` of the k-mer survival model that `aai` inverts (plan, phase 7,
step 23): real proteins keep more exact k-mers than a^k at identity a, because
substitutions cluster in variable regions. Under the regional-rate model a k-mer window
survives with probability S(a) = (1 + k (a^(-1/alpha) - 1))^(-alpha) (a^k as alpha -> inf).

`alpha` comes from protein pairs, not from a benchmark. A pair is a protein P and an
MGnify90 cluster C it aligns to (both coverages >= 0.8). Survival is the share of P's k-mers
found in the union of C's members' k-mers, over C's average member's k-mers (`pin_sum`),
which is what the query measures. Identity is DIAMOND's, P against C's representative. When
P is a member of C, it is left out of C. Pairs are spread evenly over identity 0.6-1.0.

```text
QUERIES   members table -> n random members (P) in chunks, cluster representatives
DB, ALIGN DIAMOND blastp of P against the representatives        (--chunks jobs)
PAIRS     best HSP per (P, C), coverage filter, <= per_bin pairs per 0.02 of identity
SURVIVAL  P's k-mers in C's union and pin_sum per pair, every k-mer (k, alphabet)
FIT       alpha.json (least squares), alpha_strata.tsv (by identity band, cluster size)
```

## Run

Needs Java 17+, Nextflow, Singularity or Docker (DIAMOND), and the repo's venv
(`bash workflows/setup.sh`). On Slurm, over the 1-in-100 subset used by the ablations
(`hpc/kfp-ablations/aai-alpha/run.sh` does this):

```bash
nextflow run workflows/aai-alpha -profile slurm \
    --members mgnify-proteins-subset/results/lin100-dense/members.parquet
```

`nextflow run workflows/aai-alpha -stub -profile test --members <any file>` checks the wiring.

## Use

Attach the fit to an index; its queries then estimate `aai` under it:

```bash
kmer-functional-profiler aai-model <index> results/alpha.json
```

Check `alpha_strata.tsv` first. If alpha differs much between identity bands or cluster
sizes, one parameter is not enough (plan, step 23).
