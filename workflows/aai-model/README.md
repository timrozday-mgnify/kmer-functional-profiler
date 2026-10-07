# aai-model

Fits the k-mer survival model that `aai` inverts (plan, phase 7, steps 23 and 25). Real
proteins keep more exact k-mers than a^k at identity a, because conserved sites cluster:
rates vary along the sequence in regions. The model (`kmer_functional_profiler.survival`)
draws each site's rate from `categories` gamma categories of shape `shape`, and redraws it
with probability 1/`region` at each site, a Markov chain along the sequence. Identity,
window survival and the co-survival of windows j apart then follow exactly, from small
matrix products tabulated once per model.

The parameters come from protein pairs, not from a benchmark. A pair is a protein P and an
MGnify90 cluster C it aligns to (both coverages >= 0.8). Survival is the share of P's k-mers
found in the union of C's members' k-mers, over C's average member's k-mers (`pin_sum`):
what the query measures. Identity is DIAMOND's, P against C's representative. When P is a
member of C, it is left out of C. Survival against identity says how far survival exceeds
a^k; co-survival of P's windows j apart (`both_<j>`) separates a strong rate contrast over
short regions from a weak one over long regions.

```text
QUERIES   members table -> n random members (P; standard residues only) in chunks, representatives
DB, ALIGN DIAMOND blastp of P against the representatives        (--chunks jobs)
PAIRS     best HSP per (P, C), coverage filter, <= per_bin pairs per 0.02 of identity
SURVIVAL  per pair: P's k-mers in C's union, pin_sum, co-survival by lag (every k-mer)
FIT       aai_model.json (least squares over shape, region), model_strata.tsv
          (refitted by identity band and cluster size)
```

## Run

Needs Java 17+, Nextflow, Singularity or Docker (DIAMOND), and the repo's venv
(`bash workflows/setup.sh`). On Slurm, over the
1-in-100 subset used by the ablations (`hpc/kfp-ablations/aai-model/run.sh` does this):

```bash
nextflow run workflows/aai-model -profile slurm \
    --members mgnify-proteins-subset/results/lin100-dense/members.parquet
```

`nextflow run workflows/aai-model -stub -profile test --members <any file>` checks the wiring.

## Use

Attach the fit to an index; its queries then estimate `aai` under it:

```bash
kmer-functional-profiler aai-model <index> results/aai_model.json
```

Check `model_strata.tsv` first. If shape or region differ much between identity bands or
cluster sizes, one parameter pair is not enough (plan, step 25).
