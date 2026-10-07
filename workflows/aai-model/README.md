# aai-model

Fits the k-mer survival model that `aai` inverts (plan, phase 7, steps 23, 25, 33 and 34). Real
proteins keep more exact k-mers than a^k at identity a, because conserved sites cluster.
The model (`kmer_functional_profiler.survival`, Markov-beta) splits a sequence into regions,
a Markov chain that starts a new region with probability 1/`region` per site. Each region's
identity scatters around the pair's identity as Beta(a φ, (1 − a) φ), φ = `concentration`.
`ends` is the share of a unit's windows inside the region a homolog aligns to. A unit of n
members keeps more than one member does: its survival is 1 − (1 − S)^m, with m = n^(g₀ + g₁(a − 0.8))
effective members (`union`, `union_slope`). Identity,
window survival and the co-survival of windows j apart then follow exactly, from small
matrix products tabulated once per model.

The parameters come from protein pairs, not from a benchmark. A pair is a protein P and an
MGnify90 cluster C it aligns to (both coverages >= 0.8). Survival is the share of P's
windows found in the union of C's members' k-mers, which is what the query measures.
Identity is DIAMOND's, P against its **nearest member of C**. Identity to C's representative
(step 25) made survival in large clusters far exceed what it predicts, because P has closer
members. When P is a member of C, it is left out of C. Co-survival of P's windows j apart
(`both_<j>`) separates a strong contrast over short regions from a weak one over long regions.

```text
QUERIES     members table -> n random members (P; standard residues only) in chunks, representatives
REP_ALIGN   DIAMOND blastp of P against the representatives          (--chunks jobs)
CANDIDATES  best HSP per (P, C), coverage filter, <= per_bin_candidates per 0.02 of rep identity
MEMBER_DB   the candidates' clusters' members (FASTA, member -> cluster), their P in chunks
NEAREST_ALIGN  DIAMOND blastp of P against those members             (--chunks jobs)
PAIRS       per candidate, the member of C with the highest identity (coverage filter, not P),
            <= per_bin pairs per 0.02 of that identity; identity_rep kept for comparison
SURVIVAL    per pair: P's windows in C's union, pin_sum, co-survival by lag (every k-mer)
FIT         aai_model.json, two stages on the bins' means: phi, region and ends from one-member
            pairs (bins of 0.01 identity), then union and union_slope from multi-member pairs
            (bins of identity x cluster size); model_strata.tsv (one-member model by identity
            band, union term by cluster size, and the overall model's error per size)
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

`results/pairs_stats.json` counts the candidates, those that kept a member hit (overall and
by cluster size) and the pairs kept. In pass 2, a P in a large cluster can fill
`--max-target-seqs` with that cluster's members and lose the hit to another candidate
cluster. If many candidates are lost, raise it in `--diamond_nearest_args`.

## Use

Attach the fit to an index; its queries then estimate `aai` under it:

```bash
kmer-functional-profiler aai-model <index> results/aai_model.json
```

Check `model_strata.tsv` first. The `members 1, identity …` rows should agree on
`concentration` and `region`. For the `members 2-3` … `members 101+` rows, `rmse_overall`
(the shipped model's binned error at that size) should be ≤ 0.02 (plan, step 34).
`ends` comes from pairs at >= 0.8 coverage. If the benchmark's genes cover their units
differently, compare `"ends": 1` on the same profiles (`aai-score --model`).
