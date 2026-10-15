# spike-benchmark

Strains of known identity and depth, and withheld siblings, spiked into reads (optionally a
real background), then profiled by kfp and followed up by `workflows/pileup`'s `PILEUP`.
Exact truth for three checks (plan: phase 13, steps 3-6 and its gate; phase 14, step 5 and
its gate):

- **the pileups:** read-assignment precision, and whether decoys catch a withheld sibling's
  reads;
- **`aai` / `ani`:** the profile's similarity against the strain's identity to the nearest
  member, and step 6's check (`aai` against the reads' alignment identity) with the truth
  beside it;
- **nucleotide mode:** a nucleotide index (bwa-mem2 follow-up) against KMA and against kfp's
  amino-acid query on the translated targets.

```text
FAMILIES       (without --targets) spike.py families: --families seeds, each 3 units of 3
               members (paralogs at 80-90%, members at 97%), protein or back-translated CDS
SIBLINGS       spike.py siblings: 10% of units get sib_<unit>, a member mutated to 80-90%;
               siblings join the index but are withheld from the pileup's reference
SIMULATE       per seed: 30% of units carry a strain (a member at 100/98/95/90%), each
               sibling one; lognormal depth; paired reads with 0.2% substitutions, between
               random flanks; optional background pairs appended
INDEX_TARGETS  kfp index --targets --groups (+ aai-model with --aai_model); with
               --molecule nucleotide also the translated targets' index (TRANSLATE)
PILEUP         workflows/pileup: QUERY, REFERENCE (--withhold the siblings), DIAMOND or
               BWAMEM2, CALL
QUERY_AA, KMA  nucleotide only: the amino-acid arm and KMA (1.6.18)
SCORE, SUMMARY spike.py score / summary -> <seed>/scores.tsv, <seed>/units.tsv, summary.tsv
```

`summary.tsv` has one row per arm (`kfp`, `kfp_aa`, `kma`), means over samples; the columns
are documented in `spike.py`'s docstring. The ones the gates read:

| Column | Gate |
| --- | --- |
| `read_precision` | phase 13: >= 0.95 |
| `sibling_caught` (`_parent_target`, `_parent_neighbour`, `_parent_none`) | phase 13: >= 0.8; split by the role the profile gave the sibling's parent unit (a target cannot catch them) |
| `aln_gap_<bin>`, `aln_cover` | phase 13, step 6: a gap above 0.01 in a bin with >= 100 units, or a cover below 0.8, is a phase-7 `aai` finding |
| `f1`, `l1` (`kfp` against `kma`) | phase 14: within 0.05 of KMA's |
| `sim_bias_<bin>`, `sim_sharp_bias_<bin>` | phase 14: `ani` \|median bias\| <= 0.01 over 0.9-1.0 |

`trace.tsv` gives the cost (phase 13, step 3's rule: the aligner against `PILEUP:QUERY`; phase
14: kfp against KMA), and `--query_stats true` publishes each query's `stats.json`
(`hits_per_pair`: phase 13, step 0).

```bash
nextflow run workflows/spike-benchmark -profile test,docker              # fixture, protein
nextflow run workflows/spike-benchmark -profile test_nucleotide,docker   # fixture, nucleotide
nextflow run workflows/spike-benchmark -profile local,docker             # laptop, ~5 min
nextflow run workflows/spike-benchmark -profile local_nucleotide,docker
nextflow run workflows/spike-benchmark -profile docker --targets genes.fna --molecule nucleotide \
    --background r1.fq.gz,r2.fq.gz --index_args '--k 21 --t-base 0.02'
```

The `local` profiles need the MGnify sample (`scripts/fetch_mgnify_sample.py`). HPC runs:
`hpc/kfp-ablations/spike-protein` and `spike-nucleotide`.

Limits of the simulation: substitutions are uniform over sites (no clustered variation, no
indels); nucleotide strains avoid stop codons but are not biased to synonymous sites, so at
90% nucleotide identity the protein is ~80% identical, which is hard on the amino-acid arm;
siblings are synthetic neighbours, not MGnify's.
