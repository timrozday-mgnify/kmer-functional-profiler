# species-benchmark

The species model's strain hold-out benchmark (plan: phase 11, step 10). Each sample
mixes non-representative genomes of MGnify catalogue species. Those genomes are left out
of the species index, so every sample strain is new to it, as a strain in a real sample
would be. Scores cover:

- species detection and abundance;
- unit and Pfam presence, observed against updated by the species' carriage;
- carriage on accessory units;
- the function × species table;
- the lineage model prototype (phase 11, step 12): carriage of accessory units, and
  whether two strains of one species come out as two lineages.

```text
METADATA            genomes-all_metadata.tsv of the catalogue (storeDir)
PICK                bench.py uhgg-pick -> samples.tsv (sample, genome, species, depth),
                    exclude.txt (held-out genomes), species.txt (+ distractors); samples
                    after --replicates (--two_strain_replicates) hold two genomes per species
FETCH_CATALOGUE     fetch.sh species: each species' genome/ (.faa, .fna) and pan-genome/
                    (Rtab, pan-genome.fna) in the FTP layout
FETCH_GENOMES       fetch.sh genomes: the held-out genomes' DNA (their GFF's FASTA section)
SIMULATE            InSilicoSeq --coverage_file: each genome at its depth (log-uniform)
INDEX / --index     the index profiled against, built from --members (+ --pfam) or given
PROFILE             query -> profiles/sample<N>.tsv
SPECIES_INDEX       species-index --catalogue --exclude, per variant (default, noshrink:
                    --alpha 1e-9, nocompl: --no-completeness)
SPECIES_FIT         species per sample and arm (species, species_nobg: --background-prior 0,
                    species_noshrink, species_nocompl)
LINEAGE_FIT         lineage per sample on the default species index (arm lineage)
REPS, ANNOTATE_REPS the representatives' proteins as a genome set, annotated (G1's reference)
REPS_SPECIES_INDEX  species-index --genomes on it: kfp-prior's carriage table
GENOME_FIT          genomes (G1) and kfp-prior update (B1) per sample
SYLPH_DB, SYLPH     sylph on the representatives' DNA
SCORE               bench.py uhgg-score -> species_scores.tsv (long: sample, arm, metric,
                    value), species_calibration.tsv (zero-hit units by predicted presence)
```

Metrics (`metric` column):

- *Detection:* `purity`, `completeness`, `f1`, `l1`, `spearman`. Arms are matched by
  species id (the representative); G1 reports a representative genome, which is the same id.
- *Presence:* `{unit,pfam}_{bin}_{completeness,purity,n_pred}_{observed,updated}`, with
  `bin` the carrier's depth (`all`, or a range). Unit truth is the held-out genomes'
  units, from the Panaroo families present in their `gene_presence_absence.Rtab` column,
  in `held_out.parquet`. Pfam comes from the index's labels.
- *Carriage:* `carriage_auc`, `carriage_logloss` and `prevalence_logloss`, over the true
  species' accessory units (0.1 < q < 0.9) that had no hits. This is where the prior
  carries the call: the gain is `prevalence_logloss` − `carriage_logloss`.
- *Function × species:* `ft_l1`, `ft_f1`, `ft_right`, `ft_unclassified` at rank species.
- *Lineages* (arm `lineage`): detection as above, a species' lineages summed; carriage
  from its deepest lineage; `two_strain_species` and `two_strain_resolved`, the share of
  species with two strains (both ≥ 1×) given ≥ 2 lineages.

Carriage metrics count one-strain species only: a two-strain species has no single
carriage truth.

Gate (plan, phase 11): species F1 and L1 within 0.05 of sylph's. At 0.1–1× depth, Pfam
completeness should rise by ≥ 10 points (updated against observed) at ≤ 2 points of
purity. Calibration error should be ≤ 0.05 on zero-hit units, and better than kfp-prior's.
The lineage model is ported only if its `carriage_logloss` is ≥ 10% below the species
arm's, its `l1` no worse, and `two_strain_resolved` ≥ 0.8 (plan, Species model, step 12).
The lineage arm needs the `phylo` dependency group, which `workflows/setup.sh` installs.

Locally (the mini catalogue in `tests/data/mini_uhgg`, sylph in Docker):

```bash
nextflow run workflows/species-benchmark -profile test,docker
```

On HPC: `hpc/kfp-ablations/species-benchmark/run.sh`. The defaults draw 5 samples of 30
species from the human-gut catalogue v2.0.2, plus 3 two-strain samples, with 300
distractor species. The mini catalogue has no genomes left for two-strain samples, so the
test profile draws none.
