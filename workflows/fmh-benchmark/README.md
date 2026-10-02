# fmh-benchmark

KO detection by kmer-functional-profiler indexes against fmh-funprofiler's KO sketches on
simulated metagenomes, following the fmh-funprofiler paper (Bioinformatics 2024) but with
InSilicoSeq instead of CAMISIM. This was the phase-3 gate (completeness and purity at equal
density); from phase 5 it also runs other tools on the same metagenomes (`--tools`, see
[Other tools](#other-tools)) and scores them the same way, and it scores Pfam as well as
KO (`--labels`, see [Pfam](#pfam)), Pfam being the primary label from phase 5 on.

```text
FETCH            Zenodo 10055954 (CC-BY): KEGG genomes (9.1 GB zip), KEGG proteins (3.2 GB),
                 gene -> KO table, KO sketches. Downloaded once into --data_dir.
MEMBERS          proteins grouped by KO -> members.parquet
PFAM_DB          Pfam-A.hmm from --pfam_url (or --pfam_hmm), once into --db_dir
GENOME_PROTEINS  every genome's proteins, in --pfam_chunks FASTA files
PFAM_ANNOTATE, PFAM_DOMAINS
                 hmmsearch --cut_ga per chunk ->
                 pfam/domains.parquet (gene, Pfam, envelope) and the domains as members of
                 one unit per Pfam
INDEX            kmer-functional-profiler index, one per --indexes entry and label (KO
                 units as named; Pfam units as pfam_<name>)
IMPORT_SKETCHES  the KO sketches as a sourmash-compatible index ("fmh_compat")
SAMPLE           --n_genomes random genomes per replicate (seed = replicate number)
SIMULATE         iss generate, lognormal abundances (--iss_mode perfect: no read errors)
TRUTH            reads mapped back with minimap2 (mappy); a KO is present if a read
                 overlaps one of its genes (the paper's rule, on CAMISIM's alignments); its
                 depth is aligned bases / gene length, summed over its genes. The same over
                 Pfam domains (their nt span on the genome, either strand) -> truth/seedN_pfam.csv
MGNIFY_INDEX, MGNIFY_REPS, MGNIFY_DB, MGNIFY_ANNOTATE, MGNIFY_GENES
                 per --mgnify_indexes entry: the index (built here from members, or given),
                 its cluster representatives, and every genome protein searched against them
                 (DIAMOND blastp) -> mgnify/<name>_gene_units.parquet (MGnify90-level truth)
PROFILE          kmer-functional-profiler query, every metagenome x every index; profiles of
                 --mgnify_indexes (Pfam-labelled units) are summed per Pfam
AAI_SCORE, AAI_SUMMARY
                 unit profiles of --mgnify_indexes against that truth: detection at the 90%
                 level, nearest-cluster recall beyond it, aai and aai_naive against alignment
                 identity -> aai_summary.tsv, aai_scores.tsv
AAI_CALIBRATE    per MGnify index and arm, over all seeds: an aai -> identity map fitted on
                 half the clusters and scored on the other half ->
                 calibration/<name>[~arm].json, calibration/<name>[~arm]_scores.tsv
DETECTED         per KO gather keeps, true or false: its evidence and where its hit k-mers
                 come from (holders, hits, in the sample genomes or not) -> detected/,
                 detected.tsv (all samples)
DIAMOND_DB, DIAMOND, FMH_FUNPROFILER, KMERMAID_MODEL, KMERMAID, HUMANN_DB, METAPHLAN_DB,
HUMANN, HUMANN4
                 other tools (--tools), in containers; databases built once into --db_dir
TOOL_PROFILE     each tool's output as a profile (name, evidence, abundance) -> profiles/
SCORE, SUMMARY   purity, completeness, completeness of the 25% least-covered true KOs,
                 base-weighted completeness per count (kmers_hit; kmers_unique after gather;
                 kmers_wta, kmers_ufirst after winner-take-all, uniqueness-first)
                 and --min_hits value; abundance of the detected KOs against truth depth
                 (spearman_tp, l1) -> summary.tsv (mean, sd),
                 scores.tsv; profiles/, kmers/ (tier-2 hits per KO and k-mer) and truth/ keep
                 the per-sample tables
```

`run.json` in the output directory records what produced it: a one-line description
(read model, sizes, indexes), all parameters, the code commit and whether the checkout had
uncommitted changes, the command line and the run's status. Runs with different settings
(e.g. `--iss_mode perfect`) should use their own `--outdir`: files from an earlier run are
overwritten or left in place, not removed.

The `fmh_compat` row reproduces fmh-funprofiler's KO calls: its hits equal the overlaps
`sourmash prefetch` reports (tested in `tests/python/test_compat.py`), and a KO counts as
detected at one shared hash, as with fmh-funprofiler's default `--threshold-bp 1000` at
scaled 1000. The paper's reference numbers (CAMISIM, wgsim, 64 genomes, 1 Gbp, k = 11):
purity 0.98, completeness 0.61 at 0% error; they are not directly comparable to runs
here, because the simulator differs.

## Setup

Needs Java 17+, Nextflow, [uv](https://docs.astral.sh/uv/) and a Rust toolchain; for the
other tools, Docker (local) or Singularity/Apptainer (HPC).

```bash
bash workflows/setup.sh                            # .venv with the package, InSilicoSeq and mappy
bash workflows/fmh-benchmark/containers/build.sh   # kMermaid, HUMAnN 3.9, HUMAnN 4 images (docker)
```

Our steps run on the host in `.venv`; only the tool steps use containers: public
biocontainers for DIAMOND and fmh-funprofiler, and three images built here. kMermaid has
no package. HUMAnN 3.9's biocontainer resolves bowtie2 2.2.3, whose `bowtie2-build` has no
`--threads`, so HUMAnN fails at its index build; `containers/humann3` adds bowtie2 2.5.5.
HUMAnN 4 (4.0.0.alpha.2) is only on GitHub; `containers/humann4` installs it at a fixed
commit on top of that image, which already has the MetaPhlAn (4.1.1) and DIAMOND (2.0.15)
versions it pins.

## Run

Test profile (tiny fake inputs in `tests/data/mini_fmh`, a few minutes; DIAMOND,
fmh-funprofiler and kMermaid, not HUMAnN, whose databases are ~45 GB for 3.9 and ~70 GB
for 4; `-stub --tools humann,humann4` checks their wiring):

```bash
nextflow run workflows/fmh-benchmark -profile test,docker
```

With Docker Desktop, run from a directory it shares (e.g. under your home; not `/tmp`).
Without containers, `--tools ''` runs our indexes only, as before phase 5.

### HPC

1. Environment and images (login node, internet):

   ```bash
   bash workflows/setup.sh
   bash workflows/fmh-benchmark/containers/pull.sh /shared/singularity-cache
   ```

   The built images: on a machine with Docker, `bash workflows/fmh-benchmark/containers/build.sh --sif`
   writes `kfp-kmermaid.sif`, `kfp-humann3.sif` and `kfp-humann4.sif` (or `.tar` files,
   without Singularity there, to convert on the cluster with
   `singularity build kfp-<name>.sif docker-archive://kfp-<name>.tar`); copy them to shared
   storage. Or push them to a registry with `build.sh --push <registry>`.
2. Site config: `cp workflows/fmh-benchmark/hpc.example.config workflows/fmh-benchmark/hpc.config`
   and fill in the paths (data, results, Singularity cache, the three built images),
   partition and account. If compute nodes have no internet, uncomment the line that runs
   FETCH and the HUMAnN database downloads on the head job's node.
3. Run the head job:

   ```bash
   sbatch workflows/fmh-benchmark/run_hpc.sh
   ```

   Extra arguments go to Nextflow, e.g. `sbatch workflows/fmh-benchmark/run_hpc.sh --tools diamond,fmh_funprofiler`.
   With `-resume` (always on in `run_hpc.sh`) an earlier run's indexes, samples and truth
   are reused if `data_dir` and the work directory are the same.
4. Cost per tool: `.venv/bin/python workflows/fmh-benchmark/bench.py cost <outdir>/trace.tsv`
   writes `cost.tsv` (tasks, mean wall and CPU hours, peak RSS per step and index/tool).

First run downloads: ~13 GB of Zenodo inputs into `--data_dir`; into `--db_dir`, ~45 GB of
HUMAnN 3.9 databases (ChocoPhlAn, UniRef90, utility mapping, MetaPhlAn vJun23) and ~70 GB
of HUMAnN 4 ones (ChocoPhlAn v4 alpha, 45 GB compressed; EC-filtered UniRef90; utility
mapping; MetaPhlAn vOct22_202403).
To download them ahead of (or apart from) the benchmark, e.g. on a node with internet, run
`sbatch workflows/fmh-benchmark/run_hpc.sh --dbs_only --tools humann,humann4,diamond,kmermaid --dbs pfam,hostile,host,decoy`.

| Parameter | Default | Meaning |
| --- | --- | --- |
| `--data_dir` | `fmh-benchmark-data` | Where the Zenodo inputs are stored (reused across runs) |
| `--sketches` | `KOs_sketched_scaled_1000.sig.zip` | KO sketches from the same Zenodo record |
| `--ksize` | `11` | Protein k of the sketches used |
| `--n_genomes` | `64` | Genomes per metagenome |
| `--replicates` | `10` | Metagenomes (seeds 1..N) |
| `--n_reads` | `6600000` | InSilicoSeq reads, both mates (~1 Gbp at 151 bp) |
| `--iss_model` | `novaseq` | InSilicoSeq error model: `hiseq`, `novaseq` or `miseq` |
| `--iss_mode` | `kde` | `perfect`: error-free 151 bp reads, to see which false positives come from read errors (run through `bench.py iss`, which patches iss 2.0.1's perfect model: it otherwise exits 0 without output, and makes 125 bp reads) |
| `--draws` | `100` | Posterior draws for 95% intervals on `coverage_zi` / `abundance_zi` and ambiguity groups (0: none) |
| `--min_hits` | `1,2` | Distinct k-mers for a KO to count as detected; each value is scored from the same profiles |
| `--diamond_min_hits` | `1,2,3,5,10,20,50,100` | Read pairs for a KO to count as detected by DIAMOND, in place of `--min_hits`: DIAMOND's purity/completeness curve, to compare at matched purity |
| `--indexes` | four configs (see `nextflow.config`) | `[name:, args:]` maps of `index` options |
| `--tools` | `diamond,fmh_funprofiler,kmermaid,humann,humann4` | Other tools to run and score; `''` for none |
| `--db_dir` | `--data_dir` | Where tool databases are built once (DIAMOND, kMermaid model, HUMAnN) |
| `--dbs_only` | `false` | Only fetch the Zenodo inputs and build the `--tools`' and `--dbs`' databases, then stop; later runs with the same `--data_dir`/`--db_dir` reuse them |
| `--dbs` | `''` | With `--dbs_only`, also any of `pfam` (Pfam-A HMMs), `hostile` (hostile index), `host` (host genome), `decoy` (decoy proteome) |
| `--sketch_scaled` | `1000` | Scaled of `--sketches`, for fmh-funprofiler |
| `--diamond_args` | `''` | Extra `diamond blastx` options, e.g. `--sensitive` |
| `--kmermaid_max_members` | `50` | Proteins sampled per KO to train kMermaid |
| `--kmermaid_container` | `kfp-kmermaid:edcb4ed` | kMermaid image: docker tag, `.sif` path or `docker://` URI |
| `--humann3_container` | `kfp-humann3:3.9-bt2.5.5` | HUMAnN 3.9 image, as above |
| `--humann4_container` | `kfp-humann4:e07b3a3` | HUMAnN 4 image, as above |
| `--metaphlan_index` | `mpa_vJun23_CHOCOPhlAnSGB_202403` | MetaPhlAn database for HUMAnN 3.9 (MetaPhlAn 4.1.1) |
| `--humann4_metaphlan_index` | `mpa_vOct22_CHOCOPhlAnSGB_202403` | MetaPhlAn database for HUMAnN 4 (the one it checks for) |
| `--labels` | `ko,pfam` | Labels scored (see [Pfam](#pfam)) |
| `--pfam_url`, `--pfam_hmm`, `--pfam_threshold`, `--pfam_chunks` | current release, `''`, `--cut_ga`, `64` | Pfam-A source (or a local file), hmmsearch inclusion, hmmsearch jobs |
| `--mgnify_indexes` | `[]` | `[name:, path:, members:, pfam:, args:]` maps of MGnify90 indexes (config file only; see [MGnify90-level truth](#mgnify90-level-truth-unit-resolution-and-containment-aai)) |
| `--mgnify_min_id`, `--mgnify_min_cov` | `0.9`, `0.8` | A gene is in its nearest cluster at this identity and coverage |
| `--mgnify_diamond_args` | `--sensitive --max-target-seqs 25 --id 50 --query-cover 50` | Which near hits are kept |
| `--query_arms` | one plain arm | `[name:, args:, reads:, mask:, decoy:]` maps (config file only; see [Ablations](#ablations)) |
| `--host_fractions` | `0` | Host share of read pairs, comma-separated, e.g. `0,0.5,0.9,0.99` |
| `--host_genome_url`, `--host_extra_url` | T2T-CHM13v2.0; rCRS chrM + PhiX174 (NCBI) | Host genome for spike-in reads and masks |
| `--decoy_proteome_url` | UniProt UP000005640 (human) | Proteome of the decoy index, one unit per protein |
| `--hostile_index` | `human-t2t-hla` | hostile's bowtie2 index |
| `--fastp_args` | `''` | Extra fastp options (always `--trim_poly_g`) |

Defaults compare, at k = 11: fmh-funprofiler's sketches (scaled 1000); our index at the
same base rate without and with the per-KO floor (`--n-min 8`); and our index at 10x
density (`--t-base 0.01`), the plan's "scaled = 100" baseline, alone and with a dense tier
at 1 in 10 (`--t-dense 0.1`: detection unchanged, EM abundances fitted on 10x more
k-mers). `index_*/meta.json` in the output records each index's size (`dense_bytes` for the
dense tier); `trace.tsv` records each PROFILE task's runtime.

With `-resume`, INDEX, IMPORT_SKETCHES and PROFILE rerun when the Python package changes,
and TRUTH and SCORE when `bench.py` does; Rust kernel changes need a fresh run.

Abundance is scored per row (`abundance` column) with an estimate that goes with its count:
`coverage` (hits per kept k-mer) with `kmers_hit`; `coverage_em` (EM over the units gather
keeps), `coverage_zi` (zero-inflated EM: coverage of the k-mers present) and `coverage_zib`
(zero-inflated with an empirical-Bayes prior on the present fraction) and `coverage_zip`
(zero-inflated with each k-mer's presence proportional to its in-KO frequency *p_in*) and
`abundance_zi` (`coverage_zi` times the gene copies present) with `kmers_unique`; and the hits each one-pass rule assigns per kept k-mer with its count:
`coverage_wta` with `kmers_wta` (each hit k-mer to the holding KO with the highest
containment, as sylph) and `coverage_ufirst` with `kmers_ufirst` (to the holding KO with the
highest Σ 1 / KOs-per-hit-k-mer, scaled by 1 / *t_g*). `spearman_tp` is the rank correlation with truth depth over true positives;
`l1` is the L1 distance between relative abundances over all true and detected KOs (0 is
exact, 2 is disjoint). With `--draws` > 0, `ci_cover` is the share of true
positives whose posterior interval holds truth depth (on the estimate's scale) and
`ci_width` the median log(hi / lo); for `abundance_zi`, `fp_grouped` is the share of false
positives placed in an ambiguity group (shared-evidence) with a true KO, and `group_cover`
the share of groups whose interval holds the members' true total. Also for `abundance_zi`,
`prob_tp` / `prob_fp` are the mean `present_prob` of true and false positives and
`flag_tp` / `flag_fp` the share of each below 0.5.

`detected.tsv` has one row per KO gather keeps per sample and index: `tp`, its counts,
`present_prob`, `own_evidence`, `ambiguity_group`, and over its tier-2 hit k-mers the
median `holders` (index KOs holding the k-mer), `hits_max` (most hits on one k-mer) and
`in_genome` (share present in the six-frame translation of the sample genomes; null for
`fmh_compat`, which hashes with sourmash). A false positive's k-mers that are in the
genomes are real sequence (another gene or KO); those that are not come from read errors.
Rerunning with `--iss_mode perfect` shows how many false positives errors cause.

## Pfam

With `--labels ko,pfam` (the default), every gene of every genome in the KEGG extraction is
annotated with Pfam-A (`hmmsearch --cut_ga`, envelope coordinates). Pfam truth uses the
domains, not the genes: a Pfam is present if a read overlaps one of its domains on the
genome, and its depth is the sum over those domains of aligned bases / domain length
(Benchmark labels in the plan). Two kinds of profile are scored against it:

- *Pfam units* (`pfam_<name>`, one per `--indexes` entry): each domain is a member of its
  Pfam's unit, the Pfam analogue of the KO index, which isolates the method from the reference.
- *MGnify90 indexes* (`--mgnify_indexes`, e.g. the phase-6 full build): the profile is summed
  per Pfam through the index's `unit_pfam.parquet`, each unit counting for each of its Pfams,
  the shipped configuration. Pass them in a config file:
  `params.mgnify_indexes = [[name: 'mgnify_full', path: '/shared/full-build/index']]`.

`summary.tsv` and `scores.tsv` gain a `label` column (`ko`, `pfam`). Set `--pfam_url` to the
Pfam release MGnify's `mgy_proteins_pfam` used, so the labels match (default: current
release). Other tools are scored on KO only, for now.

## MGnify90-level truth (unit resolution and containment AAI)

KO and Pfam truth check functions; they cannot check which MGnify90 cluster a read came from,
or the AAI the profile reports for it. For each `--mgnify_indexes` entry with `members`,
every protein of the KEGG genomes is searched against the representatives of the index's
clusters (`diamond blastp --sensitive`, up to 25 hits at ≥ 50% identity and ≥ 50% query
coverage: `--mgnify_diamond_args`), giving per gene its *nearest cluster* (best bitscore) and
its near hits, each with identity and coverage. A present gene (reads in the sample, from
`truth/seedN_genes.csv`) is *in* its nearest cluster at the MGnify90 level if the hit has
≥ `--mgnify_min_id` (0.9) identity and ≥ `--mgnify_min_cov` (0.8) query and subject coverage.

`aai_summary.tsv` (means over seeds, per index and arm) and `aai_scores.tsv`:

| Column | Meaning |
| --- | --- |
| `completeness_90` | Clusters holding a present gene at the 90% level that are detected (`kmers_unique` >= 1) |
| `purity_nearest`, `purity_near` | Detected units that are some present gene's nearest cluster, or any of its hits |
| `recall_<lo>` | Present genes whose nearest cluster is detected, by that hit's identity (0.95, 0.9, 0.8, 0.7, 0.5): resolution below 90% |
| `aai_bias_<lo>`, `aai_cover_<lo>`, `aai_spearman`, `aai_n` | Detected units' `aai` against the depth-weighted identity of the present genes they are nearest to; `aai_cover` = share inside `aai_lo`–`aai_hi` |
| `naive_bias_<lo>`, `naive_within05_<lo>`, `naive_spearman` | Every profiled unit a present gene hits (near hits included): `aai_naive` against the best identity of a present gene to it |

`calibration/<name>[~arm].json` is an inverse calibration of `aai` (plan, phase 7, steps
16–17): knots `aai` -> `identity`, the interval widening `widen`, the index's build
parameters it holds for, and the fit's provenance. `_scores.tsv` has the `aai_*` columns
above, raw and calibrated, per seed on the held-out half of the clusters. A map is attached
to an index with `kmer-functional-profiler calibrate-aai INDEX_DIR MAP.json` (the query then
reports calibrated `aai` and keeps `aai_raw`); it is fitted on this benchmark only, so check
it on others before attaching it to an index for real samples (plan, step 17). A test
fixture has too few units for a map (`widen` null; `calibrate-aai` refuses it).

The identity is to the cluster *representative*, which is what `aai` estimates (to the
cluster's consensus; plan, phase 7, steps 7 and 9). Entries take either a built index and the
members it came from, or members (with `pfam` and `args`) to build here:

```groovy
params.mgnify_indexes = [
    [name: 'mgnify_1in100', path: '/shared/cost-nested/1in100/index',
     members: '/shared/cost-nested/members.1in100.parquet'],
]
```

`members` may be a directory of members tables (the full build's MERGE buckets). Entries
with the same `members` (e.g. sparse and dense builds of one subset) share one annotation,
published as `mgnify/<first entry's name>_gene_units.parquet`. Annotating
against a nested subset (1 in 100 or 1 000 clusters) is cheap; against the whole release's
~1.7×10⁹ representatives it is ~10³–10⁴ CPU-hours. `ablations/aai.config` is a template for
this run (`run_ablations.sh aai OUTDIR`). The test profile builds `mgnify_mini` from a
fixture whose 24 clusters' representatives are mini-genome proteins mutated to 100-70%
identity (`tests/data/mini_fmh/mgnify_truth.tsv`); DIAMOND recovers those identities to
±0.01.

## Ablations

Phase 7 measures each change against the phase-6 method on this benchmark (plan, phase 7).
Three runs, each a config in `ablations/` on top of the site's `hpc.config`, with its own
output and work directory, so they can run at once:

```bash
sbatch workflows/fmh-benchmark/run_ablations.sh index /shared/kfp-ablations/index
sbatch workflows/fmh-benchmark/run_ablations.sh reads /shared/kfp-ablations/reads-novaseq
sbatch workflows/fmh-benchmark/run_ablations.sh reads /shared/kfp-ablations/reads-miseq --iss_model miseq
sbatch workflows/fmh-benchmark/run_ablations.sh host  /shared/kfp-ablations/host
```

| Run | Arms | Profiles (10 seeds) |
| --- | --- | --- |
| `index` | Base (k 11, 20 letters, *t_base* 1/1000, *n_min* 8) and one change each: *n_min* 0/4/16, *t_base* 1/100, dense tier 0.05/0.1/0.2, k 9/10/12, Murphy-10 k 13/15, Dayhoff k 17/20. Built per label (KO and Pfam units). | 31 indexes × 10 |
| `reads` | On the base and *t_base* 1/100 indexes: frames stop-free, edges *m* = 15/20/30, all six; quality mask *q* = 10/20/30; fastp then raw; mask × edges; fastp × edges. Once per read model. | 5 indexes × 11 arms × 10 |
| `aai` | Sparse against dense at the MGnify90 level: one nested subset (1 in 100) built sparse and with *t_dense* 0.02/0.05/0.1, against one MGnify90-level truth; plain and edges arms. Detection and `aai_naive` are tier-2 only, so the ladder isolates the dense tier's effect on `aai`, its interval and Pfam abundance | 4 indexes × 2 arms × 10 |
| `host` | 0/50/90/99% host read pairs (T2T-CHM13 + rCRS chrM + PhiX, simulated with the same read model, microbial pairs subsampled to keep depth); no handling, hostile, mask, human-proteome decoy, mask + decoy | 4 shares × (3 × 2 + 2 × 3) × 10 |

Every profile is queried with `--all-estimators`, so the estimator ablations (EM against
gather, winner-take-all and uniqueness-first; zero inflation; *p_in* weighting) are the
`abundance` rows of the same `summary.tsv`. Read the results from `summary.tsv` and
`scores.tsv`, grouped by `label`, `index` and `arm`, at `min_hits` 1 and the shipped pair
(`kmers_unique`, `abundance_zi`), Pfam first:

- *every arm:* `purity`, `completeness`, `completeness_low25`, `l1`, `spearman_tp`,
  `ci_cover`; cost per arm from `bench.py cost <outdir>/trace.tsv` (PROFILE rows are named
  `index~arm`);
- *index:* also `split_l1` (within-component share error, weighted by true depth) and
  `group_size_mean` (ambiguity groups), which the function × taxon table depends on;
- *host:* also `n_pred - tp` (false detections) and `host_like_detected` per share and arm,
  and the mask's cost in `masks/<index>.json` (masked hashes and postings).

The host samples have ids `<seed>h<percent>` (e.g. `seed3h90`); their truth is recomputed on
the mixed reads, so it counts only the microbial reads kept. Masks and decoys exist only for
the indexes built in the run (not `fmh_compat` or `--mgnify_indexes`). Locally,
`-profile test,docker -c <arms config>` runs every arm on the mini fixture (a 60 kb random
"host" and the fixture proteins as decoy) in a few minutes; hostile needs its ~4 GB index, so
leave it out there.

Not here: the divergence ladder and containment-AAI calibration (truth needs held-out
genomes with relatives at known identity; the simulation benchmark has them), and the
other tools on Pfam.

## Other tools

Each tool profiles the same 10 metagenomes and is scored against the same truth; its rows
in `summary.tsv` have `count` = `evidence` and `abundance` = `abundance`, on the tool's own
scale (Spearman and L1 only need proportionality to depth).

| Tool | Reference | Detection (`evidence`) | Abundance |
| --- | --- | --- | --- |
| `diamond` (DIAMOND 2.2.8 blastx, `-k 1`, default sensitivity) | The KEGG proteins our indexes are built from | Read pairs whose best hit is a gene of the KO | Σ aligned / subject length over those reads (as truth depth) |
| `fmh_funprofiler` (fmh-funprofiler 1.1.1, `funcprofiler`) | Its KO sketches (`--sketches`, scaled 1000, k = 11) | Hashes shared with the KO (`intersect_bp` / scaled), `--threshold_bp` = scaled | Its output (normalised `f_match_query`) |
| `kmermaid` (kMermaid at `edcb4ed`) | Retrained with KOs as clusters, `--kmermaid_max_members` random proteins per KO | Read pairs assigned to the KO (score ≥ 3) | Reads / mean member length |
| `humann` (HUMAnN 3.9, MetaPhlAn 4.1.1) | Its own ChocoPhlAn + UniRef90 databases; UniRef90 families regrouped to KOs with its mapping | 1 for every KO it reports (no read counts; `min_hits` 2 repeats 1) | Unstratified RPK |
| `humann4` (HUMAnN 4.0.0.alpha.2 at `e07b3a3`, MetaPhlAn 4.1.1) | Its v4 alpha ChocoPhlAn + EC-filtered UniRef90 (the only protein database it distributes); regrouped as above | As `humann` | Unstratified adjusted CPM |

Why each one:

- **DIAMOND** is the alignment baseline on the same reference: the accuracy a translated
  search reaches on exactly our KO units, and what it costs.
- **fmh-funprofiler** runs the released tool (sourmash sketch translate + prefetch) for its
  cost; its calls should equal `fmh_compat`'s `kmers_hit` rows (they do on the test profile).
- **kMermaid** ships a model of RefSeq protein clusters without KO labels, so it is retrained
  on the KO units with its own training code (`kmermaid_kfp.py`), as its README allows. The
  per-KO cap keeps the model near the size of the shipped one (~55 proteins per cluster);
  uncapped KEGG KOs would not fit in memory as its Python dicts.
- **HUMAnN** is the de facto standard, on its own databases. The sample genomes' proteins
  are in UniRef90, so it has no divergence handicap, but its KO calls go through the
  UniRef90 → KO mapping rather than KEGG's gene → KO table the truth uses, so part of any
  gap is annotation, not detection.
- **HUMAnN 4** is the current development release (alpha, GitHub only), run as distributed:
  its translated search covers only UniRef90 families with an EC number, so KOs without
  one are found only through the nucleotide search against the pangenomes of species
  MetaPhlAn detects. MetaPhlAn runs as a separate command before HUMAnN, whose output goes
  in with `--taxonomic-profile`. Run from inside HUMAnN 4 alpha.2, it would fail: HUMAnN
  looks for its database tag (`vOct22_CHOCOPhlAnSGB_202403`) in `metaphlan --version`, and
  MetaPhlAn 4.1.1 does not print one.

Threads differ (DIAMOND 8, HUMAnN 16, the others 1); compare CPU hours in `cost.tsv`, not
wall time.

## Results (10 metagenomes, InSilicoSeq novaseq)

Means over seeds 1..10; sd of purity and completeness <= 0.015. "Low 25%" is completeness
on the least-covered quarter of true KOs. `kmers_hit` counts every hit k-mer for every KO
holding it; `kmers_unique` counts only the k-mers gather assigns to the KO.

| Index | Count | min_hits | Purity | Completeness | Low 25% |
| --- | --- | --- | --- | --- | --- |
| fmh_compat | kmers_hit | 1 | 0.975 | 0.688 | 0.295 |
| fmh_compat | kmers_unique | 1 | 0.991 | 0.672 | 0.268 |
| kfp_s1000 | kmers_hit | 1 | 0.985 | 0.677 | 0.285 |
| kfp_s1000 | kmers_unique | 1 | 0.996 | 0.663 | 0.263 |
| kfp_s1000_floor8 | kmers_hit | 1 | 0.976 | 0.719 | 0.404 |
| kfp_s100 | kmers_hit | 1 | 0.951 | 0.960 | 0.853 |
| kfp_s100 | kmers_hit | 2 | 0.976 | 0.911 | 0.694 |
| kfp_s100 | kmers_hit | 3 | 0.985 | 0.858 | 0.533 |
| **kfp_s100** | **kmers_unique** | **1** | **0.985** | **0.955** | **0.837** |
| kfp_s100 | kmers_unique | 2 | 0.994 | 0.903 | 0.665 |

`kfp_s100` scored on `kmers_unique` at `min_hits` 1 is the recommended setting: the purity
of `kmers_hit` at 3 with 10 points more completeness, and 27 points more than fmh_compat.
Gather removes 70% of the false positives (348 to 105 per sample; none left in all 10
samples, against 29 before) for 31 true KOs per sample, 90% of them in the lowest-coverage
quarter: single-k-mer hits taken by a more abundant relative (PTS and ABC transporter
paralogs, for example). The false positives it removes are mostly modular PKS/NRPS KOs
(pks2, pks8, tyrocidine and rapamycin synthases), whose shared domains carry identical
k-mers. Collagen VII (`K16628`) is the most persistent one left, in 5 of 10 samples.

At scaled 1000, `min_hits` above 1 costs most low-coverage KOs (completeness 0.49 at 2).
The floor helps at scaled 1000 but not at scaled 100, where nearly every KO already
samples more than 8 k-mers (`n_min 8` at `t_base 0.01` scored within 0.005 of `n_min 0`),
so that configuration was dropped. Gather adds no measurable query time (10-15 s per
sample either way).

### Abundance (same run design, `kmers_unique` >= 1, true positives)

| Index | Abundance | Spearman | L1 |
| --- | --- | --- | --- |
| fmh_compat | coverage | 0.08 | 1.41 |
| kfp_s100 | coverage | 0.32 | 1.35 |
| kfp_s100 | coverage_em | 0.35 | 1.25 |
| kfp_s100 | coverage_zi | 0.58 | 0.94 |
| kfp_s100_d10 | coverage_zi | 0.63 | 0.93 |
| kfp_s1000 | abundance_zi | 0.81 | 0.57 |
| kfp_s100 | **abundance_zi** | **0.96** | **0.22** |
| kfp_s100_d10 | abundance_zi | 0.99 | 0.14 |

Truth depth sums over each KO's gene copies in the sample, while `coverage_zi` is depth per
copy; `abundance_zi` scales it by the copies present (present k-mers / the k-mers an average
member holds) and matches truth closely. Recommended: `kfp_s100`, `kmers_unique` >= 1,
`abundance_zi`. For `fmh_compat` it reduces to plain EM, since imported sketches carry no
members. The dense tier costs 15x the table size, ~6x the query time and 36.5 GB to build.
