# hmm-annotate

Annotates MGnify proteins with any HMM collection (Pfam, KOfam, ...): kfpfam's k-mer prefilter picks
candidate models, and HMMER (pyhmmer) scores only those. Or nail's seeds are the candidates, with
`--engine nail`, or HMMER scores every model, with `--engine hmmsearch`. It gives the domain rows with alignments
that `lookup-index` builds from, and a label table (protein, model) for any collection.

```text
PREP           HMM file -> kfpfam prep, tiled index and gapped bars, pressed models, thresholds.tsv
PLAN           ids of the targets -> ranges.tsv, protein_id ranges of ~--chunk_proteins   (pfam-reannotate's)
FASTA          per range: the sequences                                                   (pfam-reannotate's)
SEARCH         per range: filter join | annotate.py search -> domains/ and labels/
               (--engine hmmsearch: annotate.py search --all-pairs, no index)
NAIL_SEED      per range, with --engine nail: its seeds, for SEARCH; NAIL_MEMBERS, the same of the members
MEMBER_PAIRS   with --members: the models each representative nearly passes -> its members' pairs, by member id
MEMBER_FASTA   per --member_width ids: those members' sequences and pairs
MEMBER_SEARCH  per --member_width ids: HMMER on the pairs -> domains/ and labels/
               (--member_filter: on those and on the prefilter's, of every member)
CHECK          per range: 1 in --check_sample proteins against every model, no prefilter  (--check_sample N)
COMPARE        per range: the passing domains against the release's own Pfam table     (--compare_release true)
SUMMARY        -> summary.tsv, one total row
CHECK_SUMMARY  -> check.tsv and compare.tsv, recall against each, worst model first
```

## Setup

Needs Java 17+, Nextflow (not 26.04, whose `-resume` reuses tasks run with other parameters), [uv](https://docs.astral.sh/uv/) and kfpfam's image.

```bash
bash workflows/setup.sh   # creates .venv in the repo root
```

kfpfam's image is published nowhere. Build it where Docker is, from kfpfam's `prep-any-hmm` branch
(`63a4cf6`) or later, and convert it on the HPC. Older ones lack `top`'s windows, need `ACC` and `GA` lines, which KOfam's profiles
lack, panic on a protein over 65,535 residues, and are very slow at these settings.

```bash
docker build -f docker/filter.Dockerfile --platform linux/amd64 -t kfpfam-filter <kfpfam checkout>
```

```bash
docker save kfpfam-filter:latest -o kfpfam-filter.tar
```

```bash
singularity build kfpfam-filter.sif docker-archive://kfpfam-filter.tar
```

## Run

Test profile (six mini HMMs, 84 proteins; needs a kfpfam checkout with the filter built, and a Python
with numpy, pyhmmer and pyarrow):

```bash
nextflow run workflows/hmm-annotate -profile test --kfpfam_filter <kfpfam>/filter/target/release/filter --kfpfam_python <python> --kfpfam_prep_py <kfpfam>/measure/prep.py
```

On HPC, `hpc/kfp-ablations/hmm-annotate/run.sh` runs 1 in 1,000 proteins with the check by default.

| Option | Default | Meaning |
| --- | --- | --- |
| `--hmm` | | The collection: one HMM file, plain or gzipped |
| `--engine` | `kfpfam` | `kfpfam`: its prefilter, then HMMER on the candidates; `nail`: its seeds (an MMseqs2 profile search), then HMMER on those; `hmmsearch`: HMMER on every model |
| `--cutoff` | `gathering` | The HMM file's own thresholds: `gathering`, `trusted` or `noise` |
| `--thresholds` | | Or a TSV with a header: model, threshold, `full` or `domain`. KOfam's `ko_list` as it is |
| `--targets` | `all` | `all`: every protein; `reps`: MGnify90 representatives |
| `--members` | `false` | With `--targets reps`: their members too, so every protein. A member (90% identity) is searched with HMMER against only the models its representative passes to within `--member_slack` (8) bits |
| `--member_filter` | | With `--members`: `kfpfam` or `nail` on every member as well, and HMMER on the models of either |
| `--min_cluster` | `1` | With `--member_filter`: `--engine` has the representatives of clusters this size or more only. The prefilter has the whole of the others, whose members get no candidates from their representative |
| `--sample` | `1000` | 1 in N targets by `protein_id`; with `--members`, 1 in N clusters |
| `--check_sample` | `0` | 1 in N of those also against every model, without the prefilter; `0`: none |
| `--compare_release` | `false` | Compare with the release's `mgy_proteins_pfam.parquet`, made with Pfam 38.2 |
| `--domains` | `true` | Write the domain rows; labels are always written |
| `--weak_evalue` | `0` | Also domains under the thresholds, at this E-value (Z = 10⁶); `0`: none |
| `--chunk_proteins` | `1000000` | Proteins per `SEARCH` task, the same for every engine so that they share FASTAs |
| `--nail_seed_pvalue` | `0.01` | nail's `-S`, as the B11 arm's; nail's own default is 1e-4 |
| `--nail_proteins` | `1250000` | Proteins nail has at a time, over a chunk's most; a member task's (10 M) are split |
| `--nail_sif` | nail 0.7.1's biocontainer | nail's image, with MMseqs2 18. 0.8.0 writes its seeds differently |
| `--kfpfam_sif` | | kfpfam's image |
| `--kfpfam_mode` | `top` | `top`: as near every domain as it goes; `fast`: gapped bars per model from its threshold, 99.8% of domains. Settings: `kfpfam_modes` in `nextflow.config` |
| `--kfpfam_band`, `_diagonals` | 15, `any` | A band around every diagonal at the bar, until one reaches the gapped bar |
| `--kfpfam_block` | `80000` | Proteins the filter holds at a time |

A model in the HMM file with no threshold is an error. A `-` in the table is not: that model never
passes, and shows only with `--weak_evalue`.

## Output

- **`domains/chunkN.domains.parquet`**: one row per passing domain, in the HMM file's model order:
  `protein_id`, `model` (accession without its version, or the name), `hmm_from`, `hmm_to`,
  `ali_from`, `ali_to`, `env_from`, `env_to`, `cigar`, `score`, `seq_score`, `i_evalue`, `passes`,
  and `left`, `right`, `residues`: the aligned residues with up to 10 either side, and how many
  those are. `cigar` is `M`, `I` and `D` against the model.
- **`labels/chunkN.labels.parquet`**: `protein_id`, `model`, for proteins passing the model's
  thresholds.
- **`ranges.tsv`**, and **`ranges.members.tsv`** with `--members` (chunks `mN`): the `protein_id` range of each chunk.
  `lookup-index` reads both.
- **`thresholds.tsv`**: the thresholds used (`inf`: never passes; `-inf`: not used).
- **`hmm.sha256`**: of the HMM file, unzipped. `lookup-index` checks its models against it.
- **`logs/chunkN.log`**: the CPU, and kfpfam's log with its SIMD kernel and CPU time. `bench.py` compares
  runs of the same proteins from it, `summary.tsv` and `trace.tsv`; `hpc/kfp-ablations/hmm-annotate/benchmark.sh`
  makes the runs.
- **`summary.tsv`**: proteins, residues, proteins over 100,000 residues (HMMER's limit, not
  searched), candidate pairs, domains and proteins passing, and pyhmmer's core-seconds.
- **`check.tsv`**: per model, passing (protein, model) pairs and domains by full `hmmsearch`, how many
  `SEARCH` found, and how many it has that `hmmsearch` doesn't (`extra`, none expected). `check.members.tsv`: of the members.

- **`compare.tsv`**: the same columns with the release's own Pfam table as the truth; a domain is the same if
  its envelope is. Use it with Pfam 38.2 only, and `--cutoff gathering`, as the table was made.

`work/` holds the FASTA of every chunk until deleted: about 1.1 TB for every protein of the release.

## Measured

On the laptop (M4 Pro), 19,655 proteins of the 40-family slice (4.9 M residues) against its 40 models:
13,293 of the 13,318 domains at GA (0.9981), worst family 0.975 (PF14306), no extra domain, 64 bytes a
domain row. This is whole proteins of families that are there; the recall over all of Pfam, the recall on
KOfam (whose models are longer) and the cost per residue are for the 1-in-1,000 run to give. That was
with kfpfam's Mode A, before `--kfpfam_mode`.

Every fourth model of Pfam 38.2 (7,534), M4 Pro, 4 threads, CPU-s, native ARM, 2026-11-25. Truth is
`--engine hmmsearch`, which gives `hmmsearch`'s domains exactly; no mode found one it doesn't.

| 80,000 proteins (13.6 M residues) | filter | HMMER | against `hmmsearch` | domains missed, of 13,880 |
| --- | ---: | ---: | ---: | ---: |
| `hmmsearch` | | 404 | 1 | 0 |
| top as for reads (band 16, every diagonal) | 344 | 105 | 0.9× | 0 |
| **`top`** (band 15, `any`) | 150 | 105 | 1.6× | 0 |
| HMMER on the passing pairs alone: the floor | | 7.5 | 54× | 0 |

`fast` on 108,471 other proteins (18.0 M residues): 26 + 30 CPU-s against 530, **9.4×**, 38 of 17,667
domains missed (0.2%). Its bars were chosen on the first 80,000, and were 0.7 bits higher then, with
a query length in them.

`fast` assumes no query length. The same settings, (protein, model) pairs passing:

| | filter | HMMER | against `hmmsearch` | pairs missed |
| --- | ---: | ---: | ---: | ---: |
| 80,000 whole proteins | 21 | 26 | 8.6× of 403 | 41 of 12,867 (0.32%) |
| the same in 290,248 pieces of 50 residues | 9 | 4.5 | 32× of 438 | 13 of 16,070 (0.08%) |

Why this is not the 30× of reads:
- the same residues cut into 50-residue pieces take top's filter 47 to 54 CPU-s, against 344 (band 16)
  or 150 (band 15): a diagonal at the bar is scored over the whole protein, and each hit's window
  is full;
- top's bars are for reads, so it sends HMMER 2.2 M pairs, 170 times the passing ones;
- `hmmsearch` is at its best on this CPU: on those pieces, all of it, top is 6.8× faster. Codon's x86
  (AVX-512 for kfpfam, SSE for HMMER) should favour kfpfam. Not measured on proteins.

`--members`, 49,998 proteins in 8,953 clusters of the 40-family slice, the same models, through the
workflow; domains missed against `hmmsearch` on every protein (357 CPU-s), none extra:

| Representatives by | representatives missed, of 5,115 | members missed, of 16,143 |
| --- | ---: | ---: |
| `hmmsearch` | 0 | 12 (0.07%) |
| `kfpfam`, `top` | 1 | 12 (0.07%) |
| `kfpfam`, `fast` | 6 (0.12%) | 228 (1.4%), 196 of them in one cluster of 574 |

Members cost about 9 CPU-s, 35× less than `hmmsearch` on them, whatever the representatives' engine. A
representative missed loses its cluster's domain, so `fast`'s loss is the same rate in expectation but
comes in clusters; looser bars (`slack` 14 or 18) did not change it. The slice is of chosen families with
large clusters: 1 in 1,000 clusters and `check.members.tsv` give it for the release. Estimated for the
release at 38.2 and laptop speed: 63,600 core-h by `hmmsearch`; with `--members` about 16,000 with
`hmmsearch` on the representatives (23% of the residues) and 2,500 with `fast`. 39.0 has 1.8× the models.

With `--member_filter` and `fast`, the same 49,998 proteins; CPU-s of kfpfam and HMMER together:

| | representatives missed, of 5,115 | members missed, of 16,143 | CPU-s | against `hmmsearch` |
| --- | ---: | ---: | ---: | ---: |
| `hmmsearch` on everything | 0 | 0 | 357 | 1 |
| `fast` on everything | 6 | 249 | 41 | 8.7× |
| `fast` on everything, and representatives' hits to members | 6 | 215 | 42 | 8.5× |
| `hmmsearch` on representatives, no `--member_filter` | 0 | 12 | 85 | 4.2× |
| **`hmmsearch` on representatives**, `--member_filter` | 0 | 0 | 108 | 3.3× |
| the same, `--min_cluster 2` | 1 | 0 | 75 | 4.8× |
| the same, `--min_cluster 10` | 3 | 1 | 51 | 7.0× |

211 of the domains `fast` misses are of PF14306, 196 of them in one cluster: it misses that family in members as it does in
their representative, so its own hits do not bring them back, and `hmmsearch` on the representative does. Here
representatives are 21%, 11% and 2.6% of the residues at `--min_cluster` 1, 2 and 10. In the release (three row
groups of `mgy_cluster_seqs`, 2026-11-25) 62-74% of clusters are single proteins, and those of 10 or more
are 4-9% of clusters with 53-71% of the proteins. Estimated for the release at 38.2 and laptop speed, against 63,600 core-h: about 21,000 at
`--min_cluster` 1, 12,000 at 2 and 9,500 at 10.

nail 0.7.1 and MMseqs2 18, the same 49,998 proteins and 7,534 models, 4 threads, in one task (`bench.py`):

| | domains missed, of 21,258 | CPU-s, seeds or filter + HMMER | against `hmmsearch` |
| --- | ---: | ---: | ---: |
| `hmmsearch` | 0 | 357 | 1 |
| `kfpfam`, `top` | 1 (0.005%) | 139 + 103 | 1.5× |
| `kfpfam`, `fast` | 255 (1.2%) | 18 + 22 | 8.9× |
| `nail`, `-S 0.01` | 928 (4.4%) | 164 + 29 | 1.8× |
| `nail`, `-S 1e-4` | 963 (4.5%) | 107 + 14 | 3.0× |
| `nail` on representatives, HMMER on members | 915 (4.3%) | 117 + 15 | 2.7× |
| `hmmsearch` on representatives, `--member_filter nail` | 1 | 135 + 100 | 1.5× |
| the same, `--min_cluster 10` | 285 (1.3%) | 161 + 39 | 1.8× |

Every domain nail misses is of a pair it did not seed; the families missed most are hEGF (PF12661, 251 of 1,907), PF14306
(216 of 3,331), PF09393 and EGF. Every pair it does seed gives `hmmsearch`'s domains.

Why 1 M proteins a task (measured at 1,000, 20,000 and 80,000 proteins and half the models):
- nail spends about 65-100 CPU-s on 7,534 models whatever the proteins, so some 460-700 on 39.0's 53,275,
  and then about 2 ms a protein (14 at 39.0). So start-up is 3-5% of a task's CPU at 1 M proteins, 25-40%
  at 100,000. The table above is of one task of 50,000, where start-up is nearly half: nail should
  be faster against `hmmsearch` at 1 M (estimate: 2.5-3×).
- its temporary files are 8 KB a protein, so about 55 GB a task at 39.0 (estimate), in the task's directory and
  removed as it ends; its memory is 0.7 GB and 3 MB per 1,000 proteins, so some 10 GB with the seeds
  (estimate; it asks for 32). At 2.7 KB of seeds a protein, 1 M is 3 GB, less gzipped; all of the release's would
  be about 2 TB in `work/`.

### Why whole proteins are slower than reads, and what helps (2026-11-25)

24,000 proteins (6.2 M residues, median 191) of the slice, and the same residues cut into pieces; 7,534 models,
M4 Pro, 4 threads, CPU-s. `hmmsearch` takes 179-196 at every length.

| `fast` | join | gapped | survivors | HMMER | all | against `hmmsearch` |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 50-residue pieces | 3.3 | 0.7 | 0.3 | 3.9 | 8.0 | 25× |
| 200 | 4.2 | 3.2 | 0.7 | 9.3 | 17.0 | 11× |
| whole | 4.8 | 5.0 | 1.0 | 12.2 | 22.0 | 8.1× |

- **Seeds are not it.** Hits per residue hardly change with length (2.40 to 2.75 G), as the index is
  per residue. A 4-, 5- or 6-residue seed costs the same at any length, so indexes by query length
  would not help.
- **The gapped stage scores the whole stretch of the diagonal that meets the model**, up to the model's
  length: about 50 rows on a piece and 140 on a protein, and more pairs reach the bars by chance.
- **HMMER is most of it.** Of its 11.8 CPU-s on the whole proteins, its filters and Forward are 2.8 and
  domain definition 9.0. Domain definition on the 8,654 true pairs takes 3.0, the floor (60× `hmmsearch`).
  False candidates that pass HMMER's own filters cost 6.9: 7,000 under the threshold
  by Forward (3.0 CPU-s), and 4,100 above it that only the composition correction rejects, mostly short
  repeat models on long proteins (4.0).

Tried:

| | Change | Cost | Recall |
| --- | --- | --- | --- |
| **`--forward-slack 8`** (built, the default) | no domain definition for a candidate 8 bits under the sequence threshold by HMMER's Forward filter | HMMER −22% (`fast`), −20% (`top`) | no domain lost; 4 bits lost 2 |
| θ 4 for `fast`, not 5 | fewer seeds | filter 9.8 → 7.6 | pairs missed 69 → 89 of 8,723 |
| a spaced 5-mer, 1101011, θ 4 | 0.54× the hits of the 4-mer at θ 4 | join the same: its lists fill 32% of the SIMD lanes (74% in 4× larger blocks, then 0.85× the 4-mer's) | the same; index 12× the size |
| 6-mers | 0.34× the hits at θ 5 | | index 137× the size |
| HMMER's MSV on every pair first | | 92 CPU-s alone, more than all of `fast` | |
| a composition correction on the filter's best segment | | removes 0.8 of HMMER's 12 at no loss | |

Not tried: an X-drop gapped stage that stops when a chance extension dies (estimated up to 2 CPU-s of the 22).
So at most about 11× on whole proteins at `fast`'s recall on this CPU, against 25× on 50-residue pieces: what is
left is HMMER on true hits and on repeat-rich near misses, which `hmmsearch` pays as well. Members
(`MEMBERS=fast`) are where the large factor is.
