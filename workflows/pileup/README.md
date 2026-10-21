# pileup

Read recruitment and crude pileups downstream of a kfp profile (plan: Read recruitment and
variant pileups, phase 13; the nucleotide follow-up of phase 14, step 4). kfp profiles the
sample and shrinks the reference; an existing aligner does the alignment; `kfp-pileup` (a
separate tool, plan: Tools and scope) turns the alignments into read posteriors and pileups.

Per sample:

1. `QUERY`: `kmer-functional-profiler query` against the index (built from
   `--targets` by `INDEX` if no `--index` is given); `--aai true` adds the `aai` (`ani`)
   columns, and `--min_identity` (theta, implies `--aai`) makes relatives below theta not
   present. Neither is needed for targets.
2. `REFERENCE`: `kfp-pileup reference`: the top units by `present_prob` at a false
   discovery rate of `--max_fdr` (`present_q` <= 0.05 by default) as targets (at most
   `--max_targets`, the most hit first; 0, the default, is no cap), the other hit units of their components as neighbours (at most
   `--neighbours` per target in a component, the most-hit first); up to `--members` sequences per unit, from kfp's `fetch` (the
   index's `targets.parquet`, or `--store` members tables, e.g. the partitioned build's MERGE
   output).
3. The aligner the index's molecule calls for:
   - amino-acid index: `ROUTE` (`kfp-pileup route`; `--route false` to skip) packs the
     reference's components into blocks of about `--route_block_letters` residues and sends
     each pair to every block it shares at least `--route_min_kmers` distinct k-mers with
     (every k-mer of the reference at `--route_k`, 0 for the index's k; the reads in all six
     frames); other pairs are dropped. `--route_members` (3, as `--members`) sets the sequences per
     unit in routing's k-mer table apart from DIAMOND's reference (`REFERENCE` then writes
     `route.fa`). `DIAMOND` then runs per block
     on its pairs (`makedb`, then `blastx` per mate, e-values on the whole reference's size
     with `--dbsize`, every hit within `--diamond_top` % of the read's best,
     frameshift-aware with `-F --frameshift` (0: off), tabular with BTOP);
   - nucleotide index: `BWAMEM2` (`index`, then `mem -a`, SAM).
4. `CALL`: `kfp-pileup call`: read posteriors (leave one column out), neighbour-best reads dropped, pileups,
   variants and per-unit identity; read assignment (steps 8-9); with DIAMOND, the reads
   brought back to nucleotides along their alignments, a nucleotide pileup and a consensus
   coding sequence per unit (step 6b).
5. The hand-off (steps 6b-6c): `MAP_CONSENSUS` re-maps the recruited pairs to the consensus
   with bwa-mem2 (amino-acid index; a nucleotide index uses `BWAMEM2`'s alignments to
   `ref.fa`); `HANDOFF` (`kfp-pileup handoff`) keeps each mate's best alignment to its pair's
   unit where step 8 gives that unit at least `--min_assign`, MAPQ from the posterior; `BAM`
   sorts and indexes it.
6. Existing tools on the BAM (step 6c): `FREEBAYES` (`--pooled-continuous`), `LOFREQ`
   (`--lofreq true`), `CONSEQUENCES` (`bcftools csq`, amino-acid index: each consensus is
   one CDS), `FLORIA` (strain haplotypes; `--floria false` to skip), `CLIQUESNV` (haplotypes from linked
   SNVs, per sequence with at least `--cliquesnv_min_reads` hand-off records; `--cliquesnv
   false` to skip) and `INSTRAIN`
   (`inStrain profile`, per-gene microdiversity; `--instrain false` to skip).

Outputs under `<outdir>/<sample>/`:

| File | Contents |
| --- | --- |
| `profile.tsv` | kfp's profile |
| `ref/ref.fa`, `ref/refs.tsv` | the per-sample reference and its roles |
| `route.json` | routing: pairs, pairs routed and dropped, blocks, reference residues and k-mers |
| `routes.parquet` | routing: every (pair, block) hit, routed or not: `read_name`, `pair` (input order), `block`, `kmers` (distinct k-mers shared); another `--route_min_kmers` is a filter on `kmers` |
| `pileup/assignments.parquet` | every alignment: read, mate, reference sequence, unit, role, score (bits), posterior, `kept` (best hit a target); steps 8-9's `em_posterior`, `em_novel` (the mate's novel-relative share), `em_flag` (best source differs with the prior) |
| `pileup/pileup.parquet` | summed weights per (sequence, column, residue), with the part from sibling residues |
| `pileup/variants.tsv` | columns with depth >= `--min_depth` and an alternative residue at >= `--min_freq`, with `sibling_share` |
| `pileup/units.tsv` | per target unit: reads, alignment identity (all sequences, representative), the representative's consensus identity, and the profile's `aai` (`ani`) columns, for the `aai` check; with steps 8-9, `em_pairs`, `em_abundance`, `em_divergence` |
| `pileup/sources.tsv` | step 8: per source (unit, or a component's novel relative): pairs (likelihood-only and with the prior), divergence, effective length, abundance |
| `pileup/haplotypes.tsv` | step 9: per target sequence, haplotypes (0: the original sequence), share, reads, divergence from the original and of their reads from the unit's nearest sequence, `novel` |
| `pileup/pileup_nt.parquet` | step 6b (DIAMOND): summed weights per (sequence, nucleotide column, allele): a base, `-` (deleted), or `+` and inserted bases (before the column); a diagnostic |
| `pileup/consensus.fna`, `.gff`, `.genes.fna`, `.tsv` | step 6b: per unit, the consensus coding sequence of its best-supported sequence (substitutions only, `N` where uncovered; columns 3c..3c+2 are reference column c), as one CDS for `bcftools csq`, as prodigal genes for inStrain, and its coverage |
| `pileup/recruited_{1,2}.fq.gz` | step 6b: both mates of every pair with a target alignment |
| `variants/handoff.bam` | the hand-off: assigned reads against `consensus.fna` (or `ref.fa`), MAPQ from step 8's posterior |
| `variants/freebayes.vcf.gz`, `lofreq.vcf.gz` | step 6c: variant calls (SNVs, indels with frameshifts, MNPs) and their frequencies |
| `variants/consequences.vcf.gz`, `.tsv` | step 6c: `bcftools csq` per variant (synonymous, missense, stop gained, frameshift, in-frame indel), amino-acid positions in the reference protein's columns |
| `variants/floria/`, `variants/instrain/` | step 6c: Floria's haplotype sets per sequence; inStrain's per-gene and per-SNV tables |
| `variants/cliquesnv/` | step 6c: CliqueSNV per sequence (`<seq_id>/`: haplotype FASTA and JSON; `<seq_id>.log`), and `haplotypes.tsv` (sequence, haplotype, frequency) |

`--query_stats true` publishes `<sample>/stats.json` (the query's stage times and
`hits_per_pair`, step 0's histogram). `kfp-pileup aai-check <sample>/pileup/units.tsv` is
step 6's check on a real sample: `aai` against the reads' alignment identity by identity,
depth, length and cluster size, flagging step 6's rule.

`trace.tsv` records each task's time: the aligner's against `QUERY`'s is step 3's rule (add
routing if the aligner takes more than ~5x; routing is on by default).

```bash
nextflow run workflows/pileup -profile docker --samples samples.csv --index idx/
nextflow run workflows/pileup -profile test,docker              # fixture proteins, DIAMOND
nextflow run workflows/pileup -profile test_nucleotide,docker   # fixture CDS, bwa-mem2
```

The steps are the subworkflow `PILEUP` (reads, index, store, units to withhold from the
reference), which `workflows/read-xref` and `spike-benchmark` reuse. Options reach the
processes as inputs, not as `params` in their scripts, because Nextflow does not hash
`params`: `-resume` reruns a task when its options change.

`samples.csv` has columns `sample,r1,r2`; relative paths are read from the CSV's directory.
DIAMOND, bwa-mem2, samtools, freebayes, LoFreq, bcftools, Floria, CliqueSNV and inStrain run in biocontainers (`-profile docker` locally, `singularity` on HPC);
kfp and kfp-pileup run on the host, in the venv `workflows/setup.sh` makes.
