# read-xref

The read cross-reference pipeline (plan: Genome profile and functional update, Read
cross-reference pipeline; phase 12): a function (MGnify90 unit) × genome table built from
reads directly, each read classified twice and independently, then joined on read id. A
comparison arm and development aid for `kfp-genomes`, not truth: bowtie2 needs references
within ~95% ANI of the sample's strains, and DIAMOND's sensitivity floor carries over.

Per sample:

1. *Read → function:* `PILEUP` from `workflows/pileup` (kfp query, the recruited reference,
   DIAMOND or bwa-mem2, read posteriors with decoys): `assignments.parquet`.
2. *Read → genome:* `SYLPH` (`sylph profile`, against a sketch of `--genomes` made by
   `SKETCH`, or `--sylph_db`) picks the genomes present; `BOWTIE2` maps every pair to their
   concatenated contigs (`-k --bowtie2_k`) and keeps per-mate hits as TSV.
3. `XREF` (`python -m kfp_genomes.xref`): mates below `--min_identity` (or `--min_mapq`)
   dropped; per pair, the genomes where the most mates align with the best summed `AS`
   share it evenly; per pair, the target units by renormalised posteriors. A pair counts
   once and mates that disagree are split. Pairs with a function and no genome go to
   `unclassified`.

Outputs under `<outdir>/<sample>/` (beside `workflows/pileup`'s):

| File | Contents |
| --- | --- |
| `sylph.tsv` | sylph's profile |
| `contigs.tsv` | contig (`genome\|n`), genome, length |
| `xref/function_genome.tsv` | `unit`, `rank` (`genome`; `species`, `genus` when `--genomes` has `taxonomy`), `taxon`, `reads` (pairs), `share` of the unit's pairs, `coverage`, `copies` (coverage / the taxon's depth) |
| `xref/genomes.tsv` | per `rank` and `taxon`: `reads` (all pairs bowtie2 gave it), `function_reads`, `length`, `depth` |

```bash
nextflow run workflows/read-xref -profile docker --samples samples.csv --index idx/ --genomes genomes.tsv
nextflow run workflows/read-xref -profile test,docker   # fixture reads, two fixture genomes, no sylph
```

`genomes.tsv` has columns `genome`, `path` (nucleotide FASTA, gzipped or not, relative to
the TSV allowed, unique file names) and optional `taxonomy` (GTDB-style
`d__;...;g__;s__`). Every option of `workflows/pileup` applies (`--index` or `--targets`,
`--store`, `--min_present`, ...). The test profile maps to every genome (`--sylph false`):
sylph detects nothing from 40 read pairs.

Not done: the full-reference function arm (every unit the panel's genomes carry; needs
`kfp-genomes panel`), Pfam and KO roll-ups (the scoring code's, from the index's unit
tables), and scoring against CAMI II's read-level truth.
