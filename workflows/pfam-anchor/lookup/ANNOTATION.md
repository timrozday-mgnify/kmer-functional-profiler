# An annotated MGnify proteins index: KOfam, dbCAN and others. What size?

2026-11-24. Scripts in `annot/` (`genomes.py`, `collapse.py`, `sizes.py`), tables in `data/pfam-anchor/lookup/report/annot_*`. *Measured* on
21 gut genomes (mimicc 20HM+: 71,278 proteins with KOfam at its own thresholds, dbCAN by 2 of 3 tools, Pfam 37.4 at GA) and on the 40-family
MGnify90 slice. Release sizes are *estimated*.

## Answer

- **Labels are small.** On the residue index (80 GB at 5 bits, 123 GB at 1 byte), one scheme such as KO or CAZy adds **at most 2 GB** (2 bytes a row), and a
  sparse one like CAZy under 0.5 GB.
- **An index only for annotation can be much smaller, as it need not tell clusters apart.** Keep one row per group of rows that are 70% identical by column, and it is
  **about 26 GB (5 bits) or 40 GB (1 byte) for every Pfam-domain row**, about 15 GB for rows with a KO, and about 1.3 GB for rows with a CAZy family. At 60%: 16, 9 and 0.8 GB.
- **Collapsing to 70% costs almost nothing for annotation (estimated):** expected KO recall .987 and precision .985, against .987 and .988 uncollapsed.
- **Why use MGnify and not the KO models on reads.** A KO's threshold is for the whole protein: the median is 321 bits and 84% are over 110, more than a 50-residue
  read can usually score. dbCAN's rule needs a third of a model covered. Annotating whole proteins and passing the label to reads avoids both.
- **Reach, the main limit.** 98.8% of KO proteins and 97.4% of CAZy proteins have a Pfam GA domain, and 71% and 65% of their residues lie inside one. **But only 41% of 50-residue windows of KO proteins align at GA**
  (35% for CAZy; measured on held-out genomes), as a fragment of a domain often scores under GA. An earlier version of this page said about 70% of such reads would be reached: that was wrong. Alignments under GA recover much of it (below).

## The self-contained bundle (the user's requirement, 2026-11-24)

Sizes were accepted. The bundle must hold everything kfpfam and this tool need to turn metagenome reads into a profile of annotations, with no protein database or other file beside it.

| Part | Size | Status |
|---|---|---|
| kfpfam: index (r999), prep, Pfam HMMs | 3.4 GB on 37.4 (1.14 + 0.34 + 1.88); 6.4 GB on 39.0 | measured on Codon, B11 report |
| Residues of each kept row, inserts included, and 10 either side; row header and CIGAR (16 bytes, measured) | **32 GB at 5 bits** (rows collapsed at 70%); 20 (60%); 102 (every cluster) | estimated; the format is 5 bits |
| Labels, 2 bytes a row and scheme | 0.6 GB a scheme at 70% | estimated |
| Bits per column and residue (1 byte) | 0.1 GB | 4.75 M columns x 20 |
| Label names, a cluster's accession, recall and precision by similarity, GA thresholds | under 1.5 GB, nearly all of it accessions | estimated |
| Manifest: Pfam and MGnify releases, SHA-256 of the models | nothing | |

So **about 40 GB for KO, CAZy and a few more, rows collapsed at 70%**.

- **The flanks are needed, whole proteins are not (measured).** With the domain's residues only, a best cluster is called for .9947 of real reads (.9991 with whole proteins). With 10 residues either side it is .9991 again
  (within 5% of the best for .9997 in both). Inserts add 2.8% to the residues (measured).
  The stored reads are the aligned part only; simulated reads that run past the domain need longer flanks (.982 at 10, .994 at 20, .996 whole).
- **Unlabelled rows must stay (measured).** An index of KO proteins alone sends reads of unlabelled proteins to labelled relatives. On the 21 genomes, with a .4 similarity floor, precision falls from .959 to .899 (recall .969 to .993);
  with .6, .981 to .961. So the "rows with a KO only" and "CAZy only" columns below are not safe sizes for an index: the full collapsed column is.
- **One Pfam release for everything.** The models reads are aligned to must be the models the rows were aligned to. The manifest's SHA-256 has to be checked at start. MGnify publishes 38.2 and production runs 37.4, so the rows have to be re-aligned to 37.4, or production moved.
- **Collapsing needs another tool at scale.** The greedy all-against-all used here will not run on a 14 M-row family. DIAMOND or MMseqs2 clustering at 70% would, with column identity as a check. Untested.

## The bundle, built (2026-11-24)

`bundle.py` builds, checks and uses it (its header has the format and the commands); `rs/src/bundle.rs` stores and reads the rows; `test_bundle.py` and `bundle_check.py` are its tests.
Format 3 (2026-11-25) adds members: with `--clusters` and `--collapse 0`, each family's block holds the representatives' rows and then the distinct member rows of their clusters, stored as differences from the representative (31.4 bytes a member row, measured), and a call names the cluster (REPORT.md, "Members second pass"). An annotation bundle collapsed at 70% has no members, as collapsing would merge each into its representative.

- **The same calls as the prototype:** 6,779 of 6,779 reads have the same score as the scan with 10-residue flanks, and name one of its best clusters (`bundle_check.py`).
- **Stored size (measured):** 88.8 bytes a row on the slice's 705,166 GA rows: 73 of residues at 5 bits and 16 of everything else, not the 10 assumed at first. So 102 GB for every cluster, **32 GB at 70%**, 20 GB at 60% (estimated), and about
  **40 GB in all at 70%** with kfpfam's files and labels.
- **Refused, with a message (tested):** models other than the bundle's, a family or a column the bundle's Pfam does not have, a file cut or changed, a Parquet whose metadata names other models. The pipeline's Parquets now carry `pfam_version`, and `profile` refuses a family of another version; for
  older Parquets it warns unless `--hmm` is given. `build` trusts that its hits were aligned to its `--hmm`.
- **Self-contained (tested here):** kfpfam's `filter join` ran on the bundle's own index and prep, and its candidates gave 98.0% of the GA alignments (the k4 T7.5 index). The pipeline's pyhmmer step ran on the same models (the `.h3m` of the bundle's HMM file), natively.

**Held-out genomes (measured).** 16 of the 21 genomes make the index (95,714 domain rows, 5,917 families, KO and CAZy); 12,000 windows of 50 residues from the other 5 are the sample. 3,983 align at GA.
Of the windows that align:

| Collapse | Rows | KO recall, floor 0 / .4 | KO precision | Profile Bray-Curtis | CAZy recall | CAZy precision |
|---|---|---|---|---|---|---|
| none | 95,714 | .951 / .906 | .944 / .958 | .056 / .070 | .927 / .865 | .949 / .965 |
| 70% | 72,776 | .948 / .900 | .942 / .957 | .059 / .074 | .927 / .865 | .944 / .962 |
| 50% | 53,291 | .934 / .878 | .933 / .948 | .073 / .093 | .927 / .857 | .942 / .959 |

**Reach: all 6,737 KO windows, aligned or not (70% collapse).** E is HMMER's at Z = 10^6 sequences, a family at a time.

| Alignments used | Windows aligned | KO recall of all windows, floor 0 / .4 | Precision |
|---|---|---|---|
| GA | 3,983 | .388 / .369 | .942 / .957 |
| GA, kfpfam's candidates | 3,903 | .379 / .360 | .942 / .957 |
| E <= 0.01 | 4,747 | .449 / .423 | .933 / .953 |
| **E <= 1** | 6,256 | **.539 / .501** | .912 / **.948** |
| E <= 10 | 7,763 | .574 / .532 | .869 / .936 |

- **Weak alignments take recall from .37 to .50 for one point of precision** (floor .4, E <= 1). The label check by similarity is what makes them safe. The pipeline now reports them with `pfam_weak_evalue` (off by default, uncommitted on `pfam-orf-filter-nail`).
  **The E-value's scale matters:** E <= 1 here is at Z = 10^6, a P-value of 10^-6. The pipeline's Z is its number of models, so the same cut is `pfam_weak_evalue = 0.025` on 37.4, 0.030 on 38.2 and 0.053 on 39.0.
- **What the pipeline gives (measured, its own Nextflow subworkflow with kfpfam's k4 T7.5 index):** recall of all KO windows **.468** at precision .949 with weak alignments, .360 at .957 with `--ga-only`. The .50 above is without a prefilter, so about 3 points
  are kfpfam's candidates. The k4 T6.5 and r999 indexes were not tried.
- **The floor of .4** is worth it with weak alignments. At GA alone it costs 5 points of recall for 1.5 of precision.
- Merged fragments longer than 50 residues should align more often. Not measured here. All five genomes are gut bacteria with relatives among the 16.

## Sizes (estimated)

Rows at GA over the release: 0.98 billion, 115 billion residues (about ±15%, a guess). Row header 8 bytes; no flanks or inserts, which add about 20% (the bundle above has them). *KO rows* and *CAZy rows* use the share of Pfam-hit proteins
with a label in the 21 genomes (57% and 5.1%): an upper bound, since 80% of MGnify representatives are partial and partial proteins miss KOfam's whole-protein thresholds.

| Rows collapsed at | Rows left | 5 bits | 1 byte | KO rows only, 5 bits | CAZy rows only, 5 bits | Labels, per scheme | Scan, ms a read domain | KO recall / precision |
|---|---|---|---|---|---|---|---|---|
| none (clusters) | 985 M | 80 GB | 123 GB | 45 GB | 4.1 GB | 2.0 GB | 1.9 | .987 / .988 |
| 90% | 631 M | 53 | 82 | 30 | 2.7 | 1.3 | 1.2 | .987 / .988 |
| 80% | 449 M | 38 | 59 | 22 | 2.0 | 0.9 | 0.9 | .988 / .987 |
| **70%** | 292 M | **26** | **40** | **15** | **1.3** | 0.6 | 0.6 | **.987 / .985** |
| 60% | 171 M | 16 | 25 | 9 | 0.8 | 0.3 | 0.3 | .982 / .979 |
| 50% | 84 M | 8 | 13 | 5 | 0.4 | 0.2 | 0.2 | .970 / .962 |

- *Rows left* is measured on the slice's 705,166 GA rows (64, 46, 30, 17, 9%), and scaled up. Very large families may collapse further.
- Dense column layout pads rows about 1.5×, as in `REPORT.md`.
- **The rest of the proteins**, for reads outside Pfam domains, would be an unanchored index, not sized here beyond this: 29% of KO-protein residues, on at most 23% of clusters, is roughly 19 billion residues (12 GB at 5 bits) plus its
  own seeds.

## How accurate is passing the label? (measured)

6,000 random 50-residue windows of the 21 genomes' proteins, each searched (exhaustive Smith-Waterman) against the other 20 genomes. The call is the commonest label among the tied best hits. Windows mostly inside a Pfam domain:

| Similarity to the best hit | KO recall | KO precision |
|---|---|---|
| ≥ .95 | .988 | .992 |
| .8–.95 | .994 | .989 |
| .6–.8 | .973 | .969 |
| .4–.6 | .921 | .891 |
| < .4 | .548 | .485 |

- **A similarity floor is needed, about .4:** below it, half of the KO calls are wrong. Windows outside Pfam domains are worse (KO recall .71, precision .53, against .89 and .86 inside).
- **CAZy:** .96–.99 recall and .96–1 precision from .4 up, on only 22–123 labelled windows a bin. Below .4: .44 and .40.
- **The Pfam family alone is not enough:** taking each family's commonest KO set gives 63–65% of domain residues, in-sample, and a domain residue's family holds 10 different KO sets on average. So cluster-level
  labels do add a lot, over a Pfam-to-KO table.
- The expected recall and precision in the size table weight these bins by the B11 reads' similarity to the targets left after collapsing. They assume the kept row carries the group's label, which collapsing inside
  each label would guarantee; I did not do that. B11's two hard libraries are further from MGnify, so worse than these averages.

## What other schemes would add

- **Per-protein labels, 2 bytes each:** eggNOG orthologous groups or COG categories, MEROPS, AMR, TCDB, NCBIfam: each about the KO figure, or less when sparse. Nothing here limits the number of schemes.
- **Derived labels need no row data:** EC numbers, KEGG modules and reactions come from KOs, GO from Pfam and KO. A table of a few MB each.
- **Domain-level labels** (dbCAN coordinates) need 4 bytes more a row to say where on the protein.

## Making it (estimated, not run)

- MGnify Proteins 2026_07 publishes **Pfam 38.2 and nothing else** (no KO, no dbCAN), and its Pfam release is not production's 37.4. So the labels have to be computed.
- **KOfam:** 27,576 models, 16.25 M columns, 3.4× Pfam 38.2's 4.75 M. At the `hmmsearch` rate measured on this laptop, 14,900 core-h for Pfam on all representatives becomes about **51,000 core-h for KOfam**, or perhaps half restricted to
  representatives with a Pfam hit (a guess), as 98.8% of KO proteins have one. kfpfam in front could cut it several-fold; it is untested on whole proteins and on KOfam.
- **dbCAN:** HMMs and DIAMOND, small beside it.
- **Partial proteins:** label a cluster by its longest or full-length member, not its representative. I could not measure how many partial representatives lose their KO.

## Limits

- 21 genomes, all human gut, isolates or good assemblies: more annotated than the average MGnify protein. Use the shares as upper bounds.
- The transfer test's targets are 20 genomes, far sparser than MGnify. Nearness is compared by similarity, not by where the hit comes from.
- No null: a read from an unlabelled source may be given a label. The precision above includes that for gut genomes only.
- KOfam has a threshold for 26,254 of 28,099 KOs.
- 40 families set the collapse ratios, at most 58,000 rows each.
