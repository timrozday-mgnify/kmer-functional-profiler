# Looking up MGnify proteins 90 from Pfam alignments: options, measurements, full-scale estimates

2026-11-24. Prototype, scripts and outputs: `workflows/pfam-anchor/lookup/` (kmer-functional-profiler), data in
`data/pfam-anchor/lookup/` (git-ignored, 2.3 GB). Nothing is committed. *Measured* is on this laptop; *estimated* is
extrapolation, and is labelled so.

**Status: the residue scan was chosen on 2026-11-24.** Annotated and collapsed indexes: `ANNOTATION.md`. **One very large family (PF00005, 14 M rows) was measured on the HPC on 2026-11-25:** recall .989 end to end, flat with family size; see that section.

**Corrections (2026-11-24, after the first version).**
- *Members were overstated.* The simulated reads were cut from member proteins that stayed in the members index. With the source protein left out, the source cluster is among the top ties for **.83 (BLOSUM62) to .90
  (identity)**, not .88–.95. The two-step route (top 500 clusters on representatives, then their members) gives .82–.89. Representatives (.63–.68) are unaffected. The gain from members is about 20 points.
- *The end-to-end recall is lower than the first stage's.* The scan was built as one program (`lookup scan DIR OUT 500:200`, representatives' GA rows, quarter bits, top 500 rows, column score, 200 rows to the banded score): **.9991 real (3 misses of 3,399) and .9962 simulated**.
  With only 50 rows to the banded score, .9982 and .9929. Adding rows under GA: .9994 and .9988 (200), .9997 and 1.0000 (500).
- *Cost is now measured on that program, nothing else running, one M4 core:* the scan is 1.37 ns a row (33 GB/s), but choosing the top 500 rows adds 1.1 ns a row, so **2.5 ns a row**, plus 43 µs columns and 410 µs banded a read domain.
  At 0.77 M rows that is about **2.4 ms** (estimated), not 1.5. A library of 0.34–0.99 M domains needs 0.2–0.65 CPU-h, and 10× depth 2–6.5 CPU-h.
- *The earlier timings, and the follow-up stage times, were taken while other jobs ran* (including swapping); only the numbers above were re-timed alone.
- *A self-contained index holds domain residues, not proteins.* With 10 residues either side of the domain the calls on real reads are unchanged (.9991); see `ANNOTATION.md`.
- *MGnify's release Pfam is 38.2* (its README), so its table does not match production's 37.4.

## Bottom line

- **Recommended: score the read's residues against every cluster of its Pfam family by HMM column, weighting each match by its
  information bits, and take the best 100–500 clusters (or those within 20 bits of the best); then a banded alignment.** Measured on
  3,399 real read domains, the best target is kept for .9985 (top 100) or .9997 (top 500) after the first stage; end to end, .9991 is called (Corrections). No seed index is needed: the index is the aligned residues themselves.
- **Seeds look hard to justify at MGnify's scale (estimated).** On top of the residues they add 400–860 GB (110–215 GB if every 4th column is indexed) and leave 1.3–28% of a family to
  check. As prototyped they cost 0.7–14× the scan's CPU. If checking a candidate were 4× cheaper than measured, 10- and 14-bit seeds would cost 2–4× less than the scan, and lose 1.2–4.2% of best targets.
  All of that loss is in reads under .95 similarity to MGnify: from .95 up, the 8-, 10- and 14-bit and 4-residue Golomb seeds kept every best target (2,368 reads).
- **A global bit threshold is poor, a per-target one is better, and a per-read margin is best.** Same recall (about .999): a global
  T passes 9,000 clusters, a per-target fraction of its own bits about 1,000, a margin of 20 bits below the read's best about 220.
- **Index the representatives' GA domains, and decide about members.** Representatives alone put the source cluster among the top
  ties for only 63–68% of reads from member proteins; adding members gives 83–90% (source protein left out). Held as differences from the representative,
  members are small (4.9 differing residues per row on average).
- **Decision needed:** see the last section. The main unknown is at real scale: the typical family a read lands in is 17× (median) to 47× (mean)
  the average one tested here, and the largest about 860×.

## What exists

- **The lookup already exists** in `kmer-functional-profiler/workflows/pfam-anchor` (`anchor.py`, plan in
  `docs/pfam-anchor-plan.md`). kfpfam's `recruit/` is the step before it (weak reads to families), and writes the alignments that feed it.
  Its method: the 16 highest-entropy columns, keep the top 1% of clusters by matches, then column identity. It tested seeds, found
  anchor residues better, and replaced them. Its largest family had 3,090 clusters, and its release index was sized only roughly.
- **The pipeline stores** one row per domain in `*_pfam_alignments.parquet`: `hmm_from`, `cigar` (M/I/D) and `residues`, so each
  residue's column follows directly. The model version is not stored.
- **Release mismatch found.** The local B11 Parquets and the local Pfam HMM disagree for at least PF12715 (a read reaches column 388;
  the local model has 379 columns) and PF07580. So I re-aligned the read residues to the local HMM (3,399 of 3,463 stay at GA). On Codon, check that
  the Parquets' release, the target alignments' release and 37.4 are the same: a wrong version does not fail loudly.

## Data and truth

- **Targets.** The MGnify90 slice of 630,402 clusters (4.2 M proteins): every cluster with an annotated member of 40 Pfam families, so those
  families are complete. All members searched against the 40 models: 980,579 representative rows (705,166 at GA), and 2,022,097 distinct
  GA member rows. 6,000–58,000 representative rows per family (PF12661 has 250,000, almost all weak).
- **Queries.** 3,399 real read domains from the 11 B11 libraries (mean 47 residues, 37 families, from the pipeline's stored aligned
  residues), and 3,380 windows of 50 residues cut from non-representative members, mutated to 100–60% identity (BLOSUM62, indels). Real reads
  span .4–1 similarity to their best cluster; ERR4334731 and ERR5866877 are the hard ones (27% and 37% of reads at ≥ .95; the other nine 73–99%).
- **Truth.** Smith-Waterman (BLOSUM62, 11/1) of every query against all 630,402 representatives (187 M residues, 6×10¹³ cells, 30 min). The
  matrix and 200 scores agree with Biopython. *Best* is recall of the family's best target; *near* is within 5% of its score.
- **The best target is always within the read's family** in this slice (all 6,779), but that is by construction: the slice holds only clusters of
  these families. Proteins with no Pfam hit cannot be reached by any version of this method and are not in this test. The release has 33.5% of clusters
  with a Pfam hit on any member.

## First stage: which clusters to look at (measured, real reads)

*Candidates* is the mean number of clusters left, of a family of 16,330 on average; the share is constant over 1–16× subsampling for
seeds. *Residues* are scored on the representatives' GA and weak rows. Recall brackets are a 95% bootstrap over families. *Hard* is ERR4334731 and ERR5866877 (1,183 reads).

| First stage | Index, postings per residue | Best kept | 95% CI | Hard libraries | Simulated | Candidates | Share of family |
|---|---|---|---|---|---|---|---|
| k-mer 5 | .90 | .9950 | .991–.997 | .987 | .996 | 2,471 | 15% |
| Golomb `1100101` (4 residues) | .87 | .9997 | .999–1 | .999 | 1.000 | 4,636 | 28% |
| — every 4th column | .22 | .9915 | .987–.995 | .976 | .986 | 1,884 | 12% |
| Golomb 5 residues `110010000101` | .80 | .9953 | .991–.998 | .987 | .996 | 2,066 | 13% |
| Info-bit seed, 8 bits | .93 | .9962 | .993–.999 | .989 | .991 | 1,051 | 6.4% |
| Info-bit seed, 10 bits | .91 | .9876 | .981–.994 | .965 | .973 | 523 | 3.2% |
| Info-bit seed, 14 bits | .88 | .9582 | .938–.975 | .883 | .927 | 217 | 1.3% |
| Per-target bit threshold (θ 16) | .93 | .9985 | .997–1 | .996 | .999 | 1,758 | 11% |
| Residues, all columns, bits: **top 100** | — (matrix) | .9985 | .997–1 | .996 | .997 | 103 | 0.6% |
| — **top 500** | — | .9997 | .999–1 | .999 | 1.000 | 503 | 3% |
| — **within 20 bits of the best** | — | .9991 | .998–1 | .998 | .999 | 220 | 1.3% |
| — ≥ 0.3 of the target's own bits | — | .9988 | .998–1 | .997 | 1.000 | 1,028 | 6% |
| — global bits ≥ 20 | — | .9932 | .982–1 | .998 | .998 | 2,701 | 17% |
| — global bits ≥ 10 | — | .9979 | .993–1 | 1.000 | 1.000 | 9,015 | 55% |
| 16 highest-entropy columns, top 1,000 | — | .9991 | .998–1 | .998 | 1.000 | 1,565 | 10% |
| 16 highest-entropy columns, top 1% (existing plan) | — | .9944 | .989–.998 | .984 | .996 | 283 | 2% |

Best kept by the read's similarity to its best cluster (real reads):

| First stage | ≥ .95 (2,368) | .8–.95 (733) | .6–.8 (254) | .4–.6 (43) |
|---|---|---|---|---|
| Golomb `1100101` | 1 | 1 | 1 | 1 |
| Info-bit 8 bits | 1 | 1 | .972 | .86 |
| Info-bit 10 bits | 1 | .990 | .929 | .63 |
| Info-bit 14 bits | 1 | .950 | .713 | .28 |
| Residues, top 100 | 1 | 1 | 1 | .88 |

- **Info-bit seeds are the best seeds**: 8 bits keeps .996 at 1,051 candidates, where a 5-mer needs 2,471, a Golomb ruler of 5 residues 2,066, and two hits of the Golomb 4-residue seed 2,130 (.998).
- **Golomb rulers are only slightly better than contiguous words** of the same weight: at 5 residues .9953 and 2,066 candidates against .9950 and 2,471; at 4 residues .9997 and 4,636 against .9985 and 4,756;
  at 6 residues .982 and 907 against .988 and 1,384. Extending a ruler until it reaches T bits is worse than consecutive columns.
- **Sampling every 4th column** cuts a seed index 4× for about 1 point of recall.
- **Per-target thresholds for seeds gain nothing clear.** θ 12 keeps .988 at 747 candidates, where one 9-bit threshold keeps .993 at 711. θ 16 (.9985 at 1,758) lies on the line between 8 bits
  (.996 at 1,051) and 6 bits (1.0 at 2,607). For the residue scores they do matter (Bottom line), and a margin below the read's best does better still. A log-likelihood-ratio score of matches and mismatches
  keeps at most .984 at any threshold, as it drops diverged targets.
- **Posting lists of rare residues only** (≥ 3 bits, 37% of the postings): .981 at 524 candidates, .962 at 1,479 with a fixed 8-bit sum. Worse than scoring every residue.

## Follow-up stages (measured, on residue top 500 and 8-bit seeds)

Keep the top *n* by each stage, or everything within *d* of the best stage score; *kept* is the mean number left; then the final Smith-Waterman picks.

| After the first stage (top 500) | Rule | Best kept, real / simulated | Kept | Cost per row |
|---|---|---|---|---|
| Anchored full diagonal in columns (BLOSUM62, all columns) | best −30 | .9991 / .9982 | 59–95 | 0.08 µs |
| — | top 50 | .9985 / .9959 | 55–63 | |
| — | best only | .9603 / .9254 | 10 | |
| Window ±8 columns, on seed hits (8-bit seeds) | best −30 | .9926 / .9876 | 149–231 | 0.25 µs |
| Ungapped sequence diagonal | best only | .9962 / .9805 | 9 | 0.41 µs |
| **Banded gapped (±8) on the seed or column diagonal** | best only | **.9997 / 1.0000** | 9 | 2.2 µs |
| Smith-Waterman, whole protein | best only | .9997 / 1.0000 | 9 | 54 µs |

- **The banded gapped score is the final score**: among top-500 candidates the target with the best Smith-Waterman score always also has the best banded score, real and simulated (with 8-bit seeds, for all but 2 of 3,380 simulated reads). So whole-protein
  Smith-Waterman, 24× the cost, adds nothing on these short reads.
- **Use: columns → banded gapped.** The window stage costs 3× the column sum and keeps less. The ungapped diagonal fails on indels (simulated, .98).
- **Cost of the whole follow-up at top 500:** 40 µs for columns, then 50 rows × 2.2 µs = 110 µs. About 0.15 ms per read domain.

## Accuracy

- **Simulated reads from member proteins: is the source cluster among the top ties?** Exhaustive Smith-Waterman against representatives: .64. Column scores on representatives: .63–.68. On every member's
  row: .88 (BLOSUM62), .89 (bits), .95 (identity); the median tie is 2 clusters in each. Representatives lose reads from members, whose differences from the representative are the reason another cluster's representative can win.
  By identity (column identity), members fall from .98 at 100% to .83 at 60%, representatives are .70 and .67.
- **A two-step route works**: the top 500 clusters by bits on representatives hold the source cluster for .987 of reads, so members need only be compared there. Not built.
- **Clusters whose representative has no domain** (12% in the earlier 40-family plan) are invisible to a representatives-only index. The best Smith-Waterman target was in reach of representatives at GA
  for .9997 of real reads here, but the slice was chosen by family, so it is not a fair test.
- **Zymo.** 1,569 reads from the 6 mock runs fall in these 40 families; 1,488 have a gene of the standard (aligned residues ≥ 95% identical).
  Right, at the gene's best cluster, is .83 (per run: .77–.88; per read credit .40 because ties are many), **identical for every method, and for exhaustive search**.
  94% of these reads are ≥ .95 similar to a cluster, so Zymo measures how well a tie is resolved, not the first stage. It cannot rank the options above.
  Of 124 gene domains, 112 have a cluster over .9 similar (median .99). The log sample ERR2935805 has only 16 reads in the 40 families, so
  use SRR23926920 (358 reads with a gene; .77 right, 12 clusters tied at best) or the D6300/D6305 runs for a resolution test.
  Zymo for all families needs the genes searched against all of MGnify90 on Codon.

## Full scale: MGnify proteins 90 (estimated)

**Model.** Release rows = 1.66 B representatives × the 148,471-representative sample's domains per representative (full Pfam 38.2 at GA,
kfpfam's `rest.domtbl`), corrected for the sample being 12% longer than the release. Rows per family come from the same sample; each read domain is
scanned against its family; the family's share of the reads is B11's (6.4 M domains, 11 libraries).

| Quantity | Value | Status |
|---|---|---|
| Representative domain rows at GA | 1.0 billion (0.98) | estimated; ±15% is my guess (the sample has 38% of reps with a domain, the release 33.5% of clusters with a hit) |
| Residues in those rows | 115 billion | estimated |
| Rows a read domain meets: mean / median / 90% / 99% | 0.77 M / 0.27 M / 1.7 M / 14 M | estimated |
| Largest family | about 14 M rows | estimated from 1 family's share; PF00005 has 62.6 M domains in the release (all proteins) |
| Families unseen in the sample | 6% of read domains | estimated |

**Index and CPU per read domain, M4 core.** *Scan* is measured on a memory bench (27 GB/s, one thread) and multiplied by 0.77 M rows. *Seed* is lookups
(200 ns each, a guess) + 2 ns a posting (a guess) + checking the surviving rows at 80 ns each (measured, scalar), or in brackets at 20 ns (a guess at an optimised check).

| Option | Extra index | Candidates (rows checked) | CPU first stage, ms | Real recall (measured) |
|---|---|---|---|---|
| Residue scan, all 44 columns | none beyond the matrix | 500 | 1.3 | .9997 (top 500) |
| Residue scan, 16 columns, top 1,000 | none | 1,000–1,600 | 0.5 | .9991 |
| Info-bit seed 8 bits | 430–860 GB | 49,000 | 4.3 [1.4] | .996 |
| — every 4th column | 110–215 GB | 27,000 | 2.2 [0.6] | .985 |
| Info-bit seed 10 bits | 420–840 GB | 25,000 | 2.2 [0.7] | .988 |
| Info-bit seed 14 bits | 400–810 GB | 10,000 | 1.0 [0.4] | .958 |
| Golomb `1100101` | 400–800 GB | 218,000 | 18.6 [5.5] | .9997 |
| Per-target bit threshold (θ 16) | 430–860 GB | 83,000 | 7.2 [2.2] | .9985 |
| Golomb 5 residues | 370–730 GB | 97,000 | 8.3 [2.4] | .995 |
| k-mer 5 | 415–830 GB | 116,000 | 10.1 [3.1] | .995 |

- **Index sizes** are postings (residues × postings per residue) × 4, 6 or 8 bytes (low, used, high): short lists compress badly, so a key directory is in
  there only as a guess. The seed options also need the residues for checking, on top.
- **The matrix** (needed either way): 115 B residues at 1 byte = 123 GB with a 8-byte row header; 72 + 8 = 80 GB at 5 bits. Laid out by column for the scan, rows that cover part of a
  model pad to 1.5× on average (101 M dense cells for 68 M residues here), so about 180 GB at 1 byte or 115 GB at 5 bits if every column is stored for every row. *Estimated.*
- **Members** add about 2.5 B distinct rows (0.59 distinct of 4.2 B GA member domains, as in the slice), at about 12 bytes each if held as differences from the representative
  (mean 4.9 differing residues per non-representative row, measured; 17% of those rows equal the representative's), so **about 30 GB (20–50)**; stored in full, 3.4× the representatives' residues (measured on the slice, whose clusters are bigger than the release's).
- **Memory at search.** One family at a time, read domains batched by family: the largest family is under 2 GB of residues (14 M × 130 columns); the rest stays on disk. Serving
  from RAM would need about 130–180 GB.
- **Reading.** Every library hits nearly all large families, so one pass reads the whole matrix (2–4 min at 0.5–1 GB/s), once per library or shared by several. Random access suits seeds, and
  a seed index of 400–860 GB cannot be held in RAM or on this laptop's disk (14 GB free).
- **CPU per library.** 0.34–0.99 M domains per library at the 2 M-read cap (nail arm, measured), so **0.1–0.4 CPU-h for the scan** and 0.01–0.04 for the follow-up (M4 core). At 10× depth 1.4–4 CPU-h per
  library (estimated), against 10 for kfpfam r999 on 37.4 (114 over 11 libraries; 18 on 39.0). On Codon the per-core memory bandwidth is the unknown: assume 2–3× slower.
- **Build.** The matrix is a by-product of aligning representatives to Pfam. Searching all representatives with `hmmsearch` is 14,900 core-h (plan, laptop-measured rate), all proteins 63,600.
  Aligning only the pairs MGnify's own Pfam table holds (110 GB, scores and coordinates, no alignments) is far cheaper but untried; its Pfam version must match the pipeline's.
  A new release or Pfam version means rebuilding: 37.4, 38.2 and 39.0 differ for nearly all families.

## One very large family at full size (measured on Codon, 2026-11-25)

`hpc/kfp-ablations/pfam-anchor-large/run.sh` with `FAMILY=5`; outputs copied back to `data/pfam-anchor/lookup/large/PF00005/`.

- **Data.**
  - *Family and release:* PF00005 (Pfam 39.0, PF00005.34), against MGnify release 2026_07. 59.3 M proteins carry the family; 11.5 M of the release's 1.66 B cluster representatives do. Those give 14.3 M domain rows, 12.2 M at GA.
  - *Queries:* 2,176 real read domains at GA from the 11 B11 libraries, 200 a library (kfpfam top arm, Pfam 39.0, re-aligned to the same model), mean 41 residues. Each read is also reversed, as a null.
  - *Truth:* exhaustive Smith-Waterman against every representative, and against 1 in 4, 16 and 64 of them (the subsamples).
  - *Scan:* GA rows only, domain ±10 residues (`FLANK=10`), as a self-contained index holds.
- **Size.** 1.56 B stored residues: **1.1 GB at 5 bits, 1.7 GB at 1 byte**. This matches the estimate for the largest family (about 14 M rows, under 2 GB; Full scale).
- **Recall is about 1 point below the slice's, and does not fall with family size.** Best cluster called (*right*), and within 95% of the best score (*near*):

  | K:SECOND | Full: right | Full: near | Full: missed | 1 in 4: right | 1 in 16: right | 1 in 64: right |
  |---|---|---|---|---|---|---|
  | 100:200 | .986 | .998 | 30 | .973 | .975 | .972 |
  | **500:200** | **.989** | .997 | 25 | .985 | .988 | .988 |
  | 2000:200 | .990 | .997 | 22 | .987 | .989 | .989 |
  | 2000:1000 | .994 | .999 | 13 | .993 | .992 | .994 |

  At 1 in 64 (about 180,000 clusters and 220,000 rows, 4× the slice's largest family in rows) recall is already .988 at 500:200. So the gap from the slice's .9991 comes with this family or its reads, not with N, over a 64× range up to 14 M rows. Per library, 500:200 ranges .970–1.000.
- **Misses are near ties.**
  - *At 500:200:* 17 of the 25 misses are 1–4 score units under the best; 2 are far off (30 and 40 short).
  - *What they share:* shorter reads with lower best scores (medians 38 residues and 169, against 39 and 188), and a single best cluster (median ties 1, against 4).
  - *Cheapest fix:* more rows to the banded score. 2000:1000 halves the misses to 13.
  - *Untested cause:* the scan sees domain ±10 residues, while the truth aligns to whole proteins.
- **Null: reversed reads stay apart from real ones.** At full size the reversed reads' best score has median 63, 99th percentile 82 and maximum 93; every real read scores at least 104. The null maximum grows with N: 80, 83, 87 and 93 from 1 in 64 to full. At 1 in 64 the lowest real read (82) is only 2 above the null maximum, so the floor has to grow with family size.
- **Cost on x86** (Codon; 12.2 M rows a read domain):

  | Run | Scan | Choosing top rows | Columns + banded | A read domain | Wall |
  |---|---|---|---|---|---|
  | 1 thread, 500:200, 200 reads | 1.48 ns a row | 0.66 ns a row | 0.9 ms | **27 ms** | |
  | 16 threads, 500:200, 4,352 reads | 3.14 ns a row | 0.83 ns a row | 1.6 ms | **50 ms of CPU** | 3.1 ms |

  - *One core* is as fast as the laptop's M4 (2.1 against 2.5 ns a row).
  - *16 threads* are 1.9× slower per core. Memory bandwidth, or other jobs on the node, could explain it: the job could not be `--exclusive` there.
  - *2000:1000:* the banded score takes 4.9 ms of CPU a read domain, against 1.2 at 500:200.
  - *Worst case:* this is the largest family, about 11× the 2.4 ms estimated at the mean family size (0.77 M rows).
- **Not measured.** The truth step's runtime (the jobs shared one log name, so only the report's survived), PF00072, and Pfam 37.4.

## Uncertainties, and what I could not test

1. **Family size.** The tests stop at 58,000 rows per family; the median read domain lands at 0.27 M and the mean at 0.77 M, against a mean of 16,000 here: 17× and 47×, and 860× at the largest (14 M). Over a 16× range of N (1/16 subsamples), top-100 recall did not fall (.9991, .9982, .9985), and within-20-bits
   candidates grew as N^0.46 (61, 121, 220), so about 1,300 at 0.77 M if that continues, which is far past the tested range. The share of a family a
   seed passes is constant by construction in subsampling (it removes units at random). *Larger families are more diverse, so real shares may fall or rise.* This is why I recommend an HPC
   test of one very large family. *Tested since (2026-11-25):* in PF00005 at 14 M rows, recall is .989 at 500:200 and the same from 1 in 64 to full size; it is lower than here for that family
   or its reads, not for its size ("One very large family").
2. **Reach.** Only reads with a Pfam GA hit and only clusters with a Pfam hit. Reads outside domains, or clusters with no domain, get nothing.
3. **No null.** I have no false-assignment rate: a read whose source is absent from MGnify still gets its nearest cluster. A reversed-read null is the way to set a floor. *Measured since in PF00005:* reversed reads reach at most 93, real ones at least 104, and the null maximum grows with family size.
4. **Read truth.** Real reads are the aligned part only, so nothing here tests extension past the domain. Truth is the best Smith-Waterman score; several clusters often tie.
   The families are mid-sized, the reads are 47 residues on average, and 37 families carry the whole real set, so a family-level CI is wide.
5. **Emulation.** No x86 timing here (since measured on Codon for PF00005: 2.1 ns a row on one core, the same as the M4). The scan is on an M4 Pro, the stage times are scalar Rust in four threads with swap in use. Treat them as ±2×.
6. **Pfam.** Everything here uses the local HMMs, labelled 38.2 (the B11 Parquets are labelled 38.2 too, and differ for at least two families); production is 37.4. Sizes scale with model columns and domain counts; the choice between methods should not change.

## Reproducing

```
workflows/pfam-anchor/lookup/run.sh          # targets queries pack truth cand sub follow tables
```

Runs one step at a time on this laptop. Wall: search 13 min, pack 1, truth 8 + 30 min, options 3, follow-up 23 min (most of it full Smith-Waterman), tables 1. Peak 2.7 GB for a Rust step (Python steps not recorded);
the laptop was swapping (8.0 GB of 8.7 used at the end), so timings may be inflated. `rs/` is the prototype, `prep.py` inputs, `analyze.py`, `stages.py`, `zymo.py`, `members.py`, `sizing.py`, `extrapolate.py` the tables,
`check_truth.py` checks the truth against Biopython. Not optimised, as asked. Each step was run on its own in this order; `run.sh` has not been run from scratch in one go. The local B11 Parquets are the nail arm.

## Decision

1. **Decided: the residue scan** (representatives' GA rows, per-read margin or top 500, then columns, then a banded alignment) as the method, with no seed index. Not taken: the 16-column scan as a first
   pass (3× less CPU, the same recall within error, but the number it must keep at release family sizes is untested), or 8-bit seeds every 4th column if RAM-resident lookup must be fast and 1–2% of best targets may be lost.
2. **Representatives only, or representatives plus members as differences** (+20–50 GB estimated; the source cluster is among the top ties for .83–.90 of reads from members, against .63–.68)? *Recommend adding members as a second pass on the top 500 clusters.*
3. **Run on HPC** (`hpc/kfp-ablations/pfam-anchor-large`): one very large family at full size, to check recall and CPU past 58,000 rows. *PF00005 done 2026-11-25* (Pfam 39.0): recall holds with
   family size, at about 1 point below the slice. **Decided (2026-11-25): 2000:1000 is the default** (.994 on PF00005, 13 misses against 25; the banded score 4.9 ms of CPU a read domain there, against 1.2); `bundle.py profile` uses it. Still open: PF00072; Pfam 37.4's alignments; a wider flank; and the Zymo
   genes against all of MGnify90 if a full Zymo truth is wanted.
