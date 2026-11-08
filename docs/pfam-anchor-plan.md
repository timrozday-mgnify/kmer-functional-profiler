# From raw-read Pfam hits to MGnify90 units

Plan, 2026-11-06 (the user's request): assign a read to MGnify90 units through the HMM alignments the read and the units share, with k-mers keyed by HMM column. Includes a local test case for development. Partly tested: phase 1b ran on 2026-11-06 with windows of real proteins, and phase 0's test case was built on 2026-11-08, with simulated reads and DIAMOND's best hits, and phase 1's oracle and phase 2's seeds passed their gates the same day (Progress log). The comparison with DIAMOND and the tree index on speed, and read-level precision, are still to do.

## Verdict

Worth a short local test before building anything. The sizing favours it, and there are three known limits:

- **The search space falls sharply.** A read's Pfam hit picks a family and its HMM columns, so a k-mer only has to be unique among the residues at those columns. On the test case, 6-mers keyed by column have 1.82 units per key.
- **Short k-mers become possible.** Even k = 4 has only 2.65 units per key, so k-mers short enough to survive divergence may stay selective.
- **Reach is the first limit.** Only units with a Pfam hit can be reached: 33.5% of the release's units at GA. A read must also fall inside a domain: GA domains of the 40 test families cover 43% of the residues of proteins with a hit, though domains from other families will raise this.
- **Resolution is the second limit.** Even an error-free 40-column window is unique to one unit only 59% of the time (p90 8 units, p99 147). Many reads will tie across units, however good the k-mers are, so ties must be kept and resolved downstream.
- **Alignment agreement looks fine so far.** Column identity is only as good as the two HMM alignments. In phase 1b, 99% of a 40-column fragment's residues landed on the same column as in the full protein's alignment, and calls barely changed. Fragments had no sequencing errors or real ORF calling, so phase 1 still has to confirm this.
- **Brute-force AAI is probably enough for most families.** Comparing a read against every unit of its family should cost less than the Pfam search itself, even at release, so k-mers may only be needed for the largest families (phase 1b).

## The method

**Common template.** An HMM match column is a shared coordinate. A CIGAR maps residues of a read or protein onto those columns, so a read and a unit can be compared column by column with no read-to-unit alignment.

```text
Build   protein domains (family, hmm_from, ali_from, cigar) + sequences
          -> (family, column, k-mer) -> units

Query   read -> six-frame ORFs -> Pfam hits with alignments
          -> (family, column, k-mer) of each hit -> candidate units
          -> column score per candidate -> best units, ties kept
```

**K-mers.** A k-mer is k residues at consecutive match columns, keyed by its first column. Inserted residues are skipped and a deleted column breaks the k-mer. Matches are exact against the key, with no neighbourhood list. Seeds need not be consecutive: any set of columns can form a key, and phase 1b tests which sets work best.

**Column sampling.** Keys are positions, so the index can keep every s-th column, with no hash sampling. A query uses all of its columns, so any run of k + s − 1 matching columns is guaranteed a seed.

**Score.** Candidates are compared by plain identity over all the read's columns, so a unit missing some of them loses (phase 1; identity over shared columns only, phase 1b's choice, favoured partial domains). A rarity likelihood ratio from the HMM's emissions is the likely tie-break (phase 1b). Reads keep every family they hit, as in [recall-first pileup](kmer-functional-profiler-plan.md): candidates are not cut by a presence call, and related families are not collapsed.

## What is measured so far

Local test case below, 2026-11-06. The 40 HMMs against the mini-slice took 24 s on 4 threads.

| Measure | Value |
| --- | --- |
| Proteins, units | 294,822, 55,433 |
| Proteins and units with a GA domain from the 40 families | 183,652, 49,673 |
| GA domains, and all domains reported at `-E 10` | 200,134, 307,698 |
| Mean HMM columns per GA domain | 102 |
| Share of residues in GA domains, proteins with a hit | 43% |
| Median units per family, median HMM length | 1,043, 138 |
| Units with a GA domain whose representative has none | 6,145 (12%) |

Column k-mers, GA domains, every column (s = 1), each unit's members merged:

| k | (key, unit) postings | Keys | Units per key | Share of postings from representatives alone |
| ---: | ---: | ---: | ---: | ---: |
| 4 | 5,277,454 | 1,994,113 | 2.65 | 0.82 |
| 6 | 5,205,113 | 2,852,280 | 1.82 | 0.80 |
| 8 | 5,111,374 | 3,258,683 | 1.57 | 0.78 |

Expected chance of at least one seed, calculated with independent column matches. Real conservation is not independent, which is why the test case has real relatives as well as simulated strains.

| Aligned columns, k | 95% | 90% | 80% | 70% | 60% |
| --- | ---: | ---: | ---: | ---: | ---: |
| 30, k 4 | 1.00 | 1.00 | 1.00 | 0.97 | 0.86 |
| 30, k 6 | 1.00 | 1.00 | 0.93 | 0.72 | 0.43 |
| 30, k 8 | 1.00 | 0.96 | 0.73 | 0.40 | 0.16 |
| 40, k 4 | 1.00 | 1.00 | 1.00 | 0.99 | 0.93 |
| 40, k 6 | 1.00 | 1.00 | 0.97 | 0.83 | 0.54 |
| 40, k 8 | 1.00 | 0.99 | 0.84 | 0.51 | 0.22 |

Rough release size, assuming the 105 postings per unit seen here hold for all 556 M units with a GA hit: 5.8×10¹⁰ postings at k 6, s 1, and about 1.9×10¹⁰ at s 3. Singletons and short domains make this uncertain; the 1-in-100 `pfam-reannotate` run replaces it (phase 5).

## Local test case

Everything runs on the laptop, one job at a time, on capped threads, in about 30 minutes in total.

| Piece | Source | Notes |
| --- | --- | --- |
| Reference proteins | `data/mini-slice/members.parquet` | Family-complete, so near units are present. Git-ignored; comes from `results-hpc/tree-index-slice` through `slice.py mini` |
| Families | `results-hpc/tree-index-slice/results/families.tsv` | 40 families, cut from `~/Downloads/Pfam-A.hmm.gz` (Pfam 38.2) |
| Reference alignments | `hmmsearch`, then `workflows/pfam-reannotate/reannotate.py parse` | 24 s; the same hits format as the HPC run |
| Simulated strains | `spike.py simulate`, as in `tree-index/build.nf` | 2,000 strains, 100–70% identity, 150 bp pairs, flanks, per-pair coding truth (`reads.parquet`) |
| Real relatives | `recall.py holdout` | Proteins held out by identity to their nearest relative, removed from the reference. Real conservation, which uniform simulation lacks |
| Read ORFs | Six-frame translation, split at stops, at least 25 aa | Follows kfpfam's `measure/real_orfs.py`; 25 aa kept .998 of GA recall there |
| Read hits | `hmmsearch --notextw -E 10 --domE 1000000 -Z 1000000 --domZ 1000000` | The 40 HMMs for iteration; full Pfam-A once, for wrong-family calls |
| Null | The same reads reversed | Sets the read-level bar, as in the tree-index benchmarks |
| Baselines | DIAMOND read to all members; tree index tier 1 and 2 on the same reads | DIAMOND from the biocontainer through Docker; `build.nf` on the mini-slice |

Two local limits carry over from the tree-index work:

- **Chance hits scale with the database.** The mini-slice holds 295 k proteins rather than 5.7 B, so precision is optimistic. Local results rank methods; HPC sets the bars.
- **Families are preselected.** Reads from these families against these HMMs say nothing about reads with no family. Full Pfam-A on the reads, and the reversed reads, cover part of that.

A unit test fixture also exists: `tests/data/mini_reannotate` and `tests/data/mini_fmh/Pfam-mini.hmm` (6 HMMs, 84 proteins).

## Metrics

Reuse the tree-index definitions (`recall.py reads`): a pair is *called* when its best score clears the reversed-read bar, and *right* when any tied best is the source's unit. Report by identity bin, per strain with intervals.

Each read's outcome is split, so a failure can be traced to one stage:

1. **Hit:** the read has a hit on any family of its source protein.
2. **Candidate:** the source unit is among the candidates, given a hit.
3. **Right:** the source unit is among the tied best, given a candidate.
4. **Overall:** read recall and read precision by unit, plus median tied units per call.

Stratify by the pair's coding share and by whether the source position lies in a domain. Also report postings, query time per read and peak memory.

## Phases

Each phase records its outcome in this plan's Progress log, and only local work is planned before phase 5.

### Phase 0: test case

`workflows/pfam-anchor/prepare.py`, with its steps as subcommands:

1. Cut the 40 HMMs, annotate the mini-slice, and write `domains.parquet`.
2. Produce strains and held-out proteins, reads, and reversed reads, reusing `spike.py` and `recall.py`.
3. Translate reads and search them. Add `orf_id` to read, mate, frame and start mappings, since `reannotate.py parse` needs integer names. Alternatively, kfpfam's `recruit/pairs.py` writes the alignments with string ORF names (see "Read hits from kfpfam").
4. Run DIAMOND once for the method-independent best identity of each read.

**Done when:** one command rebuilds the case, and `hit` recall by identity bin is reported.

*Done 2026-11-08 (`workflows/pfam-anchor/prepare.py all`); see Progress log. Full Pfam-A on the reads is not run yet: it is only needed for wrong-family calls (phase 3).*

### Phase 1: oracle, no k-mers

`anchor.py oracle`: every read hit against every domain of that family, by column score. This is exhaustive, but fine at this size.

1. Rank units for each read and report `right` given a hit.
2. Report the same reads through DIAMOND, restricted to those units, to isolate alignment disagreement.
3. For pairs where both methods align the read to the same protein, report how often the two alignments place residues differently.

**Gate:** given a hit on the source family, oracle `right` is within 0.05 of DIAMOND's at each identity from 100% to 80%, on both simulated strains and real relatives. If not, the shared template is losing the information, and k-mers can't recover it. Stop, or fix the alignment step first.

*Passed 2026-11-08 (`workflows/pfam-anchor/anchor.py oracle`), with identity over the read's columns; see Progress log. The oracle ties more units than DIAMOND, so its credit is lower.*

### Phase 1b: brute-force AAI, HMM weights and discriminating columns

Added 2026-11-06 at the user's request. Runs on the same families and needs no index.

1. **Is AAI against every possible match feasible?** Measure comparisons per second, and from the release's family sizes work out the cost per read hit. If it is cheap enough, k-mers may only be needed for the largest families.
2. **Can AAI be weighted per column from the HMM?** Compare plain identity with:
   - conservation weights (information content from the HMM's match emissions);
   - the inverse, favouring variable columns;
   - rarity, where a match on a residue that is rare at that column counts for more (from HMM emissions or from the units);
   - a likelihood ratio that also penalises mismatches;
   - BLOSUM62, as the alignment-free reference.
3. **Which residues and columns discriminate?** For each column, measure how well it separates units (entropy across units) and how often it survives between a unit and its nearest relative (retention). The best seed columns need both. Report these by HMM conservation and by amino acid.
4. **Which combinations discriminate?** Measure how much columns overlap in the information they give, how that varies with distance along the model, and how many columns it takes to single out a unit.
5. **Which seed shape narrows candidates best?** A seed is a set of columns, not necessarily adjacent, compared on recall against candidates per read:
   - contiguous k-mers, as the baseline;
   - fixed spaced patterns;
   - the top k columns in a window, ranked by entropy, conservation, retention or entropy per unit of retention lost;
   - columns chosen greedily, each adding the most new information per unit of retention lost;
   - random columns, as a control.

**Tests:** reads are windows of real proteins held out of the reference: members of a known unit (truth: their own unit) and whole held-out units (truth: the nearest remaining unit). Column statistics come from the reference only. Windows use the full protein's alignment; a subset is re-aligned as fragments to measure how much that matters.

**Outcome:** picks the score for phase 1's oracle and the seed shapes for phase 2's sweep.

*Done 2026-11-06 (`workflows/pfam-anchor/columns.py`); see Progress log. Plain identity stays the score, with the likelihood ratio to test as a tie-break. Conservation weights are dropped. Phase 2 keeps contiguous k 5, 6 and 8, spaced `1101011` and `11010010011`, and greedy columns in 12-column windows at k 3, 4 and 5.*

### Phase 2: column seeds

`anchor.py build` and `anchor.py query`, in Python and NumPy.

1. Index: per family, sorted `(seed, residues, unit)`. A plain Parquet table per family is enough here.
2. Sweep the seed shapes that phase 1b kept, with s of 1, 2 and 4, and a minimum of 1 or 2 shared seeds per unit. Compare each with brute-force AAI on time and memory, and find the family size above which seeds pay for themselves.
3. Add one reduced-alphabet arm (Murphy-10, from `tree-index/alphabets.py`), since low-identity seeds are the likely weak point.
4. Score candidates by column, with phase 1's score (identity over the read's columns).

**Gate:** an arm keeps at least 95% of the oracle's right calls, with median candidates per read at most 1% of the family's units for reads from known units, and at most 10% for reads from novel units. Phase 1b found no seed near 1% for novel units (.834 kept at best), so the first version of this gate could not pass. It only has to pass in families too large for brute force.

*Passed in large families 2026-11-08 (`workflows/pfam-anchor/anchor.py seeds`); see Progress log. Phase 3 takes `s1101011` at s 2 with 2 shared seeds, and contiguous k 5 at s 2 with 2 shared seeds as the alternative. The index is held in memory: no Parquet yet.*

### Phase 3: comparison

The chosen arm against DIAMOND and tree index tiers 1 and 2, on the same reads, with the same metrics and full Pfam-A read hits.

**Gate:** recall within 0.05 of DIAMOND's on reads with a hit, and at least 10× faster than DIAMOND for the unit-assignment step. Report the Pfam-search cost separately, since kfpfam is working on reducing it. Also report what share of DIAMOND's right calls have no Pfam hit, since these reads are beyond this method's reach.

### Phase 4: coverage and fallbacks

Each of these is measured separately, only if phases 1–3 pass:

- **Weak hits:** reference and read domains below GA, with `pfam-reannotate`'s reversed arm setting the bar. Do they add right calls beyond the wrong ones? On the read side, kfpfam's recruitment already gives each weak hit a reversed-read q-value and a status (see "Read hits from kfpfam").
- **Representatives or members:** representatives alone give 80% of k 6 postings but lose 12% of units with a domain. Do they lose right calls? This decides whether the all-representatives run is enough for the reference.
- **Extending past the domain:** the anchor gives a diagonal into the protein, so reads overhanging a domain could be extended through the sequence store, with sequences held outside the index as in the tree-index design.
- **Reads with no hit:** whether they go to tier 1 or DIAMOND, and at what cost.

### Phase 5: HPC scale

Nextflow workflow (Slurm, Singularity, `test` profile on the fixture), run by the user:

1. Build from the 1-in-100 `pfam-reannotate` output, to replace the release size estimate.
2. Run on `tree-index-slice` (4.2 M proteins), which has real near-unit density.
3. Measure the largest families, which set worst-case query cost.
4. Write an index format and port lookups to Rust only after this.

## Layout and scope

- `workflows/pfam-anchor/`: `prepare.py`, `anchor.py`, `README.md`; Nextflow from phase 5.
- `tests/python/test_pfam_anchor.py`: column mapping with insertions and deletions, k-mer keys, and a read whose source unit must come first on the fixture.
- **Not in kfp's CLI:** this is a development workflow. If it succeeds, its output is a candidate unit list in `kfp-pileup reference`'s format, as tree-index tier 2 hands off.
- **kfpfam provides the read-hit step.** The contract is the `pfam-reannotate` hits schema, with a query id, so a faster filter can replace `hmmsearch` without touching this work. Built in kfpfam on 2026-11-07: see the next section.

## Read hits from kfpfam

Added 2026-11-07. kfpfam (`~/Documents/kfpfam`) now writes read hits with their alignments, and recruits reads too weak for GA to families the sample confirms, keeping sets where a read is ambiguous. The write-up is kfpfam's `docs/recruit.md`, and the measurements are in its `docs/recruit-plan.md`. All of it ran on one gut metagenome (ERR10889147, 1.07 M reads) and on simulated samples.

**What it writes.**

- **The steps.** Six-frame ORFs come from `measure/real_orfs.py`, which also writes each ORF's place in the read. The k-mer filter's candidates (Mode A) go to pyhmmer in `recruit/pairs.py`. `recruit/recruit.py` then recruits, against reversed ORFs as the null.
- **`pairs.tsv`:** every (ORF, family) pair at E ≤ 1000, Z = 20 M, with the strongest domain's alignment. This is the recall-first input: no presence cut.
- **`hits.tsv`:** recruitment's result, one row per (locus, family in its set), with these columns:
  - `status` (`ga`, `one`, `set`, `none`) and `top`, which marks one call per locus;
  - `weight`, P(true) and q from reversed reads;
  - `read_from` and `read_to`, where the alignment lies on the read.

**Matching this plan's hits schema.**

- **Same columns.** `pfam_accession` (integer), `seq_evalue`, `score`, `c_evalue`, `i_evalue`, `hmm_from`, `hmm_to`, `env_from`, `env_to`, `ali_from`, `ali_to`, `cigar` and `above_ga`. The CIGAR uses `reannotate.py`'s convention (M, I residue only, D model only), so `columns.py`'s column mapping reads it unchanged. Checked on 50,000 reads' 16,362 hits: every hit maps consistently, and residues equal the HMM consensus at .48 (GA) and .38 (recruited) of their columns.
- **Differences.**
  - The query is `query`, a string ORF name (`R<read>_<n>`), not an integer `protein_id`.
  - Files are TSV, not Parquet.
  - There is one row per (ORF, family): the strongest domain, with no `domain` or `seq_score` columns.
  - `above_ga` tests the domain score only.
- **HMMs.** Pfam 38.2, from the same file as here (kfpfam's `Pfam-A.hmm.gz` links to `~/Downloads/Pfam-A.hmm.gz`), so HMM columns agree.
- **Thresholds.** Z = 20 M and E ≤ 1000, which is wider than this plan's `-E 10` at Z = 10⁶ (equal to E ≤ 200 at Z = 20 M). On shared pairs, pyhmmer's scores and coordinates are identical to `hmmsearch` on all ORFs. The filter in front loses .003 of GA pairs and .03–.11 of the weak bins.
- **ORFs.** At least 15 aa there, against 25 aa here. Filter on `orf_len` to match.

**What was measured that bears on this plan.**

- **More read hits.** On 1.07 M reads, there are 172,619 loci at GA (0.16 per read) and 108,750 recruited at an estimated 2% per-read FDR, of which 4,179 are sets. That is 0.26 hits per read in all. Phase 1b's brute-force AAI cost scales with hits per read, so its 1–9 core-ms per read becomes about 1.6–15 with recruited hits, plus a little for set members.
- **Weak hits with a bar already exist.** Phase 4's "Weak hits" item needs a reversed-read bar for read hits below GA, and recruitment supplies one per locus.
  - On a simulated sample of MGnify proteins at 100–50% identity, GA loci were .999 right. Recruited loci with one family were .960 right, and sets held the right family .944 of the time.
  - Recruited reads found more domains than GA calls alone from 80% identity down, and more than doubled them at 70% and below.
- **The reversed-read FDR is a lower bound.** On the simulations the true FDR ran 2–4 points above the estimate: .037 at an estimated .020. Two kinds of false hit come from real sequence, which reversed reads don't imitate:
  - ORFs in other frames of real coding DNA;
  - proteins hitting a family they aren't annotated with.

  They come in runs from one region. The unit-level score here is what should remove them, so report phase 3 and 4 precision by `status`.
- **Frames compete.** Placed on the read, 346 of 350 reads with GA calls on two ORFs had both over the same bases in different frames, never within one clan, so one frame was wrong. Only 4 were two genes. Take the `top` frame, or the set, rather than every frame's hit, or wrong frames will seed chance units.
- **Short ORFs rarely hit.** Windows of 10, 15, 20 and 25 aa cut from GA reads got any pair .04, .29, .62 and .85 of the time. In a contig simulation, only 9 of 180 domain-bearing segments at gene boundaries got a pair. This sizes the "Few or weak hits" risk below.
- **Cost of the read-hit step.** Filter plus pyhmmer comes to about 1.2 core-ms per read, and the reversed reads' null costs as much again. The per-read step needs only about 10,000 null reads (2% of 1.07 M gave the same counts within 1.5%), but confirming families from weak evidence needs the whole null. Recruitment itself takes 26 s per 1.07 M reads.

**What changes in the phases.**

- **Phase 0, step 3:** `recruit/pairs.py` in its `all:ACCS.txt` mode can search the reads against the 40 HMMs and write alignments directly, with no `hmmsearch` text parsing and no integer `orf_id`.
- **Phase 3:** take the Pfam-search cost from kfpfam's figures.
- **Phase 4, weak hits:** use `hits.tsv`'s `status`, `q` and `weight` as the weak-hit arm. Report right calls by status, and `pairs.tsv` against `hits.tsv` as recall-first against recruited.

## Risks

| Risk | Effect | Where it shows |
| --- | --- | --- |
| Fragment and full-protein alignments disagree | Wrong columns, lost seeds and identity | Phase 1b: 99% agree, without sequencing errors; phase 1 confirms |
| Reads outside domains | Low reach, whatever the k-mers | `hit` stage, by domain position |
| Few or weak hits on 30–50 aa ORFs | Reads never reach the index; weak hits add chance families | `hit` stage, reversed reads; kfpfam: a pair for .29 of 15 aa and .85 of 25 aa windows; recruitment adds 0.10 hits per read at 2% estimated FDR |
| Wrong-frame and unannotated-family read hits | Chance units seeded from real but wrong sequence; the reversed null misses them | kfpfam simulation: true FDR 2–4 points above the reversed-read estimate; frames over the same read bases compete in `hits.tsv` |
| Ties in conserved columns | Calls span many units | Median tied units; downstream EM or pileup |
| Family and clan cross-talk | Read and unit annotated to related families share no columns | Full Pfam-A arm; phase 3 misses |
| Very large families | Long postings at conserved columns | Phase 5 |
| Pfam version change | Columns change, so the index must be rebuilt | Pin 38.2 for now |
| Uniform simulated mutation | Overstates low-identity seeding or misstates it | Real relatives arm |

## Decisions for the user

1. **Gate thresholds.** The numbers in phases 1–3 are proposals.
2. **Where it lives.** This plan puts the work in this repo, as units are this project's concept, with kfpfam supplying read hits. It could go in kfpfam instead.

## Progress log

- **Test case sized (2026-11-06).** The tables under "What is measured so far". Scripts were scratch only, so phase 0 still has to write them into the repo.
- **Phase 1b: brute-force AAI, HMM weights and discriminating columns (2026-11-06, the user's request).** `workflows/pfam-anchor/columns.py`, on the 40 families' GA domains in the mini-slice, 80 s in one thread, 1 GB. 8,764 reads: windows of 40 HMM columns from real proteins, held out of the reference. *Known* (5,478): a non-reference member, truth its own unit. *Novel* (3,286): a unit held out whole, truth the remaining unit nearest over the whole domain.
  - **Limits of this test.** Mid-sized families only (20 k–150 k annotated proteins each; the largest Pfam families are far bigger). The mini-slice has fewer near units than the release. Windows have no sequencing errors and at least 75% occupancy. Each unit is one domain. Novel truth is itself identity over the domain, which favours plain identity as a score, so known reads are the fair comparison of scores.
  - **Brute-force AAI is feasible.** NumPy compares 2.7 GB/s on one thread, or 14.7 ms per 1 M units at 40 columns. Rough release cost, assuming read hits follow domain counts: STATS.md's top 25 families hold 438 M domains (PF00005 62.6 M), and 2.99 B proteins at about 1.75 GA domains each is about 5.2 B. A read hit then lands in a family with about 2–13 M domains, or roughly 0.5–3.6 M units at 3.5–4 domains per unit, costing 7–53 ms. At 0.17 GA hits per read (kfpfam, `ERR10889147`) that is about 1–9 core-ms per read, against 40–70 for `hmmsearch` and 29 for kfpfam's filter, so about 3–25 core-hours per 10 M reads. Worst case is about 250 ms for a hit in PF00005. Weak hits would add to this. SIMD should be several times faster (not measured).
  - **Memory favours the matrix.** One byte per unit and column is about 100 bytes per domain here, so about 60–100 GB for representatives at release (rough), mapped per family and laid out by column so a read touches only its own columns. A k 6 index at every column has about 105 postings per unit, which is larger at any realistic bytes per posting.
  - **Scores.** "Right" is truth among the tied best; "credit" is right divided by ties, so a score can't gain by tying more. Known reads at ≥ 90% identity, difference in credit from plain identity with standard errors over families:

    | Score | Right | Credit | Unique top | Credit difference |
    | --- | ---: | ---: | ---: | ---: |
    | Likelihood ratio, unit frequencies | .827 | .629 | .632 | +.019 ± .003 |
    | Likelihood ratio, HMM emissions | .827 | .628 | .630 | +.018 ± .003 |
    | Rarity, HMM emissions | .826 | .627 | .628 | +.017 ± .003 |
    | BLOSUM62 | .826 | .625 | .623 | +.015 ± .003 |
    | Plain identity | .831 | .610 | .593 | — |
    | Entropy weights | .819 | .610 | .613 | −.000 ± .001 |
    | Inverse conservation weights | .818 | .609 | .614 | −.001 ± .001 |
    | Conservation weights | .817 | .608 | .615 | −.002 ± .002 |

    On novel reads, identity has the highest right (.690, against .595–.646) and credit (.523), with HMM ratio level (−.000 ± .005) and conservation weights worst (−.055 ± .008), though the truth favours identity there. So conservation weighting does not help. Rarity and likelihood ratios help slightly, by breaking ties. Identity is more accurate when it gives a single best unit (.855 against .848), so a ratio as a tie-break after identity is the natural next test.
  - **Resolution per read is limited.** On known reads the true unit is among the tied best for only 83%, with about 5 ties on average. The called unit is close, though: for 93% it is within 2 points of the best whole-domain identity, on all reads of that scenario. Abundance from many reads has to resolve the rest.
  - **Fragment alignments.** The same windows as sequences, re-aligned with `hmmsearch`: 78% hit their source family at `-E 10`, Z = 10⁶ (no null yet), at about 36 columns. Of their residues, 99% are on the same column as in the full protein's alignment, and 88–95% of fragments agree completely, at every identity. Right on known reads is unchanged at .802, and on novel reads .690 becomes .680.
  - **Columns.** Across 6,243 columns with at least 30% occupancy, separating units and surviving divergence pull against each other: entropy against retention, Spearman −.77; entropy against HMM conservation −.68. Median entropy is 2.63 bits and median retention .80 between a unit and its nearest relative.

    | HMM conservation quintile | Conservation (bits) | Entropy across units (bits) | Retention | Entropy per bit of retention lost |
    | ---: | ---: | ---: | ---: | ---: |
    | 1 (least) | 0.19 | 3.16 | .708 | 6.4 |
    | 2 | 0.30 | 2.97 | .712 | 6.1 |
    | 3 | 0.41 | 2.65 | .757 | 6.8 |
    | 4 | 0.68 | 2.14 | .820 | 7.9 |
    | 5 (most) | 1.58 | 1.35 | .889 | 8.1 |

  - **Residues.** Value is a shared residue's rarity at its column, per bit of retention lost. Best: W 14.6 (retained .89), C 9.7 (.81), Y 9.6 (.86), F 8.9 (.83), M 8.1 (.75), P 7.9 (.85), G 7.6 (.85), H 7.6 (.78). Worst: A 5.8, E 5.7, S 5.7, I 5.9, V 5.9. These are averages over columns, and not yet used in any seed.
  - **Combinations.** Adjacent columns are the most redundant: mutual information over the smaller entropy falls steadily from .46 at distance 1 to .41 at 4, .33 at 11 and .23 at 30, with no helical periodicity at 3–4. This also reflects shared ancestry, so it isn't specific to column pairs. Spacing columns out gains about 0.5 bits per pair at most.
  - **Seeds.** Kept is the share of oracle-right reads with the true unit among candidates; share is median candidates as a share of the family; seeds per column sets index size.

    | Seed | Hits needed | Novel: kept | Novel: share | Known: kept | Known: share | Seeds per column |
    | --- | ---: | ---: | ---: | ---: | ---: | ---: |
    | Contiguous k 4 | 1 | .979 | 16.7% | 1.00 | 17.3% | .98 |
    | Contiguous k 5 | 1 | .940 | 7.3% | 1.00 | 7.9% | .97 |
    | Spaced `1101011` (k 5) | 1 | .959 | 7.8% | 1.00 | 8.2% | .95 |
    | Contiguous k 6 | 1 | .903 | 3.2% | 1.00 | 3.9% | .96 |
    | Spaced `11010010011` (k 6) | 1 | .908 | 2.7% | 1.00 | 3.5% | .92 |
    | Contiguous k 8 | 1 | .805 | 0.7% | .999 | 1.3% | .94 |
    | Contiguous k 8 | 2 | .770 | 0.4% | .996 | 1.0% | .94 |
    | Greedy, 12 columns, k 3 | 1 | .960 | 17.7% | .999 | 16.2% | .15 |
    | Greedy, 12 columns, k 4 | 1 | .915 | 6.3% | .996 | 6.5% | .15 |
    | Greedy, 12 columns, k 5 | 1 | .868 | 2.4% | .992 | 2.9% | .15 |
    | Most variable, 12 columns, k 4 | 1 | .827 | 1.1% | .998 | 1.7% | .51 |
    | Most conserved, 12 columns, k 4 | 1 | .965 | 29.0% | 1.00 | 29.3% | .51 |
    | Random, 12 columns, k 4 | 1 | .973 | 15.0% | 1.00 | 14.8% | .91 |

    - **Non-contiguous seeds:** spacing adds about 2 points of recall at the same k and a similar share (.940 to .959 at k 5). Random columns in a window do about as well as fixed spacing.
    - **Choosing columns:** this mainly shrinks the index, not the candidates. Greedy columns use 6–7× fewer seeds than contiguous k-mers, losing about 2–4 points of novel recall at similar shares. Conserved columns alone are poor at narrowing, and variable columns alone miss diverged relatives.
    - **Budgets:** known reads narrow to 1% with .997 kept (most variable, 12 columns, k 4, 2 hits). Novel reads keep at best .77 at 0.5%, .834 at 1%, .898 at 2%, .938 at 5%, .963 at 10% and .983 at 20%.
    - **Family size:** families of at least 1,500 units narrow about twice as far on share (k 6 novel 1.8% against 3.8%), with novel recall about 2 points lower. Whether that continues up to the largest families is unmeasured.
    - **Real fragment alignments:** results are similar: novel k 5 .958 and spaced `1101011` .976, at about 10% share.
  - **Reproducing.** Scratch steps until phase 0: cut the HMMs in `results-hpc/tree-index-slice/results/families.tsv`, run `hmmsearch` with `pfam-reannotate`'s thresholds on `data/mini-slice`, then `reannotate.py parse`, then `columns.py --width 40 --per-family 150`.
- **Read hits from kfpfam (2026-11-07, the user's request).** kfpfam's read-hit step and weak-read recruitment are built (kfpfam `docs/recruit.md`), and summarised in "Read hits from kfpfam" above. That section covers the schema mapping and the measurements that bear on phases 0, 3 and 4. Two of those measurements add risks, now in the table: the reversed-read FDR is 2–4 points optimistic, and wrong frames compete with real ones over the same read bases. Nothing in this plan's phases has run on kfpfam's output yet.
- **Phase 0: test case (2026-11-08, the user's request).** `workflows/pfam-anchor/prepare.py`, steps `annotate`, `reads`, `search`, `diamond` and `report` (`all` runs them in order), into `data/pfam-anchor/`. About 8 minutes in all, peak 2.1 GB (holdout's relative search).
  - **What it builds.**
    - `hits.parquet`: the 40 HMMs on all 294,822 mini-slice proteins, the same counts as the sizing.
    - `holdout/`: 1,686 real relatives held out of the reference (`recall.py holdout`, 250 per bin). The 85–70% bins only filled to 87–132.
    - `domains.parquet`: the reference's hits.
    - `strains/`: 2,000 strains of reference proteins with a GA domain, 100–70%, 11,953 pairs.
    - `relatives/`: the held-out proteins unchanged, 10,146 pairs.
    - Each read set gets `truth.parquet` (source protein; truth unit, the nearest remaining relative's for relatives), ORFs of at least 25 aa with read, mate, strand and start, hits for the reads and for the reads reversed, and DIAMOND's best hit per mate.
    - `spike.py simulate`'s `reads.parquet` now also records each pair's fragment (start from the CDS, length, strand), so a pair can be placed on its source.
  - **Hit recall.** Pairs with a mate overlapping a GA domain of the source by at least 25 residues; *any* is any hit at `-E 10`, Z = 10⁶, on a family of the source:

    | Identity bin | Strains: GA | Strains: any | Relatives: GA | Relatives: any |
    | ---: | ---: | ---: | ---: | ---: |
    | 1.00 | .517 | .805 | .572 | .850 |
    | 0.98 | .550 | .842 | .556 | .841 |
    | 0.95 | .492 | .754 | .511 | .811 |
    | 0.90 | .410 | .707 | .551 | .805 |
    | 0.85 | .407 | .726 | .651 | .866 |
    | 0.80 | .295 | .603 | .517 | .857 |
    | 0.75 | .189 | .538 | .551 | .812 |
    | 0.70 | .172 | .488 | .531 | .817 |
    | 0.60 | | | .491 | .790 |

    - **Relatives keep their hit rate.** Their reads are the held-out protein itself, so the bin is the distance to the reference, not to the HMM. Strains lose hits with identity, as uniform mutation erodes the HMM match.
    - **Reach is low over a whole protein.** Over all pairs, `any` is .32–.58 for strains and .26–.57 for relatives. Pairs outside a domain rarely hit (under .15), and only 29–66% of pairs lie in one. Relatives' sources are drawn from all members, so some have no domain from the 40 families.
    - **The null is quiet.** Reversed reads hit any family for .001 of pairs and none at GA.
  - **DIAMOND.** Best hit per pair against the reference (`blastx --sensitive -k 1`): .98–1.00 of pairs have one, and their median identity follows the bin (strains 100 to 72%, relatives 100 to 77%). This is the method-independent identity for phase 1.
  - **Crash and fix.** The first `reads` run crashed the laptop: `recall.py holdout`'s relative search joined every pool k-mer with all 295 k members at once. It now works through the pool in chunks of 500, and each `prepare.py` step runs under a watchdog that kills its process group above `--max-gb` (6 GB by default). DIAMOND's container is capped the same way, since Colima's VM can take 20 GB.
- **Phase 1: oracle, no k-mers (2026-11-08, the user's request).** `workflows/pfam-anchor/anchor.py oracle`, on phase 0's case: 16 s, 2.9 GB peak.
  - **Method.** Every hit of a read pair, at `-E 10` and in any family, is compared with every reference GA domain of that family. Per hit, each unit takes its best domain, which must share at least 60% of the hit's columns. Per pair and family, matches are summed over the hits, and the units with the best identity are the call, with ties kept.
  - **DIAMOND.** `prepare.py diamond` now keeps every hit within 10% of each mate's best score (`--top 10`), with coordinates and CIGAR. Per unit, DIAMOND's score is the best bitscore summed over the mates, with ties kept. *Restricted* keeps only units with a GA domain in a family the pair hits, which is the oracle's universe. It differs from all units by at most .013.
  - **The score.** Identity over shared columns, phase 1b's choice, failed the gate on strains, 6–15 points below DIAMOND. It favours units whose domain covers only part of the read. For example, a read at 100% scored 46/48 and 47/48 on its own unit, but 48/48 and 38/38 on a partial domain. Identity over all the read's columns fixes this. A unit's missing columns then count against it, much as alignment length raises DIAMOND's bitscore.
  - **Results.** Pairs with a hit on a family of the source:

    | Identity bin | Strains: oracle right | DIAMOND | Oracle credit | DIAMOND credit | Relatives: oracle right | DIAMOND | Oracle credit | DIAMOND credit |
    | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
    | 1.00 | .950 | .990 | .715 | .843 | .961 | .956 | .714 | .776 |
    | 0.98 | .941 | .968 | .738 | .828 | .974 | .964 | .699 | .789 |
    | 0.95 | .948 | .956 | .708 | .829 | .856 | .796 | .656 | .709 |
    | 0.90 | .942 | .929 | .671 | .789 | .696 | .654 | .520 | .578 |
    | 0.85 | .946 | .935 | .687 | .814 | .634 | .591 | .520 | .556 |
    | 0.80 | .946 | .912 | .676 | .803 | .487 | .487 | .407 | .463 |
    | 0.75 | .931 | .912 | .668 | .791 | .485 | .447 | .427 | .424 |
    | 0.70 | .927 | .894 | .707 | .818 | .310 | .337 | .238 | .330 |
    | 0.60 | | | | | .351 | .312 | .286 | .294 |

    - **The gate passes on `right`.** Strains are at most .040 below DIAMOND (at 100%). Relatives are level with it or above, at every bin.
    - **The oracle ties more.** It averages 2–5 tied units against DIAMOND's 1.1–2.1, since it sees only the domain's columns, not the rest of the read. So its credit is .09–.13 lower on strains and up to .09 lower on relatives. Phase 1b's likelihood-ratio tie-break, or extending past the domain (phase 4), are the ways to recover it.
    - **Candidates are almost complete.** The truth unit qualifies on .96–1.00 of pairs. The pairs lost are those where its domain doesn't cover the read's columns.
    - **Real relatives are the hard case.** Right falls from .96 at 100% to .49 at 80% and .31–.35 below that, for DIAMOND as much as for the oracle. A single read pair rarely separates a novel unit's nearest relative from its neighbours. Reads at GA alone change little (`right_ga` in `oracle.tsv`).
  - **Alignment agreement (step 3).** For 10,427 ORF hits that DIAMOND also aligns to the truth protein, the residue pairs the two HMM alignments imply match DIAMOND's for 97.6–99.2% of residues. In each bin, 80–91% of hits agree completely, with no trend in identity, so the shared template loses little.
  - **Other findings.**
    - Short repeat families, such as PF12661 (22 columns), lose the most. Reads carry several overlapping domain hits, some of only 2–7 columns with scores below zero, and the copies of the repeat in the reference are often below GA.
    - With the shared-column score, 48 of 633 misses at 80% and above (8%) fell only on a reference domain below GA. Phase 4's weak-domain arm covers them.
- **Phase 2: column seeds (2026-11-08, the user's request).** `workflows/pfam-anchor/anchor.py seeds`, on the pairs with a hit on a family of their source: 30 s, 1.8 GB peak.
  - **Method.**
    - **Index:** per family and shape, sorted postings of (seed, residues, unit), merged by unit and held in memory. Smaller s keeps only seeds whose first column is a multiple of s.
    - **Query:** a pair uses every seed of all its hits in the family. A unit is a candidate with at least h distinct matching seeds.
    - **Scoring:** candidates are scored as in the oracle, so *kept* is the share of the oracle's right pairs that stay right, and *share* is candidates over the family's units (median per pair).
    - **Sets:** strains are the known-unit reads, relatives the novel-unit ones. *Large* families have at least 1,500 units: 9 of the 40, the largest 3,090.
  - **Results.** Selected arms; all 48 are in `seeds.tsv`.

    | Arm | s | h | Postings per unit | Strains: kept | Strains: share | Large: kept | Large: share | Relatives: kept | Relatives: share | Large: kept | Large: share |
    | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
    | k 5 | 1 | 1 | 105 | .996 | 8.0% | .997 | 6.0% | .996 | 11.8% | .995 | 8.8% |
    | k 5 | 2 | 2 | 53 | .971 | 2.0% | .965 | 0.9% | .979 | 2.6% | .981 | 1.5% |
    | k 6 | 1 | 2 | 105 | .968 | 2.1% | .960 | 1.0% | .980 | 2.9% | .982 | 1.6% |
    | k 6 | 4 | 1 | 26 | .947 | 1.9% | .940 | 0.9% | .971 | 2.7% | .969 | 1.6% |
    | k 8 | 1 | 1 | 103 | .925 | 1.0% | .916 | 0.4% | .962 | 1.5% | .961 | 0.7% |
    | `1101011` | 2 | 2 | 51 | .965 | 1.4% | .958 | 0.6% | .978 | 1.9% | .977 | 1.0% |
    | `11010010011` | 1 | 2 | 97 | .965 | 1.1% | .964 | 0.4% | .969 | 1.5% | .967 | 0.7% |
    | Greedy k 4 | 1 | 1 | 15 | .933 | 6.1% | .930 | 5.1% | .983 | 9.7% | .986 | 8.6% |
    | Murphy-10 k 8 | 2 | 2 | 49 | .903 | 1.4% | .883 | 0.4% | .968 | 2.3% | .967 | 1.0% |

    - **The gate passes in large families, not over all of them.** In large families, contiguous k 5 or k 6, and both spaced shapes, keep at least .95 at no more than 1% for strains, with relatives far inside their 10%. Over all families, the best strains arm is just above 1% (`11010010011`, s 1, h 2: .965 at 1.1%). The gate only has to hold where brute force is too slow, so it passes.
    - **Novel reads narrow as well as known ones.** Relatives keep more than strains at a similar share. This differs from phase 1b, where novel windows were the hard case. Here, relatives' reads are the held-out protein itself, so they have no sequencing divergence, while strains are mutated to 70%. Strains at 70–80% are the weak point.
    - **Two shared seeds beat a longer k.** At half the postings, k 5 or `1101011` at s 2 and h 2 keep 4–5 points more on strains than k 8 at h 1 (.965–.971 against .925), for a little more share (1.4–2.0% against 1.0%).
    - **Greedy columns and Murphy-10 don't help.** Greedy seeds are 7× smaller, but at .93 kept on strains they narrow only to 6%. Murphy-10 is worse than the full alphabet at the same k, s and h.
  - **Speed.** Brute force costs about 63 ns per reference domain per hit, so 0.25 ms per hit in the median family here. A seed lookup costs about 0.02 ms per hit, mostly Python overhead.
    - Scoring the candidates costs brute force times their share, so seeds pay off above about 330 domain rows (about 100–150 units).
    - With the chosen arms, unit assignment is 7.5–8× faster than brute force overall, and 12–14× faster in large families. The gain grows with family size, since the share falls (phase 1b). At release-sized families of millions of units, brute force is about 60 ms per million domains per hit, against a lookup plus 1% of that.
  - **Limits.** Families here have at most 3,090 units, so the share at release sizes is extrapolated. The postings are per unit, not per protein. Reads have the simulator's 0.2% sequencing errors, but no real ORF calling.
- **Anchoring ablation (2026-11-08, the user's question).** `anchor.py seeds --unanchored`: the same seeds, read off the same HMM alignments, but keyed by residues only, so a seed matches at any column of its family. Recall (kept) is unchanged to .001 for every arm. Candidate share is unchanged for k ≥ 6 and spaced `11010010011`, and for every arm at h 2. It grows slightly for short seeds at h 1: k 5 s 1, 8.0% to 8.8%; `1101011` s 1, 9.2% to 10.1%. Postings fall by under 0.2%. So within a family, the column in the key adds almost nothing. A 5-residue seed rarely recurs at another column of the same family, so the narrowing comes from restricting to the read's family, and from reading seeds off match columns (insertions dropped). Still unmeasured: a pooled index across families, which would show what the family restriction itself is worth, and short seeds (k ≤ 4) or the largest families, where the column may matter more.
