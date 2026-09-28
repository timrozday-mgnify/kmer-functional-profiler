# Protein k-mer functional profiler — assessment & implementation plan

Sep 26, 2026 · @Tim

## Verdict

The core idea has been done: fmh-funprofiler already does FracMinHash functional profiling on translated amino-acid k-mers against KEGG orthologs, and kMermaid and UProC are older/parallel k-mer approaches. What is not done is your specific combination: a guaranteed per-gene sampling floor, frame-aware read translation, and a proper abundance model. That is a real but narrow niche, and it targets fmh-funprofiler's documented weakness (poor completeness for low-abundance functions).

Keep, change, drop:

- **Keep: translation to amino acids.** Sound and established. Silent-site tolerance is real; the cost is weaker specificity per residue, so k must rise.
- **Keep, but downgrade: frame detection.** Mostly a speed and noise reduction (about 6 frames to \~1.5), not a sensitivity gain. A stop-codon filter gets most of it for free.
- **Keep, but reformulate: length-scaled density.** Do it as a per-gene (or per-function) FracMinHash threshold with nested thresholds. The post-hoc correction is then exact, not approximate.
- **Drop: start/stop anchoring.** It either breaks sketch consistency (position-based selection that reads cannot reproduce) or adds almost nothing (content-based anchors hit few reads). The per-gene floor already delivers what anchoring was meant to.
- **Drop: stacking minimizers + hash threshold + anchors.** Use one context-free sampler. Three layered samplers make the statistics hard and the reads-vs-reference consistency fragile.

The harder problems are not sparsity. They are sensitivity to divergent homologs (exact amino-acid k-mers fail below roughly 80% identity), k-mers shared across functions (domains, paralogs), and database licensing (KEGG is not freely redistributable). The plan below addresses those first.

## Prior art

The closest tool is fmh-funprofiler; everything in your plan except the density control and frame handling has a precedent.

| Tool | What it does | Relevance to your plan |
| --- | --- | --- |
| [fmh-funprofiler](https://github.com/KoslickiLab/fmh-funprofiler) (Koslicki lab, Bioinformatics 2024) | Sketches each KEGG ortholog group and the translated metagenome with sourmash protein FracMinHash, then uses `sourmash prefetch` overlaps as abundances. Default scaled = 1000. | Direct predecessor. Reported 39–99x faster and up to 40–55x less memory than DIAMOND, with comparable completeness and better purity. Its low completeness on rare KOs comes from downsampling, which is exactly your target. |
| [kMermaid](https://pmc.ncbi.nlm.nih.gov/articles/PMC12507277/) (PLoS Comp Biol 2025) | Assigns reads to taxa-agnostic clusters of homologous proteins using amino-acid k-mer frequency profiles in a nested hash map. | Unsparsified k-mer approach; a speed/accuracy baseline and an example of clustering proteins before indexing to reduce multimapping. |
| [UProC](https://academic.oup.com/bioinformatics/article/31/9/1382/200454) (2015) | "Mosaic matching" of amino-acid words against Pfam or KEGG families, including ORF selection from DNA reads. | Early proof that word-based protein classification beats profile HMMs on 100 bp reads. |
| [Orpheum](https://github.com/czbiohub-sf/orpheum) (CZ Biohub) | Picks the coding frame of each read by six-frame translation and the highest reduced-alphabet k-mer Jaccard against a reference proteome. | Your frame-detection idea, done with database containment rather than codon statistics. Evaluated on bacteria in the [metapangenome paper](https://dib-lab.github.io/2021-paper-metapangenomes/). |
| [Metabuli](https://github.com/steineggerlab/Metabuli) (Nat Methods 2024) and its 2026 update | "Metamers": 8-codon k-mers storing both the amino-acid translation and codon identity; exact amino-acid matching, DNA Hamming distance for specificity. The update adds spaced metamers, reduced alphabets and syncmers. | Taxonomic, not functional, but the best template for your k-mer encoding and sampling. Syncmers halved database size and doubled speed at a small recall cost. |
| [FragGeneScanRs](https://github.com/unipept/FragGeneScanRs) | Rust HMM gene prediction for short, error-prone reads. GPL-3.0. | Option for frame calling, especially with indel errors. The licence matters if you link it. |
| sylph | FracMinHash containment with a statistical model for low coverage (taxonomic). | Template for coverage-aware abundance and containment correction. |
| HUMAnN 3/4, DIAMOND, MMseqs2 | Translated alignment. | Accuracy baselines; HUMAnN is the de facto standard for gene-family profiles. |

Nothing I found anchors k-mer selection to start/stop codons or scales density by gene length. There may be a reason: see the next section.

## Critique of each element

### Sketch consistency is the constraint everything must respect

A sketch works only if a read selects the same k-mers the reference selected, using information the read has. Any rule that depends on position within a gene fails this: a 150 bp read usually does not know where the gene starts. Selection must be a function of the k-mer itself (FracMinHash, syncmers) or of a local window the read also contains (minimizers, with edge losses at read ends).

### Start/stop anchoring

- **Position-based version** ("keep k-mers within x residues of the start") is not reproducible from reads. Broken.
- **Content-based version** ("always keep k-mers containing a start/stop") is consistent, but weak. Internal Met is \~2% of residues, so "contains M" is just a biased subsample. Stops appear only if you append `*` to reference proteins; then only reads crossing the stop produce that anchor, a few percent of reads per gene.
- True starts are ambiguous (GTG/TTG, annotation errors), and termini are usually less conserved than core domains.
- In protein space, the reference is already a set of proteins with defined ends. ORF finding on the reference is unnecessary; it only matters on the read side.

### Translation and frame detection

Translation is the right call. The cost: an exact amino-acid k-mer of length 11 survives at identity *p* with probability *p*^11, so 0.31 at 90%, 0.09 at 80%, 0.02 at 70% identity. Functional orthologs are often below 70%. Reduced alphabets and spaced seeds (Metabuli's route) are how you recover sensitivity.

Frame detection by stop codons is cheap and effective at moderate GC. A random 50-codon frame has no stop with probability (61/64)^50 ≈ 0.09, so a coding read keeps \~1.5 of 6 frames. At 70% GC, stop codons are rarer (\~1.9% per codon), about 0.38 of off-frames survive and you keep \~2.9 frames. Frame filtering is a speed and false-positive optimisation, not a sensitivity gain: off-frame k-mers rarely match a protein database anyway.

### Minimizers vs hash threshold

Random minimizers do nest: with a fixed order, a (w=20)-minimizer is also a (w=10)-minimizer, since it is the minimum of every sub-window containing it. But they depend on neighbouring k-mers, so short reads lose picks at their edges, and the density is not an unbiased fraction. FracMinHash is context-free and gives a clean binomial model. Closed syncmers are the middle ground: context-free with a spacing guarantee (Edgar 2021). Use one sampler, not three.

### Length-scaled density and post-hoc correction

The problem is real. At scaled = 1000, a single 300-aa gene has \~0.29 expected sketched k-mers, so it is detected with probability \~25%. fmh-funprofiler hides this by pooling all genes of a KO, and pays in completeness for rare KOs.

The fix works if thresholds are nested. Give each unit *g* (protein cluster or function) its own threshold *t\_g* = max(*t\_base*, *n\_min* / (*L\_g* − *k* + 1)). The query keeps every k-mer with hash < max(*t\_g*). A hit counts for *g* only if hash < *t\_g*. Because {hash < *t\_g*} ⊆ {hash < *t\_max*}, per-unit containment and coverage are unbiased, with *t\_g* as the known sampling rate. No approximation needed.

The catch: the query side runs at the densest rate any unit needs. So the saving is in index size and lookups, not hashing, which you do for every k-mer regardless (SIMD hashing makes that cheap). Stream the query and count hits; do not store a query sketch.

### What the plan misses

- **Shared k-mers.** Multi-domain proteins and paralogs share k-mers across functions. Summing hits double-counts. Needs EM over read-level equivalence classes (salmon/kallisto style) or a gather-style greedy assignment.
- **Normalisation.** Hits per sketch k-mer ≈ k-mer coverage; converting to copies per genome needs a single-copy marker or average-genome-size normaliser.
- **Database.** KEGG requires a licence for bulk data; eggNOG, Pfam, UniRef and KOfam HMM-derived annotations are open alternatives.

## Revised design

One context-free sampler with nested per-unit thresholds, frame-filtered translated reads, and an EM abundance model.

**Reference (index build)**

1. Input: MGnify Proteins 2026\_07 FASTA plus its cluster membership: 90% clusters (MGnify90) as units, Pfam as the functional label (see Database below). No ORF calling.
2. Optional reduced alphabet (e.g. Murphy-10, Dayhoff) and optional spaced-seed mask; amino-acid k of 8–12, packed into u64.
3. Per 90% cluster *g*: *t\_g* = min(1, max(*t\_base*, *n\_min* / *n\_kmers(g)*)) for non-singleton clusters, *t\_g* = *t\_base* for singletons. Keep k-mers with hash < *t\_g*.
4. Store hash → list of (unit id); store *t\_g* and *n\_kmers(g)* per unit. Sorted array or minimal perfect hash, memory-mapped.

**Query (streaming, no stored sketch)**

1. Read FASTQ; for each read and mate, translate the frames that survive the stop filter (default) or all six (`--frames all`); split at stops; drop segments shorter than k.
2. Hash every k-mer; skip if hash ≥ *t\_max*; look up survivors.
3. Per read (pair), collect the set of units hit, counting a hit for *g* only if hash < *t\_g*. Emit equivalence class counts.

**Abundance**

- Per unit, hits ≈ Poisson(λ\_g × *m\_g*), where *m\_g* is the sketch size and λ\_g the k-mer coverage.
- Two-stage assignment: uniqueness-weighted detection, then EM per connected component of clusters that share kept k-mers (see Assignment below).
- Report: coverage λ\_g, containment (fraction of sketch k-mers seen) with a sylph-style low-coverage correction, a copies-per-genome value using a single-copy marker set, and relative abundance.

**Knobs to ablate:** *t\_base*, *n\_min*, k, alphabet, frame mode, sampler (FracMinHash vs closed syncmers), unit granularity.

A `--sourmash-compat` mode using sourmash's protein encoding and hash lets you check your containment numbers against fmh-funprofiler exactly before diverging.

## Database: MGnify Proteins

Use MGnify's own 90% clusters as units and index the k-mers of all members, not just representatives; skip your own 95% dereplication. MGnify90 clusters are the only grouping level: the floor applies per cluster but only to non-singleton clusters, or the index will not fit in memory.

**What the release provides.** The latest release is 2026\_07: 5.74 billion sequences in 1.66 billion 90% clusters (73% singletons), CC0-licensed, so pre-built indexes can be redistributed. Pfam hits are provided for all proteins (52% have one). There is no 30% membership table (MGnify30 ships only as representative FASTAs of three subsets), so this plan uses MGnify90 clusters only; grouping them into 30% families is future work. Details, file sizes and Parquet layout: [mgnify-2026\_07.md](mgnify-2026_07.md).

**95% vs 90%.**

- No 95% clustering is published. Making one means clustering 5.7 billion sequences yourself.
- The identity threshold only changes k-mer content if you index representatives alone. A member at 90% identity shares only \~0.9^11 ≈ 31% of its 11-mers with the representative, so a rep-only index loses most k-mers of divergent members.
- If member k-mers are indexed and labelled with their cluster, the threshold only sets unit granularity. 90% is already finer than any functional question needs; 95% would add units and shared k-mers with no gain.
- If memory forces a cut, subsample member k-mers (lower *t\_base*), not the identity threshold.

**Index size (rough, to verify on real data).**

| Component | Estimate | Assumption |
| --- | --- | --- |
| Base sketch at *t\_base* = 1/1000 | \~5×10^8 hashes, \~6–10 GB | \~5.7×10^9 proteins × \~200 aa, \~half distinct |
| Floor at 90%-cluster level, *n\_min* = 8 | ≥ 1.3×10^10 hashes, \~160 GB | 1.66×10^9 clusters in 2026\_07 (0.45×10^9 non-singletons) |
| Floor at 90% non-singleton clusters only, *n\_min* = 8 | ≥ 3.6×10^9 hashes, \~45 GB | 0.45×10^9 non-singleton clusters |

So: apply the floor to non-singleton 90% clusters only, and give singletons the base rate. If that is still too large, lower *n\_min* or restrict the floor to clusters with a full-length member.

**Functional labels.** Pfam hits are published for every protein, so each cluster can take the Pfam labels of its members. KO or eggNOG labels would need your own annotation run on cluster representatives (0.45×10^9 non-singletons: expensive).

**Practical options.**

- Development index from a biome subset (e.g. human gut) using MGnify's biome counts per sequence: small enough for a laptop.
- A two-stage screen like sylph's new `.syl2db` format: a very sparse first pass picks candidate clusters, a denser second pass quantifies only those ([sylph releases](https://github.com/bluenote-1577/sylph/releases)).

## Assignment of shared k-mers

Uniqueness-first is a good rule for deciding which units are present, but a poor one for quantifying them; use it as the detection stage and EM for abundances.

| Method | Rule | Weakness here |
| --- | --- | --- |
| sourmash gather | Greedy: take the target with the largest remaining overlap, remove its k-mers, repeat. | Order-dependent; favours big units. |
| sylph profile | Winner-take-all: each shared k-mer goes entirely to the genome with the highest containment ANI, then ANI is recomputed on the reassigned k-mers ([sylph paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC12339375/)). | Works at species level (>95% ANI, few close neighbours). Dense 90% protein clusters have many near-equal competitors. |
| Uniqueness-first greedy (your idea) | Rank units by unique or least-shared k-mer support; assign those first. | See below. |
| EM on equivalence classes | Shared hits split in proportion to current abundance estimates; unique hits pin those estimates. | More code; needs candidate pruning at MGnify scale. |

**Why pure uniqueness-first greedy is risky for 90% clusters**

- Unique k-mers are not a random sample of a unit. They sit at the variable positions, which is exactly where a sample strain tends to differ from every reference. Abundance from them is biased low for divergent strains.
- Units in dense neighbourhoods (many near-identical clusters) have few unique k-mers, so estimates are noisy just where competition is highest.
- A sample protein sitting between two references (93% to each) should split roughly evenly; any greedy rule gives it all to one.
- With per-cluster thresholds, raw unique counts must be normalised by each unit's sampling rate *t\_g* before ranking.

**Recommended two stages**

1. **Detection (your idea, formalised).** At build time, store per unit its sketch size *m\_g* and unique sketch k-mers *u\_g* (kept k-mers found in no other cluster). At query time, test unique hits against a Poisson background of false hits, and rank by an IDF-weighted score Σ hits(x) / |units(x)|. Keep units that pass.
2. **Quantification.** Within each connected component, EM over the surviving units: hits(x) \~ Poisson(Σ\_{g ∋ x} λ\_g), zero-inflated for divergent k-mers as in sylph. Components (clusters linked by shared kept k-mers, computed at build time after dropping promiscuous k-mers) are independent, so this runs in parallel.

At the sketch densities planned, most reads carry zero or one hit, so read-level equivalence classes collapse to k-mer-level classes. Run the EM on k-mer hit counts; it is simpler and loses nothing.

Pfam totals do not depend on how hits are split between clusters that carry the same Pfam label, so most errors in the within-component split do not reach the level most users will report.

**Evaluation:** all four rules run on the same hit table, so implement gather, winner-take-all and uniqueness-first as cheap baselines in phase 4 and compare against EM in phase 5.

## Homologous groups and k-mer weighting

Treat each 90% group as one unit whose k-mer set is the union of its members, and store how often each k-mer occurs inside the group. That in-group frequency, set against the out-of-group frequency, is both the selection criterion and the likelihood weight.

**How fmh-funprofiler does it.** It builds one FracMinHash sketch per KEGG ortholog from the amino-acid k-mers of all its genes (k = 7, 11, 15; scaled = 1000; abundances tracked). The authors tried per-gene sketches and dropped them because genes are too short to keep enough k-mers. They run `sourmash prefetch`, which reports each KO's overlap independently, and derive abundances from those overlaps ([paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC11373326/)). So there is no within-group weighting and no reassignment of k-mers shared between KOs: a shared k-mer counts for every KO containing it.

**In-group vs out-of-group frequency.** For k-mer *x* and group *G*:

- *p\_in(x)*: fraction of G's members containing *x*. High means core: a read from an unseen member of G probably carries it. Private k-mers (one member) are weak evidence.
- *p\_out(x)*: number of other groups containing *x*. Low means specific.
- Score: log(*p\_in* / *p\_out*), a naive-Bayes log-likelihood ratio. Precedents: CLARK keeps only target-discriminative k-mers; MetaPhlAn picks markers by in-clade core plus out-of-clade uniqueness.
- Compute *p\_in* over full-length members only, with sequence weights (as in profile HMMs) so over-sampled lineages and partial predictions do not skew it. Quantise it to 2–4 bits per posting.

**Score-based selection is free on the query side.** The query checks every k-mer with hash < *t\_max*. So the index can keep any score-chosen subset of those k-mers, and queries stay consistent. Rules:

1. Among k-mers passing the hash threshold, keep the highest-scoring ones per group until the floor *n\_min* is met.
2. Drop k-mers present in more than N clusters outright (low-complexity runs, common motifs such as Walker A). This also removes the longest posting lists.
3. Record which k-mers were kept, so each group's expected hit count uses the kept set, not a nominal rate.

**Likelihood with *p\_in*.** hits(*x*) \~ Poisson(λ\_G · *p\_in(x)*) turns strain variation into part of the model rather than zero-inflation noise. It also yields a per-group divergence estimate: a containment-AAI analogous to sylph's containment ANI. That is a useful novelty signal per function.

**Tiered lookup ("test the most unique first").**

1. A small tier-1 index of high-score k-mers (high *p\_in*, low *p\_out*), a few per group, stays in RAM. Every query k-mer probes it.
2. Groups with tier-1 support above background become candidates.
3. Only candidates' full k-mer sets (tier 2, memory-mapped) are probed, via a second pass over the reads or buffered query hashes.

This is the same pattern as sylph's two-stage database. At MGnify scale it is probably what makes the index usable on a normal node.

## Index encoding

Store \~12–16 bits per hash, with the high bits implicit in the layout, and spend the real effort on compressing unit IDs, which will dominate the index.

**Fingerprints.** A bare 16-bit key fails: \~5×10^8 hashes over 2^16 values puts thousands of k-mers on every value, and no model recovers the signal. What sets the false-hit rate is the effective hash length:

- Partition by the high bits (e.g. 2^24 buckets addressed directly) or use Elias–Fano over sorted hashes; entries store only the low bits.
- Hashes below the FracMinHash threshold share log₂(1/*t*) leading zero bits for free.
- Stored bits per entry ≈ log₂(1/ε) + 2, for false-hit rate ε per lookup.

| Stored bits per entry | False-hit rate ε per lookup |
| --- | --- |
| 8 | \~1.6×10^-2 |
| 12 | \~10^-3 |
| 16 | \~6×10^-5 |

**Statistics.** A spurious hit lands on a random entry, so group *g* gets an expected Q·ε·*m\_g*/n false hits, where Q is the number of query lookups, *m\_g* the group's kept k-mers and n the index size. Add this as a known Poisson background in detection and EM, and verify it with shuffled or decoy reads.

**The floor that matters more.** Genuine chance matches already occur. With 20 amino acids and k = 11 there are \~2^47.5 possible k-mers; MGnify might hold \~2^39 distinct ones (rough), so \~0.3% of off-target query k-mers hit something by coincidence. A 16-bit fingerprint adds far less. Reduced alphabets shrink the space (Murphy-10 at k = 11 gives \~2^36.5, near saturation), so they need a larger k.

**Unit IDs.** With \~10^9 units an ID costs \~30 bits, and conserved k-mers carry many. Map each k-mer to the ID of its distinct unit set and store each set once, as kallisto, Themisto and Fulgor do; sets repeat heavily. Store component IDs in tier 1 and resolve units only in tier 2.

**Rules.** Compare thresholds on the full 64-bit hash before truncating, so nested per-cluster thresholds stay exact. Keep the fingerprint width a build parameter and measure ε against index size.

## Error vs variation model

Adopt sylph's separation of sequencing error from divergence, with two changes: a mixture-aware likelihood, and a dense second tier so each group has enough k-mers to fit.

**What transfers.** Errors thin every k-mer's count roughly equally and lower the effective coverage λ. Divergence removes a k-mer from all reads and shows up as a structural zero. Estimating λ from the counts of present k-mers gives:

- abundance that stays correct when the sample variant is only \~90% identical to the reference (total hits divided by sketch size understates it);
- coverage-corrected containment, and from it a containment-AAI per group.

In protein space errors are gentler than in DNA: \~25% of random substitutions are synonymous, and a stop-creating error just splits the read. For Illumina this all folds into λ and cancels in relative abundances.

**Model.** For group *G* and kept k-mer *x*:

- *z\_x* \~ Bernoulli(*π\_x*), with *π\_x* informed by *p\_in(x)*;
- count(*x*) | *z\_x* = 1 \~ NegBin(mean λ\_G · *r*, dispersion φ\_G), where *r* is the error-thinning factor; count = 0 if *z\_x* = 0;
- empirical-Bayes priors on λ\_G and on divergence, shared across groups, to stabilise sparse groups.

**Mixtures.** Groups are often covered by several strains at different coverages, so core k-mers get summed coverage and variant k-mers partial coverage. The negative binomial absorbs mild mixing. A dispersion test per group flags groups whose counts don't fit a single Poisson; the flag is reported as "multi-variant". A small latent-variant mixture is a later option if the flag turns out to be common.

**Enough k-mers.** A sketched group has tens of k-mers, too few below \~1x where nearly all counts are 1. Tier 2 is therefore dense (or near-dense) for candidate groups only, giving hundreds of k-mers per fit at bounded cost. Report confidence intervals, not point estimates, at low coverage.

**Absolute coverage (optional).** Estimate *r* from read quality scores, or from near-invariant single-copy marker groups where divergence ≈ 0.

**Out of scope for now.** Long-read indels cause frameshifts, which remove every downstream k-mer in that frame. That is not uniform thinning and needs a frameshift-aware model.

## Development strategy: Python first, Rust at the end

Prototype the algorithm in Python, with stable hot loops in Rust from day one via PyO3, then port the rest to Rust once evaluation has frozen the design. The Python version stays in the repo as the reference oracle for the Rust port.

| Component | Language during development | Why |
| --- | --- | --- |
| FASTQ streaming, translation, stop-filter frames, amino-acid k-mer packing, hashing, threshold filter | Rust from phase 1 (PyO3 module) | Hot, well-specified, unlikely to change. Pure Python would make every experiment on real data too slow. |
| Index lookup (sorted hashes, fingerprints) | Rust once the layout settles (end of phase 3) | Memory and speed matter on MGnify-scale subsets. |
| Index build: clustering tables, p\_in/p\_out scoring, floors, tiers | Python (polars, DuckDB) | Where most design changes happen. |
| Assignment, EM, zero-inflated negative-binomial model, dispersion test | Python (numpy, scipy) | Statistical design still open; easiest to iterate and inspect. |
| Evaluation, plots, ablations | Python + Nextflow | Stays Python permanently. |
| Final CLI and full pipeline | Rust (phase 6) | Release product. |

**Rules for the split**

- The Rust kernels expose batch APIs that return numpy arrays (e.g. read index, frame, hash per kept k-mer), never per-k-mer Python calls.
- Every Rust kernel has a slow pure-Python twin under `reference/`, used only in tests; property tests check they agree.
- Before porting (end of phase 5), write a short algorithm spec: parameters, formulas, file formats. The port implements the spec; the Python outputs on fixtures become golden files.
- A Python component moves to Rust early only if profiling shows it blocks experiments, and its interface has not changed for a phase.

## Implementation plan

Eight phases, each with a go/no-go gate; phases 1–5 are the Python prototype with Rust kernels, phase 6 is the Rust port. Phase 3 is the first point where the prototype must at least match fmh-funprofiler, or the project should stop.

| Phase | Deliverable | Language | Gate |
| --- | --- | --- | --- |
| 0. Skeleton | Mixed repo (uv + maturin + Cargo workspace), CI, pre-commit, stub Python CLI. Fixtures, the parity job, coverage and `bench.yml` move to phase 1, when there are kernels to test. | Both | CI green on Linux and macOS arm64 for both languages. |
| 1. Rust kernels (done) | PyO3 module: FASTQ streaming, codon tables (11, 4), six-frame translation, stop-filter frames, reduced alphabets, amino-acid k-mer packing, hashing, FracMinHash filter; batch numpy outputs. Pure-Python reference twins. | Rust + Python tests | Property tests pass (frame symmetry, threshold nesting, synonymous invariance); Rust matches reference; ≥ 1 M reads/min/thread. |
| 2. Index prototype | Build from a MGnify biome subset (DuckDB/Parquet): member k-mers per 90% group, p\_in/p\_out scores, per-cluster floor (non-singletons), connected components, tier 1/tier 2, fingerprints; stored as Parquet + numpy. | Python | Sizes match the analytical estimates; scoring behaves on hand-checked clusters. |
| 3. Query + naive counts | Hit counting and containment; `--sourmash-compat` using the sourmash Python API. Lookup moves to Rust once the layout settles. | Python (+ Rust lookup) | Matches fmh-funprofiler containment in compat mode; ≥ parity in completeness/purity at equal density. **Stop here if not.** |
| 4. Model | Uniqueness-weighted detection, EM per connected component, zero-inflated negative-binomial model, dispersion flag, genome normalisation, dense tier 2. | Python | Clear completeness gain for low-abundance groups over phase 3 at ≤ 2x index size; calibrated intervals on simulations. |
| 5. Evaluation and freeze | Benchmarks vs fmh-funprofiler, HUMAnN, kMermaid, DIAMOND; ablations; divergence ladder. Algorithm spec written; golden outputs recorded. | Python + Nextflow | Defaults chosen; spec reviewed. |
| 6. Rust port | Index build, query, model and CLI in Rust, implementing the spec. Differential tests against the Python golden outputs (exact for counts, tolerance for EM). | Rust | All golden tests pass; ≥ 10x Python end to end; MGnify-scale index builds on one node. |
| 7. Release | Rust binary via cargo-dist, bioconda recipe, Nextflow module for the hybrid profiling pipeline; optional Python wheel of the bindings. | Rust | Tagged release reproduces phase 5 results. |

Out of scope initially: long reads (indels break frames; would need FragGeneScan-style frameshift handling), eukaryotic genes, metatranscriptomes.

## Libraries

Most of the plumbing exists; the amino-acid k-mer hashing and the translation LUT are small enough to write yourself.

| Need | Crate | Note |
| --- | --- | --- |
| FASTA/FASTQ parsing | `paraseq` or `needletail` | paraseq has parallel paired-end record sets; needletail is the mature default. |
| Compression | `niffler` + `flate2` (zlib-rs backend), `zstd` | Transparent gz/zst input. |
| DNA minimizers/hashing | [`simd-minimizers`](https://github.com/rust-seq/simd-minimizers), `packed-seq`, `seq-hash` | Needs AVX2 or NEON; also accepts ASCII text, so could sample amino-acid strings. |
| General hashing | `xxhash-rust` or `wyhash` | For u64-packed amino-acid k-mers. |
| sourmash compatibility | `sourmash` (Rust core) | Protein/Dayhoff/HP encodings and its hash, for the compat mode. |
| Parallelism | `rayon`, `crossbeam-channel` | Reader → worker → reducer pipeline. |
| Index storage | `memmap2`, `rkyv` or `bincode`; `ptr_hash` or sorted arrays | Zero-copy, memory-mapped index. |
| Stats | `statrs` | Poisson/binomial for coverage correction. |
| CLI, errors, logs | `clap`, `anyhow`/`thiserror`, `tracing`, `indicatif` | Standard. |
| Testing | `proptest`, `insta`, `assert_cmd`, `cargo-nextest` | Invariants, golden outputs, CLI tests. |
| Benchmarks | `criterion` or `divan` | Micro-benchmarks of hashing and translation. |
| Frame calling (optional) | [FragGeneScanRs](https://github.com/unipept/FragGeneScanRs) | GPL-3.0: call as a subprocess or keep your crate GPL. |

For binaries using the SIMD crates, pin `target-cpu=x86-64-v3` in release builds; Groot Koerkamp's [post on distributing SIMD binaries](https://curiouscoding.nl/posts/distributing-rust-simd-binaries/) covers bioconda and the `ensure_simd` check.

**Bridge and Python side (phases 0–5)**

| Need | Library | Note |
| --- | --- | --- |
| Rust ↔ Python | `pyo3`, `maturin`, `numpy` (the crate `rust-numpy`) | Zero-copy numpy outputs from batch kernels; release the GIL inside kernels. |
| Environment and builds | `uv` with the maturin build backend | `uv sync` rebuilds the extension; `maturin develop --release` for fast iteration. |
| MGnify subsets | `duckdb` over MGnify's Parquet files | Remote querying of MGnify Proteins via Parquet/DuckDB is documented by MGnify. |
| Tables and storage | `polars`, `pyarrow`, `numpy` | Parquet for index tables during prototyping. |
| Statistics | `scipy` (optimize, stats), `numpy` | EM, zero-inflated negative binomial, dispersion tests. |
| sourmash compat | `sourmash` Python package | Direct comparison with fmh-funprofiler sketches. |
| CLI | `typer` | Prototype CLI mirroring the planned Rust one. |
| Testing | `pytest`, `hypothesis`, `pytest-benchmark` | Property tests mirror the Rust `proptest` ones. |
| Lint, format, types | `ruff`, `mypy` | Run in pre-commit and CI. |
| Plots, notebooks | `matplotlib`, Jupyter or marimo | Evaluation only; kept out of the package. |

## Test and benchmark data

Three tiers: tiny hand-built fixtures in the repo, a small simulated community for CI, and large benchmarks run outside CI.

| Tier | Data | Tests | Where |
| --- | --- | --- | --- |
| Unit fixtures (< 1 MB) | \~20 proteins with their CDS; hand-built reads: in-frame, each of 6 frames, reverse strand, spanning a stop, with synonymous and non-synonymous changes, with N, with an indel. | Translation, frame filter, silent-change invariance, threshold nesting. | `tests/data/`, generated by a script in `scripts/` so they are reproducible. |
| CI community (\~5–10 MB) | 5 complete genomes spanning GC: *S. aureus* (low GC), *E. coli* K-12, *B. subtilis* 168, *P. aeruginosa* PAO1 (high GC), one archaeon. Reads from InSilicoSeq with known origins; truth = per-gene coverage from read origins mapped through annotations. | End-to-end golden output (insta snapshot), GC effect on frame filter, abundance error bounds. | Repo via Git LFS, or Zenodo with a checksum-verified download step. |
| fmh-funprofiler datasets | Their simulated metagenomes and KO ground truth. | Head-to-head at equal density; compat-mode parity. | Zenodo record from their repo. |
| Divergence ladder | Hold out the source genomes' proteins from the index, keeping relatives at \~95/90/80/70% amino-acid identity. | Sensitivity decay vs alphabet, k, spaced seeds. The key experiment for any amino-acid k-mer method. | Built by the eval pipeline. |
| CAMI II (marine, strain madness, plant-associated) | Public simulated metagenomes with genome-level truth. | Functional truth derived by annotating the source genomes (eggNOG-mapper or Bakta) and projecting read origins. | CAMI data portal. |
| Mock communities | ZymoBIOMICS standards (e.g. the gut standard) with published reference genomes. | Real sequencing error and library bias with known composition. | ENA/SRA runs. |
| Real cohort | A subset of HMP2/IBDMDB with published HUMAnN outputs. | Concordance with HUMAnN; runtime at scale. | IBDMDB portal. |
| Negative controls | Shuffled reads, human reads, intergenic-only simulated reads. | False-positive rate. | Generated. |

**Reference databases:** develop on a human-gut biome subset of MGnify Proteins 2026\_07; release on the full set (CC0, so pre-built indexes can be shared). Ground-truth genomes for the CI community and CAMI must be annotated against the same unit definitions, e.g. by mapping their predicted proteins to MGnify clusters.

## Repo setup

One GitHub repo holding a uv-managed Python package built with maturin and a Cargo workspace; pre-commit runs fast checks for both languages, and CI tests both plus Rust-vs-Python parity.

**Create**

Done in phase 0 (repo `timrozday-mgnify/kmer-functional-profiler`, package `kmer_functional_profiler`, GPL-3.0-or-later). `uv init --build-backend maturin` was not used: it puts a single crate at the root, which conflicts with the workspace layout below. Instead `pyproject.toml` points maturin at `crates/py` (`manifest-path`, `python-source = "python"`, `module-name = "kmer_functional_profiler._core"`). Dependencies are added in the phase that first uses them (polars, pyarrow, duckdb in phase 2; scipy in phase 4; sourmash in phase 3; pytest-benchmark and rust-numpy in phase 1).

**Layout**

```text
pyproject.toml        # maturin backend; module <name>._core
uv.lock
Cargo.toml            # workspace, shared deps, lints, release profile
rust-toolchain.toml   # pinned stable
.cargo/config.toml    # target-cpu=x86-64-v3 for release
crates/core/          # kernels: fastx, translation, k-mers, hashing, samplers, lookup
crates/py/            # PyO3 bindings over core -> <name>._core
crates/cli/           # phase 6: final Rust binary (index, query, model)
python/<name>/        # prototype: index build, query, model, typer CLI
python/<name>/reference/  # slow pure-Python twins of the kernels (test oracles)
tests/python/  tests/data/  tests/golden/   # golden outputs from the prototype
scripts/              # fixture generation, data download
eval/                 # Nextflow benchmarks, notebooks
docs/spec.md          # algorithm spec, written at end of phase 5
```

Keeping `crates/core` free of PyO3 means phase 6 reuses the kernels unchanged; only `crates/py` knows about Python.

**Pre-commit** (the `pre-commit` framework, or its Rust port `prek`):

- `pre-commit-hooks`: trailing-whitespace, end-of-file-fixer, check-yaml, check-toml, check-added-large-files (limit 1 MB, so data goes to LFS/Zenodo).
- Python: `ruff check --fix` and `ruff format`; `mypy` on `python/` at the pre-push stage.
- Rust: `cargo fmt --all -- --check`; `cargo clippy --workspace --all-targets -- -D warnings` at the pre-push stage.
- `uv lock --check` so the lockfile never drifts from `pyproject.toml`.
- `typos` for spelling, `taplo` for TOML formatting.

**CI (GitHub Actions)**

- `ci.yml` on push and pull request, three jobs:
  - **rust**: fmt, clippy with `-D warnings`, `cargo nextest`, `cargo doc`, MSRV build, `cargo deny`.
  - **python**: `astral-sh/setup-uv`, `uv sync` (builds the extension in release mode), ruff, mypy, pytest with hypothesis, on Python 3.11–3.13.
  - **parity**: Rust kernels vs `reference/` twins on fixtures and the CI community; from phase 6, the Rust CLI vs `tests/golden/`.
  - Matrix: ubuntu-latest and macos-latest (arm64, covers NEON). Coverage from `pytest-cov` and `cargo-llvm-cov` to Codecov.
- Use `dtolnay/rust-toolchain`, `Swatinem/rust-cache`, and `taiki-e/install-action` for nextest, cargo-deny and llvm-cov.
- `bench.yml`: manual or weekly `pytest-benchmark` and criterion runs; not a merge gate.
- `release.yml` (phase 7): `cargo-dist` binaries on tags, optional wheels via `PyO3/maturin-action`; bioconda recipe after the first tagged release.
- Dependabot for Cargo, uv and Actions; branch protection on `main` requiring `ci.yml`.

## Risks and open questions

The biggest risk is that the gain over fmh-funprofiler with a lower scaled value is too small to justify a new tool.

- **MGnify scale.** A floor on all 1.66×10^9 clusters would need \~160 GB of index; on the 0.45×10^9 non-singletons, \~45 GB. Confirm the full index fits a 64 GB node before phase 4.
- **Component size.** Promiscuous k-mers can chain clusters into one giant component, which makes the EM serial. Measure component sizes on the development subset in phase 2 and tune the N-clusters cut-off.

* **Marginal novelty.** Setting sourmash to scaled = 100 may recover most of the completeness gap at modest cost. Run that baseline in phase 3 before building phase 4.
* **Query density set by the smallest units.** If *n\_min* forces *t\_max* near 1, the query does almost no sparsification. Measure the distribution of *t\_g* on the real database early.
* **Divergence.** Exact amino-acid k-mers miss distant homologs regardless of sampling. If the ladder shows steep decay below 80% identity, reduced alphabets or spaced seeds become mandatory, not optional.
* **Shared k-mers and hierarchy.** EM at protein-cluster level, then aggregation to function, is likely better than EM directly on functions. Untested.
* **Normalisation.** Which single-copy marker set, and whether to report per-genome copies by default.
* **Frame filter at high GC.** Keeps \~3 frames at 70% GC; acceptable, but check false positives there specifically.
* **Decided (phase 1):** k-mer hash is the splitmix64 finalizer of the bit-packed k-mer (5/4/3 bits per residue for protein/Murphy-10/Dayhoff, so k ≤ 12/16/21), a bijection on u64; threshold rule is keep iff hash ≤ max\_hash, with max\_hash = ⌊*t*·2^64⌋ − 1. The sourmash hash is used only in `--sourmash-compat`.
* **Decided:** MGnify90 clusters are the only grouping level; no 30% families.
* **Future work:** 30% families by mapping MGnify90 representatives onto the 128.7 M MGnify30-C2 representatives, as a coarser level for floors, EM partitions and annotation.
* **Decided:** package name `kmer_functional_profiler` (tool name may still change before release); licence GPL-3.0-or-later, so FragGeneScanRs can be linked.
* **Open:** whether KO/eggNOG labels are worth the annotation run, or Pfam suffices.
* **Full-scale runs:** anything over the whole release (subset extraction, index build) ships as a Nextflow pipeline with README and setup scripts for HPC; local work uses samples only.

## Sources

- [Metagenomic functional profiling: to sketch or not to sketch? (Bioinformatics 2024)](https://academic.oup.com/bioinformatics/article/40/Supplement_2/ii165/7749078)
- [fmh-funprofiler repository](https://github.com/KoslickiLab/fmh-funprofiler)
- [fmh-funprofiler preprint v2 (completeness discussion)](https://www.biorxiv.org/content/10.1101/2023.11.06.565843v2.full.pdf)
- [kMermaid (PLoS Comp Biol)](https://pmc.ncbi.nlm.nih.gov/articles/PMC12507277/)
- [UProC (Bioinformatics 2015)](https://academic.oup.com/bioinformatics/article/31/9/1382/200454)
- [Orpheum repository](https://github.com/czbiohub-sf/orpheum)
- [Protein k-mers enable assembly-free microbial metapangenomics](https://dib-lab.github.io/2021-paper-metapangenomes/)
- [Metabuli repository](https://github.com/steineggerlab/Metabuli)
- [Metabuli update: spaced metamers, reduced alphabets, syncmers (bioRxiv 2026)](https://www.biorxiv.org/content/10.64898/2026.03.13.711249v2)
- [simd-minimizers](https://github.com/rust-seq/simd-minimizers) and [docs](https://docs.rs/simd-minimizers/latest/simd_minimizers/)
- [Distributing Rust SIMD binaries](https://curiouscoding.nl/posts/distributing-rust-simd-binaries/)
- [FragGeneScanRs (BMC Bioinformatics 2022)](https://pmc.ncbi.nlm.nih.gov/articles/PMC9148508/)

* [MGnify Proteins resource (2026\_07)](https://docs.mgnify.org/src/docs/mgnify-proteins.html)
* [MGnify: the microbiome sequence data analysis resource in 2023 (NAR)](https://academic.oup.com/nar/article/51/D1/D753/6880769)
* [sylph (Nature Biotechnology 2024)](https://pmc.ncbi.nlm.nih.gov/articles/PMC12339375/)
* [sylph releases (two-stage .syl2db)](https://github.com/bluenote-1577/sylph/releases)
