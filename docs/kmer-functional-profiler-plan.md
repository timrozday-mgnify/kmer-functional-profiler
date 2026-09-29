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
| 2. Index prototype (built; gate needs a real subset) | Build from a MGnify biome subset (DuckDB/Parquet): member k-mers per 90% group, p\_in/p\_out scores, per-cluster floor (non-singletons), connected components, tier 1/tier 2, fingerprints; stored as Parquet + numpy. | Python | Sizes match the analytical estimates; scoring behaves on hand-checked clusters. |
| 3. Query + naive counts (done; gate passed, see Progress log) | Hit counting and containment; `--sourmash-compat` using the sourmash Python API. Lookup moves to Rust once the layout settles. | Python (+ Rust lookup) | Matches fmh-funprofiler containment in compat mode; ≥ parity in completeness/purity at equal density. **Stop here if not.** |
| 4. Model (in progress: baselines, EM, zero-inflated EM adopted, dense tier, *p\_in*-weighted presence, copies-scaled abundance, posterior intervals, ambiguity groups, presence probability, copies error calibrated on real data; simulation benchmark; see Progress log) | Uniqueness-weighted detection, EM per connected component, zero-inflated negative-binomial model, dispersion flag, genome normalisation, dense tier 2. | Python | Clear completeness gain for low-abundance groups over phase 3 at ≤ 2x index size; calibrated intervals on simulations. |
| 5. Evaluation and freeze (in progress: tool benchmark with DIAMOND, fmh-funprofiler, kMermaid, HUMAnN 3.9 and 4 built, results pending on HPC; see Progress log) | Benchmarks vs fmh-funprofiler, HUMAnN, kMermaid, DIAMOND; ablations; divergence ladder. Algorithm spec written; golden outputs recorded. | Python + Nextflow | Defaults chosen; spec reviewed. |
| 6. Rust port and full-scale tuning | Index build, query, model and CLI in Rust, implementing the spec. Differential tests against the Python golden outputs (exact for counts, tolerance for EM). Full-scale cost study: nested all-biome MGnify subsets (1 in 10⁴ to 1 in 10 clusters) and a full-release statistics pass (per-cluster k-mer counts, predicted postings) to fit how storage, build and query cost scale, then the full build. Tune for cost: *t\_base*, *n\_min*, *t\_cap*, dense-tier rate and scope (e.g. non-singletons only), fingerprint width, unit-ID encoding, memory-mapped lookup; re-run the fmh and simulation benchmarks at each candidate to choose defaults on cost vs accuracy. Expect many iterations and some accuracy given up for cost. | Rust | All golden tests pass; ≥ 10x Python end to end; MGnify-scale index builds on one node at a chosen cost/accuracy point, with the accuracy given up versus phase 5 recorded. |
| 7. Release | Rust binary via cargo-dist, bioconda recipe, Nextflow module for the hybrid profiling pipeline; optional Python wheel of the bindings. | Rust | Tagged release reproduces phase 5 results. |

Out of scope initially: long reads (indels break frames; would need FragGeneScan-style frameshift handling), eukaryotic genes, metatranscriptomes.

## Progress log

What each step did, and the choices, results and interpretations behind it, newest phase last. Add an entry with every step; keep superseded choices and say what replaced them.

### Project-wide

* **Decided:** MGnify90 clusters are the only grouping level; no 30% families.
* **Decided:** package name `kmer_functional_profiler` (tool name may still change before release); licence GPL-3.0-or-later, so FragGeneScanRs can be linked.
* **Full-scale runs:** anything over the whole release (subset extraction, index build) ships as a Nextflow pipeline with README and setup scripts for HPC; local work uses samples only.

* **Decided (2026-09-28): full-scale check in phase 6, not phase 4.** No test or benchmark so far uses the full MGnify release: unit tests use fixtures, the simulation 300 synthetic units from 100 MGnify seed proteins, the fmh benchmark KEGG KO indexes (tier 2 72 MB, dense tier 1.09 GB at 1/10), and the only real MGnify build is the pre-fix gut 1-in-10,000 subset (213,645 proteins, about 4×10⁻⁵ of the release). The Python build holds the members table and k-mer tables in memory, so a full-release build needs the Rust port (or a partitioned build) anyway. Phase 6 therefore carries the full-scale cost study and tuning (see the phase table). Expect a lot of tuning there to bring cost down, and some compromises on accuracy; defaults chosen in phase 5 are provisional until phase 6 has priced them, and the phase-4/5 accuracy numbers are the reference that any cost saving is measured against.

### Phase 1

* **Decided (phase 1):** k-mer hash is the splitmix64 finalizer of the bit-packed k-mer (5/4/3 bits per residue for protein/Murphy-10/Dayhoff, so k ≤ 12/16/21), a bijection on u64; threshold rule is keep iff hash ≤ max\_hash, with max\_hash = ⌊*t*·2^64⌋ − 1. The sourmash hash is used only in `--sourmash-compat`.

### Phase 2

* **Decided (phase 2):** the build hashes members three times (distinct k-mers per unit; candidate hashes at *t\_max*; presence of candidates in every unit) so no dense k-mer table is held. *p\_in* is the unweighted fraction of full-length members (all members if a unit has none); sequence weights are deferred until over-sampled lineages are seen to skew it. Score = log2(max(*p\_in*, 0.5/*n*) / *n\_groups*), *n\_groups* counting the unit itself. *p\_in* is quantised to 4 bits in tier 2 (value = unit << 4 | *p\_in* level). Every candidate is kept (no score-based trimming beyond the promiscuity cut, default 64 groups); tier 1 is the top 4 by score per unit, mapped to component ids. Fingerprint layout: leading zero bits implicit, ⌈log2 keys⌉ bucket bits addressing an offsets array, 16 fingerprint bits stored; colliding keys share the union of their unit sets.
* **Phase 2 subsets:** member sequences are spread over the whole 1 TB sequence file, so no real subset fits a laptop download; `workflows/mgnify-subset` extracts one on HPC (default: human gut, 1 in 1000 clusters).
* **Phase 2 gate (gut, 1 in 10,000 clusters: 213,645 proteins, 24,117 clusters):** candidates matched Σ *t\_g*·*n\_kmers* (67,652 vs 67,493); shared k-mers link clusters with a common Pfam (562 of 574 labelled). Four problems, fixed on `index-fixes`:
  1. *Query density.* 90% clusters have few distinct k-mers (median 114), so the floor applied to 8,136 of 8,141 non-singletons and *t\_max* reached 0.57. **Decided:** *t\_g* = max(*t\_base*, min(*t\_cap*, oversample·*n\_min*/*n\_kmers*)) with *t\_cap* = 0.2, so the query samples 20%; clusters with fewer than 40 distinct k-mers get fewer than *n\_min* (430 of 8,141 non-singletons, 4 with none). A 0.05 cap was tried first: 5% sampling but 3,322 clusters below the floor.
  2. *Layout.* 19.4 B per hash, mostly one offsets bucket per key. Now ~4 keys per bucket and the smallest unsigned dtype per array: 7.7 B per hash on the subset (small unit ids; larger ids at full scale need wider values).
  3. *Floor k-mers chosen by hash only.* 30% had *p\_in* = 0 (fragments only). **Decided:** floored clusters sample at 4× (`oversample`) and keep their *n\_min* best-scoring candidates; median *p\_in* rose from 0.13 to 0.87, *p\_in* = 0 fell to 1%.
  4. *Adapter artefacts.* The most shared k-mers were six-frame translations of Illumina TruSeq/Nextera adapters and P5/P7 ends in MGnify proteins (4.1% of members). **Decided:** mask every 6-mer of those translations with X before hashing (shuffled-protein control: 0.08% hit). After masking no k-mer exceeded 64 clusters and the largest component fell from 59 to 6 clusters.

  Rebuilt subset (*t\_cap* = 0.2): 66,337 postings (candidates 210,658 vs 210,940 expected), about 3.7×10⁹ for the full release, so roughly 30–45 GB depending on value widths at full scale.

### Phase 3

* **Decided (phase 3):** the query samples at *t\_max*, looks hashes up in tier 2 and counts a hit for a unit only if hash ≤ *t\_g* of that unit, which also drops most fingerprint false hits. Output per unit: hits, distinct k-mers hit, reads hit, containment (k-mers hit / *m\_g*) and coverage (hits / *m\_g*).
* **Done (phase 3, compat):** `import-sourmash` loads protein signatures (e.g. fmh-funprofiler's KO sketches) as an index; queries on it hash reads with sourmash. On fmh-funprofiler's demo reads (k = 11, scaled = 1000) it reports the same 211 KOs with the same overlaps as `sourmash prefetch`, so the same abundances (fmh-funprofiler's abundance is `f_match_query` = overlap / query hashes, renormalised). Our containment is raw (overlap / *m\_g*); sourmash's `f_query_match` divides it by 1 − (1 − 1/scaled)^(*m\_g*·scaled).
* **Benchmark data (phase 3 gate):** fmh-funprofiler's paper does not publish its simulated metagenomes. It publishes 64-genome-scale inputs to regenerate them with CAMISIM (Zenodo 10055954, CC-BY: KEGG genomes 9.1 GB, `protein_ref_db_giant.faa` 3.2 GB, `present_genes_and_koids.csv` 0.26 GB) and its mean metrics (sourmash k = 11: purity 0.98, completeness 0.61 at 0% error).
* **Benchmark pipeline:** `workflows/fmh-benchmark` builds our KO indexes from the same KEGG proteins as the sketches, simulates metagenomes with InSilicoSeq (64 genomes, ~1 Gbp, lognormal abundances) instead of CAMISIM, derives KO truth by mapping reads back with minimap2 (a KO is present if a read overlaps one of its genes, the paper's rule) and scores purity and completeness for fmh-funprofiler (compat) and our indexes. Runs on HPC; results under Phase 3 below and in `workflows/fmh-benchmark/README.md`.
* **Phase 3 gate: passed (10 InSilicoSeq novaseq metagenomes, 64 genomes, ~1 Gbp, k = 11).** Means over seeds 1..10, sd ≤ 0.015; full table in `workflows/fmh-benchmark/README.md`.

  | Index | Count | min\_hits | Purity | Completeness | Low 25% |
  | --- | --- | --- | --- | --- | --- |
  | fmh\_compat (fmh-funprofiler, scaled 1000) | kmers\_hit | 1 | 0.975 | 0.688 | 0.295 |
  | kfp\_s1000 (ours, same density) | kmers\_hit | 1 | 0.985 | 0.677 | 0.285 |
  | kfp\_s1000\_floor8 | kmers\_hit | 1 | 0.976 | 0.719 | 0.404 |
  | kfp\_s100 (*t\_base* 0.01) | kmers\_hit | 1 | 0.951 | 0.960 | 0.853 |
  | **kfp\_s100** | **kmers\_unique** | **1** | **0.985** | **0.955** | **0.837** |

  *Interpretation.* At equal density we match fmh-funprofiler (parity, as the gate requires); the floor adds 4 points of completeness and 11 on the least-covered quarter. The real gain is density: scaled 100 lifts completeness from 0.69 to 0.96. Its false positives were mostly modular PKS/NRPS KOs whose shared domains carry identical k-mers, so a hit on a shared k-mer is weak evidence.

  **Decided:** add gather-style reassignment (`query.gather`: repeatedly take the unit with the most unassigned hit k-mers, scaled by 1/*t\_g*) as the phase-4 detection baseline; profiles report `kmers_unique` beside `kmers_hit`. Gather removes 70% of false positives (348 → 105 per sample) for 31 true KOs lost, and adds no measurable query time. Recommended setting: `kfp_s100`, `kmers_unique` ≥ 1. `min_hits` > 1 trades too much low-coverage completeness (0.49 at 2 on scaled 1000). The floor does not help at scaled 100 (within 0.005 of *n\_min* 0), so `kfp_s100_floor8` was dropped from the defaults.

  *Caveats.* InSilicoSeq, not CAMISIM, so not directly comparable with the paper's 0.98 / 0.61; KO units from KEGG proteins, not MGnify90 clusters. The phase-2 gate on a real MGnify subset (index size at scale) is still open.

### Phase 4

* **Phase 4, step 1 — EM quantification (`query.em`).** On the units gather keeps, each hit k-mer's count is Poisson with mean Σ λ\_g over the units holding it, and every kept k-mer of a unit (*m\_g*, hit or not) counts in its expectation; the multiplicative EM update gives the MLE. Profiles gain `coverage_em` (0 for units gather drops). Choices: EM on k-mer counts, not read equivalence classes (see Assignment); all components in one sparse product, since they share no k-mers and so update independently (per-component parallelism only matters in the Rust port); stop at max relative change 10⁻⁶ or 1,000 iterations. Tests: a two-unit case against the closed-form MLE, and Σ λ\_g·*m\_g* = observed hits on the fixture reads.

  *Not yet:* zero-inflation and *p\_in* weights, fingerprint false-hit background, uniqueness-first and winner-take-all baselines, dispersion flag, genome normalisation, dense tier 2. *Result pending:* abundance accuracy has not been benchmarked; the fmh benchmark scores detection only (truth has per-KO read and base counts to compare against).

* **Phase 4, step 2 — abundance scoring in the fmh benchmark.** TRUTH adds per-KO `depth` = Σ over the KO's genes of aligned bases / gene length (read depth × copies, the quantity λ should track up to a constant). SCORE pairs each detection count with an abundance estimate: `coverage` (hits / *m\_g*) with `kmers_hit`, `coverage_em` with `kmers_unique`. On the detected KOs it reports `spearman_tp` (rank correlation with depth over true positives, so it measures quantification separately from detection) and `l1` (L1 between relative abundances over true ∪ detected KOs, so misses and false positives count). TRUTH now takes `bench.py` as an input, so `-resume` recomputes the truth tables (one minimap2 pass per sample). Checked on the test profile only; *result pending* on HPC.

  *Expectation to test.* A KO unit is the union of k-mers of every KEGG gene with that KO, and a sample carries only a few of those genes, so λ = hits / *m\_g* is diluted by the share of the unit's k-mers the sample holds, which varies widely between KOs. Plain EM does not correct this; the zero-inflated model (λ from present k-mers, as in sylph) should. If `spearman_tp` is low for both estimates, that is the case for doing zero-inflation next.

* **Phase 4, step 3 — winner-take-all and uniqueness-first baselines (`query.assign_best`).** With a fixed score per unit, both rules reduce to one pass: each hit k-mer goes to the holding unit with the highest score (ties to the lowest id). Gather differs only in that its scores fall as k-mers are taken.
  - *Winner-take-all* (sylph): score = containment, k-mers hit / *m\_g*. One pass; sylph's recomputation of containment on the reassigned k-mers is not repeated.
  - *Uniqueness-first*: score = Σ over the unit's hit k-mers of 1 / (units holding the k-mer), an IDF weight, divided by *t\_g* as this plan requires. The plan's Poisson test of unique hits against a false-hit background waits for that background model.

  Profiles gain `kmers_wta`, `coverage_wta`, `kmers_ufirst`, `coverage_ufirst` (coverage = hits on the assigned k-mers / *m\_g*). The benchmark scores each count at every `min_hits` for detection and pairs it with its coverage for abundance, so all four rules (gather + EM, winner-take-all, uniqueness-first, none) run on the same hit table as planned. Tests: hand-worked cases where the two rules disagree (a fully contained small unit wins under winner-take-all, loses under uniqueness-first unless it is sampled more sparsely), and on fixture reads every hit k-mer and hit is assigned exactly once. *Result pending* on HPC.

* **Phase 4, step 4 — simulated abundance benchmark (`workflows/sim-benchmark/sim.py`).** The fmh benchmark needs HPC, so a local simulation with exact truth guides model work between HPC runs (about 3 s per seed). 100 families from the MGnify sample, each with 3 paralogous units at 85–95% identity to the seed (shared k-mers), 4 members per unit at 97%; half the units present, each as one strain at 100/95/90/85% identity to its centroid with lognormal depth; random-codon CDSs, 150 bp single-end reads, 0.2% substitutions, 20% random decoy reads. Index configs: `dense` (*t\_base* 1), `s10` (0.1), `floor` (MGnify defaults). Metrics: unit purity and completeness; Spearman over true positives; L1 over units and over families (the Pfam-level view); sd of log(estimate / depth); `bias_<identity>` = median estimate / depth at that strain identity over the median at 100%.

  ```bash
  uv run python workflows/sim-benchmark/sim.py --seeds 5
  ```

* **Phase 4, step 5 — zero-inflated EM (`em(zero_inflated=True)`, `coverage_zi`, `present_zi`).** Each unit gets a present fraction π\_g beside λ\_g: only π\_g·*m\_g* of its kept k-mers occur in the sample, the rest are structural zeros. Hits are shared in proportion to λ·π; π\_g = (unit's share of hit k-mers) / (*m\_g*·(1 − e^−λ)), capped at 1. With no shared k-mers the fixed point is sylph's zero-truncated Poisson MLE; with no excess zeros it reduces to plain EM (both tested).

  Simulation, 5 seeds (min\_hits 1):

  | Config | Rule | Purity | Spearman | L1 | L1 family | log-ratio sd | bias 95% | bias 90% | bias 85% |
  | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
  | dense | none (all hits) | 0.58 | 0.75 | 0.55 | 0.41 | 0.74 | 0.58 | 0.34 | 0.22 |
  | dense | gather + EM | 0.91 | 0.79 | 0.50 | 0.41 | 0.78 | 0.56 | 0.30 | 0.16 |
  | dense | winner-take-all | 0.90 | 0.76 | 0.52 | 0.41 | 0.84 | 0.55 | 0.28 | 0.14 |
  | dense | uniqueness-first | 0.90 | 0.76 | 0.52 | 0.41 | 0.84 | 0.55 | 0.28 | 0.14 |
  | **dense** | **gather + ZI EM** | 0.91 | **0.95** | **0.23** | **0.18** | **0.52** | **0.94** | **0.94** | **0.91** |
  | s10 | gather + EM | 0.97 | 0.75 | 0.54 | 0.45 | 0.83 | 0.56 | 0.29 | 0.16 |
  | s10 | gather + ZI EM | 0.97 | 0.92 | 0.21 | 0.17 | 0.65 | 0.95 | 0.98 | 0.89 |
  | floor | gather + EM | 0.96 | 0.71 | 0.60 | 0.49 | 0.77 | 0.54 | 0.32 | 0.25 |
  | floor | gather + ZI EM | 0.96 | 0.84 | 0.39 | 0.32 | 0.72 | 0.94 | 0.92 | 0.84 |

  *Interpretation.* Divergence, not shared k-mers, dominates abundance error. Without zero-inflation, estimates fall roughly as identity^k (0.95^11 = 0.57, 0.90^11 = 0.31, 0.85^11 = 0.17), because a divergent strain carries only that share of its unit's k-mers; members' own variation dilutes further (median π ≈ 0.25 even at 100% identity). Zero-inflation removes the bias at every config. The three simple rules and plain EM differ little from each other (EM best by 0.03 Spearman); all gain purity over counting every hit (0.58 → 0.91 dense), and gather's detection is what does that. Residuals by depth (dense): median log-ratio within ±0.07 above depth 1 (−0.14 at ≤ 1), with sd 1.15 at depth ≤ 1, 0.45 at 1–3, 0.14 at 3–10, 0.09 above 10. With the floor index units have 2–4 hit k-mers and sd stays 0.25–1.0 at every depth. So the remaining error is thin data: low coverage, few k-mers.
  **Decided:** `coverage_zi` is the recommended abundance estimate; gather stays the detection rule. `coverage_em`, `coverage_wta` and `coverage_ufirst` remain as baselines for the fmh benchmark.

* **Phase 4, step 6 — empirical-Bayes prior on the present fraction (`fit_present_prior`, `coverage_zib`, experimental).** A Beta prior on π, matched by moments to units with λ ≥ 3 (where π is measured), with π updated by EM over unhit k-mers' presence. Result (5 seeds): dense log-ratio sd 0.52 → 0.42 and L1 0.23 → 0.21 at a small bias cost (85%: 0.91 → 0.88); s10 sd 0.65 → 0.46 but bias at 85% 0.89 → 0.76; floor worse everywhere (Spearman 0.84 → 0.74, bias at 85% 0.84 → 0.39), because with ~8 kept k-mers the prior outweighs the data and pulls divergent strains' π up, so λ down. **Decided:** not the default. It helps only where units have many k-mers, which is what the dense tier 2 will give, so revisit it there; `coverage_zib` stays in profiles so the HPC fmh run measures it on real data.

  *Next, by this evidence:* (1) the fmh benchmark on HPC, to confirm ZI on real genomes (KO units are far larger unions than MGnify90 clusters, so π is smaller); (2) dense tier 2 for candidate units, since k-mer count per unit limits precision more than the model does; (3) *p\_in*-weighted presence (core k-mers more likely present); (4) confidence intervals, which low-coverage units need.

* **Phase 4, step 7 — dense tier (`IndexParams.t_dense`, `--t-dense`).** The plan's two-stage lookup: the sparse tier 2 detects, a denser table quantifies.
  - *Build:* a fourth hashing pass keeps each unit's k-mers with hash ≤ max(*t\_dense*, *t\_g*) (`max_hash_dense`), minus promiscuous ones (same `max_groups` cut, counted over all units), with quantised *p\_in*, in a `dense` table of the tier-2 layout; `m_dense` per unit. The per-unit max makes each unit's dense set contain its tier-2 set. A uniform *t\_dense* was tried first: at 1 in 50 it gave floored units (tier 2 up to 0.2) fewer k-mers than tier 2, and ZI Spearman fell below the floor alone (0.83 vs 0.84), with completeness 0.78.
  - *Query:* after gather on tier-2 hits, the reads are streamed again at the dense rate and only the detected units' hits are kept; `coverage_em`, `_zi`, `_zib` are then fitted on those over `m_dense`, and `kmers_dense` reports the dense k-mers hit. Detection and the other baselines are unchanged. A second pass rather than buffering hashes, since the kernel runs at ~54 M reads/min and buffered hashes at dense rates would not fit in memory.
  - *EM fix found on the way:* with more shared k-mers visible, EM explains some gather-detected units away entirely; their weight decays towards 0, and zero-inflated EM divided 0 by 0 (NaN). Divisions are now guarded and units left with < 10⁻³ expected hits (`EXPLAINED_AWAY`) get coverage 0. Up to 1.5% of true units are explained away this way (at 1 in 50), so they drop out of true positives (floor 0.879 → 0.867–0.879 completeness).
  - *Test:* sparse tier 2 plus a fully dense tier gives exactly the EM and ZI estimates of an index that is dense throughout.

  Simulation, 5 seeds, zero-inflated EM (`gather_zi`); tables = size of the lookup tables relative to `floor`:

  | Config | Tables | Spearman | L1 | L1 family | log-ratio sd | bias 85% |
  | --- | --- | --- | --- | --- | --- | --- |
  | floor | 1× | 0.835 | 0.391 | 0.317 | 0.719 | 0.84 |
  | floor + dense 1/50 | 3.5× | 0.891 | 0.325 | 0.256 | 0.715 | 0.90 |
  | floor + dense 1/20 | 3.9× | 0.897 | 0.313 | 0.241 | 0.698 | 0.89 |
  | floor + dense 1/5 | 11× | 0.922 | 0.251 | 0.198 | 0.519 | 0.92 |
  | floor + dense 1/1 | 55× | 0.940 | 0.231 | 0.185 | 0.443 | 0.93 |
  | dense index (reference) | 54× | 0.948 | 0.228 | 0.178 | 0.517 | 0.91 |

  *Interpretation.* The dense tier recovers most of the gap between the floor index and a fully dense one: at 1 in 5, Spearman 0.84 → 0.92 and L1 0.39 → 0.25, with detection unchanged. With the dense tier, the empirical-Bayes prior no longer hurts (`gather_zib` at 1 in 5: Spearman 0.925, sd 0.42, but bias at 85% 0.81), as step 6 predicted; it still trades divergence bias for precision. The size ratios here overstate the cost at MGnify scale in one way and understate it in another: the simulated units are all non-singleton and floored (8 tier-2 k-mers each), while MGnify's 73% singletons sit at *t\_base* = 0.001 in tier 2 and would be raised to *t\_dense* in the dense table. The dense table is only probed for detected units and can stay on disk (memory-mapped), so its size bounds disk, not RAM, but the phase-4 gate (≤ 2× index size) is not met at any useful rate here.
  **Decided:** keep the dense tier as an option (default off) and measure it on real data before choosing a rate: `kfp_s100_d10` (tier 2 at 0.01, dense at 0.1) added to the fmh benchmark. Open: restrict the dense tier to non-singleton clusters, or to a rate that keeps it ≤ 1× tier 2 at MGnify scale; decide once the MGnify subset build reports `dense_bytes`.

* **Fix (fmh benchmark):** `index_*/meta.json` (index sizes, incl. `dense_bytes`) was never published: Nextflow folds a file inside a declared directory output into that directory, which did not match the publish pattern. INDEX and IMPORT\_SKETCHES now copy `meta.json` out and publish it to `<outdir>/index_<name>/meta.json`.

* **Phase 4, step 8 — *p\_in*-weighted presence (`em_pin`, `coverage_zip`, `present_zip`, experimental).** Zero-inflated EM with a per-k-mer presence probability from *p\_in(x)* (the plan's *z\_x* ~ Bernoulli(*π\_x*), *π\_x* informed by *p\_in*). The index stores per unit a histogram of kept k-mers over the 16 quantised *p\_in* levels (`pin_hist`, `pin_hist_dense`), so unhit k-mers are summed by level; the query now carries each hit's `pin_q` from tier 2 / the dense table. *p\_in* is the level midpoint, so no k-mer is certainly absent. Hit k-mers count as present by their share, as in `em`. Tests: all k-mers at one level reproduces zero-inflated EM; a k-mer shared between a unit where it is core and one where it is private goes mostly to the first.
  - *First link, logistic:* presence = σ(α\_g + logit *p\_in*), α\_g fitted by a Newton step per iteration. Worse than plain ZI everywhere (5 seeds): dense Spearman 0.948 → 0.936, bias at 85% 0.91 → 0.82, sd 0.52 → 0.55. Divergence scales presence multiplicatively (each k-mer survives with probability ≈ identity^k whatever its *p\_in*), and a logit shift cannot do that: bringing core k-mers down drives low-*p\_in* ones to 0.
  - *Adopted link, multiplicative:* presence = *s\_g*·*p\_in(x)* (*s\_g* = 1: the strain carries k-mers like a random member; < 1: divergent), capped so presence ≤ 1. Closed-form M-step: *s\_g* = expected present k-mers / Σ *p\_in*.

  Simulation, 5 seeds, `gather_zi` → `gather_zip`:

  | Config | Spearman | L1 | log-ratio sd | bias 85% |
  | --- | --- | --- | --- | --- |
  | dense | 0.948 → 0.947 | 0.228 → 0.229 | 0.517 → 0.411 | 0.91 → 0.88 |
  | s10 | 0.919 → 0.918 | 0.213 → 0.209 | 0.653 → 0.481 | 0.89 → 0.84 |
  | floor | 0.835 → 0.829 | 0.391 → 0.394 | 0.719 → 0.672 | 0.84 → 0.82 |
  | floor + dense 1/50 | 0.891 → 0.892 | 0.325 → 0.315 | 0.715 → 0.496 | 0.90 → 0.86 |
  | floor + dense 1/5 | 0.922 → 0.923 | 0.251 → 0.244 | 0.519 → 0.405 | 0.92 → 0.87 |
  | floor + dense 1/1 | 0.940 → 0.940 | 0.231 → 0.227 | 0.443 → 0.348 | 0.93 → 0.89 |

  *Interpretation.* *p\_in* tells the model which zeros are expected (private k-mers) and which are informative (core k-mers missed), so per-unit estimates scatter 20–30% less; ranks and L1 barely move, because they are dominated by the large depth range. The residual bias (≈ 4% at 85% identity) fits the simulation's cluster model, where presence given *p\_in* is not proportional to *p\_in*: members are star-like mutants of a centroid, so a k-mer at *p\_in* = 1/4 is almost never in the centroid, while the model expects it in a quarter of strains. Against the empirical-Bayes prior (step 6), it gives a similar precision gain at a quarter of the bias cost and does not break the floor index. **Decided:** not the default yet: `coverage_zip` joins the fmh benchmark, where KO units are unions over many species and *p\_in* is low for most k-mers, a very different regime. Adopt it if it holds there.

* **Fix (fmh benchmark):** INDEX and IMPORT\_SKETCHES now take the package sources as an input, like PROFILE, so `-resume` rebuilds indexes when the build code changes; cached indexes without `pin_hist` made every PROFILE fail. Indexes from runs before this change are rebuilt once.

* **Phase 4, step 9 — fmh benchmark with abundance (HPC, PR #13 code, 10 InSilicoSeq metagenomes).** Detection reproduced the phase-3 numbers exactly. Abundance against truth depth (Σ over the KO's genes of aligned bases / gene length), true positives, `kmers_unique` ≥ 1:

  | Index | Abundance | Spearman | L1 |
  | --- | --- | --- | --- |
  | fmh\_compat | coverage (hits / *m\_g*) | 0.08 | 1.41 |
  | fmh\_compat | ZI EM | 0.40 | 0.97 |
  | kfp\_s100 | coverage | 0.32 | 1.35 |
  | kfp\_s100 | EM / winner-take-all / uniqueness-first | 0.35 / 0.34 / 0.37 | 1.25 / 1.32 / 1.23 |
  | kfp\_s100 | ZI EM | 0.58 | 0.94 |
  | kfp\_s100 | ZI EM + prior | 0.53 | 0.98 |
  | kfp\_s100\_d10 | ZI EM (dense tier 1/10) | 0.63 | 0.93 |

  As in the simulation, zero-inflation is the largest gain and the three simple rules and plain EM are indistinguishable; the prior does not help. But the level is far below the simulation's 0.95. *Why:* a KO unit is the union of its genes from many species, and a sample holds several of them at different depths; truth sums depth over those copies, while ZI coverage is depth per present k-mer, about one copy's. Scaling by copies present (offline, from the profiles, with copies ≈ `present_zi` × `n_members`, valid when members share few k-mers as KO members do) gives Spearman **0.935** (kfp\_s100) and **0.964** (d10), log-ratio sd 0.63 and 0.45. Hits / *t\_g* / (k-mers per member), with no model at all, gives 0.92: most of the gap was the estimand, not the model.
  **Decided:** profiles report `copies_zi` = present k-mers / `pin_sum` and `abundance_zi` = `coverage_zi` × `copies_zi`. `pin_sum` (`pin_sum_dense`) is Σ *p\_in* over a unit's kept k-mers, computed at build from unquantised *p\_in*: how many kept k-mers an average member holds, so present / `pin_sum` counts member-equivalents present. The 4-bit *p\_in* levels cannot give this for KOs (most KO k-mers have *p\_in* < 1/15). The two estimates answer different questions: `coverage_zi` is depth per copy and stays the per-cluster estimate for MGnify90 units (one strain each; the simulation shows `abundance_zi` there brings back the divergence bias, 0.16 at 85%, because a divergent strain looks like fewer copies); `abundance_zi` is total depth over copies, which is what the KO-unit benchmark measures. For the real tool, function-level abundance is Σ over clusters of `coverage_zi`, per the plan (EM at cluster level, then aggregate). The fmh benchmark now scores `abundance_zi` to confirm the exact form (needs a rerun; these indexes predate `pin_sum`).

  *Dense tier cost on real data* (kfp\_s100\_d10 vs kfp\_s100): tables 1.09 GB vs 72 MB (15×; 7.9 B vs 5.2 B per posting), index build 36.5 GB vs 7.8 GB peak memory, query 2m16–2m46 vs 23 s and 5 GB vs 0.58 GB (the prototype materialises 8 B per key for lookup). For +0.05 Spearman (ZI) or +0.03 (copies-scaled), and a 30% lower log-ratio sd. **Decided:** the dense tier stays off by default; revisit after the Rust lookup (memory-mapped, no materialised keys) and with a rate or singleton restriction chosen from MGnify-subset `dense_bytes`.

  *Also:* `em_pin` read quantised *p\_in* levels as (*l* + 0.5)/16, but the build quantises as round(15·*p\_in*); now *l*/15 with level 0 at 1/30. Simulation after the fix (`gather_zip` vs `gather_zi`): dense sd 0.52 → 0.42 as before, floor Spearman 0.835 → 0.806 (was 0.829). The `gut-lin10000` archive from the same HPC batch is the pre-fix phase-2 build (*t\_max* 0.57, 19.4 B per hash), already recorded under Phase 2.

* **Phase 4, step 10 — confidence intervals (`bootstrap_zi`, `query --bootstrap B`).** 95% percentile intervals for `coverage_zi` and `abundance_zi` (`*_lo`, `*_hi`) from a Poisson bootstrap over reads: each replicate weights every read pair by a Poisson(1) draw, rebuilds the k-mer counts and refits zero-inflated EM (sharing between units included). Reads, not k-mers, are resampled because one 150 bp read hits ~40 consecutive amino-acid k-mers of a unit, so k-mer counts are correlated and a Poisson/Fisher interval would be too narrow. To do this the query now keeps hits per (unit, k-mer, read) for detected units, in both passes. Off by default (B = 0); tests: intervals bracket the estimate, leave every other column unchanged, and are reproducible (fixed seed).
  - *Simulation calibration* (5 seeds, B = 100, `coverage_zi` against true k-mer coverage = depth × a per-sample scale, the median estimate / depth over strains at 100% identity and depth > 5):

  | Config | Coverage (95% nominal) | Coverage, depth ≤ 2 | Median width, log(hi / lo) |
  | --- | --- | --- | --- |
  | dense | 0.964 | 0.962 | 1.17 (×3.2) |
  | s10 | 0.973 | 0.979 | 1.42 (×4.1) |
  | floor | 0.968 | 0.983 | 2.18 (×8.8) |
  | floor + dense 1/50 | 0.958 | 0.966 | 1.63 (×5.1) |
  | floor + dense 1/5 | 0.950 | 0.925 | 1.15 (×3.2) |

  *Interpretation.* Calibrated to slightly conservative everywhere, including at low depth where the point estimate is least reliable; widths track the information per unit (the floor index's ~8 k-mers per unit give ×9 intervals). This meets the phase-4 gate's "calibrated intervals on simulations" for `coverage_zi`. Caveat: the scale is estimated from the same sample; the ±1–2% conservatism is within what that allows.
  - *fmh benchmark:* PROFILE runs with `--bootstrap` (`params.bootstrap`, default 100) and SCORE reports `ci_cover` and `ci_width` for estimates with intervals, against truth depth on the estimate's scale (median estimate / depth over true positives). Result pending on HPC; the dense-tier profile will be the slowest (100 EM refits over its larger hit table).

* **Phase 4, step 11 — fmh benchmark with copies-scaled abundance (HPC, PR #14 code).** Detection and every earlier estimate reproduced the previous run exactly. `abundance_zi` against truth depth (`kmers_unique` ≥ 1, true positives, 10 seeds):

  | Index | `coverage_zi` Spearman / L1 | `abundance_zi` Spearman (sd) / L1 |
  | --- | --- | --- |
  | kfp\_s1000 | 0.373 / 0.97 | 0.807 (0.011) / 0.57 |
  | kfp\_s1000\_floor8 | 0.369 / 1.02 | 0.827 (0.014) / 0.56 |
  | kfp\_s100 | 0.576 / 0.94 | **0.957 (0.003) / 0.22** |
  | kfp\_s100\_d10 | 0.627 / 0.93 | **0.988 (0.002) / 0.14** |

  *Interpretation.* The exact form (present k-mers / Σ *p\_in*) beats the offline approximation (0.935), and on KO units abundance is now as good as detection: kfp\_s100 detects 95.5% of KOs at 98.5% purity and ranks their abundance at 0.96. The dense tier's gain is larger here than for per-copy coverage (L1 0.22 → 0.14, 36% lower), but its cost is unchanged (15× tables, ~6× query time); still off by default, and a candidate for the phase-6 cost study. At scaled 1000 the floor helps abundance a little (0.81 → 0.83).
  - `fmh_compat`: `abundance_zi` equals plain EM (0.09), as expected: imported sketches carry no members, so *p\_in* = 1 and `pin_sum` = *m\_g*, and copies reduce to the present fraction. fmh-funprofiler's own sketches cannot give copy-aware abundance.
  - `coverage_zip` equals `coverage_zi` on our KO indexes (0.576 vs 0.576): nearly all KO k-mers have *p\_in* < 1/15 and fall in quantised level 0, so presence cannot vary between them. On `fmh_compat` it equalled plain EM, which exposed a bug: with level 15 read as exactly 1, a k-mer present with probability 1 stays so, because EM then reads every unhit copy as present. **Fixed:** level probabilities are clipped to [1/30, 1 − 1/30]; a test covers all-core units. The simulation is unchanged by the fix (its `coverage_zip` bias is structural, not this bug). `coverage_zip` stays experimental: it cannot help KO-like units at 4-bit *p\_in*, and its intended case, MGnify90 clusters with varied *p\_in*, has no truth-bearing benchmark until phase 6.
  **Decided:** recommended setting for KO-like units is `kfp_s100`, detection `kmers_unique` ≥ 1, abundance `abundance_zi`; for MGnify90 clusters `coverage_zi`, summed over clusters for function-level abundance. Intervals (PR #15) have not been run on HPC yet.

* **Phase 4, step 12 — bootstrap intervals on real data (HPC, PR #15 code, B = 100).** Detection and all point estimates reproduced the previous run; the `em_pin` fix shows as expected (`fmh_compat` `coverage_zip` 0.09 → 0.395, in line with `coverage_zi`). Interval coverage of truth depth on the estimate's scale (nominal 95%, true positives):

  | Index | `abundance_zi` coverage | median width, log(hi / lo) | `coverage_zi` coverage |
  | --- | --- | --- | --- |
  | kfp\_s1000 | 0.556 | 1.28 | 0.529 |
  | kfp\_s100 | 0.607 | 0.62 | 0.333 |
  | kfp\_s100\_d10 | 0.717 | 0.38 | 0.163 |

  (`coverage_zi` estimates depth per copy, not the KO total, so its low coverage is expected.) *Why `abundance_zi` intervals are too narrow:* the estimate is unbiased at every depth (median log residual within ±0.03 per depth quartile on kfp\_s100), but the bootstrap half-width shrinks with depth (0.88 → 0.12 from the lowest to the highest quartile) while the residual sd does not (0.69 → 0.34). Coverage therefore falls with depth, 0.80 → 0.43 on kfp\_s100 and 0.88 → 0.49 with the dense tier. The bootstrap measures read-sampling noise only; what it misses is a per-unit model error. Adding a constant σ on the log scale in quadrature restores calibration offline: σ = 0.2 gives 0.95 (d10) and 0.88 (s100), σ = 0.3 gives 0.97 and 0.93 (s1000 needs more: 0.82 at 0.3). The simulation had no such error (single-strain units, so copies are exact), which is why it was calibrated there. The likeliest source is the copies estimate: it divides present k-mers by an *average* member's kept k-mers, but the KO genes in a sample are specific members whose lengths vary between species.
  *Cost:* the bootstrap multiplies PROFILE time by 4–6 on kfp\_s100 (23 s → 1m31–2m09) and 8–11 with the dense tier (2m16–2m46 → 18–30 min; peak memory 6 GB).
  **Decided:** intervals are not calibrated on real data yet; `--bootstrap` stays off by default in the tool, and results must not be read as 95% intervals. Next: add a per-unit model-error term, from the spread of member lengths (kept k-mers per member, measurable at build) divided by the square root of the copies present, plus a small calibrated floor, and test it on the next HPC run; reduce the benchmark's default replicates (or bootstrap only non-dense indexes) to cut cost.

* **Phase 4, step 13 — posterior intervals and ambiguity groups (`posterior_zi`, `query --draws D`; replaces the read bootstrap).** *Intent (clarified):* intervals should carry the uncertainty of k-mer matches, i.e. of splitting hits between similar units, including off-target matches from sequencing errors, which mostly land on close relatives. An LCA over a similarity tree was considered and not needed: the likelihood already shows which units the data cannot separate. Aggregating to functions is out of scope for the tool (downstream analysis); only the benchmarks do it.
  - *Method:* each of D draws reweights reads (Poisson(1) per read pair), refits zero-inflated EM (read-sampling uncertainty), then runs 10 Gibbs sweeps on those counts and keeps the last state (split uncertainty). A sweep splits each hit k-mer's count among its holders multinomially ∝ coverage × present fraction; updates coverage from the counts on the k-mers each unit was given, a zero-truncated Poisson likelihood free of the present fraction (Metropolis on log coverage, 5 steps, Gamma(1, 0.01) prior); and draws the present fraction exactly on a 512-point grid from hit k-mers ~ Binomial(*m*, π(1 − e^−λ)). A unit given no hits in a draw is absent in it. *Ambiguity groups:* holders of a shared hit k-mer whose coverage draws correlate below −0.5 are linked; connected units form a group with `ambiguity_group`, `group_size` and interval on the group total (`group_coverage_zi_*`, `group_abundance_zi_*`).
  - *Tried first, and why they failed:* (1) Gibbs with data augmentation of presence and conjugate Gamma/Beta draws, reads reweighted every sweep: intervals too narrow (dense coverage 0.85), because coverage and present fraction trade off along a ridge the chain crawls along. (2) The same with the zero-truncated coverage update: still too narrow (0.81), because a single chain only partly follows each sweep's reweighted counts and so averages over them. Nesting (refit per draw, short chain from it) fixed both.
  - *Simulation with near-identical twins* (5 seeds, D = 100, `--twins 0.3`: 30% of families gain a unit that is a 99% copy of another, each present independently):

  | Config | Coverage | depth ≤ 2 | twins | width | TP grouped | group cover | twin pairs grouped |
  | --- | --- | --- | --- | --- | --- | --- | --- |
  | dense | 0.933 | 0.964 | 0.682 | 1.21 | 0.015 | 0.53 | 0.10 |
  | s10 | 0.959 | 0.982 | 0.682 | 1.71 | 0.014 | 0.88 | 0.22 |
  | floor | 0.969 | 0.966 | 0.883 | 2.66 | 0.000 | – | 0.00 |
  | floor + dense 1/5 | 0.946 | 0.954 | 0.863 | 1.35 | 0.037 | 0.80 | 0.21 |

  *Interpretation.* Intervals are calibrated overall and at low depth. Most twins are resolved by gather (the unsupported twin gets no gather k-mers and is not reported), which is correct. The hard case is both twins present at different depths: the split is biased, the true twin's interval covers 68–88%, and only 10–22% of such pairs form a group, because their draws are not anti-correlated enough. The likely cause is the split step's approximation (weights ∝ coverage × present fraction, ignoring the holders' joint presence on each shared k-mer).
  - *Removed:* `bootstrap_zi` and `--bootstrap` (`--draws` replaces them; the fmh benchmark's `params.bootstrap` is now `params.draws`). *fmh benchmark:* SCORE adds `fp_grouped` (share of false-positive KOs grouped with a true KO) and `group_cover` for `abundance_zi`; result pending on HPC. The copies model error found in step 12 is separate and not addressed here.
  - *Next:* model the holders' joint presence exactly in the split (most shared k-mers have two holders, so enumeration is cheap), and check whether twin coverage and grouping improve; tune the anti-correlation threshold on the twin simulation.

* **Fix (fmh benchmark):** SUMMARY failed when a metric was empty in one sample's score file (e.g. `group_cover` when no ambiguity group formed): Polars reads an all-empty column as String and cannot stack it with Float64. Metrics are now cast to Float64 on read; a regression test stacks a file with an empty metric.

* **Phase 4, step 14 — posterior intervals and ambiguity groups on real data (HPC, PR #17 + #18 code, D = 100).** Point estimates unchanged. `abundance_zi` interval coverage of truth depth (nominal 95%; the step-12 bootstrap in brackets): kfp\_s1000 0.69 (0.56), kfp\_s100 0.72 (0.61), kfp\_s100\_d10 0.75 (0.72), fmh\_compat 0.43 (0.36); widths grew by 10–30%. Better, still not calibrated: the per-KO model error of the copies estimate (step 12) is untouched by sampling. PROFILE time: kfp\_s100 2m17–2m39 (23 s without draws), dense tier 20–25 min, kfp\_s1000 39–44 s.
  - *Ambiguity groups did not form:* none on kfp\_s100 or kfp\_s100\_d10 in any of 10 samples; groups in 2–4 samples on the scaled-1000 and compat indexes; no false-positive KO was ever grouped with a true one (`fp_grouped` empty). *Why:* the false positives left after gather are small units with a median of 1 gather-assigned k-mer (true KOs: 15), borrowing hits from a well-supported true KO. The true KO's coverage barely moves when the small unit takes or loses those hits, so their draws are not anti-correlated: correlation is the wrong signal for lopsided pairs.
  - *Absent in some draws* (lower bound 0) as a false-positive flag: catches 52% of false positives on kfp\_s100 (21% with the dense tier) but also 8% of true KOs, so only 9% of flagged KOs are false. Weak.
  - *Next (to decide):* group by shared evidence instead of correlation: link a unit to another when most of its allocated hits (across draws) lie on k-mers it shares with that unit, which catches lopsided pairs; and report each unit's share of evidence that is its own. Separately, the copies model error still keeps `abundance_zi` intervals from calibrating.

* **Phase 4, step 15 — ambiguity groups by shared evidence (replaces anti-correlation).** Over all posterior draws, each unit's allocated hits are accumulated per k-mer; a unit is linked to another holder when at least half (`shared_evidence` = 0.5) of its allocated hits lie on k-mers the two share, and linked units form a group. A unit never given hits (explained away in every draw) links to every co-holder of its hit k-mers. Profiles gain `own_evidence`: a unit's share of allocated hits on k-mers no other detected unit holds. Test: a small unit with one stray hit of its own and three k-mers shared with a well-supported unit (the lopsided pair correlation missed) is grouped with it, and so is one with no own k-mers; an independent unit is not.
  - *Simulation* (5 seeds, D = 100, `--twins 0.3`), anti-correlation (step 13) → shared evidence: twin pairs grouped dense 0.10 → **1.0**, s10 0.22 → **1.0**, floor + dense 1/5 0.21 → **1.0**, floor 0.00 → 0.17 (≈8 k-mers per unit leave little shared evidence to count); group intervals holding the true total 0.53 → 0.83 (dense), 0.88 → 0.78 (s10), 0.80 → 0.90 (floor + dense 1/5), – → 0.85 (floor); true positives in a group 3–23%. Per-unit interval coverage unchanged (0.93–0.97 overall; twins both present still 0.70–0.88).
  - *Interpretation:* grouping now finds the pairs the data cannot separate, and group totals are closer to calibrated than members, as intended; group intervals still undercover a little (0.78–0.90), for the same split-approximation reason as twins' member intervals. The 0.5 threshold is a default, not tuned; the fmh run (`fp_grouped`) will show how many real false positives it catches. Result pending on HPC.

* **Phase 4, step 16 — shared-evidence groups on real data (HPC, PR #19 code, D = 100).** Point estimates and intervals unchanged from step 14. Per index, over 10 samples (true positives TP, false positives FP, `kmers_unique` ≥ 1):

  | Index | FPs grouped | TPs grouped | groups / sample | groups of true KOs only | FP median `own_evidence` | group cover |
  | --- | --- | --- | --- | --- | --- | --- |
  | kfp\_s1000 | 1.0% | 2.0% | 27 | 99% | 1.0 | 0.58 |
  | kfp\_s100 | 5.6% | 2.3% | 53 | 91% | 1.0 | 0.57 |
  | kfp\_s100\_d10 | 5.0% | 1.3% | 34 | 87% | 0.98 | 0.57 |
  | fmh\_compat | 1.9% | 2.1% | 38 | 98% | 1.0 | 0.36 |

  Of the FPs that are grouped, 84–100% share the group with a true KO (the benchmark's `fp_grouped` reported this conditional share by mistake; see fix below). Groups have a median of 2 KOs (up to 60 on kfp\_s100).
  - *Interpretation.* The FPs left after gather are mostly not sequence-similarity confusions: their hits sit on k-mers no other detected KO holds (median `own_evidence` 1.0). Gather already removed the similarity-driven FPs in phase 3 (348 → 105 per sample), so what remains looks like independent evidence: chance or error matches, or genes the truth rule does not count (collagen VII was the most persistent in step 3). Groups mostly flag pairs of true KOs present together (paralogous KOs), which is the intended output; their totals are not calibrated (0.57 cover), for the same copies-model-error reason as the members' intervals. The per-unit calibration gap on real data is therefore the copies estimate, not ambiguity.
  - *Fix:* `fp_grouped` averaged `is_in` over FPs, but FPs without a group gave null and were skipped, so it reported the share of *grouped* FPs that sit with a true KO (0.84–1.0) instead of the share of all FPs grouped with one (≈ 1–5%). Nulls now count as not grouped; regression test added.
  - *Next options:* model the copies error (spread of member lengths per unit, at build) so KO-level intervals calibrate; examine the ungrouped FPs directly (which KOs, which k-mers, error-rate dependence) to decide whether a background false-hit term is needed; exact joint presence in the posterior split (step 13) only matters for near-identical units, which are rare in this benchmark.

* **Phase 4, step 17 — presence probability, copies error, false-positive diagnostic.** Three changes, from step 16's next options plus a new requirement: uncertainty must also reflect weak evidence, i.e. a unit hit on few k-mers, especially k-mers many units share.
  - *Presence probability (`present_prob`, always computed).* Per unit gather keeps: present vs. detected by background. Evidence is the unit's own k-mers (the tier-2 hit k-mers gather gave it), so hits other units explain do not count. An absent unit is detected by background (off-target homologs, error k-mers) with probability 1 − exp(−μ), μ = `BACKGROUND` × reads × *t\_g*, then has *h* own k-mers ~ Geometric(`CLUMP`). Each k-mer's odds are multiplied by its index *holders*, since conserved motifs shared by many units are what unindexed genes carry too. A present unit's *h* follows a distribution fitted to the sample. With A absent index units, P(present | *h*) = *w\_h* / (*w\_h* + A · P(*h* | absent)), fitted by fixed-point iteration; at the fixed point, *w\_h* ≈ units seen with *h* minus expected background. Posterior draws drop a unit with probability 1 − `present_prob`; group totals keep a dropped member's share. Defaults `BACKGROUND` = 1e-6 per read pair per unit of *t\_g* and `CLUMP` = 0.3 are calibrated on the step-16 fmh run (below).
  - *Tried and dropped (on the step-16 fmh profiles, offline):*
    1. Background fitted per sample (β per posting, presence prior ρ, per-unit hit rate *h*/*m*): **not identifiable**. A free per-unit hit rate explains any few-hit unit as well as background does, so fits went to one extreme or the other. In 6 of 10 s100 samples β ≈ 7e-8 and nothing was flagged; in 4, β ≈ 0.017 and 79–90% of all KOs were flagged. The simulation had not shown this, because its false positives are explained by other detected holders.
    2. A beta-binomial population prior on the hit rate, or a free present-unit *h* distribution with ρ free: ρ → 1 and every unit ≈ 1.
    3. Background ∝ *m\_g* (β per posting): calibratable, but within one-k-mer units it ranks the wrong way (AUC 0.47 at s100, 0.19 at s1000), because false positives are small KOs. ∝ *t\_g* ranks neutrally at fixed *t\_g* and well where *t\_g* varies (floor index: AUC 0.87 within one-k-mer units, where false positives sit on floored KOs).
    4. A Poisson count given background: overconfident at 2 own k-mers (1.00 vs 0.96 true at s100), because background comes in clumps (37% of s100 false positives have ≥ 2 own k-mers). Replaced by the geometric.
    5. (Simulation) a floor of one copy's k-mers at the fitted coverage: flagged 54% of true units, because divergent strains keep 0.85¹¹ = 17% of k-mers at 85% identity.
  - *Calibration (step-16 fmh profiles, offline, 3.3 M read pairs, holders = 1).* P(true | own k-mers = 1, 2, 3), predicted vs. observed: kfp\_s100 0.80/0.85, 0.94/0.96, 0.98/0.98; kfp\_s1000 0.99/0.99, 0.99/1.00, 1.00/1.00; kfp\_s1000\_floor8 0.98/0.95, 0.99/0.99, 0.99/0.99; fmh\_compat 0.98/0.97, 0.99/1.00, 1.00/1.00. Among single-k-mer KOs, no per-unit feature in the profile separates true from false (AUC 0.42–0.53 for hits, reads per k-mer, *m\_g*, hit fraction, coverage, gather rank), so the calibrated probability of the count is what the data support; *t\_g* adds separation only where it varies. Holders cannot be checked on KO indexes (*u\_g* / *m\_g* ≈ 0.99: KO k-mers are almost never shared); it will matter for MGnify90 clusters. False positives are spread over many KOs (694 KOs for 1,048 s100 false positives, 29% of rows from KOs false in ≥ 3 samples), except with the floor, where some are false in all 10 samples (e.g. K03044, an RNA polymerase subunit homologous to bacterial rpoB).
  - *Copies error.* The index stores `len_cv` (`len_cv_dense`): the coefficient of variation of kept k-mers per counting member; its mean is `pin_sum` (tested). Each posterior draw scales copies by exp(N(0, *s*²)), *s*² = log(1 + `len_cv`² / copies) + `copies_error`². *Offline check on the step-16 run:* the excess error in `abundance_zi`, √(residual sd² − posterior sd²), falls with copies (dense tier 0.30 → 0.19, s100 0.54 → 0.33 from the lowest to the fourth copies quintile) but levels off. Fit σ² ≈ a / copies + floor²: a ≈ 0.07 and floor ≈ 0.17 with the dense tier; a ≈ 0.25 and floor ≈ 0.30 at s100. So `len_cv` should cover the low-copies end, but a floor of ≈ 0.2 will be needed. **Decided:** `copies_error` stays 0 until the HPC run measures `len_cv` alone, then the floor is fitted.
  - *Simulation* (5 seeds, D = 100, `--twins 0.3`; step 13/15 values in brackets). With ~10⁴ reads, background is negligible and `present_prob` ≈ 1, so intervals are as before: coverage dense 0.936 (0.933), s10 0.964 (0.959), floor 0.970 (0.969), floor + dense 1/5 0.953 (0.946); twins 0.71, 0.76, 0.88, 0.86; group cover 0.81, 0.84, 0.75, 0.91. Its false positives are paralogs keeping own k-mers after gather; the ambiguity groups handle those, not presence.
  - *False-positive diagnostic (fmh benchmark).* `query --kmers` writes the tier-2 hits per (unit, k-mer) with `hits` and `holders`. A new DETECTED step writes one row per KO gather keeps (`detected.tsv`) with `tp`, counts, `present_prob`, `own_evidence`, group and, over its hit k-mers, median holders, most hits on one k-mer and `in_genome`: the share found in the six-frame translation of the sample genomes. A false positive's k-mers in the genomes are real sequence (another gene or KO); those not in them come from read errors. `--iss_mode perfect` reruns with error-free reads. SCORE adds `prob_tp`, `prob_fp`, `flag_tp`, `flag_fp` for `abundance_zi`. Result pending on HPC (default run, then `--iss_mode perfect`).
  - *Next:* from the HPC run, (1) whether `present_prob` stays calibrated in the tool (own k-mers from gather, reads counted, holders from lookups) and how its interval dropping changes `abundance_zi` coverage; (2) fit `copies_error` with `len_cv` in place; (3) what false-positive k-mers are (in genome or not, error dependence), which decides whether `BACKGROUND` should depend on read errors.

* **Phase 4, step 18 — fmh benchmark with presence, `len_cv` and error-free reads (HPC, PR #22 code, D = 100).** Reads simulated with `--iss_mode perfect`: inferred from fresh SIMULATE tasks and every hit k-mer being in the sample genomes; the run's parameters were not in the results. Detection matches step 16 (read errors) almost exactly: kfp\_s100 completeness 0.956 (0.955), purity 0.985 (0.985), 103 false positives per sample (105).
  - *False positives come from real sequence, not read errors.* Every hit k-mer of every false-positive KO is in the six-frame translation of the sample genomes (`in_genome` = 1.0 on all kfp indexes), and removing read errors removes almost no false positives. They have a median of 1 holder (the k-mers are unique among KOs), so they are sequence in the sample genomes that the truth rule does not count for that KO: another gene, an unannotated gene, or another frame. `BACKGROUND` therefore need not depend on the error rate.
  - *`present_prob` in the tool* reproduces the offline calibration (step 17): P(true | own k-mers = 1, 2, 3), predicted vs. observed, kfp\_s100 0.80/0.85, 0.94/0.96, 0.98/0.98; kfp\_s1000 0.99/0.99, 0.99/1.00, 1.00/1.00; floor8 0.98/0.95, 0.99/0.99, 0.99/0.99; fmh\_compat 0.98/0.97, 0.99/1.00, 1.00/1.00 (corrected: first read from `detected.tsv`, whose `fmh_compat` rows were shifted one column, see step 19). Mean over false positives 0.85 vs. 0.98 over true ones at s100. Within single-k-mer KOs it separates only on the floored index (AUC 0.87; 0.50–0.52 elsewhere), as expected.
  - *`abundance_zi` intervals with `len_cv` (`copies_error` = 0):* coverage of truth depth (step 16 in brackets) kfp\_s100 0.896 (0.72), kfp\_s100\_d10 0.857 (0.75), kfp\_s1000 0.900 (0.69), floor8 0.903 (0.70), fmh\_compat 0.45 (0.43, no members so no `len_cv`); median width log(hi / lo) 1.29 (0.79), 0.58 (0.42), 2.77 (1.64). Group cover 0.80 (0.57), 0.69 (0.57), 0.87 (0.58), 0.82 (0.58). The gain is `len_cv`: coverage among intervals with a positive lower bound is about the same as overall (0.88 vs 0.89 at s100), so presence-dropping only widens already-wide intervals. The lowest two copies quintiles now cover 0.92–0.95; the highest still undercover, with excess error ≈ 0.15 (dense tier) to 0.3–0.45 (s100, s1000), the floor step 17 predicted.
  - **Decided:** `copies_error` = 0.2. Widening the intervals offline by 0.2 in quadrature gives kfp\_s100 0.95, kfp\_s100\_d10 0.97, kfp\_s1000 0.93, floor8 0.92 (0.15: 0.94, 0.96, 0.92, 0.91). Fitted on one intact profile per index (see fix), so the next run must confirm it over 10 samples.
  - *Fix:* PROFILE published `profile.tsv` and `kmers.parquet` under the same name, so 42 of 50 published profiles were k-mer tables. SCORE, DETECTED and the summaries read the task outputs and are unaffected. The tables now go to `profiles/` and `kmers/`.
  - *Next:* rerun the default (read-error) benchmark to confirm `copies_error` over all samples. To find what the false-positive sequence is, locate each false positive's k-mers in the sample genomes against the gene coordinates: inside a gene of another KO, a gene with no KO, or outside genes. That decides whether these are errors of the truth rule or genuine off-target homology.

* **Step 19 — reproducible posterior draws; `detected.tsv` column fix.**
  - `posterior_zi` is seeded, but its draws depended on input row order, which Polars `group_by` does not fix: two identical fmh test runs differed in `own_evidence`. It now sorts its hits by unit, hash and read first; a test checks that shuffled input gives identical output, and two fresh test-profile runs give identical `detected.tsv` and `scores.tsv`.
  - `DETECTED` wrote only the columns a profile has, so `fmh_compat` rows (no `n_members`) had 15 fields under a 16-field header in the stacked `detected.tsv`. Their columns were shifted one place from `n_members` on, or the file failed to parse, depending on which file came first. Every column is now written, null where missing. The step-18 `fmh_compat` presence calibration was recomputed from the shifted rows and corrected; SCORE's `prob_*` metrics read the profiles directly and were right.

* **Phase 4, step 20 — fmh benchmark with `copies_error` = 0.2 and read errors (HPC, PR #24 code, iss novaseq, D = 100; described by the run's `run.json`).** Detection and point estimates as in step 16.
  - *`abundance_zi` intervals are calibrated on real data.* Coverage of truth depth over 10 samples (step 18, `len_cv` only, in brackets): kfp\_s100 **0.950** (sd 0.005; 0.90), kfp\_s100\_d10 0.962 (0.86), kfp\_s1000 0.926 (0.90), floor8 0.926 (0.90), as the offline fit predicted (0.95, 0.97, 0.93, 0.92). By copies quintile, kfp\_s100 covers 0.945–0.956 and kfp\_s100\_d10 0.956–0.971 in every quintile; kfp\_s1000 and floor8 fall to 0.90 in the top two quintiles, so the sparse indexes could take a slightly larger floor. Median width log(hi / lo) 1.50 at s100 (1.29) and 0.95 with the dense tier (0.58). Group cover 0.93, 0.93, 0.88, 0.87 (0.80, 0.69, 0.87, 0.82). fmh\_compat stays at 0.50: its sketches carry no members, so no `len_cv`, and its copies are not meaningful.
  - *Read errors cause few false positives.* 96% of s100 false positives have every hit k-mer in the sample genomes (mean share 0.977); about 4% involve a k-mer from a read error, consistent with 105 false positives per sample here vs. 103 with error-free reads (step 18). 14% of true KOs have some error-derived hit, which does no harm.
  - *`present_prob`* calibration is the same as with error-free reads (kfp\_s100 0.80/0.85, 0.94/0.96, 0.98/0.98 at 1/2/3 own k-mers).
  - The run predates #25 and #26: 8 floor8 DETECTED tasks were again killed at 8 GB before passing on retry, and `detected.tsv` has the column shift (fmh\_compat file first here, so the other indexes' rows carry an extra field). Both are fixed on `main`.
  - **Decided:** `copies_error` = 0.2 stays. With the phase-4 gate's "calibrated intervals on simulations" met and intervals now calibrated on real data at the recommended setting (kfp\_s100), the remaining phase-4 item is where false-positive sequence comes from (step 18, *Next*).

### Phase 5

* **Phase 5, step 1 — other tools in the fmh benchmark (`--tools`).** The fmh benchmark now runs DIAMOND, fmh-funprofiler, kMermaid and HUMAnN on the same 10 metagenomes and scores them with the same truth and metrics as our indexes, so every tool sits in one `summary.tsv`. Each tool's output becomes a profile of `name`, `evidence` (its detection count) and `abundance` (`bench.py tool-profile`); SCORE treats (`evidence`, `abundance`) as one more rule. Choices:
  - *DIAMOND 2.2.8* blastx, best hit per read (`-k 1`), default sensitivity, against the same KEGG proteins our KO indexes are built from: the alignment baseline without a reference handicap. Evidence: read pairs; abundance: Σ aligned / subject length, the quantity truth depth sums.
  - *fmh-funprofiler 1.1.1*, the released `funcprofiler` on its KO sketches (scaled 1000), for its cost; `fmh_compat` already reproduces its calls. `--threshold_bp` = scaled (one shared hash): its default at scaled 1000, and on the test profile (scaled 10) it gives the same 12 KOs as `fmh_compat`, where the default 1000 bp (100 hashes) gave 2.
  - *kMermaid* (commit `edcb4ed`) ships a RefSeq-cluster model with no KO labels, so it is retrained with KOs as clusters through its own `get_kmer_dict_cluster`, as its README allows (`kmermaid_kfp.py`). Clusters are capped at 50 random proteins per KO (`--kmermaid_max_members`), near the shipped model's ~55 per cluster; uncapped KEGG KOs would not fit in memory as its Python dicts. Evidence: read pairs; abundance: reads / mean member length.
  - *HUMAnN 3.9* (MetaPhlAn 4.1.1, `mpa_vJun23_CHOCOPhlAnSGB_202403`) on its own ChocoPhlAn and UniRef90 databases, UniRef90 families regrouped to KOs with its mapping; unstratified RPK. It reports no read counts, so every KO it reports has evidence 1. The sample genomes are in UniRef90, but its KO calls go through UniRef90 → KO rather than KEGG's gene → KO table, so part of any gap will be annotation, not detection.
  - *Cost.* The trace is now raw (ms, bytes) with `%cpu` and `peak_rss`; `bench.py cost` summarises wall hours, CPU hours and peak RSS per step and index or tool. Threads differ (DIAMOND 8, HUMAnN 16, the rest 1), so CPU hours are the comparison. One-off database builds (DIAMOND_DB, KMERMAID_MODEL, HUMANN_DB) are separate rows, like our INDEX.
  - *Containers.* Tool steps run in pinned biocontainers (the kMermaid image is built by `containers/build.sh`, since kMermaid has no package); `-profile docker` locally, `-profile singularity` on HPC. Our own steps still run on the host in `.venv`. HPC kit: `hpc.example.config`, `run_hpc.sh` (Slurm head job), `containers/pull.sh` (pre-pulls images into the Singularity cache).
  - *Checked* on the test profile with Docker (DIAMOND, fmh-funprofiler, kMermaid; HUMAnN only on HPC, its databases are ~45 GB): all steps run and score, and native fmh-funprofiler matches `fmh_compat`. `tests/python/test_fmh_benchmark.py` covers each tool's parsing and `cost`. Also fixed: under `-profile test` the trace was written to `results/`, since `trace.file` resolved before the profile set `outdir`.
  - *Result pending on HPC.* Expected: DIAMOND near-complete detection at a much higher CPU cost (the fmh-funprofiler paper reports 39–99x); kMermaid's k = 5 whole-read scores trade purity for speed; HUMAnN lower completeness on KEGG truth from the UniRef90 → KO mapping. Our query is still the Python prototype, so its cost is an upper bound until the phase-6 port.
  - *Next:* the HPC run; then the divergence ladder and ablations (phase 5's other items).

* **Phase 5, step 2 — HUMAnN 4 (`humann4`); HUMAnN images fixed.** HUMAnN 4 is the current development release, 4.0.0.alpha.2 (GitHub `master`, commit `e07b3a3`, July 2026). It is not on bioconda or PyPI, so `containers/humann4` installs it from GitHub. It runs as distributed: v4 alpha ChocoPhlAn and the EC-filtered UniRef90, the only protein database it distributes, so KOs without an EC number rely on the nucleotide search. MetaPhlAn 4.1.1 uses `mpa_vOct22_CHOCOPhlAnSGB_202403`, the database HUMAnN 4 checks for. KOs come from the v4 `map_ko_uniref90`, scored like `humann` (unstratified adjusted CPM).
  - *Workaround:* HUMAnN 4 alpha.2 checks the MetaPhlAn database by looking for its tag in `metaphlan --version`. MetaPhlAn 4.1.1 prints no tag there, so HUMAnN exits whenever it runs MetaPhlAn itself. HUMANN4 therefore runs MetaPhlAn first (same `-t rel_ab_w_read_stats`) and passes its profile with `--taxonomic-profile`, whose header check accepts it. Its CPU time counts in the same task.
  - *Fix, found while building it:* the HUMAnN 3.9 biocontainer resolves bowtie2 2.2.3, whose `bowtie2-build` rejects `--threads`, so step 1's HUMANN would have failed at its index build on HPC. `containers/humann3` adds bowtie2 2.5.5 to that image, and the HUMAnN 4 image builds on it. Both images ran their bundled demos to gene families, and HUMAnN 4 through `humann_regroup_table`, locally (x86 emulation). A stub run checks the pipeline wiring with both versions. The full databases (~45 GB for 3.9, ~70 GB for 4) are only downloaded on HPC.
  - *Result pending on HPC.* Expected: HUMAnN 4 lower completeness than 3.9 on KOs without an EC number (translated search limited to EC-annotated families). If so, a second run pairing HUMAnN 4 with the 3.9 full UniRef90 (UniRef 2019 on both) would separate the database effect from the algorithm.

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

- **MGnify scale.** A floor on all 1.66×10^9 clusters would need \~160 GB of index; on the 0.45×10^9 non-singletons, \~45 GB. Deferred to phase 6 (see Progress log): the Python prototype cannot build the full release on one node, so full-scale storage and compute are measured and tuned there.
- **Component size.** Promiscuous k-mers can chain clusters into one giant component, which makes the EM serial. Measure component sizes on the development subset in phase 2 and tune the N-clusters cut-off.

* **Marginal novelty.** Setting sourmash to scaled = 100 may recover most of the completeness gap at modest cost. Run that baseline in phase 3 before building phase 4.
* **Query density set by the smallest units.** If *n\_min* forces *t\_max* near 1, the query does almost no sparsification. Measure the distribution of *t\_g* on the real database early.
* **Divergence.** Exact amino-acid k-mers miss distant homologs regardless of sampling. If the ladder shows steep decay below 80% identity, reduced alphabets or spaced seeds become mandatory, not optional.
* **Shared k-mers and hierarchy.** EM at protein-cluster level, then aggregation to function, is likely better than EM directly on functions. Untested.
* **Normalisation.** Which single-copy marker set, and whether to report per-genome copies by default.
* **Frame filter at high GC.** Keeps \~3 frames at 70% GC; acceptable, but check false positives there specifically.
* **Future work:** 30% families by mapping MGnify90 representatives onto the 128.7 M MGnify30-C2 representatives, as a coarser level for floors, EM partitions and annotation.
* **Open:** whether KO/eggNOG labels are worth the annotation run, or Pfam suffices.

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
