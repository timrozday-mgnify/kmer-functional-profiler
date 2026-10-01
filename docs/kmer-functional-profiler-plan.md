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

The whole-frame filter has a cost at gene ends. A read covering a gene's stop codon has a `*` in its coding frame and is dropped, though most of it is coding; so is a read crossing the start, whenever the upstream sequence has an in-frame stop (\~80% within 100 bp of intergenic sequence). With 150 bp reads and k = 11, coverage of the last \~40 residues ramps from full to zero, and the first \~40 lose part of theirs. That is \~10–12% of the k-mer coverage of a 330 aa protein and \~30% of a 100 aa one: a length bias that falls on short clusters, which the floor exists for. Splitting at stops (`--frames all`) recovers it but hashes \~4 frames instead of \~1.5 and keeps off-frame segments. A minimum segment length separates them poorly: off-frame stop-free runs average \~21 codons, so \~38% are ≥ 20 aa. A middle mode, *stop-free + edges*, keeps stop-free frames whole and, in frames with stops, only the terminal segments (read edge to first or last stop) of ≥ *m* aa: a read crossing a gene boundary always has its coding part as a terminal segment, and an internal segment between two stops is almost never coding in a 150 bp read. Compared in the phase-7 ablations (see Progress log, project-wide, 2026-10-01).

### Read QC

The query reads raw FASTQ, and most of what read QC does is already covered or harmless for exact protein k-mers. A separate QC pass (fastp or similar) costs one more decompress, recompress and copy of the reads. What each QC step does here:

- *Ns:* a codon containing a non-ACGT base translates to `X`, and k-mers never span `X`. Already handled.
- *Adapters:* translated adapter 6-mers are masked in the index (phase 2), so adapter k-mers do not hit. Read-through into an adapter usually puts a stop in the coding frame and drops the read under the stop-free filter; `--frames edges` keeps the insert part as a terminal segment.
- *Low-quality bases:* an error replaces the ≤ k k-mers spanning its codon with k-mers that almost never hit. That thins λ about evenly, which the model absorbs. Errors also cost more than that, though: about 4% of substitutions in a sense codon create a stop (23 of 549), which drops the read's coding frame under the stop-free filter. At 0.5% error per base, that is about 3% of 150 bp coding reads, and more on reads with poor tails. Error k-mers are also sampled and looked up for nothing.
- *Poly-G tails (two-colour chemistry), low complexity:* translate to Gly/Pro runs. Such k-mers sit in many clusters, so the promiscuity cut (64 groups) should have removed them; to be checked by counting hits on poly-G reads.
- *Host and PhiX removal, deduplication:* these are not quality steps. Host removal stays outside the tool, and deduplication is a separate ablation (Q6).

**Quality mask.** `--min-qual q` replaces each base with Phred < *q* by `N` before translation. It costs one byte comparison per base on a record that has already been parsed. It reuses the `X` rule above, so a masked base removes only the k-mers that span its codon, can never create a stop, and saves the lookups of its error k-mers. Selection is still a function of the k-mer alone, so sketch consistency holds; masking thins λ, as errors do. FASTA input is unaffected. With a mask, the remaining reason for external QC is adapter read-through, which `--frames edges` covers.

Rejected:

- *Trimming low-quality tails:* a special case of the mask, and it also discards good bases after one bad base.
- *Per-read expected-error filter:* drops whole reads, and so the good k-mers in them.
- *Weighting k-mers by base accuracy:* gives non-integer counts that break the count model.

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

1. Read FASTQ, raw or QC'd; optionally mask bases below Phred *q* as `N` (`--min-qual`, planned; see Read QC); for each read and mate, translate the frames that survive the stop filter (default), those plus the terminal segments of frames with stops (`--frames edges`, planned; see Translation and frame detection), or all six (`--frames all`); split at stops; drop segments shorter than k.
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

Pfam totals do not depend on how hits are split between clusters that carry the same Pfam label, so most errors in the within-component split do not reach the unstratified Pfam profile. They do reach the function × taxon table (Genome mode), where sibling clusters of one Pfam usually belong to different taxa.

**Evaluation:** all four rules run on the same hit table, so implement gather, winner-take-all and uniqueness-first as cheap baselines in phase 4 and compare against EM in the phase-7 ablations.

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

## Query cost

Querying at the densest per-unit rate (*t\_max* = *t\_cap* = 0.2) is affordable if lookups are sorted and merged instead of random, and if the query keeps per-(unit, k-mer) counts instead of per-read rows. The density raises the size of the sample sketch, not the cost of hashing, which is paid for every k-mer at any rate.

**Why the query is dense.** A read cannot tell which unit it came from, so it must test every k-mer that any unit could have kept: hash ≤ max *t\_g*. 43% of non-singleton clusters reach *t\_cap* (phase 6, step 3), so *t\_max* = 0.2: 200× fmh-funprofiler's rate (scaled 1000), and 200× what the 73% singletons need. No context-free sampler avoids this; score-based selection keeps the index small but not the query.

**Rough cost per 10 Gbp sample** (6.7×10⁷ reads of 150 bp, ~1.5 stop-free frames, ~60 11-mers per read, full index at the defaults; to be measured):

| Stage | Scales with | Estimate | Now (Python prototype) |
| --- | --- | --- | --- |
| Decompress, translate, hash | reads × frames | ~4×10⁹ hashes, ~1–2 min on 12 cores (kernel ~54 M reads/min); gzip decompression single-threaded unless bgzf | Rust, fine |
| Sampled k-mers | × *t\_max* | ~8×10⁸ at 0.2; distinct perhaps half (error k-mers and low coverage are mostly distinct) | streamed per batch |
| Lookup in tier 2 | sampled k-mers, and index residency | ~31 GB tier 2 + unit thresholds. Resident: 2–4 cache misses per lookup, ~1 min on 12 cores. Not resident (cold mmap on a network file system): a page fault per lookup, hours | numpy `searchsorted` on 8 B per key materialised (~27 GB) |
| Hit rows | hits × set size | hit fraction unknown on MGnify (10–50%?): 10⁸–4×10⁸ rows | one polars row per (unit, hash, read), ~25 B each: 3–10 GB |
| Gather | hit (unit, k-mer) pairs | 10⁶–10⁷ units, ~10⁸ pairs | Python sets and a heap: slow and tens of GB at this size |
| Model fits | detected (unit, k-mer) pairs × iterations | em, zi, zib, zip all fitted, one reported | 4 fits |
| Posterior draws (*D* = 100) | pairs × *D* × sweeps | 5–10× all the rest (fmh benchmark: 23 s → 2m17) | per-read rows needed for the read resampling |
| Dense tier | second read pass + dense table | 15× tables, ~6× query time (phase 4, step 9) | hash-major table, same lookup |

So at full scale the cost is set by three things, in order: per-read hit rows and Python model code (memory), index residency (time), and the posterior (time). Hashing and the sampling rate come after.

**How sylph keeps profiling cheap** (source read at v1.0.0):

- *Sketch the sample first.* Reads become a sample sketch (`.sylsp`): distinct sampled k-mers with counts, at c = 200. Everything downstream works on distinct k-mers, not on read occurrences, and a sketch is reused across databases.
- *Two stages (`.syl2db`, v1.0).* Stage 1 is a pooled minimal perfect hash (boomphf) over all genomes' sparse k-mers (screen c = 3000), with a multi-owner CSR from k-mer to genomes: loaded in full, probed once per distinct sample k-mer, and permissive (85% ANI). Stage 2 stores each genome's dense k-mers as its own Golomb-Rice block at a known offset; only screen survivors are decoded, by positional reads (`pread`) rather than mmap, and decoded blocks are not cached, so RSS stays bounded by the survivors. The screen's k-mer floor is scaled by c / screen c so small genomes are not dropped before the dense stage.
- *Cuckoo filters* are used for read deduplication (a scalable cuckoo filter over (k-mer, read-pair markers), approximate at a false-positive rate of 10⁻⁴), to save memory against an exact hash set. It is an accuracy feature for PCR duplicates, not a lookup structure.

**What transfers, and what does not.**

- *Sketch-first transfers directly*, with one change: the posterior resamples reads, so pass 1 keeps counts per distinct hash only and read ids come back in an optional second pass limited to the detected units' k-mers (as the dense tier's second pass already is). Pass 1 memory is then the sketch, not the hits.
- *Sorted merge instead of random lookup.* Tier 2 is hash-ordered, so a batch of sampled hashes sorted by hash (radix sort) is resolved by one forward walk through the tier-2 arrays: every page read once, in order, whether or not the index is resident. Radix-partitioning the sampled hashes by the build's pack hash ranges (spilled to local disk when a memory budget is set) makes the query an external merge-join: one sequential pass over the index per sample, with memory = one range of the sketch plus one range of tier 2. With the index in the page cache (node-local disk or `/dev/shm`), the same walk is just cache-friendly.
- *Cohorts.* The merge-join takes the k-way merge of many samples' sorted sketches, so one index pass serves a batch of samples. On HPC this is the cheapest mode for large studies; separately, jobs on one node share one page-cached copy of an mmapped index.
- *Two-stage screening transfers only partly.* Sylph's screen works because genomes are long: at c = 3000 a 2 Mbp genome still has ~700 screen k-mers. A floored 90% cluster has ~8 k-mers at *t\_g* up to 0.2, so any sparser screen cannot see it on its own; screening at a coarser unit (connected component, Pfam) would pass a rare short cluster only when relatives are present too. That gives up exactly the low-abundance completeness the floor exists for, so a screen is an opt-in fast mode for ablation, not the default.
- *Sylph's stage-2 layout is the right one for the dense tier.* The dense tier is only read for detected units, so it should be unit-major like `.syl2db`: per-unit blocks of sorted hash fingerprints (Golomb-Rice or Elias-Fano, with *p\_in*), read with `pread` for detected units only, loaded into a small hash map that pass 2's reads probe at the dense rate. Its size then costs disk, not RAM, and its query cost scales with detected units, not the table: the reason the dense tier is off by default (15× tables, 6× time) mostly goes away.
- *Cuckoo or Bloom prefilter.* Tier 2 is already a quotient-style filter (implicit bucket bits plus fingerprints), so a separate cuckoo filter only helps if misses dominate and the index is not resident. The blocked Bloom filter from the partitioned build (`_core.bloom_contains`) can be tried as an in-RAM prefilter (~4 GB at 10 bits per key) in that regime; with sorted merges it is unnecessary. Read deduplication with a cuckoo filter is an accuracy option to ablate, not a cost saving.

**Low memory (a 12 GB laptop).** The partitioned merge-join is the laptop mode. `--memory` sets the partition count, and the cost of a small budget is scratch disk and one sequential index read, not many passes. There is one code path: a node with the index in RAM just runs it with one partition.

- *Shard passes vs disk lookups with caching.* Random lookups against an index on disk lose. After sketch dedup, a 10 Gbp sample leaves ~4×10⁸ distinct lookups, spread uniformly over ~50 GB of index. A page cache of ~8 GB then hits ~15% of them, and an NVMe drive at 10⁵–5×10⁵ random 4 KB reads/s needs 15 min to 1 h+ (macOS and thread count dependent). Caching does little beyond what dedup already does: repeated k-mers of high-coverage genes are merged in the sketch, and what remains has no reuse. Reading the index shard by shard, in order, reads every page once: 50 GB at 2–3 GB/s is under a minute.
- *Spill the sketch, not the reads.* Pass 1 radix-partitions sampled hashes into *R* scratch files by the index's pack hash ranges: 8 B per sampled k-mer, ~6.4 GB for 10 Gbp (16 B, ~13 GB, when read ids are kept for the posterior). Each partition is then sorted, counted and merged with its slice of tier 2, which is mmapped and read sequentially (`MADV_SEQUENTIAL`, pages dropped after), so it sits in the page cache and not in RSS. With *R* = 8, a partition is under 1 GB. Re-reading the FASTQ once per index shard instead needs no scratch space, but costs a decompress-and-hash pass per shard (~2–3 min each on a laptop); keep that only as the fallback when scratch disk is short.
- *Everything after the lookup must also stream.*
  - Hits come out per hash range but the model works per unit, so per-(unit, hash) counts are spilled by unit range.
  - The `max_hash_g` check moves after aggregation and becomes a sorted join with the unit table. That table is mmapped by column and read in unit order, never randomly.
  - Gather, EM and the posterior run one component at a time, reading only that component's units. Peak memory is then set by the largest component, which Q1 has to measure: promiscuous k-mers are cut at 64 clusters, but a giant component is possible.
- *The dense tier needs no second read pass* when *t\_dense* ≤ *t\_max*, because the spilled sketch already holds every dense-rate hash. Detected units' dense blocks are read, re-partitioned by hash range and merged with the same sketch partitions. A second read pass is needed only for *t\_dense* > *t\_max* or for read ids.
- *Rough laptop budget, 10 Gbp:*
  - hashing ~2–3 min (8 cores);
  - spill write and read ~6–13 GB (seconds on NVMe);
  - sorting ~10–20 s;
  - index read ~50 GB, ~0.5–2 min;
  - model: to be measured.

  That is about 1.5–2× the in-RAM time, for ≤ 8 GB RSS. Needed on disk: the index (~50 GB after phase 6's free cuts) plus scratch.
- *Smaller indexes are the other lever.* These trade accuracy for size, unlike the merge, which does not:
  - biome subsets (e.g. human gut), which may fit in RAM whole;
  - a lighter laptop index (*t\_cap* 0.1, *n\_min* 4, singletons out).

**Rules.**

- The shipped query has no per-read state in pass 1; read ids exist only in pass 2, for detected units' k-mers, and only when the posterior or the dense tier needs them.
- One estimator is fitted by default (the one phase 7 picks); the others and the baselines (em, zib, zip, winner-take-all, uniqueness-first, `kmers_out`) are benchmark options.
- Components are independent throughout (gather, EM, posterior), so they run in parallel, and a component with one unit and no shared k-mers needs no Gibbs sweeps: nothing is split, so its interval comes from a 2-D grid over (λ, π) or the read resampling alone.
- Every query records time and peak RSS per stage, so ablations report cost beside accuracy.

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
| Index lookup (sorted hashes, fingerprints) | Rust once the layout settles (end of phase 3); memory-mapped in phase 6 | Memory and speed matter on MGnify-scale subsets. |
| Index build: clustering tables, p\_in/p\_out scoring, floors, tiers | Python (polars, DuckDB) through phase 5; Rust or a partitioned build in phase 6 | Where most design changes happen, but the Python build cannot hold the full release in memory, so it is the component that blocks full-scale work. |
| Assignment, EM, zero-inflated negative-binomial model, dispersion test | Python (numpy, scipy) | Statistical design still open; easiest to iterate and inspect. |
| Evaluation, plots, ablations | Python + Nextflow | Stays Python permanently. |
| Final CLI and full pipeline | Rust (phase 8) | Release product. |

**Rules for the split**

- The Rust kernels expose batch APIs that return numpy arrays (e.g. read index, frame, hash per kept k-mer), never per-k-mer Python calls.
- Every Rust kernel has a slow pure-Python twin under `reference/`, used only in tests; property tests check they agree.
- Before porting (end of phase 7), write a short algorithm spec: parameters, formulas, file formats. The port implements the spec; the Python outputs on fixtures become golden files.
- A Python component moves to Rust early only if profiling shows it blocks experiments, and its interface has not changed for a phase.

## Implementation plan

Ten phases, each with a go/no-go gate, plus two optional ones; phases 1–5 are the Python prototype with Rust kernels, phase 6 makes the method run on all of MGnify Proteins at a reasonable cost, phase 7 runs the ablations against that full-scale method and freezes it, phase 8 ports the rest to Rust, and phase 10 (desirable, not essential) lets users add their own proteins to a released index, and phase 11 (desirable, not essential) predicts which genomes are present from the function profile. A separate companion tool (see Genome-informed unit presence) feeds those genomes back as a prior on which units are present. Phase 3 is the first point where the prototype must at least match fmh-funprofiler, or the project should stop.

| Phase | Deliverable | Language | Gate |
| --- | --- | --- | --- |
| 0. Skeleton | Mixed repo (uv + maturin + Cargo workspace), CI, pre-commit, stub Python CLI. Fixtures, the parity job, coverage and `bench.yml` move to phase 1, when there are kernels to test. | Both | CI green on Linux and macOS arm64 for both languages. |
| 1. Rust kernels (done) | PyO3 module: FASTQ streaming, codon tables (11, 4), six-frame translation, stop-filter frames, reduced alphabets, amino-acid k-mer packing, hashing, FracMinHash filter; batch numpy outputs. Pure-Python reference twins. | Rust + Python tests | Property tests pass (frame symmetry, threshold nesting, synonymous invariance); Rust matches reference; ≥ 1 M reads/min/thread. |
| 2. Index prototype (built; gate needs a real subset) | Build from a MGnify biome subset (DuckDB/Parquet): member k-mers per 90% group, p\_in/p\_out scores, per-cluster floor (non-singletons), connected components, tier 1/tier 2, fingerprints; stored as Parquet + numpy. | Python | Sizes match the analytical estimates; scoring behaves on hand-checked clusters. |
| 3. Query + naive counts (done; gate passed, see Progress log) | Hit counting and containment; `--sourmash-compat` using the sourmash Python API. Lookup moves to Rust once the layout settles. | Python (+ Rust lookup) | Matches fmh-funprofiler containment in compat mode; ≥ parity in completeness/purity at equal density. **Stop here if not.** |
| 4. Model (in progress: baselines, EM, zero-inflated EM adopted, dense tier, *p\_in*-weighted presence, copies-scaled abundance, posterior intervals, ambiguity groups, presence probability, copies error calibrated on real data; simulation benchmark; see Progress log) | Uniqueness-weighted detection, EM per connected component, zero-inflated negative-binomial model, dispersion flag, genome normalisation, dense tier 2. | Python | Clear completeness gain for low-abundance groups over phase 3 at ≤ 2x index size; calibrated intervals on simulations. |
| 5. Tool benchmarks (in progress: DIAMOND, fmh-funprofiler, kMermaid, HUMAnN 3.9 and 4 scored on KO; cost and Pfam pending; see Progress log) | Benchmarks vs fmh-funprofiler, HUMAnN, kMermaid, DIAMOND on the fmh benchmark, with CPU time and peak memory per tool. Pfam truth and scoring on the same metagenomes (the primary benchmark, see Benchmark labels); KO scoring kept for the direct comparison with fmh-funprofiler. | Python + Nextflow | Every tool scored with the same truth and metrics, on Pfam where it can report Pfam and on KO otherwise; our accuracy and cost relative to each recorded. |
| 6. Full-scale method (in progress: nested builds done and cost model checked; free index cuts, Rust pass-1 kernel and partitioned build (tier 2 packed by hash range, sets deduplicated by set hash) done; whole-release statistics pending on HPC; query cost planned (see Query cost), Q1 ladder measured to 40 M pairs (step 17) and re-run after the memory steps (step 31: 207 s, 5.8 GB anonymous; next: unit-column readahead, read prefetch, `presence`, posterior on real data; ladder to 200 M reads on ERR7738575 and pooled runs (step 32); EM and posterior on giant components: exact methods first (step 33); several inputs per query process (step 35); step 32's run to 100 M pairs: 287 s, 6.3 GB anonymous, EM at `max_iter` in every cell (step 36); E1–E3 and A1's measurement in code, step 32 to re-run with the pooled sample and posterior cells (step 37)); Q3's Rust tier-2 lookup built; unit table memory-mapped by column; Q2 done (`coverage_em` default, EM per component); see Progress log) | Make the method run on all of MGnify Proteins at a reasonable compute cost, before any ablation. Full-scale cost study: nested all-biome MGnify subsets (1 in 10⁴ to 1 in 10 clusters) and a full-release statistics pass (per-cluster k-mer counts, predicted postings) to fit how storage, build and query cost scale. Query cost: per-stage query time and memory against the nested and full indexes, a Rust sorted-merge lookup, a unit-major dense tier and cohort mode (see Query cost). Scalable build: index build and lookup in Rust (ported ahead of the spec, as the rules above allow for a component that blocks experiments) or a hash-partitioned Nextflow build; memory-mapped lookup. Tune for cost: *t\_base*, *n\_min*, *t\_cap*, dense-tier rate and scope (e.g. non-singletons only), fingerprint width, unit-ID encoding; re-run the Pfam and simulation benchmarks at each candidate to choose defaults on cost vs accuracy (KO alongside, for comparison). Then the full build and a query of the fmh-benchmark metagenomes against it, scored on Pfam truth (the full MGnify index carries Pfam labels only). Expect many iterations and some accuracy given up for cost. | Rust + Python + Nextflow | Full MGnify index builds on one HPC node at a chosen cost/accuracy point (build time, index size, query memory and time recorded), with the accuracy given up versus phases 4–5 recorded. |
| 7. Ablations and freeze | Ablations with the phase-6 method as the baseline, so each measures a change against what will ship: EM vs gather/winner-take-all/uniqueness-first, zero inflation, *p\_in* weighting, dense tier, floors, alphabet and k (scored also on within-component split error and ambiguity-group size, which the function × taxon table depends on; see Risks), frame mode (stop-free, stop-free + edges, all six), read QC (raw, quality mask, fastp); divergence ladder (with containment-AAI calibration, see Sequence similarity); host spike-in ladder with the joint query (`--extra-index`) and human mask sidecar (`mask`) built for it (Additional references, steps 1–2), deciding whether the release ships the human mask. Where an ablation changes the index, it is run on a nested subset whose accuracy phase 6 has tied to the full build. Algorithm spec written; golden outputs recorded. | Python + Nextflow | Defaults chosen on Pfam metrics at full-scale cost, including host handling; spec reviewed and covering joint queries and masks. If the divergence ladder forces a different alphabet or k, phase 6's cost study is repeated for it. |
| 8. Rust port | Query, model and CLI in Rust, implementing the spec (the index build and lookup already ported in phase 6 are brought in line with it). Differential tests against the Python golden outputs (exact for counts, tolerance for EM). | Rust | All golden tests pass; ≥ 10x Python end to end; full-scale results of phase 7 reproduced. |
| 9. Release | Rust binary via cargo-dist, bioconda recipe, Nextflow module for the hybrid profiling pipeline; optional Python wheel of the bindings. | Rust | Tagged release reproduces phase 7 results. |
| 10. Additional references (desirable, not essential) | Decoy role for contaminant proteomes; a Nextflow recipe for study indexes (contigs or MAGs → pyrodigal → linclust 90% → `index`) queried jointly with the base; `extend` overlays and `compact` only if joint queries prove too approximate. Joint query and host mask come earlier, in phase 7. See Additional references below. | Rust (+ Nextflow module) | Base + study index recovers ≥ 90% of the completeness gain of a rebuild that includes the study's proteins, with no change to unrelated base units; extending with one study (~10⁶ proteins) takes minutes on one node. |
| 11. Genome mode (desirable, not essential) | `annotate-genomes`: a set of reference genomes (protein FASTA; nucleotide through pyrodigal) streamed through the query's first pass with read id = genome, giving each genome's raw hits per unit. `genomes`: genome abundances from a profile's per-unit hits by gather and weighted, zero-inflated EM over those contents; unexplained fraction; the per-sample function × taxon table (Pfam × genome, rolled up to species and genus, with an unclassified remainder), a primary output. See Genome mode below. | Python, then Rust (+ Nextflow module) | On the fmh benchmark's 64 genomes plus distractors, genome detection and abundance within 0.05 (F1, L1) of sylph on the same genomes; function × taxon abundance L1 and (taxon, function) F1 no worse than HUMAnN 3.9's stratified output on KO; 10⁵ genomes annotated in ≤ 1 day on one node; the genome fit ≤ 10% of the query's time. |

Out of scope initially: long reads (indels break frames; would need FragGeneScan-style frameshift handling), eukaryotic genes, metatranscriptomes.

## Progress log

What each step did, and the choices, results and interpretations behind it, newest phase last. Add an entry with every step; keep superseded choices and say what replaced them.

### Project-wide

* **Decided:** MGnify90 clusters are the only grouping level; no 30% families.
* **Decided:** package name `kmer_functional_profiler` (tool name may still change before release); licence GPL-3.0-or-later, so FragGeneScanRs can be linked.
* **Full-scale runs:** anything over the whole release (subset extraction, index build) ships as a Nextflow pipeline with README and setup scripts for HPC; local work uses samples only.

* **Decided (2026-09-28): full-scale check in phase 6, not phase 4.** No test or benchmark so far uses the full MGnify release: unit tests use fixtures, the simulation 300 synthetic units from 100 MGnify seed proteins, the fmh benchmark KEGG KO indexes (tier 2 72 MB, dense tier 1.09 GB at 1/10), and the only real MGnify build is the pre-fix gut 1-in-10,000 subset (213,645 proteins, about 4×10⁻⁵ of the release). The Python build holds the members table and k-mer tables in memory, so a full-release build needs the Rust port (or a partitioned build) anyway. Phase 6 therefore carries the full-scale cost study and tuning (see the phase table). Expect a lot of tuning there to bring cost down, and some compromises on accuracy; defaults chosen in phase 5 are provisional until phase 6 has priced them, and the phase-4/5 accuracy numbers are the reference that any cost saving is measured against.
* **Decided (2026-09-29): full-scale method before ablations; phases re-ordered.** Supersedes the phase-6 placement above. Ablations measured on the prototype at small scale could pick defaults that the full-scale cost tuning then overturns, so the full-scale method now comes first. Phase 5 keeps only the tool benchmarks; the new phase 6 (full-scale method) takes the cost study, scalable build, tuning and full build from the old phase 6; ablations, divergence ladder, spec and golden outputs move to phase 7 and run against the phase-6 method; the Rust port of query, model and CLI becomes phase 8 and release phase 9. The index build and lookup move to Rust (or a partitioned build) in phase 6, ahead of the spec, because they block full-scale work; phase 8 brings them in line with the spec.
* **Decided (2026-09-29): reference extension is a desirable, non-essential objective (phase 10).** Users should be able to add their own proteins (a study's MAGs and assembly functions) to the MGnify90 index. It does not gate the release, but earlier phases should not rule it out: see the constraints in Extending the reference (superseded 2026-10-01: see below and Additional references).
* **Decided (2026-09-29): Pfam is the primary benchmark label; KO is secondary.** The tool labels MGnify90 clusters with Pfam, so defaults and gates from phase 5 on are chosen on Pfam truth. The KO benchmark stays, for direct comparison with fmh-funprofiler (which profiles KOs) and with the phase-3 to phase-5 KO results. The fmh benchmark gains Pfam truth from the same genomes and Pfam scoring for our indexes, DIAMOND and HUMAnN; fmh-funprofiler and kMermaid stay on KO. Details: Benchmark labels. Earlier KO-based decisions (phases 3–4) stand until the Pfam benchmark has rerun them.
* **Decided (2026-10-01): frame mode is a phase-7 ablation, with a new *stop-free + edges* mode.** The default stop-free filter drops every read whose coding frame holds a stop, so reads crossing a gene's stop codon (and often its start) are lost whole, under-covering the \~40 residues at each gene end: a length bias against short clusters (estimate in Translation and frame detection; not yet measured). Plan:
  - *Mode:* `FrameMode::Edges` (`--frames edges`) in `DnaScanner::scan`: stop-free frames whole; in frames with stops, the segment before the first stop and the segment after the last, each kept if ≥ *m* aa (*m* a developer-only constant, default 20, swept in the benchmark). Internal segments are dropped. Reference twin in `reference/`; tests: a read crossing a stop keeps its coding terminal segment, internal segments never hash, stop-free reads hash as in `StopFree`.
  - *Benchmark:* a `frames` list parameter in the fmh-benchmark and simulation workflows (stopfree, edges at *m* = 15/20/30, all), scored side by side. Beyond the usual metrics (completeness, purity, L1, Spearman; Pfam, KO alongside): abundance error and completeness binned by protein length, and per-residue coverage near gene ends in the simulation (where read positions are known), which shows the bias directly. Cost from the query ladder: hashes and sampled k-mers per read, lookup and total time.
  - *Decision rule:* change the default from stop-free only if edges (or all) cuts the short-protein abundance error without losing purity, at a hashing cost the Q6 front accepts.
* **Decided (2026-10-01): read QC becomes optional, through a built-in quality mask; ablated in phase 7 with frame mode.** The aim is to profile raw FASTQ in one pass, without a separate QC pass that rewrites the reads. Reasoning in Read QC. Plan:
  - *Mask:* in `FastxHits::read_next`, bases whose quality byte is below 33 + *q* become `N` as they are copied into the batch (needletail already parses `qual()`). `--min-qual` (CLI and `query`; 0 = off, the default until phase 7). Reference twin: a `mask_quality` helper in `reference/` applied before `hash_dna`. Tests: a low-quality base removes exactly the k-mers that span its codon; a low-quality base that would make a stop leaves the frame stop-free; *q* = 0 and FASTA input are unchanged.
  - *Benchmark:* a `qc` list parameter in the fmh-benchmark and simulation workflows: raw, mask at *q* = 10, 20, 30, and fastp (defaults, with poly-G trimming) followed by raw. Use InSilicoSeq's novaseq model (binned qualities) and its miseq model (poor tails). Score the usual metrics (Pfam, KO alongside), plus coding reads dropped by the stop filter (simulation, where frames are known), sampled k-mers per read, and end-to-end time including fastp. On the Q1 run (ERR7738575), record the fraction of bases masked and the change in sampled k-mers and hits, for raw vs mask. Run crossed with the frame-mode arms (mask × edges), since both recover reads the stop filter drops.
  - *Decision rule:* make the mask the default at the *q* that matches fastp-then-raw on Pfam completeness, purity and abundance error. Then document external QC as optional. If fastp still wins, identify which step accounts for the gap (adapters, poly-G) before adding any of it to the tool.

* **Decided (2026-10-01): additional references by joint query and mask sidecar; overlays only if needed.** This supersedes the overlay-first design of 2026-09-29. Gather, EM and the posterior take components from the hit table, so several indexes with one hash scheme can be queried jointly with offset unit ids. A study index or a host-proteome decoy then competes with MGnify units with no change to the model. In-place updates are I/O-heavy (any added unit touches almost every tier-2 hash range, and existing clusters' *p\_in* needs their members), so the base stays immutable and extras are layers. Host reads are handled upstream (hostile) for host-rich samples, plus a mask sidecar from the six-frame translated host genome (≈10⁻⁵ of postings lost at 20 letters, k = 11), which also flags host-derived MGnify clusters. The joint query and mask move to phase 7 (they decide whether the release ships a human mask); decoy role, study recipe and conditional overlays stay in phase 10. Details: Additional references.

* **Decided (2026-10-01): genome mode is a desirable, non-essential objective (phase 11); batching by input count, not by shared lookups.** Users want to know which genomes are present, predicted from the function profile. Genomes are annotated by the query's own first pass (read id = genome), so their content is defined in the same sampled k-mer space as the reads' hits, and the genome fit reuses gather and the zero-inflated EM one level up (units play the role of k-mers, genomes the role of units). Batching: the measured per-sample fixed costs are per job (index staging 125–280 s, step 31) and per process (start-up, `Index.load`), not per lookup, so several inputs per process (phase 6, step 35) captures most of the gain; a shared, deduplicated lookup across samples (the cohort merge-join of Query cost) is built only where its deduplication is measured to pay. Details: Genome mode.

* **Decided (2026-10-01): genome-informed unit presence is a separate companion tool (working name `kfp-prior`), as Bracken is to Kraken; desirable, not essential.** Units missed in low-coverage genomes that almost always carry them should be reported as probably present, and units missed in high-coverage genomes as confident absences. Detected genomes (phase 11) set a per-unit prior from clade carriage frequencies, which each unit's own hits then update. It is imputation from taxonomy: it raises completeness but adds no information beyond the taxonomic profile, so imputed values are reported apart from observed ones. A separate tool keeps the profiler's output evidence-only, and lets the prior take taxa from any source (genome mode or sylph). It talks to the profiler only through files, and needs its own release cycle and genome-set data, not the index. Details: Genome-informed unit presence.

* **Decided (2026-10-01): the function × taxon table is a primary output of genome mode, not an optional `--stratify`.** Users need to know which taxa carry each function in each sample (HUMAnN's stratified output). Genome mode already fits genome depths to unit hits, so each unit's hits split among the genomes that carry it in proportion to the fit's responsibilities, λ\_G *c\_{G,u}* / Σ\_H λ\_H *c\_{H,u}*; summing units by Pfam and genomes by GTDB rank gives the table. Hits no genome explains go to `unclassified`; ambiguity groups of near-identical genomes are reported at their lowest common rank, not split. Phase 11 stays desirable, not essential, but its gate now includes the table. Two consequences for earlier phases: the within-component unit split now matters (sibling units of one Pfam can sit in different taxa), so phase 7 scores it; and the profile keeps per-unit hits after the EM split as well as raw. Details: Genome mode, Function × taxon table.

* **Decided (2026-10-01): containment AAI per unit, sylph-style; member-level placement only if the data call for it.** Users asked for sylph's ANI outputs in protein space, both for the matches and as a blanket over references, to place sample functions more finely than MGnify90. The plan is in Sequence similarity (containment AAI). In short: `aai` = min(1, `copies_zi`)^(1/k) on detected units, with intervals from the posterior; `aai_naive` (per-unit zero-truncated Poisson correction, no reassignment) on every hit unit, which the profile already lists, so the blanket mode costs one column; and a `component` column so neighbour AAIs can be read from the profile. Rejected: finer units (95% or member-level) and in-tool alignment. Deferred: a member sidecar for nearest-member placement, built only if calibration data show enough close, diverse-cluster hits. Calibration and its gate go with the phase-7 divergence ladder. Status: planned, no code.

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

* **Phase 5, step 3 — tool comparison on the fmh benchmark (HPC, commit `02109ec`, iss novaseq, D = 100; `results/tools-v1`, described by its `run.json`).** Our indexes reproduce step 20 (kfp\_s100 completeness 0.955, purity 0.985). KO scoring, 10 seeds, detection at `min_hits` 1 (ours: `kmers_unique`; tools: `evidence`), abundance ours `abundance_zi`, tools their own:

  | Index / tool | Predicted | Purity | Completeness | Low-25% completeness | Spearman (TP) | L1 |
  | --- | --- | --- | --- | --- | --- | --- |
  | DIAMOND | 8110 | 0.872 | 0.999 | 0.996 | 0.973 | 0.139 |
  | kfp\_s100\_d10 | 6866 | 0.985 | 0.955 | 0.837 | 0.988 | 0.137 |
  | kfp\_s100 | 6866 | 0.985 | 0.955 | 0.837 | 0.957 | 0.223 |
  | kMermaid | 8234 | 0.732 | 0.852 | 0.814 | 0.440 | 1.261 |
  | HUMAnN 4 | 5379 | 0.986 | 0.749 | 0.428 | 0.866 | 0.452 |
  | HUMAnN 3.9 | 5298 | 0.984 | 0.737 | 0.416 | 0.863 | 0.473 |
  | kfp\_s1000\_floor8 | 5101 | 0.980 | 0.706 | 0.388 | 0.827 | 0.559 |
  | fmh-funprofiler | 4987 | 0.976 | 0.688 | 0.295 | 0.707 | 0.664 |
  | fmh\_compat | 4801 | 0.991 | 0.672 | 0.268 | 0.090 | 1.357 |
  | kfp\_s1000 | 4712 | 0.996 | 0.663 | 0.263 | 0.807 | 0.574 |

  - *DIAMOND* is the accuracy ceiling for detection (completeness 0.999) but pays ~1000 false-positive KOs per sample (purity 0.872; 0.889 at `min_hits` 2), likely best-hit reads landing on paralogous KOs (not yet checked). kfp\_s100\_d10 matches its abundance accuracy (L1 0.137 vs 0.139, Spearman 0.988 vs 0.973) at much higher purity, losing 4.5 points of completeness, nearly all in the low-coverage quartile (0.837 vs 0.996), at DIAMOND's lowest threshold.
  - *DIAMOND purity/completeness curve.* A threshold on read pairs per KO (`--diamond_min_hits`, default 1,2,3,5,10,20,50,100, scored in place of `--min_hits` for DIAMOND only) gives the curve, for comparing at matched purity. Recomputed offline from this run's profiles (purity / completeness / low-25%): 1: 0.872 / 0.999 / 0.996; 2: 0.889 / 0.995 / 0.980; 3: 0.901 / 0.989 / 0.959; 5: 0.916 / 0.977 / 0.910; 10: 0.937 / 0.941 / 0.766; 20: 0.958 / 0.881 / 0.529; 50: 0.979 / 0.750 / 0.139; 100: 0.989 / 0.608 / 0.045. No threshold reaches kfp\_s100's point (0.985 / 0.955): at 10 read pairs DIAMOND is below it on both axes, and matching its purity costs DIAMOND a third of its completeness. Many false positives carry many reads, so they are not sampling noise; a per-read filter (bitscore or identity) or the multi-KO genes that `diamond` counts for every KO are the next things to check.
  - *HUMAnN* sits between scaled 1000 and scaled 100: completeness 0.74–0.75 with purity 0.985, quantification better than any scaled-1000 index. **HUMAnN 4 did not lose completeness to 3.9** (0.749 vs 0.737), so the expected EC-filter loss is not visible at KO level and the paired-database run is not needed. `min_hits` 2 rows are empty by construction (evidence is always 1).
  - *kMermaid* over-calls (purity 0.73) and quantifies poorly (Spearman 0.44), as expected from k = 5 whole-read scores over KO-level clusters.
  - *fmh-funprofiler* (native) matches `fmh_compat` detection within 0.016 completeness; its own abundance (Spearman 0.71) beats `fmh_compat`'s `abundance_zi` (0.09), since the compat index has no copies and holders to model. At equal density our `kfp_s1000` is at parity on detection and better on abundance (0.81); `kfp_s1000_floor8` is +0.02 completeness.
  - *Cost: pending.* The archive has no `trace.tsv`, so `bench.py cost` cannot run; CPU hours and peak RSS per tool need the trace from the HPC run directory.
  - *Interpretation:* at scaled 100 the prototype beats every non-alignment tool on both detection and abundance, and with the dense tier matches DIAMOND's abundance accuracy; its remaining gap to DIAMOND is low-coverage completeness. The phase-5 gate (every tool scored on the same truth) is met for KO; Pfam scoring and cost remain.

### Phase 6

* **Phase 6, step 1 — cost-study pipelines (`workflows/mgnify-subset`, `kmer_functional_profiler.cost`).** Both halves of the cost study run from the subset workflow, so the full-release statistics and the nested builds share one extraction path.
  - *Nested subsets.* `--sample` takes a list (e.g. `100,1000,10000`): MEMBERSHIP and EXTRACT run once at the densest sample, then SUBSET keeps `cluster_rep % N == 0` for each N, so every subset contains the sparser ones, and INDEX builds each under `1inN/index/`. The raw Nextflow trace (as in the fmh benchmark) records each build's time and peak RSS.
  - *Whole-release statistics.* Per-cluster distinct k-mer counts need a cluster's members together, but members are spread over the whole `protein_id` range. EXTRACT therefore tags every member with a bucket (DuckDB `hash(cluster_rep) % --buckets`) and MERGE writes one members table per bucket, so each bucket holds whole clusters: a one-off 1 TB shuffle whose output the full build will reuse. STATS runs pass 1 of the build (`index._n_kmers`, same masking, k and alphabet) per bucket for exact `n_kmers`, plus `n_members`, `n_full_length`, `sum_len`, `max_len`, and samples (hash, cluster) pairs at 1/1000 (`--stats_rate`). COMBINE counts clusters per sampled hash (`groups.tsv`, the promiscuity profile) and predicts sizes (`cost.tsv`).
  - *Cost model (`cost.predict_cost`).* Applies the build's rules to `n_members` and `n_kmers`: *t\_g* from the same expression as the build (`index._t_g`, now shared), candidates *t\_g*·*n\_kmers*·*keep*, where *keep* is the fraction of sampled pairs on k-mers in at most `max_groups` clusters; floored units keep min(candidates, *n\_min*); dense postings max(*t\_dense*, *t\_g*)·*n\_kmers*·*keep*. Bytes use the dtypes `PackedTable.build` picks, with one key per posting and `--sets` value sets per posting (both upper bounds; the build stores each (unit, *p\_in* level) set once, and the nested builds give the real `tier2_sets / postings`). Default grid: *n\_min* 4/8/16 × *t\_cap* 0.05/0.2 × *t\_dense* 0/0.02/0.1, for 1-in-1 to 1-in-10⁴.
  - *Checked:* on random clusters `cluster_stats` gives the build's `n_kmers` exactly, and the predicted postings, dense postings and tier-2 bytes (with the build's set ratio) are within 5% of the built index (`tests/python/test_cost.py`); the workflow's test profile runs nested builds, statistics and COMBINE end to end. Also changed: MEMBERSHIP writes membership sorted by `protein_id` (so each EXTRACT range reads only its row groups instead of the whole table) and checks `--max_protein_id` itself; EXTRACT checks that every member in its range has a sequence, replacing MERGE's global check, which a bucketed MERGE cannot do cheaply.
  - *Expected cost:* pass 1 runs at ~16 M residues/s on a laptop (Rust hashing plus polars distinct), so the 1.07×10¹² residues take ~20 h of one such node, spread over 256 STATS jobs. At full scale the unit id needs 35 bits with the 4 *p\_in* bits, so each set value becomes a u64, and ~4×10⁹ postings put set ids near the u32 limit: the byte estimate will show how much unit-ID encoding (the phase-6 tuning item) has to recover.
  - *Fix after the first HPC run:* MEMBERSHIP on the whole release was killed at 96 GB (third attempt). DuckDB's `memory_limit` was the whole allocation, `--biome root` still joined against all 1.66×10⁹ cluster ids, and the output was sorted globally (5.7×10⁹ rows). Now DuckDB gets 75% of the task's memory, `root` skips the join, and membership is written partitioned by EXTRACT's `protein_id` ranges (`membership/shard=i/`) instead of sorted, so each range still reads only its own files.
  - *Result pending on HPC* (commands in the workflow README). *Next:* check the predictions for the nested subsets against their builds, fit build time and memory against subset size, then choose between the Rust port of the build and a hash-partitioned build.

* **Phase 6, step 2 — nested all-biome builds (`cost-nested`), cost model checked.** The nested run (`--biome root --sample 100,1000,10000`, default build parameters: *t\_base* 0.001, *n\_min* 8, *t\_cap* 0.2, no dense tier, 16-bit fingerprints) completed on HPC. The whole-release statistics run (`cost-full`) is still to come.
  - *Builds* (INDEX rows of the trace; 8 CPUs requested, 3.5–5 used):

    | Subset | Proteins | Residues | Units | Postings | Tier 1 + 2 | Index dir | Build | Peak RSS |
    | --- | --- | --- | --- | --- | --- | --- | --- | --- |
    | 1 in 10⁴ | 0.53 M | 9.8×10⁷ | 0.15 M | 0.34 M | 4.5 MB | 12 MB | 7 s | 4.9 GB |
    | 1 in 10³ | 5.35 M | 1.0×10⁹ | 1.54 M | 3.38 M | 48 MB | 126 MB | 72 s | 6.5 GB |
    | 1 in 10² | 53.2 M | 9.9×10⁹ | 15.4 M | 33.8 M | 514 MB | 1.3 GB | 16.4 min | 39.2 GB |

    Units, postings and table bytes scale linearly (×10.0 per step); ~2.2 postings and ~33 bytes of tier 1 + 2 per unit, tier 2 ~9–10 bytes per posting. Promiscuous candidates dropped (*max\_groups* 64) grow with density: 0, 0.04%, 0.47% of candidates. `tier2_sets / postings` is 0.239, 0.239, 0.243. `postings.parquet` and `units.parquet` make up about 60% of the index directory. The other steps were small: MEMBERSHIP 3.3 min at 3.5 GB, EXTRACT ≤ 73 s per shard, MERGE 6 s at 20 GB.
  - *Build scaling.* Time is ~75 ns per residue at 1 in 10⁴ and 10³ but 99 ns at 1 in 10²; memory is ~4.5 GB of fixed cost plus ~3.5 bytes per residue. Extrapolated to the full release (1.07×10¹² residues, 100× the 1-in-100 subset): ≥ 27 h and ~3.7 TB for one process. The single-process Python build cannot do the full release on one node; the full build has to be partitioned (by k-mer hash, with a global union-find over postings for components) whatever language it is in.
  - *Fix: uint32 overflow in residue sums.* polars' `str.len_bytes()` is UInt32, so `n_residues` in `meta.json` wrapped at 1 in 100 (1.33×10⁹ reported; 9.92×10⁹ true, from 1.33×10⁹ + 2·2³²), and so did the cumulative sum `_batches` splits on: batch ids restarted after each wrap, merging batches into ~2.3× the intended 2×10⁷ residues. Units were never split, so the index is correct, but batch memory was higher than set, which may account for part of the 1-in-100 superlinearity. At a full-release bucket (~4×10⁹ residues) the same would have merged batches too. Both now cast to UInt64 before summing, as does `cluster_stats`' `sum_len` (which cast after).
  - *Cost model against the builds.* STATS and COMBINE run locally on the 1-in-1000 members (47 s, 4.5 GB), predicting 1 in 10³ and 10⁴ with the build's parameters and `--sets 0.24`: postings +0.57% / +0.52%, tier-1 hashes +0.56% / +0.19%, tier-2 bytes +0.63% / +0.52%. Tier-1 bytes were over-predicted (+49% / +99%): the model assumed one tier-1 set per unit, but only units with a posting have one (519 K of 1.54 M), and at 1 in 10⁴ the set ids then fit in u16. `predict_cost` now counts units with a posting as Σ(1 − e^(−*m*)) (Poisson on the expected postings *m*; 525 K predicted vs 519 K built, 52.5 K vs 52.5 K), which puts tier-1 bytes at +0.65% / +0.08%; `tests/python/test_cost.py` now checks tier-1 bytes too. The model is good enough to size the full index from `cost-full` without building it.
  - *Interpretation.* Linear extrapolation gives ~1.5×10⁹ units, ~3.4×10⁹ postings and ~50 GB of tier 1 + 2 for the full release before the u64 value width (35-bit unit ids) and growing promiscuity are counted; `cost-full` will give the real figures. Postings stay under the u32 limit (4.3×10⁹) but not by much, so offsets and set ids need a width check at full scale.
  - *Next:* the `cost-full` run (now with the overflow fix, which matters for its per-bucket batching); then design the hash-partitioned build, and profile where the 1-in-100 build spends its time and memory to decide how much of it moves to Rust.

* **Phase 6, step 3 — feasibility of size and speed reductions (analysis only, no code changed).** Measured on the 1-in-1000 build and its members (local), the cost model (COMBINE grid on the 1-in-1000 statistics, ×1000 for the full release), and the simulation (5 seeds, `gather_zi`). Full-release figures are linear extrapolations until `cost-full` reports.
  - *Where the size is.* At full scale the query-resident index would be ~290 GB, and the tiers are the smaller part. The unit table loads at 154 B per unit (1.54×10⁹ units → ~240 GB): 16 columns, of which `pin_hist` (16 × u32) is 64 B, `name` duplicates `cluster_rep`, and `t_g`/`max_hash_g` are 16 B derivable from `n_kmers`. 66% of units have no posting (90% of singletons), so they can never be hit. Tier 2 is ~31 GB (~34 GB with u64 values), tier 1 ~17 GB. On disk, `postings.parquet` and `units.parquet` add ~75 GB.
  - *Free cuts (no change to query results).* (1) Tier 1 and the build's components are not read by the query, which recomputes components on the hit subgraph (`query.py`); `u_g` is unused too. Dropping them saves ~17 GB and removes the only global graph step (union-find over ~3.4×10⁹ edges) from a partitioned build. Tier 1 comes back with the two-stage lookup, if the Rust query needs it. (2) Keep unit rows only for units with a posting, renumbered (5.3×10⁸ at full scale, 29 bits) and slimmed (`cluster_rep`, `n_kmers`, singleton flag, `m_g`, `pin_sum` f32, `len_cv` f16, `pin_hist` as u8 while `m_g` ≤ 255 (it is ≤ 10 at the defaults); dense-tier columns only with a dense tier): ~36 B per unit, ~19 GB, memory-mapped by column. (3) `postings.parquet` becomes optional (inspection only). Together: ~290 GB → ~53 GB, with no accuracy change.
  - *Precision.* Fingerprints 16 → 8 bits saves ~1 B per posting (11%) but raises ε from ~6×10⁻⁵ to ~1.6×10⁻² per lookup: ~0.3 false hits per 150-bp read at *t\_max* 0.2 (~20 lookups), several times the genuine chance-match background the plan estimates (~0.3% of off-target k-mers, ~0.06 per read). 12 bits saves only if bit-packed (~0.5 B, 5%). Neither is worth it. *p\_in* bits: re-packing the 1-in-1000 postings with 4/3/2/1/0 bits gives 808/783/725/575/522 K sets and 30.9/30.7/30.3/29.1/28.6 MB, so dropping *p\_in* saves ≤ 7% and loses `em_pin`'s 20–30% lower scatter. What *p\_in* bits do change is the value width at full scale: renumbered units (29 bits) plus 4 bits is 33 bits (u64); 3 bits, or *p\_in* in a separate u8 array beside u32 unit values, keeps u32 (~1 B per posting, ~10%).
  - *Density.* *n\_min* sets the size: floored non-singletons hold 96% of postings (412 K × 8), singletons 3.7%. Tier 2 at full scale by *n\_min* 2/4/8/16 (*t\_cap* 0.2): 8.6/16/31/59 GB. Simulation (units only): completeness 0.49/0.73/0.85/0.94, Spearman 0.74/0.77/0.84/0.88. So halving *n\_min* costs 12 points of unit-level completeness; the Pfam-level cost should be smaller (a Pfam present is usually carried by several clusters) and needs the Pfam benchmark. *t\_base* changes only singletons: 0.0001/0.001/0.01 give 30/31/45 GB. The dense tier is the most expensive feature: +120 GB at *t\_dense* 0.02 and +250 GB at 0.1 (*n\_min* 8), against +0.05 Spearman on KO units; it stays off for the first full build. Dropping singletons altogether saves 3.7% of postings and 22% of the (slimmed) unit rows; 90% of singletons already have no posting.
  - *Query rate.* The query looks up every read k-mer with hash ≤ *t\_max* = *t\_cap*, since 43% of non-singletons (fewer than 160 distinct k-mers) reach the cap. *t\_cap* 0.1 halves lookups for 5% fewer postings but leaves 19% of floored units below *n\_min* (3.6% at 0.2; 43% at 0.05); the simulation's units are long enough that 0.1 changes nothing and 0.05 costs 1 point of purity, so short MGnify clusters need the real benchmark. Lookup cost at full scale is not measured yet: the prototype materialises 8 B per key for lookup (~27 GB at full scale), so the memory-mapped Rust lookup comes first anyway.
  - *Build time and memory.* Profile of the 1-in-1000 build (81 s, 4.9 GB locally): pass 1 (exact distinct k-mers per unit) 47 s (58%), almost all in polars `unique` over every k-mer (the Rust hashing is 10.7 s over all passes); `_len_cv` (a fourth hashing pass) 6.3 s; packing and writing 5.7 s. A Rust distinct-count kernel over unit-sorted members would roughly halve build time. Peak memory holds about three copies of the sequences (original, masked, and `partition_by` batches, with `members` kept alive for later joins), which fits 39 GB at 1 in 100 (9.9 GB of sequence); holding one copy is a small Python change.
  - *Partitioned build.* With components gone, everything except *n\_groups* (*p\_out*) is local to a unit and so to a `cost-full` cluster bucket: *n\_kmers*, *t\_g*, candidates, *p\_in*, `len_cv`, the floor's top-*n\_min*. *n\_groups* feeds the promiscuity cut and the score. Counting it from candidates only (a hash-partitioned reduce of each bucket's candidates) matches the exact count for 99.45% of postings at 1 in 1000 and keeps 354 promiscuous candidates above 64 where the exact count drops 4,603. But sharing grows about 4× per decade of density (shared postings 0.25%, 1.1%, 4.0% at 1 in 10⁴, 10³, 10²), so the approximation will be much worse at full scale. The exact alternative broadcasts a Bloom filter of candidate hashes (~10¹⁰ at full scale, ~14 GB at 10 bits each) to every bucket job, which emits (hash, unit) presence for matching k-mers only, then reduces by hash range. `cost-full`'s `groups.tsv` will show how much sharing there is to get right.
  - *Recommended order.* (1) Free cuts: drop tier 1, build components and `u_g`; slim and trim the unit table; make `postings.parquet` optional; hold one copy of the sequences. (2) Rust distinct-count kernel for pass 1. (3) Partitioned build over the `cost-full` buckets, with exact *n\_groups* via the Bloom-filter broadcast. (4) Only then trade accuracy, on the Pfam benchmark: *n\_min* 4 vs 8, *t\_cap* 0.1 vs 0.2, singletons in or out, dense tier rate. Fingerprint and *p\_in* precision stay as they are.

* **Phase 6, step 4 — free cuts (index format 2).** Item (1) of step 3's order, none of which changes what the query computes.
  - *Dropped:* tier 1, the build's connected components and `u_g`, and the `tier1_per_unit` parameter (supersedes phase 2's "tier 1 is the top 4 by score per unit, mapped to component ids"; tier 1 returns with the two-stage lookup if the Rust query needs it). The query already finds components on the hit subgraph. The cost model no longer predicts tier-1 bytes, and sizes set values by units with a posting instead of all clusters.
  - *Unit table:* only units with a posting get a row, renumbered from 0 in `cluster_rep` order (`unit_pfam.parquet` follows). Columns: `unit`, `cluster_rep`, `n_members`, `n_kmers`, `t_g`, `max_hash_g`, `m_g`, `pin_hist` (smallest unsigned type holding the largest `m_g`), `pin_sum` and `len_cv` as f32, plus the dense-tier columns with a dense tier. Dropped `name` (the query derives it from `cluster_rep`; sourmash imports keep theirs), `n_full_length`, `n_counting`, `n_candidates`, `n_promiscuous` (still summed into `promiscuous_dropped`). Stats gain `n_clusters`; `n_units` now counts indexed units. The query reads format-1 indexes too (it ignores `tier1_per_unit`).
  - *Other:* `postings.parquet` is written only with `--postings` (`build_index(..., postings_parquet=True)`); the members table drops its sequences once batched, so the batches hold the only copy.
  - *Result, 1 in 1000 (local):* 525,764 of 1,543,842 units indexed (the cost model predicted 525 K); unit table 238 MB → 34 MB in memory (154 → 64 B per unit); tier 1 (16.8 MB) gone; tier 2 unchanged byte for byte (30.9 MB, same postings and sets); index directory 126 MB → 39 MB; build 81 s / 4.9 GB → 76 s / 4.4 GB peak RSS. Extrapolated ×1000: ~290 GB query-resident → ~31 GB tier 2 + ~34 GB units.
  - *Not done:* deriving `t_g`/`max_hash_g` from `n_kmers` and `n_members` at load (16 B per unit; sourmash imports have no `n_kmers`, so it waits for the Rust loader), f16 `len_cv` (polars has no f16), and memory-mapping the unit table by column (with the Rust lookup). Together they would reach step 3's ~36 B per unit.
  - *Next:* (2) the Rust distinct-count kernel for pass 1.

* **Phase 6, step 5 — Rust distinct-count kernel for pass 1 (`_core.distinct_kmers`).** Item (2) of step 3's order. The kernel takes a batch's sequences and their unit ids (non-decreasing, which `load_members` and `_batches` guarantee, and which it checks), hashes every k-mer of each unit's run of members into one buffer, sorts and dedups it, and returns the count per unit. Runs are split over threads (`std::thread::scope`, slices of about equal residues; no new dependency). It replaces the polars `unique` over every (member, hash) in `index._n_kmers`, which `cost.cluster_stats` (the `cost-full` STATS step) also calls.
  - *Checked:* identical `n_kmers` for all 1,543,842 units at 1 in 1000 (the kernel also returns 0 for the 850 units without a k-mer, which callers filled with 0 before); rebuilt index identical to step 4's (units table, stats and every `.npy`); Rust test on hand-built runs, including unsorted ids and a length mismatch.
  - *Result, 1 in 1000 (local, 12 cores):* pass 1 46.9 s → 1.4 s (~700 M residues/s against ~21 M); whole build 76 s → 30 s, peak RSS unchanged (4.4 GB). STATS on `cost-full` should now be bound by reading its bucket and by the sampled-pairs pass, not by pass 1 (step 1 estimated ~20 node-hours for pass 1 over the release).
  - *Where the time is now (profile, 28 s):* `_kmers` 15.8 s, the passes that keep (unit, member, hash) rows (passes 2 and 3 at *t\_max* 0.2, and `_len_cv`): `hash_proteins` 6.5 s (single-threaded) and polars `unique`/joins the rest; `_len_cv` 5.1 s of it; packing and writing 2.5 s. Peak memory is still set by holding the sequences and the pass-2/3 tables, which the partitioned build (item 3) splits by bucket.
  - *Next:* (3) the partitioned build over the `cost-full` buckets, with exact *n\_groups* via the Bloom-filter broadcast. The HPC environment needs the rebuilt extension before `cost-full` runs.

* **Phase 6, step 6 — partitioned build (`kmer_functional_profiler.partition`, `mgnify-subset --build`).** Item (3) of step 3's order. The build runs over the MERGE buckets (whole clusters each), with units keyed by `cluster_rep` until the last step: CANDIDATES (per bucket: pass 1, unit table, candidate hashes), BLOOM (once: Bloom filter of all candidate hashes, 7 probes by double hashing, `--bloom_bits` per hash; plus 1,025 quantiles of the candidate hashes), PRESENCE (per bucket: every distinct (unit, k-mer) up to the largest candidate hash that passes the filter, with its counting members and whether it is the unit's own candidate, split into `--ranges` hash ranges), GROUPS (per range: exact *n\_groups*, then candidate rows only, split by bucket), POSTINGS (per bucket: score, promiscuity cut, floor, `len_cv`, Pfam) and PACK (once: renumber and write). `build_index` was split into shared helpers (`_prepare`, `_select_postings`, `_unit_pfam`, `finish_index`) that both paths use, so the rules exist once.
  - *Choices.* *n\_groups* is exact, not the candidates-only approximation step 3 rejected: the filter only thins the rows each bucket emits, and GROUPS drops its false hits because they have no candidate row. Hash ranges are cut at candidate-hash quantiles rather than at equal widths: candidates sit below each unit's own *t\_g*, so equal widths put most of them in the lowest range (seen in the test: 1,130 real rows per range against several times that in range 0 before). PRESENCE stops at the largest candidate hash, since no k-mer above it can be a candidate. Stage outputs are new files, never rewritten, because Nextflow stages inputs as links into upstream work directories. No dense tier yet (it stays off for the first full build).
  - *Checked:* on 300 random clusters in 3 buckets at *k* 4 (so k-mers are shared across buckets and the promiscuity cut drops some), with a 2-bit filter so ~40% of presence rows are false hits, the partitioned index equals `build_index`'s file for file (units, postings, unit Pfam, every tier-2 array, stats); the same holds at 1 in 1000 (8 buckets, 8 ranges, 10-bit filter) against step 5's index, and for the workflow's `-profile test,test_build` run (two buckets) against `build_index` on the same members.
  - *Result, 1 in 1000 (local, stages run one after another):* CANDIDATES 11.3 s over 8 buckets (≤ 1.5 s each), BLOOM 0.8 s, PRESENCE 19.1 s (≤ 2.7 s each), GROUPS 0.8 s, POSTINGS 6.4 s, PACK 2.3 s; 42 s in all against 30 s single-process, as every per-bucket stage re-reads, masks and batches its members. 10.8 M candidate keys (counted per bucket), 11.7 M presence rows (1.08 per candidate), 2.6% of them false hits (more in the upper ranges, where fewer k-mers are candidates), ~11.5 B per row on disk, ranges within 10% of each other.
  - *Full-scale expectation (×1000):* Bloom filter ~14 GB, read by every PRESENCE job; ~1.2×10¹⁰ presence rows (~140 GB) over `--ranges` 256, ~47 M rows per GROUPS job; per-bucket steps as STATS plus one pass at *t\_max*. BLOOM's probes are numpy (`bitwise_or.at` over ~7×10¹⁰ positions, likely an hour or more); a Rust kernel if it matters. PACK still holds every posting (~3.4×10⁹) and the string-keyed set dedup of `PackedTable.build` in one process: it is the remaining single-node step, and the next to split (pack tier 2 by hash range, whose keys are already contiguous, with set ids made global afterwards).
  - *Next:* split PACK; then run `--build` on the whole release (it needs `--pfam true`, unlike `cost-full`).

* **Phase 6, step 7 — tier 2 packed by hash range (PACK split into UNITS, PACK\_RANGE, CONCAT).** PACK held every posting and ran `PackedTable.build` over them in one process. Now POSTINGS also computes each unit's index columns (`m_g`, `pin_hist`, `pin_sum`, f32 `len_cv`; `write_index`'s `unit_columns`) and writes its postings sorted by hash. UNITS (once) concatenates the unit tables of units with a posting, numbers them in `cluster_rep` order, casts `pin_hist` to the width the largest `m_g` needs, and fixes the tier-2 layout from the total postings and *t\_max*; pack ranges are the candidate quantiles rounded down to key boundaries, so hashes that share a key (and so pool their value sets) never fall in two ranges. PACK\_RANGE (per range) reads that range's postings from every bucket (row-group pushdown on the hash-sorted files), maps them to unit ids and builds a `PackedPart`: sorted keys, a local set per key, and the range's distinct sets with a 128-bit content hash each. CONCAT (once) deduplicates the sets globally, concatenates fingerprints, set ids and offsets bucket counts, and writes `meta.json` with `write_meta`.
  - *Changed in `PackedTable.build`* (both builds, so they stay identical): value sets are numbered in order of their content hash (two splitmix64 sums over the set's values, 128 bits; a collision would merge two sets, ~10⁻²⁰ at 10⁹ sets) instead of by string rank, which also replaces the string-keyed dedup step 2's ponytail note flagged; and the offsets bucket count comes from the number of (hash, value) rows instead of distinct hashes, because the partitioned layout is fixed before any range has seen its hashes. Postings and distinct hashes differ by 0.2% at 1 in 1000, the same bit width, so the layout did not change: the rebuilt 1-in-1000 tier 2 has the same offsets, fingerprints, bytes (30,920,774) and value set per key as step 5's; only set numbering differs. The cost model already sized buckets by postings. Format stays 2 (the layout fields in `meta.json` mean the same).
  - *Checked:* `PackedTable` built in three ranges with 4-bit fingerprints (many pooled keys, shared sets) and concatenated equals the table built at once, dtypes included; the random-cluster test, the 1-in-1000 run (8 buckets, 8 ranges) and `-profile test,test_build` all equal `build_index` again (`candidates_expected` to float rounding).
  - *Result, 1 in 1000:* UNITS 0.1 s, PACK\_RANGE 1.2 s over 8 ranges, CONCAT 0.7 s (the old PACK took 2.3 s). The 8 ranges hold 2.04 M local sets against 0.81 M distinct ones (2.5×), and 262 K to 750 K postings each: ranges are cut at candidate quantiles, and postings sit lower in hash than candidates.
  - *Full-scale limits:* no step holds every posting now, but CONCAT holds tier 2 (~31 GB) and every range's sets before dedup (more duplicated at 256 ranges than at 8; up to ~2×10⁹ at 16 B of hash each plus contents), and UNITS the ~34 GB unit table. If CONCAT outgrows a 128 GB node, dedup sets in a reduce keyed by set hash first. Pack ranges could be balanced on posting-hash quantiles recorded by POSTINGS.
  - *Next:* run `--build` on the whole release (it needs `--pfam true`, unlike `cost-full`), after `cost-full` has checked the predicted sizes.

* **Phase 6, step 8 — Bloom filter in Rust (`_core.bloom_insert`, `_core.bloom_contains`).** BLOOM's numpy probes (`bitwise_or.at` over 7 positions per key) were the slow part of that step at scale, and PRESENCE's numpy lookups the largest part of its time. Both are now Rust, with a blocked filter: each key's 7 bits lie in one 64-byte block (block from a remix of the hash by multiply-shift, bits from a second remix, 9 bits each), so an insert or lookup touches one cache line instead of up to 7, which matters for a ~14 GB filter that no cache holds. Inserts are threaded without atomics: each thread owns a run of blocks and scans every hash, setting only its own bits. Lookups split the hashes over threads. The filter length is rounded up to whole blocks. PRESENCE memory-maps the filter, so jobs on one node share its pages.
  - *Cost of blocking:* 0.97% false hits at 10 bits per key against 0.82% for the unblocked filter (10.8 M random keys); at 1 in 1000, 3.1% of presence rows are false hits (was 2.6%), all dropped by GROUPS as before.
  - *Checked:* Rust tests (no false negatives, false-hit rate between 0.5% and 2%, threaded insert equal to a sequential one on a block count that does not divide by the threads, partial blocks rejected); the Python filter test; the 1-in-1000 partitioned build and `-profile test,test_build` still equal `build_index`.
  - *Result:* 10.8 M keys insert in 0.03 s against 0.88 s for numpy (29×) and look up in 0.03 s against 0.48 s (16×); 10⁸ keys insert in 0.42 s and look up in 0.25 s (12 cores). At 1 in 1000 BLOOM takes 0.2 s (was 1.0 s; the rest is reading candidates and the quantiles) and PRESENCE 9.9 s (was 19.8 s). At full scale (~10¹⁰ keys) filling the filter should take about a minute, against ~14 min extrapolated for numpy; BLOOM is then bound by reading the buckets' candidate files (~80 GB).

* **Phase 6, step 9 — set dedup reduce before CONCAT (`dedup_sets`, `SetSlice`, `mgnify-subset` DEDUP).** Step 7 left CONCAT holding every pack range's value sets before deduplication (2.5× the distinct sets at 8 ranges, more at 256). Each `PackedPart` keeps its sets in content-hash order, so DEDUP job *s* of `--ranges` takes, from every part, the run of sets whose first hash word lies in [*s*·2⁶⁴/*n*, (*s*+1)·2⁶⁴/*n*), deduplicates them, and writes the distinct sets in hash order with each part's local-to-job id map. Global set ids are then a job's base (the running sum of earlier jobs' counts) plus its job id, so CONCAT only concatenates. `PackedTable.concat` takes the slices; without them it deduplicates in one slice, as `build_index` does, so there is still one code path.
  - *Checked:* `PackedTable` from parts with sets deduplicated in three set-hash ranges equals the table packed at once (dtypes included); the random-cluster test, the 1-in-1000 run (8 buckets, 8 ranges) and `-profile test,test_build` equal `build_index`.
  - *Result, 1 in 1000:* DEDUP 0.7 s over 8 jobs, CONCAT 0.2 s (was 0.6 s); each job holds 100,713–101,380 distinct sets (±0.4%, the content hashes being uniform). CONCAT now holds tier 2 (~31 GB at full scale) and one part's id map at a time; each DEDUP job holds 1/*n* of all local sets.

* **Phase 6, step 10 — query cost plan (analysis only, no code changed).** With the build partitioned, the other half of "reasonable compute" is profiling a metagenome against the full index. Design in Query cost above; sylph v1.0.0's source read for its sample sketch, two-stage `.syl2db` and cuckoo-filter dedup. All figures are estimates until step Q1 measures them; Q1 decides which of Q2–Q7 are built, and in what order.
  - *Where the cost is.* Hashing is fixed per read. At full scale the query is limited by (1) per-read hit rows in polars and the Python gather/model code, (2) random lookups into a ~31 GB tier 2 that may not be resident, (3) posterior draws (5–10× the rest). The dense rate *t\_max* = 0.2 matters through the sketch size (~8×10⁸ sampled k-mers per 10 Gbp), which sorting and spilling bound.
  - *Order.*
    - **Q1. Measure first, on HPC, with a real metagenome.** Development of Q2–Q7 is guided by these numbers, and a step is built only if Q1 shows its stage matters at full scale (see *Rule* below), to avoid optimising what does not cost. Run the current query unchanged, apart from logging, as a grid:
      - *Reads:* one public human gut metagenome of ~40 M read pairs (150 bp, ~12 Gbp; run to be chosen), plus a nested read ladder drawn from it with a fixed seed so each subset contains the smaller ones: 0.4, 1.2, 4, 12 and 40 M pairs.
      - *Indexes:* the nested all-biome builds (1 in 10⁴, 10³, 10²; 1 in 10 if built), default parameters. Lookups do not depend on the index subset and hits grow with it, so the read × index grid gives both scaling exponents and extrapolates to 40 M pairs against the full index.
      - *Settings:* `draws` 0 over the whole grid; `draws` 100 on the 1-in-100 column only. Cells that run out of memory or time are kept as results: the ladder makes sure smaller cells succeed and shows where the prototype breaks.
      - *Recorded per cell:* wall time, CPU time and peak RSS from the Nextflow trace; per-stage time and RSS from timers in `profile` (read/hash, lookup, hit aggregation, gather, each fit, presence, posterior); and the counts the design choices hinge on: sampled k-mers, distinct sampled k-mers (the value of sketch dedup), hits and hit rows (hit fraction), (unit, hash) pairs, detected units, and component sizes (largest component, which bounds the low-memory mode).
      - *Lookup at full size:* a synthetic tier 2 (random hashes, 3.4×10⁹ postings, the real layout), resident vs cold, until the full build exists.
      - *Pipeline:* a query-cost entry in `mgnify-subset` (read subsampling, then query per read subset × index), results in a `query_cost.tsv` next to the build costs; the query prediction in `cost.py` (lookups = reads × k-mers per read × frames × *t\_max*; hits by the fitted exponent) is checked against it, as the build model was.
      - *Rule:* after Q1, a later step goes ahead only if its stage is a large share (tentatively ≥ 20%) of time or peak memory at the extrapolated full-scale cell, or blocks the laptop gate; the ladder is re-run after each step, so priorities are re-read from new data rather than from this plan.
    - **Q2. Free cuts (no change to shipped estimates).** Fit only the shipped estimator by default; em, zib, zip, winner-take-all, uniqueness-first and `kmers_out` behind a benchmark flag. Pass 1 aggregates to per-(unit, hash) counts; per-read rows only when `draws` > 0, and only for detected units (second read pass, in-RAM set of their hashes). Skip Gibbs for one-unit components. Parallelise gather and EM by component.
    - **Q3. Rust lookup and aggregation.** A `_core` query kernel over the memory-mapped tier-2 arrays: sampled hashes in batches, radix-sorted and counted, resolved by one forward walk over `offsets`/`fingerprints`/`set_ids` (no materialised keys), set expansion and the `max_hash_g` check, output per-(unit, hash) counts with `pin_q` and holders as numpy. `madvise` sequential; `--memory` sets the partition count: the sketch is radix-partitioned by the pack hash ranges and spilled to scratch when it exceeds the budget, so a sample needs one pass over the index at any budget (see Low memory). Hit counts are spilled by unit range, the `max_hash_g` check becomes a sorted join with the mmapped unit table, and gather, EM and posterior stream one component at a time. Checked against the Python `unit_hits` path on the fmh fixtures (exact counts).
    - **Q4. Unit-major dense tier.** Rebuild the dense tier as per-unit blocks (sylph's stage 2): sorted hash fingerprints plus *p\_in*, Golomb-Rice or Elias-Fano, read by `pread` for detected units only, probed by pass 2. Then re-run the dense-tier comparison of phase 4, step 9 for cost (disk, RSS, time) and accuracy; the phase-4 gate's "≤ 2× index size" becomes a disk bound.
    - **Q5. Deployment and cohorts.** Index staged once per node (local SSD or `/dev/shm`), several samples per node sharing its pages; a cohort mode that merges many samples' sorted sketches in one index pass; Nextflow module for both. Record per-sample cost at cohort sizes 1, 10, 100.
    - **Q6. Cost ablations on the Pfam benchmark** (with the phase-7 accuracy ablations, same runs; KO alongside), each reporting time, peak RSS and index bytes beside Pfam completeness, purity and abundance error, to draw the cost/accuracy front:

      | Knob | Values | Cost it moves | Expected accuracy risk |
      | --- | --- | --- | --- |
      | *t\_cap* (sets *t\_max*) | 0.2, 0.1, 0.05 | sampled k-mers ∝ *t\_cap*; sketch memory | short clusters below *n\_min* (19% at 0.1, 43% at 0.05) |
      | Frames | stop-free, stop-free + edges (*m* = 15, 20, 30), all six, frame caller (one) | hashing and sampled k-mers ~1.5 : ? : 4 : 1 | high-GC reads, reads crossing stops; abundance error by protein length |
      | Posterior | *D* = 0, 20, 100; Gibbs only on multi-unit components | model time 1–10× | interval calibration (phase 4, step 20) |
      | Estimators | shipped one only vs all | model time ~3–4× | none for the shipped column |
      | Dense tier (unit-major) | off, 0.02, 0.1, full | disk; pass-2 time ∝ detected units | none; accuracy gain measured |
      | Singletons | in at *t\_base*, out | unit table −22%, postings −4% | rare singleton functions |
      | Lookup mode / `--memory` | random resident, sorted resident, partitioned spill at 64 / 12 / 8 GB, re-read reads per shard, cohort | time vs RSS vs scratch disk | none (identical counts) |
      | Two-stage screen | off, component-level sparse screen | lookups for dense stage only | loses rare short clusters: the floor's gain |
      | Read QC | raw, `--min-qual` 10 / 20 / 30, fastp then raw | no QC pass (one less read-and-rewrite of the FASTQ); fewer lookups of error k-mers | adapter read-through, poly-G tails |
      | Read dedup (cuckoo, as sylph) | off, on | sketch memory +, per-read state | abundance on duplicate-heavy libraries |
    - **Q7. Only if Q1–Q6 leave the query too costly:** the two-stage screen as an opt-in fast mode, a Bloom prefilter for non-resident random lookups, 12-bit packed fingerprints.
  - *Gate for this part of phase 6:* a 10 Gbp metagenome profiles against the full index on one HPC node in ≤ 1 h and ≤ 64 GB RSS with the index not preloaded (≤ 16 GB with the partitioned spill), with counts identical across lookup modes; and on a laptop (12 GB RAM, 8 cores, NVMe) at `--memory 8G`, peak RSS ≤ 8 GB and wall time ≤ 2× the HPC node's, with identical results. Targets to confirm after Q1.

* **Phase 6, step 11 — Q1 query-cost workflow and stage timers (`profile(timer=)`, `query --stats`, `mgnify-subset --query`).** Q1's measurement, built; the run is pending on HPC.
  - *Reads:* ERR7746321 (PRJEB49206, Hadza gut, NovaSeq 6000, ~140 bp paired, 29 Gbp, ~106 M pairs by ENA's counts; to be confirmed by the LADDER count). Deeper than the planned ~40 M pairs, so every ladder rung, 40 M included, is a proper random subset, and deeper rungs can be added later without a new sample.
  - *Stage timers.* `Timer` records wall time, CPU time and the process's peak RSS per stage: load, keys (tier 2's materialised keys, built once on first lookup), hash (the Rust read kernel, timed per batch), lookup (`unit_hits`: tier-2 lookup, set expansion, `max_hash_g`), aggregate (per-read rows and per-(unit, hash) counts), components, gather, presence, dense, fit\_em, fit\_zi, fit\_zib, fit\_zip, result, posterior, baselines, total. Peak RSS is the high-water mark at each stage's last exit (`getrusage`), so the stage that raises it shows as a step. Counts: reads, sampled k-mers, distinct sampled k-mers (estimated on 1 in 256 of hash space, so the count costs no memory), hit k-mers (sampled occurrences with ≥ 1 unit), hit rows, per-read rows, (unit, hash) pairs, hit and detected units, and the connected components of hit units linked by shared k-mers (number, largest in units and pairs). The query is otherwise unchanged.
  - *Workflow* (`--query true`): FETCH\_READS (ENA, md5-checked; or `--query_reads` local files), LADDER (nested subsets by rank in one seeded permutation, so each rung holds the smaller ones, exact sizes, mates in step), QUERY per rung × index (`draws` 0 everywhere; `--query_draws` 100 on `--query_draws_on` 1in100 only), QUERY\_COST → `query_cost.tsv` (one row per cell, stage columns `{stage}_wall_s`, `_cpu_s`, `_peak_rss`). QUERY runs at `--query_memory` (128 GB) without retries and failures are ignored, so cells that break the prototype show as FAILED in `trace.tsv` and the rest still run.
  - *Fixed:* `posterior_zi` crashed when no unit was detected (an empty `bincount` bound); it now returns an empty table.
  - *Checked:* the `test,test_query` profile (fixture reads, ladder 10 and 40 pairs, both test indexes, draws 3 on 1in1) and a stub run of the download path; the CLI test checks the stats against the profile (hit rows = total hits, hit units = rows); a test checks that ladder subsets nest and pair.
  - *Not built yet:* the synthetic full-size tier 2 for resident vs cold lookups, and the check of `cost.py`'s query prediction against `query_cost.tsv`; both come after the run.

* **Phase 6, step 12 — first Q1 run: every cell out of memory; stats written per stage.** Against `full-build`, all five rungs (0.4–40 M pairs) were OOM-killed at 128 GB within minutes, the 0.4 M rung included, so the cost does not come from the reads. `--stats` was written only on success, so the failed cells recorded nothing.
  - *Suspect:* the `keys` stage. `stored_keys` materialises tier 2's full keys as uint64 (~27 GB at ~3.4×10⁹ keys) through temporaries of the same size (bucket repeat, shift, fingerprint cast).
  - *Changed:* `Timer` rewrites the stats JSON at every stage start and end, with `running` naming the stage in progress, and logs each stage's start and end with current and peak RSS to stderr (`.command.err`), so a killed cell shows where it died. The ladder gains 10 k and 100 k rungs, where the read-side cost is negligible.
  - *Re-run (10 k to 40 M pairs):* every cell killed in 30–70 s; the 10 k cell's stats show `load` done at 43.5 GB peak RSS (20.6 s wall, the unit table read into memory) and `running: keys`.
  - *Fixed:* `stored_keys` shifts and fills the key array in place (`<<=`, `|=` with a buffered cast) instead of building the bucket, shifted, cast and result arrays side by side. Measured with `tracemalloc` on 2×10⁸ keys at ~4 keys per bucket: peak 4.8 → 2.4 GB, identical keys; ×17 for the full tier 2 (~3.4×10⁹ keys), ~82 → ~41 GB, so `load` + `keys` ≈ 85 GB, under 128 GB.
  - *Re-run:* memory fits (10 k cell: `load` 43.5 GB, `keys` 87.9 GB peak, no later stage higher), but every cell exited 1 after `gather`: `index.units["pin_hist"].to_numpy()` panics with Polars' maximum length, the full index's 16 levels per unit (~8.5×10⁹ values) exceeding 2³². `rt64` stays rejected (see the Pfam join fix), so `em_pin` takes the column as a Series and gathers the hit units' rows only.
  - *Lookup is I/O-bound:* 70.7 s wall against 3.7 s CPU for 239 K sampled k-mers (~0.3 ms per lookup): random page faults into the memory-mapped tier 2 on `/hps`. Q1's resident vs cold question, measured on real data; QUERY now reads tier 2 once into the page cache before the query (`--query_preload`, on by default; its time goes to `.command.err`), so the next run shows the resident case. It also copies the files the query reads (meta, units, tiers) to node-local scratch first (`--query_scratch`, default `$TMPDIR`; in place if the copy fails), so pages evicted by other jobs or by the query's own growth are re-read locally, not from `/hps` (which may be GPFS, whose cache the preload does not fill). `--query_in_memory` (`query --in-memory`, off by default) reads the tiers into RSS instead of memory-mapping them: no page faults, ~31 GB more peak RSS.
  - *First complete cell (10 k pairs, before the scratch copy):* 127 s wall, 98.6 GB peak RSS. `load` 49 s and 43.5 GB, `keys` 42 s and 87.9 GB, `lookup` 28 s wall against 3.4 s CPU; every per-read stage under 0.5 s. `load` and `keys` were ~3× slower than the run before, likely from the ladder's cells reading `/hps` at once. `fit_zip` and `result` raised peak RSS by 10.7 GB at 10 k reads: `em_pin`'s `pin_hist` rows, the unit columns joined into the result and the baselines' `m_g`/`t_g` all went through the full unit table.
  - *Fixed:* a `hit_units` stage filters the unit table to the hit units once (~32 K rows at 10 k pairs); `em_pin` takes those rows (unit, `pin_hist`), and the result and baselines join them. Measured on a 5×10⁷-unit table (50 chunks, 16 × u8 `pin_hist`): +0.78 → +0.02 GB RSS, identical output; ~8 GB at full scale.
  - *100 k pairs (scratch copy, preload):* 68 s wall, 111 GB peak RSS. The copy made `load` 10 s, `keys` 11 s and `lookup` CPU-bound (1.2 s wall = CPU for 2.4 M sampled k-mers); `lookup` raised RSS by 14 GB of mapped tier-2 pages (file-backed, reclaimable). `gather` and `presence` each added ~4.3 GB: both converted the full `t_g` column (~5.3×10⁸ × 8 B) to index it at hit units. The benchmark fits (`fit_em`, `fit_zib`, `fit_zip`, baselines) took 36.7 s of the 68 s, `fit_zip` alone 30.6 s (11× for 10× reads); cutting them waits until the benchmark has settled which estimator ships (Q2).
  - *Fixed:* after `hit_units`, units are renumbered 0.. in the hit units' order, so `t_g`, `m_g`, `pin_hist`, `pin_sum` and `len_cv` (and their dense forms) come from the hit units' rows; the output maps back to the index's ids. `em_pin` takes a plain array again. The only full-column conversions left are `max_hash_g`/`max_hash_dense`, which the lookup indexes per batch (Q3). Profiles on the fixtures (18 of 20 units hit, with and without a dense tier, 5 draws) are byte-identical to before.
  - *Next:* Q3.
  - *Next:* re-run the ladder. Walking the memory-mapped `offsets`/`fingerprints` without materialising keys (Q3) and memory-mapping the unit table remain the full-scale fixes.


* **Phase 6, step 13 — Q3, first part: Rust tier-2 lookup over the memory-mapped arrays (`_core.unit_hits`, `_core.packed_lookup`).** The 10 k and 100 k cells showed materialising tier 2's keys (`keys`: ~44 GB of the 88 GB peak, 11–42 s) as the largest fixed cost of a query, so the lookup moves to Rust first.
  - *Kernel* (`crates/core/src/packed.rs`): per sampled hash, the bucket from its top bits addresses `offsets` directly, the bucket's few fingerprints are scanned, and the set is expanded with the `max_hash_g` check; rows (unit, hash, read, `pin_q`, holders) in the same order as the numpy path. Arrays of any unsigned width are read in place, so a memory-mapped table stays on disk except for the buckets and sets touched. `stored_keys` and the `keys` stage are gone; `PackedTable.lookup` uses the same kernel.
  - *Checked:* equal to the replaced numpy path (materialised keys + `searchsorted`), set ids and hit rows in order, on a random table with 4-bit fingerprints (shared, merged sets), random `max_hash_g`, hashes above `max_hash`, in memory and memory-mapped; fixture profiles (with and without a dense tier, 5 draws) byte-identical.
  - *Result, 5×10⁷ keys (memory-mapped, 407 MB), 2.4 M lookups (a 100 k-pair batch):* 0.62 s and +0.52 GB RSS (the table's touched pages) against 0.85 s and +2.6 GB for the numpy path with its keys. At full scale this removes ~44 GB of anonymous memory and the `keys` stage.
  - *Not built (Q1 rule: no measured cost yet):* batch radix sort and a forward walk (the per-hash lookup is O(1) via `offsets`; sorting only helps a non-resident index, which the scratch copy and preload now avoid); aggregation to per-(unit, hash) counts in Rust (`aggregate` took 0.03 s at 100 k pairs); `madvise`, `--memory` partitioning and spill, component streaming; the `max_hash_g` check as a join with a memory-mapped unit table (the lookup still takes `max_hash_g` as a full ~4 GB array). Threads: the kernel is single-threaded (~0.26 µs per lookup here).
  - *Next:* re-run the ladder; then memory-map the unit table by column (`load`'s 43.5 GB), which also serves `max_hash_g` to the lookup without a copy.
* **Phase 6, step 14 — unit table memory-mapped by column (`UnitTable`, `units.<column>.npy`, `unit-columns`).** `load` read the whole unit table into memory (43.5 GB at full scale, 20–49 s), and the query then used ~32 K of its ~5.3×10⁸ rows plus `max_hash_g` for the lookup.
  - *Format:* beside `units.parquet` (kept, for inspection and string columns), `write_index` and `concat` write each numeric column but `unit` (the row number) as `units.<column>.npy`, fixed-size arrays (`pin_hist`) as 2-D; one column in memory at a time. `Index.units` is a `UnitTable`: `units[column]` memory-maps the column (read into memory with `--in-memory`), `rows(units)` gathers whole rows for the hit units (strings, e.g. sourmash `name`, by a filtered parquet scan), `frame()` reads the whole table for tests. Without the `.npy` files (indexes built before this step, e.g. `full-build`) a column is read from the parquet on first use, so old indexes still query, at one column's memory rather than the table's; `kmer-functional-profiler unit-columns <index>` adds the files in place. QUERY's scratch copy includes them. `max_hash_g`/`max_hash_dense` go to `_core.unit_hits` as the mapped arrays, no copy.
  - *Checked:* rows gathered from the columns equal the parquet's; fixture profiles (sparse + dense tier, 20 draws) identical with the column files and with the parquet fallback; the sourmash import (string `name`) through the scan.
  - *Result, 5×10⁷ units (full-build schema, 16 × u8 `pin_hist`; 1.8 GB parquet, 2.8 GB of `.npy`), 32 K hit units:* the unit columns the query uses plus the hit units' rows in 0.08 s against 1.6–2.1 s for reading the parquet and filtering, with no anonymous copy of the table. Peak RSS still reads ~2.1 GB (vs 3.5–5.0 GB): 32 K random rows touch most pages of each column at this size; those pages are file-backed and reclaimable, and at full scale (10× the units, same hits) a far smaller share. Writing the columns took 3.1 s. At full scale ~32 GB of `.npy` beside the parquet (60 B per unit at the current dtypes); the `load` stage's 43.5 GB should go.
  - *Baseline, 10 k pairs with Q3 and the unit table still read whole (`full-build`, scratch copy, preload):* 15.1 s wall, 53.9 GB peak RSS (was 127 s and 98.6 GB before Q3): the `keys` stage is gone. `load` 10.3 s (80 s CPU) and 43.5 GB, 81% of peak: the unit table this step removes. `lookup` 0.39 s CPU-bound for 239 K sampled k-mers (~1.6 µs each), +9.8 GB of touched tier-2 pages (file-backed). `hit_units` 1.27 s (filtering the full table; ~0.1 s expected from the columns). `fit_zip` 2.1 s is the largest stage after `load` (14% of wall, ~40% once `load` goes): `em_pin` iterates all 17.9 K detected units jointly until the slowest converges, with an `np.add.at` per iteration. Every other stage under 0.3 s. 32 K hit units, 17.9 K detected, largest component 58 units.
  - *Fixed, first run against `full-build` (no `.npy` files yet):* `hit_units` panicked with Polars' maximum length: the parquet fallback read `pin_hist` whole (16 × ~5.3×10⁸ values > 2³²), as `unit-columns` would have. Columns are now read in slices of 2²⁶ rows (`COLUMN_SLICE`), into a preallocated array (fallback) or straight into the `.npy` via `open_memmap` (`unit-columns`, so one slice in memory). Checked: slices of 3 rows on the fixture index give the same arrays both ways; 5×10⁷ units in 12 slices, identical `pin_hist`, 9.4 s against 3.9 s unsliced (~8 slices at full scale).
  - *Next:* `unit-columns` on `full-build`, then re-run the ladder. If `fit_zip` stays a large share: `bincount` over `col * 16 + level` instead of `np.add.at` (same output), or Q2's fit-only-the-shipped-estimator.

* **Phase 6, step 15 — 100 k pairs with memory-mapped unit columns; `em_pin` without `np.add.at`.**
  - *Result, 100 k pairs against `full-build` (scratch copy, preload):* 58.5 s wall, 65.7 GB peak RSS. `load` 0.02 s and 0.19 GB (was 10.3 s and 43.5 GB). `lookup` 1.46 s for 2.4 M sampled k-mers (~0.6 µs each), RSS to 36.6 GB: nearly all of tier 2's mapped pages touched (file-backed). `hit_units` 14.7 s (82.6 s CPU) and +29 GB: the unit columns read whole from the parquet (~60 B × 5.3×10⁸ units), i.e. the queried index had no `units.*.npy`, so the parquet fallback ran; `unit-columns` on the index removes it. `fit_zip` 30.6 s (52% of wall), then `fit_zi` 3.1, `baselines` 2.7, `fit_zib` 2.2, `gather` 1.3, `fit_em` 1.1 s. 292 K hit units, 164 K detected, largest component 71 units, 527 K (unit, hash) pairs.
  - *`fit_zip`:* `em_pin` runs all detected units jointly until the slowest converges, so it hits `max_iter` (1000); `np.add.at` into the (unit, level) table was 37% of it. Replaced by a `bincount` over `unit * 16 + level`: identical output; on a synthetic table at this size (164 K units, 330 K pairs) 18.1 → 12.4 s. Fewer iterations is not a fix: at 50/100/200/1000 iterations, 18.7 K/9.0 K/4.7 K/530 units are > 0.1% off the 5000-iteration coverage.
  - *Next:* `unit-columns` on `full-build`, re-run. Then `em_pin` per component (components are independent; the largest is 71 units): iterate only components not yet converged, which leaves the per-iteration cost with the slow units only. Stopping becomes per component, so estimates change within `tol`; check against the joint fit on the fixtures and the fmh benchmark. `fit_zi`/`em` share the pattern.
  - *`em_pin` per component:* components share no k-mer, so each is fitted on its own and stops at its own convergence (max change ≤ `tol` × the component's largest coverage, scale within `tol`); an iteration only touches components still moving, compacted when half of the current set has converged. Components come from the same unit–k-mer graph as `components()`. Synthetic data shaped like the 100 k cell (164 K units, 470 K pairs, 20% of k-mers shared within groups of 4; 113 K components): 13.7–14.1 s joint → 0.75–0.95 s. Against a 20 000-iteration joint fit, units off by > 0.1% (coverage / presence): joint at `tol` 1e-6 194 / 204; per component at 1e-6 254 / 315, 1e-7 200 / 208, 1e-8 194 / 204 (0.89 s), 1e-9 194 / 204. **Chosen:** `tol` 1e-8, same accuracy as the joint fit. The remaining ~200 units are slow EM components that stop at `max_iter`; 3000 iterations leaves 5 / 7 at 1.1 s, so raising `max_iter` is now cheap if accuracy calls for it. A test checks that fitting two components together equals fitting each alone (to 1e-12).

* **Phase 6, step 16 — Q2: `coverage_em` by default, benchmark estimators behind `--all-estimators`; EM per component.**
  - **Decided:** `coverage_em` (plain EM over the units gather keeps) is the default estimate. Only it is fitted by default; `coverage_zi` (with `copies_zi`/`abundance_zi`), `_zib`, `_zip`, winner-take-all and uniqueness-first need `query --all-estimators` (`profile(all_estimators=True)`). `draws` > 0 still fits `_zi`, whose intervals the posterior gives; `present_prob` (0.03 s) stays on. The fmh and sim benchmarks pass the flag; `mgnify-subset` QUERY does not, so Q1 cells now measure the default path. The Pfam benchmark's "shipped configuration" (Σ `coverage_zi` over a Pfam's clusters, phase 5) needs the flag or a switch to `coverage_em`.
  - *Saving at 100 k pairs (from step 15's stages):* `fit_zi` 3.1 + `fit_zib` 2.2 + `fit_zip` 30.6 + `baselines` 2.7 = 38.6 of 58.5 s. With the unit columns in place (`hit_units` < 1 s), the default cell should be ~6 s: `lookup` 1.5, `gather` 1.3, `fit_em` 1.1 s, the rest small.
  - *EM per component:* `em` (plain, zero-inflated, with prior) now runs on the same per-component driver as `em_pin` (`_fit_components`: components compacted as they converge, each stopping on its own), with `bincount` sums in place of the sparse products. Pairs are sorted by (unit, hash) first, so results do not depend on the order polars hands them over (the sparse product did not, and the fixture tests caught the difference). Synthetic 100 k-shaped data (164 K units, 470 K pairs), units > 0.1% off a converged fit (coverage / presence): plain joint `tol` 1e-6 1.82 s, 54 / 0; per component 1e-6 0.16 s, 130 / 0; 1e-7 0.18 s, 50 / 0; 1e-8 0.19 s, 50 / 0. Zero-inflated: joint 3.67 s, 26 / 46; per component 1e-6 0.26 s, 219 / 390; 1e-7 0.30 s, 31 / 50; 1e-8 0.34 s, 26 / 46. **Chosen:** `tol` 1e-8, as for `em_pin`: the joint fit's accuracy at ~10× less time. The posterior's per-draw zero-inflated fits gain the same.
  - *Checked:* the default profile equals the `--all-estimators` profile on its columns; fitting two components together equals fitting each alone (to 1e-12) for `em`, zero-inflated `em` and `em_pin`.
  - *Not built (Q1 rule, stage below 20% or not measured):* per-read rows only for detected units with `draws` > 0 (`aggregate` 0.03 s at 100 k); skipping Gibbs for one-unit components (`draws` 0 in Q1); `gather` by component (~1 s at 100 k, a Python heap loop: a Rust port if it grows with reads).
  - *Next:* re-run the ladder with the unit columns (`unit-columns` on `full-build`).

* **Phase 6, step 17 — Q1 ladder complete on the default path (Q2 + Q3, `full-build`, 10 k to 40 M pairs, `draws` 0, scratch copy + preload, 8 CPUs).** Every cell completed. 40 M pairs (~11 Gbp, the gate's sample size): **488.6 s query wall (570 s job, ~80 s of copy and preload), 96.7 GB peak RSS**.

  | Pairs | hash | lookup | aggregate | hit\_units | gather | presence | fit\_em | total (s) | peak RSS (GB) |
  | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
  | 10 k | 0.0 | 0.6 | 0.0 | 14.7 | 0.1 | 0.0 | 0.0 | 16.4 | 47.0 |
  | 100 k | 0.3 | 1.5 | 0.0 | 16.9 | 1.4 | 0.0 | 0.2 | 21.5 | 65.7 |
  | 400 k | 1.1 | 2.8 | 0.1 | 14.7 | 5.8 | 0.1 | 0.6 | 26.8 | 71.1 |
  | 1.2 M | 3.2 | 6.2 | 0.3 | 14.2 | 14.1 | 0.2 | 1.8 | 42.4 | 72.5 |
  | 4 M | 10.7 | 17.0 | 1.1 | 14.6 | 32.9 | 0.6 | 5.4 | 87.0 | 75.4 |
  | 12 M | 31.8 | 75.2 | 3.6 | 15.7 | 73.1 | 1.5 | 12.2 | 224.7 | 81.0 |
  | 40 M | 107.1 | 174.1 | 9.8 | 14.9 | 122.9 | 13.7 | 24.1 | 488.6 | 96.7 |

  - *Unit columns not in place:* `hit_units` is a fixed ~15 s (~83 s CPU) and +29 GB at 10 k: the parquet fallback, so `unit-columns` had not been run on the queried copy of `full-build`. Every figure below includes it.
  - *Counts per pair (40 M; exponent in reads, 4 → 40 M):* 24 sampled k-mers (1.00); distinct 4.5 (0.58), so **only 19% of sampled k-mers are distinct** and sketch dedup cuts lookups 5.4× at this depth, more with depth; hit k-mers 3.3 (1.00; hit fraction 13.8%, the low end of the 10–50% guess); hit rows 5.7 (1.00); (unit, hash) pairs 0.77 (0.48); hit units 0.30 (0.41); detected 0.19 (0.44), i.e. 7.4 M detected units.
  - *Time at 40 M:* lookup 36%, gather 25%, hash 22%, fit\_em 5%, everything else ≤ 3%. Lookup and hash scale linearly with reads and are single-threaded (CPU = wall; `fastx.rs` and `packed.rs` use one thread; job CPU/wall 1.2 on 8 CPUs). Lookup is 0.18 µs per sampled k-mer, resident. `gather` (the Python heap loop) grows as (unit, hash) pairs (0.57), 123 s. `fit_em` per component 24 s for 7.4 M detected units (0.65). `presence` jumps 1.5 → 13.7 s from 12 M to 40 M (1.37): not yet explained.
  - *Memory at 40 M (high-water steps):* lookup 48.7 GB (mostly tier 2's file-backed pages, 18 GB already at 10 k), aggregate +23 GB (224 M per-read rows in polars), hit\_units +7 GB (parquet fallback), gather +13 GB, the rest +2 GB. `getrusage` does not split file-backed from anonymous, so the anonymous share is estimated: ~60 GB, of which per-read rows ~23 GB.
  - *A giant component is forming:* largest component 58 → 110 K units (289 K pairs) across the ladder, exponent 2.1 from 4 M to 40 M, while components overall grow at 0.5. At 40 M it is still 1.5% of detected units, so per-component streaming holds, but it will grow faster than the rest with depth or cohort pooling.
  - *Against the gate* (10 Gbp, ≤ 1 h, ≤ 64 GB, index not preloaded): time passes 7× over (preloaded); memory fails at 96.7 GB, ~66 GB without the parquet fallback, of which ~35 GB is reclaimable mapped index.
  - *Q1 rule applied (≥ 20% of time or memory at the 40 M cell):* (1) `unit-columns` on `full-build` (no code). (2) Lookup: dedup sampled hashes before the lookup (sort + count per batch: 5.4× fewer probes, and sorted probes walk tier 2 in order), then threads. (3) Gather: per component, in Rust or with the EM's component driver. (4) Hash: threads (reader → workers). (5) Per-read rows only for `draws` > 0 (Q2, deferred at step 16 on time; it now qualifies on memory, ~23 GB). Also explain `presence`'s jump, and measure anonymous RSS separately (`/proc/self/status` `RssAnon`) so the gate's memory figure is not inflated by page cache.

* **Phase 6, step 18 — threaded tier-2 lookup (`_core.unit_hits`); dedup deferred.**
  - *Threads:* `unit_hits` splits a batch's hashes into contiguous runs, one per available CPU (at least 2¹⁴ hashes each), looks them up in parallel and joins the rows in order, so output is unchanged. Lookups are latency-bound random reads, so threads overlap the misses. Synthetic table (5×10⁷ keys, 5×10⁶ units, 16-bit fingerprints, in memory), 2.4 M lookups at 14% hits (a 100 k-pair batch): 0.17–0.18 s on one thread → 0.025 s on 12 CPUs (~7×), identical rows. On the 8-CPU QUERY job, expect the 40 M cell's 174 s to fall to ~30 s. Checked: a Rust test (threaded equals serial over 5 runs of hashes), the Python suite.
  - *Dedup not built:* duplicates are spread across the sample, not within a batch. The 100 k-pair rung is one `batch_reads` batch and 92% of its sampled k-mers are distinct, so per-batch dedup saves ~8% of lookups, about the cost of the sort. The 5.4× (at 40 M) needs dedup across the whole sample, i.e. sketch-first: count distinct sampled hashes over all batches, then look each up once. That needs per-(unit, hash) counts only, so it goes with (5), per-read rows only for `draws` > 0.
  - *Next:* re-run the ladder (with `unit-columns`); then sketch-first dedup with (5).

* **Phase 6, step 19 — gather in Rust (`_core.gather`, `crates/core/src/gather.rs`).** Step 17's (3): `gather` was 123 s at 40 M pairs (25% of wall), a Python heap loop over per-unit sets.
  - *Kernel:* the same lazy-heap greedy, over dense k-mer ids (pairs sorted by hash once) and a per-unit CSR of untaken k-mers compacted in place, with a flat taken bitmap. Heap order is the Python one (highest score, then lowest unit; a stale top is pushed back only if strictly below the next), so ranks match on ties too. Single-threaded; not split by component, since the Python cost was the constant factor, not the algorithm, and a global heap keeps the ranks without a merge.
  - *Checked:* Rust tests (explaining away, the `t_g` weighting, ties, repeated pairs, bad input); synthetic pairs with ~30% shared k-mers and three `t_g` values (many ties), 3 K and 1 M pairs: frames identical to the Python version's; the Python suite.
  - *Result:* 1 M pairs (272 K units taken) 1.38 → 0.06 s (23×); 31 M pairs and 12 M units, the 40 M cell's size, 5.4 s including the conversion from polars (was 123 s on HPC).

* **Phase 6, step 20 — threaded read hashing (`FastxHits`).** Step 17's (4): `hash` was 107 s at 40 M pairs (22% of wall), on one thread.
  - *Split measured first:* 1 M synthetic pairs (140 bp, gzip level 6): the hash stage 2.80 s, `gzip -dc` of both files 0.29 s. Translation and hashing are ~90%, decompression ~10%.
  - *Change:* a batch's records are read on one thread (decompress, parse, copy into one buffer), then contiguous runs of reads (≥ 4096 each) are scanned on all CPUs, each with its own `DnaScanner`, and the hits joined in order: output unchanged at any thread count.
  - *Checked:* a Python test with 12 K random pairs, batches split over up to 4 threads, equal to `hash_dna` over the interleaved mates (read, mate, frame, hash, in order); the existing FASTX tests; both suites.
  - *Result:* 2.80 → 1.06 s per 1 M pairs on 12 CPUs (2.6×). The serial read is now most of what is left. Expect the 40 M cell's 107 s to fall to ~40 s on 8 CPUs (estimate). If it stays a large share, the next cut is a reader thread that prefetches the next batch while the current one is scanned and looked up.

* **Phase 6, step 21 — per-read rows only for `draws` > 0; sketch-first dedup dropped.** Step 17's (5), and (2)'s dedup half, deferred at step 18.
  - *Dedup measured, not built:* 24 M sampled hashes (1 M pairs' worth, 19% distinct as at 40 M) against the synthetic table of step 18: threaded lookup of all 24 M 0.21 s; `np.unique` alone 3.0 s; lookup of the 4.6 M distinct 0.02 s. With the lookup threaded, a probe costs about what the dedup's own memory access does, so dedup at best breaks even in time and adds the sketch's memory (7.7 GB of sampled hashes at 40 M before dedup). It pays only when the index is not resident, where a sorted walk replaces page faults: the partitioned merge-join (laptop mode), not the HPC query.
  - *Change:* each batch's hits are reduced at once to per-(unit, hash) hits and per-unit distinct reads (read ids never span batches, so per-unit reads add up exactly), and summed into running totals (`_Summed`: batches wait until they outgrow the total, then re-aggregate, so memory stays ~2× the total). Per-read rows are kept only when `draws` > 0, for the posterior. The result's `hits`, `kmers_hit` and `reads` come from the totals. `read_rows` in `--stats` is 0 when none are kept. The dense tier's second pass still builds per-read rows (detected units only; off by default).
  - *Checked:* profiles equal the previous code's on the fixtures for sparse and dense indexes, `draws` 0 and 5, with and without `--all-estimators`, batches of 7 and 100 000 reads, merging at every batch (16 combinations, all columns); a test that summed batches equal one aggregation; both suites.
  - *Expected at 40 M pairs:* the `aggregate` stage's +23 GB of per-read rows (224 M rows) replaced by ~31 M (unit, hash) rows, ~1–2 GB with the waiting batches. To measure on the next ladder run.

* **Phase 6, step 22 — peak anonymous RSS per stage (`peak_anon`).** Step 17's (6): `peak_rss` (`getrusage`) counts the resident pages of the memory-mapped tier 2 and unit columns, which the kernel reclaims under pressure, so it overstates what the query needs against the gate's memory bound (step 17 estimated ~35 GB of the 96.7 GB as mapped index).
  - *Change:* `Timer` reads `RssAnon` from `/proc/self/status` (anonymous memory: arrays, polars frames, heap), and a sampler thread polls it every 50 ms while an outermost stage runs, as Linux keeps no anonymous peak. Each stage records `peak_anon`, the high-water mark at its end, next to `peak_rss`; `--log` prints it. Without `/proc` (macOS) the key is left out. `query_cost.tsv` gains `{stage}_peak_anon` columns without changes to the workflow.
  - *Checked:* a test with a faked `RssAnon` that spikes only between stage start and end (the sampler catches it), and the no-`/proc` case; both suites.
  - *Limit:* a spike shorter than 50 ms can be missed. `peak_rss` stays the hard bound (what Slurm's cgroup limit counts, with page cache charged but reclaimable).

* **Phase 6, step 23 — `unit-columns` without memory-mapped writes.** `unit-columns` on `full-build` (index on `/hps`) segfaulted on the HPC.
  - *Diagnosis:* `faulthandler` put the main thread in native code inside the call that fills the `open_memmap` `.npy` (no Python frame below it). The same slice reader into an in-memory array had worked in the query's fallback, reading the node-local copy, so the memory-mapped write to `/hps` is the suspect. Not confirmed by a native backtrace.
  - *Second problem:* `open_memmap` creates each file full-size and zero-filled before writing, so an interrupted run leaves columns that load as valid zeros, and the query would use them.
  - *Fixed:* `write_unit_columns` appends each slice to a plain file (`.npy` v2 header, then `tofile`) under `units.<column>.npy.partial`, renamed when complete.
  - *Checked:* the slice test (3-row slices) loads identical arrays and leaves no `.partial`; both suites. Any `units.*.npy` from the crashed run must be deleted before re-running.

* **Phase 6, step 24 — query memory for a 16 GB laptop (target < 12 GB peak), results unchanged.** Measured with a synthetic harness that fakes read hashing and the tier-2 lookup at step 17's 40 M-pair counts (24 sampled k-mers and 5.0 hit rows per pair, 25 M (unit, hash) pairs, 10.4 M hit units, small components) and runs `profile` unchanged from there, so everything after the lookup allocates as in the real query. Peaks are macOS `phys_footprint` (the harness itself holds ~0.8 GB). Each change was checked against the committed code: result frames identical, floats to 1e-12 (`present_prob` already differs between identical runs in the last bit: polars' threaded sum of `log_holders`).
  - *`peak_anon` on macOS:* `anon_rss` reads `phys_footprint` (`proc_pid_rusage`), which leaves out clean pages of memory-mapped files, so the laptop target can be checked locally.
  - *`_Summed` in 16 parts* by the first key (mod 16), each merged on its own: a merge's polars group-by covers 1/16 of the table. The read loop's peak was 8.2 GB (footprint steady at ~2 GB, with spikes at each merge); now 3.3 GB. `aggregate` is also faster (24.8 → 18.9 s).
  - *int32 indices* (`_ids`) in `em`, `em_pin`, `components` and `_fit_components`, and the sorted input freed before the fit: `em` alone on 12.6 M pairs +2.71 → +1.93 GB, same time.
  - *Posterior (`draws` > 0):* the present-fraction grid was three (detected units × 512) float64 arrays per sweep, 8 GB each at 2 M units (a 1-in-4 harness run passed 45 GB and was killed); now a chunk of units at a time, with the uniform draws taken first, at the same point in the stream. The four (draws × units) matrices (4 × 6.3 GB at *D* = 100 and 7.9 M units) are memory-mapped from an unlinked file in `TMPDIR` above `SCRATCH_BYTES` (256 MB), with group totals formed in place and quantiles taken per chunk of units. At 1 in 20 with 10 draws: posterior 257 → 215 s, peak RSS 8.7 → 3.9 GB. `test_posterior_scratch_file_changes_nothing` forces the file and one-unit chunks.
  - *Result at 40 M pairs, `draws` 0:* footprint peak 8.95 → 6.5 GB (≈ 5.7 GB without the harness), wall 216 → ~199 s. Peaks are now spread: `components` (+1.4, stats only), `gather` (+0.8), `fit_em` (+0.7).
  - *Runtime cost:* none measured; merges and the grid got faster. The one trade-off is disk: at *D* = 100 the draw matrices write ~25 GB to `TMPDIR`, written once per draw row and read once. No low-memory switch is needed for these changes.
  - *Not covered by the harness:* the real index's page cache (file-backed, reclaimable, not in the footprint, but it sets the laptop's lookup time; the partitioned merge-join of Query cost is the fix), the giant component (110 K units at 40 M), and per-read rows for the posterior (~224 M rows, ~6 GB at 40 M; next: keep reads per hash instead of per (unit, hash), since holders share them, with u32 read ids).

* **Phase 6, step 25 — posterior reads stored per k-mer, not per (unit, k-mer).** Step 24's next item. Every unit holding a k-mer sees the same reads, so the per-read rows kept for `draws` > 0 were repeated once per holder, and carried `unit`, `pin_q` and `holders` as well.
  - *Change:* each batch keeps (`hash`, `read`, `n`), with `n` = the batch's rows for that pair divided by `holders` (one row per holding unit per hit). `posterior_zi` takes the detected (unit, hash) pairs and those rows. A key's resampled count is its k-mer's (`per_hash[row]`): the same reads summed in the same order, so the draws are bit-identical. The read ids and k-mer indices are int32 where they fit. The dense pass keeps its per-unit rows for the fit and collapses them the same way for the posterior. u32 read ids were not taken up, since a sample could exceed 2³² reads.
  - *Checked:* fixture profiles at `draws` 5, sparse and dense indexes, batches of 7 and 100 000 reads, equal to the step-24 code; harness at 1 in 20 with 10 draws, equal; test helper `by_hash` builds the new inputs from the old per-read rows.
  - *Result (harness, 1 in 20, 10 draws):* read rows 9.9 M → 6.6 M (5.0 → 3.3 per pair, the hit k-mers per pair) at 20 → 16 B each; footprint before the posterior 1.58 → 0.90 GB; posterior peak 3.8 → 2.95 GB; time unchanged (211 s). At 40 M pairs, ~6.5 GB of per-read rows becomes ~2.1 GB. The rest of the posterior's peak (~2 GB at 1 in 20) is in the per-draw EM refits and the Gibbs sweeps.
  - *1 in 4, 10 draws* (the size that passed 45 GB and was killed before step 24): completes, posterior 1102 s, footprint peak 7.4 GB (3.1 GB before the posterior; 2.0 M detected units, 33 M read rows). The posterior's own ~4.3 GB grows with detected units, so at 40 M pairs with draws it is expected at ~20 GB: over the laptop target, and the next thing to cut if draws are to run on a laptop.

* **Phase 6, step 26 — `--scratch-in-memory` (workflow `--query_scratch_in_memory`).** Some HPC nodes have no usable local disk, so step 24's scratch file for the posterior's draw matrices (~25 GB at *D* = 100 on the 40 M cell) can be turned off: the matrices are then plain arrays in RAM. It is one branch in `_scratch` (`np.zeros` vs `np.memmap`); everything downstream is the same code, and the tests check both give the same frame. It is kept apart from `--in-memory` (the index tiers), since the two trade-offs are independent. `TMPDIR=/dev/shm` does the same with no flag but is capped by the node's tmpfs size and is charged as shared memory.

* **Phase 6, step 27 — fits in batches of components; cheaper component labels.** First step of dividing the post-lookup work by connected component (units linked by shared k-mers), which never interact in gather, the fits or the posterior.
  - *Labels:* `_unit_components` links each unit to the next holder of the same k-mer (a unit graph with pairs − distinct k-mers edges) instead of the bipartite unit-k-mer graph. Components are still numbered by their smallest unit, so labels are identical (checked on 200 random graphs).
  - *Batches:* `component_batches` sorts the fitted (unit, hash) rows into batches of whole components, ~`MAX_BATCH_PAIRS` (2 M) rows each, and `em`/`em_pin` run per batch, with outputs re-sorted by unit. Every fit stops per component and `max_iter` counts per component, so results are identical. The batches are taken from the pairs being fitted, since dense-tier k-mers can link units tier 2 does not. `detected` is kept as the batches' rows and `own` freed after `presence`, so no second copy is held.
  - *Checked:* `test_component_batches_change_nothing`: one component per batch, sparse and dense indexes built with each protein twice (a mutated copy as its own cluster, so components have two units), all estimators and `draws` 3, equal to one batch; batches that split components change the result (mutation check). 40 M-pair harness: identical to step 23's code.
  - *Result (40 M pairs, `draws` 0):* footprint peak 6.46 → ~5.4 GB (runs vary by ~0.4 GB), now set by `gather`, not the fits. Labels +1.4 → +0.9 GB. Time +14 s (`components` 9 → 14 s, `fit_em` 19 → 28 s: the sort by k-mer and the labelling).
  - *Next:* the posterior per component batch (per-read Poisson weights from a counter-based generator keyed by read id, one stream per component). That changes the posterior's draws once, so it waits for approval; gather per component (+0.8 GB) only if it becomes the peak on real data.

* **Phase 6, step 28 — posterior per batch of components; scratch file removed. Changes the posterior's draws (approved).** Step 27's next item.
  - *Change:* `posterior_zi` splits the detected pairs into batches of whole components (`component_batches`, sized to `POSTERIOR_BATCH_BYTES` = 512 MB at ~1.6 KB + 16 B × draws per (unit, hash) row, measured), routes each k-mer's read rows to its batch, and runs the draws per batch (`_posterior_batch`, the old body) with the draw matrices in RAM. Each batch has its own generator, seeded by (seed, the batch's smallest unit). A read's Poisson(1) weight is a SplitMix64 hash of (seed, draw, read id) through the Poisson inverse CDF (`_poisson1`), so a read hitting several batches weighs the same in each, as in one global draw.
  - *Removed:* step 24's scratch file (`_scratch`) and step 26's `--scratch-in-memory` / `--query_scratch_in_memory`: a batch's matrices are bounded, so neither disk nor a switch is needed.
  - *Grid chunks:* the present-fraction grid's chunks were 256 MB per array, ~1 GB per step whatever the batch; now 32 MB (`CHUNK_BYTES`), same time.
  - *Results:* the draws differ from step 25's, not their distribution. Harness, 1 in 20, 10 draws, old / new: estimate inside its interval 0.805 / 0.804 (coverage), 0.905 / 0.905 (abundance); median interval width over estimate 1.404 / 1.403 and 1.677 / 1.676; mean group size 1.266 / 1.267; same ambiguity group for 96.5% of units; everything outside the posterior identical. **The fmh benchmark calibration (`COPIES_ERROR`, 95% intervals at kfp\_s100) is to be re-checked on HPC.**
  - *Checked:* a component's intervals are the same run alone or batched next to one sharing its reads (one component per batch); `_poisson1` has mean and variance 1 and gives a read its own weight in any order; the fixture tests (interval coverage ≥ 0.9, seeded reproducibility, order independence, lopsided pairs) pass unchanged.
  - *Memory (harness, footprint):* 1 in 20, 10 draws: posterior peak 2.95 → 1.47 GB at the same time (213 s). 1 in 4, 10 draws: 7.4 GB (step 25), 6.0 GB with the batches before the grid fix, **5.1 GB** after it (posterior 1068 s, was 1102 s). The posterior now adds 2.3 GB over the 2.8 GB before it, mostly held for the whole posterior rather than per batch: the read rows routed to batches (a sorted copy of the 33 M (hash, read) rows next to the caller's) and the batch assignment. Routing without the copy is the next cut (step 30).

* **Phase 6, step 29 — `--low-memory` (workflow `--query_low_memory`): re-read the reads for the posterior.** For runs where memory must stay capped rather than be raised (a laptop, not HPC).
  - *Change:* with `draws` > 0 the first pass normally keeps (hash, read, n) rows of every hit unit through the query (~2 GB at 40 M pairs). With `low_memory` it keeps none, and after the fits a second pass hashes the reads again and counts, per read, the sampled hashes that are detected k-mers. A k-mer's rows over its holders are its occurrences in the read, so no lookup is needed and the rows are the first pass's for those k-mers: results identical. With a dense tier the per-read rows already come from its own pass, so the first pass no longer keeps tier-2 rows then, with or without the flag.
  - *Checked:* fixture with shared k-mers (every protein also as a mutated copy), `draws` 5, batches of 7 reads: profile equal to the default; harness, 1 in 20, 10 draws: identical output.
  - *Cost:* one more read-hashing pass (~40–100 s at 40 M pairs on 8 CPUs by step 17/20's rates; 2.8 s in the harness, whose reads are synthetic). Saving: no per-read rows through gather and the fits (0.13 GB at 1 in 20, ~2 GB expected at 40 M), and none for undetected units. The posterior still holds the detected k-mers' rows, most of them, so its own peak barely moves (1.47 → 1.43 GB at 1 in 20); the flag lowers the query's peak only where that is set before the posterior. Rows per posterior batch (a read pass per batch) would go further, at a pass per batch.

* **Phase 6, step 30 — posterior read rows routed without a copy.** Step 28 joined every (hash, read) row to its batch and sorted the result, a second copy of all read rows held through the posterior. Each batch now takes its k-mers' rows with a semi-join on its own hashes: one scan of the rows per batch, nothing the size of all rows copied.
  - *Checked:* harness, 1 in 20, 10 draws: output identical to step 28's; both suites.
  - *Result:* 1 in 20, 10 draws: posterior peak 1.47 → 1.32 GB, 213 s both. **1 in 4, 10 draws: 5.1 → 3.06 GB** (posterior 1075 s); the posterior now adds 0.3 GB over the 2.8 GB before it, so the query's peak is set before the posterior again, where `--low-memory` applies.
  - *Scaling to 40 M pairs:* per-batch memory is bounded by `POSTERIOR_BATCH_BYTES`; what grows is the per-read rows (~2 GB, or none held until the posterior with `--low-memory`) and the stages before the posterior (~5.4 GB at `draws` 0, step 27). Expected peak with draws: ~6–8 GB, under the 12 GB target (to confirm with `peak_anon` on HPC).

* **Phase 6, step 31 — Q1 ladder re-run on `main` (#55: steps 18–26; not yet steps 27–30), with `unit-columns` (HPC, `full-build`, 10 k to 40 M pairs, `draws` 0, scratch copy + preload, 8 CPUs).** Every cell completed. 40 M pairs: **206.5 s query wall (was 488.6 s, step 17), 5.76 GB peak anonymous memory (`peak_anon`, ~60 GB estimated at step 17), 71.2 GB peak RSS (96.7 GB).**

  | Pairs | hash | lookup | aggregate | hit\_units | components | gather | presence | fit\_em | total (s) | peak anon (GB) | peak RSS (GB) |
  | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
  | 10 k | 0.0 | 0.1 | 0.0 | 0.2 | 0.0 | 0.0 | 0.0 | 0.0 | 0.4 | 0.24 | 21.6 |
  | 100 k | 0.2 | 0.3 | 0.1 | 0.5 | 0.1 | 0.1 | 0.0 | 0.2 | 1.6 | 0.50 | 52.2 |
  | 400 k | 0.8 | 0.7 | 0.3 | 0.9 | 0.3 | 0.3 | 0.1 | 1.0 | 4.7 | 1.03 | 60.6 |
  | 1.2 M | 1.8 | 0.8 | 0.6 | 0.9 | 0.9 | 1.0 | 0.2 | 2.4 | 9.2 | 1.60 | 63.0 |
  | 4 M | 6.8 | 2.7 | 2.4 | 1.3 | 2.9 | 2.8 | 0.6 | 7.4 | 28.4 | 2.70 | 66.0 |
  | 12 M | 20.2 | 12.5 | 7.8 | 1.3 | 5.0 | 5.0 | 1.4 | 13.1 | 69.1 | 4.11 | 68.9 |
  | 40 M | 71.1 | 39.6 | 27.1 | 1.7 | 9.4 | 12.1 | 13.7 | 25.0 | 206.5 | 5.76 | 71.2 |

  - *Time at 40 M:* hash 34% (107 → 71 s; CPU/wall only 1.5 on 8 CPUs, so the serial read and decompression now bound it), lookup 19% (174 → 40 s; CPU/wall 4.4), aggregate 13% (9.8 → 27 s: step 21's per-batch sums replace per-read rows, the trade that removed ~23 GB), fit\_em 12%, presence 7%, gather 6% (123 → 12 s), hit\_units 1.7 s (15 s with the parquet fallback).
  - *`presence` still jumps* 1.4 → 13.7 s from 12 M to 40 M pairs (detected units only 4.6 → 7.4 M). Not explained; likely its fixed-point iteration reaching `max_iter` (500) at 40 M.
  - *Memory:* peak anonymous 5.76 GB, matching the synthetic harness at the same code (~5.7 GB without its own 0.8 GB, step 24), so the harness is representative. It passes the gate's 64 GB with room, and the 12 GB laptop target. Peak RSS is file-backed index pages: tier 2 (41 GB mapped by the end of lookup) and **hit\_units +25 GB with anonymous memory flat** (+8 GB at 10 k pairs for 32 K units): random reads of the memory-mapped unit columns pull in readahead around every unit. Reclaimable, but charged to the job's cgroup.
  - *Components:* largest 110 K units (289 K pairs) at 40 M, as at step 17; within the 2 M-pair batches of steps 27–28.
  - *Job overhead:* the scratch copy and preload take 125–280 s per job (on the shared file system's load), more than the 40 M query itself, so they set the time of every cell below 40 M.
  - *Against the gate* (10 Gbp, ≤ 1 h, ≤ 64 GB, index not preloaded): memory passes (anonymous); time passes preloaded (207 s); not-preloaded time still to measure.
  - **Next, in order:**
    1. *Unit columns read without readahead:* `madvise(MADV_RANDOM)` on the memory-mapped unit columns (and check tier 2's lookup the same way), so hit\_units maps the pages it reads, not +25 GB around them. Results unchanged; measured by peak RSS at hit\_units.
    2. *Prefetching read batches:* a reader thread decompresses and parses the next batch while the current one is hashed and looked up (step 20's next cut), targeting hash's 71 s at CPU/wall 1.5.
    3. *Explain `presence`'s jump:* record its iterations and convergence at 12 M and 40 M; if it hits `max_iter`, find why before changing `tol` or the iteration.
    4. *Re-run the ladder with steps 27–30 and `draws` > 0* on the 40 M cell: the posterior's memory and time on real data (expected ~6–8 GB peak), and the fmh benchmark's interval calibration (`COPIES_ERROR`) re-checked after step 28's change to the draws.

* **Phase 6, step 32 — the ladder to 200 M reads (ERR7738575), pooled runs beyond it, EM convergence counts.** A deep metagenome is ~200 M reads (100 M pairs), 2.5× step 31's top cell, and the EM has to hold up in more complex samples (more taxa and functions) than that.
  - *Workflow:* `--query_run` defaults to ERR7738575 (human gut, NovaSeq, 111.5 M pairs = 223 M reads, 31 Gbp; ENA's `read_count` counts reads, ~139 bp each), and `--query_ladder` gains a 100 M-pair cell. Comma-separated runs are fetched and pooled (`POOL`: R1s and R2s concatenated in the given order) into one deeper, more diverse sample, e.g. ERR7738575 + ERR7746321 (two people's gut, 217.6 M pairs = 435 M reads) with a 200 M-pair cell. Checked by stub runs: one run, two runs (with `POOL`) and local reads (no fetch).
  - *Counts:* `--stats` (and `query_cost.tsv`) gain `fit_batches`, `fit_largest_batch_pairs` (a component above `MAX_BATCH_PAIRS` is a batch of its own, so this is the largest component when it exceeds 2 M pairs), `em_iterations` (the most any component took) and `em_unconverged_units` (units of components stopped at `max_iter` = 1000, not converged) for the shipped EM. Tests: a component cut off at 3 iterations is reported; the fixture's fit converges in one batch.
  - *What depth does to the EM (from step 31):* the largest component grew 5.4 K → 110 K units from 12 M to 40 M pairs (exponent ~2.5), while components overall grow ~0.4. Extrapolated: ~1 M units and ~3 M pairs at 100 M pairs, more in pooled or more diverse samples: a giant component may be forming (percolation). EM results do not depend on how components are batched, so this is cost and convergence, not correctness:
    - *EM memory:* ~110–160 B per pair, so even a 10 M-pair component is ~1.5 GB. Safe.
    - *EM convergence:* a giant component couples many units and may converge slowly; at `max_iter` its units keep their last values silently. `em_unconverged_units` now measures this.
    - *Posterior:* ~1.6 KB per pair, and a component is never split, so a 3 M-pair component is ~5 GB in one batch, beyond `POSTERIOR_BATCH_BYTES`. The posterior is the stage depth puts at risk.
  - *Run:* `--query_run ERR7738575` (default) to 100 M pairs, `draws` 0 on `full-build`; then the pooled sample to 200 M pairs. Read `largest_component_*`, `fit_largest_batch_pairs`, `em_iterations`, `em_unconverged_units`, `peak_anon` and stage times per cell.
  - *If the giant component appears:* first find what links it. Fingerprint false hits and sequencing-error k-mers link unrelated units at random, which would grow with depth as observed; conserved motifs below the 64-cluster promiscuity cut would not. Count the linking k-mers by hits and holders before changing anything.
  - *If EM stops at `max_iter`:* accelerate the fixed point (SQUAREM), which reaches the same fixed point in fewer iterations, rather than raising `max_iter`.
  - *If the posterior's giant batch is too large:* its Gibbs sweeps need only one k-mer's holders at a time, so the per-pair arrays can be streamed in k-mer blocks inside the batch.

* **Phase 6, step 33 — plan for EM and posterior robustness on giant components (no code yet).** Step 32's run decides which of these are needed; exact methods come first, and the approximate one only if they are not enough.
  - **E1. SQUAREM acceleration of the per-component EM** (exact: the same fixed point, to `tol`). Step 32's `em_iterations` and `em_unconverged_units` say whether slow convergence is real. Accelerate the fixed-point map (Varadhan & Roland 2008: two EM steps, a squared extrapolation, a stabilising EM step, falling back to plain EM if the likelihood drops) per component, inside `_fit_components`, so converged components still drop out. Checked: coverages within `tol` of plain EM on the fixtures and the harness; fewer iterations on the largest component; `em_unconverged_units` 0 at 100 M pairs. Not raising `max_iter`: that hides slow convergence and costs time on every component.
  - **E2. Posterior streamed in k-mer blocks inside a batch.** A giant component is one batch whatever `POSTERIOR_BATCH_BYTES` says (~1.6 KB per pair, ~5 GB at 3 M pairs). A Gibbs sweep splits each k-mer's count among its holders, so it needs one k-mer's holders at a time: sweep the batch's per-pair arrays in k-mer blocks, keeping only per-unit totals across blocks. Same model and distribution; the draws change unless the random numbers are keyed per k-mer (as the read weights already are, step 28), which would make them independent of the block size.
  - **E3. Block-wise EM for components too big to fit at once** (exact at convergence). Split a component into blocks, hold the shares of k-mers shared across blocks fixed, run EM within each block, update the shared k-mers' shares, and repeat until they stop changing. Memory is per block, at the cost of outer rounds. Needed only if a component outgrows memory in the fits (~110–160 B per pair: ~1.5 GB at 10 M pairs), so after E1 and E2, and only if step 32's pooled run shows one that large.
  - **A1. Weak-link threshold (approximate; only if E1–E3 are not enough).**
    - *Strength* of the link between detected units A and B: *s*(A, B) = shared hit k-mers / min(hit k-mers of A, of B), roughly the most a cut can move either unit's coverage.
    - *Cut:* keep links with *s* ≥ τ, re-form components from them, give each k-mer spanning a cut to the component of its gather owner, and drop it from the other holders' hits and `m_g` (neither evidence nor a zero). EM and the posterior then run per component as now.
    - *The parameter:* a module constant like `COPIES_ERROR`, default 0 (exactly today's results), with a developer-only override for the benchmarks; not a user setting.
    - *Measure first:* a histogram of link strength and shared hits per link in the largest component, recorded with `--stats` only. If the giant component is held together by links of strength ~0.01 (one k-mer, one or two hits: fingerprint false hits or sequencing-error k-mers), a small τ breaks it at little cost; if by strong links (multi-domain families), cutting costs accuracy and E1–E3 are the route.
    - *Sweep (development only, as the step-18 calibrations):* τ ∈ {0, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2} as a list parameter in the fmh-benchmark and query-cost workflows (like `--diamond_min_hits`), scored side by side. Cost on the query ladder (largest component, EM iterations, posterior batch, time, memory at 100 M pairs and on the pooled sample); accuracy on the fmh benchmark against τ = 0 (L1, Spearman, completeness, purity, interval coverage; KO, then Pfam), also at depth, where weak links multiply. Choose the largest τ with negligible accuracy loss, record the evidence here, and make it the default.

* **Phase 6, step 34 — step 31's next items 1–4 in code (#59, review fixes after it); measurements pending on HPC.**
  - *1. `madvise(MADV_RANDOM)`* on every memory-mapped unit column and packed table (tier 2 and dense), through `_madvise_random`; a test checks a loaded array's base is still the `mmap`, so the hook is not skipped silently. Effect to measure: peak RSS at hit\_units, not preloaded (preloaded pages are resident whatever the advice).
  - *2. Prefetching read batches:* `FastxHits` reads and decompresses on its own thread, up to two batches ahead (a bounded channel), while the current batch is scanned. Errors arrive after the batches read before them (tested); a reader panic is re-raised when the channel closes, not taken as end of input. Effect to measure: hash time and CPU/wall at 40 M.
  - *3. `presence` convergence:* `--stats` (and `query_cost.tsv`) gain `presence_iterations` and `presence_converged` (0 at `max_iter` = 500). Tested on the fixture, including a cut-off at one iteration. To read at 12 M and 40 M before changing `tol` or the iteration.
  - *4. Draws on chosen cells:* `--query_draws_at` (default `40000000`; `''` = every cell) picks the ladder cells that also get the posterior on `--query_draws_on`. Stub-checked both ways. The run itself, the posterior's memory and time on real data and the `COPIES_ERROR` re-check, is still to do.

* **Phase 6, step 35 — plan for several inputs per query process (no code yet).** Where a batch of samples saves work, from step 31's costs:

  | Cost | Per | Batching saves it? |
  | --- | --- | --- |
  | Index staging (scratch copy, preload) | job | Yes: 125–280 s per job, more than the whole 40 M-pair query (207 s), so it sets the time of every smaller sample. |
  | Interpreter start-up, imports, `Index.load` (memory-mapping, `meta.json`) | process | Yes; seconds, which matters only for thousands of tiny inputs. |
  | Page faults on tier 2 and unit columns | node | Partly: the page cache is shared by every process on a node, so co-locating jobs already shares it. |
  | Hashing, translation | read | No. |
  | Lookup | sampled k-mer | Only through k-mers repeated across samples (one lookup per distinct hash in the batch) and a sequential index walk. Lookup is 19% of wall at 40 M pairs; a 2× cross-sample dedup would save ~10% with the index resident. With a cold index (laptop, network file system) the sequential walk is what matters: that is Query cost's cohort merge-join. |
  | Gather, EM, presence, posterior | sample | No: samples are fitted independently. Small samples can run side by side to fill the cores their small components leave idle. |

  - **B1 (build): several inputs in one process.** `query --samples SAMPLES.tsv --out-dir DIR` (columns `sample`, `r1`, optional `r2`): load the index once and loop over `profile()`, one `profile.tsv` (and `--stats` JSON) per sample, named by `sample`. Results identical to separate runs (tested). The query workflows group small samples into one job up to a total read budget (e.g. 40 M pairs per job, a parameter), so staging is paid once per job.
  - **B2 (conditional): one lookup per distinct hash across a batch.** Hash every sample of the batch with a sample tag, sort and deduplicate the sampled hashes, look each up once, and join the hits back per sample. Memory grows with the batch's sampled hashes (8–16 B each), so the batch is capped by the scratch budget, as in the laptop mode. Build it only if a measurement says it pays: the distinct-hash ratio across a batch (e.g. 10 gut samples of one study; a set of strains of one species). Expected: small for metagenomes with a resident index, large for redundant genome sets and for a cold index, where it is the cohort merge-join already planned in Query cost.
  - **Genome annotation needs neither:** many genomes stream as one input with read id = genome (Genome mode), so one process and one lookup pass serve the whole set.
  - *Measure:* per-sample wall time and peak memory for batches of 1, 10 and 100 samples of 1 M pairs, and the staging time saved per job; the distinct-hash ratio for B2's decision.

* **Phase 6, step 36 — step 32's run: ERR7738575 to 100 M pairs (HPC, `full-build`, `draws` 0, scratch copy + preload, 8 CPUs, step 34's code).** Every cell completed. **100 M pairs (200 M reads): 287 s query wall, 6.30 GB peak anonymous, 76.5 GB peak RSS.** The posterior cells did not run: `--query_draws_on` defaulted to `1in100`, which matched no index of this run (fixed in step 37).

  | Pairs | hash | lookup | aggregate | hit\_units | components | gather | presence | fit\_em | total (s) | peak anon (GB) | peak RSS (GB) | largest component (units / pairs) | EM unconverged units (of detected) | presence iterations |
  | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
  | 10 k | 0.0 | 0.1 | 0.1 | 0.2 | 0.0 | 0.0 | 0.0 | 0.0 | 0.5 | 0.27 | 24.2 | 50 / 50 | 65 of 19 K | 13 |
  | 100 k | 0.1 | 0.2 | 0.1 | 0.5 | 0.1 | 0.1 | 0.0 | 0.3 | 1.6 | 0.55 | 56.2 | 80 / 165 | 574 of 157 K | 16 |
  | 400 k | 0.3 | 0.7 | 0.3 | 0.7 | 0.3 | 0.2 | 0.1 | 0.8 | 3.8 | 0.97 | 64.5 | 205 / 425 | 1.8 K of 456 K | 23 |
  | 1.2 M | 0.5 | 1.1 | 0.6 | 0.8 | 0.6 | 0.6 | 0.2 | 2.0 | 6.9 | 1.43 | 66.8 | 346 / 980 | 4.0 K of 938 K | 22 |
  | 4 M | 2.2 | 5.8 | 3.4 | 1.5 | 1.7 | 1.7 | 0.4 | 4.9 | 23.3 | 2.21 | 69.3 | 2.1 K / 4.6 K | 7.6 K of 1.8 M | 21 |
  | 12 M | 5.3 | 15.1 | 8.9 | 1.3 | 3.8 | 5.0 | 0.7 | 9.8 | 52.6 | 3.22 | 72.0 | 7.4 K / 17 K | 13.9 K of 3.1 M | 22 |
  | 40 M | 31.6 | 74.4 | 27.6 | 1.4 | 6.2 | 6.7 | 1.7 | 16.1 | 171.0 | 5.16 | 74.7 | 40 K / 102 K | 24.0 K of 5.1 M | 34 |
  | 100 M | 51.3 | 109.3 | 67.6 | 1.7 | 8.9 | 9.8 | 4.7 | 24.3 | 286.9 | 6.30 | 76.5 | 398 K / 1.08 M | 39.9 K of 6.9 M | 83 |

  - *Against the gate* (10 Gbp, ≤ 1 h, ≤ 64 GB): 100 M pairs is 28 Gbp, so preloaded time and anonymous memory pass with room. Not-preloaded time is still unmeasured.
  - **EM stops at `max_iter` in every cell**, from 10 k pairs up (`em_iterations` 1000 everywhere), for 0.3–0.6% of detected units. Step 32's question is answered: slow convergence is real and not only a depth effect. So E1 is needed, and step 37 builds it.
  - *Largest component:* 40 K → 398 K units (102 K → 1.08 M pairs) from 40 M to 100 M pairs. That is exponent ~2.5 in depth, as from 12 M to 40 M, so the component is still growing super-linearly (percolation-like). At 1.08 M pairs it is under the EM's 2 M-pair batches and about 1.7 GB in the posterior (~1.6 KB per pair), so E2 and E3 are not yet forced at one run's depth. The pooled 200 M-pair sample may force them.
  - *`presence` is explained:* it converges in every cell (13–83 iterations, `presence_converged` 1), so step 31's 13.7 s at 40 M was not `max_iter`. It now takes 1.7 s at 40 M and 4.7 s at 100 M. Step 31's jump was probably ERR7746321's own structure. Item 3 of step 31 is closed.
  - *Hash (prefetch, step 34):* 40 M 71 → 31.6 s, CPU/wall 1.5 → 2.85. Prefetching the read batches works.
  - *Lookup regressed:* 40 M 39.6 → 74.4 s, CPU/wall 4.4 → 2.7, with the same work as step 17 (24 sampled k-mers per pair). The likely cause is step 34's `MADV_RANDOM` on tier 2: on preloaded pages it turns off fault-around, so each 4 KB page is its own minor fault, and the threads contend on them. Step 37 removes the advice from the packed tables (unit columns keep it). The next run is a direct A/B against this one.
  - *hit\_units still +27 GB RSS* (43.9 → 71.0 GB at 40 M) with `MADV_RANDOM` on the unit columns. This is expected preloaded, because the advice cannot shrink pages already in the page cache (step 34). It is file-backed and reclaimable; the effect is to measure not preloaded.

* **Phase 6, step 37 — step 33's E1, E2 and E3 in code, A1's measurement, step 32 ready to re-run (pooled sample, posterior cells).** Developed on synthetic components: families of 1–7 units sharing 5–60 core k-mers, Poisson hits, ~40% of units absent. These reproduce step 36: 179 of 7.5 K units at `max_iter`, and plain EM needs 14 K steps to converge.
  - **E1, SQUAREM (`_squarem`, inside `_fit_components`, so `em` in all its forms and `em_pin`).** One cycle is two EM steps, the SqS3 step length per component (−‖r‖/‖v‖, capped by a per-component `step_max` that grows 4× when reached), and a stabilising EM step from the extrapolated point. Coverages are kept ≥ 10⁻³ × the second step's, because EM is multiplicative and a coverage of 0 could never grow back. `other` is kept within [0, `other_max`].
    - *Safeguard:* residual-based, not likelihood-based, because the zero-inflated fits have no simple objective. A component whose EM residual grows over the cycle falls back to its two plain steps, and its cap resets to 1.
    - *Iteration count and convergence test:* `max_iter` counts EM steps, three per cycle. Convergence is still tested on a plain EM step, so it means what it did. `step_max` is kept per component through the block compactions, so results do not depend on batching (`test_component_batches_change_nothing` caught this).
    - *Result (synthetic, 7.5 K units):* plain EM converges in 3.9 K steps (was 33.5 K) and agrees with plain EM run to convergence within 6×10⁻⁶. Error at `max_iter` = 1000 fell from 0.044 to 6×10⁻⁴ in coverage. Zero-inflated EM converges in 18 K steps, where plain EM had not converged at 300 K; error at 1000 steps fell from 8.0 to 3×10⁻⁴.
    - *What is left of `em_unconverged_units`:* units that are not converged at `tol` 10⁻⁸ but are already within ~10⁻⁴ of the fixed point. Only 11 of 7.5 K differ by more than 10⁻⁷. Tested: same coverages as plain EM to convergence, in fewer than half the steps (`test_squarem_reaches_plain_em_fixed_point`).
  - **E3, block-wise EM (`_fit_blockwise`).** It applies to components above `MAX_FIT_PAIRS` = 10 M pairs (~1.5 GB of fit memory), so it changes no result at the sizes seen so far. The units are cut into runs of ≤ 10 M pairs. Each round gives every block `BLOCK_STEPS` = 30 EM steps, holding its k-mers' other holders fixed as a per-k-mer `offset` in the expected hits (`_Block.offset`). Rounds repeat until no unit moves by more than the usual test (at most `MAX_BLOCK_ROUNDS` = 1000). `--stats` gains `em_block_rounds`.
    - *First version superseded:* the first version fitted each block to convergence per round. It reported convergence at a wrong point (one unit at 0 against 0.107): a block fitted alone pushes a unit towards 0 that it needs back once the others move, and multiplicative EM brings it back very slowly. Few steps per block fixed it.
    - *Result (one chained component of 170–1.1 K units, blocks of 300–3000 pairs):* within 5×10⁻⁷ of the whole fit for plain EM, at 3–20× its time.
    - *Zero-inflated fits:* the zero-inflated `em` and `em_pin` have several fixed points, and block-wise reaches a different one. Both points have joint EM residual ~10⁻⁸; coverage differs by up to 1.2. So E3 is exact only for the shipped `coverage_em`, and the posterior's zero-inflated refits (a 10 M-pair component) would be path-dependent. Tested: plain EM matches the whole fit (`test_blockwise_em_matches_whole`).
  - **E2, posterior streamed inside a batch.** Each draw's counts are rebuilt from the read rows in chunks of `READ_CHUNK` = 2²² rows. These are integer sums, so the output is identical (tested). Each Gibbs sweep runs over whole k-mers in blocks of about `SWEEP_BLOCK_PAIRS` = 2²⁰ entries, keeping only per-unit totals across blocks, plus each entry's evidence on the last sweep.
    - *Effect on results:* a batch under one block draws the same random numbers as before, so every result to date is unchanged. A bigger batch gets the same model with other random numbers. Keying draws per k-mer would make it block-independent, but `rng.binomial` cannot be vectorised that way, and it was not worth it.
    - *Result (one synthetic component, 1.15 M pairs, 10 draws, blocks of 2¹⁷):* numpy peak 206 → 171 B per pair, same time (82 s). Mean upper bound 1.968 vs 1.970; ambiguity groups 5058 vs 5057.
    - *What this says about the ~1.6 KB per pair (step 28):* most of it is not the sweep's arrays. It is presumably the polars frames and per-read rows (`hash_reads`: several reads per k-mer, sorted and joined per batch). Step 32's posterior cells will measure it on real data before more is streamed.
  - **A1's measurement (`link_cuts`, `--stats` only, stage `links`).** On the largest component of detected units it records `links` (unit pairs sharing hit k-mers) and two weak-link counts: `links_one_kmer` (one shared k-mer) and `links_le2_hits` (≤ 2 hits on the shared k-mers). For each τ in step 33's sweep (0.005–0.2) it records `cut_links_{τ}` and `cut_largest_units_{τ}`, the largest component left after cutting links weaker than τ (strength = shared k-mers / the smaller unit's hit k-mers). The run then shows directly whether a small τ breaks the giant component (fingerprint false hits or error k-mers) or only a large one does (real multi-domain links), which decides whether A1 is needed. Tested on a hand-built case.
  - **Workflow (step 32's runs):**
    - `--query_draws_on` defaults to `''`, meaning every index (`''`, `true` or `false` from the command line all mean this). `--query_draws_at` defaults to `40000000,100000000`. Stub-checked with defaults, with `''`, and with one named index at every cell.
    - The README gives both commands: ERR7738575 to 100 M pairs, and ERR7738575 + ERR7746321 pooled at 40, 100 and 200 M pairs with the posterior at 100 and 200 M.
    - *Fixed after the first attempt:* both FETCH_READS jobs failed with HTTP 403 from `ftp.sra.ebi.ac.uk` over HTTPS. The same download had worked in step 36's run, and the URLs served normally from elsewhere, so the cause is likely throttling or blocking of the HPC host. FETCH_READS now downloads one run at a time, tries HTTPS then FTP with back-off, and retries after waits of 2, 4 and 8 min. `--query_reads` takes several pairs (`R1,R2,R1b,R2b`) and pools them, so runs downloaded on a login node can be pooled too. Stub-checked: one and two local pairs, an odd file count (error), two runs.
    - *Fixed after the second attempt:* Slurm rejected POOL, which had no time limit (the site requires one). POOL now has 6 h and 2 GB, and every process defaults to 4 h and 2 GB, so a process added without resources still submits.
  - **To read from the run:**
    1. `em_unconverged_units` and `em_iterations` per cell (E1 on real data; expect a few hundred units, all near their fixed point).
    2. Lookup time at 40 M and 100 M against step 36 (the tier-2 advice).
    3. On the posterior cells, `posterior_peak_anon` and time at 1–3 M-pair components (E2), then the `COPIES_ERROR` re-check of step 34.
    4. The pooled sample's largest component and `fit_largest_batch_pairs` (E3 engages above 10 M pairs; `em_block_rounds`).
    5. The `cut_*` columns: if a τ ≤ 0.02 cuts the largest component to the size of the second largest, run A1's τ sweep on the fmh benchmark; if only τ ≥ 0.1 does, A1 costs accuracy and E1–E3 stay the route.

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

## Benchmark labels

Pfam is the primary label for benchmarking; KO is kept as a secondary label.

- **Why Pfam.** The released index labels MGnify90 clusters with MGnify's own Pfam hits (no KO or eggNOG annotation run), so Pfam is what users will report and what defaults should be chosen on. Pfam is also open, so truth sets and indexes can be shared.
- **Why keep KO.** fmh-funprofiler profiles KEGG orthologs and is the direct predecessor, so a KO benchmark is the only like-for-like comparison with it (the phase-3 gate was passed on KO). KO results are reported beside Pfam results, but decisions and gates from phase 5 on use Pfam.
- **Pfam truth.** Annotate the benchmark genomes' proteins with Pfam-A (`hmmsearch --cut_ga`, the Pfam release of MGnify's `mgy_proteins_pfam`). A Pfam is present if a read overlaps one of its domain hits (the domain's coordinates, not the whole gene), and its depth is Σ over those domains of aligned bases / domain length, as for KOs. Multi-domain proteins count towards each of their Pfams.
- **Our profiles.** Two index types: Pfam units built from the benchmark's reference proteins (the analogue of the KO units, isolating the method from the reference), and MGnify90 cluster indexes where each Pfam's abundance is Σ over clusters carrying it of `coverage_zi` (the shipped configuration).
- **Other tools on Pfam.** DIAMOND: a read counts towards the Pfams whose domains its alignment overlaps on the subject. HUMAnN: UniRef90 families regrouped to Pfam (`humann_regroup_table -g uniref90_pfam`). fmh-funprofiler and kMermaid stay on KO only: fmh-funprofiler ships KO sketches, and kMermaid's disjoint clusters do not fit multi-domain proteins.

## Test and benchmark data

Three tiers: tiny hand-built fixtures in the repo, a small simulated community for CI, and large benchmarks run outside CI.

| Tier | Data | Tests | Where |
| --- | --- | --- | --- |
| Unit fixtures (< 1 MB) | \~20 proteins with their CDS; hand-built reads: in-frame, each of 6 frames, reverse strand, spanning a stop, with synonymous and non-synonymous changes, with N, with an indel. | Translation, frame filter, silent-change invariance, threshold nesting. | `tests/data/`, generated by a script in `scripts/` so they are reproducible. |
| CI community (\~5–10 MB) | 5 complete genomes spanning GC: *S. aureus* (low GC), *E. coli* K-12, *B. subtilis* 168, *P. aeruginosa* PAO1 (high GC), one archaeon. Reads from InSilicoSeq with known origins; truth = per-gene coverage from read origins mapped through annotations. | End-to-end golden output (insta snapshot), GC effect on frame filter, abundance error bounds. | Repo via Git LFS, or Zenodo with a checksum-verified download step. |
| fmh-funprofiler datasets | Their simulated metagenomes and KO ground truth, plus Pfam truth from the same genomes. | Head-to-head at equal density; compat-mode parity (KO); primary accuracy (Pfam). | Zenodo record from their repo. |
| Divergence ladder | Hold out the source genomes' proteins from the index, keeping relatives at \~95/90/80/70% amino-acid identity. | Sensitivity decay vs alphabet, k, spaced seeds. The key experiment for any amino-acid k-mer method. | Built by the eval pipeline. |
| CAMI II (marine, strain madness, plant-associated) | Public simulated metagenomes with genome-level truth. | Functional truth derived by annotating the source genomes with Pfam (primary) and KofamScan KOs (for fmh-funprofiler) and projecting read origins. | CAMI data portal. |
| Mock communities | ZymoBIOMICS standards (e.g. the gut standard) with published reference genomes. | Real sequencing error and library bias with known composition. | ENA/SRA runs. |
| Real cohort | A subset of HMP2/IBDMDB with published HUMAnN outputs. | Concordance with HUMAnN; runtime at scale. | IBDMDB portal. |
| Negative controls | Shuffled reads, human reads, intergenic-only simulated reads. | False-positive rate; human reads also with and without the human mask. | Generated. |
| Host spike-in ladder | fmh benchmark metagenomes plus simulated T2T-CHM13 reads at 0/50/90/99% host. | Host handling: none, upstream removal, mask, decoy (Additional references). | Generated (InSilicoSeq). |

**Reference databases:** develop on a human-gut biome subset of MGnify Proteins 2026\_07; release on the full set (CC0, so pre-built indexes can be shared). Ground-truth genomes for the CI community and CAMI must be annotated against the same unit definitions, e.g. by mapping their predicted proteins to MGnify clusters, and with the same Pfam release as MGnify's labels.

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
crates/cli/           # phase 8: final Rust binary (index, query, model)
python/<name>/        # prototype: index build, query, model, typer CLI
python/<name>/reference/  # slow pure-Python twins of the kernels (test oracles)
tests/python/  tests/data/  tests/golden/   # golden outputs from the prototype
scripts/              # fixture generation, data download
eval/                 # Nextflow benchmarks, notebooks
docs/spec.md          # algorithm spec, written at end of phase 7
```

Keeping `crates/core` free of PyO3 means phases 6 and 8 reuse the kernels unchanged; only `crates/py` knows about Python.

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
  - **parity**: Rust kernels vs `reference/` twins on fixtures and the CI community; from phase 8, the Rust CLI vs `tests/golden/`.
  - Matrix: ubuntu-latest and macos-latest (arm64, covers NEON). Coverage from `pytest-cov` and `cargo-llvm-cov` to Codecov.
- Use `dtolnay/rust-toolchain`, `Swatinem/rust-cache`, and `taiki-e/install-action` for nextest, cargo-deny and llvm-cov.
- `bench.yml`: manual or weekly `pytest-benchmark` and criterion runs; not a merge gate.
- `release.yml` (phase 9): `cargo-dist` binaries on tags, optional wheels via `PyO3/maturin-action`; bioconda recipe after the first tagged release.
- Dependabot for Cargo, uv and Actions; branch protection on `main` requiring `ci.yml`.

## Additional references (study proteins, hosts, contaminants)

Users need three kinds of extra reference beside MGnify90. **Augment:** a study's own proteins (MAG gene predictions, assembly gene calls and their annotations) are the closest references its reads will ever have, and much of a new study's protein content is not yet in MGnify. **Subtract:** host reads (human, mouse, plant), mitochondria and PhiX spike-in should not register as microbial functions. **Compete:** contaminants that only exist as proteins (diet plants, chloroplasts, fungi) should absorb the reads that match them rather than leave those reads to their nearest microbial relatives. The plan handles all three with two query-time mechanisms: querying several indexes jointly, and a mask sidecar. Both work on a released index without rebuilding it. Index overlays and compaction (the earlier "Extending the reference" design) become a conditional last step, built only if joint queries prove too approximate.

### Can several indexes be queried together?

Yes, and almost nothing breaks, because the model links units at query time. Gather, EM and the posterior find components from the k-mers that actually got hits (`_unit_components` over the hit table), not from build-time components. If the hit tables of several indexes are joined with offset unit ids, a study unit and an MGnify unit that share hit k-mers fall into one component and split those hits in the EM. Units of different indexes that share no hits never interact, so the base profile is unchanged.

| Concern | Effect | Handling |
| --- | --- | --- |
| Hash scheme (k, alphabet, hash, seed) | Units are comparable only if their k-mers are hashed the same way. | Required to match, checked from `meta.json`. One read pass serves every index; the query samples at the largest *t\_max*, and each index looks up only hashes ≤ its own *t\_max*. |
| `n_groups` (score, promiscuity cut, floor choice) | Each build counts only its own units. A k-mer held by 40 base units and 40 study units passes both cuts although it is promiscuous across the union. | Accepted as an approximation. The effect is somewhat larger components and slightly wrong scores; the evaluation measures it. An overlay build (step 5) fixes it by probing the base. |
| Duplicate units (a study protein already in an MGnify90 cluster) | Two units with near-identical k-mer sets split the coverage, so each shows about half. Pfam totals are unaffected when the labels agree. | The existing ambiguity groups (posterior) report such pairs. Assigning study proteins to existing clusters (step 5) removes the duplicates. |
| Fingerprint false hits | Each index adds its own background, Q·ε·*m\_g*/n. | Background computed per index, as now. |
| Unit ids and labels | Unit ids clash across indexes. | Offset ids at load; the output gains an `index` column; each index carries its own label tables. |
| Normalisation (single-copy markers) | Markers come from the base. | Computed from the base only. |
| Cost | Lookups and residency add up across indexes; hashing does not. | A study index (~10⁶ proteins) is MB to GB, negligible next to the base. |

Indexes built with different k or alphabets cannot be joined: their units share no k-mer space and so cannot compete. They can still be queried in one read pass, but their profiles stay separate and may double-count reads.

### Can an index be updated in place?

It can, but it is I/O-heavy. Adding units costs about as much I/O as compaction, and much less than a rebuild; changing existing units is the expensive case.

- *Adding units* (new clusters): their postings spread uniformly over hash space. Tier 2 is packed by hash range, so almost every range changes, which amounts to one sequential rewrite of tier 2 and the dense tier (~50 GB: minutes on HPC). Set ids deduplicated by set hash must be re-resolved for the k-mers that change.
- *Changing existing units* (members added to an MGnify90 cluster): *p\_in* and the floor's choice of kept k-mers need every member's k-mers, and so the cluster's member sequences from the 1 TB release. That cost is high, which is why overlays keep the old *p\_in*.
- *Global ripple:* a new k-mer raises `n_groups` for every unit that holds it, which changes scores, the promiscuity cut, and other units' floors. A strictly exact update recomputes all affected units, which in practice means a partitioned rebuild of those components.
- *Removing* (masks): dropping postings needs the same tier-2 rewrite. A sidecar of masked hashes avoids it.

So the index is treated like an LSM tree: an immutable base, small layers next to it (study indexes, decoy indexes, masks), and an occasional `compact` that merges them into a new base in one sequential pass. Layers come first because they cost nothing to the base and are easy to drop.

### Options compared

| Option | Kind | Handles | Cost | Weakness |
| --- | --- | --- | --- | --- |
| Upstream read removal (hostile or KneadData: bowtie2 or minimap2 against T2T-CHM13) | Subtract, per read | Host, PhiX | One alignment pass, but host-rich samples then hash 10–100× fewer reads | Needs the host genome. A sparse sketch cannot classify single reads, so read-level removal must happen upstream; this tool cannot replace it. |
| Mask sidecar (host genome k-mers ∩ index) | Subtract, per k-mer | Host coding and non-coding reads, mitochondria, PhiX; also flags host-derived MGnify clusters | Build: one six-frame pass over the genome and one merge with tier 2. Query: one `isin` on sampled hashes | Drops microbial k-mers that coincide with host ones: ~10⁻⁵ of postings by chance at 20 letters and k = 11, more where proteins are truly conserved. At Murphy-10 and k = 11, ~2% by chance: use a proteome-only mask there. |
| Decoy index (host or contaminant proteins) queried jointly | Compete | Coding reads of anything with a proteome (diet, chloroplast, fungi, a host without a genome) | A small index plus joint query | Catches only coding reads; host reads are ~98% non-coding. |
| Study index queried jointly | Augment | MAGs, assembly gene calls | A small index plus joint query; no new code beyond step 1 | Duplicate units and `n_groups` approximations (table above) |
| Overlay plus `compact` (old phase-10 design) | Augment | Same, exactly | Cluster assignment by alignment; compaction rewrites tier 2 | Most code; build only if joint queries fall short |
| Rebuild a subset index including the study proteins | Augment | Biome subsets, development | A subset build (laptop to one node) | Not for the full release |
| Post-hoc decontamination with negative controls (e.g. decontam on unit tables) | Subtract, statistical | Kit and reagent bacteria, which have no separate reference because they are in MGnify | None for the tool | Needs control samples; out of scope, but the output tables should stay usable for it |

Recommended defaults:

- **Host-rich samples** (biopsy, skin, oral, sputum): upstream removal in the Nextflow module, plus the human mask as a safety net.
- **Everything else:** the human mask alone. It costs ~10⁻⁵ of postings and removes the ~0.3% per-k-mer chance hits of residual human reads.
- **Decoys and study indexes:** opt-in.

**Why the mask beats a host decoy for genomes.**

- *Coverage of host reads.* The mask is built from the six-frame translation of the whole host genome, so it covers non-coding reads too: repeats, and the low-complexity peptides from translated LINE/Alu.
- *Exactness.* It removes host-derived hits outright, instead of relying on the EM to give them to a decoy. Host variation leaves ~3% of host 11-mers off-reference (heterozygosity ~10⁻³/bp over 33 bp). Those k-mers hit the index at the ~0.3% chance rate, so ~10⁻⁴ of host k-mers still hit by chance.
- *Flagging contaminated clusters.* The per-unit masked fraction flags MGnify90 clusters that come from host contigs gene-called as prokaryotes, a contamination of MGnify itself that no decoy would expose. Such units are reported as `host_like` rather than dropped.

### Implementation plan

Steps 1–2 are small Python changes needed for the phase-7 host ablation, because the result decides whether the release ships a human mask by default. The phase-8 spec and the Rust query must then include them. Steps 3–5 stay in phase 10.

1. **Joint query (`query --extra-index DIR`, repeatable).**
   - `Index.load` each index; error unless k, alphabet, hash scheme and seed match.
   - Sample at the largest *t\_max*; look each hash up in every index whose *t\_max* admits it.
   - Offset unit ids by the cumulative unit counts and join the unit tables (`m`, `pin_sum`, `m_dense`, `pin_sum_dense`, `max_hash_g`, labels) with an `index` column.
   - Leave gather, EM, presence and posterior untouched: they already work on the joined hit table. The dense second pass runs per index.
   - Tests:
     - querying [base] and [base, unrelated] gives identical rows for base units;
     - [A, B] equals one index built from A ∪ B when A and B share no k-mers;
     - a scheme mismatch errors.
2. **Mask sidecar (`mask GENOME.fa INDEX OUT`, `query --mask OUT`).**
   - Cut the genome into overlapping windows (10 kb, overlap 3k − 1 nt) and run them through the existing read kernel with `frames="all"`. Keep hashes ≤ the index's *t\_max*, then sort and deduplicate (human: ~7×10⁸ sampled, ~3 GB of u64 while building, not kept).
   - Look them up in tier 2 and the dense tier, and write `mask.npy` (masked hashes found in the index; expected MB).
   - Write the per-unit decrements of `m`, `pin_sum`, `m_dense` and `pin_sum_dense` (from the masked postings' *p\_in*), and the masked fraction.
   - At query, drop sampled hashes in the mask before lookup, subtract the decrements, and report units with masked fraction > 0.5 as `host_like`. Record the mask's source genome and checksum in the output `meta`.
   - Release one human mask: T2T-CHM13 + rCRS chrM + PhiX174.
   - Tests:
     - an index built from proteins encoded in a fixture "genome" is fully masked by that genome;
     - containment of unmasked units is unchanged;
     - the decrements equal a rebuild that drops the masked hashes (when the floor is not involved).
   - Fallback if the evaluation shows the mask costs too much: a build-time drop rule beside the promiscuity cut in `_select_postings` and the partitioned build, so floored units refill from other candidates.
3. **Decoy role.** Add `role: decoy` in an index's `meta.json` (`index --role decoy`). Decoy units join gather and the EM but are left out of the unit and Pfam tables; one summary row per decoy index reports its hits and coverage. Build decoy indexes with the existing `index` command from a proteome clustered at 90% (or one unit per protein).
4. **Study-index recipe** (Nextflow module, no tool code): contigs or MAGs → pyrodigal (meta mode) → MMseqs2 linclust at 90% → members Parquet, plus Pfam by `hmmsearch --cut_ga` with MGnify's Pfam release → `index` with the base's k, alphabet and seed → `query --extra-index`. Protein FASTA input skips the gene calling.
5. **Overlay and `compact` (conditional).** Build only if step 4's evaluation shows duplicate units or `n_groups` approximations costing more than 2 points of completeness or purity against a rebuild. The design is unchanged from before:
   - assign study proteins to existing MGnify90 clusters (dense-tier containment for candidates, then alignment to their representatives) and cluster the rest at 90%;
   - build the new units under the base's rules, with `n_groups` probed from the base;
   - keep *p\_in* of existing units, and use appendable unit ids;
   - `compact` partition-rebuilds the affected components.

**Constraints on earlier phases** (keep these possible; do not build for them):

- Keep the per-unit build rules local to the unit (*t\_g*, floor, *p\_in*), with only `n_groups` and components global.
- Record the base release, k, alphabet, hash scheme and seed in `meta.json`.
- Keep unit ids appendable and the unit table separate from the hash tables.
- Let the phase-6 Rust lookup probe more than one table and merge their hits.
- Keep per-unit expected counts (`m`, `pin_sum` and the dense versions) as plain columns that a layer can adjust, not baked into the hash tables.

### Evaluation

- **Host spike-in ladder** (phase 7):
  - *Samples.* The fmh benchmark metagenomes mixed with simulated human reads (T2T-CHM13 plus chrM, InSilicoSeq with the same error model) at 0, 50, 90 and 99% host.
  - *Arms.* No handling; upstream hostile; mask; human-proteome decoy; mask + decoy.
  - *Metrics.* Pfam purity and completeness, false detections, `host_like` units detected, and query time.
  - *Negative control.* The existing human-reads control should detect ~0 units with the mask.
  - *Mask cost.* The fraction of postings masked, and the units losing > 10% of their kept k-mers (expected: conserved proteins such as EF-Tu, DnaK, GroEL, ATP synthase).
- **Study ladder** (phase 10): on the divergence ladder, hold genomes out of the index and compare four arms: base alone; base + study index (joint query); overlay, if built; and a full rebuild that includes them. Also confirm that a joint query with an unrelated index leaves every base unit unchanged.
- **Gate:**
  - the mask removes ≥ 99% of false detections from 90% host reads at ≤ 0.1% loss of base-unit completeness;
  - the joint query with a study index recovers ≥ 90% of the completeness gain of a rebuild;
  - extending with one study (~10⁶ proteins) takes minutes on one node.

## Genome mode (predicting genomes from functions)

Given a sample's function profile, predict which reference genomes are present and at what depth. Annotate each genome with the index by running its proteins through the query's first pass, then explain the sample's per-unit hits as a mixture of genome contents with the same gather and zero-inflated EM the unit level uses. Desirable, not essential: phase 11.

**Why from functions.** One index and one read pass give both profiles; genomes are seen in protein space, so synonymous variation does not cost sensitivity; and the EM's allocation gives a stratified table (which genome carries which function's hits), as HUMAnN reports.

**Prior art.**

| Tool | Approach | Relevance |
| --- | --- | --- |
| sylph | DNA FracMinHash containment per genome, zero-truncated Poisson coverage, winner-take-all on shared k-mers. | The accuracy baseline at species level, run on the same genomes. |
| sourmash gather | Greedy minimum set cover of the sample sketch by genome sketches; also works on protein sketches of proteomes. | The detection rule we already use for units (`gather`), one level up. |
| PanPhlAn | Matches a species' gene-family presence/absence profile in a metagenome to reference strains. | Genome identity from gene content, as here; per species only. |
| MetaPhlAn, mOTUs | Clade-specific or single-copy marker genes. | Markers are a restriction to ablate (core units only), and give the genome-equivalents normaliser. |

### Annotating genomes: what differs from reads

The query's first pass already gives what the model needs. With read id = genome, `unit_hits` returns per (genome, unit) the raw hits on the unit's kept k-mers (each k-mer counted for every holder, before any splitting) and the distinct k-mers hit. That table is the genome's content *c\_{G,u}*.

- **Same sampler, same space.** Content must be measured with the tier and thresholds the sample's hits use (tier 2 and *t\_g*; the dense tier only if the fit uses dense hits). Units the sampler cannot see, e.g. a singleton at *t\_base* with ~0.3 kept k-mers, are invisible in both the genome and the sample, so they drop out consistently. This is the sketch-consistency rule of the critique applied to genomes.
- **No splitting at annotation.** Raw hits add across genomes, so a sample of genomes at depths λ\_G has expected raw unit hits *h\_u* = Σ\_G λ\_G *c\_{G,u}* exactly, shared k-mers included. Gather, EM and the posterior are not run per genome: at depth 1 with no error their statistics mean nothing.
- **No containment threshold.** A genome protein that shares one k-mer with an unrelated unit really does put hits there in the sample, so low-containment links stay in *c* with their small weight. A threshold is applied only to the human-readable annotation (best unit per protein).
- **Input.** Long records break the read kernel's assumptions: `stopfree` drops every frame of a contig, since each contains stops. Default input is protein FASTA (genome sets ship `.faa`; otherwise pyrodigal, as in the study-index recipe), hashed by a protein mode of `FastxHits` (the `hash_proteins` scanner behind the same streaming reader). Six-frame windows (the mask sidecar's path) are a fallback: ~8×10⁶ off-frame k-mers per 4 Mbp genome, sampled at 0.2 and hitting at the ~0.3% chance rate, give ~5×10³ chance hits per genome, which the model then has to explain.
- **Batching is built in.** All genomes stream as one input (a list of files, each record tagged with its genome), so there is one process and one lookup pass for the whole set.

**Size.** A 4 Mbp genome has ~1.1×10⁶ amino-acid k-mers, ~2.2×10⁵ sampled at *t\_max* = 0.2. 10⁵ genomes (order of GTDB's species representatives) are ~2.2×10¹⁰ sampled k-mers, about 30 deep metagenomes' worth: ~30 min of hashing and lookup on 8 CPUs at step 31's rates. The content table is ~5×10³ units per genome, ~5×10⁸ rows for 10⁵ genomes (~4 GB as unit-major CSR of (genome u32, hits f32)).

### Genome model: options

| Option | Model | Strength | Weakness |
| --- | --- | --- | --- |
| G0. Genome gather | Greedy cover of hit units by genome contents (`_core.gather` with genome as unit and unit as item). | Existing code, order of minutes. Detection only. | Order-dependent; favours large genomes; no abundances. |
| **G1. Two-stage weighted ZI EM (default)** | *h\_u* \~ Poisson(Σ\_G λ\_G *c\_{G,u}*) over the units any candidate carries, zero-inflated: only a fraction π\_G of G's content is present (strain divergence, incomplete MAGs), as in the unit-level ZI EM. Gather (G0) first for detection. | Reuses `profile.tsv` (`hits`); small (detected units × genomes holding them); the same fixed point as the unit EM. | Raw unit hits are correlated (a shared k-mer counts in each holder), so the Poisson likelihood is composite: point estimates are fine, intervals need care. |
| G2. Joint EM on k-mer hits | Holders of a k-mer are the genomes carrying any unit that holds it; the existing `em` with genome as unit. | Exact likelihood; no correlation issue. | Holder lists blow up for core k-mers (every strain of a species); needs per-(genome, hash) content, ~8× the G1 table. Ablation on small genome sets only. |
| G3. NNLS / sparse regression | Unit coverages on the content matrix, L1 or Dirichlet (α < 1) penalty to pick few genomes among near-identical ones. | Standard; the sparsity prior is a cheap answer to near-identical strains. | Gaussian on counts; same information as G1 with a worse noise model. Sparsity prior can be tried inside G1 instead. |
| G4. Marker units only | G1 restricted to units single-copy and core for each species cluster. | Robust to mobile elements and accessory genes; gives genome-equivalents for Normalisation. | Throws away accessory content, the part that separates strains. |

Recommended: G0 detection, G1 abundances, G4 as an ablation and the normaliser, G2 as an exactness check on the fmh benchmark. Code changes:

- A per-pair `weight` column (default 1, identical results) in `em`, its zero-inflated form and `_fit_components`: the expected count of item *u* from holder G is λ\_G *w\_{G,u}*, and the ZI present fraction counts an item as hit with probability 1 − e^(−λ *w*).
- Candidate screen before the fit: genomes with content-weighted containment (Σ *c* over hit units / Σ *c*) ≥ a floor, so the EM never sees the 10⁵ reference genomes at once. A unit-major CSR (unit → genomes) is read for hit units only.
- Ambiguity groups and intervals from the existing posterior, with the overdispersion of unit hits handled by drawing unit hits from the unit-level posterior rather than Poisson.

**What it reports.** Per genome: λ\_G (depth), π\_G (present fraction), relative abundance, genome-equivalents (Σ λ\_G, the normaliser of Risks). Per sample: the fraction of unit hits explained by genomes (the rest are organisms not in the set), and the function × taxon table below.

### Function × taxon table

A primary output: per sample, abundance of each Pfam (and unit) attributed to each taxon, as HUMAnN's stratified table.

- **Attribution.** After the G1 fit, unit *u*'s split hits (the unit-level EM's allocation, so a k-mer shared across units is not counted twice) go to genome G with responsibility *r\_{G,u}* = λ\_G *c\_{G,u}* / Σ\_H λ\_H *c\_{H,u}*. Hits above what the fit predicts (*h\_u* − Σ\_G λ\_G *c\_{G,u}*, floored at 0), and all hits on units no candidate genome carries, go to `unclassified`. This is the exact E-step of G1, so it costs nothing beyond the fit.
- **Ranks.** Genomes roll up to species, genus and family by the `taxonomy` column; rows are `Pfam|rank:taxon` plus `Pfam|unclassified`, with the unstratified total kept so rows sum to it. Units roll up to Pfam as the unstratified profile does (multi-label units count for each label).
- **Ambiguity.** Genomes in one posterior ambiguity group (near-identical content) are not split; their share is reported at the group's lowest common rank. Intervals per row come from the genome-level posterior draws.
- **Normalisation.** Same units as the profile (hits, coverage, copies per genome-equivalent); a per-taxon option divides by that taxon's genome-equivalents (copies of the function per cell of that taxon).
- **Limits.** Taxa come only from the reference genome set: a biome the set covers poorly leaves most hits `unclassified`. Within a species, protein k-mers cannot see synonymous differences between strains, so strain attribution is weak; species and genus are the ranks to read. A unit-level fallback (lowest common taxon of a MGnify90 cluster's member contigs, if the release gives contig taxonomy; to check) could label part of `unclassified` and is built only if the unclassified fraction on real samples is large.


**Reference sets.** GTDB species representatives (CC BY-SA 4.0: check share-alike before shipping a precomputed table), MGnify genome catalogues per biome, and a study's own MAGs through the study-index recipe's pyrodigal step. Units carried by no candidate genome are left out of the fit and counted as unexplained.

### Implementation plan

1. **`annotate-genomes INDEX GENOMES.tsv OUT`** (columns `genome`, `path`, optional `taxonomy`): protein mode for `FastxHits` (Rust, small), read id = genome; writes `genome_units` (unit-major CSR: unit → (genome, hits, kmers)) and `genomes` (genome, proteins, Σ *c*, taxonomy), with the index's `meta.json` checksum so a mismatched index errors. Tests: a fixture "genome" made of an index's member proteins carries those units; two genomes annotated together equal each annotated alone; reads simulated from a genome at depth λ give unit hits ≈ λ·*c*.
2. **`genomes PROFILE.tsv GENOME_INDEX OUT`**: screen, G0, G1 per component, outputs above. Tests: the weighted EM with all weights 1 equals today's `em`; two genomes with disjoint content recover their depths; a genome held out with a 95%-identical relative in the set is reported as the relative with π < 1.
3. **Function × taxon table and posterior** (responsibilities, rank roll-up, ambiguity groups at their common rank, `unclassified`). Tests: two genomes with disjoint content get all of their own units' hits; a unit shared by two genomes at depths 1 and 3 splits 1:3; hits on a unit no genome carries are all `unclassified`; rows sum to the unstratified total. Then the G2 and G4 ablations.
4. **Nextflow:** genome-set annotation as a one-off HPC job published beside the index; the profile pipeline gains an optional genome step.

**Constraints on earlier phases** (keep these possible; do not build for them):

- Keep `unit_hits` callable with any per-record id, not just read ids.
- Keep the raw per-unit `hits` (tier 2, before splitting) in `profile.tsv`.
- Keep the per-unit hits after the EM split in `profile.tsv` too (the function × taxon table splits those).
- Keep gather and the EM generic over (holder, item) pairs, as they are now.

### Evaluation

- **fmh benchmark** (phase 11): its 64 genomes have exact simulated abundances. Annotate them plus distractors (other KEGG genomes from the same Zenodo record), against the full MGnify index and the KO index. Metrics: genome purity, completeness, F1, abundance L1 and Spearman. Baselines on the same genomes: sylph (DNA) and sourmash gather on protein sketches.
- **Function × taxon table** (phase 11): truth per (genome, KO or Pfam) from the 64 genomes' annotations and simulated abundances, rolled up to species and genus. Metrics: stratified abundance L1, (taxon, function) F1, and the fraction of hits attributed to the right taxon; baseline HUMAnN 3.9's stratified output (KO), already run in phase 5.
- **Hold-out ladder:** remove each source genome from the set, keeping relatives at ~95/90% ANI, as the divergence ladder does for proteins: the relative should be reported, with π\_G falling as identity falls.
- **Cost:** annotation time per genome and per 10⁵ genomes; fit time against the query's.
- **Gate:** see phase 11 in the table.

## Genome-informed unit presence (Bayesian update from taxa)

Use the genomes genome mode finds to set a prior on which units are present, then update it with each unit's own hits. A unit missed in a low-coverage genome that always carries it is then reported as probably present; a unit missed in a high-coverage genome becomes a confident absence. Desirable, not essential; it needs genome mode (phase 11) or an external taxonomic profile.

**A companion tool, not a profiler stage** (working name `kfp-prior`), as Bracken is to Kraken. It reads the profiler's output files and writes its own; the profiler never imports it.

- *The profiler's output stays evidence-only.* Imputed presence is a model of taxonomy, not of the reads, and should not be mixed into `profile.tsv` by default.
- *Taxa from any source.* Genome mode's `genomes.tsv` or a sylph profile (option B2) are interchangeable inputs, so the companion does not depend on phase 11 shipping.
- *Different data and cadence.* Its reference is a carriage table built from a genome set with taxonomy (GTDB, MGnify catalogues), not from the index; it can be rebuilt when a catalogue updates, without touching the index.
- *Interface (the file contract):*
  - **From the profiler:** `profile.tsv` with per-unit own hits *h*, `present_prob` and its likelihood-ratio term, expected hits per unit of coverage (from `m`, `pin_sum`), component id and Pfam labels, and the index checksum.
  - **From genome mode:** `genomes.tsv` (genome, *P\_G*, λ\_G, π\_G) and the best-unit annotation per genome protein. **From sylph:** genome, abundance, ANI, mapped to annotated genome ids.
  - **Written by the companion:** `presence.tsv` (unit, `prior`, `present_prob_updated`, `expected_hits`, `imputed`) and `pfam_presence.tsv` (observed, and observed + imputed).

**Is it well formed?** Yes, with four qualifications.

- *It is imputation, so it adds no new information.* The prior is a function of the taxonomic profile. Imputed units raise completeness and give pathway-level answers ("this sample can do X"). They cannot show functional differences beyond what taxonomy already shows (the PICRUSt critique). So `observed` and `imputed` are reported separately, and differential analyses should use observed values only.
- *Unit is not the same as function.* A MGnify90 cluster tracks lineage closely, often species or genus. The sample's strain may carry the same function as a sibling cluster at < 90% identity, which is then hit in place of the expected unit. A zero for the expected unit is then right about the unit and wrong about the function. Siblings in the same connected component that are detected explain the missing unit away (step 2). Pfam-level presence is the summary users should read.
- *The taxa come from the same reads.* Genome mode infers genomes from unit hits, and the update then reuses those hits. A single unit counts for little in its genome's fit (~5×10³ units per genome), and for a unit with zero hits the double counting only pulls its genome down. So the cavity correction (the genome fit without unit *u*) is left out, with an ablation to measure it. An external taxonomic profile (sylph on the same reads) can be used instead and does not have this problem.
- *Genome false positives become function false positives.* One falsely detected genome imports thousands of units at its core frequencies. The prior is weighted by the genome's presence probability, and detection must be calibrated before this step can be.

**Prior art.** PICRUSt2, Tax4Fun2 and Piphillin predict functions from 16S taxa alone, with no read evidence. PICRUSt2 uses hidden-state prediction on a tree (castor). HUMAnN aligns to the pangenomes of detected species before the translated search: it is taxon-informed, but only through the search space, with no prior. CheckM uses lineage-specific marker sets, the same "always present in this clade" signal, applied to MAG completeness. MGnify genome catalogues (e.g. UHGG) ship per-species pangenomes with core/accessory calls (Panaroo presence/absence).

### Model

For unit *u* and detected genome *G* (presence probability *P\_G*, depth λ\_G, present fraction π\_G from G1):

- *Carriage frequency* *q\_{G,u}*: the probability that a genome related to *G* carries *u*. Phylogenetic similarity is done by taxonomy-rank shrinkage, a beta-binomial down the GTDB ranks: *q\_{s,u}* = (*n\_{s,u}* + α *q\_{genus,u}*) / (*N\_s* + α), and likewise genus toward family. *n\_{s,u}* is the number of species *s*'s genomes that carry *u* (best-unit annotation from `annotate-genomes`), *N\_s* the species' genome count, and α is fitted by leave-one-genome-out log loss. A tree-based alternative (castor-style hidden-state prediction, or a phylogenetic kernel) is used only if rank shrinkage loses on held-out genomes.
- *Prior* (noisy-OR over genomes): carriage via *G*, *z\_{G,u}* \~ Bernoulli(*ρ\_{G,u}*), *ρ\_{G,u}* = *P\_G* · *q\_{G,u}*, independent across genomes. A unit that no detected genome carries keeps today's global prior (the *w\_h* / A of `present_prob`).
- *Likelihood.* Own hits *h\_u* \~ Poisson(background + Σ\_G *z\_{G,u}* λ\_G π\_G *e\_u*), where *e\_u* is the expected hits per genome copy on *u*'s kept k-mers (*c\_{G,u}* from the content table). π\_G discounts for strain divergence; it mixes gene loss with sequence divergence, which the ablation checks.
- *Zero hits (the case this exists for) has a closed form.* Background is ≈ 0 at *h* = 0, so genomes factorise:
  P(*z\_{G,u}* = 1 | *h* = 0) = *ρ* e^(−λπ*e*) / (1 − *ρ* + *ρ* e^(−λπ*e*)), and P(*u* present | 0) = 1 − Π\_G (1 − that).
  Worked example: *q* = 0.95, *P\_G* = 1. At λπ*e* = 0.5 expected hits, a zero gives P = 0.92 (missing is unsurprising). At 5 expected hits it gives 0.11 (a confident absence).
- *Hit units (h > 0).* Multiply the odds of the existing `present_prob` by the per-unit prior odds over the global prior odds. Nothing else changes, so a unit with no genome-informed prior gets today's value exactly.
- *Abundance (optional).* An imputed unit's coverage is Σ\_G P(*z\_{G,u}* = 1 | 0) λ\_G × its expected copies. Report it, but leave it out of the EM.

### Options

| Option | What | Strength | Weakness |
| --- | --- | --- | --- |
| **B1. Post-hoc closed form (default)** | The model above, run on the profile and genome outputs; per-unit prior, closed-form update | A few hundred lines; reuses `present_prob` and genome-mode outputs; components independent | Double counting (cavity left out); independence across genomes |
| B2. External taxa | B1 with genomes and depths from sylph (or MetaPhlAn) instead of G1 | No circularity; sylph is the better species detector | Needs sylph genome ids mapped to annotated genomes; π\_G must come from sylph's ANI |
| B3. Joint hierarchical model | Gibbs over genome presence, *z\_{G,u}* and unit coverage together | Exact: no double counting, correct joint uncertainty | Much more code and time; only worth it if B1 is miscalibrated |
| B4. Function-level prior | Prior and update at Pfam (or component) level, not unit | Avoids the sibling-cluster problem | Loses unit resolution; Pfam presence from units is an aggregation anyway |

Recommended: B1, reported at both unit and Pfam level, with B2 as an ablation arm. B3 only if B1's calibration fails.

### Implementation plan

1. **Carriage table (`kfp-prior build GENOME_UNITS TAXONOMY OUT`).**
   - Input: the genome-mode content (best unit per protein above the annotation threshold) and GTDB taxonomy per genome.
   - Output: `carriage` (clade, unit, *n*, *N*) at species, genus and family, and fitted α per rank.
   - Fit α on held-out genomes, then write *q* per (species, unit) for units with *q* above a floor (e.g. 0.05). The table is sparse; the floor is a parameter.
   - Tests: a species whose genomes all carry *u* gives *q* → 1 as *N* grows; a one-genome species shrinks toward its genus; leave-one-out log loss beats both no shrinkage and a genus-only prior on a fixture.
2. **Update (`kfp-prior update PROFILE.tsv GENOMES.tsv CARRIAGE OUT`).**
   - Per detected genome, read its species' *q* row; noisy-OR into per-unit priors; closed form for undetected units and odds rescaling for detected ones.
   - Explain-away: when a unit has *h* = 0 and a sibling in its component is detected, the unit keeps its unit-level value but the function is marked `sibling_hit`. The sample's strain probably carries the function as the sibling's variant.
   - Output: `presence.tsv` and `pfam_presence.tsv` (above); `imputed` is true when *h* = 0 and `present_prob_updated` ≥ 0.5. `profile.tsv` is not modified. Errors if the profile's and the carriage table's index checksums differ.
   - Tests: with no genomes detected, every value equals the profile's `present_prob`; the worked example's numbers; a unit at *q* = 1 in a genome with λπ*e* → ∞ and *h* = 0 gives P → 0.
3. **Taxa adapters and ablations:** a sylph adapter (B2: genome ids mapped to annotated genomes, π\_G from ANI); the cavity correction (refit G1 per component without unit *u*, small sets only); π\_G on vs off; species-only *q* vs shrinkage.
4. **Packaging and Nextflow:** a separate Python package and CLI in this repository (no import of `kmer_functional_profiler`; tests use fixture files in the contract's format), split into its own repository only if its release cycle diverges. Nextflow: the carriage table is built once per genome set; the profile pipeline gains an optional module after the genome step.

**Constraints on the profiler** (keep these possible; do not build for them):

- `profile.tsv` columns in the file contract are stable and documented, with the index checksum in its header or `meta`.
- Genome mode keeps the per-protein best-unit annotation, not only the raw content.
- Genome input carries taxonomy (the `taxonomy` column already planned).
- `present_prob` keeps its prior as a separable term (global odds × likelihood ratio).

### Evaluation

- **Depth ladder (fmh benchmark):** subsample reads so the 64 genomes span 0.05–5× coverage. Truth is the units carried by each source genome. Metrics per coverage bin: unit and Pfam completeness and purity, observed-only vs updated; and calibration of `present_prob_genome` (reliability curve, with the observed fraction among units with *h* = 0).
- **Hold-out arm:** remove each source genome from the reference set, keeping same-species relatives. The prior then comes from relatives, so accessory units are where precision should fall; report precision by *q* bin.
- **False-genome control:** add a distractor genome forced as detected at *P\_G* = 0.5 and count its imputed units.
- **Expected outcome:** large completeness gains for core units (*q* > 0.9) of genomes at 0.1–1×, where floored units still miss often (8 kept k-mers at 0.2× is ~1.6 expected hits, P(0) ≈ 0.2). Little change for accessory units, where the evidence dominates anyway. Nothing for organisms absent from the genome set, which limits the benefit to well-catalogued biomes (human gut more than soil).
- **Gate:** at 0.1–1× coverage, Pfam completeness rises by ≥ 10 points at ≤ 2 points of purity loss; calibration error ≤ 0.05 on *h* = 0 units; and no change when no genomes are detected.

## Sequence similarity (containment AAI)

sylph reports a containment ANI per genome: in `profile` for the genomes it assigns abundance to, and in `query` for every genome above a minimum ANI, with no reassignment of shared k-mers. The same is possible here in protein space, as a containment AAI per MGnify90 unit. It is cheap because the quantities it needs are already in the profile. It adds granularity of one kind, how far the sample's variant is from the reference and which references bracket it, rather than finer functional units. Desirable, not essential; steps 1–2 are small enough to build whenever phase 6 allows, and calibration belongs to the phase-7 divergence ladder.

**Estimator.** Under independent substitutions at amino-acid identity *a*, an exact k-mer survives with probability *a*^k, so *â* = *c*^(1/k), where *c* is the fraction of the reference's kept k-mers present in the sample. This is sylph's formula with protein k (11, so *c* = 0.57 at 95%, 0.31 at 90%, 0.17 at 85%).

- *What c is measured against.* A unit's kept k-mers are the union of its members' k-mers, so raw containment is diluted even for a strain identical to one member (median π ≈ 0.25 in the phase-4 simulation). `copies_zi` = present k-mers / `pin_sum` already divides by an average member's kept k-mers. So *c* = min(1, `copies_zi`), and *â* is the AAI to a typical member of the cluster.
- *Coverage.* At low coverage, k-mers the strain carries are missed by chance. The zero-inflated EM already separates the present fraction from coverage, as sylph's zero-truncated Poisson does, so `copies_zi` is coverage-corrected. Sequencing errors thin λ and do not reach *c*.
- *Shared k-mers.* After EM, a k-mer's presence is credited to the units that explain it, so the profile-mode AAI is not inflated by a close neighbour's reads.

**Two modes, as in sylph.**

| Mode | Units | *c* from | Meaning |
| --- | --- | --- | --- |
| Profile (`aai`) | Units gather keeps | `copies_zi` after zero-inflated EM; interval from posterior draws of `present` | How far the sample's variant of this function is from its cluster |
| Blanket (`aai_naive`) | Every hit unit (the profile already lists them all) | Per unit, independently: k-mers hit / `pin_sum`, divided by 1 − e^(−λ̂), with λ̂ the zero-truncated Poisson MLE from that unit's hits per hit k-mer (sylph's correction) | Similarity of the reads to each reference, shared k-mers counted for every holder |

Reporting for "all references" is feasible here because a unit with no hits has no measurable AAI. Like sylph's minimum ANI, it is below detection. The profile already has one row per hit unit, so the blanket mode adds a column, not a table. A `--min-aai` filter controls output size. When λ̂ is too small to correct (hits per hit k-mer ≈ 1), `aai_naive` is the uncorrected value and is flagged as a lower bound.

**Measurable range and precision** (k = 11, high coverage, binomial error on *c* only):

| Kept k-mers per unit | Lowest AAI one k-mer can show | SE of *â* at 97% (*c* ≈ 0.72) |
| --- | --- | --- |
| ~8 (floor, tier 2) | 0.83 | ~0.02 |
| ~30 | 0.73 | ~0.01 |
| ~100s (dense tier) | ≤ 0.66 | ≤ 0.006 |

So the useful band is ~85–100% on the floor index and ~75–100% with the dense tier, close to sylph's 90–100% ANI band. Precision at the floor is enough to tell 99% from 93% but not 97% from 95%. That is a further argument for the dense tier on detected units.

**Known biases** (to measure, then calibrate):

1. *Clustered variation.* Substitutions concentrate in variable regions, so more k-mers survive than *a*^k predicts, and *â* runs high. Mash and skani show the same effect. If the ladder shows a bias above ~1 AAI point, fit a monotone correction per (alphabet, k) on the ladder and store it in the index metadata.
2. *Union mosaic.* A sample k-mer missing from its nearest member may sit in another member, so *â* lies between the AAI to the nearest member and the AAI to the cluster's mosaic. Members longer or shorter than average shift `pin_sum`; `len_cv` measures that spread. Values above 1 are capped.
3. *Gene-end loss.* Under the stop-free frame filter, the last ~40 residues lose coverage structurally. That looks like divergence: about −1 AAI point for a 330 aa protein and −3 for a 100 aa one. `--frames edges` removes most of it; otherwise correct `pin_sum` for expected end loss at the read length.
4. *Reduced alphabets* measure identity in the reduced alphabet, which is higher than AAI. Calibrate per alphabet and name the column accordingly.
5. *Strain mixtures.* Several variants of one unit pool their k-mers, so *c* exceeds any single variant's and *â* runs high. This is sylph's limitation too. The planned dispersion flag ("multi-variant") marks these units.

**Granularity beyond MGnify90.** MGnify90 clusters are already finer than function, so extra resolution here is about variants and origin, not function. Options, cheapest first:

| Option | What it gives | Cost | Verdict |
| --- | --- | --- | --- |
| A. Per-unit AAI (`aai`) | Novelty per detected function: "present, as a variant ~93% to its cluster". Aggregated per Pfam (coverage-weighted), it gives how novel that function's carriers in the sample are. | A column; one more EM fit when the shipped estimator is not ZI | Build (step 1) |
| B. Neighbour AAIs | For each detected unit, the other hit units in its component with `aai_naive` ≥ threshold: the references that bracket the sample variant (e.g. 96% to A, 91% to B). This is complementary to ambiguity groups, which report split uncertainty, not similarity. | A `component` column, so users group the profile by it; no new table | Build (step 2) |
| C. Member placement | The nearest MGYP member within a unit and its AAI, with the runner-up. Member ids link to biomes and assemblies (`mgy_biomes`, `mgy_clusters`); the protein release has no taxonomy. | The index keeps only *p\_in* per k-mer, so this needs a member sidecar: per kept dense-tier k-mer, a bitset over at most *M* members per unit (16–64, chosen to span the cluster's diversity), unit-major and read by `pread` for detected units, like the dense tier. Disk ~*M* bits per posting. | Only if step 3's data show it matters (see below) |
| D. Alignment refinement | Exact per-read identity and best member: reads of detected units aligned (DIAMOND) to those units' member sequences | Needs a unit-ordered member sequence store; members are scattered over the 1 TB sequence file | External recipe at most (phase 9's hybrid pipeline); not in the tool |
| E. Finer units (95% or member-level) | Finer EM units | More than 10⁹ units and many more shared k-mers; already rejected (see Database: 95% vs 90%) | No |

Within-unit placement (C) helps only when a cluster is diverse and the sample variant sits close to one member. The decision rule: build C if, on the real fmh-benchmark and MGnify-subset queries, at least ~20% of detected abundance is in units with `aai` ≥ 0.97 whose members are diverse (mean *p\_in* of kept k-mers ≤ 0.5). Otherwise *â* to the cluster says nearly as much.

### Implementation plan

1. **`aai` in profiles.** `aai` = min(1, `copies_zi`)^(1/k), plus `aai_lo`/`aai_hi` from the posterior's draws of `present` (`posterior_zi` already draws them). Null when present k-mers < 3. A `--aai` flag fits ZI when it is not already fitted. Tests: on the simulation's identity ladder (dense), the median *â* is within 0.01 of the true identity at 95% and 90%; a strain identical to a member gives *â* ≥ 0.99.
2. **`aai_naive` and `component` on every hit unit.** Closed-form zero-truncated Poisson λ̂ per unit (vectorised Newton), a lower-bound flag, and `--min-aai` to drop rows. `component` uses the existing `_unit_components` labels. Tests: with no shared k-mers, `aai_naive` equals `aai`; a unit sharing all its k-mers with a well-covered neighbour gets `aai_naive` ≈ 1 but `aai` ≈ null (explained away).
3. **Calibration (phase 7, divergence ladder and simulation).** Truth is the alignment identity of each source protein to the nearest member of its unit and to the representative. Report bias and interval coverage by identity bin, protein length, frame mode, alphabet and tier. Fit the correction in bias 1 only if needed. Also measure the decision rule for C.
4. **Member sidecar (conditional on step 3):** option C above.

**Gate.** With the dense tier at a rate phase 6 accepts: |median bias| ≤ 0.01 AAI over 0.85–1.0 after any calibration; `aai` interval coverage ≥ 0.9; Spearman(*â*, true identity) ≥ 0.9 on the ladder. On the floor index alone, report the same metrics without a gate.

## Risks and open questions

The biggest risk is that the gain over fmh-funprofiler with a lower scaled value is too small to justify a new tool.

- **MGnify scale.** A floor on all 1.66×10^9 clusters would need \~160 GB of index; on the 0.45×10^9 non-singletons, \~45 GB. Phase 6 (see Progress log): the Python prototype cannot build the full release on one node, so full-scale storage and compute are measured and tuned there, before the phase-7 ablations.
- **Component size.** Promiscuous k-mers can chain clusters into one giant component, which makes the EM serial. Measure component sizes on the development subset in phase 2 and tune the N-clusters cut-off.

* **Marginal novelty.** Setting sourmash to scaled = 100 may recover most of the completeness gap at modest cost. Run that baseline in phase 3 before building phase 4.
* **Query density set by the smallest units.** If *n\_min* forces *t\_max* near 1, the query does almost no sparsification. Measure the distribution of *t\_g* on the real database early.
* **k and the within-component split.** Longer k makes more of a unit's k-mers unique, which helps the EM separate sibling units, but fewer k-mers survive divergence. For siblings at 90% identity and a sample strain at 95% to one of them, the fraction of sampled k-mers that hit that unit uniquely is ≈ 0.95^k (1 − 0.9^k): 0.36, 0.39, 0.39, 0.38, 0.37 at k = 7, 9, 11, 13, 15. Separation is flat near k = 11, and rises with k only for strains almost identical to a reference (1 − 0.9^k: 0.52 at k = 7, 0.69 at 11, 0.79 at 15). So k is not a strong lever for the split; phase 7 measures it on split error and ambiguity-group size alongside sensitivity. If it matters, a targeted fix is a longer-k disambiguation tier for large components only, not a longer k everywhere. Near-identical genomes in the function × taxon table share the same proteins, so no protein k separates them.
* **Divergence.** Exact amino-acid k-mers miss distant homologs regardless of sampling. If the ladder shows steep decay below 80% identity, reduced alphabets or spaced seeds become mandatory, not optional.
* **Shared k-mers and hierarchy.** EM at protein-cluster level, then aggregation to function, is likely better than EM directly on functions. Untested.
* **Normalisation.** Which single-copy marker set, and whether to report per-genome copies by default.
* **Frame filter at high GC.** Keeps \~3 frames at 70% GC; acceptable, but check false positives there specifically.
* **Additional references.**
  * *Joint queries.* They inherit per-index `n_groups`, so promiscuity cuts and scores are only approximate across indexes, and duplicate units split coverage. The step-5 overlay is the fix, built only if measured to matter.
  * *Human mask.* It removes truly conserved microbial k-mers, and with reduced alphabets ~2% of postings by chance; with reduced alphabets, mask with the host proteome only.
  * *Host contamination in MGnify itself.* MGnify90 clusters from host contigs are flagged `host_like`, not removed; whether to drop them from the release is open.
* **Genome mode.**
  * *Visible content.* Genomes from biomes MGnify covers poorly have few units the sampler sees, so little power; report Σ *c* per genome so users can tell.
  * *Shared mobile content.* Plasmid and phage units carried by many genomes can prop up genomes whose own content is absent; the zero inflation (π\_G) and G4's core-only fit are the checks.
  * *Near-identical strains.* Strains sharing almost all units are not identifiable; the posterior's ambiguity groups report them rather than a forced split.
  * *Two-stage correlation.* G1's composite likelihood may give intervals that are too narrow; check coverage on the fmh benchmark before reporting them.
* **Genome-informed presence.**
  * *Imputation is not evidence.* Imputed units follow from taxa; downstream statistics on them re-test taxonomy.
  * *Sibling clusters.* A strain carrying the function as a < 90% variant makes the expected unit look absent; read function-level presence, not unit-level.
  * *Double use of the reads.* Genomes come from the same unit hits; the cavity ablation and the sylph arm (B2) measure the effect.
* **Containment AAI.** Clustered variation and the union-of-members k-mer set both push *â* above the true identity, and gene-end loss pulls it below. At the floor (~8 kept k-mers) its precision is ±2 AAI points, so it is only informative with the dense tier. Read it as similarity to a cluster, not to a member.
* **Future work:** 30% families by mapping MGnify90 representatives onto the 128.7 M MGnify30-C2 representatives, as a coarser level for floors, EM partitions and annotation.
* **Open:** whether KO/eggNOG labels are worth the annotation run for users, or Pfam suffices. Benchmarking uses Pfam (see Benchmark labels), with KO only for the fmh-funprofiler comparison.

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
* [PanPhlAn (strain identification from gene-family profiles)](https://github.com/SegataLab/panphlan)
