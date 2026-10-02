"""Generate the unit-test fixtures in tests/data/ (deterministic; rerun to regenerate).

Twenty random proteins are back-translated into CDSs followed by random flanking DNA.
Each protein yields two read pairs (300 bp fragment, 2 x 150 bp); R1 gets one of the
test cases below and R2 is the reverse complement of the fragment's far end, so the
true frames cover all six. ``truth.tsv`` records each read's source protein, true
frame and amino-acid span.

``mini_release/`` holds a few MGnify90 clusters in the release's Parquet schemas, for
testing ``workflows/mgnify-subset`` without the real release. ``mini_fmh/`` mimics the
fmh-funprofiler benchmark inputs (Zenodo 10055954) for ``workflows/fmh-benchmark``: three
genomes with gene mapping tables, their proteins, gene-to-KO table and KO sketches, plus
``Pfam-mini.hmm``: Pfam-style HMMs (with GA cut-offs) built by ``hmmbuild`` from segments of
some of those proteins, one of them carrying two, so Pfam truth has sub-gene domains.
"""

import gzip
import random
import subprocess
import tempfile
from pathlib import Path

import polars as pl
import sourmash
from sourmash.signature import save_signatures_to_json

from kmer_functional_profiler.reference import BASES, CODE_11, reverse_complement

OUT = Path(__file__).resolve().parent.parent / "tests" / "data"
AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"
CASES = ["clean", "synonymous", "nonsynonymous", "n", "indel", "stop"]
READ_LEN, FRAG_LEN, FLANK_LEN = 150, 300, 300
GUT, MARINE = "root:Host-associated:Human:Digestive system", "root:Environmental:Aquatic:Marine"

CODONS = [a + b + c for a in BASES for b in BASES for c in BASES]
SYNONYMS = {aa: [c for c, t in zip(CODONS, CODE_11, strict=True) if t == aa] for aa in CODE_11}


def substitute_codon(read: str, offset: int, rng: random.Random, synonymous: bool) -> str:
    """Replace one in-frame codon of ``read`` with a synonymous or missense codon."""
    starts = [
        i
        for i in range(offset, len(read) - 2, 3)
        if len(SYNONYMS[CODE_11[CODONS.index(read[i : i + 3])]]) > 1
    ]
    i = rng.choice(starts)
    aa = CODE_11[CODONS.index(read[i : i + 3])]
    if synonymous:
        choices = [c for c in SYNONYMS[aa] if c != read[i : i + 3]]
    else:
        choices = [c for c, t in zip(CODONS, CODE_11, strict=True) if t not in (aa, "*")]
    return read[:i] + rng.choice(choices) + read[i + 3 :]


def truth_row(pair: int, mate: int, protein: int, cds_start: int, case: str) -> str:
    """TSV row for a read covering the CDS from ``cds_start`` (reverse strand for mate 1)."""
    offset = -cds_start % 3
    aa_start = (cds_start + offset) // 3
    fields = (pair, mate, f"prot{protein}", 3 * mate + offset, aa_start)
    return "\t".join(map(str, (*fields, aa_start + (READ_LEN - offset) // 3, case)))


def fastq(name: str, seq: str) -> str:
    return f"@{name}\n{seq}\n+\n{'I' * len(seq)}\n"


def main() -> None:
    rng = random.Random(20260928)
    OUT.mkdir(parents=True, exist_ok=True)
    proteins, cdss, r1s, r2s = [], [], [], []
    truth = ["read\tmate\tprotein\tframe\taa_start\taa_end\tcase"]
    for p in range(20):
        protein = "M" + "".join(rng.choice(AMINO_ACIDS) for _ in range(rng.randint(120, 400)))
        cds = "".join(rng.choice(SYNONYMS[aa]) for aa in protein) + rng.choice(SYNONYMS["*"])
        context = cds + "".join(rng.choice("ACGT") for _ in range(FLANK_LEN))
        proteins.append(f">prot{p}\n{protein}\n")
        cdss.append(f">prot{p}\n{cds}\n")
        for j in range(2):
            pair = 2 * p + j
            case = CASES[pair % len(CASES)]
            # A stop-case R1 crosses the stop codon; its R2 lies in the flank.
            start = len(cds) - 90 if case == "stop" else rng.randrange(len(cds) - 3 - FRAG_LEN)
            r1 = context[start : start + READ_LEN]
            r2 = reverse_complement(
                context[start + FRAG_LEN - READ_LEN : start + FRAG_LEN].encode()
            ).decode()
            offset = -start % 3
            if case in ("synonymous", "nonsynonymous"):
                r1 = substitute_codon(r1, offset, rng, synonymous=case == "synonymous")
            elif case == "n":
                i = rng.randrange(READ_LEN)
                r1 = r1[:i] + "N" + r1[i + 1 :]
            elif case == "indel":
                i = rng.randrange(30, READ_LEN - 30)
                r1 = r1[:i] + r1[i + 1 :]
            r2_start = start + FRAG_LEN - READ_LEN
            truth.append(truth_row(pair, 0, p, start, case))
            truth.append(truth_row(pair, 1, p, r2_start, "flank" if case == "stop" else "clean"))
            r1s.append(fastq(f"pair{pair}/1", r1))
            r2s.append(fastq(f"pair{pair}/2", r2))

    (OUT / "proteins.faa").write_text("".join(proteins))
    (OUT / "cds.fna").write_text("".join(cdss))
    (OUT / "truth.tsv").write_text("\n".join(truth) + "\n")
    for name, records in (("reads_1.fastq.gz", r1s), ("reads_2.fastq.gz", r2s)):
        with (
            (OUT / name).open("wb") as f,
            gzip.GzipFile(filename="", mode="wb", fileobj=f, mtime=0) as gz,
        ):
            gz.write("".join(records).encode())


def mini_release(rng: random.Random) -> None:
    """Twelve clusters of 1-5 mutated copies, with member ids scattered over 1..999."""
    out = OUT / "mini_release"
    out.mkdir(exist_ok=True)
    ids = iter(rng.sample(range(1, 1000), 60))
    clusters, members, proteins, pfam = [], [], [], []
    for c in range(12):
        base = "M" + "".join(rng.choices(AMINO_ACIDS, k=rng.randint(40, 150)))
        size = rng.choice([1, 1, 2, 3, 5])
        member_ids = sorted(next(ids) for _ in range(size))
        rep = member_ids[0]
        biomes = f"root;{GUT}:Large intestine" if c % 2 else f"root;{MARINE}"
        clusters.append((rep, size, size, biomes, biomes))
        for m in member_ids:
            seq = "".join(rng.choice(AMINO_ACIDS) if rng.random() < 0.05 else a for a in base)
            members.append((rep, m))
            proteins.append((m, rng.random() < 0.7, seq))
            if c % 3:
                pfam.append((m, 1000 + c, 1e-10, 50.0, 1, 30, 2, 31))
    tables = {
        "mgy_clusters": (
            clusters,
            [
                "cluster_rep",
                "cluster_size",
                "cluster_assembly_count",
                "cluster_rep_biomes",
                "cluster_members_biomes",
            ],
        ),
        "mgy_cluster_seqs": (members, ["cluster_rep", "cluster_member"]),
        "mgy_protein_sequences": (proteins, ["protein_id", "full_length", "sequence"]),
        "mgy_proteins_pfam": (
            pfam,
            [
                "protein_id",
                "pfam_accession",
                "i_evalue",
                "score",
                "hmm_from",
                "hmm_to",
                "env_from",
                "env_to",
            ],
        ),
    }
    for name, (rows, schema) in tables.items():
        frame = pl.DataFrame(rows, schema=schema, orient="row")
        frame.sort(schema[0]).write_parquet(out / f"{name}.parquet")


def mini_fmh(rng: random.Random) -> None:
    """Three one-contig genomes of 25 genes each on both strands, 12 KOs, 10% genes unlabelled."""
    out = OUT / "mini_fmh"
    genomes_dir = out / "genomes_extracted_from_kegg"
    kos = [f"ko:K{i:05d}" for i in range(1, 13)]
    ko_proteins: dict[str, list[str]] = {ko: [] for ko in kos}
    faa, ko_rows = [], []
    for g in ("gaa", "gbb", "gcc"):
        contig, genome, rows = f"CP{rng.randrange(10**6):06d}.1", [], []
        for i in range(25):
            genome.append("".join(rng.choices("ACGT", k=rng.randint(50, 300))))
            protein = "M" + "".join(rng.choices(AMINO_ACIDS, k=rng.randint(100, 300)))
            cds = "".join(rng.choice(SYNONYMS[aa]) for aa in protein) + rng.choice(SYNONYMS["*"])
            strand = rng.choice("+-")
            start = sum(map(len, genome)) + 1
            genome.append(cds if strand == "+" else reverse_complement(cds.encode()).decode())
            gene = f"{g}:G{i:04d}"
            rows.append((g, f"{g}_{contig}", gene, gene, contig, start, start + len(cds) - 1,
                         strand, protein, cds))  # fmt: skip
            faa.append(f">{gene}|{gene}|{g}_{contig}|{contig}|{i}|{len(protein)}\n{protein}\n")
            if rng.random() < 0.9:
                ko = rng.choice(kos)
                ko_rows.append((gene, ko))
                ko_proteins[ko].append(protein)
        (genomes_dir / g).mkdir(parents=True, exist_ok=True)
        seq = "".join(genome) + "".join(rng.choices("ACGT", k=300))
        lines = "\n".join(seq[i : i + 80] for i in range(0, len(seq), 80))
        (genomes_dir / g / f"{g}.fasta").write_text(f"> {contig} Mini genome {g}\n{lines}\n")
        pl.DataFrame(
            rows,
            schema=["genome_name", "assembly_id", "gene_name", "protein_id", "contig_id",
                    "start_position", "end_position", "strand", "aa_sequence", "nt_sequence"],
            orient="row",
        ).with_row_index("").write_csv(genomes_dir / g / f"{g}_mapping.csv")  # fmt: skip
    (out / "protein_ref_db_giant.faa").write_text("".join(faa))
    pl.DataFrame(ko_rows, schema=["gene_id", "ko_id"], orient="row").with_row_index("").write_csv(
        out / "present_genes_and_koids.csv"
    )
    sigs = []
    for ko, proteins in ko_proteins.items():
        mh = sourmash.MinHash(n=0, ksize=11, is_protein=True, scaled=10, track_abundance=True)
        for protein in proteins:
            mh.add_protein(protein)
        sigs.append(sourmash.SourmashSignature(mh, name=ko))
    with (out / "KOs_mini.sig").open("w") as f:
        save_signatures_to_json(sigs, f)
        f.write("\n")


def mini_pfam(rng: random.Random) -> None:
    """Six single-sequence HMMs (``hmmbuild``, HMMER 3.4) from 40-80 aa segments of mini_fmh
    proteins: two from one protein, the rest from five others (genes in several genomes).

    Kept if present (delete the file to rebuild): hmmbuild's last digits differ between
    platforms (macOS arm64 vs Linux x86_64), so the file cannot be byte-reproducible.
    """
    if (OUT / "mini_fmh" / "Pfam-mini.hmm").exists():
        return
    genes = pl.concat(
        pl.read_csv(p, columns=["gene_name", "aa_sequence"])
        for p in sorted((OUT / "mini_fmh" / "genomes_extracted_from_kegg").glob("*/*_mapping.csv"))
    )
    picked = genes.sample(6, seed=rng.randrange(2**31))["aa_sequence"].to_list()
    segments = [(picked[0], 5, 55), (picked[0], 70, 120)]
    for protein in picked[1:5]:
        start = rng.randrange(len(protein) - 80)
        segments.append((protein, start, start + rng.randint(40, 80)))
    with tempfile.TemporaryDirectory() as tmp:
        hmms = []
        for i, (protein, start, end) in enumerate(segments, start=1):
            sto = Path(tmp) / f"{i}.sto"
            sto.write_text(
                f"# STOCKHOLM 1.0\n#=GF ID Mini{i}\n#=GF AC PF9{i:04d}.1\n"
                f"#=GF GA 25.00 25.00;\n#=GF TC 25.00 25.00;\n#=GF NC 24.00 24.00;\n"
                f"seq{i} {protein[start:end]}\n//\n"
            )
            hmm = Path(tmp) / f"{i}.hmm"
            subprocess.run(["hmmbuild", "--amino", "--seed", "1", str(hmm), str(sto)],
                           check=True, capture_output=True)  # fmt: skip
            # hmmbuild stamps the build date; drop it so reruns give the same bytes
            hmms.append("".join(line for line in hmm.read_text().splitlines(keepends=True)
                                if not line.startswith("DATE")))  # fmt: skip
    (OUT / "mini_fmh" / "Pfam-mini.hmm").write_text("".join(hmms))


def mini_host(rng: random.Random) -> None:
    """A stand-in host for the fmh-benchmark test profile: 60 kb of random DNA (two
    chromosomes) gzipped, the extra records (``chrM``, PhiX) plain, and the fixture proteins
    gzipped as a decoy proteome."""
    out = OUT / "mini_fmh"
    chroms = "".join(f">chr{i}\n{''.join(rng.choices('ACGT', k=30_000))}\n" for i in (1, 2))
    extra = f">NC_012920.1 mito\n{''.join(rng.choices('ACGT', k=2_000))}\n" + (
        f">NC_001422.1 phiX\n{''.join(rng.choices('ACGT', k=1_000))}\n"
    )
    for name, text in (
        ("host_mini.fa.gz", chroms),
        ("decoy_mini.faa.gz", (OUT / "proteins.faa").read_text()),
    ):
        with (
            (out / name).open("wb") as f,
            gzip.GzipFile(filename="", mode="wb", fileobj=f, mtime=0) as gz,
        ):
            gz.write(text.encode())
    (out / "host_extra.fa").write_text(extra)


def mini_mgnify(rng: random.Random) -> None:
    """A stand-in MGnify90 release for the mini_fmh genomes: for 24 of their proteins
    (including every one with a Pfam-mini domain), a cluster whose representative is the
    protein mutated to a known identity (100% to 70%) and whose two other members are the
    representative mutated to 97%; plus 6 unrelated clusters. Members and Pfam tables in the
    release's schemas (integer Pfam accessions: each cluster takes its source gene's
    domains, found by ``hmmsearch --cut_ga``). ``mgnify_truth.tsv`` records each
    representative's source gene and identity."""
    out = OUT / "mini_fmh"
    genes = pl.concat(
        pl.read_csv(p, columns=["gene_name", "aa_sequence"])
        for p in sorted((out / "genomes_extracted_from_kegg").glob("*/*_mapping.csv"))
    )
    with tempfile.TemporaryDirectory() as tmp:
        faa = Path(tmp) / "p.faa"
        faa.write_text("".join(f">{g}\n{s}\n" for g, s in genes.iter_rows()))
        tbl = Path(tmp) / "d.tbl"
        subprocess.run(["hmmsearch", "--cut_ga", "--domtblout", tbl, "-o", "/dev/null",
                        out / "Pfam-mini.hmm", faa], check=True)  # fmt: skip
        domains: dict[str, set[int]] = {}
        for line in tbl.read_text().splitlines():
            if not line.startswith("#"):
                f = line.split()
                domains.setdefault(f[0], set()).add(int(f[4][2:].split(".")[0]))
    others = genes.filter(~pl.col("gene_name").is_in(list(domains)))
    chosen = (
        genes.filter(pl.col("gene_name").is_in(list(domains))).rows()
        + others.sample(24 - len(domains), seed=rng.randrange(2**31)).rows()
    )
    identities = [1.0, 0.95, 0.9, 0.85, 0.8, 0.7]

    def mutate(seq: str, identity: float) -> str:
        out_seq = list(seq)
        for i in rng.sample(range(len(seq)), round((1 - identity) * len(seq))):
            out_seq[i] = rng.choice(AMINO_ACIDS.replace(seq[i], ""))
        return "".join(out_seq)

    members, truth, pfam, pid = [], [], [], 1
    for n, (gene, protein) in enumerate(chosen):
        identity = identities[n % len(identities)]
        rep = mutate(protein, identity)
        cluster = [(pid, pid, True, rep)] + [
            (pid + i, pid, True, mutate(rep, 0.97)) for i in (1, 2)
        ]
        members += cluster
        truth.append((pid, gene, identity))
        pfam += [(m[0], acc) for m in cluster for acc in sorted(domains.get(gene, ()))]
        pid += 3
    for _ in range(6):
        rep = "M" + "".join(rng.choices(AMINO_ACIDS, k=rng.randint(150, 300)))
        members += [(pid, pid, True, rep), (pid + 1, pid, True, mutate(rep, 0.97))]
        pid += 2
    schema = ["protein_id", "cluster_rep", "full_length", "sequence"]
    pl.DataFrame(members, schema=schema, orient="row").write_parquet(out / "mgnify_members.parquet")
    pl.DataFrame(truth, schema=["cluster_rep", "gene_name", "identity"], orient="row").write_csv(
        out / "mgnify_truth.tsv", separator="\t"
    )
    pl.DataFrame(pfam, schema=["protein_id", "pfam_accession"], orient="row").write_parquet(
        out / "mgnify_pfam.parquet"
    )


if __name__ == "__main__":
    main()
    mini_release(random.Random(20260929))
    mini_fmh(random.Random(20260930))
    mini_pfam(random.Random(20261001))
    mini_host(random.Random(20261002))
    mini_mgnify(random.Random(20261003))
