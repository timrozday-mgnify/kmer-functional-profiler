"""Spike-in benchmark: exact truth for nucleotide mode, the pileups and ``aai``/``ani``
(plan: phase 13, steps 3-6 and its gate; phase 14, step 5 and its gate).

Subcommands (run by ``main.nf``; each prints what it wrote):

- ``families``: targets with MGnify-like structure from seed proteins: per seed, ``--paralogs``
  units (the seed mutated to ``--paralog-identity``), each with ``--members`` members (its
  centroid mutated to ``--member-identity``, as in a 90% cluster). For a nucleotide molecule
  each seed is back-translated once and mutated at nucleotide level. Writes ``targets.fa``
  and ``groups.tsv`` (member, unit; ``index --groups``).
- ``siblings``: ``--siblings`` of the units get a sibling unit ``sib_<unit>``, one of its
  members mutated to ``--sibling-identity`` (MGnify90 neighbours are below 90%). Siblings
  join the index (so the profile sees them, and its EM credits their k-mers to them) but are
  *withheld* from the pileup's reference (``withhold.txt``), the plan's test of the decoys:
  a withheld sibling's reads land on its relatives, which are neighbour decoys unless the
  profile called them present. Writes ``index_targets.fa``, ``index_groups.tsv`` (every
  target's unit), ``siblings.tsv`` (sibling, parent unit) and ``withhold.txt``.
- ``simulate``: one sample. ``--present`` of the (non-sibling) units carry a strain: one of
  its members mutated to an identity drawn from ``--identities``, at a lognormal depth;
  each sibling carries one with probability ``--sibling-present``. Protein strains are
  back-translated with random codons; every strain sits between random flanks. Paired
  reads (``--read-len``, fragments ``--fragment`` +- 10%, either strand) get ``--error``
  substitutions; reads are named ``s<strain>_<n>``, background pairs (``--background``,
  appended) ``bg<n>``. Writes ``reads_1.fq.gz``, ``reads_2.fq.gz`` and ``truth.tsv``.
- ``translate``: a nucleotide targets FASTA translated (code 11, internal stops as X), for the
  amino-acid arm on the same targets (phase 14, step 5).
- ``score``: one sample's arms against its truth -> ``scores.tsv`` (one row per arm) and
  ``units.tsv`` (per unit, every arm's call, abundance and similarity beside the truth).
- ``summary``: means over samples per arm -> ``summary.tsv``.

Truth: a unit is present when it carries a strain (a sibling's included: it is in the
index); its depth is the strain's read bases in the CDS over the CDS length; its identity is
the strain's to the nearest member of the unit (the convention ``aai`` is fitted at), and
``identity_bin`` the ladder value drawn.

Scores per arm (``kfp``: the profile, ``present_prob`` >= 0.5; ``kfp_aa``: the profile
against the translated targets; ``kma``: KMA's ``.res``, similarity its ``Query_Identity``,
over the aligned part):

- ``purity``, ``completeness``, ``f1``; ``l1`` between relative abundances over the
  union of true and called units; ``spearman_tp`` over true positives;
- ``sim_bias_<bin>``, ``sim_cover_<bin>``, ``sim_n_<bin>``: median similarity (``aai``,
  ``ani``, KMA's query identity) minus truth by identity bin, and the share of truths
  inside ``*_lo``-``*_hi``; ``sim_sharp_bias_<bin>``, ``sim_sharp_n_<bin>``: the same bias over
  intervals at most ``INFORMATIVE`` wide (a degenerate interval marks an unidentified ``aai``);
- ``kfp`` only, from the pileup: ``read_precision`` (posterior-weighted share of target
  assignments that go to the read's own unit, over reads from strains), ``sibling_caught``
  (share of withheld siblings' mate weight not assigned to a target; the gate wants >= 0.8),
  ``sibling_caught_parent_target``, ``_neighbour`` and ``_none`` (split by the role the
  profile gave the sibling's parent unit: a target cannot catch them),
  ``background_assigned`` (mates of background reads kept on a target);
- ``kfp`` only, step 6's ``aai`` check: ``aln_gap_<bin>`` (median ``aai`` minus the reads'
  alignment identity, units.tsv's ``identity``), ``aln_truth_<bin>`` (median alignment
  identity minus truth) and ``aln_cover`` (share of alignment identities inside the ``aai``
  interval). Step 6's rule flags a gap above 0.01 in a bin with >= 100 units, or a cover
  below 0.8.
- ``kfp`` only, step 6b (amino-acid indexes): ``consensus_identity_<bin>``, the median
  identity of a unit's consensus coding sequence to its strain's CDS over covered columns
  (``consensus_n_<bin>`` strains).
"""

import argparse
import gzip
import math
import random
from collections.abc import Iterator
from pathlib import Path
from typing import Final

import numpy as np
import polars as pl

from kmer_functional_profiler.index import read_fasta
from kmer_functional_profiler.reference import BASES, CODE_11, reverse_complement

AMINO: Final = "ACDEFGHIKLMNPQRSTVWY"
CODONS: Final = [a + b + c for a in BASES for b in BASES for c in BASES]
SYNONYMS: Final = {
    aa: [c for c, t in zip(CODONS, CODE_11, strict=True) if t == aa] for aa in AMINO + "*"
}
TRANSLATE: Final = dict(zip(CODONS, CODE_11, strict=True))
FLANK: Final = 150
PRESENT_PROB: Final = 0.5
INFORMATIVE: Final = 0.1  # sim_sharp_*: intervals at most this wide


def mutate(
    seq: str,
    identity: float,
    rng: random.Random,
    alphabet: str,
    weights: dict[str, tuple[str, list[float]]] | None = None,
) -> str:
    """Substitute a (1 - identity) share of the sites, uniformly over sites. The new residue
    is uniform over the alphabet, or drawn by ``weights`` (residue -> (choices, weights);
    see :func:`blosum_weights`). Nucleotide sequences are CDS in frame: a substitution that
    would make a stop codon takes another base (or none), as real genes carry no internal
    stops."""
    out = list(seq)
    for i in rng.sample(range(len(seq)), round((1 - identity) * len(seq))):
        if weights and seq[i] in weights:
            out[i] = rng.choices(*weights[seq[i]])[0]
            continue
        choices = alphabet.replace(seq[i], "")
        if alphabet == BASES:
            c = i - i % 3
            choices = "".join(
                b
                for b in choices
                if TRANSLATE.get("".join(out[c:i] + [b] + out[i + 1 : c + 3])) != "*"
            )
        if choices:
            out[i] = rng.choice(choices)
    return "".join(out)


def blosum_weights() -> dict[str, tuple[str, list[float]]]:
    """Per residue, the other residues weighted 2^(S/2) by BLOSUM62 (half-bit scores, so
    the target-to-background odds): substitutions as conservative as real homologs'."""
    from Bio.Align import substitution_matrices

    m = substitution_matrices.load("BLOSUM62")
    out = {}
    for a in AMINO:
        others = AMINO.replace(a, "")
        out[a] = (others, [2.0 ** (m[a][b] / 2) for b in others])
    return out


SENSE: Final = [c for c, t in zip(CODONS, CODE_11, strict=True) if t != "*"]


def back_translate(protein: str, rng: random.Random) -> str:
    """A CDS for the protein; a residue outside the 20 (X, B, Z, J, U, O in real units) takes
    any sense codon."""
    return "".join(rng.choice(SYNONYMS.get(aa, SENSE)) for aa in protein) + rng.choice(
        SYNONYMS["*"]
    )


def translate(cds: str) -> str:
    """Code 11; a final stop dropped, internal ones (synthetic substitutions make them) as X."""
    protein = "".join(TRANSLATE.get(cds[i : i + 3], "X") for i in range(0, len(cds) - 2, 3))
    return protein.removesuffix("*").replace("*", "X")


def identity(a: str, b: str) -> float:
    return sum(x == y for x, y in zip(a, b, strict=True)) / len(a)


def read_groups(path: Path | None) -> dict[str, str]:
    if path is None:
        return {}
    return dict(line.split("\t")[:2] for line in path.read_text().splitlines() if line)


def families(args: argparse.Namespace) -> None:
    rng = random.Random(args.seed)
    nucleotide = args.molecule == "nucleotide"
    alphabet = BASES if nucleotide else AMINO
    seeds = [
        s for _, s in read_fasta(args.seed_fasta) if 200 <= len(s) <= 600 and set(s) <= set(AMINO)
    ][: args.families]
    targets, groups = [], []
    for f, seed in enumerate(seeds):
        base = back_translate(seed, rng)[:-3] if nucleotide else seed
        for p in range(args.paralogs):
            centroid = mutate(base, rng.uniform(*args.paralog_identity), rng, alphabet)
            unit = f"f{f}p{p}"
            for m in range(args.members):
                targets.append(
                    (f"{unit}m{m}", mutate(centroid, args.member_identity, rng, alphabet))
                )
                groups.append((f"{unit}m{m}", unit))
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    (out / "targets.fa").write_text("".join(f">{i}\n{s}\n" for i, s in targets))
    (out / "groups.tsv").write_text("".join(f"{t}\t{u}\n" for t, u in groups))
    print(f"{len(seeds)} families, {len(targets)} targets -> {out}")


def siblings(args: argparse.Namespace) -> None:
    rng = random.Random(args.seed)
    alphabet = BASES if args.molecule == "nucleotide" else AMINO
    targets = dict(read_fasta(args.targets))
    groups = read_groups(args.groups)
    units: dict[str, list[str]] = {}
    for t in targets:
        units.setdefault(groups.get(t, t), []).append(t)
    picked = rng.sample(sorted(units), round(args.siblings * len(units)))
    sibs = {}
    for unit in picked:
        parent = targets[rng.choice(units[unit])]
        sibs[f"sib_{unit}"] = (
            unit,
            mutate(parent, rng.uniform(*args.sibling_identity), rng, alphabet),
        )
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    records = list(targets.items()) + [(s, seq) for s, (_, seq) in sibs.items()]
    (out / "index_targets.fa").write_text("".join(f">{i}\n{q}\n" for i, q in records))
    (out / "index_groups.tsv").write_text(
        "".join(f"{t}\t{groups.get(t, t)}\n" for t in targets)
        + "".join(f"{s}\t{s}\n" for s in sibs)
    )
    (out / "siblings.tsv").write_text(
        "sibling\tparent\n" + "".join(f"{s}\t{u}\n" for s, (u, _) in sibs.items())
    )
    (out / "withhold.txt").write_text("".join(f"{s}\n" for s in sibs))
    print(f"{len(sibs)} siblings of {len(units)} units -> {out}")


def _fastq_pairs(r1: Path, r2: Path) -> Iterator[tuple[str, str, str, str]]:
    with gzip.open(r1, "rt") as f1, gzip.open(r2, "rt") as f2:
        while True:
            a, b = [f1.readline() for _ in range(4)], [f2.readline() for _ in range(4)]
            if not a[0] or not b[0]:
                return
            yield a[1].strip(), a[3].strip(), b[1].strip(), b[3].strip()


def simulate(args: argparse.Namespace) -> None:
    rng = random.Random(args.seed)
    nprng = np.random.default_rng(args.seed)
    nucleotide = args.molecule == "nucleotide"
    alphabet = BASES if nucleotide else AMINO
    targets = dict(read_fasta(args.targets))
    groups = read_groups(args.groups)
    units: dict[str, list[str]] = {}
    for t in targets:
        units.setdefault(groups.get(t, t), []).append(t)
    sibling_parent = (
        dict(line.split("\t") for line in args.siblings.read_text().splitlines()[1:] if line)
        if args.siblings
        else {}
    )
    names = sorted(u for u in units if u not in sibling_parent)
    present = set(rng.sample(names, round(args.present * len(names))))
    present |= {s for s in sorted(sibling_parent) if rng.random() < args.sibling_present}

    weights = blosum_weights() if args.substitution == "blosum62" and not nucleotide else None
    strains = []  # (kind, unit, parent unit, member, identity_bin, sequence)
    for unit in sorted(present):
        member = rng.choice(units[unit])
        ident = rng.choice(args.identities)
        kind = "sibling" if unit in sibling_parent else "strain"
        seq = mutate(targets[member], ident, rng, alphabet, weights)
        strains.append((kind, unit, sibling_parent.get(unit, unit), member, ident, seq))

    qual = chr(33 + round(-10 * math.log10(max(args.error, 1e-4))))
    rows = []
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    with (
        gzip.open(out / "reads_1.fq.gz", "wt", compresslevel=1) as f1,
        gzip.open(out / "reads_2.fq.gz", "wt", compresslevel=1) as f2,
    ):

        def write(name: str, r1: str, r2: str) -> None:
            for f, r in ((f1, r1), (f2, r2)):
                read = bytearray(r.encode())
                for i in nprng.choice(
                    len(read), nprng.binomial(len(read), args.error), replace=False
                ):
                    read[i] = ord(rng.choice(BASES.replace(chr(read[i]), "")))
                f.write(f"@{name}\n{read.decode()}\n+\n{qual * len(read)}\n")

        for k, (kind, unit, parent, member, ident, seq) in enumerate(strains):
            cds = seq if nucleotide else back_translate(seq, rng)
            flank = lambda: "".join(rng.choices(BASES, k=FLANK))  # noqa: E731
            dna = flank() + cds + flank()
            depth = min(math.exp(rng.gauss(args.depth_mu, args.depth_sigma)), args.max_depth)
            n_pairs = round(depth * len(dna) / (2 * args.read_len))
            cds_bases = 0
            for i in range(n_pairs):
                frag = max(args.read_len, round(rng.gauss(args.fragment, args.fragment / 10)))
                frag = min(frag, len(dna))
                s = rng.randrange(len(dna) - frag + 1)
                for a in (s, s + frag - args.read_len):
                    cds_bases += max(0, min(a + args.read_len, FLANK + len(cds)) - max(a, FLANK))
                fwd = dna[s : s + frag]
                if rng.random() < 0.5:
                    fwd = reverse_complement(fwd.encode()).decode()
                r2 = reverse_complement(fwd[-args.read_len :].encode()).decode()
                write(f"s{k}_{i}", fwd[: args.read_len], r2)
            nearest = max(
                identity(seq, targets[m]) for m in units[unit] if len(targets[m]) == len(seq)
            )
            rows.append(
                (k, kind, unit, parent, member, ident, nearest, depth, len(cds), cds_bases,
                 n_pairs, cds)
            )  # fmt: skip
        n_bg = 0
        if args.background:
            for n_bg, (a, _, b, _) in enumerate(_fastq_pairs(*args.background), 1):
                if n_bg > args.background_pairs > 0:
                    n_bg -= 1
                    break
                f1.write(f"@bg{n_bg}\n{a}\n+\n{'I' * len(a)}\n")
                f2.write(f"@bg{n_bg}\n{b}\n+\n{'I' * len(b)}\n")
    truth = pl.DataFrame(
        rows,
        schema=[
            "strain",
            "kind",
            "unit",
            "parent",
            "member",
            "identity_bin",
            "identity",
            "depth_drawn",
            "cds_len",
            "cds_bases",
            "pairs",
            "cds",
        ],  # fmt: skip
        orient="row",
    ).with_columns(depth=pl.col("cds_bases") / pl.col("cds_len"))
    truth.write_csv(out / "truth.tsv", separator="\t")
    print(f"{len(strains)} strains, {truth['pairs'].sum()} pairs, {n_bg} background -> {out}")


def translate_targets(args: argparse.Namespace) -> None:
    records = [(i, translate(s)) for i, s in read_fasta(args.targets)]
    Path(args.out).write_text("".join(f">{i}\n{p}\n" for i, p in records if p))
    print(f"{len(records)} targets translated -> {args.out}")


def _l1(truth: pl.DataFrame, called: pl.DataFrame) -> float:
    both = truth.join(called, on="unit", how="full", coalesce=True).fill_null(0.0)
    t, e = both["depth"], both["estimate"]
    if t.sum() == 0 or e.sum() == 0:
        return float("nan")
    return float((t / t.sum() - e / e.sum()).abs().sum())


def _arm_scores(
    truth: pl.DataFrame, called: pl.DataFrame, bins: list[float]
) -> dict[str, float | None]:
    """``called``: unit, estimate, and optionally sim, sim_lo, sim_hi."""
    tp = truth.join(called, on="unit")
    n_true, n_called = truth.height, called.height
    purity = tp.height / n_called if n_called else float("nan")
    completeness = tp.height / n_true if n_true else float("nan")
    out: dict[str, float | None] = {
        "true": n_true,
        "called": n_called,
        "purity": purity,
        "completeness": completeness,
        "f1": 2 * purity * completeness / (purity + completeness)
        if purity + completeness > 0
        else 0.0,
        "l1": _l1(truth.select("unit", "depth"), called.select("unit", "estimate")),
        "spearman_tp": tp.select(pl.corr("depth", "estimate", method="spearman")).item()
        if tp.height > 2
        else None,
    }
    if "sim" in called.columns:
        for b in bins:
            at = tp.filter(pl.col("identity_bin") == b, pl.col("sim").is_not_null())
            out[f"sim_n_{b}"] = at.height
            out[f"sim_bias_{b}"] = (
                at.select((pl.col("sim") - pl.col("identity")).median()).item()
                if at.height
                else None
            )
            if "sim_lo" in at.columns:
                sharp = at.filter(pl.col("sim_hi") - pl.col("sim_lo") <= INFORMATIVE)
                out[f"sim_sharp_n_{b}"] = sharp.height
                out[f"sim_sharp_bias_{b}"] = (
                    sharp.select((pl.col("sim") - pl.col("identity")).median()).item()
                    if sharp.height
                    else None
                )
                out[f"sim_cover_{b}"] = (
                    at.select(
                        pl.col("identity").is_between(pl.col("sim_lo"), pl.col("sim_hi")).mean()
                    ).item()
                    if at.height
                    else None
                )
    return out


def _pileup_scores(
    truth: pl.DataFrame,
    assignments: pl.DataFrame,
    units: pl.DataFrame,
    refs: pl.DataFrame,
    sim: str,
    bins: list[float],
) -> dict[str, float | None]:
    a = assignments.with_columns(
        strain=pl.col("read_name").str.extract(r"^s(\d+)_", 1).cast(pl.Int64),
        background=pl.col("read_name").str.starts_with("bg"),
    ).join(truth.select("strain", "kind", source="unit"), on="strain", how="left")
    on_target = a.filter(pl.col("kept"), pl.col("role") == "target")
    from_strains = on_target.filter(pl.col("kind") == "strain")
    w = from_strains["posterior"].sum()
    out: dict[str, float | None] = {
        "read_precision": from_strains.filter(pl.col("cluster_rep") == pl.col("source"))[
            "posterior"
        ].sum()
        / w
        if w
        else None,
        "background_assigned": on_target.filter("background")
        .select("read_name", "mate")
        .n_unique(),
    }
    # siblings: per mate, the weight a target took (0 when the mate aligned nowhere)
    sib = truth.filter(pl.col("kind") == "sibling")
    taken = (
        on_target.filter(pl.col("kind") == "sibling")
        .group_by("strain")
        .agg(taken=pl.col("posterior").sum())
    )
    # the role the profile gave the sibling's parent unit: target, neighbour or none
    roles = refs.select(parent="cluster_rep", parent_role="role").unique()
    sib = (
        sib.join(taken, on="strain", how="left")
        .join(roles, on="parent", how="left")
        .with_columns(
            mates=2 * pl.col("pairs"),
            taken=pl.col("taken").fill_null(0.0),
            parent_role=pl.col("parent_role").fill_null("none"),
        )
    )
    for name, part in (
        ("", sib),
        *(
            (f"_parent_{r}", sib.filter(pl.col("parent_role") == r))
            for r in ("target", "neighbour", "none")
        ),
    ):
        mates = part["mates"].sum()
        out[f"sibling_caught{name}"] = 1 - part["taken"].sum() / mates if mates else None
    # step 6's aai check: alignment identity against aai and against truth
    u = units.join(
        truth.filter(pl.col("kind") == "strain").select(
            "unit", "identity_bin", truth_identity="identity"
        ),
        left_on="cluster_rep",
        right_on="unit",
    ).filter(pl.col("identity").is_not_null(), pl.col("reads") > 0)
    for b in bins:
        at = u.filter(pl.col("identity_bin") == b)
        out[f"aln_n_{b}"] = at.height
        out[f"aln_gap_{b}"] = (
            at.select((pl.col(sim) - pl.col("identity")).median()).item()
            if at.height and sim in at.columns
            else None
        )
        out[f"aln_truth_{b}"] = (
            at.select((pl.col("identity") - pl.col("truth_identity")).median()).item()
            if at.height
            else None
        )
    with_interval = u.filter(pl.col(f"{sim}_lo").is_not_null())
    out["aln_cover"] = (
        with_interval.select(
            pl.col("identity").is_between(pl.col(f"{sim}_lo"), pl.col(f"{sim}_hi")).mean()
        ).item()
        if with_interval.height
        else None
    )
    if "em_posterior" in a.columns:
        out |= _em_scores(truth, a, sib, units)
    return out


def _em_scores(
    truth: pl.DataFrame, a: pl.DataFrame, sib: pl.DataFrame, units: pl.DataFrame
) -> dict[str, float | None]:
    """Steps 7-9 (the joint EM, consensus updates, haplotypes): the crude pileup scores again
    with ``em_posterior`` (no ``kept`` filter: the EM drops reads itself), the share of
    strain mates given to the novel class, and abundance from ``em_abundance``."""
    on_target = a.filter(pl.col("role") == "target")
    strains = on_target.filter(pl.col("kind") == "strain")
    w = strains["em_posterior"].sum()
    out: dict[str, float | None] = {
        "read_precision_em": strains.filter(pl.col("cluster_rep") == pl.col("source"))[
            "em_posterior"
        ].sum()
        / w
        if w
        else None,
    }
    mates = a.group_by("read_name", "mate").agg(
        pl.col("kind").first(), novel=pl.col("em_novel").first()
    )
    out["strain_novel_em"] = mates.filter(pl.col("kind") == "strain")["novel"].mean()
    taken = (
        on_target.filter(pl.col("kind") == "sibling")
        .group_by("strain")
        .agg(taken_em=pl.col("em_posterior").sum())
    )
    sib = sib.join(taken, on="strain", how="left").with_columns(pl.col("taken_em").fill_null(0.0))
    for name, part in (
        ("", sib),
        *(
            (f"_parent_{r}", sib.filter(pl.col("parent_role") == r))
            for r in ("target", "neighbour", "none")
        ),
    ):
        n = part["mates"].sum()
        out[f"sibling_caught_em{name}"] = 1 - part["taken_em"].sum() / n if n else None
    # abundance over the reference's targets and the strains (withheld siblings have no
    # source here): the EM's against the profile's coverage_zi on the same units
    t = truth.filter(pl.col("kind") == "strain").select("unit", "depth")
    for name, col in (("em", "em_abundance"), ("profile", "coverage_zi")):
        called = units.filter(pl.col("em_pairs").is_not_null()).select(
            unit="cluster_rep", estimate=pl.col(col).fill_null(0.0)
        )
        tp = t.join(called, on="unit")
        out[f"l1_targets_{name}"] = _l1(t, called)
        out[f"spearman_targets_{name}"] = (
            tp.select(pl.corr("depth", "estimate", method="spearman")).item()
            if tp.height > 2
            else None
        )
    return out


def _consensus_scores(
    truth: pl.DataFrame, pileup: Path, bins: list[float]
) -> dict[str, float | None]:
    """Step 6b: the consensus coding sequence's identity to its unit's strain over covered
    (non-N) columns, median by identity bin (``consensus_identity_<bin>``, ``_n_<bin>``);
    strains whose consensus is from a member of another length are left out."""
    cons = dict(read_fasta(pileup / "consensus.fna"))
    info = pl.read_csv(pileup / "consensus.tsv", separator="\t", infer_schema=False)
    seq_of = dict(zip(info["cluster_rep"], info["seq_id"], strict=True))
    rows = []
    for unit, b, cds in (
        truth.filter(pl.col("kind") == "strain").select("unit", "identity_bin", "cds").iter_rows()
    ):
        c = cons.get(seq_of.get(str(unit), ""), "")
        # the strain's CDS carries a stop codon past the protein's columns
        cols = [i for i in range(len(c)) if c[i] != "N"] if len(cds) - len(c) in (0, 3) else []
        if cols:
            rows.append((b, sum(c[i] == cds[i] for i in cols) / len(cols)))
    t = pl.DataFrame(rows, schema=["bin", "identity"], orient="row")
    out: dict[str, float | None] = {}
    for b in bins:
        at = t.filter(pl.col("bin") == b)
        out[f"consensus_n_{b}"] = at.height
        out[f"consensus_identity_{b}"] = (
            at.select(pl.col("identity").median()).item() if at.height else None
        )
    return out


def score(args: argparse.Namespace) -> None:
    truth_all = pl.read_csv(args.truth, separator="\t")
    truth = truth_all.select("unit", "depth", "identity", "identity_bin")
    bins = sorted(truth["identity_bin"].unique().to_list(), reverse=True)
    rows, per_unit = [], [truth]

    def kfp_called(path: Path) -> tuple[pl.DataFrame, str]:
        p = pl.read_csv(path, separator="\t", infer_schema_length=None)
        sim = "ani" if "ani_naive" in p.columns else "aai"
        cols = [sim, f"{sim}_lo", f"{sim}_hi"] if sim in p.columns else []
        called = p.filter(pl.col("present_prob") >= PRESENT_PROB).select(
            unit=pl.col("cluster_rep").cast(pl.String),
            estimate="coverage_zi",
            **dict(zip(["sim", "sim_lo", "sim_hi"], cols, strict=False)),
        )
        return called, sim

    called, sim = kfp_called(args.profile)
    row = {"arm": "kfp", **_arm_scores(truth, called, bins)}
    if args.pileup:
        row |= _pileup_scores(
            truth_all,
            pl.read_parquet(args.pileup / "assignments.parquet").join(
                pl.read_csv(args.ref / "refs.tsv", separator="\t", infer_schema=False).select(
                    "seq_id", "cluster_rep"
                ),
                on="seq_id",
            ),
            pl.read_csv(
                args.pileup / "units.tsv", separator="\t", infer_schema_length=None
            ).with_columns(pl.col("cluster_rep").cast(pl.String)),
            pl.read_csv(args.ref / "refs.tsv", separator="\t", infer_schema=False),
            sim,
            bins,
        )
        if (args.pileup / "consensus.fna").exists():
            row |= _consensus_scores(truth_all, args.pileup, bins)
    rows.append(row)
    per_unit.append(called.rename(lambda c: c if c == "unit" else f"kfp_{c}"))
    if args.aa_profile:
        called, _ = kfp_called(args.aa_profile)
        rows.append(
            {
                "arm": "kfp_aa",
                **_arm_scores(truth, called.drop("sim", "sim_lo", "sim_hi", strict=False), bins),
            }
        )
        per_unit.append(called.rename(lambda c: c if c == "unit" else f"kfp_aa_{c}"))
    if args.kma:
        groups = read_groups(args.groups)
        res = pl.read_csv(args.kma, separator="\t", infer_schema=False)
        res = res.rename({res.columns[0]: "template"}).with_columns(
            pl.col("Depth", "Query_Identity").str.strip_chars().cast(pl.Float64)
        )
        called = (
            res.with_columns(
                target=pl.col("template").str.strip_chars().str.split(" ").list.first()
            )
            .with_columns(unit=pl.col("target").replace(groups))
            .group_by("unit")
            .agg(estimate=pl.col("Depth").sum(), sim=pl.col("Query_Identity").max() / 100)
        )
        rows.append({"arm": "kma", **_arm_scores(truth, called, bins)})
        per_unit.append(called.rename(lambda c: c if c == "unit" else f"kma_{c}"))
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    scores = pl.DataFrame(rows, infer_schema_length=None).with_columns(sample=pl.lit(args.sample))
    scores.write_csv(out / "scores.tsv", separator="\t")
    table = per_unit[0]
    for t in per_unit[1:]:
        table = table.join(t, on="unit", how="full", coalesce=True)
    table.write_csv(out / "units.tsv", separator="\t")
    with pl.Config(tbl_cols=-1, tbl_width_chars=250):
        print(scores)


def summary(args: argparse.Namespace) -> None:
    scores = pl.concat(
        [pl.read_csv(p, separator="\t", infer_schema_length=None) for p in args.scores],
        how="diagonal_relaxed",
    )
    out = (
        scores.group_by("arm", maintain_order=True)
        .agg(pl.len().alias("samples"), pl.exclude("sample").mean())
        .with_columns(pl.selectors.float().round(4))
    )
    out.write_csv(args.out, separator="\t")
    with pl.Config(tbl_cols=-1, tbl_rows=-1, tbl_width_chars=250, fmt_str_lengths=40):
        print(out.transpose(include_header=True, column_names="arm"))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = parser.add_subparsers(dest="command", required=True)
    p = sub.add_parser("families")
    p.add_argument("--seed-fasta", type=Path, required=True)
    p.add_argument("--families", type=int, default=50)
    p.add_argument("--paralogs", type=int, default=3)
    p.add_argument("--members", type=int, default=3)
    p.add_argument("--paralog-identity", type=float, nargs=2, default=(0.8, 0.9))
    p.add_argument("--member-identity", type=float, default=0.97)
    p.add_argument("--molecule", choices=["protein", "nucleotide"], default="protein")
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--out", type=Path, default=Path("."))
    p = sub.add_parser("siblings")
    p.add_argument("--targets", type=Path, required=True)
    p.add_argument("--groups", type=Path)
    p.add_argument("--molecule", choices=["protein", "nucleotide"], default="protein")
    p.add_argument("--siblings", type=float, default=0.1)
    p.add_argument("--sibling-identity", type=float, nargs=2, default=(0.8, 0.9))
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--out", type=Path, default=Path("."))
    p = sub.add_parser("simulate")
    p.add_argument("--targets", type=Path, required=True)
    p.add_argument("--groups", type=Path)
    p.add_argument("--molecule", choices=["protein", "nucleotide"], default="protein")
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--present", type=float, default=0.3)
    p.add_argument("--identities", type=float, nargs="+", default=[1.0, 0.98, 0.95, 0.9])
    p.add_argument("--siblings", type=Path, help="siblings.tsv of `siblings`")
    p.add_argument("--sibling-present", type=float, default=1.0)
    p.add_argument(
        "--substitution", choices=["uniform", "blosum62"], default="uniform",
        help="protein strains: new residues uniform, or weighted by BLOSUM62",
    )  # fmt: skip
    p.add_argument("--depth-mu", type=float, default=1.5)
    p.add_argument("--depth-sigma", type=float, default=1.0)
    p.add_argument("--max-depth", type=float, default=100.0)
    p.add_argument("--read-len", type=int, default=150)
    p.add_argument("--fragment", type=int, default=350)
    p.add_argument("--error", type=float, default=0.002)
    p.add_argument("--background", type=Path, nargs=2, metavar=("R1", "R2"))
    p.add_argument("--background-pairs", type=int, default=0, help="0 = all")
    p.add_argument("--out", type=Path, default=Path("."))
    p = sub.add_parser("translate")
    p.add_argument("--targets", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    p = sub.add_parser("score")
    p.add_argument("--sample", required=True)
    p.add_argument("--truth", type=Path, required=True)
    p.add_argument("--profile", type=Path, required=True)
    p.add_argument("--pileup", type=Path)
    p.add_argument("--ref", type=Path)
    p.add_argument("--aa-profile", type=Path)
    p.add_argument("--kma", type=Path)
    p.add_argument("--groups", type=Path)
    p.add_argument("--out", type=Path, default=Path("."))
    p = sub.add_parser("summary")
    p.add_argument("scores", type=Path, nargs="+")
    p.add_argument("--out", type=Path, default=Path("summary.tsv"))
    args = parser.parse_args()
    {
        "families": families,
        "siblings": siblings,
        "simulate": simulate,
        "translate": translate_targets,
        "score": score,
        "summary": summary,
    }[args.command](args)


if __name__ == "__main__":
    main()
