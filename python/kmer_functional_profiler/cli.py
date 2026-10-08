"""Prototype CLI; mirrors the planned Rust one."""

import json
from dataclasses import fields, replace
from pathlib import Path
from typing import Annotated

import polars as pl
import typer

from kmer_functional_profiler import __version__
from kmer_functional_profiler.compat import import_signatures
from kmer_functional_profiler.genomes import (
    MIN_CONTAINMENT,
    MIN_UNITS,
    GenomeIndex,
    genome_profile,
)
from kmer_functional_profiler.genomes import annotate_genomes as annotate
from kmer_functional_profiler.index import (
    AAI_CALIBRATION,
    AAI_MODEL,
    Index,
    IndexParams,
    build_index,
    write_unit_columns,
)
from kmer_functional_profiler.mask import Mask, build_mask
from kmer_functional_profiler.query import MIN_AAI_KMERS, Timer, check_aai_calibration, profile
from kmer_functional_profiler.species import (
    BG_PRIOR,
    PRIOR_PRESENT,
    SpeciesIndex,
    species_index_from_catalogue,
    species_index_from_genomes,
    species_profile,
)
from kmer_functional_profiler.survival import check_aai_model

app = typer.Typer(no_args_is_help=True)
DEFAULTS = IndexParams()


@app.command()
def version() -> None:
    """Print the version."""
    typer.echo(__version__)


@app.command()
def index(
    members: Annotated[Path, typer.Argument(help="Parquet: protein_id, cluster_rep, ...")],
    out_dir: Path,
    pfam: Annotated[Path | None, typer.Option(help="Parquet: protein_id, pfam_accession")] = None,
    k: int = DEFAULTS.k,
    alphabet: str = DEFAULTS.alphabet,
    t_base: float = DEFAULTS.t_base,
    t_base_singleton: Annotated[
        float | None, typer.Option(help="t_base for singleton clusters (default: t_base)")
    ] = None,
    n_min: int = DEFAULTS.n_min,
    t_cap: float = DEFAULTS.t_cap,
    oversample: float = DEFAULTS.oversample,
    mask_adapters: bool = DEFAULTS.mask_adapters,
    max_groups: int = DEFAULTS.max_groups,
    fp_bits: int = DEFAULTS.fp_bits,
    t_dense: float = DEFAULTS.t_dense,
    batch_residues: int = DEFAULTS.batch_residues,
    postings: Annotated[bool, typer.Option(help="Also write postings.parquet")] = False,
    like: Annotated[
        Path | None,
        typer.Option(
            help="An index whose build parameters to copy (all but --batch-residues), so "
            "the new one can be queried jointly with it"
        ),
    ] = None,
    role: Annotated[
        str | None,
        typer.Option(
            help="decoy: when queried as an --extra-index, its units compete but are "
            "reported as one row (host or contaminant proteomes)"
        ),
    ] = None,
) -> None:
    """Build an index from a members table; print its stats as JSON."""
    if role not in (None, "decoy"):
        raise typer.BadParameter("role must be decoy")
    if like is not None:
        base = json.loads((like / "meta.json").read_text())
        if base.get("hash") != "kfp":
            raise typer.BadParameter("--like needs an index built by `index` (kfp hash)")
    params = (
        IndexParams(
            **{
                f.name: base["params"][f.name]
                for f in fields(IndexParams)
                if f.name in base["params"]
            }
        )
        if like is not None
        else IndexParams(
            k=k,
            alphabet=alphabet,
            t_base=t_base,
            t_base_singleton=t_base_singleton,
            n_min=n_min,
            t_cap=t_cap,
            oversample=oversample,
            mask_adapters=mask_adapters,
            max_groups=max_groups,
            fp_bits=fp_bits,
            t_dense=t_dense,
        )
    )
    params = replace(params, batch_residues=batch_residues)
    stats = build_index(members, out_dir, params, pfam, postings)
    if role is not None:
        meta = json.loads((out_dir / "meta.json").read_text())
        (out_dir / "meta.json").write_text(json.dumps(meta | {"role": role}, indent=2) + "\n")
    typer.echo(json.dumps(stats, indent=2))


@app.command()
def import_sourmash(
    signatures: Annotated[Path, typer.Argument(help="sourmash signatures, e.g. .sig.zip")],
    out_dir: Path,
    ksize: int = 11,
) -> None:
    """Build an index from protein FracMinHash signatures (queried with sourmash hashing)."""
    typer.echo(json.dumps(import_signatures(signatures, out_dir, ksize), indent=2))


@app.command()
def unit_columns(index_dir: Path) -> None:
    """Write the unit table's numeric columns as .npy files, for an index built without them."""
    write_unit_columns(index_dir)


@app.command()
def calibrate_aai(
    index_dir: Path,
    calibration: Annotated[
        Path | None,
        typer.Argument(help="JSON from the benchmark's aai-calibrate; omit to remove the map"),
    ] = None,
) -> None:
    """Attach an aai calibration to an index (``aai_calibration.json``): its queries then
    report ``aai`` as alignment identity and keep the raw estimate as ``aai_raw``."""
    target = index_dir / AAI_CALIBRATION
    if calibration is None:
        target.unlink(missing_ok=True)
        return
    cal = json.loads(calibration.read_text())
    meta = json.loads((index_dir / "meta.json").read_text())
    try:
        check_aai_calibration(cal, meta["params"])
    except ValueError as e:
        raise typer.BadParameter(str(e)) from e
    target.write_text(json.dumps(cal, indent=1) + "\n")


@app.command()
def aai_model(
    index_dir: Path,
    model: Annotated[
        Path | None,
        typer.Argument(help="JSON from the aai-model workflow; omit to remove the model"),
    ] = None,
) -> None:
    """Attach a k-mer survival model to an index (``aai_model.json``): its queries then
    estimate ``aai`` under it (Markov-modulated regions, ``survival.SurvivalModel``)
    instead of a^k."""
    target = index_dir / AAI_MODEL
    if model is None:
        target.unlink(missing_ok=True)
        return
    fitted = json.loads(model.read_text())
    meta = json.loads((index_dir / "meta.json").read_text())
    try:
        check_aai_model(fitted, meta["params"])
    except ValueError as e:
        raise typer.BadParameter(str(e)) from e
    target.write_text(json.dumps(fitted, indent=1) + "\n")


@app.command(name="mask")
def mask_command(
    genome: Annotated[Path, typer.Argument(help="Host genome FASTA (plain or gzip)")],
    index_dir: Path,
    out_dir: Path,
) -> None:
    """Build a mask sidecar: the index's k-mers in the six-frame translated genome."""
    typer.echo(json.dumps(build_mask(genome, index_dir, out_dir), indent=2))


@app.command()
def annotate_genomes(
    index_dir: Path,
    genomes: Annotated[
        Path,
        typer.Argument(
            help="TSV: genome (name), path (protein FASTA, relative to the TSV), optional "
            "taxonomy (GTDB-style d__;p__;...;s__)"
        ),
    ],
    out_dir: Path,
) -> None:
    """Annotate reference genomes with the index: each genome's raw tier-2 hits per unit,
    the content `genomes` fits a profile with."""
    typer.echo(json.dumps(annotate(index_dir, genomes, out_dir), indent=2))


@app.command(name="genomes")
def genomes_command(
    profile_tsv: Annotated[Path, typer.Argument(help="profile.tsv from `query`")],
    genome_index: Annotated[Path, typer.Argument(help="Output of `annotate-genomes`")],
    out: Annotated[Path, typer.Argument(help="TSV of detected genomes")] = Path("genomes.tsv"),
    summary: Annotated[
        Path | None, typer.Option(help="JSON: explained fraction, genome-equivalents, counts")
    ] = None,
    function_taxon: Annotated[
        Path,
        typer.Option(
            help="TSV of each function's split hits (hits_em) per taxon and rank, with "
            "unclassified and the total (Pfam if the genome index has labels, else unit name)"
        ),
    ] = Path("function_taxon.tsv"),
    index: Annotated[
        Path | None,
        typer.Option(help="The profile's index, checked against the one annotated with"),
    ] = None,
    min_containment: float = MIN_CONTAINMENT,
    min_units: int = MIN_UNITS,
) -> None:
    """Genomes present and their depths, from a profile's per-unit hits."""
    gi = GenomeIndex(genome_index)
    if index is not None:
        try:
            gi.check_index(index)
        except ValueError as e:
            raise typer.BadParameter(str(e)) from e
    prof = pl.read_csv(profile_tsv, separator="\t")
    table, sample, stratified = genome_profile(
        prof, gi, min_containment=min_containment, min_units=min_units
    )
    if stratified is None:
        typer.echo("profile has no hits_em: no function x taxon table", err=True)
    else:
        stratified.write_csv(function_taxon, separator="\t")
    table.write_csv(out, separator="\t")
    if summary is not None:
        summary.write_text(json.dumps(sample, indent=2) + "\n")
    typer.echo(f"{table.height} genomes detected -> {out}", err=True)


@app.command()
def species_index(
    index_dir: Path,
    out_dir: Path,
    catalogue: Annotated[
        Path | None,
        typer.Option(
            help="MGnify genome catalogue directory (genomes-all_metadata.tsv and "
            "species_catalogue/ as on the FTP site): pangenome families annotated, carriage "
            "from gene_presence_absence.Rtab"
        ),
    ] = None,
    genomes: Annotated[
        Path | None,
        typer.Option(help="Output of `annotate-genomes` with taxonomy (species from s__), instead"),
    ] = None,
    exclude: Annotated[
        Path | None,
        typer.Option(
            help="Genome names (one per line) left out of every count; their carried units "
            "go to held_out.parquet (benchmark truth)"
        ),
    ] = None,
    species: Annotated[
        Path | None,
        typer.Option(help="With --catalogue: species representatives (one per line) to keep"),
    ] = None,
    alpha: Annotated[
        float | None,
        typer.Option(help="Fix the shrinkage α at every rank (ablation; ~0: no shrinkage)"),
    ] = None,
    completeness: Annotated[
        bool, typer.Option(help="Weight genomes by completeness (off: an ablation)")
    ] = True,
) -> None:
    """Build a species index: per (species, unit) the prevalence prior and the expected
    hits per genome copy, from a genome catalogue's pangenomes or an annotated genome set."""
    if (catalogue is None) == (genomes is None):
        raise typer.BadParameter("give one of --catalogue and --genomes")
    names = frozenset(exclude.read_text().split()) if exclude is not None else frozenset()
    if catalogue is not None:
        keep = set(species.read_text().split()) if species is not None else None
        meta = species_index_from_catalogue(
            index_dir, catalogue, out_dir, names, keep, alpha, completeness
        )
    else:
        assert genomes is not None
        meta = species_index_from_genomes(index_dir, genomes, out_dir, names, alpha, completeness)
    typer.echo(json.dumps(meta, indent=2))


@app.command(name="species")
def species_command(
    profile_tsv: Annotated[Path, typer.Argument(help="profile.tsv from `query`")],
    species_index: Annotated[Path, typer.Argument(help="Output of `species-index`")],
    out: Annotated[Path, typer.Argument(help="TSV of detected species")] = Path("species.tsv"),
    units: Annotated[
        Path, typer.Option(help="TSV of the detected species' units: prevalence, carriage")
    ] = Path("species_units.tsv"),
    presence: Annotated[
        Path, typer.Option(help="TSV of unit presence updated by the species (kfp-prior format)")
    ] = Path("presence.tsv"),
    pfam_presence: Annotated[
        Path, typer.Option(help="TSV of Pfam presence, observed and updated")
    ] = Path("pfam_presence.tsv"),
    function_taxon: Annotated[
        Path,
        typer.Option(
            help="TSV of each function's split hits (hits_em) per species, genus and family, "
            "with unclassified and the total"
        ),
    ] = Path("function_species.tsv"),
    summary: Annotated[
        Path | None, typer.Option(help="JSON: explained fraction, genome-equivalents, counts")
    ] = None,
    index: Annotated[
        Path | None,
        typer.Option(help="The profile's index, checked against the species index's"),
    ] = None,
    min_containment: float = MIN_CONTAINMENT,
    min_units: int = MIN_UNITS,
    prior: Annotated[
        float, typer.Option(help="Prior probability that a screened species is present")
    ] = PRIOR_PRESENT,
    background_prior: Annotated[
        float,
        typer.Option(help="Prior probability that a unit has background hits (0: none)"),
    ] = BG_PRIOR,
) -> None:
    """Species present and their depths, fitted jointly with the units each carries in the
    sample (prevalence as prior, the profile's hits as evidence); updated unit and Pfam
    presence; the function x species table. The profile is not modified."""
    si = SpeciesIndex(species_index)
    if index is not None:
        try:
            si.check_index(index)
        except ValueError as e:
            raise typer.BadParameter(str(e)) from e
    prof = pl.read_csv(profile_tsv, separator="\t")
    result = species_profile(
        prof, si, min_containment=min_containment, min_units=min_units, prior=prior,
        background_prior=background_prior,
    )  # fmt: skip
    outputs = {"species": out, "units": units, "presence": presence,
               "pfam_presence": pfam_presence, "function_species": function_taxon}  # fmt: skip
    for key, path in outputs.items():
        table = result[key]
        if isinstance(table, pl.DataFrame):
            table.write_csv(path, separator="\t")
    if summary is not None:
        summary.write_text(json.dumps(result["summary"], indent=2) + "\n")
    detected = result["species"]
    assert isinstance(detected, pl.DataFrame)
    typer.echo(f"{detected.height} species detected -> {out}", err=True)


@app.command()
def lineage(
    profile_tsv: Annotated[Path, typer.Argument(help="profile.tsv from `query`")],
    species_index: Annotated[Path, typer.Argument(help="Output of `species-index`")],
    out: Annotated[Path, typer.Argument(help="TSV of lineages")] = Path("lineages.tsv"),
    units: Annotated[
        Path, typer.Option(help="TSV of each lineage's units: prevalence, carriage, hits")
    ] = Path("lineage_units.tsv"),
    summary: Annotated[Path | None, typer.Option(help="JSON: the species fit's report")] = None,
    lineages: Annotated[int, typer.Option(help="Lineages fitted per species (K)")] = 2,
    dim: Annotated[int, typer.Option(help="Dimensions of the lineage coordinates")] = 2,
    starts: Annotated[int, typer.Option(help="MAP starts per species")] = 4,
    seed: int = 0,
    min_containment: float = MIN_CONTAINMENT,
    min_units: int = MIN_UNITS,
    prior: Annotated[
        float, typer.Option(help="Prior probability that a screened species is present")
    ] = PRIOR_PRESENT,
    background_prior: Annotated[
        float,
        typer.Option(help="Prior probability that a unit has background hits (0: none)"),
    ] = BG_PRIOR,
) -> None:
    """Prototype (phase 11, step 12): lineages within each detected species, placed among
    its reference genomes, with their depths and carriage. Needs the `phylo` dependency
    group (JAX, NumPyro)."""
    try:
        from kmer_functional_profiler.lineage import lineage_profile
    except ImportError as e:
        raise typer.BadParameter(f"needs the phylo dependency group ({e})") from e
    prof = pl.read_csv(profile_tsv, separator="\t")
    result = lineage_profile(
        prof, SpeciesIndex(species_index), dim=dim, k=lineages, starts=starts, seed=seed,
        prior=prior, background_prior=background_prior, min_containment=min_containment,
        min_units=min_units,
    )  # fmt: skip
    result["lineages"].write_csv(out, separator="\t")
    result["units"].write_csv(units, separator="\t")
    if summary is not None:
        summary.write_text(json.dumps(result["summary"], indent=2) + "\n")
    typer.echo(f"{result['lineages'].height} lineages -> {out}", err=True)


@app.command()
def query(
    index_dir: Path,
    r1: Path,
    r2: Annotated[Path | None, typer.Argument()] = None,
    out: Annotated[Path, typer.Option(help="TSV of per-unit counts")] = Path("profile.tsv"),
    frames: Annotated[
        str,
        typer.Option(
            help="Frames hashed: stopfree; edges (also the terminal segments, >= 30 aa, of "
            "frames with stops; edges:M sets the length); all"
        ),
    ] = "edges:30",
    min_qual: Annotated[
        int, typer.Option(help="Mask bases below this Phred quality as N (0 = off)")
    ] = 0,
    genetic_code: int = 11,
    draws: Annotated[
        int, typer.Option(help="Posterior draws for 95% intervals and ambiguity groups")
    ] = 0,
    kmers: Annotated[
        Path | None, typer.Option(help="Parquet of tier-2 hits per unit and k-mer (diagnostics)")
    ] = None,
    own_hist: Annotated[
        Path | None,
        typer.Option(
            help="Parquet of each unit's hit histogram over the k-mers only it holds (unit, "
            "hits, kmers): the strain-mixture evidence of kfp-genomes place"
        ),
    ] = None,
    stats: Annotated[
        Path | None,
        typer.Option(
            help="JSON of time and peak RSS per stage, and hit counts; rewritten after each "
            "stage, and each stage's start and end are logged to stderr with its RSS"
        ),
    ] = None,
    in_memory: Annotated[
        bool,
        typer.Option(
            help="Read the index tiers into memory instead of memory-mapping them: no page "
            "faults during lookup, at the cost of their full size in RSS"
        ),
    ] = False,
    low_memory: Annotated[
        bool,
        typer.Option(
            help="Trade time for memory where results are unchanged: with --draws, re-read "
            "the reads for the detected units' k-mers instead of keeping per-read rows of "
            "every hit unit through the query"
        ),
    ] = False,
    aai: Annotated[
        bool,
        typer.Option(
            help="Containment AAI of the detected units (aai, with a closed-form 95% interval "
            "aai_lo/aai_hi; no --draws needed); fits the zero-inflated EM. aai_naive is "
            "reported for every hit unit regardless"
        ),
    ] = False,
    min_aai: Annotated[
        float,
        typer.Option(help="Drop units with aai_naive below this, or null (0 = keep all)"),
    ] = 0.0,
    min_aai_kmers: Annotated[
        float,
        typer.Option(
            help="aai, aai_lo, aai_hi (and aai_naive) null below this many hit k-mers "
            "(aai_kmers, kmers_hit), which are still reported"
        ),
    ] = MIN_AAI_KMERS,
    extra_index: Annotated[
        list[Path] | None,
        typer.Option(
            help="Another index (same k, alphabet, hash) queried jointly; repeatable. Its "
            "units compete with the first index's; output ids are offset, with a source column"
        ),
    ] = None,
    mask: Annotated[
        Path | None,
        typer.Option(help="Mask sidecar (from `mask`) of the first index: host k-mers dropped"),
    ] = None,
    summary: Annotated[
        Path | None,
        typer.Option(
            help="JSON of sample-level explained and unknown fractions: model-based "
            "(explained_fraction) and a census of known k-mers (census_containment); fits "
            "the zero-inflated EM"
        ),
    ] = None,
    all_estimators: Annotated[
        bool,
        typer.Option(
            help="Also fit the benchmark estimators (zero-inflated, empirical-Bayes, p_in-"
            "weighted EM; winner-take-all, uniqueness-first) and the experimental strain-mix "
            "coverage_mix; coverage_zi is the default"
        ),
    ] = False,
    em_start: Annotated[
        int,
        typer.Option(
            help="Seed for random EM starting points (0 = the usual start): profiles from "
            "several seeds agree where the EM converged, a check needing no truth"
        ),
    ] = 0,
) -> None:
    """Profile reads (FASTA/FASTQ, optionally paired) against an index."""
    timer = Timer(stats, log=stats is not None)
    sample: dict[str, float | int | None] | None = None if summary is None else {}
    with timer("total"):
        with timer("load"):
            loaded = Index.load(index_dir, mmap=not in_memory)
            extra = [Index.load(d, mmap=not in_memory) for d in extra_index or []]
            masked = None if mask is None else Mask(mask)
        result = profile(
            loaded,
            r1,
            r2,
            genetic_code=genetic_code,
            frames=frames,
            min_qual=min_qual,
            draws=draws,
            kmers_out=kmers,
            own_hist_out=own_hist,
            timer=timer if stats else None,
            all_estimators=all_estimators,
            low_memory=low_memory,
            with_aai=aai,
            min_aai=min_aai,
            min_aai_kmers=min_aai_kmers,
            extra=extra,
            mask=masked,
            summary=sample,
            em_start=em_start,
        )
        result.write_csv(out, separator="\t")
        if summary is not None:
            summary.write_text(json.dumps(sample, indent=2) + "\n")
    timer.write()
    typer.echo(f"{result.height} units hit -> {out}", err=True)
