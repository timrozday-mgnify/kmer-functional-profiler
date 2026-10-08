"""kfp-genomes: genome profile and functional update, a companion of kmer-functional-profiler
(plan: Genome profile and functional update, phase 12)."""

import typer

app = typer.Typer(no_args_is_help=True)


@app.callback()
def main() -> None:
    """Place the sample's species and strains among reference genomes from the functional
    profile (`place`, on a `panel`), then update the functional profile from them (`update`).
    The profile is never modified."""
