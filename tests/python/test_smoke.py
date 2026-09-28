from typer.testing import CliRunner

import kmer_functional_profiler
from kmer_functional_profiler.cli import app


def test_version_matches_package() -> None:
    result = CliRunner().invoke(app, ["version"])
    assert result.exit_code == 0
    assert result.stdout.strip() == kmer_functional_profiler.__version__ == "0.1.0"
