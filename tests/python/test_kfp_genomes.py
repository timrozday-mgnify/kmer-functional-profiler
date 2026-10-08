"""kfp-genomes (phase 12): evidence layer, panel, place and update."""

# ruff: noqa: F811  (index_dir: the fixture imported from test_genomes)

from pathlib import Path
from typing import Any

import numpy as np
import polars as pl
import pytest
from genome_sim import Strain, simulate
from scipy.stats import chi2
from test_genomes import index_dir  # noqa: F401  (fixture)
from typer.testing import CliRunner

from kfp_genomes.cli import app
from kfp_genomes.evidence import (
    Evidence,
    Histogram,
    UnitTerms,
    log_likelihood,
    mixture_loglik,
    zero_hit_present,
)
from kfp_genomes.panel import (
    Panel,
    SpeciesPanel,
    panel_from_catalogue,
    panel_from_genomes,
    write_panel,
)
from kfp_genomes.place import NOVEL, OneStrain, place_species, screen
from kmer_functional_profiler.query import coverage_interval, ztnb_one_rate

# --- Evidence layer (step 1) ---


def fit_terms(m: float, c: float, pi: float) -> tuple[float, float]:
    """k and H of a unit fitted at coverage_zi c and present_zi pi."""
    return pi * m * -np.expm1(-c), pi * m * c


def test_likelihood_peaks_at_coverage_zi_and_present_zi() -> None:
    m, c, pi = 200.0, 2.0, 0.6
    k, h = fit_terms(m, c, pi)
    d = np.linspace(0.5, 8, 1501)[:, None]
    f = np.linspace(0.05, 1, 951)[None, :]
    ll = log_likelihood(k, h, m, d, f)
    i, j = np.unravel_index(np.argmax(ll), ll.shape)
    assert d[i, 0] == pytest.approx(c, abs=0.01) and f[0, j] == pytest.approx(pi, abs=0.002)


def test_profile_interval_over_survival_is_coverage_interval() -> None:
    m, c, pi, w = 200.0, 3.0, 0.3, 2.5
    k, h = fit_terms(m, c, pi)
    d = np.geomspace(0.5, 20, 40001)
    f_hat = np.minimum(1.0, k / (m * -np.expm1(-d)))  # the maximum over f at each D
    drop = log_likelihood(k, h, m, c, pi, w) - log_likelihood(k, h, m, d, f_hat, w)
    inside = d[drop <= chi2.ppf(0.95, 1) / 2]
    lo, hi = coverage_interval(np.array([c]), np.array([k]), np.array([w]))
    assert inside.min() == pytest.approx(lo[0], rel=1e-3)
    assert inside.max() == pytest.approx(hi[0], rel=1e-3)


def test_zero_hit_term_matches_simulation() -> None:
    rng = np.random.default_rng(1)
    m, f, d, n = 40, 0.4, 0.1, 200_000
    present = rng.binomial(m, f, n)
    p_zero = (rng.poisson(d * present) == 0).mean()
    expect = np.exp(log_likelihood(0, 0, m, d, f))
    assert p_zero == pytest.approx(expect, abs=4 * np.sqrt(expect * (1 - expect) / n))
    assert np.exp(log_likelihood(0, 0, m, d, 1.0)) == pytest.approx(np.exp(-d * m))  # f = 1: e^-Dm
    eps = zero_hit_present(np.array([5.0, 500.0]))
    assert 1 > eps[0] > eps[1] > 0  # more k-mers, less chance to be missed when present


def test_histogram_with_one_rate_is_rate_mixtures_one_rate_fit() -> None:
    rng = np.random.default_rng(2)
    rows = []
    for unit, rate in ((7, 3.0), (9, 8.0)):
        hits = rng.poisson(rate, 400)
        hits = hits[hits > 0][:60]
        values, counts = np.unique(hits, return_counts=True)
        rows += [(unit, int(v), int(c)) for v, c in zip(values, counts, strict=True)]
    hist = Histogram.from_frame(pl.DataFrame(rows, schema=["unit", "hits", "kmers"], orient="row"))
    mu, ll = ztnb_one_rate(hist.row, hist.h, hist.c, hist.v)
    assert mixture_loglik(hist, mu[:, None]) == pytest.approx(ll)
    # and the fit's rate maximises it
    for scale in (0.9, 1.1):
        assert (mixture_loglik(hist, scale * mu[:, None]) < ll).all()


def toy_profile() -> pl.DataFrame:
    """Two hit units fitted at coverage 2 (one explained away), the profile's columns."""
    m = np.array([100.0, 80.0, 60.0])
    return pl.DataFrame({
        "unit": [3, 5, 8], "hits": [150, 90, 4], "reads": [100, 60, 4], "m_g": m,
        "coverage_zi": [2.0, 1.5, 0.0], "present_zi": [0.8, 0.9, 0.0],
        "coverage_zi_dispersion": [1.0, 1.2, 1.0], "present_prob": [0.999, 0.99, 0.0],
        "present_llr": [12.0, 10.0, None],
    })  # fmt: skip


def test_unit_ratios_favour_the_fitted_depth_and_doubt_missed_units() -> None:
    ev = Evidence.from_profile(toy_profile())
    assert ev.index(np.array([5, 4, 8])).tolist() == [1, -1, 2]
    assert ev.informative.tolist() == [True, True, False]
    terms = UnitTerms.build(ev, np.array([3, 5, 8, 11]), np.array([100.0, 80.0, 60.0, 50.0]))
    assert terms.hit.tolist() == [True, True, False, False] and terms.blank[2]
    d = np.array([[2.0], [0.2], [2.0], [2.0]])  # unit 3 at its depth and well below
    lr = terms.log_ratio(np.c_[d, d * 5], np.full((4, 2), 0.8), beta=0.01)
    assert lr[0, 0] > 0 > lr[0, 1]  # at its fitted depth, carriage beats not carried
    assert lr[2].tolist() == [0.0, 0.0]  # explained away: no evidence either way
    # a zero-hit unit: carried at depth 2 is unlikely, against no hits from elsewhere
    expect = log_likelihood(0, 0, 50, 2.0, 0.8) - np.log1p(-0.01 * (1 - terms.eps[3]))
    assert lr[3, 0] == pytest.approx(expect) and lr[3, 0] < -30


# --- Panel (step 2) ---

TAXONOMY = "d__B;p__P;c__C;o__O;f__F;g__G;s__"


def synthetic_panel(
    tmp: Path,
    carried: dict[str, set[int]],
    species: dict[str, str],
    completeness: dict[str, float] | None = None,
    f: float = 0.8,
    **kwargs: Any,
) -> Panel:
    """A panel of genomes carrying the units listed (every unit at survival ``f``, one copy),
    each in the species named; index ``meta.json`` faked; 100 kept k-mers per unit."""
    names = list(carried)
    genomes = pl.DataFrame({
        "genome": range(len(names)), "name": names,
        "taxonomy": [TAXONOMY + species[n] for n in names],
        "completeness": [(completeness or {}).get(n, 1.0) for n in names],
    }).cast({"genome": pl.UInt32})  # fmt: skip
    rows = pl.DataFrame(
        [(g, u, 1, f) for g, n in enumerate(names) for u in sorted(carried[n])],
        schema={"genome": pl.UInt32, "unit": pl.UInt32, "n": pl.UInt16, "f": pl.Float64},
        orient="row",
    )
    index_dir = tmp / "index"
    index_dir.mkdir(exist_ok=True)
    (index_dir / "meta.json").write_text("{}")
    m_g = np.full(1 + max(max(c) for c in carried.values()), 100)
    write_panel(tmp / "panel", index_dir, m_g, genomes, rows, **kwargs)
    return Panel(tmp / "panel")


def test_panel_carriage_at_a_reference_and_at_the_species_average(tmp_path: Path) -> None:
    rng = np.random.default_rng(3)
    core = set(range(40))
    carried = {f"g{i}": core | set((40 + rng.choice(60, 20, replace=False)).tolist())
               for i in range(8)}  # fmt: skip
    panel = synthetic_panel(tmp_path, carried, dict.fromkeys(carried, "A"))
    sp = panel.species_panel(0)
    assert sp.carried.shape == (8, sp.units.size) and (sp.edges.shape[0] == 8 * 7)
    for g in range(8):  # a complete reference with t = 0, ell = 0: its own carriage
        assert np.array_equal(sp.carriage(g, sp.edges[g * 7, 1], 0.0, 0.0), sp.carried[g])
    # ell = 1: the rank-shrunk prevalence (species.py's q, from all genomes)
    genomes = panel.genomes.join(pl.DataFrame({"name": list(carried)}), on="name")
    assert genomes.height == 8
    q = dict(panel.species_units.select("unit", "q").iter_rows())
    assert np.allclose(sp.carriage(0, 1, 0.3, 1.0), [q[u] for u in sp.units])
    assert all(q[u] > 0.95 for u in core)
    # half-way between two references: the mean of their carriage
    g, h = sp.edges[0]
    assert np.allclose(sp.carriage(g, h, 0.5, 0.0), (sp.carried[g] * 1.0 + sp.carried[h]) / 2)


def test_incomplete_mag_missing_unit_is_shrunk_by_its_completeness(tmp_path: Path) -> None:
    carried = {f"g{i}": set(range(10)) for i in range(6)} | {"mag": set(range(8))}
    panel = synthetic_panel(tmp_path, carried, dict.fromkeys(carried, "A"), {"mag": 0.6})
    sp = panel.species_panel(0)
    mag = sp.names.index("mag")
    q = sp.q[8]
    assert sp.xt[mag, 8] == pytest.approx(q * 0.4 / (1 - q * 0.6))
    assert sp.xt[mag, 0] == 1.0 and 0 < sp.xt[mag, 8] < q


def test_farthest_point_selection_keeps_both_clades(tmp_path: Path) -> None:
    rng = np.random.default_rng(4)
    core = set(range(50))
    a_set, b_set = set(range(50, 80)), set(range(80, 110))
    carried = {f"a{i}": core | {u for u in a_set if rng.random() < 0.9} for i in range(30)}
    carried |= {f"b{i}": core | {u for u in b_set if rng.random() < 0.9} for i in range(5)}
    panel = synthetic_panel(tmp_path, carried, dict.fromkeys(carried, "A"), max_per_species=3,
                            neighbours=2)  # fmt: skip
    chosen = panel.genomes["name"].to_list()
    assert len(chosen) == 3 and {c[0] for c in chosen} == {"a", "b"}
    assert panel.species["genomes"].to_list() == [35] and panel.species["panel"].to_list() == [3]
    assert panel.neighbours.height == 3 * 2
    # q still counts every genome: clade a's accessory units are more common
    q = dict(panel.species_units.select("unit", "q").iter_rows())
    assert np.mean([q[u] for u in a_set]) > np.mean([q[u] for u in b_set])


def test_panel_from_catalogue_and_genome_set(tmp_path: Path, index_dir: Path) -> None:
    from test_species import A_GENOMES, annotate_genomes_for, fake_catalogue

    panel_from_catalogue(index_dir, fake_catalogue(tmp_path), tmp_path / "cat")
    cat = Panel(tmp_path / "cat")
    cat.check_index(index_dir)
    assert cat.species["id"].to_list() == ["MGYG000000001", "MGYG000000004"]
    # FPS starts at the representative; every genome fits under the cap
    a = cat.genomes.filter(pl.col("species") == 0).sort("fps_rank")
    assert a["name"][0] == "MGYG000000001" and a.height == 3
    rows = cat.carriage_rows.join(cat.genomes.select("genome", "name"), on="genome")
    assert (
        rows.filter(pl.col("name") == "MGYG000000002").height
        < rows.filter(pl.col("name") == "MGYG000000001").height
    )
    assert rows["f"].is_between(0.0, 1.0).all() and (rows["n"] >= 1).all()
    genomes = {g: list(m) for g, m in A_GENOMES.items()} | {"MGYG000000004": list(range(10, 20))}
    gi = annotate_genomes_for(tmp_path, index_dir, genomes)
    panel_from_genomes(index_dir, gi, tmp_path / "set")
    gset = Panel(tmp_path / "set")
    # the same proteins either way: the same carriage and the same q
    key = ["name", "unit"]
    a_rows = rows.select(*key, "n").sort(key)
    b_rows = (gset.carriage_rows.join(gset.genomes.select("genome", "name"), on="genome")
              .select(*key, "n").sort(key))  # fmt: skip
    assert a_rows.equals(b_rows)
    assert np.allclose(cat.species_units.sort("species", "unit")["q"],
                       gset.species_units.sort("species", "unit")["q"])  # fmt: skip


def test_panel_cli(tmp_path: Path, index_dir: Path) -> None:
    from test_species import fake_catalogue

    (tmp_path / "exclude.txt").write_text("MGYG000000002\n")
    args = ["panel", str(index_dir), str(tmp_path / "p"), "--catalogue",
            str(fake_catalogue(tmp_path)), "--exclude", str(tmp_path / "exclude.txt"),
            "--max-per-species", "1"]  # fmt: skip
    result = CliRunner().invoke(app, args)
    assert result.exit_code == 0, result.output
    panel = Panel(tmp_path / "p")
    assert panel.genomes["name"].to_list() == ["MGYG000000001", "MGYG000000004"]
    assert panel.neighbours.height == 0
    assert set(pl.read_parquet(tmp_path / "p" / "held_out.parquet")["name"]) == {"MGYG000000002"}
    assert panel.species_panel(0).edges.tolist() == [[0, 0]]


# --- place, one strain (step 3) ---


def two_clade_species(rng: np.random.Generator, n: int = 6) -> dict[str, set[int]]:
    """Clades a and b of species A: a core of 300 units, each clade 100 accessory units
    carried at 85% by its genomes."""
    core = set(range(300))
    pools = {"a": range(300, 400), "b": range(400, 500)}
    return {f"{c}{i}": core | {u for u in pools[c] if rng.random() < 0.85}
            for c in "ab" for i in range(n)}  # fmt: skip


def placed_on(fit: OneStrain, sp: SpeciesPanel, prefix: str) -> float:
    """Posterior mass of placements whose references are all in the clade ``prefix``
    (ℓ < 1)."""
    post = fit.point_posterior()
    names = np.array(sp.names)
    ok = (np.char.startswith(names[fit.grid.g], prefix)
          & np.char.startswith(names[fit.grid.h], prefix) & (fit.grid.ell < 1))  # fmt: skip
    return float(post[ok].sum())


def test_held_out_strain_is_placed_in_its_clade_at_its_depth(tmp_path: Path) -> None:
    rng = np.random.default_rng(5)
    refs = two_clade_species(rng)
    held = next(iter(two_clade_species(rng, 1).values()))  # a new clade-a genome
    panel = synthetic_panel(tmp_path, refs, dict.fromkeys(refs, "A"))
    for depth in (1.0, 3.0):
        prof = simulate([Strain(held, depth)], seed=int(depth))
        placed = place_species(prof, panel, species=[0])
        fit = placed.fits[0]
        assert fit.present_prob > 0.99
        assert placed_on(fit, placed.panels[0], "a") > 0.9
        mean, lo, hi = fit.depth_summary()
        assert mean == pytest.approx(depth, rel=0.1) and lo < depth < hi
        # accessory units of clade a the strain lacks, with no hits: doubted, not imputed
        sp = placed.panels[0]
        missing = [i for i, u in enumerate(sp.units) if 300 <= u < 400 and u not in held]
        assert fit.unit_carriage[missing].max() < 0.5


def test_absent_species_is_absent(tmp_path: Path) -> None:
    rng = np.random.default_rng(6)
    a = {f"a{i}": set(range(200)) | {u for u in range(200, 260) if rng.random() < 0.5}
         for i in range(5)}  # fmt: skip
    b = {f"b{i}": set(range(1000, 1200)) for i in range(5)}
    panel = synthetic_panel(tmp_path, a | b, dict.fromkeys(a, "A") | dict.fromkeys(b, "B"))
    prof = simulate([Strain(set(range(1000, 1200)), 2.0)])
    placed = place_species(prof, panel, species=[0, 1])
    assert placed.fits[0].present_prob < 0.01 and placed.fits[1].present_prob > 0.99
    assert screen(Evidence.from_profile(prof), panel)["species"].to_list() == [1]


def test_outside_organism_sharing_part_of_a_core_is_not_the_species(tmp_path: Path) -> None:
    """Phase 11's known failure: an organism outside the panel carrying 20-40% of an absent
    species' core units at about 3 hits each. No heuristic: the species' background β."""
    rng = np.random.default_rng(7)
    core = set(range(400))
    a = {f"a{i}": core | {u for u in range(400, 500) if rng.random() < 0.5} for i in range(6)}
    b = {f"b{i}": set(range(1000, 1300)) for i in range(4)}
    panel = synthetic_panel(tmp_path, a | b, dict.fromkeys(a, "A") | dict.fromkeys(b, "B"))
    for share in (0.2, 0.4):
        outside = set(rng.choice(400, int(400 * share), replace=False).tolist())
        outside |= set(range(2000, 2600))  # its own units, in no panel species
        prof = simulate([Strain(outside, 3 / (0.8 * 100)), Strain(set(range(1000, 1300)), 2.0)])
        placed = place_species(prof, panel, species=[0, 1])
        assert placed.fits[0].present_prob < 0.5, (share, placed.fits[0].log_bf)
        assert placed.fits[1].present_prob > 0.99


def test_strain_unlike_every_reference_is_novel(tmp_path: Path) -> None:
    rng = np.random.default_rng(8)
    core = set(range(200))
    pool = range(200, 400)
    genomes = {f"g{i}": core | {u for u in pool if rng.random() < 0.3} for i in range(40)}
    panel = synthetic_panel(tmp_path, genomes, dict.fromkeys(genomes, "A"), max_per_species=6)
    novel = core | {u for u in pool if rng.random() < 0.3}
    near = genomes[panel.genomes["name"][2]]  # a panel genome itself
    for strain, expect in ((novel, True), (near, False)):
        placed = place_species(simulate([Strain(strain, 2.0)]), panel, species=[0])
        fit = placed.fits[0]
        post = fit.point_posterior()
        assert (post[fit.grid.ell > NOVEL].sum() > 0.5) == expect


def test_species_sharing_units_are_fitted_together(tmp_path: Path) -> None:
    shared = set(range(100))  # MGnify90 clusters both species' genomes fall in
    a = {f"a{i}": shared | set(range(100, 400)) for i in range(4)}
    b = {f"b{i}": shared | set(range(1000, 1300)) for i in range(4)}
    panel = synthetic_panel(tmp_path, a | b, dict.fromkeys(a, "A") | dict.fromkeys(b, "B"))
    prof = simulate([Strain(a["a0"], 2.0), Strain(b["b0"], 1.0)])
    placed = place_species(prof, panel)
    assert placed.candidates["species"].to_list() == [0, 1] and placed.report["components"] == 1
    for s, depth in ((0, 2.0), (1, 1.0)):
        fit = placed.fits[s]
        assert fit.present_prob > 0.99
        assert fit.depth_summary()[0] == pytest.approx(depth, rel=0.1)
        sp = placed.panels[s]
        assert fit.unit_carriage[np.isin(sp.units, list(shared))].min() > 0.9
