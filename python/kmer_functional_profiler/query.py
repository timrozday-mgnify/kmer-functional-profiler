"""Query (phase 3 prototype): naive hit counts and containment per unit.

Reads are streamed through the Rust kernel at the index's ``t_max``; every sampled hash
is looked up in tier 1 and counts for each unit in its set whose own threshold it passes
(``hash <= max_hash_g``), which also discards most fingerprint false hits. No query
sketch is stored: only per-(unit, hash) counts are kept. Indexes imported from sourmash
signatures (``meta["hash"] == "sourmash"``) hash reads with sourmash instead; ``frames``
and ``genetic_code`` then do not apply.

Shared k-mers count for every unit that holds them in ``kmers_hit``; ``kmers_unique`` is
what a gather-style greedy assignment leaves each unit (the phase-4 detection baseline).
``coverage_em`` re-splits the hits of the units gather keeps by EM (phase-4 quantification);
``coverage_zi`` fits coverage to the k-mers present only (zero-inflated EM), the shipped
coverage: a diverged gene hits only part of a unit's k-mers, so ``coverage_em`` reads its depth
low (phase 7, step 35).
Two more baselines give each hit k-mer to one unit (:func:`assign_best`): winner-take-all
(``kmers_wta``) and uniqueness-first (``kmers_ufirst``), each with the coverage its k-mers'
hits give.
"""

import ctypes
import itertools
import json
import multiprocessing
import os
import resource
import sys
import threading
import time
from collections.abc import Callable, Iterable, Iterator, Sequence
from concurrent.futures import FIRST_COMPLETED, Future, ProcessPoolExecutor, wait
from contextlib import contextmanager
from dataclasses import fields
from functools import partial
from pathlib import Path
from typing import TYPE_CHECKING, Any, Final, NamedTuple

import numpy as np
import polars as pl
from scipy.sparse import coo_array
from scipy.sparse.csgraph import connected_components, reverse_cuthill_mckee
from scipy.special import betainc, gammaincc, gammaln
from scipy.stats import chi2, norm

from kmer_functional_profiler import _core
from kmer_functional_profiler.compat import sourmash_hits
from kmer_functional_profiler.index import (
    PIN_BITS,
    Index,
    IndexParams,
    PackedTable,
    base_rate,
    molecule,
    packed_layout,
    read_rows,
    span_nt,
    take,
)
from kmer_functional_profiler.survival import SurvivalModel

if TYPE_CHECKING:
    from kmer_functional_profiler.mask import Mask

DISTINCT_SAMPLE: Final = 256  # distinct sampled k-mers are counted on 1 in this of hash space


def peak_rss() -> int:
    """Peak resident set size of this process so far, in bytes."""
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return peak if sys.platform == "darwin" else peak * 1024  # Linux reports KiB


def rss() -> int:
    """Current resident set size in bytes (Linux); the peak elsewhere."""
    try:
        return int(Path("/proc/self/statm").read_text().split()[1]) * resource.getpagesize()
    except OSError:
        return peak_rss()


def _footprint() -> int | None:
    """macOS ``phys_footprint`` (what Activity Monitor calls Memory): anonymous and
    compressed memory, without clean pages of memory-mapped files."""
    buf = ctypes.create_string_buffer(512)  # rusage_info_v2; phys_footprint at byte 72
    if _libc is None or _libc.proc_pid_rusage(os.getpid(), 2, buf) != 0:
        return None
    return int.from_bytes(buf.raw[72:80], "little")


_libc = ctypes.CDLL(None) if sys.platform == "darwin" else None


def anon_rss() -> int | None:
    """Current anonymous RSS in bytes (Linux ``RssAnon``; macOS ``phys_footprint``): memory
    the process owns, without the page cache of memory-mapped files (the index), which the
    kernel can reclaim. None elsewhere."""
    if _libc is not None:
        return _footprint()
    try:
        with open("/proc/self/status") as f:
            return next(int(line.split()[1]) * 1024 for line in f if line.startswith("RssAnon:"))
    except (OSError, StopIteration):
        return None


ANON_SAMPLE_S: Final = 0.05  # anonymous RSS sampling interval while a stage runs
CGROUP_ROOT: Final = Path("/sys/fs/cgroup")
# cgroup memory.stat keys (v2, then v1) for the charge's parts: anonymous, file pages (page
# cache, including mapped index files), mapped file pages, shmem (tmpfs, e.g. a RAM $TMPDIR)
CGROUP_STAT: Final = {
    "anon": ("anon", "total_rss"),
    "file": ("file", "total_cache"),
    "mapped": ("file_mapped", "total_mapped_file"),
    "shmem": ("shmem", "total_shmem"),
}


def memory_cgroup(proc: Path = Path("/proc/self/cgroup"), root: Path = CGROUP_ROOT) -> Path | None:
    """The memory cgroup whose limit applies to this process (Linux, v2 or v1): the nearest
    of its own cgroup and their parents with a limit set (Slurm's job or step cgroup), else
    its own. None without cgroup memory accounting (macOS, containers without it)."""
    try:
        lines = proc.read_text().splitlines()
    except OSError:
        return None
    for line in lines:
        _, controllers, rel = line.split(":", 2)
        v2 = controllers == ""
        if not v2 and "memory" not in controllers.split(","):
            continue
        base = root if v2 else root / "memory"
        own = base / rel.lstrip("/")
        limit, usage = (
            ("memory.max", "memory.current")
            if v2
            else ("memory.limit_in_bytes", "memory.usage_in_bytes")
        )
        if not (own / usage).exists():
            continue
        for d in (own, *own.parents):
            if not d.is_relative_to(base):
                break
            try:
                text = (d / limit).read_text().strip()
            except OSError:
                continue
            if text != "max" and int(text) < 2**60:  # v1 writes "no limit" as ~2**63
                return d
        return own
    return None


def cgroup_memory(cgroup: Path) -> dict[str, int]:
    """The cgroup's memory charge now (what Slurm's limit is checked against), its kernel
    high-water mark where the kernel keeps one, its limit, its OOM kills and
    :data:`CGROUP_STAT`'s parts, in bytes; empty if unreadable."""
    v2 = (cgroup / "memory.current").exists()
    files = (
        {"cgroup": "memory.current", "cgroup_peak": "memory.peak", "cgroup_limit": "memory.max"}
        if v2
        else {
            "cgroup": "memory.usage_in_bytes",
            "cgroup_peak": "memory.max_usage_in_bytes",
            "cgroup_limit": "memory.limit_in_bytes",
        }
    )
    out = {}
    for key, name in files.items():
        try:
            text = (cgroup / name).read_text().split()[0]
        except (OSError, IndexError):
            continue
        if text != "max" and int(text) < 2**60:  # v1 writes "no limit" as ~2**63
            out[key] = int(text)
    # processes the kernel's OOM killer killed in the cgroup (v2 memory.events; v1
    # memory.oom_control from Linux 4.13)
    try:
        events = (cgroup / ("memory.events" if v2 else "memory.oom_control")).read_text()
        kills = dict(line.split() for line in events.splitlines())
        if "oom_kill" in kills:
            out["cgroup_oom_kill"] = int(kills["oom_kill"])
    except (OSError, ValueError):
        pass
    try:
        stat = dict(line.split() for line in (cgroup / "memory.stat").read_text().splitlines())
    except (OSError, ValueError):
        return out
    for part, keys in CGROUP_STAT.items():
        key = keys[0] if v2 else keys[1]
        if key in stat:
            out[f"cgroup_{part}"] = int(stat[key])
    return out


class Timer:
    """Wall time, CPU time and peak RSS per query stage, and counts (the query cost study).

    A stage entered once per read batch accumulates its times. ``peak_rss`` is the
    process's peak at the stage's last exit, so the stage that raises it shows as a step.
    It counts resident pages of memory-mapped files too, so where the OS reports it (Linux,
    macOS) ``peak_anon`` is the same high-water mark for anonymous memory only, sampled every
    ``ANON_SAMPLE_S`` while a stage runs (the kernel keeps no anonymous peak).

    Under a memory cgroup (Linux; :func:`memory_cgroup`, e.g. a Slurm job), each stage also
    keeps ``peak_cgroup``, the sampled high-water mark of the cgroup's charge, which is what
    Slurm's memory limit kills on: anonymous memory plus page cache, mapped index pages
    included, plus shmem (a ``$TMPDIR`` in RAM). ``peak_cgroup_{anon,file,mapped,shmem}``
    are its parts' own sampled peaks, ``cgroup_peak`` the kernel's high-water mark at the
    stage's end where it keeps one (v1, v2 from Linux 5.19; it counts the whole cgroup's
    life, other processes of the job included) and ``cgroup_limit`` the limit.

    With ``path``, the stats are rewritten there after every stage, with ``running`` set to
    the stage in progress, so a process killed mid-query (out of memory) leaves the stages
    it finished and the one it died in. With ``log``, each stage's start and end are printed
    to stderr with the current and peak RSS.
    """

    def __init__(self, path: Path | None = None, *, log: bool = False) -> None:
        self.stages: dict[str, dict[str, float]] = {}
        self.counts: dict[str, Any] = {}  # numbers, and the hits_per_pair histogram
        self.path, self.log, self.running = path, log, ""
        self.start = time.perf_counter()
        self.peak_anon = anon_rss()
        self.cgroup = memory_cgroup()
        self.peak_cgroup: dict[str, int] = {}
        self._update_cgroup()

    def _sample_anon(self, stop: threading.Event) -> None:
        while not stop.wait(ANON_SAMPLE_S):
            self._update_anon()
            self._update_cgroup()

    def _update_anon(self) -> None:
        now = anon_rss()
        if now is not None and self.peak_anon is not None:
            self.peak_anon = max(self.peak_anon, now)

    def _update_cgroup(self) -> None:
        if self.cgroup is None:
            return
        for key, value in cgroup_memory(self.cgroup).items():
            if key in ("cgroup_peak", "cgroup_limit", "cgroup_oom_kill"):  # read at the end
                continue
            name = "peak_" + key
            self.peak_cgroup[name] = max(self.peak_cgroup.get(name, 0), value)

    def _log(self, event: str, stage: str) -> None:
        if self.log:
            gib = 2**30
            print(
                f"[{time.perf_counter() - self.start:9.1f}s] {event:5} {stage:<11}"
                f" rss {rss() / gib:7.2f} GiB  peak {peak_rss() / gib:7.2f} GiB"
                + (
                    "" if self.peak_anon is None else f"  anon peak {self.peak_anon / gib:7.2f} GiB"
                ),
                file=sys.stderr,
                flush=True,
            )

    def write(self) -> None:
        if self.path is not None:
            self.path.write_text(json.dumps(self.as_dict(), indent=2))

    @contextmanager
    def __call__(self, stage: str) -> Iterator[None]:
        outer, self.running = self.running, stage
        self._log("start", stage)
        self.write()
        stop = threading.Event()
        sampled = self.peak_anon is not None or self.cgroup is not None
        if sampled and not outer:  # one sampler, for the outermost stage
            threading.Thread(target=self._sample_anon, args=(stop,), daemon=True).start()
        wall, cpu = time.perf_counter(), time.process_time()
        try:
            yield
        finally:
            stop.set()
            self._update_anon()
            self._update_cgroup()
            s = self.stages.setdefault(stage, {"wall_s": 0.0, "cpu_s": 0.0})
            s["wall_s"] += time.perf_counter() - wall
            s["cpu_s"] += time.process_time() - cpu
            s["peak_rss"] = peak_rss()
            if self.peak_anon is not None:
                s["peak_anon"] = self.peak_anon
            s |= self.peak_cgroup
            if self.cgroup is not None:
                now = cgroup_memory(self.cgroup)
                s |= {k: now[k] for k in ("cgroup_peak", "cgroup_limit") if k in now}
            self.running = outer
            self._log("end", stage)
            self.write()

    def as_dict(self) -> dict[str, object]:
        return {"running": self.running, "stages": self.stages, "counts": self.counts}


def _ids(x: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Distinct values of ``x`` and each element's index into them; the index is int32 when
    it fits (half the memory of ``np.unique``'s int64), which changes no result."""
    distinct, inverse = np.unique(x, return_inverse=True)
    return distinct, inverse.astype(np.int32 if len(x) < 2**30 else np.int64, copy=False)


def components(kmers: pl.DataFrame) -> tuple[int, int, int]:
    """Connected components of hit units linked by shared k-mers.

    Returns the number of components, the largest's units and its (unit, hash) pairs.
    """
    units, col = _ids(kmers["unit"].to_numpy())
    _, row = _ids(kmers["hash"].to_numpy())
    n = len(units)
    if n == 0:
        return 0, 0, 0
    label = _unit_components(col, row, n)
    per_unit = np.bincount(label)
    largest = int(per_unit.argmax())
    n_units = int(per_unit[largest])
    return len(per_unit), n_units, int((label[col] == largest).sum())


def component_labels(kmers: pl.DataFrame) -> pl.DataFrame:
    """``unit`` and its ``component`` label (:func:`_unit_components`) for every unit in
    ``kmers`` (``unit``, ``hash``)."""
    units, col = _ids(kmers["unit"].to_numpy())
    _, row = _ids(kmers["hash"].to_numpy())
    label = _unit_components(col, row, len(units)) if len(units) else np.empty(0, np.int64)
    return pl.DataFrame(
        {"unit": units, "component": label}, schema_overrides={"component": pl.UInt32}
    ).with_columns(pl.col("unit").cast(kmers.schema["unit"]))


def _unit_components(col: np.ndarray, row: np.ndarray, n: int) -> np.ndarray:
    """Component label of each of ``n`` units, linked by pairs (unit ``col``, k-mer ``row``);
    components are numbered in order of their smallest unit."""
    # Units only, each linked to the next holder of the same k-mer: a unit graph with fewer
    # edges than pairs, instead of the bipartite unit-k-mer graph.
    order = np.argsort(row, kind="stable")
    by_kmer = col[order]
    same = row[order][1:] == row[order][:-1]
    graph = coo_array(
        (np.ones(int(same.sum()), dtype=np.int8), (by_kmer[:-1][same], by_kmer[1:][same])),
        shape=(n, n),
    )
    label: np.ndarray = connected_components(graph, directed=False)[1]
    return label


LINK_TAUS: Final = (0.005, 0.01, 0.02, 0.05, 0.1, 0.2)  # A1's sweep of link-strength cuts
# a stricter promiscuity cut (build's max_groups, 16) for the largest component (step 40)
HOLDER_CUTS: Final = (2, 4, 8, 16, 32)


def link_cuts(kmers: pl.DataFrame) -> dict[str, int]:
    """What weak-link cuts would do to the largest component (phase 6, step 33, A1).

    ``kmers`` has ``unit``, ``hash`` and ``hits``. In the largest component of units linked
    by shared k-mers, a link's strength is its shared k-mers over the smaller unit's k-mers.
    Returns its size (``detected_largest_component_units``, ``_pairs``: the largest
    component the fits see, unlike :func:`components` over all hit units), its ``links``,
    the weak ones (``links_one_kmer``: one shared k-mer;
    ``links_le2_hits``: at most 2 hits on the shared k-mers), and per cut tau the links below
    it (``cut_links_{tau}``) and the largest component left (``cut_largest_units_{tau}``).

    With ``holders`` (units holding the k-mer in the index), per h in ``HOLDER_CUTS`` what
    dropping k-mers held by more than h units would leave, as a build's ``max_groups`` h
    would (bar the floored units' refill from other candidates, which the build does):
    ``holders_cut_pairs_{h}`` (pairs of the component dropped) and the largest component
    left, ``holders_cut_largest_units_{h}`` and ``_pairs_{h}``.
    """
    units, col = _ids(kmers["unit"].to_numpy())
    _, row = _ids(kmers["hash"].to_numpy())
    if len(units) == 0:
        return {}
    label = _unit_components(col, row, len(units))
    keep = label[col] == np.bincount(label).argmax()
    own = pl.DataFrame({"u": col[keep], "k": row[keep], "hits": kmers["hits"].to_numpy()[keep]})
    n_kmers = np.bincount(col[keep], minlength=len(units))
    links = (
        own.join(own.select("u", "k"), on="k", suffix="2")
        .filter(pl.col("u") < pl.col("u2"))
        .group_by("u", "u2")
        .agg(shared=pl.len(), shared_hits=pl.col("hits").sum())
    )
    a, b = links["u"].to_numpy(), links["u2"].to_numpy()
    strength = links["shared"].to_numpy() / np.minimum(n_kmers[a], n_kmers[b])
    out = {
        "detected_largest_component_units": int(n_kmers.astype(bool).sum()),
        "detected_largest_component_pairs": int(keep.sum()),
        "links": links.height,
        "links_one_kmer": int((links["shared"] == 1).sum()),
        "links_le2_hits": int((links["shared_hits"] <= 2).sum()),
    }
    for tau in LINK_TAUS:
        kept = strength >= tau
        graph = coo_array(
            (np.ones(int(kept.sum()), dtype=np.int8), (a[kept], b[kept])),
            shape=(len(units), len(units)),
        )
        left = connected_components(graph, directed=False)[1][np.unique(col[keep])]
        out[f"cut_links_{tau}"] = int((~kept).sum())
        out[f"cut_largest_units_{tau}"] = int(np.bincount(left).max())
    if "holders" in kmers.columns:
        holders = kmers["holders"].to_numpy()[keep]
        c, r = col[keep], row[keep]
        for h in HOLDER_CUTS:
            kept = holders <= h
            label = _unit_components(c[kept], r[kept], len(units))[c[kept]]
            pairs = np.bincount(label) if kept.any() else np.zeros(1, dtype=np.intp)
            out[f"holders_cut_pairs_{h}"] = int((~kept).sum())
            out[f"holders_cut_largest_pairs_{h}"] = int(pairs.max())
            out[f"holders_cut_largest_units_{h}"] = len(np.unique(c[kept][label == pairs.argmax()]))
    return out


MAX_BATCH_PAIRS: Final = 2_000_000  # (unit, hash) rows per batch of components fitted at once


def component_batches(kmers: pl.DataFrame, max_pairs: int | None = None) -> list[pl.DataFrame]:
    """``kmers`` (one row per ``unit``, ``hash``) split into batches of whole components of
    units linked by shared k-mers, about ``max_pairs`` (default ``MAX_BATCH_PAIRS``) rows
    each (a larger component is a batch of its own). Fits that stop per component
    (:func:`em`, :func:`em_pin`) give the same result per batch as on the whole, with the
    working memory of one batch."""
    if kmers.height == 0:
        return [kmers]
    units, col = _ids(kmers["unit"].to_numpy())
    _, row = _ids(kmers["hash"].to_numpy())
    label = _unit_components(col, row, len(units))
    size = np.bincount(label[col])  # rows per component
    batch = ((np.cumsum(size) - size) // (max_pairs or MAX_BATCH_PAIRS))[
        label[col]
    ]  # by each component's start
    del units, row, label
    order = np.argsort(batch, kind="stable")
    rows = kmers[order]
    ends = np.cumsum(np.unique(batch, return_counts=True)[1])
    return [rows.slice(lo, hi - lo) for lo, hi in zip(np.r_[0, ends[:-1]], ends, strict=True)]


def _ranges(starts: np.ndarray, lengths: np.ndarray) -> np.ndarray:
    """Concatenated ``arange(start, start + length)`` for each pair."""
    ends = np.cumsum(lengths)
    total = int(ends[-1]) if len(ends) else 0
    out: np.ndarray = np.arange(total) + np.repeat(starts - ends + lengths, lengths)
    return out


class _Block(NamedTuple):
    """The components still being fitted: their units (grouped by component), the pairs of
    those units, each pair's local unit ``c`` and local k-mer ``r``, the k-mers, and each
    k-mer's weight from holders fitted elsewhere (``offset``: 0, or per k-mer in a block of a
    component fitted block-wise, added to the k-mer's expected hits)."""

    unit: np.ndarray
    pairs: np.ndarray
    c: np.ndarray
    r: np.ndarray
    kmer: np.ndarray
    offset: np.ndarray | float = 0.0


# (coverage, second parameter) -> (new coverage, new second parameter, expected hits)
Step = Callable[[np.ndarray, np.ndarray], tuple[np.ndarray, np.ndarray, np.ndarray]]
# pairs -> each pair's weight in its k-mer's expected hits, at the current lam and other
PairWeight = Callable[[np.ndarray], np.ndarray]

# a component with more pairs is fitted block-wise (query --max-fit-pairs); a whole fit takes
# ~65 bytes per pair, block-wise ~90 in two blocks (step 41, synthetic): it saves memory only
# with many blocks, and on real components whole fits ran 25-32x faster at lower peak (step 42)
MAX_FIT_PAIRS = 100_000_000
# the zero-inflated EM's step cap (query --zi-max-iter), the plain EM's: 3,000 and 10,000 left
# the same units unconverged at 3-10x the fit's time, with no change in accuracy (step 45)
ZI_MAX_ITER = 1000
BLOCK_STEPS: Final = 30  # EM steps per block per round of a block-wise fit
MAX_BLOCK_ROUNDS: Final = 1000  # rounds over a block-wise component's blocks


def _squarem(
    step: Step,
    la: np.ndarray,
    se: np.ndarray,
    starts: np.ndarray,
    seg: np.ndarray,
    step_max: np.ndarray,
    other_max: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """One SQUAREM cycle per component (Varadhan & Roland 2008, SqS3): two EM steps, a
    squared extrapolation along them, and an EM step from there. Returns the point the
    last step started from, that step's result and expected hits, as one EM step would.

    The step length is capped by ``step_max`` per component, which grows 4x whenever the cap
    is reached (updated in place). Extrapolated coverages are kept above 1e-3 x the second
    step's (EM is multiplicative: a coverage of 0 could never grow back) and ``other`` inside
    [0, ``other_max``]. A component whose EM residual grows over the cycle falls back to its
    two plain steps and its cap to 1, so a cycle never does worse than plain EM's progress.
    """
    l1, s1, _ = step(la, se)
    l2, s2, att2 = step(l1, s1)

    def per_component(x: np.ndarray) -> np.ndarray:
        return np.add.reduceat(x, starts)

    r_l, r_s = l1 - la, s1 - se
    v_l, v_s = l2 - l1 - r_l, s2 - s1 - r_s
    rr = per_component(r_l**2 + r_s**2)
    vv = per_component(v_l**2 + v_s**2)
    alpha = -np.sqrt(_div(rr, vv))
    alpha = np.clip(np.where(vv > 0, alpha, -1.0), -step_max, -1.0)  # -1: the plain two steps
    step_max[alpha <= -step_max] *= 4
    a = np.repeat(alpha, seg)
    lx = np.maximum(la - 2 * a * r_l + a * a * v_l, 1e-3 * l2)
    sx = np.clip(se - 2 * a * r_s + a * a * v_s, 0.0, other_max)
    l3, s3, att3 = step(lx, sx)
    worse = per_component((l3 - lx) ** 2 + (s3 - sx) ** 2) > rr
    step_max[worse] = 1.0
    back = np.repeat(worse, seg)
    return (
        np.where(back, l1, lx),
        np.where(back, s1, sx),
        np.where(back, l2, l3),
        np.where(back, s2, s3),
        np.where(back, att2, att3),
    )


def _fit_components(
    col: np.ndarray,
    row: np.ndarray,
    lam: np.ndarray,
    other: np.ndarray,
    make_step: Callable[[_Block], Step],
    tol: float,
    max_iter: int,
    report: dict[str, float] | None = None,
    other_max: float = 1.0,
    pair_weight: PairWeight | None = None,
    offset: np.ndarray | None = None,
    pair_ids: np.ndarray | None = None,
) -> np.ndarray:
    """Iterate an EM per component of units linked by shared k-mers (pairs of unit ``col``
    and k-mer ``row``); return each unit's expected hits at its last step.

    Components share no k-mer, so each is its own fit and stops at its own convergence:
    its coverages ``lam`` change by at most ``tol`` x its largest, and ``other`` (presence
    or scale) by at most ``tol``. An iteration only touches components still moving (most
    converge in a few dozen, a few run to ``max_iter``); they are compacted into a new
    ``_Block``, and ``make_step`` called on it, once half of the current block has
    converged. ``lam`` and ``other`` are updated in place.

    Steps are accelerated by SQUAREM (:func:`_squarem`; ``other`` kept in [0, ``other_max``]),
    which reaches the same fixed point in fewer steps; ``max_iter`` counts EM steps, three
    per cycle. Convergence is still tested on a plain EM step, so it means the same.

    With ``pair_weight``, a component of more than ``MAX_FIT_PAIRS`` pairs is fitted
    block-wise (:func:`_fit_blockwise`). Only units with pairs are fitted; ``offset`` (per
    k-mer ``row``) and ``pair_ids`` (each pair's index in the caller's pairs) serve the blocks.
    ``report`` gains ``em_iterations``, ``em_unconverged_units`` and ``em_max_change``: how
    far the furthest component stopped at ``max_iter`` was from the test, as its last
    change (the larger of the coverage change over its largest coverage and the ``other``
    change; ``tol`` passes; 0 if none stopped). Components that stop there are mostly
    units sharing all their hit k-mers, whose coverage drifts slowly between them along a
    nearly flat likelihood (phase 7, step 29).
    """
    n = len(lam)
    label = _unit_components(col, row, n)
    if pair_weight is not None:
        size = np.bincount(label[col])  # pairs per component
        if (size > MAX_FIT_PAIRS).any():
            in_big = (size > MAX_FIT_PAIRS)[label[col]]
            small = np.flatnonzero(~in_big)
            attributed = _fit_components(
                col[small],
                row[small],
                lam,
                other,
                make_step,
                tol,
                max_iter,
                report,
                other_max,
                pair_ids=small,
            )
            for c in np.flatnonzero(size > MAX_FIT_PAIRS):
                pairs = np.flatnonzero(label[col] == c)
                args = (col, row, lam, other, make_step, tol, report, other_max)
                attributed += _fit_blockwise(pairs, pair_weight, *args)
            return attributed
    attributed = np.zeros(n)
    active = np.argsort(label, kind="stable")  # units grouped by component
    active = active[np.bincount(col, minlength=n)[active] > 0]
    seg = np.diff(np.flatnonzero(np.r_[True, np.diff(label[active]) != 0, True]))
    by_unit = np.argsort(col, kind="stable").astype(col.dtype)
    pair_start = np.r_[0, np.cumsum(np.bincount(col, minlength=n))]
    local = np.empty(n, dtype=col.dtype)
    step_max = np.ones(len(seg))  # per component, kept through compactions
    change = np.zeros(len(seg))  # per component, its last change against the test
    it = 0
    while len(active) and it < max_iter:
        pairs = by_unit[_ranges(pair_start[active], np.diff(pair_start)[active])]
        local[active] = np.arange(len(active))
        kmer, r = _ids(row[pairs])
        block = _Block(
            active,
            pairs if pair_ids is None else pair_ids[pairs],
            local[col[pairs]],
            r,
            kmer,
            0.0 if offset is None else offset[kmer],
        )
        step = make_step(block)
        starts = np.r_[0, np.cumsum(seg)[:-1]]
        la, se, done = lam[active], other[active], np.zeros(len(active), dtype=bool)
        while it < max_iter and done.sum() * 2 <= len(active):
            if it + 3 <= max_iter:
                la, se, new, new_se, att = _squarem(step, la, se, starts, seg, step_max, other_max)
                it += 3
            else:
                it += 1
                new, new_se, att = step(la, se)
            change = np.maximum(
                _div(
                    np.maximum.reduceat(np.abs(new - la), starts),
                    np.maximum.reduceat(new, starts),
                ),
                np.maximum.reduceat(np.abs(new_se - se), starts),
            )
            converged = change <= tol
            la, se = new, new_se
            now = np.repeat(converged, seg) & ~done
            lam[active[now]], other[active[now]], attributed[active[now]] = (
                la[now],
                se[now],
                att[now],
            )
            done |= now
        rest = ~done  # at max_iter, units still moving keep their last values
        lam[active[rest]], other[active[rest]], attributed[active[rest]] = (
            la[rest],
            se[rest],
            att[rest],
        )
        moving = ~done[starts]
        active, seg, step_max = active[np.repeat(moving, seg)], seg[moving], step_max[moving]
        change = change[moving]
    if report is not None:  # units still moving here stopped at max_iter, not converged
        report["em_iterations"] = max(report.get("em_iterations", 0), it)
        report["em_unconverged_units"] = report.get("em_unconverged_units", 0) + len(active)
        stopped = float(change.max()) if len(active) else 0.0
        report["em_max_change"] = max(report.get("em_max_change", 0.0), stopped)
    return attributed


def _unit_order(col: np.ndarray, row: np.ndarray, units: np.ndarray) -> np.ndarray:
    """``units`` (sorted) in reverse Cuthill-McKee order of their graph (linked by shared
    k-mers ``row``), so units that share k-mers sit near each other in the order."""
    local = np.searchsorted(units, col)
    order = np.argsort(row, kind="stable")
    by_kmer = local[order]
    same = row[order][1:] == row[order][:-1]  # each unit linked to the k-mer's next holder
    n = len(units)
    graph = coo_array(
        (np.ones(int(same.sum()), dtype=np.int32), (by_kmer[:-1][same], by_kmer[1:][same])),
        shape=(n, n),
    ).tocsr()
    ordered: np.ndarray = units[reverse_cuthill_mckee(graph, symmetric_mode=False)]
    return ordered


def _fit_blockwise(
    pairs: np.ndarray,
    pair_weight: PairWeight,
    col: np.ndarray,
    row: np.ndarray,
    lam: np.ndarray,
    other: np.ndarray,
    make_step: Callable[[_Block], Step],
    tol: float,
    report: dict[str, float] | None,
    other_max: float,
) -> np.ndarray:
    """Fit one component (its ``pairs``) in blocks of units of at most ``MAX_FIT_PAIRS``
    pairs, so a step's working memory is one block's. A round takes ``BLOCK_STEPS`` EM steps
    on each block in turn, with the weight of its k-mers' other holders held fixed (an
    ``offset`` per k-mer). Few steps per block, not a fit to convergence: a block fitted
    alone can drive a unit towards 0 that the others' next round would need back, and EM is
    slow to bring it back. At convergence that is the whole component's fixed point.

    Two things keep the rounds few (phase 7, step 29). Blocks are runs of units in reverse
    Cuthill-McKee order (:func:`_unit_order`), so units sharing k-mers mostly fall in one
    block and fewer k-mers are held fixed across blocks. And rounds are accelerated as
    :func:`_squarem` accelerates steps: two rounds, a squared extrapolation along them
    (coverages kept above 1e-3 x the second round's, ``other`` in [0, ``other_max``]), and a
    round from there, falling back to the two plain rounds if that round moves more than
    the first did. Without both, block-wise EM crept along slow directions and its round
    test could pass ~0.1 short of the fixed point. It stops when the last round changes no
    unit's parameters by more than :func:`_fit_components`' test, or at ``MAX_BLOCK_ROUNDS``.

    ``report`` gains ``em_block_rounds`` (the most rounds a component took) and counts
    ``em_iterations`` as rounds x ``BLOCK_STEPS``; a component still moving after the last
    round counts as unconverged, with its last change in ``em_max_change``.
    """
    units, per_unit = np.unique(col[pairs], return_counts=True)
    ordered = _unit_order(col[pairs], row[pairs], units)
    count = per_unit[np.searchsorted(units, ordered)]
    chunk = np.zeros(len(lam), dtype=np.int64)
    chunk[ordered] = (np.cumsum(count) - count) // MAX_FIT_PAIRS
    pair_chunk = chunk[col[pairs]]
    order = np.argsort(pair_chunk, kind="stable")
    blocks = np.split(pairs[order], np.flatnonzero(np.diff(pair_chunk[order])) + 1)
    n_rows = int(row[pairs].max()) + 1
    attributed = np.zeros(len(lam))

    def one_round() -> tuple[np.ndarray, np.ndarray]:
        total = np.bincount(row[pairs], weights=pair_weight(pairs), minlength=n_rows)
        for b in blocks:
            total -= np.bincount(row[b], weights=pair_weight(b), minlength=n_rows)
            att = _fit_components(
                col[b],
                row[b],
                lam,
                other,
                make_step,
                tol,
                BLOCK_STEPS,
                None,
                other_max,
                offset=total,
                pair_ids=b,
            )
            own = np.unique(col[b])
            attributed[own] = att[own]
            total += np.bincount(row[b], weights=pair_weight(b), minlength=n_rows)
        return lam[units].copy(), other[units].copy()

    def moved(la: np.ndarray, se: np.ndarray, la0: np.ndarray, se0: np.ndarray) -> float:
        return max(
            float(np.abs(la - la0).max() / max(la.max(), 1e-300)), float(np.abs(se - se0).max())
        )

    step_max, change, rounds = 1.0, np.inf, 0
    while change > tol and rounds + 3 <= MAX_BLOCK_ROUNDS:
        l0, s0 = lam[units].copy(), other[units].copy()
        l1, s1 = one_round()
        l2, s2 = one_round()
        att2 = attributed[units].copy()
        r_l, r_s = l1 - l0, s1 - s0
        v_l, v_s = l2 - l1 - r_l, s2 - s1 - r_s
        rr = float((r_l**2).sum() + (r_s**2).sum())
        vv = float((v_l**2).sum() + (v_s**2).sum())
        alpha = min(max(-np.sqrt(rr / vv) if vv > 0 else -1.0, -step_max), -1.0)
        if alpha <= -step_max:
            step_max *= 4
        lx = np.maximum(l0 - 2 * alpha * r_l + alpha**2 * v_l, 1e-3 * l2)
        sx = np.clip(s0 - 2 * alpha * r_s + alpha**2 * v_s, 0.0, other_max)
        lam[units], other[units] = lx, sx
        l3, s3 = one_round()
        rounds += 3
        if float(((l3 - lx) ** 2).sum() + ((s3 - sx) ** 2).sum()) > rr:
            step_max = 1.0  # back to the two plain rounds
            lam[units], other[units], attributed[units] = l2, s2, att2
            change = moved(l2, s2, l1, s1)
        else:
            change = moved(l3, s3, lx, sx)
    if report is not None:
        report["em_iterations"] = max(report.get("em_iterations", 0), rounds * BLOCK_STEPS)
        report["em_block_rounds"] = max(report.get("em_block_rounds", 0), rounds)
        stopped = change > tol
        report["em_unconverged_units"] = report.get("em_unconverged_units", 0) + (
            len(units) if stopped else 0
        )
        report["em_max_change"] = max(report.get("em_max_change", 0.0), change if stopped else 0.0)
    return attributed


def unit_hits(
    table: PackedTable, max_hash_g: np.ndarray, hashes: np.ndarray, reads: np.ndarray
) -> pl.DataFrame:
    """Expand sampled query hashes to (unit, hash, read, pin_q, holders) hits through ``table``.

    A hit counts for a unit only if the hash passes that unit's ``max_hash_g``; ``pin_q`` is
    the k-mer's quantised ``p_in`` in that unit and ``holders`` the number of index units
    it counts for. Rows follow the input order, then each set's order. The table's arrays
    are read in place (``_core.unit_hits``): one bucket and one set per hash, no keys.
    """
    u64 = np.uint64
    return pl.DataFrame(
        _core.unit_hits(
            table,
            np.asarray(max_hash_g, dtype=u64),
            np.asarray(hashes, dtype=u64),
            np.asarray(reads, dtype=u64),
        )
    )


# Tier-1 buckets per read of its key arrays when walking it in hash order (~16 M keys).
WALK_BUCKETS: Final = 2**22
WALK_HASHES: Final = 2**22  # hit hashes expanded to their sets at a time


def _key_sets(table: PackedTable, hashes: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """The sorted ``hashes`` that ``table`` holds and their set ids, as
    :meth:`PackedTable.lookup` finds them, with its key arrays read ``WALK_BUCKETS``
    buckets at a time in hash order (``read_rows``)."""
    hashes = hashes[hashes <= np.uint64(table.max_hash)]
    below = table.shift + table.fp_bits  # hash bits below the bucket
    block = (hashes >> np.uint64(below)) // np.uint64(WALK_BUCKETS)
    cuts = [0, *(np.flatnonzero(np.diff(block)) + 1), len(hashes)]
    none = np.zeros(1, np.uint8)
    found, sets = [np.empty(0, np.uint64)], [np.empty(0, np.int64)]
    for a, b in itertools.pairwise(cuts):
        if a == b:
            continue
        b0 = int(block[a]) * WALK_BUCKETS
        b1 = min(b0 + WALK_BUCKETS, len(table.offsets) - 1)
        offsets = read_rows(table.offsets, b0, b1 + 1, 2 * (b - a)).astype(np.uint64)
        k0, k1 = int(offsets[0]), int(offsets[-1])
        # The block as a table of its own: its hashes shifted down to bucket 0.
        part = PackedTable(
            table.max_hash,
            table.lead_bits,
            table.bucket_bits,
            table.fp_bits,
            offsets - offsets[0],
            read_rows(table.fingerprints, k0, k1, b - a),
            read_rows(table.set_ids, k0, k1, b - a),
            none,
            none,
        )
        ids = part.lookup(hashes[a:b] - np.uint64(b0 << below))
        hit = ids >= 0
        found.append(hashes[a:b][hit])
        sets.append(ids[hit])
    return np.concatenate(found), np.concatenate(sets)


class _TierHits(NamedTuple):
    """The part of a tier a sample's hashes hit, in memory: exact keys (no fingerprint
    bits dropped), the units passing their ``max_hash_g`` only, numbered ``0..`` in
    ``units`` (index unit ids) order, with their ``max_hash_g``."""

    table: PackedTable
    max_hash_g: np.ndarray
    units: np.ndarray

    def hits(self, hashes: np.ndarray, reads: np.ndarray) -> pl.DataFrame:
        hits = unit_hits(self.table, self.max_hash_g, hashes, reads)
        return hits.with_columns(unit=pl.Series(self.units[hits["unit"].to_numpy()]))


def tier_hits(
    table: PackedTable,
    max_hash_g: np.ndarray,
    hashes: np.ndarray | Iterable[np.ndarray],
) -> _TierHits:
    """What :func:`unit_hits` would give the sorted, distinct sampled ``hashes`` (or their
    consecutive hash ranges, each sorted), with ``table`` and ``max_hash_g`` read once in
    order (plan, phase 6, step 40 (c)): keys by bucket range (:func:`_key_sets`), then
    :func:`_resolve`."""
    ranges = [hashes] if isinstance(hashes, np.ndarray) else hashes
    found = [_key_sets(table, part) for part in ranges]
    return _resolve(
        table,
        max_hash_g,
        np.concatenate([f for f, _ in found] or [np.empty(0, np.uint64)]),
        np.concatenate([s for _, s in found] or [np.empty(0, np.int64)]),
    )


def _resolve(
    table: PackedTable, max_hash_g: np.ndarray, found: np.ndarray, set_ids: np.ndarray
) -> _TierHits:
    """The rows of the sorted hit hashes ``found`` (their ``set_ids`` in ``table``): the sets
    hit read in set order and the units' ``max_hash_g`` in unit order (``take``), so memory
    is a block of the tier plus the result, at any index size; ``unit_hits`` on the result
    gives the rows of ``table`` for any of the hashes walked, in the same order."""
    sets, of_set = np.unique(set_ids, return_inverse=True)
    del set_ids
    ends = take(table.set_offsets, np.stack([sets, sets + 1], axis=1).ravel()).astype(np.int64)
    del sets
    starts, lengths = ends[0::2], ends[1::2] - ends[0::2]
    values = take(table.set_values, _ranges(starts, lengths)).astype(np.uint64)
    del ends, starts
    units, local = np.unique(values >> np.uint64(PIN_BITS), return_inverse=True)
    limit = take(max_hash_g, units).astype(np.uint64)
    # unit << PIN_BITS | pin_q, in the smallest width that holds it
    value_type = np.uint32 if len(units) < 2 ** (32 - PIN_BITS) else np.uint64
    local_values = (local.astype(value_type) << value_type(PIN_BITS)) | (
        values & np.uint64(2**PIN_BITS - 1)
    ).astype(value_type)
    del values
    first = np.cumsum(lengths) - lengths
    # Per chunk of hashes, each's set in set order, kept where it passes the unit's limit:
    # the hashes with a row kept, their row counts and the rows' values.
    keys, counts, kept = [np.empty(0, np.uint64)], [np.empty(0, np.int64)], []
    for a in range(0, len(found), WALK_HASHES):
        s = of_set[a : a + WALK_HASHES]
        at = _ranges(first[s], lengths[s])
        key = np.repeat(np.arange(a, a + len(s)), lengths[s])
        keep = found[key] <= limit[local[at]]
        n = np.bincount(key[keep] - a, minlength=len(s))
        keys.append(found[a : a + len(s)][n > 0])
        counts.append(n[n > 0])
        kept.append(local_values[at[keep]])
    del local, local_values, of_set, first, lengths
    distinct, n = np.concatenate(keys), np.concatenate(counts)
    del keys, counts
    value = np.concatenate(kept or [np.empty(0, value_type)])
    del kept
    layout = packed_layout(table.max_hash, 0, len(distinct))
    fp_bits = 64 - layout["lead_bits"] - layout["bucket_bits"]  # the whole key: no false hits
    buckets = np.bincount(
        (distinct >> np.uint64(fp_bits)).astype(np.int64), minlength=2 ** layout["bucket_bits"]
    )
    exact = PackedTable(
        table.max_hash,
        layout["lead_bits"],
        layout["bucket_bits"],
        fp_bits,
        np.concatenate([[0], np.cumsum(buckets)]).astype(np.uint64),
        distinct & np.uint64(2**fp_bits - 1),
        np.arange(len(distinct), dtype=np.uint32 if len(distinct) < 2**32 else np.uint64),
        np.concatenate([[0], np.cumsum(n)]).astype(np.uint64),
        value,
    )
    return _TierHits(exact, limit, units.astype(np.uint32))


class _Joint:
    """Several indexes with one hash scheme queried as one: unit ids are offset by the units
    of the indexes before, so gather, EM and the posterior see one hit table (Additional
    references, step 1). Each index looks up only hashes up to its own threshold. A single
    index goes straight through, so its results are unchanged."""

    def __init__(self, indexes: Sequence[Index]) -> None:
        base = indexes[0]
        kind = molecule(base.meta["params"]["alphabet"])
        for other in indexes[1:]:
            if molecule(other.meta["params"]["alphabet"]) != kind:
                raise ValueError("nucleotide and protein indexes cannot be queried jointly")
            for key in ("k", "alphabet"):
                if other.meta["params"][key] != base.meta["params"][key]:
                    raise ValueError(f"indexes differ in {key}: cannot be queried jointly")
            if other.meta.get("hash") != base.meta.get("hash"):
                raise ValueError("indexes differ in hash scheme: cannot be queried jointly")
            if (other.dense is None) != (base.dense is None):
                raise ValueError("either every index or none must have a dense tier")
        self.indexes = list(indexes)
        heights = [i.units.height for i in indexes]
        self.offsets = np.cumsum([0, *heights])
        self.height = int(self.offsets[-1])
        self.dense = base.dense is not None
        self.walked: dict[str, list[_TierHits]] = {}

    def walk(self, tier: str, ranges: Iterable[np.ndarray]) -> None:
        """Read every index's ``tier`` once for the distinct sampled hashes, given as
        consecutive hash ranges, each sorted (:func:`tier_hits`, a range in memory at a
        time); :meth:`hits` then looks them up in memory."""
        column = "max_hash_g" if tier == "tier1" else "max_hash_dense"
        found: list[list[tuple[np.ndarray, np.ndarray]]] = [[] for _ in self.indexes]
        for part in ranges:
            for got, index in zip(found, self.indexes, strict=True):
                got.append(_key_sets(getattr(index, tier), part))
        self.walked[tier] = []
        for got, index in zip(found, self.indexes, strict=True):
            hit = np.concatenate([np.empty(0, np.uint64), *(f for f, _ in got)])
            sets = np.concatenate([np.empty(0, np.int64), *(s for _, s in got)])
            got.clear()
            self.walked[tier].append(_resolve(getattr(index, tier), index.units[column], hit, sets))

    def _hits(self, i: int, tier: str, hashes: np.ndarray, reads: np.ndarray) -> pl.DataFrame:
        if tier in self.walked:
            return self.walked[tier][i].hits(hashes, reads)
        index = self.indexes[i]
        column = "max_hash_g" if tier == "tier1" else "max_hash_dense"
        return unit_hits(getattr(index, tier), index.units[column], hashes, reads)

    def max_hash(self, tier: str) -> int:
        return max(int(getattr(i, tier).max_hash) for i in self.indexes)

    def hits(self, tier: str, hashes: np.ndarray, reads: np.ndarray) -> pl.DataFrame:
        """:func:`unit_hits` in every index's ``tier`` ("tier1" or "dense"), units offset;
        ``holders`` counts the units holding the k-mer over all indexes."""
        if len(self.indexes) == 1:
            return self._hits(0, tier, hashes, reads)
        parts = []
        for i, (index, offset) in enumerate(zip(self.indexes, self.offsets[:-1], strict=True)):
            ok = hashes <= np.uint64(getattr(index, tier).max_hash)
            got = self._hits(i, tier, hashes[ok], reads[ok])
            parts.append(got.with_columns(pl.col("unit").cast(pl.UInt64) + int(offset)))
        hits = pl.concat(parts)
        per_index = hits.group_by("hash", "unit").agg(pl.col("holders").first())
        holders = (
            hits.select("hash", "unit")
            .unique()
            .join(per_index, on=["hash", "unit"])
            .group_by("hash")
            .agg(total=pl.len())
        )
        return (
            hits.join(holders, on="hash", maintain_order="left")
            .with_columns(holders=pl.col("total").cast(hits.schema["holders"]))
            .drop("total")
        )

    def rows(self, units: np.ndarray) -> pl.DataFrame:
        """Unit table rows of (sorted) joint ``units``; with several indexes, ``source`` is
        each unit's index (0 = the first)."""
        if len(self.indexes) == 1:
            return self.indexes[0].units.rows(units)
        at = np.searchsorted(units, self.offsets)
        parts = [
            index.units.rows(units[lo:hi] - offset).with_columns(
                pl.col("unit").cast(pl.UInt64) + int(offset), source=pl.lit(i, pl.UInt32)
            )
            for i, (index, offset, lo, hi) in enumerate(
                zip(self.indexes, self.offsets[:-1], at[:-1], at[1:], strict=True)
            )
        ]
        return pl.concat(parts, how="diagonal_relaxed")


def gather(kmers: pl.DataFrame, t_g: np.ndarray) -> tuple[pl.DataFrame, np.ndarray]:
    """Assign each hit k-mer to one unit, greedily, as ``sourmash gather`` does.

    Repeatedly take the unit with the most unassigned hit k-mers, scaled by ``1 / t_g`` so
    units sampled at different rates compare as estimated k-mer counts, and give it those
    k-mers. Units left with none are explained away and get no row. Returns ``unit``,
    ``kmers_unique`` (k-mers assigned) and ``gather_rank`` (0 = taken first), and per row of
    ``kmers`` whether its unit took its k-mer (the k-mer's first holder in gather order).
    """
    got = _core.gather(
        kmers["unit"].cast(pl.UInt32).to_numpy(),
        kmers["hash"].cast(pl.UInt64).to_numpy(),
        np.asarray(t_g, dtype=np.float64),
    )
    owned = got.pop("owned")
    taken = pl.DataFrame(got).with_columns(
        gather_rank=pl.int_range(len(got["unit"]), dtype=pl.UInt32)
    )
    return taken, owned


def assign_best(kmers: pl.DataFrame, score: pl.Expr) -> pl.DataFrame:
    """Give each hit k-mer to the unit holding it with the highest ``score`` (ties: lowest).

    ``score`` is evaluated per row of ``kmers`` (``unit``, ``hash`` and any other columns)
    and must be constant within a unit. Static unit scores make this one pass, unlike
    :func:`gather`, whose scores fall as k-mers are taken. The rules it serves:

    - winner-take-all, as sylph: score = containment (k-mers hit / ``m_g``);
    - uniqueness-first: score = sum of 1 / (units holding the k-mer) over the unit's hit
      k-mers, divided by ``t_g`` so units sampled at different rates compare fairly.

    Returns the rows of ``kmers`` that won.
    """
    return (
        kmers.with_columns(_score=score)
        .sort(["hash", "_score", "unit"], descending=[False, True, False])
        .unique("hash", keep="first", maintain_order=True)
        .drop("_score")
    )


WTA_SCORE: Final = pl.col("hash").n_unique().over("unit") / pl.col("m_g")
UFIRST_SCORE: Final = (1 / pl.len().over("hash")).sum().over("unit") / pl.col("t_g")


EXPLAINED_AWAY: Final = 1e-3  # expected hits below which EM reports coverage 0


def _div(num: np.ndarray, den: np.ndarray) -> np.ndarray:
    """num / den, 0 where den is 0 (units or k-mers whose weight has vanished)."""
    out = np.zeros_like(num, dtype=np.float64)
    np.divide(num, den, out=out, where=den > 0)
    return out


# Random EM starts (``query --em-start``), to check convergence where there is no truth:
# coverages times a log-normal factor of this sd (x 0.37-2.7 at one sd), and presence or
# scale uniform on [START_PRESENT_MIN, 1].
START_SPREAD: Final = 1.0
START_PRESENT_MIN: Final = 0.3
# The zero-inflated fit starts from the plain EM's coverages (``query --zi-cold-start``: from
# the usual start, as before): the same fixed point in 5-60% fewer steps on synthetic families.
ZI_WARM_START = True


def em(
    kmers: pl.DataFrame,
    m_g: np.ndarray,
    *,
    zero_inflated: bool = False,
    prior: tuple[float, float] | None = None,
    tol: float = 1e-8,
    max_iter: int | None = None,
    report: dict[str, float] | None = None,
    start: np.random.Generator | None = None,
    init: np.ndarray | None = None,
) -> pl.DataFrame:
    """Per-unit k-mer ``coverage`` (and ``present`` fraction) by EM over k-mer hit counts.

    ``kmers`` has one row per (``unit``, ``hash``) with the k-mer's ``hits``. Each k-mer's
    hits are Poisson with mean the sum of ``coverage`` over the units holding it, and each
    unit's ``m_g`` kept k-mers (hit or not) all count in its expectation; EM finds the
    maximum-likelihood coverages. Components share no k-mers, so each is fitted and
    stops on its own (:func:`_fit_components`); ``tol`` 1e-8 per component matches the
    accuracy the joint fit had at 1e-6.

    ``zero_inflated`` lets only a fraction ``present`` of a unit's kept k-mers occur in the
    sample (the rest are structural zeros: the sample's strain differs there), so coverage
    is fitted to the k-mers present. Hits are shared in proportion to coverage x present,
    and ``present`` is the unit's share of hit k-mers over the kept k-mers expected to be
    hit at that coverage, capped at 1. Without shared k-mers the fixed point is the
    zero-truncated Poisson MLE sylph uses: mean hits per hit k-mer = c / (1 - exp(-c)).
    With no excess zeros ``present`` stays 1 and the result equals plain EM. Units left
    with fewer than ``EXPLAINED_AWAY`` expected hits (their k-mers explained by other units)
    get coverage 0. ``max_iter`` defaults to 1,000 steps, ``ZI_MAX_ITER`` when zero-inflated.

    ``prior`` (a, b >= 1) puts a Beta prior on ``present`` (see :func:`fit_present_prior`),
    which pulls it towards the prior mean when few k-mers could be hit, i.e. at low coverage
    or small ``m_g``; ``present`` is then updated by EM over the unhit k-mers' presence.

    ``start`` draws the starting coverages and ``present`` at random instead (see
    :data:`START_SPREAD`): fits from several starts agree where the fit converged.
    ``init`` (per unit, indexed as ``m_g``) replaces the starting coverages, e.g. a plain
    fit's to start the zero-inflated one; each is kept above 1e-3 x its usual start, since a
    coverage of 0 could never grow (a unit plain EM explains away may be present here).

    ``report`` collects ``em_iterations`` (the most any component took),
    ``em_unconverged_units`` (units of components stopped at ``max_iter``), summed over
    calls, and ``em_max_change`` (how far the furthest of those was from converging: see
    :func:`_fit_components`).
    """
    kmers = kmers.sort("unit", "hash")  # sums in a fixed order: results independent of input order
    units, col = _ids(kmers["unit"].to_numpy())
    hashes, row = _ids(kmers["hash"].to_numpy())
    hits = np.zeros(len(hashes))
    hits[row] = kmers["hits"].to_numpy()
    del kmers, hashes  # the sorted copy is not needed through the fit
    m = m_g[units].astype(np.float64)
    lam = np.bincount(col, weights=hits[row], minlength=len(units)) / m
    if init is not None:
        lam = np.maximum(init[units], 1e-3 * lam)
    pi = np.ones(len(units))
    if start is not None:  # a random starting point (:data:`START_SPREAD`)
        lam *= np.exp(start.normal(0.0, START_SPREAD, len(units)))
        if zero_inflated:
            pi = start.uniform(START_PRESENT_MIN, 1.0, len(units))

    def make_step(b: _Block) -> Step:
        h, mb, per_unit = hits[b.kmer], m[b.unit], len(b.unit)

        def over_units(x: np.ndarray) -> np.ndarray:  # per unit, sum of x over its k-mers
            return np.bincount(b.c, weights=x[b.r], minlength=per_unit)

        def step(la: np.ndarray, p: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
            w = la * p
            mu = np.bincount(b.r, weights=w[b.c], minlength=len(b.kmer)) + b.offset
            attributed = w * over_units(_div(h, mu))  # expected hits from each unit
            if zero_inflated:
                kmers_hit = w * over_units(_div(np.ones_like(mu), mu))  # expected hit k-mers
                seen = -np.expm1(-la)  # chance a present k-mer is hit
                if prior is None:
                    new_p = np.minimum(1.0, _div(kmers_hit, mb * seen))
                else:
                    # Unhit k-mers are present with odds p (1 - seen) : (1 - p).
                    odds = p * (1 - seen) / np.maximum(1 - p * seen, 1e-300)
                    unhit_present = np.maximum(mb - kmers_hit, 0) * odds
                    new_p = (kmers_hit + unhit_present + prior[0] - 1) / (mb + sum(prior) - 2)
            else:
                new_p = p
            return _div(attributed, mb * new_p), new_p, attributed

        return step

    if max_iter is None:
        max_iter = ZI_MAX_ITER if zero_inflated else 1000
    attributed = _fit_components(
        col,
        row,
        lam,
        pi,
        make_step,
        tol,
        max_iter,
        report,
        pair_weight=lambda p: lam[col[p]] * pi[col[p]],
    )
    # Units explained away by others converge towards 0 without reaching it.
    lam[attributed < EXPLAINED_AWAY] = 0.0
    return pl.DataFrame(
        {"unit": units, "coverage": lam, "present": pi},
        schema={"unit": pl.UInt32, "coverage": pl.Float64, "present": pl.Float64},
    )


# p_in of each quantised level (pin_q = round(15 p_in)), kept inside (0, 1): a k-mer certain to
# be present (p_in = 1 at s_g = 1) would stay so forever, since EM then reads every unhit copy
# of it as present; one certain to be absent could never be hit.
PIN_P: Final = np.clip(
    np.arange(2**PIN_BITS) / (2**PIN_BITS - 1), 0.5 / (2**PIN_BITS - 1), 1 - 0.5 / (2**PIN_BITS - 1)
)


def em_pin(
    kmers: pl.DataFrame,
    pin_hist: np.ndarray,
    *,
    tol: float = 1e-8,
    max_iter: int = 1000,
    start: np.random.Generator | None = None,
) -> pl.DataFrame:
    """Zero-inflated EM where each k-mer's presence follows its ``p_in``.

    k-mer x of unit g is present in the sample with probability s_g * p_in(x): a strain
    like a random member carries each k-mer with probability p_in (s_g = 1), and a more
    divergent one loses k-mers uniformly (s_g < 1). So core k-mers are expected present
    more often than private ones, which changes both how shared hits are split (by
    coverage x presence of that k-mer in each unit) and how much an unhit k-mer argues for
    low coverage rather than absence. ``p_in`` is read from its quantised level (``PIN_P``);
    s_g is capped so presence stays <= 1.

    ``kmers`` has ``unit``, ``hash``, ``hits`` and ``pin_q``; ``pin_hist[unit]`` counts the
    unit's kept k-mers per level (its sum is ``m_g``). Coverage is attributed hits
    over expected present k-mers (hit ones by their share, unhit ones by their posterior
    presence), and s_g is expected present k-mers over the sum of p_in. Like :func:`em`, a
    hit k-mer counts as present for a unit by its share of the hits. Returns ``unit``,
    ``coverage`` and ``present`` (expected fraction of kept k-mers present); units explained
    away get coverage 0. ``start`` draws a random starting point, as in :func:`em`.
    """
    kmers = kmers.sort("unit", "hash")  # sums in a fixed order: results independent of input order
    units, col = _ids(kmers["unit"].to_numpy())
    hashes, row = _ids(kmers["hash"].to_numpy())
    level = kmers["pin_q"].to_numpy().astype(np.intp)
    hits = np.zeros(len(hashes))
    hits[row] = kmers["hits"].to_numpy()
    del kmers, hashes  # the sorted copy is not needed through the fit
    n = len(units)
    hist = pin_hist[units].astype(np.float64)
    levels = hist.shape[1]
    m, expected = hist.sum(axis=1), hist @ PIN_P  # kept k-mers; present ones at s_g = 1
    lam, scale = np.bincount(col, weights=hits[row], minlength=n) / m, np.ones(n)
    if start is not None:  # a random starting point, as in :func:`em`
        lam *= np.exp(start.normal(0.0, START_SPREAD, n))
        scale = start.uniform(START_PRESENT_MIN, 1.0, n)

    def make_step(b: _Block) -> Step:
        h, lv, hs, ex = hits[b.kmer], level[b.pairs], hist[b.unit], expected[b.unit]
        cell = b.c * levels + lv

        def step(la: np.ndarray, sc: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
            weight = la[b.c] * sc[b.c] * PIN_P[lv]
            mu = np.bincount(b.r, weights=weight, minlength=len(b.kmer)) + b.offset
            share = _div(weight, mu[b.r])
            att = np.bincount(b.c, weights=share * h[b.r], minlength=len(b.unit))
            hit = np.bincount(cell, weights=share, minlength=hs.size).reshape(hs.shape)
            seen = -np.expm1(-la)[:, None]  # chance a present k-mer is hit
            pi = sc[:, None] * PIN_P
            odds = pi * (1 - seen) / np.maximum(1 - pi * seen, 1e-300)
            present = hit.sum(axis=1) + (np.maximum(hs - hit, 0) * odds).sum(axis=1)
            return _div(att, present), np.clip(_div(present, ex), 0.0, 1 / PIN_P[-1]), att

        return step

    # tol 1e-8 per component matches the joint fit's accuracy at 1e-6 (phase 6, step 15).
    attributed = _fit_components(
        col,
        row,
        lam,
        scale,
        make_step,
        tol,
        max_iter,
        other_max=1 / PIN_P[-1],
        pair_weight=lambda p: lam[col[p]] * scale[col[p]] * PIN_P[level[p]],
    )
    lam[attributed < EXPLAINED_AWAY] = 0.0
    return pl.DataFrame(
        {"unit": units, "coverage": lam, "present": scale * expected / m},
        schema={"unit": pl.UInt32, "coverage": pl.Float64, "present": pl.Float64},
    )


def fit_present_prior(
    fit: pl.DataFrame, min_coverage: float = 3.0, min_units: int = 10
) -> tuple[float, float] | None:
    """Beta (a, b) matching the mean and variance of ``present`` over well-covered units.

    ``fit`` is a zero-inflated :func:`em` result; at ``coverage >= min_coverage`` nearly all
    present k-mers are hit, so ``present`` is measured, not guessed. None if fewer than
    ``min_units`` qualify. a, b are floored at 1 (a unimodal prior).
    """
    well = fit.filter(pl.col("coverage") >= min_coverage)["present"]
    if well.len() < min_units:
        return None
    mean, var = float(well.mean()), float(well.var())  # type: ignore[arg-type]
    strength = mean * (1 - mean) / max(var, 1e-9) - 1
    return max(mean * strength, 1.0), max((1 - mean) * strength, 1.0)


# Background detections per unit per read (pair) per unit of sampling rate t_g, and the
# geometric ratio of own k-mers per background detection (a homolog gives a few at once);
# calibrated on the fmh benchmark (step 18): P(true | own k-mers) within 0.05 on all indexes.
BACKGROUND: Final = 1e-6
CLUMP: Final = 0.3
# The most of the units gather keeps that background may claim, in expectation
# (query --background-share): BACKGROUND x every index unit predicts 2e8-1e9 background
# detections on the MGnify90 standard tier, against 3.2 M units kept at 10 M pairs (plan,
# phase 13, 2026-10-20); on KO-level indexes the cap does not bind. 0.1 is the user's
# choice (2026-10-20), to calibrate against spikes on the full standard tier.
BACKGROUND_SHARE: Final = 0.1
H_CAP: Final = 20  # hit counts pooled at and above this in the present-unit distribution


def presence(
    own: pl.DataFrame,
    t_g: np.ndarray,
    n_reads: int,
    n_index_units: int,
    *,
    background: float = BACKGROUND,
    clump: float = CLUMP,
    share: float = BACKGROUND_SHARE,
    max_iter: int = 500,
    tol: float = 1e-9,
    report: dict[str, float] | None = None,
) -> pl.DataFrame:
    """Probability that each unit gather keeps is present rather than hit by background.

    ``own`` has one row per k-mer gather gave a unit (``unit``, ``hash``, ``holders``), so
    hits other units explain are not evidence. An absent unit is still detected by
    background (off-target homologs, error k-mers) with probability 1 - exp(-mu_g), mu_g =
    ``background`` x ``n_reads`` x ``t_g``, then with h own k-mers ~ Geometric (ratio
    ``clump``), each hit k-mer's odds scaled by its index ``holders`` (k-mers many units
    share are conserved motifs, which unindexed genes carry too). A present unit has h
    drawn from a distribution f fitted to the sample.
    With A absent units (``n_index_units`` minus the expected present), P(present | h) =
    w_h / (w_h + A x P(h | absent)), w_h the expected present units with h own k-mers; w
    and A are fitted by fixed-point iteration, at which w_h is about the units seen with h
    own k-mers less the background expected to give h.

    The background rate is calibrated, not fitted: fitted per sample it is not
    identifiable (a free f, or a free per-unit hit rate, absorbs any unit with few hits,
    so fits collapse to no background or to everything being background). Scaling by
    t_g rather than m_g matches the fmh benchmark, where false positives are small,
    often floored, KOs.

    The rate is per index unit, so the background it predicts grows with the index; it is
    capped at ``share`` of the units seen (``own``'s): the absent term is scaled down until
    A x mean(1 - exp(-mu)) over the units seen is at most ``share`` x their number (``inf``:
    no cap). ponytail: the units seen stand in for the index's t_g distribution; store the
    index's mean of 1 - exp(-mu) per read depth at build if the two differ.

    ``report`` gets ``presence_iterations``, ``presence_converged`` (0 when stopped at
    ``max_iter``), ``presence_background_scale`` (the cap's factor, 1 when it does not bind)
    and ``presence_background_units`` (the expected background units, sum of 1 - P).
    """
    per_unit = (
        own.group_by("unit")
        .agg(h=pl.len(), log_holders=pl.col("holders").cast(pl.Float64).log().sum())
        .sort("unit")
    )
    units = per_unit["unit"].to_numpy()
    h = per_unit["h"].to_numpy().astype(np.float64)
    mu = background * n_reads * t_g[units]
    capped = np.minimum(h, H_CAP).astype(np.intp)
    # P(own hits | absent); units at the cap are never taken for background.
    absent = np.where(
        capped < H_CAP,
        -np.expm1(-mu)
        * (1 - clump)
        * clump ** (h - 1)
        * np.exp(per_unit["log_holders"].to_numpy()),
        0.0,
    )
    detect = float(np.mean(-np.expm1(-mu))) if len(units) else 0.0
    prob = np.ones(len(units))
    it, done, scale = 0, False, 1.0
    while not done and it < max_iter:
        w = np.bincount(capped, weights=prob, minlength=H_CAP + 1)
        n_absent = n_index_units - prob.sum()
        predicted = n_absent * detect  # background units the calibrated rate predicts
        scale = min(1.0, share * len(units) / predicted) if predicted > 0 else 1.0
        new = w[capped] / (w[capped] + scale * n_absent * absent)
        done = np.abs(new - prob).max(initial=0) <= tol
        prob = new
        it += 1
    if report is not None:
        report["presence_iterations"] = it
        report["presence_converged"] = int(done)
        report["presence_background_scale"] = scale
        report["presence_background_units"] = float((1 - prob).sum())
    return pl.DataFrame(
        {"unit": units, "present_prob": prob},
        schema={"unit": pl.UInt32, "present_prob": pl.Float64},
    )


def fdr_q(prob: np.ndarray) -> np.ndarray:
    """Each unit's q-value: the least expected false share, mean(1 - P), among the units
    ranked by presence probability down to it (ties in input order), over the units given.
    Keeping units with q <= x keeps the largest top set whose expected false share is <= x."""
    order = np.argsort(-prob, kind="stable")
    fdr = np.cumsum(1 - prob[order]) / np.arange(1, len(prob) + 1)
    q = np.empty(len(prob))
    q[order] = np.minimum.accumulate(fdr[::-1])[::-1]
    return q


# Floor of the copies error (sd on the log scale), calibrated on the fmh benchmark (step 18):
# brings abundance_zi intervals to 95% at kfp_s100.
COPIES_ERROR: Final = 0.2
MH_STEPS: Final = 5  # Metropolis steps on coverage per sweep
PI_GRID: Final = (np.arange(512) + 0.5) / 512  # present-fraction grid for exact draws


def _ztp_log_post(lam: np.ndarray, hits: np.ndarray, hit: np.ndarray) -> np.ndarray:
    """Log posterior of coverage from ``hits`` over ``hit`` k-mers, zero-truncated Poisson."""
    out: np.ndarray = hits * np.log(lam) - hit * (lam + np.log(-np.expm1(-lam))) - 0.01 * lam
    return out


def posterior_zi(
    pairs: pl.DataFrame,
    hash_reads: pl.DataFrame,
    m_g: np.ndarray,
    pin_sum: np.ndarray,
    draws: int,
    *,
    present_prob: np.ndarray | None = None,
    len_cv: np.ndarray | None = None,
    copies_error: float = COPIES_ERROR,
    sweeps: int = 10,
    level: float = 0.95,
    shared_evidence: float = 0.5,
    seed: int = 0,
    workers: int | None = None,
) -> pl.DataFrame:
    """Intervals and ambiguity groups for the zero-inflated model: resampled reads + Gibbs.

    ``pairs`` has one row per hit (``unit``, ``hash``) and ``hash_reads`` one per (``hash``,
    ``read``) with the ``n`` times the read hit the k-mer (every unit holding a k-mer sees the
    same reads, so they are stored once per k-mer). Each of ``draws`` replicates

    1. reweights every read (pair) by a Poisson(1) draw and rebuilds the k-mer counts (one
       read hits many neighbouring k-mers, so counts are not independent), and refits
       zero-inflated :func:`em`: this carries the read-sampling uncertainty;
    2. runs ``sweeps`` Gibbs sweeps from that fit on those counts and keeps the last state,
       which adds the uncertainty of splitting shared k-mers between similar units. A sweep
       splits each hit k-mer's count among its holders multinomially, in proportion to
       coverage x present fraction (EM's shares, approximating the holders' joint presence);
       updates coverage from the counts on the k-mers each unit was given, a zero-truncated
       Poisson likelihood free of the present fraction (Metropolis on log coverage,
       Gamma(1, 0.01) prior); and draws the present fraction exactly on a grid from hit
       k-mers ~ Binomial(m, pi (1 - e^-lambda)), plus the present unhit k-mers for copies.

    A unit given no hits in a draw is absent in it (coverage and abundance 0), and so is one
    dropped with probability 1 - ``present_prob[unit]`` (:func:`presence`), so weak or
    background-like evidence widens intervals down to 0. Copies (present k-mers over
    ``pin_sum``) assume an average member; each draw scales them by exp(N(0, s^2)), s^2 =
    log(1 + ``len_cv[unit]``^2 / copies) + ``copies_error``^2, the spread of the mean
    kept k-mers of that many members plus a calibrated floor. Where units
    cannot be told apart, their shared counts move between them across draws, so their
    intervals widen. Ambiguity groups come from shared evidence: over all draws, a unit is
    linked to another when at least ``shared_evidence`` of the hits allocated to it lie on
    k-mers both hold, and linked units form a group. This catches lopsided pairs, e.g. a
    small unit living on hits it shares with a well-supported one, whose draws barely
    correlate. ``own_evidence`` is a unit's share of allocated hits on k-mers no other
    detected unit holds. A group's total is usually far better determined than its members.
    Returns per unit ``coverage_zi_lo``/``_hi``, ``abundance_zi_lo``/``_hi`` at ``level``,
    ``own_evidence``, and ``ambiguity_group`` (the group's smallest unit id; null when
    alone), ``group_size`` and the group totals ``group_coverage_zi_lo``/``_hi`` and
    ``group_abundance_zi_lo``/``_hi``.

    Units in different components (linked by shared k-mers) never interact, so the draws run
    per batch of components (:func:`component_batches`) sized to ``POSTERIOR_BATCH_BYTES``,
    each with its own generator seeded by ``seed`` and the batch's smallest unit. A read's
    Poisson weight is a function of ``seed``, the draw and its id (:func:`_poisson1`), so
    reads hitting several batches weigh the same in each. So the batches run in parallel on
    ``workers`` processes (default ``pl.thread_pool_size()``, i.e. ``POLARS_MAX_THREADS``)
    with the same result as one after another.
    """
    per_pair = 1600 + 16 * draws  # bytes per (unit, hash) row at peak, measured (step 28)
    parts = component_batches(
        pairs.select("unit", "hash"), max(1, POSTERIOR_BATCH_BYTES // per_pair)
    )
    batch = partial(
        _posterior_batch,
        m_g=m_g,
        pin_sum=pin_sum,
        draws=draws,
        present_prob=present_prob,
        len_cv=len_cv,
        copies_error=copies_error,
        sweeps=sweeps,
        level=level,
        shared_evidence=shared_evidence,
        seed=seed,
    )

    def rows(part: pl.DataFrame) -> pl.DataFrame:
        """This batch's k-mers' rows: a scan per batch, no copy of all rows."""
        return hash_reads.join(part.select("hash").unique(), on="hash", how="semi")

    workers = min(workers or pl.thread_pool_size(), len(parts))
    if workers <= 1:
        return pl.concat([batch(part, rows(part)) for part in parts]).sort("unit")
    # Spawn, not fork: Polars' thread pool does not survive a fork. Each worker is one core,
    # so its Polars runs one thread (children read the environment when they start).
    # ponytail: the stage's peak_rss/peak_anon are the parent's only; a worker adds about
    # one batch (POSTERIOR_BATCH_BYTES, more for a larger component), which Slurm counts.
    out: list[pl.DataFrame] = []
    pending: set[Future[pl.DataFrame]] = set()
    polars_threads = os.environ.get("POLARS_MAX_THREADS")
    os.environ["POLARS_MAX_THREADS"] = "1"
    try:
        with ProcessPoolExecutor(
            workers,
            mp_context=multiprocessing.get_context("spawn"),
            initializer=_posterior_worker_init,
            initargs=(batch,),
        ) as pool:
            # Largest first: a giant component's batch is the critical path, so start it early.
            for part in sorted(parts, key=lambda part: -part.height):
                if len(pending) >= 2 * workers:  # a few batches' rows built ahead, not all
                    done, pending = wait(pending, return_when=FIRST_COMPLETED)
                    out += [f.result() for f in done]
                pending.add(pool.submit(_posterior_worker_batch, part, rows(part)))
            out += [f.result() for f in pending]
    finally:
        if polars_threads is None:
            del os.environ["POLARS_MAX_THREADS"]
        else:
            os.environ["POLARS_MAX_THREADS"] = polars_threads
    return pl.concat(out).sort("unit")


_posterior_worker: Callable[[pl.DataFrame, pl.DataFrame], pl.DataFrame] | None = None


def _posterior_worker_init(batch: Callable[[pl.DataFrame, pl.DataFrame], pl.DataFrame]) -> None:
    """A worker's :func:`_posterior_batch` with posterior_zi's arguments, sent once."""
    global _posterior_worker
    _posterior_worker = batch


def _posterior_worker_batch(part: pl.DataFrame, rows: pl.DataFrame) -> pl.DataFrame:
    assert _posterior_worker is not None
    return _posterior_worker(part, rows)


POSTERIOR_BATCH_BYTES: Final = 2**29  # working memory of one batch of the posterior's draws
SWEEP_BLOCK_PAIRS: Final = 2**20  # (unit, hash) entries per block of a Gibbs sweep
READ_CHUNK: Final = 2**22  # (hash, read) rows reweighted at once per draw

# Poisson(1) CDF at 0..24; beyond, a probability below 1e-24.
POISSON1_CDF: Final = np.cumsum(np.exp(-1.0) / np.cumprod(np.r_[1.0, np.arange(1.0, 25.0)]))


def _poisson1(seed: int, draw: int, reads: np.ndarray) -> np.ndarray:
    """Poisson(1) weight of each read (id) in a draw: a SplitMix64 hash of (``seed``,
    ``draw``, read) to a uniform, then the inverse CDF, so the same read has the same weight
    wherever it is drawn."""
    with np.errstate(over="ignore"):
        x = np.asarray(reads, dtype=np.uint64) + np.uint64(
            (seed * 1_000_003 + draw) % 2**64
        ) * np.uint64(0x9E3779B97F4A7C15)
        for shift, mult in ((30, 0xBF58476D1CE4E5B9), (27, 0x94D049BB133111EB)):
            x = (x ^ (x >> np.uint64(shift))) * np.uint64(mult)
        x ^= x >> np.uint64(31)
    uniform = (x >> np.uint64(11)).astype(np.float64) * 2.0**-53
    return np.searchsorted(POISSON1_CDF, uniform, side="right")


def _posterior_batch(
    pairs: pl.DataFrame,
    hash_reads: pl.DataFrame,
    m_g: np.ndarray,
    pin_sum: np.ndarray,
    draws: int,
    *,
    present_prob: np.ndarray | None,
    len_cv: np.ndarray | None,
    copies_error: float,
    sweeps: int,
    level: float,
    shared_evidence: float,
    seed: int,
) -> pl.DataFrame:
    """:func:`posterior_zi` on one batch of whole components."""
    # Sorted, so the seeded draws do not depend on row order (Polars group_by does not fix it).
    keyed = pairs.select("unit", "hash").sort("unit", "hash")
    units, col = np.unique(keyed["unit"].to_numpy(), return_inverse=True)
    rng = np.random.default_rng([seed, int(units[0])])
    hashes, row = np.unique(keyed["hash"].to_numpy(), return_inverse=True)
    n_units, n_rows = len(units), len(hashes)
    by_hash = hash_reads.join(pl.DataFrame({"hash": hashes}), on="hash", how="semi").sort(
        "hash", "read"
    )
    hash_of_hit: np.ndarray = np.searchsorted(hashes, by_hash["hash"].to_numpy()).astype(
        np.int32 if n_rows < 2**31 else np.int64
    )
    reads, read_of = _ids(by_hash["read"].to_numpy())
    n = by_hash["n"].to_numpy().astype(np.int64)
    del by_hash
    # Entries sorted by k-mer, with their position within it, for sequential binomial splits.
    order = np.argsort(row, kind="stable")
    row_s, col_s = row[order], col[order]
    starts = np.searchsorted(row_s, np.arange(n_rows))
    position = np.arange(len(row_s)) - starts[row_s]
    holders = np.bincount(row_s, minlength=n_rows)
    m = m_g[units].astype(np.float64)
    keep_prob = np.ones(n_units) if present_prob is None else present_prob[units]
    cv2 = np.zeros(n_units) if len_cv is None else len_cv[units] ** 2
    coverage = np.zeros((draws, n_units))
    abundance, group_coverage, group_abundance = (np.zeros_like(coverage) for _ in range(3))
    evidence = np.zeros(len(row_s))  # hits allocated per (k-mer, unit) entry, all draws
    # Whole k-mers in blocks of about SWEEP_BLOCK_PAIRS entries, swept one block at a time.
    cuts = np.unique(np.r_[starts[row_s[::SWEEP_BLOCK_PAIRS]], len(row_s)])
    for b in range(draws):
        per_hash = np.zeros(n_rows)
        for lo in range(0, len(read_of), READ_CHUNK):  # integer sums: exact in any order
            part = slice(lo, lo + READ_CHUNK)
            weight = n[part] * _poisson1(seed, b, reads[read_of[part]])
            per_hash += np.bincount(hash_of_hit[part], weights=weight, minlength=n_rows)
        per_hash = np.rint(per_hash)
        per_key = per_hash[row]
        table = keyed.select("unit", "hash").with_columns(hits=per_key).filter(pl.col("hits") > 0)
        if table.height == 0:
            continue
        fit = em(table, m_g, zero_inflated=True)
        at = np.searchsorted(units, fit["unit"].to_numpy())
        lam, pi = np.full(n_units, 1e-6), np.full(n_units, 1e-6)
        lam[at] = np.maximum(fit["coverage"].to_numpy(), 1e-6)
        pi[at] = np.clip(fit["present"].to_numpy(), 1e-6, 1 - 1e-6)
        # Holders of a k-mer see the same reads, so any holder's count is the k-mer's.
        count = per_hash.astype(np.int64)
        for sweep in range(sweeps):
            unit_rate = lam * pi
            hits, hit = np.zeros(n_units), np.zeros(n_units)  # hits and k-mers given
            for lo, hi in itertools.pairwise(cuts):
                first = row_s[lo]
                r_b, c_b, pos_b = row_s[lo:hi] - first, col_s[lo:hi], position[lo:hi]
                held = holders[first : row_s[hi - 1] + 1]
                rate = unit_rate[c_b]
                rest = np.bincount(r_b, weights=rate, minlength=len(held))  # rate not yet visited
                left = count[first : first + len(held)].copy()
                given = np.zeros(hi - lo, dtype=np.int64)
                for j in range(int(held.max())):
                    here = pos_b == j
                    r = r_b[here]
                    share = np.where(held[r] == j + 1, 1.0, _div(rate[here], rest[r]))
                    given[here] = rng.binomial(left[r], np.clip(share, 0.0, 1.0))
                    left[r] -= given[here]
                    rest[r] -= rate[here]
                hits += np.bincount(c_b, weights=given, minlength=n_units)
                hit += np.bincount(c_b, weights=given > 0, minlength=n_units)
                if sweep == sweeps - 1:
                    evidence[lo:hi] += given
            alive = hit > 0
            for _ in range(MH_STEPS):
                step = rng.normal(0.0, 1.5 / np.sqrt(hits + 1))
                new = lam * np.exp(step)
                gain = _ztp_log_post(new, hits, hit) - _ztp_log_post(lam, hits, hit) + step
                lam = np.where(alive & (np.log(rng.random(n_units)) < gain), new, lam)
            seen = -np.expm1(-lam)
            u = rng.random(n_units)
            pick = np.empty(n_units, dtype=np.intp)
            # Units x grid arrays, a chunk of units at a time (8 GB per 2 M units at once).
            chunk = max(1, CHUNK_BYTES // (8 * len(PI_GRID)))
            for lo in range(0, n_units, chunk):
                part = slice(lo, lo + chunk)
                log_post = hit[part, None] * np.log(PI_GRID) + (m - hit)[part, None] * np.log1p(
                    -PI_GRID * seen[part, None]
                )
                cdf = np.cumsum(np.exp(log_post - log_post.max(axis=1, keepdims=True)), axis=1)
                pick[part] = (cdf < u[part, None] * cdf[:, -1:]).sum(axis=1)
            pi = PI_GRID[np.minimum(pick, len(PI_GRID) - 1)]
        unhit_odds = pi * (1 - seen) / (1 - pi * seen)
        present = hit + rng.binomial(np.maximum(m - hit, 0).astype(np.int64), unhit_odds)
        copies = present / pin_sum[units]
        sd = np.sqrt(np.log1p(cv2 / np.maximum(copies, 1.0)) + copies_error**2)
        scaled = lam * copies * np.exp(rng.normal(0.0, sd))
        # Group totals keep a dropped member's share: under absence its hits go to others.
        group_coverage[b] = np.where(alive, lam, 0.0)
        group_abundance[b] = np.where(alive, scaled, 0.0)
        alive &= rng.random(n_units) < keep_prob
        coverage[b] = np.where(alive, lam, 0.0)
        abundance[b] = np.where(alive, scaled, 0.0)

    # Ambiguity groups: link a unit to another holder of its hit k-mers when most of its
    # allocated hits lie on k-mers the two share.
    total_evidence = np.bincount(col_s, weights=evidence, minlength=n_units)
    own = _div(np.bincount(col_s, weights=evidence * (holders[row_s] == 1), minlength=n_units),
               total_evidence)  # fmt: skip
    entries = pl.DataFrame({"row": row_s, "unit": col_s, "evidence": evidence})
    shared = (
        entries.join(entries.select("row", other="unit"), on="row")
        .filter(pl.col("unit") != pl.col("other"))
        .group_by("unit", "other")
        .agg(pl.col("evidence").sum())
    )
    unit_a, unit_b = shared["unit"].to_numpy(), shared["other"].to_numpy()
    # A unit never given hits (explained away in every draw) has total 0, so it links to
    # every co-holder of its hit k-mers.
    link = shared["evidence"].to_numpy() >= shared_evidence * total_evidence[unit_a]
    graph = coo_array(
        (np.ones(int(link.sum())), (unit_a[link], unit_b[link])),
        shape=(n_units, n_units),
    )
    _, group = connected_components(graph, directed=False)
    size = np.bincount(group)[group]
    group_id = np.full(n_units, len(units), dtype=np.int64)
    np.minimum.at(group_id, group, np.arange(n_units))  # smallest member index per group
    group_unit = units[group_id[group]]
    total = [np.stack([np.bincount(group, w, minlength=n_units) for w in x])[:, group]
             for x in (group_coverage, group_abundance)]  # fmt: skip
    total = [np.where(size > 1, t, x) for t, x in zip(total, (coverage, abundance), strict=True)]
    tails = [(1 - level) / 2, (1 + level) / 2]
    columns = {"unit": units.astype(np.uint32)}
    for name, x in (("coverage_zi", coverage), ("abundance_zi", abundance),
                    ("group_coverage_zi", total[0]), ("group_abundance_zi", total[1])):  # fmt: skip
        columns[f"{name}_lo"], columns[f"{name}_hi"] = np.quantile(x, tails, axis=0)
    return pl.DataFrame(columns).with_columns(
        ambiguity_group=pl.when(pl.Series(size) > 1).then(pl.Series(group_unit.astype(np.uint32))),
        group_size=pl.Series(size.astype(np.uint32)),
        own_evidence=pl.Series(own),
    )


MIN_AAI_KMERS: Final = 5  # default hit k-mers below which aai is null (--min-aai-kmers)
AAI_CHUNK: Final = 2**20  # units per aai_naive estimate
MIN_AAI_LAMBDA: Final = 0.1  # zero-truncated coverage below which aai_naive is not corrected


def count_cdf(k: np.ndarray, mean: np.ndarray, var: np.ndarray) -> np.ndarray:
    """P(X <= k) for a count X with this ``mean`` and ``var``: binomial (n = mean^2 /
    (mean - var)) when var < mean, negative binomial (r = mean^2 / (var - mean)) when
    var > mean, Poisson between. ``k`` may be fractional (the incomplete beta and gamma
    functions interpolate the CDF between counts); 0 for k < 0."""
    k = np.asarray(k, dtype=np.float64)
    mean = np.maximum(mean, 1e-300)
    gap = var - mean
    with np.errstate(divide="ignore", invalid="ignore"):
        n = mean**2 / np.where(gap < 0, -gap, 1.0)
        binom = np.where(k >= n, 1.0, betainc(np.maximum(n - k, 1e-300), k + 1, 1 - mean / n))
        r = mean**2 / np.where(gap > 0, gap, 1.0)
        negbin = betainc(r, k + 1, r / (r + mean))
    exact = np.abs(gap) <= 1e-9 * mean
    cdf = np.where(exact, gammaincc(k + 1, mean), np.where(gap < 0, binom, negbin))
    return np.where(k < 0, 0.0, cdf)


def survival_quantiles(
    h: np.ndarray,
    pin_sum: np.ndarray,
    seen: np.ndarray,
    keep: np.ndarray,
    windows: np.ndarray,
    model: SurvivalModel,
    clumping: np.ndarray | float = 1.0,
    level: float = 0.95,
    selected: float = 1.0,
    seen_range: tuple[np.ndarray, np.ndarray] | None = None,
    n_members: np.ndarray | float = 1.0,
    bounds: bool = True,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """(point, lower, upper) identity from ``h`` hit k-mers of an average member's
    ``pin_sum`` kept windows (sampled at rate ``keep`` out of ``windows``), each window
    surviving with probability S(a) under ``model`` and hit with probability ``seen``.
    S is the union's of the unit's ``n_members`` (:meth:`SurvivalModel.union_at`; one
    member's without a union term), identity to the nearest member (step 34).

    At identity a (a row of the model's grid), h = sum_i x_i y_i over the n = ``pin_sum``
    windows, x_i survival (correlated: variance inflation d, :meth:`SurvivalModel.overlap_at`;
    neighbouring windows share sites, and under clustered rates windows a region apart
    survive together) and y_i a hit (Bernoulli(seen), in clumps of ``clumping``: a read hits
    a run of windows). So h is a count (:func:`count_cdf`) with mean n seen S and variance
    n seen S (1 - seen) + n seen^2 S (1 - S) d + (clumping - 1) n S^2 seen (1 - seen), plus
    (n seen S v)^2 for a union's scatter v (:meth:`SurvivalModel.union_scatter_at`). The
    model's survival mean and variance match simulation of it, skew included (step 25).

    An estimate is reported only if h >= ``selected`` (a unit is in the profile only if
    h >= 1, and ``aai`` is null below ``--min-aai-kmers``), so h's distribution is
    conditioned on that. The point is the median-unbiased estimate (mid-p: the a at which h
    is the conditional median); the bounds are where h is at the ``level`` tails' quantiles
    (a Neyman interval; with ``bounds`` False, the point only, three times). seen is
    estimated (from hits per hit k-mer), so the upper bound is
    taken at the low end of ``seen_range`` and the lower bound at its high end. Each is a
    bisection over the grid's rows (~11 steps). Without the conditioning, estimates read
    high near the detection limit (step 23).
    """
    tables, last = model.tables, len(model.tables["S"]) - 1
    n = np.maximum(pin_sum, 1e-9)
    n_members = np.asarray(n_members, dtype=np.float64)
    clump = np.maximum(np.asarray(clumping, dtype=np.float64), 1.0)

    def cdfs(row: np.ndarray, p: np.ndarray, k: np.ndarray = h) -> tuple[np.ndarray, np.ndarray]:
        """P(h' <= k | h' >= ``selected``) and P(h' < k | ...) at ``row``, hit chance ``p``."""
        surv = model.union_at(row, n_members)
        mean = n * p * surv
        # ponytail: one member's window correlation for the union's; a union of members
        # likely evens survival out (less correlated), so intervals err wide, not narrow
        d = model.overlap_at(row, keep, windows)
        var = n * p * surv * (1 - p) + n * p**2 * surv * (1 - surv) * d
        var += (clump - 1) * n * surv**2 * p * (1 - p)
        var += (n * p * surv * model.union_scatter_at(row, n_members)) ** 2  # union scatter
        zero = count_cdf(np.full_like(h, selected - 1), mean, var)
        scale = 1 - zero
        # selection all but impossible: given it, h' sits at the threshold, below any h
        rare = scale < 1e-12
        scale = np.where(rare, 1.0, scale)
        at = np.clip((count_cdf(k, mean, var) - zero) / scale, 0.0, 1.0)
        before = np.clip((count_cdf(k - 1, mean, var) - zero) / scale, 0.0, 1.0)
        return np.where(rare, 1.0, at), np.where(rare, (k >= selected + 1).astype(float), before)

    def first(test: Callable[[np.ndarray], np.ndarray]) -> np.ndarray:
        """The first row where ``test`` (false, then true as S falls) holds; ``last`` if none."""
        lo, hi = np.zeros(len(h), dtype=int), np.full(len(h), last)
        while np.any(lo < hi):
            mid = (lo + hi) // 2
            ok = test(mid) | (lo == hi)
            lo, hi = np.where(ok, lo, mid + 1), np.where(ok, mid, hi)
        return lo

    point = first(lambda row: np.add(*cdfs(row, seen)) / 2 >= 0.5)
    if not bounds:  # the point only, as the bounds
        return tables["a"][point], tables["a"][point], tables["a"][point]
    low_seen, high_seen = (seen, seen) if seen_range is None else seen_range
    tail = (1 - level) / 2
    upper = first(lambda row: cdfs(row, low_seen)[0] >= tail)
    # beyond the model's reach (more hits than identity 1 gives on average, as in large
    # unions, where hits over one member's pin_sum can exceed 1), the data say only "at the
    # top": the lower bound is taken as for h at that mean, not collapsed onto 1 (step 34)
    top = n * high_seen * model.union_at(np.zeros(len(h), dtype=int), n_members)
    h_low = np.minimum(h, np.maximum(top, selected))
    # the lowest S whose upper tail P(h' >= h) still holds tail: one row before it fails
    lower = np.maximum(first(lambda row: 1 - cdfs(row, high_seen, h_low)[1] < tail) - 1, 0)
    # ponytail: rows of the grid (steps ~0.003 in a), not interpolated between them
    a = tables["a"]
    return a[point], np.minimum(a[lower], a[point]), np.maximum(a[upper], a[point])


# build parameters an aai calibration depends on: a map fitted under others does not apply
AAI_CALIBRATION_PARAMS: Final = (
    "k", "alphabet", "t_base", "t_base_singleton", "n_min", "t_cap", "t_dense"
)  # fmt: skip


def _index_params(meta: dict[str, Any]) -> IndexParams:
    """Build parameters from ``meta.json``; fields older indexes lack take their defaults
    (format-1 indexes also record ``tier1_per_unit``, which is ignored)."""
    return IndexParams(
        **{f.name: meta["params"].get(f.name, f.default) for f in fields(IndexParams)}
    )


def check_aai_calibration(cal: dict[str, Any], params: dict[str, Any]) -> None:
    """Raise ValueError unless ``cal`` is a usable map for an index built with ``params``:
    ``aai`` strictly increasing knots, ``identity`` non-decreasing in [0, 1] of the same
    length, ``widen`` >= 0, and the build parameters it records (``params``, if any) equal
    to the index's in ``AAI_CALIBRATION_PARAMS``."""
    xs, ys = np.asarray(cal.get("aai", []), float), np.asarray(cal.get("identity", []), float)
    if len(xs) < 2 or len(xs) != len(ys):
        raise ValueError("aai and identity must be knots of equal length, at least 2")
    if np.any(np.diff(xs) <= 0) or np.any(np.diff(ys) < 0) or ys.min() < 0 or ys.max() > 1:
        raise ValueError("aai must increase strictly, identity must not decrease, in [0, 1]")
    if not cal.get("widen", 0.0) >= 0:
        raise ValueError("widen must be >= 0")
    fitted = cal.get("params", {})
    differ = [p for p in AAI_CALIBRATION_PARAMS if p in fitted and fitted[p] != params.get(p)]
    if differ:
        raise ValueError(f"calibration fitted on an index with other {', '.join(differ)}")


def calibrate_aai(profile: pl.DataFrame, cal: dict[str, Any]) -> pl.DataFrame:
    """``aai``, ``aai_lo`` and ``aai_hi`` through an inverse calibration map (knots
    ``aai`` -> ``identity``, linear between them, constant beyond; fitted by the benchmark's
    ``aai-calibrate``), the interval widened by ``widen`` and clipped to [0, 1]; a null
    ``aai`` stays null."""
    xs, ys = np.asarray(cal["aai"], float), np.asarray(cal["identity"], float)
    widen = {"aai": 0.0, "aai_lo": -cal["widen"], "aai_hi": cal["widen"]}
    return profile.with_columns(
        pl.col(c)
        .map_batches(lambda s: pl.Series(np.interp(s.to_numpy(), xs, ys)), return_dtype=pl.Float64)
        .add(d)
        .clip(0.0, 1.0)
        .fill_nan(None)
        for c, d in widen.items()
    )


def _calibrated(result: pl.DataFrame, calibrations: list[dict[str, Any] | None]) -> pl.DataFrame:
    """Each index's calibration on its own units (``source`` with joint indexes); the raw
    estimate kept as ``aai_raw``, ``aai_raw_lo``, ``aai_raw_hi``."""
    if "aai" not in result.columns or all(c is None for c in calibrations):
        return result
    result = result.with_columns(aai_raw="aai", aai_raw_lo="aai_lo", aai_raw_hi="aai_hi")
    if len(calibrations) == 1:
        return calibrate_aai(result, calibrations[0])  # type: ignore[arg-type]
    return pl.concat(
        calibrate_aai(part, cal) if cal is not None else part
        for cal, part in (
            (calibrations[i], result.filter(pl.col("source") == i))
            for i in range(len(calibrations))
        )
    ).sort("unit")


MIN_DISPERSION_KMERS: Final = 10  # unique hit k-mers below which a unit's dispersion is 1


def kmer_dispersion(kmers: pl.DataFrame, coverage: pl.DataFrame) -> pl.DataFrame:
    """Per unit, how much its k-mers' hit counts vary beyond a zero-truncated Poisson at its
    fitted ``coverage`` (``unit``, ``coverage``): the observed variance of hits over the hit
    k-mers only it holds in ``kmers`` (``unit``, ``hash``, ``hits``), over the zero-truncated
    variance E (1 + c - E), E = c / (1 - e^-c); floored at 1, and 1 with fewer than
    ``MIN_DISPERSION_KMERS`` such k-mers. Uneven read coverage along a gene, and strains at
    different depths in one unit, raise it (step 35); without it intervals of deep units
    were too narrow."""
    own = kmers.filter(pl.len().over("hash") == 1)
    moments = own.group_by("unit").agg(
        n=pl.len(), mean=pl.col("hits").mean(), var=pl.col("hits").var()
    )
    c = pl.col("coverage").clip(lower_bound=1e-9)
    expect = c / (1 - (-c).exp())
    ztp_var = expect * (1 + c - expect)
    return (
        coverage.select("unit", "coverage")
        .join(moments, on="unit", how="left")
        .select(
            "unit",
            dispersion=pl.when(pl.col("n") >= MIN_DISPERSION_KMERS)
            .then((pl.col("var") / ztp_var).clip(lower_bound=1.0))
            .otherwise(1.0)
            .fill_null(1.0)
            .fill_nan(1.0),
        )
    )


MIX_SHAPES: Final = (0.0, 0.02, 0.05, 0.1, 0.2, 0.35, 0.6)  # v = cv^2 of uneven coverage
MIX_SHAPE_DEPTH: Final = 3.0  # units the sample's shape is read from: mean hits >= this
MIX_ITERATIONS: Final = 50  # p99 change in coverage_mix vs 150: 0.05 log2


def _ztnb_logpmf(h: np.ndarray, mu: np.ndarray, v: float) -> np.ndarray:
    """log P(h | h > 0) of a negative binomial of mean ``mu`` and squared CV ``v``
    (Poisson at v = 0), without the log h! term (constant in the fit)."""
    if v == 0:
        return np.asarray(h * np.log(mu) - mu - np.log(-np.expm1(-mu)))
    r = 1 / v
    p0 = r * np.log(r / (r + mu))
    return np.asarray(
        gammaln(h + r) - gammaln(r) + p0 + h * np.log(mu / (r + mu)) - np.log(-np.expm1(p0))
    )


def _ztnb_mean(mean: np.ndarray, v: float, iterations: int = 60) -> np.ndarray:
    """mu whose zero-truncated mean mu / (1 - P(0)) is ``mean``, by fixed point (to ~0 where
    ``mean`` <= 1)."""
    mean = np.asarray(mean, dtype=np.float64)
    mu = np.maximum(mean, 1e-3)
    for _ in range(iterations):
        p0 = np.exp(-mu) if v == 0 else (1 + v * mu) ** (-1 / v)
        mu = np.maximum(mean * (1 - p0), 1e-3)
    return np.asarray(mu)


def rate_mixture(kmers: pl.DataFrame, coverage: pl.DataFrame) -> pl.DataFrame:
    """Per unit, ``coverage`` (``unit``, ``coverage``) corrected for a strain mix:
    ``coverage_mix`` and ``mix_rates`` (step 35, open item 3).

    Strains at different depths in one unit give its k-mers several rates: one strain's own
    k-mers its depth, the k-mers they share the sum. The hits of the hit k-mers only the unit
    holds in ``kmers`` (``unit``, ``hash``, ``hits``) are fitted as a mixture of 1-3
    zero-truncated negative binomials, chosen by BIC; the top rate is the summed depth of the
    strains sharing its k-mers. ``coverage_mix`` = ``coverage`` x top rate / one-rate fit
    (``coverage`` itself with one rate). The shape (squared CV of uneven coverage along a
    gene) is the sample's: the median of each unit's best one-rate shape over units with
    ``MIN_DISPERSION_KMERS`` such k-mers and mean hits >= ``MIX_SHAPE_DEPTH``. A per-unit
    shape absorbs a second rate into a wide one; a Poisson mixture splits uneven coverage
    into rates (simulation: top rate +0.45 log2 on one strain at cv 0.3, depth 10)."""
    own = kmers.filter(pl.len().over("hash") == 1)
    ids = own.group_by("unit").agg(n=pl.len()).filter(pl.col("n") >= MIN_DISPERSION_KMERS)
    rows = own.join(ids.with_row_index("i"), on="unit").group_by("i", "hits").agg(c=pl.len())
    if ids.height == 0:
        return coverage.select("unit", coverage_mix="coverage", mix_rates=pl.lit(1, pl.UInt8))
    n = ids.height
    u = rows["i"].to_numpy().astype(np.int64)
    h = rows["hits"].to_numpy().astype(np.float64)
    c = rows["c"].to_numpy().astype(np.float64)
    cnt = np.bincount(u, c, n)
    mean = np.bincount(u, c * h, n) / cnt

    def one(v: float) -> tuple[np.ndarray, np.ndarray]:  # one rate: closed form
        mu = _ztnb_mean(mean, v)
        return mu, np.bincount(u, c * _ztnb_logpmf(h, mu[u], v), n)

    ll = np.stack([one(v)[1] for v in MIX_SHAPES], 1)
    deep = mean >= MIX_SHAPE_DEPTH
    # ponytail: one shape per sample; per-gene-length shapes if short genes' ends dominate
    v = float(np.median(np.asarray(MIX_SHAPES)[ll[deep].argmax(1)])) if deep.any() else 0.0
    mu1, ll1 = one(v)
    best, top, rates = -2 * ll1 + np.log(cnt), mu1.copy(), np.ones(n, dtype=np.uint8)
    for k in (2, 3):
        mu = np.maximum(mean[:, None] * np.exp(np.linspace(-1, 1, k))[None], 0.5)
        w, ll = np.full((n, k), 1 / k), np.zeros(n)
        for _ in range(MIX_ITERATIONS):
            lp = np.log(w[u]) + _ztnb_logpmf(h[:, None], mu[u], v)
            mx = lp.max(1, keepdims=True)
            r = np.exp(lp - mx)
            s = r.sum(1, keepdims=True)
            r /= s
            ll = np.bincount(u, c * (mx + np.log(s))[:, 0], n)  # at the mixture before this M step
            rs = np.stack([np.bincount(u, c * r[:, j], n) for j in range(k)], 1)
            rh = np.stack([np.bincount(u, c * r[:, j] * h, n) for j in range(k)], 1)
            w = np.maximum(rs / cnt[:, None], 1e-12)
            mu = _ztnb_mean(rh / np.maximum(rs, 1e-12), v)
        bic = -2 * ll + (2 * k - 1) * np.log(cnt)
        better = bic < best
        best, top = np.where(better, bic, best), np.where(better, mu.max(1), top)
        rates[better] = k
    ids = ids.select("unit", ratio=pl.Series(top / mu1), mix_rates=pl.Series(rates))
    return (
        coverage.select("unit", "coverage")
        .join(ids, on="unit", how="left")
        .select(
            "unit",
            coverage_mix=pl.col("coverage") * pl.col("ratio").fill_null(1.0),
            mix_rates=pl.col("mix_rates").fill_null(1).cast(pl.UInt8),
        )
    )


def coverage_interval(
    coverage: np.ndarray,
    hit_kmers: np.ndarray,
    inflation: np.ndarray,
    level: float = 0.95,
) -> tuple[np.ndarray, np.ndarray]:
    """(lower, upper) of a ``level`` likelihood-ratio interval for zero-truncated Poisson
    coverage, per unit: ``hit_kmers`` k-mers hit (present x m x seen), their hits
    ``hit_kmers`` x E at the fitted ``coverage`` (the MLE), the log-likelihood
    H log c - K c - K log(1 - e^-c) divided by ``inflation`` (dispersion x clumping, a
    quasi-likelihood). Where every hit k-mer was hit about once the likelihood is flat
    towards 0 and the lower bound is 0: depth below ~1 is not identifiable from hit counts
    (step 35). Bisection on log c, 60 steps each side."""
    c = np.maximum(np.asarray(coverage, dtype=np.float64), 1e-9)
    k = np.asarray(hit_kmers, dtype=np.float64)
    h = k * c / -np.expm1(-c)
    w = np.maximum(np.asarray(inflation, dtype=np.float64), 1.0)
    cut = float(chi2.ppf(level, 1)) / 2

    def drop(x: np.ndarray) -> np.ndarray:  # log-likelihood below the maximum, scaled
        def ll(y: np.ndarray) -> np.ndarray:
            return np.asarray(h * np.log(y) - k * y - k * np.log(-np.expm1(-y)))

        return np.asarray((ll(c) - ll(x)) / w)

    def bound(far: np.ndarray) -> np.ndarray:
        near, out = np.log(c), np.log(far)  # drop(near) = 0 <= cut < drop(out), or none
        reach = drop(far) > cut
        for _ in range(60):
            mid = (near + out) / 2
            beyond = drop(np.exp(mid)) > cut
            near, out = np.where(beyond, near, mid), np.where(beyond, mid, out)
        return np.where(reach, np.exp((near + out) / 2), far)

    lower = bound(np.full_like(c, 1e-6))
    upper = bound(np.maximum(c * 1e3, 1e3))
    return np.where(lower <= 1e-6 * 1.0001, 0.0, lower), upper


def aai_fit(
    coverage: np.ndarray,
    present: np.ndarray,
    m: np.ndarray,
    pin_sum: np.ndarray,
    n_kmers: np.ndarray,
    windows: np.ndarray,
    model: SurvivalModel,
    clumping: np.ndarray | float = 1.0,
    level: float = 0.95,
    selected: float = 1.0,
    n_members: np.ndarray | float = 1.0,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Containment AAI from a zero-inflated fit, with a ``level`` interval, no draws:
    (``aai``, ``aai_lo``, ``aai_hi``, h). Arrays are per unit.

    The hit k-mers h = present x m x seen, seen = 1 - exp(-coverage), of an average
    member's ``pin_sum`` kept windows are turned into identity by :func:`survival_quantiles`
    under the index's survival model (``model``, :mod:`~kmer_functional_profiler.survival`),
    given h >= ``selected``, for units of ``n_members``; t = m / ``n_kmers`` is the windows'
    sampling rate and ``windows`` an average member's windows (length - k + 1). ``clumping``
    is 1 + the unit's expected kept k-mers per covering read: at low depth the sample size
    is the reads, not the k-mers (without it a strain at 100% identity and depth ~1 got
    intervals around 0.88-0.95, phase 7, step 8). The coverage is a zero-truncated Poisson MLE with
    variance lambda (1 - e^-l)^2 / (h (1 - e^-l - l e^-l)) (delta method, inflated by
    ``clumping``); the bounds take seen at the ends of its ``level`` interval.
    """
    lam = np.maximum(np.asarray(coverage, dtype=np.float64), 1e-9)
    m = np.asarray(m, dtype=np.float64)
    pin = np.maximum(np.asarray(pin_sum, dtype=np.float64), 1e-9)
    seen = -np.expm1(-lam)
    h = np.asarray(present, dtype=np.float64) * m * seen
    keep = np.clip(_div(m, np.asarray(n_kmers, dtype=np.float64)), 0.0, 1.0)
    clump = np.maximum(np.asarray(clumping, dtype=np.float64), 1.0)
    info = _div(h * (seen - lam * np.exp(-lam)), lam * seen**2)  # Fisher information on lambda
    spread = norm.ppf((1 + level) / 2) * np.sqrt(clump * _div(np.ones_like(info), info))
    seen_range = (-np.expm1(-np.maximum(lam - spread, 1e-9)), -np.expm1(-(lam + spread)))
    point, lo, hi = survival_quantiles(
        h, pin, seen, keep, windows, model, clump, level, max(selected, 1.0), seen_range,
        n_members,
    )  # fmt: skip
    return point, lo, hi, h


def aai_columns(
    coverage: np.ndarray,
    present: np.ndarray,
    m: np.ndarray,
    pin_sum: np.ndarray,
    n_kmers: np.ndarray,
    windows: np.ndarray,
    model: SurvivalModel,
    clumping: np.ndarray | float,
    min_kmers: float,
    n_members: np.ndarray | float = 1.0,
) -> pl.DataFrame:
    """:func:`aai_fit` as ``aai``, ``aai_lo``, ``aai_hi`` and ``aai_kmers``; the first three
    null below ``min_kmers`` hit k-mers (an estimate on a few k-mers is mostly the luck that
    put the unit in the profile) or ``MIN_AAI_LAMBDA`` coverage (where every hit k-mer has
    about one hit and the zero-inflated fit cannot tell coverage from presence), and
    estimated given h >= ``min_kmers``, the selection they are reported under."""
    *estimate, h = aai_fit(
        coverage, present, m, pin_sum, n_kmers, windows, model, clumping, selected=min_kmers,
        n_members=n_members,
    )  # fmt: skip
    # h is EM-fitted: the tolerance keeps a unit at exactly min_kmers hit k-mers
    keep = (h >= min_kmers - 1e-3) & (np.asarray(coverage) >= MIN_AAI_LAMBDA)
    names = ("aai", "aai_lo", "aai_hi")
    columns = {c: np.where(keep, x, np.nan) for c, x in zip(names, estimate, strict=True)}
    return pl.DataFrame(columns).fill_nan(None).with_columns(aai_kmers=pl.Series(h))


def unit_windows(units: pl.DataFrame, k: int) -> np.ndarray:
    """An average member's k-mer windows, ``len_mean`` - k + 1 (``n_kmers`` for indexes
    without ``len_mean``, e.g. sourmash imports, then ``m_g``)."""
    for column, shift in (("len_mean", 1 - k), ("n_kmers", 0), ("m_g", 0)):
        if column in units.columns:
            return np.asarray(np.maximum(units[column].to_numpy().astype(np.float64) + shift, 1.0))
    raise ValueError("units have no length column")


def unit_members(units: pl.DataFrame) -> np.ndarray:
    """Each unit's members (``n_members``; 1 for indexes without it, e.g. sourmash
    imports): the union the survival model's union term scales with."""
    if "n_members" not in units.columns:
        return np.ones(units.height)
    return np.asarray(np.maximum(units["n_members"].to_numpy().astype(np.float64), 1.0))


def greedy_coverage(own: pl.DataFrame, kmers_hit: np.ndarray, m: np.ndarray) -> pl.DataFrame:
    """Coverage without EM, from gather's assignment (as sylph's, on its reassigned k-mers).

    Per unit, its free k-mers are its ``m`` kept ones but those hit and taken by an earlier
    unit (``kmers_hit`` - own); of those, ``own`` (``unit``, ``hits``) are hit. The
    zero-inflated Poisson MLE over them: the zero-truncated rate of hits per own k-mer, so
    missing k-mers (divergence) do not lower it, unless that rate leaves more k-mers hit than
    are present, then the Poisson rate with all present (hits over free k-mers), as at low
    coverage, where hit k-mers are seen once and the truncated rate is 0.
    """
    per_unit = own.group_by("unit").agg(hits=pl.col("hits").sum(), kmers=pl.len()).sort("unit")
    unit = per_unit["unit"].to_numpy()
    hits, kmers = per_unit["hits"].to_numpy(), per_unit["kmers"].to_numpy()
    free = np.maximum(m[unit] - (kmers_hit[unit] - kmers), kmers).astype(np.float64)
    lam = ztp_lambda(hits / kmers)
    all_present = kmers >= free * -np.expm1(-lam)
    lam = np.where(all_present, hits / free, lam)
    return per_unit.select("unit", coverage_greedy=pl.Series(lam))


def ztp_lambda(mean: np.ndarray, iterations: int = 50) -> np.ndarray:
    """Zero-truncated Poisson MLE of the rate from the mean count of the non-zero items:
    lambda / (1 - exp(-lambda)) = mean, by Newton from lambda = mean; 0 where mean <= 1."""
    mean = np.asarray(mean, dtype=np.float64)
    lam = mean.copy()
    fit = mean > 1
    for _ in range(iterations):
        e = np.exp(-lam[fit])
        seen = -np.expm1(-lam[fit])
        f = lam[fit] / seen - mean[fit]
        df = (seen - lam[fit] * e) / seen**2
        lam[fit] = np.maximum(lam[fit] - f / df, 1e-9)
    lam[~fit] = 0.0
    return lam


def aai_naive(
    hits: np.ndarray,
    kmers_hit: np.ndarray,
    pin_sum: np.ndarray,
    keep: np.ndarray,
    windows: np.ndarray,
    model: SurvivalModel,
    min_kmers: float = MIN_AAI_KMERS,
    clumping: np.ndarray | float = 1.0,
    n_members: np.ndarray | float = 1.0,
) -> pl.DataFrame:
    """Per unit, independently (shared k-mers count for every holder, as sylph's ``query``):
    the point of :func:`survival_quantiles` (given at least ``min_kmers`` hit) for the
    k-mers hit of ``pin_sum``, seen = 1 - exp(-lambda), lambda the zero-truncated Poisson
    MLE from hits per hit k-mer, reads hitting ``clumping`` windows each. Below
    ``MIN_AAI_LAMBDA`` the coverage correction is not
    identifiable: seen is taken as 1 and ``aai_naive_lower_bound`` set. Null below
    ``min_kmers`` k-mers hit, as ``aai``."""
    lam = ztp_lambda(_div(hits.astype(np.float64), kmers_hit.astype(np.float64)))
    low = lam < MIN_AAI_LAMBDA
    seen = np.where(low, 1.0, -np.expm1(-np.maximum(lam, MIN_AAI_LAMBDA)))
    # Units are independent: only those reported are estimated, point only, AAI_CHUNK at a
    # time (the bisection's arrays are per unit; step 40).
    point = np.full(len(hits), np.nan)
    ok = np.flatnonzero(kmers_hit >= min_kmers)
    for start in range(0, len(ok), AAI_CHUNK):
        at = ok[start : start + AAI_CHUNK]

        def on(x: Any, at: np.ndarray = at) -> Any:
            return x[at] if np.ndim(x) else x

        point[at], _, _ = survival_quantiles(
            on(kmers_hit).astype(np.float64), on(pin_sum).astype(np.float64), on(seen), on(keep),
            on(windows), model, on(clumping), selected=max(min_kmers, 1.0),
            n_members=on(n_members), bounds=False,
        )  # fmt: skip
    return pl.DataFrame({"aai_naive": point, "aai_naive_lower_bound": low}).with_columns(
        pl.col("aai_naive").fill_nan(None)
    )


CHUNK_BYTES: Final = 2**25  # units x grid arrays of the posterior are built in chunks of this


class _Summed:
    """Per-batch aggregates summed into one table, re-aggregated each time the batches
    waiting outgrow the running total, so memory stays about twice the total's.

    The table is split into ``parts`` by the first key (mod ``parts``), each merged on its
    own, so a merge's transient memory (polars' group-by over the concatenation) is that of
    one part, not the whole table."""

    min_rows = 1_000_000  # rows waiting before the first merge, over all parts

    def __init__(self, keys: list[str], aggs: dict[str, pl.Expr], parts: int = 16) -> None:
        self.keys, self.aggs, self.parts = keys, aggs, parts
        self.totals: list[pl.DataFrame | None] = [None] * parts
        self.pending: list[list[pl.DataFrame]] = [[] for _ in range(parts)]
        self.rows = [0] * parts

    def add(self, df: pl.DataFrame) -> None:
        self.empty = df.clear()
        split = df.with_columns(_part=pl.col(self.keys[0]) % self.parts).partition_by(
            "_part", as_dict=True, include_key=False
        )
        for key, part in split.items():
            p = int(key[0])
            self.pending[p].append(part)
            self.rows[p] += part.height
            total = self.totals[p]
            if self.rows[p] > max(
                self.min_rows // self.parts, 0 if total is None else total.height
            ):
                self.totals[p], self.pending[p], self.rows[p] = self._merged(p), [], 0

    def _merged(self, p: int) -> pl.DataFrame:
        total = self.totals[p]
        parts = [total, *self.pending[p]] if total is not None else self.pending[p]
        return pl.concat([self.empty, *parts]).group_by(self.keys).agg(**self.aggs)

    def total(self) -> pl.DataFrame:
        """The summed table. It releases the parts, so it is called once, at the end."""
        merged = [self._merged(p) for p in range(self.parts)]
        self.totals, self.pending = [None] * self.parts, [[] for _ in range(self.parts)]
        return pl.concat(merged, rechunk=False)


SENSE_CHANGE: Final = 0.75  # share of random base substitutions that change the amino acid


def census_counts(
    batch: dict[str, np.ndarray],
    hit_hashes: np.ndarray,
    census_max: int,
    promiscuous: np.ndarray | None = None,
) -> tuple[int, int]:
    """Census k-mers of a batch of sampled k-mers, and how many of them are in the index.

    A census k-mer has hash <= ``census_max`` (a rate below which the index holds every
    unit's k-mers) and comes from a mate with one hashed frame: if the mate is coding, that
    frame is the coding one, so off-frame k-mers, which never hit, stay out. Frames are
    counted from the batch's sampled k-mers, so a kept frame with none sampled is missed
    (``census_frame_miss`` in :func:`sample_summary`). Counted per occurrence, as reads are.
    ``hit_hashes`` holds the batch's hashes found in the index; ``promiscuous`` (sorted)
    known hashes with no posting, which count as found too.
    """
    mate = batch["read"].astype(np.uint64) * 2 + batch["mate"]
    frames = np.unique(mate * 8 + batch["frame"])
    mates, n_frames = np.unique(frames // 8, return_counts=True)
    take = (batch["hash"] <= np.uint64(census_max)) & np.isin(mate, mates[n_frames == 1])
    census = batch["hash"][take]
    known = np.isin(census, hit_hashes)
    if promiscuous is not None and len(promiscuous):
        at = np.searchsorted(promiscuous, census).clip(max=len(promiscuous) - 1)
        known |= promiscuous[at] == census
    return len(census), int(known.sum())


def sample_summary(
    k: int,
    reads: int,
    mates: int,
    base_stats: dict[str, float] | None,
    census: tuple[int, int] | None,
    t_max: float,
    alphabet: str = "protein",
) -> dict[str, float | int | None]:
    """How much of the sample the index knows (Phase 7, step 14), as sylph's
    ``--estimate-unknown`` does for genomes. ``error_thinning`` r = (1 - d)^s, with s the
    nucleotides a k-mer spans (3k translated, k for a nucleotide index) and d the per-base
    chance of losing a k-mer window: bases masked or not A/C/G/T, plus the mean Phred error
    probability, times ``SENSE_CHANGE`` when reads are translated. Errors clustered in read
    tails make r too low, and the containment too high. (A model-based
    ``explained_fraction``, summed over detected units, was dropped in step 40: a gene hit by
    several units counted once per unit, up to 0.59 of read bases against 0.41 truth on the
    fmh benchmark.)

    ``census_containment`` (model-free): the fraction of census k-mers (:func:`census_counts`)
    found in the index, over r and capped at 1: includes units below detection, but exact
    k-mers only hit at ~80% AAI or more. ``census_frame_miss`` is the chance that a kept frame
    of a mate has no sampled k-mer, so the mate passes as one-frame and its off-frame
    k-mers lower the census.

    ``base_stats`` is :attr:`FastxHits.base_stats` (None for sourmash-hashed queries);
    ``census`` the summed (census k-mers, hits), None when the indexes do not support one.
    Values that cannot be computed are None. ``census_unknown`` is 1 minus the containment.
    """
    out: dict[str, float | int | None] = dict.fromkeys(
        (
            "bases", "lost_bases", "expected_errors", "mean_read_length", "error_thinning",
            "census_kmers", "census_hits", "census_containment", "census_unknown",
            "census_frame_miss",
        )
    )  # fmt: skip
    out["reads"], out["mates"] = reads, mates
    r = None
    if base_stats is not None and base_stats["bases"] > 0 and reads > 0:
        bases = base_stats["bases"]
        mate_nt = bases / (reads * mates)
        nucleotide = molecule(alphabet) == "nucleotide"
        sense = 1.0 if nucleotide else SENSE_CHANGE
        lost = (base_stats["lost"] + sense * base_stats["expected_errors"]) / bases
        r = (1 - min(lost, 1.0)) ** span_nt(k, alphabet)
        # k-mer windows per frame; a nucleotide read is one frame, so no frame is missed
        windows = float("inf") if nucleotide else (mate_nt - 2) / 3 - k + 1
        out |= {
            "bases": int(bases),
            "lost_bases": int(base_stats["lost"]),
            "expected_errors": float(base_stats["expected_errors"]),
            "mean_read_length": mate_nt,
            "error_thinning": r,
            "census_frame_miss": float((1 - t_max) ** max(windows, 0)),
        }
    if census is not None:
        out["census_kmers"], out["census_hits"] = census
        if census[0] > 0 and r is not None and r > 0:
            contained = min(1.0, census[1] / census[0] / r)
            out["census_containment"] = contained
            out["census_unknown"] = 1 - contained
    return out


def profile(
    index: Index,
    r1: str | Path,
    r2: str | Path | None = None,
    *,
    genetic_code: int = 11,
    frames: str = "edges:30",  # plan, phase 7, step 31
    min_qual: int = 0,
    batch_reads: int = 100_000,
    draws: int = 0,
    kmers_out: str | Path | None = None,
    timer: Timer | None = None,
    all_estimators: bool = False,
    low_memory: bool = False,
    with_aai: bool = False,
    min_aai: float = 0.0,
    min_aai_kmers: float = MIN_AAI_KMERS,
    background_share: float = BACKGROUND_SHARE,
    min_identity: float = 0.0,
    extra: Sequence[Index] = (),
    mask: "Mask | None" = None,
    summary: dict[str, float | int | None] | None = None,
    em_start: int = 0,
    greedy: bool = False,
    link_stats: bool = False,
) -> pl.DataFrame:
    """Per-unit hits, distinct k-mers hit, reads hit, containment and mean coverage.

    ``containment`` is the fraction of the unit's kept k-mers seen at least once;
    ``coverage`` is hits per kept k-mer; ``kmers_unique`` and ``gather_rank`` come from
    :func:`gather` (0 and null for units explained away) and ``coverage_em`` from :func:`em`
    over the units gather keeps (0 for the rest; ``hits_em`` its hits after the split), with
    ``coverage_zi`` and ``present_zi`` from its zero-inflated form, and
    ``coverage_zip``/``present_zip`` from :func:`em_pin`.
    ``copies_zi`` is present k-mers over an average member's kept k-mers (``pin_sum``), the
    member-equivalents present, and ``abundance_zi`` = ``coverage_zi`` x ``copies_zi``: for a
    unit whose members come from many genomes (a KO), total depth over its gene copies.
    ``kmers_wta``/``coverage_wta`` and ``kmers_ufirst``/``coverage_ufirst`` are the k-mers
    :func:`assign_best` gives each unit and their hits per kept k-mer. Units without hits
    are omitted. ``coverage_zi`` is the shipped coverage, with a likelihood-ratio interval
    ``coverage_zi_lo``/``_hi`` (:func:`coverage_interval`, widened by read clumping and by
    ``coverage_zi_dispersion``, :func:`kmer_dispersion`; with ``draws`` > 0 the posterior's
    instead); ``coverage_em`` and ``hits_em`` stay. ``coverage_mix`` and ``mix_rates``
    (:func:`rate_mixture`, strains at several depths) are experimental, with
    ``all_estimators`` only.
    ``aai`` is fitted with ``with_aai``,
    ``draws``, ``summary`` or ``all_estimators``; the other estimators (``_zib``, ``_zip``,
    ``_wta``, ``_ufirst``) only with ``all_estimators`` (benchmarks).

    With a dense tier, the reads are streamed a second time at its rate and the EM
    estimates (``coverage_em``, ``_zi``, ``_zib``) are fitted on the dense hits of the units
    gather keeps, over their ``m_dense`` k-mers; ``kmers_dense`` counts those hit.

    ``present_prob`` (:func:`presence`) is the probability that a unit gather keeps is
    present rather than hit by background, from the tier-1 k-mers gather gave it, with
    background at most ``background_share`` of those units; units gather drops get 0.
    ``present_q`` is its q-value (:func:`fdr_q`, over every row before ``min_aai`` drops
    any): keep ``present_q`` <= x for an expected false share of at most x. With
    ``min_identity`` (theta) > 0, implying ``with_aai``, "present" means a member within
    theta identity is in the sample: units whose ``aai_hi`` (the 95% upper bound, calibrated
    where the index is) is below theta are relatives, not the unit, and get
    ``present_prob`` 0; units without an ``aai`` (fewer than ``min_aai_kmers``) are not
    tested.

    ``draws`` > 0 adds 95% posterior intervals for ``coverage_zi`` and ``abundance_zi`` and
    ambiguity groups from that many Gibbs sweeps (:func:`posterior_zi`). Its per-read rows
    are kept from the first pass over the reads, or with ``low_memory`` rebuilt by a second
    pass for the detected units' k-mers only (same result; a read pass more, the rows of
    all hit units never held). ``kmers_out`` writes
    the tier-1 hits per (``unit``, ``hash``) with ``hits`` and ``holders`` to Parquet.

    Without ``with_aai`` (or ``all_estimators``, ``min_aai`` > 0) no AAI is estimated and
    the profile holds the units gather keeps only (step 40); with it, every hit unit.
    Containment AAI (sylph's containment ANI in protein space): ``aai_naive`` on every hit
    unit from its own tier-1 hits (:func:`aai_naive`), and ``aai`` on the units gather keeps
    (with the zero-inflated fit), from the hit k-mers after the EM split
    (``aai_kmers``), with a closed-form 95% interval ``aai_lo``/``_hi`` (:func:`aai_fit`; no
    draws needed). Both invert the survival model of the index's ``aai_model.json``
    (:func:`survival`; a^k without one) and are null below ``min_aai_kmers`` hit k-mers
    (``aai_kmers``, or ``kmers_hit`` for ``aai_naive``, still reported). An index with an
    ``aai_calibration.json`` maps them to alignment identity (:func:`calibrate_aai`) and
    keeps the raw ones as ``aai_raw``, ``aai_raw_lo``/``_hi``.
    ``component`` labels the hit units linked by shared k-mers by their smallest unit id in
    the profile, so
    the units that bracket a sample variant can be read together. ``min_aai`` > 0 drops rows
    with ``aai_naive`` below it or null.

    ``extra`` indexes (same k, alphabet and hash scheme) are queried jointly with ``index``:
    their units compete with its units in gather, EM and the posterior, with ids offset by
    the units of the indexes before, and a ``source`` column (0 = ``index``) in the output.

    ``mask`` (:class:`~kmer_functional_profiler.mask.Mask`, built against ``index``) drops
    masked sampled hashes before any lookup and subtracts the masked postings from the unit
    rows' expected counts; ``masked_fraction`` and ``host_like`` are added.

    ``summary``, if given, is filled with :func:`sample_summary`'s sample-level stats and
    known-k-mer census. The census needs every index to have a complete
    base stratum (``base_stratum_complete`` in ``meta.json``) and the kfp hash; an index
    built without ``promiscuous.npy`` counts its promiscuous k-mers as unknown.

    The index's ``alphabet`` sets how reads are hashed: amino-acid alphabets translate them
    (``frames`` and ``genetic_code`` apply), ``dna`` hashes them untranslated as canonical
    k-mers (``frames`` and ``genetic_code`` are ignored), and the ``aai`` columns are named
    ``ani`` (``min_aai`` then filters on ``ani_naive``).

    ``greedy`` replaces the EM fits by gather's assignment (:func:`greedy_coverage`, as
    sylph; ``coverage_greedy``): no ``_em``/``_zi`` columns and no ``aai``, so it excludes
    ``draws`` and ``all_estimators``; with ``with_aai``, ``aai_naive`` is reported.

    ``link_stats`` (with ``timer``) adds the largest component's link and promiscuity cuts to
    the counts (:func:`link_cuts`, ablations only: GBs at depth).

    ``em_start`` > 0 seeds random starting points for the EM fits (:func:`em`,
    :func:`em_pin`; not the posterior's): profiles from several starts agree where the fits
    converged, a check that needs no truth (phase 7, step 30).

    ``timer`` records each stage's time and peak RSS and the counts the query's cost
    hinges on (sampled k-mers, distinct ones estimated on 1 in ``DISTINCT_SAMPLE`` of hash
    space, hit k-mers, hit rows, (unit, hash) pairs, detected units, component sizes).
    """
    if greedy and (draws > 0 or all_estimators):
        raise ValueError("greedy has no zero-inflated fit: no draws or other estimators")
    # aai (aai_naive, and aai from the zero-inflated fit) only when asked for (step 40);
    # without it the profile holds the detected units only.
    aai = with_aai or all_estimators or min_aai > 0 or min_identity > 0
    record = timer is not None
    timer = timer or Timer()
    counts = timer.counts
    params = _index_params(index.meta)
    # ponytail: the first index's survival model for every unit; extra indexes share its k
    # and alphabet, which the model is fitted for
    model = SurvivalModel.from_json(index.aai_model, params.k)
    joint = _Joint([index, *extra])
    sources: list[_core.FastxHits] = []  # the read streams, for their base counts
    census_max = (
        min(
            min(_core.max_hash(base_rate(_index_params(i.meta))), int(i.tier1.max_hash))
            for i in joint.indexes
        )
        if summary is not None
        and index.meta.get("hash") != "sourmash"
        and all(i.meta.get("base_stratum_complete") for i in joint.indexes)
        else None
    )
    census = [0, 0]
    # Promiscuous base-stratum hashes are known but never hit (index.promiscuous_base).
    promiscuous = (
        None
        if census_max is None
        else np.unique(
            np.concatenate(
                [np.empty(0, np.uint64)]
                + [i.promiscuous for i in joint.indexes if i.promiscuous is not None]
            )
        )
    )

    def stream(max_hash: int) -> Iterable[dict[str, np.ndarray]]:
        batches = unmasked(max_hash)
        if mask is None:
            return batches
        return (
            {c: v[keep] for c, v in b.items()} for b in batches for keep in [mask.keep(b["hash"])]
        )

    def unmasked(max_hash: int) -> Iterable[dict[str, np.ndarray]]:
        if index.meta.get("hash") == "sourmash":
            return sourmash_hits(r1, r2, params.k, max_hash, batch_reads)
        source = _core.FastxHits(
            r1,
            r2,
            k=params.k,
            alphabet=params.alphabet,
            genetic_code=genetic_code,
            frames=frames,
            max_hash=max_hash,
            batch_reads=batch_reads,
            min_qual=min_qual,
        )
        sources.append(source)
        return source

    def by_read(hits: pl.DataFrame) -> pl.DataFrame:
        return hits.group_by("unit", "hash", "read").agg(
            n=pl.len(), pin_q=pl.col("pin_q").first(), holders=pl.col("holders").first()
        )

    def reread(kmers: np.ndarray) -> pl.DataFrame:
        """(hash, read, n) rows for the sampled ``kmers``: times each read has each, which is
        also a hit's rows over its holders, so these are the first pass's rows for them."""
        wanted = np.sort(kmers).astype(np.uint64)
        parts = [pl.DataFrame(schema={"hash": pl.UInt64, "read": pl.UInt64, "n": pl.UInt32})]
        for b in stream(joint.max_hash("tier1")):
            at = np.minimum(np.searchsorted(wanted, b["hash"]), max(len(wanted) - 1, 0))
            keep = wanted[at] == b["hash"] if len(wanted) else np.zeros(len(at), dtype=bool)
            found = pl.DataFrame({"hash": b["hash"][keep], "read": b["read"][keep]})
            parts.append(found.group_by("hash", "read").agg(n=pl.len().cast(pl.UInt32)))
        return pl.concat(parts)

    def per_kmer(per_read: pl.DataFrame) -> pl.DataFrame:
        return per_read.group_by("unit", "hash").agg(
            hits=pl.col("n").sum(), pin_q=pl.col("pin_q").first(), holders=pl.col("holders").first()
        )

    # Per-(unit, hash) hits and per-unit reads are summed batch by batch (read ids never
    # span batches), so per-read rows are kept only for the posterior.
    first = {"pin_q": pl.col("pin_q").first(), "holders": pl.col("holders").first()}
    pairs = _Summed(["unit", "hash"], {"hits": pl.col("hits").sum(), **first})
    reads_summed = _Summed(["unit"], {"reads": pl.col("reads").sum()})
    empty = np.empty(0, dtype=np.uint64)
    batches, n_reads, subsample = [], 0, []
    # The posterior's per-read rows come from this pass, unless a dense tier's pass gives
    # them or ``low_memory`` re-reads the reads for the detected k-mers only.
    keep_reads = draws > 0 and not joint.dense and not low_memory
    sampled = hit_kmers = hit_rows = 0
    per_pair_hist = np.zeros(1, dtype=np.int64)  # pairs by hit k-mers (phase 13, step 0)
    # The sample's distinct sampled hashes first, a read pass more, so that tier 1 is read
    # once in hash order and dropped (plan, phase 6, step 40 (c)), not paged in at random.
    # Kept sorted in hash ranges in Rust (~1.25x the distinct hashes; step 40).
    distinct = _core.DistinctHashes(joint.max_hash("tier1"))
    with timer("sketch"):
        for batch in stream(joint.max_hash("tier1")):
            distinct.add(np.ascontiguousarray(batch["hash"], dtype=np.uint64))
    with timer("walk"):
        joint.walk("tier1", (distinct.take(i) for i in range(len(distinct))))
    del distinct
    reads = iter(stream(joint.max_hash("tier1")))
    b: dict[str, np.ndarray] | None = {"hash": empty, "read": empty}
    while b is not None:
        with timer("lookup"):
            hits = joint.hits("tier1", b["hash"], b["read"])
        with timer("aggregate"):
            pairs.add(hits.group_by("unit", "hash").agg(hits=pl.len(), **first))
            reads_summed.add(
                hits.select("unit", "read").unique().group_by("unit").agg(reads=pl.len())
            )
            if keep_reads:
                # Once per (hash, read): each hit has a row per unit holding the k-mer.
                batches.append(
                    hits.group_by("hash", "read").agg(n=pl.len() // pl.col("holders").first())
                )
        if len(b["read"]):  # reads are numbered in input order; the last has sampled hashes
            n_reads = max(n_reads, int(b["read"].max()) + 1)
            if census_max is not None:
                got = census_counts(b, hits["hash"].unique().to_numpy(), census_max, promiscuous)
                census = [census[0] + got[0], census[1] + got[1]]
        sampled += len(b["hash"])
        subsample.append(b["hash"][b["hash"] <= joint.max_hash("tier1") // DISTINCT_SAMPLE])
        hit_rows += hits.height
        hit_kmers += round((1 / hits["holders"]).sum()) if hits.height else 0  # rows per hit
        if record and hits.height:  # occurrences per pair, holders collapsed
            per_pair = hits.group_by("read").agg((1 / pl.col("holders")).sum().round())
            pair_counts = np.bincount(per_pair[:, 1].cast(pl.Int64).to_numpy())
            grow = max(0, len(pair_counts) - len(per_pair_hist))
            per_pair_hist = np.pad(per_pair_hist, (0, grow))
            per_pair_hist[: len(pair_counts)] += pair_counts
        with timer("hash"):
            b = next(reads, None)
    joint.walked.pop("tier1", None)  # the walk's table: every batch is looked up
    with timer("aggregate"):
        per_read = pl.concat(batches) if batches else None
        kmer_hits, unit_reads = pairs.total(), reads_summed.total()
    counts |= {
        "reads": n_reads,
        "sampled_kmers": sampled,
        "distinct_sampled_kmers_est": len(np.unique(np.concatenate(subsample))) * DISTINCT_SAMPLE,
        "hit_kmers": hit_kmers,
        "hit_rows": hit_rows,
        "read_rows": 0 if per_read is None else per_read.height,
        "unit_kmer_pairs": kmer_hits.height,
        "hit_units": kmer_hits["unit"].n_unique(),
    }
    if record:  # index i: pairs with i hit k-mer occurrences; 0 = pairs with no hit
        per_pair_hist[0] = n_reads - per_pair_hist[1:].sum()
        counts["hits_per_pair"] = per_pair_hist.tolist()
    with timer("hit_units"):  # rows gathered from the unit columns, never the whole table
        hit_info = joint.rows(np.sort(kmer_hits["unit"].unique().to_numpy()))
        if mask is not None:
            hit_info = mask.adjust(hit_info)
    if kmers_out is not None:
        kmer_hits.write_parquet(kmers_out)
    # From here units are numbered 0.. in the hit units' order (``index`` keeps the index's
    # ids), so per-unit arrays come from their rows only: a full column is ~4 GB per 8 bytes.
    ids = hit_info.select(
        pl.col("unit").cast(kmer_hits.schema["unit"]),
        index=pl.int_range(pl.len(), dtype=pl.UInt32),
    )

    def renumber(df: pl.DataFrame) -> pl.DataFrame:
        return df.join(ids, on="unit").drop("unit").rename({"index": "unit"})

    kmer_hits, unit_reads = map(renumber, (kmer_hits, unit_reads))
    hit_info = hit_info.with_columns(index=pl.col("unit"), unit=ids["index"])
    with timer("components"):  # a union-find on the local ids (step 40)
        pair_unit = kmer_hits["unit"].cast(pl.UInt32).to_numpy()
        label = _core.components(pair_unit, kmer_hits["hash"].to_numpy(), hit_info.height)
        component = pl.DataFrame(
            {"unit": np.arange(len(label), dtype=np.uint32), "component": label}
        ).with_columns(pl.col("unit").cast(kmer_hits.schema["unit"]))
        if record and kmer_hits.height:  # a sample with no hit has no component
            sizes = np.bincount(label)
            largest = int(np.argmax(sizes))  # the first of the largest
            counts["components"] = len(sizes)
            counts["largest_component_units"] = int(sizes[largest])
            counts["largest_component_pairs"] = int(np.count_nonzero(label[pair_unit] == largest))
        del pair_unit
    t_g = hit_info["t_g"].to_numpy()
    with timer("gather"):
        assigned, owned = gather(kmer_hits.select("unit", "hash"), t_g)
        # Each hit k-mer is gather's: the first unit in gather order holding it.
        own = kmer_hits.filter(pl.Series(owned)).select("unit", "hash", "holders", "hits")
        del owned
    counts["detected_units"] = assigned.height
    detected = kmer_hits.join(assigned.select("unit"), on="unit")
    m_g = hit_info["m_g"].to_numpy()
    pin_hist = hit_info["pin_hist"].to_numpy()
    pin_sum = hit_info["pin_sum"]
    len_cv = hit_info["len_cv"].to_numpy()
    with timer("presence"):
        present_prob = presence(
            own, t_g, n_reads, joint.height, share=background_share, report=counts
        )
    greedy_fit = (
        greedy_coverage(own, np.bincount(kmer_hits["unit"].to_numpy(), minlength=len(m_g)), m_g)
        if greedy
        else None
    )
    del own
    if not aai:  # the undetected units are dropped here, units renumbered 0.. over the rest
        keep = assigned.select("unit").sort("unit")
        new = keep.with_columns(new=pl.int_range(pl.len(), dtype=pl.UInt32))

        def detected_only(df: pl.DataFrame) -> pl.DataFrame:
            return df.join(new, on="unit").drop("unit").rename({"new": "unit"})

        kmer_hits, unit_reads, component, detected, assigned, present_prob = map(
            detected_only, (kmer_hits, unit_reads, component, detected, assigned, present_prob)
        )
        if greedy_fit is not None:
            greedy_fit = detected_only(greedy_fit)
        hit_info = hit_info[keep["unit"].to_numpy()].with_columns(unit=new["new"])
        ids = ids.join(new.rename({"unit": "index"}), on="index").select("unit", index="new")
        t_g, m_g = hit_info["t_g"].to_numpy(), hit_info["m_g"].to_numpy()
        pin_hist, pin_sum = hit_info["pin_hist"].to_numpy(), hit_info["pin_sum"]
        len_cv = hit_info["len_cv"].to_numpy()
        del keep, new
    if joint.dense:
        # Second pass: every k-mer at the dense rate, for the detected units only.
        keep = assigned.select("unit")  # detected units are hit units, so renumber keeps all
        with timer("dense"):
            detected_reads = pl.concat(
                [
                    by_read(
                        renumber(joint.hits("dense", b["hash"], b["read"])).join(
                            keep, on="unit", how="semi"
                        )
                    )
                    for b in [{"hash": empty, "read": empty}, *stream(joint.max_hash("dense"))]
                ]
            )
            detected = per_kmer(detected_reads)
            if draws > 0:
                per_read = detected_reads.group_by("hash", "read").agg(pl.col("n").first())
        m_g = hit_info["m_dense"].to_numpy()
        pin_hist = hit_info["pin_hist_dense"].to_numpy()
        pin_sum = hit_info["pin_sum_dense"]
        len_cv = hit_info["len_cv_dense"].to_numpy()
    if link_stats:
        with timer("links"):
            counts |= link_cuts(detected)
    # Clumping: a read covering a unit hits ~ Poisson(mu) of its kept k-mers, which inflates
    # the variance of hit counts by 1 + mu. mu from tier-1 hits per hitting read (a
    # zero-truncated mean).
    ratio = (
        kmer_hits.group_by("unit")
        .agg(hits=pl.col("hits").sum())
        .join(unit_reads, on="unit")
        .select("unit", per_read=pl.col("hits") / pl.col("reads"))
    )
    mu_tier1 = np.zeros(len(m_g))
    mu_tier1[ratio["unit"].to_numpy()] = ztp_lambda(ratio["per_read"].to_numpy())
    if greedy:  # gather's assignment instead of the EM fits, as sylph (no fit_em, fit_zi)
        assert greedy_fit is not None
        fits = [greedy_fit]
    else:
        start = np.random.default_rng(em_start) if em_start else None
        with timer("fit_em"):
            # Fits stop per component, so batches of components give the same units' results.
            parts = component_batches(detected)
            detected = pl.concat(parts, rechunk=False)  # the batches' rows, not a second copy

            def per_batch(fit: Callable[[pl.DataFrame], pl.DataFrame]) -> pl.DataFrame:
                return pl.concat([fit(part) for part in parts]).sort("unit")

            counts["fit_batches"] = len(parts)
            counts["fit_largest_batch_pairs"] = max(part.height for part in parts)
            plain = per_batch(lambda part: em(part, m_g, report=counts, start=start)).select(
                "unit", coverage_em="coverage"
            )
            # Hits after the EM split (what the function x taxon table splits): coverage x m.
            plain = plain.with_columns(hits_em=plain["coverage_em"] * m_g[plain["unit"].to_numpy()])
        fits = [plain]
        # coverage_zi is the shipped coverage (step 35): the zero-inflated fit always runs
        with timer("fit_zi"):
            zi_report: dict[str, float] = {}  # its em_* (step 40: 9,059 s block-wise, unreported)
            warm = None
            if ZI_WARM_START:
                warm = np.zeros(len(m_g))
                warm[plain["unit"].to_numpy()] = plain["coverage_em"].to_numpy()
            inflated = per_batch(
                lambda part: em(
                    part, m_g, zero_inflated=True, report=zi_report, start=start, init=warm
                )
            )
            counts |= {f"zi_{key}": value for key, value in zi_report.items()}
        # Present k-mers over an average member's kept k-mers: member-equivalents present.
        copies = pl.col("present") * pl.col("m") / pl.col("pin_sum")
        zi = inflated.join(
            pl.DataFrame(
                {
                    "unit": np.arange(len(m_g), dtype=np.uint32),
                    "m": m_g,
                    "pin_sum": pin_sum,
                    # sourmash imports have no n_kmers: windows taken as all kept (t = 1)
                    "n_kmers": hit_info["n_kmers"] if "n_kmers" in hit_info.columns else m_g,
                    "windows": unit_windows(hit_info, params.k),
                    "n_members": unit_members(hit_info),
                }
            ),
            on="unit",
        )
        # scaled to the fitted tier's k-mers
        mu = mu_tier1 * _div(m_g.astype(np.float64), hit_info["m_g"].to_numpy().astype(np.float64))
        estimate = aai_columns(
            zi["coverage"].to_numpy(),
            zi["present"].to_numpy(),
            zi["m"].to_numpy(),
            zi["pin_sum"].to_numpy(),
            zi["n_kmers"].to_numpy(),
            zi["windows"].to_numpy(),
            model,
            1 + mu[zi["unit"].to_numpy()],
            min_aai_kmers,
            zi["n_members"].to_numpy(),
        )
        fitted_zi = zi.select(
            "unit",
            coverage_zi="coverage",
            present_zi="present",
            copies_zi=copies,
            abundance_zi=pl.col("coverage") * copies,
        )
        if aai:
            fitted_zi = fitted_zi.hstack(estimate)
        unit_kmers = detected.select("unit", "hash", "hits")
        if all_estimators:  # experimental: rarely splits real units (step 35, HPC results)
            fitted_zi = fitted_zi.join(rate_mixture(unit_kmers, zi), on="unit", how="left")
        if draws == 0:  # with draws, the posterior gives coverage_zi's interval
            spread = kmer_dispersion(unit_kmers, zi)
            zi = zi.join(spread, on="unit", how="left")
            lam = zi["coverage"].to_numpy()
            hit_kmers = zi["present"].to_numpy() * zi["m"].to_numpy() * -np.expm1(-lam)
            inflation = zi["dispersion"].fill_null(1.0).to_numpy() * (1 + mu[zi["unit"].to_numpy()])
            lo, hi = coverage_interval(lam, hit_kmers, inflation)
            fitted_zi = fitted_zi.with_columns(
                coverage_zi_lo=pl.Series(np.where(lam > 0, lo, 0.0)),
                coverage_zi_hi=pl.Series(np.where(lam > 0, hi, 0.0)),
                coverage_zi_dispersion=zi["dispersion"].fill_null(1.0),
            )
        fits.append(fitted_zi)
        if all_estimators:
            with timer("fit_zib"):
                prior = fit_present_prior(inflated)
                fitted = (
                    per_batch(
                        lambda part: em(part, m_g, zero_inflated=True, prior=prior, start=start)
                    )
                    if prior
                    else inflated
                )
                fits.append(fitted.select("unit", coverage_zib="coverage", present_zib="present"))
            with timer("fit_zip"):
                fits.append(
                    per_batch(lambda part: em_pin(part, pin_hist, start=start)).select(
                        "unit", coverage_zip="coverage", present_zip="present"
                    )
                )
    naive = None
    if aai:  # aai_naive of every hit unit
        with timer("aai"):
            per_unit = (
                kmer_hits.group_by("unit")
                .agg(hits=pl.col("hits").sum(), kmers=pl.len())
                .sort("unit")
            )
            at = per_unit["unit"].to_numpy()
            tier1 = hit_info["m_g"].to_numpy().astype(np.float64)
            naive = aai_naive(
                per_unit["hits"].to_numpy(),
                per_unit["kmers"].to_numpy(),
                hit_info["pin_sum"].to_numpy()[at],
                np.clip(_div(tier1, hit_info["n_kmers"].to_numpy().astype(np.float64)), 0, 1)[at]
                if "n_kmers" in hit_info.columns
                else np.ones(len(at)),
                unit_windows(hit_info, params.k)[at],
                model,
                min_aai_kmers,
                1 + mu_tier1[at],
                unit_members(hit_info)[at],
            ).with_columns(unit=per_unit["unit"])
    unit_info = hit_info.select(pl.exclude("^(pin_(hist|sum)|len_cv).*$"))
    if "name" not in unit_info.columns:  # sourmash imports name units, builds by cluster_rep
        unit_info = unit_info.with_columns(name=pl.col("cluster_rep").cast(pl.String))
    with timer("result"):
        result = (
            kmer_hits.group_by("unit")
            .agg(hits=pl.col("hits").sum().cast(pl.UInt64), kmers_hit=pl.len().cast(pl.UInt32))
            .join(unit_reads.with_columns(pl.col("reads").cast(pl.UInt64)), on="unit")
            .join(unit_info, on="unit")
            .join(assigned, on="unit", how="left")
            .join(present_prob, on="unit", how="left")
        )
        if naive is not None:
            result = result.join(naive, on="unit")
        result = result.join(component, on="unit")
        for fit in fits:
            result = result.join(fit, on="unit", how="left")
    if detected.height == 0:  # nothing detected (e.g. host-only reads): no posterior to draw
        draws = 0
    if draws > 0 and per_read is None:
        with timer("reread"):
            per_read = reread(detected["hash"].unique().to_numpy())
    if draws > 0:
        prob = np.zeros(len(m_g))
        prob[present_prob["unit"].to_numpy()] = present_prob["present_prob"].to_numpy()
        assert per_read is not None  # per-read rows are kept when draws > 0
        with timer("posterior"):
            intervals = posterior_zi(
                detected.select("unit", "hash"),
                per_read,
                m_g,
                pin_sum.to_numpy(),
                draws,
                present_prob=prob,
                len_cv=len_cv,
            )
        result = result.join(intervals, on="unit", how="left")
    if joint.dense:
        result = result.join(
            detected.group_by("unit").agg(kmers_dense=pl.len().cast(pl.UInt32)),
            on="unit",
            how="left",
        ).with_columns(pl.col("kmers_dense").fill_null(0))
    if all_estimators:
        rated = kmer_hits.join(hit_info.select("unit", "m_g", "t_g"), on="unit")
        with timer("baselines"):
            for rule, score in (("wta", WTA_SCORE), ("ufirst", UFIRST_SCORE)):
                won = (
                    assign_best(rated, score)
                    .group_by("unit")
                    .agg(
                        pl.len().cast(pl.UInt32).alias(f"kmers_{rule}"),
                        pl.col("hits").sum().alias(f"_{rule}"),
                    )
                )
                result = result.join(won, on="unit", how="left").with_columns(
                    pl.col(f"kmers_{rule}").fill_null(0),
                    (pl.col(f"_{rule}").fill_null(0) / pl.col("m_g")).alias(f"coverage_{rule}"),
                )
            result = result.drop("_wta", "_ufirst")
    # The hit tables are done with: freed before the result's last copies (step 40).
    del kmer_hits, detected, unit_reads, component, hit_info, unit_info, assigned, naive
    del present_prob, fits
    with timer("finish"):
        result = _finish(result, joint, min_aai, params.alphabet, min_identity)
    if summary is not None:  # over every unit, before min_aai drops rows
        # The first stream is the full first pass over the reads.
        summary |= sample_summary(
            params.k,
            n_reads if not sources else sources[0].n_reads,
            2 if r2 is not None else 1,
            sources[0].base_stats if sources else None,
            None if census_max is None else (census[0], census[1]),
            joint.max_hash("tier1") / 2**64,
            params.alphabet,
        )
    return result


def _finish(
    result: pl.DataFrame, joint: _Joint, min_aai: float, alphabet: str, min_identity: float = 0.0
) -> pl.DataFrame:
    """The profile's last columns, its rows in unit order, aai calibrated, presence tested
    against ``min_identity`` and given q-values, rows filtered by ``min_aai``."""
    result = (
        result.with_columns(
            pl.col(
                "^(coverage|present|copies|abundance|hits)_(em|zi|zib|zip|greedy)(_lo|_hi)?$"
            ).fill_null(0.0),
            pl.col("present_prob").fill_null(0.0),
            pl.col("^aai_kmers$").fill_null(0.0),
            containment=pl.col("kmers_hit") / pl.col("m_g"),
            coverage=pl.col("hits") / pl.col("m_g"),
            kmers_unique=pl.col("kmers_unique").fill_null(0),
        )
        .with_columns(unit=pl.col("index"))  # back to the index's unit ids
        .with_columns(component=pl.col("unit").min().over("component"))
        .drop("index")
        .sort("unit")
    )
    result = _calibrated(result, [i.aai_calibration for i in joint.indexes])
    if min_identity > 0:  # confidently below theta: a relative's sequence, not the unit
        result = result.with_columns(
            present_prob=pl.when(pl.col("aai_hi") < min_identity)
            .then(0.0)
            .otherwise(pl.col("present_prob"))
        )
    result = result.with_columns(present_q=pl.Series(fdr_q(result["present_prob"].to_numpy())))
    if min_aai > 0:
        result = result.filter(pl.col("aai_naive") >= min_aai)
    if molecule(alphabet) == "nucleotide":  # identity of nucleotides: ANI, as sylph's
        result = result.rename({c: "ani" + c[3:] for c in result.columns if c.startswith("aai")})
    return result
