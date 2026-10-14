"""The shared execution path behind ``anqchor match|cis|trans run``, ``anqchor run``,
``anqchor report`` and ``anqchor steps``."""

from __future__ import annotations

import contextlib
import copy
import errno
import fcntl
import os
import sys
from collections.abc import Iterator
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import typer

from anqchor.config.load import (
    ConfigError,
    config_to_dict,
    load_config,
    read_yaml,
    write_resolved_config,
)
from anqchor.config.model import Config
from anqchor.engine.minimake import changed_signatures, log_path, write_pipeline_makefile
from anqchor.engine.provenance import now_iso, write_run_json
from anqchor.engine.runner import StepStatus, invalidate_missing_outputs, run_make, step_status
from anqchor.engine.step import StepGraph
from anqchor.pipelines import PIPELINES, BuildContext, Layout, build_graph
from anqchor.qpgen.commands import THREADS_ENV
from anqchor.qpgen.locate import QpgenNotFound, locate
from anqchor.util.files import sha256_file
from anqchor.util.log import get_logger

log = get_logger("cli")

EXIT_CONFIG = 2
EXIT_STEP_FAILED = 3
EXIT_QPGEN = 4


@dataclass
class RunOptions:
    """Options common to every command that executes pipelines."""

    config: Path | None = None
    out: Path | None = None
    jobs: int = 1
    threads: int | None = None
    dataset: list[str] = field(default_factory=list)
    strategy: list[str] = field(default_factory=list)
    dry_run: bool = False
    until: str | None = None
    only: str | None = None
    force: bool = False
    keep_going: bool = False
    qpgen: Path | None = None
    plotlyjs: str | None = None
    overrides: dict[str, Any] = field(default_factory=dict)
    reports: bool = True
    """Include the HTML report steps (``all``) or stop at the summaries (``run``)."""


def fail(msg: str, code: int) -> typer.Exit:
    """Print ``error: msg`` in red to stderr and return the ``typer.Exit`` to raise."""
    typer.secho(f"error: {msg}", fg=typer.colors.RED, err=True)
    return typer.Exit(code)


def build_config(opts: RunOptions, pipelines: list[str]) -> tuple[Config, list[str]]:
    """Load the config file and flag overrides, then apply ``--dataset``/``--strategy``.

    Returns the configuration and the pipelines to run. ``--dataset`` and ``--strategy``
    restrict each pipeline to the named datasets/strategies it is configured with. When
    several pipelines run (``anqchor run``), a pipeline left with no dataset or no
    strategy is skipped (with a note); it is an error only when no pipeline is left. A
    single explicitly named pipeline (``anqchor cis run --dataset X``) still accepts a
    dataset outside its ``datasets`` list.
    """
    overrides = copy.deepcopy(opts.overrides)  # never mutate the caller's options
    if opts.out is not None:
        overrides["out_dir"] = str(opts.out)
    if opts.threads is not None:
        overrides.setdefault("qpgen", {})["threads"] = opts.threads
    if opts.qpgen is not None:
        overrides.setdefault("qpgen", {})["binary"] = str(opts.qpgen)
    if opts.plotlyjs is not None:
        overrides.setdefault("report", {})["plotlyjs"] = opts.plotlyjs
    if opts.config is None and not overrides.get("genotypes"):
        raise ConfigError("give --config FILE, or the input flags (see --help)")
    if opts.config is None:
        # Flag-only invocation: run just the requested pipelines.
        pl = overrides.setdefault("pipelines", {})
        for p in PIPELINES:
            pl.setdefault(p, {})["enabled"] = p in pipelines
    else:
        retarget_flag_names(overrides, read_yaml(Path(opts.config)))
    cfg = load_config(opts.config, overrides or None)
    if opts.dataset:
        unknown = [d for d in opts.dataset if d not in cfg.datasets]
        if unknown:
            raise ConfigError(f"--dataset: unknown dataset(s) {', '.join(unknown)}")
    single = len(pipelines) == 1
    active: list[str] = []
    for p in pipelines:
        pc = getattr(cfg.pipelines, p)
        if opts.dataset:
            in_pipeline = cfg.datasets_for(p)  # type: ignore[arg-type]
            keep_ds = [d for d in in_pipeline if d in opts.dataset]
            if not keep_ds and single:
                keep_ds = list(opts.dataset)  # the user named this pipeline explicitly
            if not keep_ds:
                log.info(
                    "--dataset %s: not in pipeline %s; skipping it", ", ".join(opts.dataset), p
                )
                continue
            pc.datasets = keep_ds
        if opts.strategy and p in ("cis", "trans"):
            keep = [s for s in pc.strategies if s.name in opts.strategy]
            if not keep:
                if single:
                    raise ConfigError(
                        f"--strategy: no {p} strategy named {', '.join(opts.strategy)}"
                    )
                log.info(
                    "--strategy %s: no %s strategy of that name; skipping %s",
                    ", ".join(opts.strategy),
                    p,
                    p,
                )
                continue
            pc.strategies = keep
        active.append(p)
    if not active:
        what = []
        if opts.dataset:
            what.append(f"--dataset {', '.join(opts.dataset)}")
        if opts.strategy:
            what.append(f"--strategy {', '.join(opts.strategy)}")
        raise ConfigError(
            f"{' and '.join(what)}: matches none of the pipelines {', '.join(pipelines)}"
        )
    if opts.strategy and not any(p in ("cis", "trans") for p in active):
        raise ConfigError(f"--strategy: no cis/trans strategy named {', '.join(opts.strategy)}")
    return cfg, active


def retarget_flag_names(overrides: dict[str, Any], raw: dict[str, Any]) -> None:
    """Apply input flags to the configuration file's dataset/reference when it has one.

    Flags such as ``--pheno`` describe a dataset called ``--dataset-name`` (default
    ``dataset``). When a configuration file defines exactly one dataset under another
    name, the flags override that dataset instead of adding a second one; likewise for
    ``--pairs`` and the reference. An explicit ``--dataset-name``/``--reference-name``
    is respected. A reference keeps the kind given in the file.
    """
    for key in ("datasets", "references"):
        ov = overrides.get(key)
        file_entries = raw.get(key)
        if not isinstance(ov, dict) or len(ov) != 1 or not isinstance(file_entries, dict):
            continue
        ((name, value),) = ov.items()
        default_name = "dataset" if key == "datasets" else "reference"
        if name != default_name or name in file_entries or len(file_entries) != 1:
            continue
        (file_name,) = file_entries
        value = dict(value)
        if key == "references":
            value.pop("kind", None)
        overrides[key] = {file_name: value}


def check_inputs_exist(cfg: Config, pipelines: list[str]) -> None:
    """Fail early (exit 2) when a configured input file is missing."""
    missing: list[str] = []

    def chk(label: str, p: Path | None) -> None:
        if p is not None and not Path(p).exists():
            missing.append(f"{label}: {p}")

    g = cfg.genotypes
    chk("genotypes.pgen_list", g.pgen_list)
    if g.pfile is not None:
        for suf in (".pgen", ".psam", ".pvar.idx.gz"):
            chk("genotypes.pfile", Path(str(g.pfile) + suf))
    refs = {cfg.reference_for(p) for p in pipelines}  # type: ignore[arg-type]
    for r in sorted(refs):
        rc = cfg.references[r]
        chk(f"references.{r}.pairs", rc.pairs)
        chk(f"references.{r}.variants", rc.variants)
    dss = {d for p in pipelines for d in cfg.datasets_for(p)}  # type: ignore[arg-type]
    for d in sorted(dss):
        dc = cfg.datasets[d]
        for f in ("pheno", "cov", "trait_annot", "sample_map", "trait_map"):
            chk(f"datasets.{d}.{f}", getattr(dc, f))
    if "match" in pipelines:
        chk("pipelines.match.prs_samples", cfg.pipelines.match.prs_samples)
    for p in ("cis", "trans"):
        if p in pipelines:
            for s in getattr(cfg.pipelines, p).strategies:
                chk(f"pipelines.{p}.strategies.{s.name}.pheno", s.pheno)
                if s.cov not in (None, "default") and s.cov.file is not None:
                    chk(f"pipelines.{p}.strategies.{s.name}.cov.file", s.cov.file)
    if "trans" in pipelines:
        t = cfg.pipelines.trans
        if isinstance(t.variants, Path):
            chk("pipelines.trans.variants", t.variants)
        if isinstance(t.traits, Path):
            chk("pipelines.trans.traits", t.traits)
    if missing:
        raise ConfigError("input files not found:\n  " + "\n  ".join(missing))


def prepare_run(
    opts: RunOptions,
    pipelines: list[str],
    *,
    lock: contextlib.ExitStack | None = None,
    dry_run: bool = False,
    record_settings: bool = True,
) -> tuple[Config, Path, StepGraph, dict[str, list[str]], Path]:
    """Resolve config and qpgen, write the resolved config and Makefile.

    With ``lock``, the run directory's lock (:func:`run_lock`) is taken before anything in
    it is written and held until ``lock`` closes. With ``dry_run`` existing step
    signatures and stamps are left untouched (see :func:`write_pipeline_makefile`).
    ``anqchor.resolved.yaml`` records the settings the run's outputs were made with, so
    without ``record_settings`` (a dry run, ``anqchor steps``) an existing one is kept.

    Returns (config, run_dir, graph, groups, qpgen binary); ``groups`` has one entry per
    pipeline that runs (``--dataset``/``--strategy`` may skip some)."""
    try:
        cfg, pipelines = build_config(opts, pipelines)
        check_inputs_exist(cfg, pipelines)
    except ConfigError as e:
        raise fail(str(e), EXIT_CONFIG) from e
    try:
        located = locate(cfg.qpgen.binary)
    except QpgenNotFound as e:
        raise fail(str(e), EXIT_QPGEN) from e
    run_dir = Path(cfg.out_dir)
    run_dir.mkdir(parents=True, exist_ok=True)
    if lock is not None:
        lock.enter_context(run_lock(run_dir))
    resolved = run_dir / Layout.resolved_config
    if record_settings or not resolved.exists():
        write_resolved_config(cfg, resolved)
    ctx = BuildContext(cfg=cfg, qpgen=located.path, qpgen_sha256=sha256_file(located.path))
    try:
        graph, groups = build_graph(ctx, pipelines, reports=opts.reports, index=len(pipelines) > 1)
    except ValueError as e:
        raise fail(str(e), EXIT_CONFIG) from e
    write_pipeline_makefile(graph, run_dir, groups, dry_run=dry_run)
    return cfg, run_dir, graph, groups, located.path


def execute(opts: RunOptions, pipelines: list[str], targets: list[str] | None = None) -> None:
    """Build and run the requested pipelines; exits with the documented exit codes."""
    with contextlib.ExitStack() as lock:
        _execute(opts, pipelines, targets, lock)


def _execute(
    opts: RunOptions, pipelines: list[str], targets: list[str] | None, lock: contextlib.ExitStack
) -> None:
    cfg, run_dir, graph, groups, qpgen = prepare_run(
        opts, pipelines, lock=lock, dry_run=opts.dry_run, record_settings=not opts.dry_run
    )
    pipelines = [p for p in pipelines if p in groups]  # --dataset/--strategy may skip some
    if targets is None:
        targets = list(pipelines)
    if opts.until and opts.only:
        raise fail("--until and --only cannot be combined", EXIT_CONFIG)
    invalidated = invalidate_missing_outputs(graph, run_dir, dry_run=opts.dry_run)
    for n in invalidated:
        log.info("output of %s is missing; it will run again", n)
    # A dry run changes no signature or stamp; make is told which signatures would have
    # changed (-W) so that it reports the same steps a real run would execute.
    what_if = (
        sorted(set(changed_signatures(graph, run_dir)) | set(invalidated)) if opts.dry_run else []
    )
    pending = set(what_if)
    for n in what_if:
        pending.update(graph.downstream(n))
    assume_old: list[str] = []
    force = opts.force
    if opts.until:
        if opts.until not in graph:
            raise fail(f"--until: unknown step {opts.until!r} (see 'anqchor steps')", EXIT_CONFIG)
        targets = [opts.until]
    if opts.only:
        if opts.only not in graph:
            raise fail(f"--only: unknown step {opts.only!r} (see 'anqchor steps')", EXIT_CONFIG)
        upstream = graph.closure([opts.only])[:-1]
        up_status = step_status(graph, run_dir, upstream)
        missing_up = [u for u in upstream if up_status[u] is not StepStatus.DONE or u in pending]
        if missing_up:
            raise fail(
                f"--only {opts.only}: upstream steps not run (or out of date): "
                f"{', '.join(missing_up)}",
                EXIT_CONFIG,
            )
        # 'make -B <step>' with every upstream stamp marked old (-o): exactly this step
        # runs, nothing upstream does, and no stamp is touched before make starts (so a
        # dry run changes nothing).
        assume_old = upstream
        force = True
        targets = [opts.only]
    started = now_iso()
    config_text = (run_dir / Layout.resolved_config).read_text(encoding="utf-8")
    result = run_make(
        run_dir,
        targets,
        jobs=opts.jobs,
        dry_run=opts.dry_run,
        force=force,
        keep_going=opts.keep_going,
        assume_old=assume_old,
        what_if=what_if,
        env={THREADS_ENV: str(cfg.qpgen.threads)},
    )
    if not opts.dry_run:
        write_run_json(
            run_dir,
            config_text=config_text,
            qpgen_binary=qpgen,
            argv=sys.argv,
            started=started,
            finished=now_iso(),
            returncode=result.returncode,
            graph=graph,
        )
    if not result.ok:
        if result.failed_step:
            where = f" in step {result.failed_step}; log: {result.log}"
        elif result.message:
            where = f": {result.message}"
        else:
            where = ""
        raise fail(f"run failed{where}", EXIT_STEP_FAILED)
    if opts.dry_run:
        would = (
            ", ".join(result.would_run) if result.would_run else "nothing (all steps up to date)"
        )
        typer.echo(f"dry run: would run: {would}")
        typer.echo(f"dry run: Makefile written to {run_dir / Layout.makefile}")
        return
    ran = ", ".join(result.ran) if result.ran else "nothing (all steps up to date)"
    typer.echo(f"ran: {ran}")
    for p in pipelines:
        rep = run_dir / Layout.report(p)
        if not rep.exists():
            continue
        if f"report.{p}" in graph or not results_newer_than(run_dir, graph, rep):
            typer.echo(f"{p} report: {rep}")
        else:  # 'match run'/'cis run' leave an earlier report as it was
            typer.echo(
                f"{p} report: {rep} is older than the results; update it with "
                f"'anqchor {p} report --out {run_dir}'"
            )
    idx = run_dir / Layout.report("index")
    if idx.exists():
        typer.echo(f"index: {idx}")


def results_newer_than(run_dir: Path, graph: StepGraph, path: Path) -> bool:
    """True when a step of ``graph`` finished (its stamp was touched) after ``path`` was
    last written, i.e. a report rendered earlier no longer shows the current results."""
    from anqchor.engine.minimake import stamp_path

    mtime = path.stat().st_mtime_ns
    for name in graph.names():
        st = run_dir / stamp_path(name)
        if st.exists() and st.stat().st_mtime_ns > mtime:
            return True
    return False


LOCK_FILE = ".anqchor.lock"


@contextlib.contextmanager
def run_lock(run_dir: Path) -> Iterator[None]:
    """Hold an exclusive lock on ``<run dir>/.anqchor.lock`` while a command rewrites the
    Makefile, signatures and stamps and runs make. A second invocation on the same run
    directory fails at once (exit 2) instead of running shared steps concurrently."""
    fd = os.open(run_dir / LOCK_FILE, os.O_RDWR | os.O_CREAT, 0o644)
    try:
        try:
            fcntl.flock(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except OSError as e:
            if e.errno in (errno.EWOULDBLOCK, errno.EAGAIN, errno.EACCES):
                raise fail(
                    f"another anqchor command is running in {run_dir} (lock {run_dir / LOCK_FILE}); "
                    "wait for it to finish",
                    EXIT_CONFIG,
                ) from e
            raise
        with contextlib.suppress(OSError):
            os.ftruncate(fd, 0)
            os.write(fd, f"{os.getpid()}\n".encode())
        yield
    finally:
        os.close(fd)


def show_steps(opts: RunOptions, pipelines: list[str]) -> None:
    """Print every step of the run with its status (``done``, ``stale``, ``missing``)."""
    with contextlib.ExitStack() as lock:
        _cfg, run_dir, graph, _groups, _ = prepare_run(
            opts, pipelines, lock=lock, record_settings=False
        )
        _show_steps(run_dir, graph)


def _show_steps(run_dir: Path, graph: StepGraph) -> None:
    status = step_status(graph, run_dir)
    width = max((len(n) for n in graph.names()), default=10)
    colors = {"done": typer.colors.GREEN, "stale": typer.colors.YELLOW, "missing": typer.colors.RED}
    for step in graph:
        st = status[step.name].value
        typer.echo(
            f"{step.name:<{width}}  "
            + typer.style(f"{st:<7}", fg=colors[st])
            + f"  {step.description}"
        )
    typer.echo(f"\nrun directory: {run_dir}   logs: {run_dir / log_path('<step>')}")


def config_dict(cfg: Config) -> dict[str, Any]:
    """The configuration as a plain dict (paths as strings)."""
    return config_to_dict(cfg)


def render_report(run_dir: Path, pipeline: str, plotlyjs: str | None = None) -> Path:
    """Re-render ``report/<pipeline>.html`` of an existing run from its summaries.

    Reads the run's resolved configuration, rebuilds the report step's definition and runs
    it directly (no make, no qpgen), so it is fast enough to iterate on reports. Fails
    with exit code 2 when the run or its summaries do not exist."""
    from anqchor.config.load import load_resolved_config
    from anqchor.pipelines import pipeline_module
    from anqchor.steps import StepError, run_step_now

    run_dir = Path(run_dir)
    if not (run_dir / Layout.resolved_config).exists():
        raise fail(
            f"{run_dir} is not an anqchor run directory (no {Layout.resolved_config}); "
            f"run 'anqchor {pipeline} run --out {run_dir}' first",
            EXIT_CONFIG,
        )
    try:
        cfg = load_resolved_config(run_dir)
    except ConfigError as e:
        raise fail(str(e), EXIT_CONFIG) from e
    if plotlyjs is not None:
        cfg.report.plotlyjs = plotlyjs  # type: ignore[assignment]
    qpgen = cfg.qpgen.binary or Path("qpgentools")  # report steps never run qpgen
    module = pipeline_module(pipeline)
    steps = [
        s
        for s in module.build(BuildContext(cfg=cfg, qpgen=Path(qpgen)))
        if s.name == module.REPORT_STEP
    ]
    if not steps:
        raise fail(f"pipeline {pipeline} has no report step", EXIT_CONFIG)
    with contextlib.ExitStack() as lock:
        lock.enter_context(run_lock(run_dir))
        try:
            run_step_now(steps[0], run_dir)
        except StepError as e:
            raise fail(str(e), EXIT_CONFIG) from e
    return run_dir / Layout.report(pipeline)


def render_cis_comparison(compare_dir: Path, plotlyjs: str | None = None) -> Path:
    """Render ``report/cis.html`` of a ``cis compare`` directory from its tables."""
    from anqchor.engine.step import Step
    from anqchor.steps import StepError, run_step_now
    from anqchor.steps.cis_compare import COMPARISON_DATASET, SUMMARY_FILES
    from anqchor.util.files import read_json

    compare_dir = Path(compare_dir)
    summary = read_json(compare_dir / "summary.json")
    step = Step(
        name="report.cis",
        kind="report_cis",
        inputs=tuple(Path(n) for n in (*SUMMARY_FILES, "runs.tsv")),
        outputs=(Layout.report("cis"),),
        params={
            "datasets": [COMPARISON_DATASET],
            "summary_dirs": {COMPARISON_DATASET: "."},
            "comparison": True,
            "reference": summary.get("reference"),
            "engine": summary.get("engine"),
            "plotlyjs": plotlyjs or "embed",
        },
        pipeline="cis",
    )
    try:
        run_step_now(step, compare_dir)
    except StepError as e:
        raise fail(str(e), EXIT_CONFIG) from e
    return compare_dir / Layout.report("cis")


def render_trans_comparison(compare_dir: Path, plotlyjs: str | None = None) -> Path:
    """Render ``report/trans.html`` of a ``trans compare`` directory from its tables."""
    from anqchor.engine.step import Step
    from anqchor.steps import StepError, run_step_now
    from anqchor.steps.trans_compare import TABLES

    compare_dir = Path(compare_dir)
    step = Step(
        name="report.trans",
        kind="report_trans",
        inputs=tuple(Layout.compare_table(n) for n in TABLES),
        outputs=(Layout.report("trans"),),
        params={"comparison": True, "plotlyjs": plotlyjs or "embed"},
        pipeline="trans",
    )
    try:
        run_step_now(step, compare_dir)
    except StepError as e:
        raise fail(str(e), EXIT_CONFIG) from e
    return compare_dir / Layout.report("trans")
