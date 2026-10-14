"""``anqchor`` command-line entry point."""

from __future__ import annotations

from pathlib import Path
from typing import Annotated

import typer

from anqchor import __version__
from anqchor.cli import options as O
from anqchor.cli.execute import RunOptions, execute, show_steps
from anqchor.cli.options import flag_overrides, set_nested
from anqchor.util.log import setup_logging

app = typer.Typer(
    name="anqchor",
    help="Genetically anchored quality control for large-scale multi-omics xQTL analysis.",
    no_args_is_help=True,
    add_completion=False,
    pretty_exceptions_enable=False,
)

_state: dict[str, object] = {"qpgen": None}


def _version(value: bool) -> None:
    if value:
        typer.echo(f"anqchor {__version__}")
        raise typer.Exit()


@app.callback()
def root(
    version: Annotated[
        bool,
        typer.Option(
            "--version", callback=_version, is_eager=True, help="Show the version and exit."
        ),
    ] = False,
    verbose: Annotated[bool, typer.Option("--verbose", "-v", help="Debug logging.")] = False,
    qpgen: Annotated[
        Path | None, typer.Option("--qpgen", help="Path to qpgentools (overrides discovery).")
    ] = None,
) -> None:
    """anqchor: sample identity matching and xQTL strategy evaluation anchored on known xQTLs."""
    setup_logging(verbose)
    _state["qpgen"] = qpgen


def _opts(**kw: object) -> RunOptions:
    opts = RunOptions(**{k: v for k, v in kw.items() if k in RunOptions.__dataclass_fields__})  # type: ignore[arg-type]
    opts.dataset = list(kw.get("dataset") or [])  # type: ignore[call-overload]
    opts.strategy = list(kw.get("strategy") or [])  # type: ignore[call-overload]
    if opts.qpgen is None and _state["qpgen"] is not None:
        opts.qpgen = _state["qpgen"]  # type: ignore[assignment]
    return opts


def _parse_lambda(text: str) -> str | float:
    """``--lambda`` value: ``auto`` or a number (range checked by the config model)."""
    if text == "auto":
        return text
    try:
        return float(text)
    except ValueError:
        raise typer.BadParameter(
            f"expected 'auto' or a number, got {text!r}", param_hint="--lambda"
        ) from None


# --------------------------------------------------------------------------- match
match_app = typer.Typer(
    help="Sample identity matching between genotypes and omics data.", no_args_is_help=True
)
app.add_typer(match_app, name="match")


def _match_command(reports: bool):  # type: ignore[no-untyped-def]
    """Build the ``match run`` (summaries only) or ``match all`` (summaries + report) command."""

    def match_cmd(
        config: O.Config = None,
        out: O.Out = None,
        jobs: O.Jobs = 1,
        threads: O.Threads = None,
        dataset: O.Dataset = None,
        dry_run: O.DryRun = False,
        until: O.Until = None,
        only: O.Only = None,
        force: O.Force = False,
        keep_going: O.KeepGoing = False,
        plotlyjs: O.PlotlyJs = None,
        geno_list: O.GenoList = None,
        pfile: O.PFile = None,
        pheno: O.Pheno = None,
        pheno_format: O.PhenoFormat = None,
        cov: O.Cov = None,
        cov_format: O.CovFormat = None,
        cat_cov_names: O.CatCovNames = None,
        trait_id_rule: O.TraitIdRule = None,
        sample_map: O.SampleMap = None,
        pairs: O.Pairs = None,
        dataset_name: O.DatasetName = "dataset",
        reference_name: O.RefName = "reference",
        mahalanobis: Annotated[
            bool | None,
            typer.Option(
                "--mahalanobis/--independent",
                help="Whiten correlated traits (slower) instead of scoring them as independent (default).",
            ),
        ] = None,
        lam: Annotated[
            str | None,
            typer.Option("--lambda", help="With --mahalanobis: 'auto' or a value in [0, 1]."),
        ] = None,
        rint: Annotated[
            bool | None,
            typer.Option(
                "--rint/--no-rint",
                help="Rank-based inverse normal transform after covariate adjustment.",
            ),
        ] = None,
        use_cov: Annotated[
            bool | None,
            typer.Option(
                "--use-cov/--no-use-cov", help="Adjust phenotypes for covariates before matching."
            ),
        ] = None,
        z_threshold: Annotated[
            float | None, typer.Option("--z-threshold", help="Z-score for a lenient self match.")
        ] = None,
        z_diff: Annotated[
            float | None,
            typer.Option("--z-diff", help="Z-score gap that makes a new best match clear."),
        ] = None,
        prs_samples: Annotated[
            Path | None,
            typer.Option(
                "--prs-samples",
                exists=True,
                dir_okay=False,
                help="Genotype IDs to compute PRS for (default: all).",
            ),
        ] = None,
    ) -> None:
        ov = flag_overrides(
            geno_list=geno_list,
            pfile=pfile,
            pheno=pheno,
            pheno_format=pheno_format,
            cov=cov,
            cov_format=cov_format,
            cat_cov_names=cat_cov_names,
            trait_id_rule=trait_id_rule,
            sample_map=sample_map,
            pairs=pairs,
            dataset_name=dataset_name,
            reference_name=reference_name,
        )
        m = "pipelines.match."
        set_nested(ov, m + "mahalanobis", mahalanobis)
        if lam is not None:
            set_nested(ov, m + "lambda", _parse_lambda(lam))
        set_nested(ov, m + "rint", rint)
        set_nested(ov, m + "use_cov", use_cov)
        set_nested(ov, m + "z_threshold", z_threshold)
        set_nested(ov, m + "z_diff", z_diff)
        set_nested(ov, m + "prs_samples", str(prs_samples) if prs_samples else None)
        execute(
            _opts(
                config=config,
                out=out,
                jobs=jobs,
                threads=threads,
                dataset=dataset,
                dry_run=dry_run,
                until=until,
                only=only,
                force=force,
                keep_going=keep_going,
                plotlyjs=plotlyjs,
                overrides=ov,
                reports=reports,
            ),
            ["match"],
        )

    return match_cmd


match_app.command(
    "run",
    help="Run sample identity matching up to the summary statistics (pair-prs, match-prs-pheno, "
    "summary tables). No HTML; render it with 'anqchor match report --out DIR'.",
)(_match_command(reports=False))
match_app.command(
    "all", help="'match run' followed by 'match report': analysis, summaries and the HTML report."
)(_match_command(reports=True))


@match_app.command("report")
def match_report(
    out: Annotated[
        Path, typer.Option("--out", "-o", help="Run directory written by 'anqchor match run'.")
    ],
    plotlyjs: O.PlotlyJs = None,
) -> None:
    """(Re-)render report/match.html from an existing run's summaries. Runs no analysis."""
    from anqchor.cli.execute import render_report

    path = render_report(out, "match", plotlyjs=plotlyjs)
    typer.echo(f"match report: {path}")


# ----------------------------------------------------------------------------- cis
cis_app = typer.Typer(
    help="cis-xQTL power evaluation by replication of known pairs.", no_args_is_help=True
)
app.add_typer(cis_app, name="cis")


def _cis_command(reports: bool):  # type: ignore[no-untyped-def]
    """Build ``cis run`` (summaries only) or ``cis all`` (summaries + report)."""

    def cis_cmd(
        config: O.Config = None,
        out: O.Out = None,
        jobs: O.Jobs = 1,
        threads: O.Threads = None,
        dataset: O.Dataset = None,
        strategy: O.Strategy = None,
        dry_run: O.DryRun = False,
        until: O.Until = None,
        only: O.Only = None,
        force: O.Force = False,
        keep_going: O.KeepGoing = False,
        plotlyjs: O.PlotlyJs = None,
        geno_list: O.GenoList = None,
        pfile: O.PFile = None,
        pheno: O.Pheno = None,
        pheno_format: O.PhenoFormat = None,
        cov: O.Cov = None,
        cov_format: O.CovFormat = None,
        cat_cov_names: O.CatCovNames = None,
        cov_include: Annotated[
            str | None,
            typer.Option(
                "--cov-include",
                help="Use only these covariates of --cov: comma-separated names; ranges such as "
                "PEER1-10 expand. Default: every covariate in --cov.",
            ),
        ] = None,
        trait_annot: O.TraitAnnot = None,
        trait_id_rule: O.TraitIdRule = None,
        pairs: O.Pairs = None,
        dataset_name: O.DatasetName = "dataset",
        reference_name: O.RefName = "reference",
        rint: Annotated[
            bool | None,
            typer.Option(
                "--rint/--no-rint",
                help="Apply a rank-based inverse normal transform to the phenotypes after "
                "covariate adjustment. Default: test the phenotypes as given (one transform per run; compare "
                "runs with 'cis compare').",
            ),
        ] = None,
        assoc_engine: Annotated[
            str | None,
            typer.Option("--assoc-engine", help="'pair' (pair-assoc) or 'rect' (rect-assoc)."),
        ] = None,
    ) -> None:
        ov = flag_overrides(
            geno_list=geno_list,
            pfile=pfile,
            pheno=pheno,
            pheno_format=pheno_format,
            cov=cov,
            cov_format=cov_format,
            cat_cov_names=cat_cov_names,
            trait_annot=trait_annot,
            trait_id_rule=trait_id_rule,
            pairs=pairs,
            dataset_name=dataset_name,
            reference_name=reference_name,
        )
        set_nested(ov, "pipelines.cis.assoc_engine", assoc_engine)
        set_nested(ov, "pipelines.cis.rint", rint)
        if cov_include is not None:
            names = [c.strip() for c in cov_include.split(",") if c.strip()]
            if not names:
                raise typer.BadParameter(
                    "give at least one covariate name", param_hint="--cov-include"
                )
            set_nested(
                ov, "pipelines.cis.strategies", [{"name": "default", "cov": {"include": names}}]
            )
        execute(
            _opts(
                config=config,
                out=out,
                jobs=jobs,
                threads=threads,
                dataset=dataset,
                strategy=strategy,
                dry_run=dry_run,
                until=until,
                only=only,
                force=force,
                keep_going=keep_going,
                plotlyjs=plotlyjs,
                overrides=ov,
                reports=reports,
            ),
            ["cis"],
        )

    return cis_cmd


cis_app.command(
    "run",
    help="Evaluate replication of known cis pairs for one set of inputs (genotypes, phenotypes, "
    "covariates, phenotype transform, reference), up to the summary tables. No HTML: use 'cis report' "
    "or 'cis compare'.",
)(_cis_command(reports=False))
cis_app.command("all", help="'cis run' followed by 'cis report' on the same run directory.")(
    _cis_command(reports=True)
)


@cis_app.command("compare")
def cis_compare(
    runs: Annotated[
        list[Path],
        typer.Argument(
            help="Run directories written by 'anqchor cis run'.", exists=True, file_okay=False
        ),
    ],
    out: Annotated[Path, typer.Option("--out", "-o", help="Comparison directory to write.")],
    label: Annotated[
        list[str] | None,
        typer.Option(
            "--label",
            help="Label for each run, in order (repeat once per run). Default: the run directory names (<parent>_<name> when two are the same).",
        ),
    ] = None,
    report: Annotated[
        bool, typer.Option("--report/--no-report", help="Also write report/cis.html.")
    ] = True,
    plotlyjs: O.PlotlyJs = None,
) -> None:
    """Compare finished cis runs (any inputs) and write comparison tables and the HTML report."""
    import sys

    from anqchor.cli.execute import EXIT_CONFIG, fail, render_cis_comparison
    from anqchor.steps import StepError
    from anqchor.steps.cis_compare import compare_runs

    try:
        summary = compare_runs(
            list(runs), out, labels=list(label) if label else None, argv=sys.argv
        )
    except StepError as e:
        raise fail(str(e), EXIT_CONFIG) from e
    for w in summary.get("warnings", []):
        typer.secho(f"warning: {w}", fg=typer.colors.YELLOW, err=True)
    typer.echo(
        f"compared {len(summary['runs'])} run(s): {', '.join(r['label'] for r in summary['runs'])}"
    )
    typer.echo(f"tables: {out}")
    if report:
        typer.echo(f"cis report: {render_cis_comparison(out, plotlyjs)}")


@cis_app.command("report")
def cis_report(
    out: Annotated[
        Path, typer.Option("--out", "-o", help="A cis run directory or a 'cis compare' directory.")
    ],
    plotlyjs: O.PlotlyJs = None,
) -> None:
    """(Re-)render report/cis.html from existing tables. Runs no analysis."""
    from anqchor.cli.execute import EXIT_CONFIG, fail, render_cis_comparison, render_report
    from anqchor.steps.compare import comparison_kind

    kind = comparison_kind(out)
    if kind == "trans":
        raise fail(
            f"{out} is a trans comparison; use 'anqchor trans report --out {out}'", EXIT_CONFIG
        )
    path = render_cis_comparison(out, plotlyjs) if kind else render_report(out, "cis", plotlyjs)
    typer.echo(f"cis report: {path}")


# --------------------------------------------------------------------------- trans
trans_app = typer.Typer(
    help="trans-xQTL evaluation: replicated reference pairs, hotspots and calibration.",
    no_args_is_help=True,
)
app.add_typer(trans_app, name="trans")


def _trans_command(reports: bool):  # type: ignore[no-untyped-def]
    """Build ``trans run`` (summaries only) or ``trans all`` (summaries + report)."""

    def trans_cmd(
        config: O.Config = None,
        out: O.Out = None,
        jobs: O.Jobs = 1,
        threads: O.Threads = None,
        dataset: O.Dataset = None,
        strategy: O.Strategy = None,
        dry_run: O.DryRun = False,
        until: O.Until = None,
        only: O.Only = None,
        force: O.Force = False,
        keep_going: O.KeepGoing = False,
        plotlyjs: O.PlotlyJs = None,
        geno_list: O.GenoList = None,
        pfile: O.PFile = None,
        pheno: O.Pheno = None,
        pheno_format: O.PhenoFormat = None,
        cov: O.Cov = None,
        cov_format: O.CovFormat = None,
        cat_cov_names: O.CatCovNames = None,
        cov_include: Annotated[
            str | None,
            typer.Option(
                "--cov-include",
                help="Use only these covariates of --cov: comma-separated names; ranges such as "
                "PEER1-10 expand. Default: every covariate in --cov.",
            ),
        ] = None,
        trait_annot: O.TraitAnnot = None,
        trait_id_rule: O.TraitIdRule = None,
        pairs: O.Pairs = None,
        variants: Annotated[
            Path | None,
            typer.Option(
                "--variants",
                exists=True,
                dir_okay=False,
                help="Variant list (CHROM:POS:REF:ALT per line) instead of the reference's variants.",
            ),
        ] = None,
        traits: Annotated[
            str | None,
            typer.Option(
                "--traits",
                help="'reference' (default), 'all', or a file of trait IDs ('file:' prefix optional).",
            ),
        ] = None,
        dataset_name: O.DatasetName = "dataset",
        reference_name: O.RefName = "reference",
        rint: Annotated[
            bool | None,
            typer.Option(
                "--rint/--no-rint",
                help="Apply a rank-based inverse normal transform to the phenotypes after covariate "
                "adjustment. Default: test the phenotypes as given (one transform per run).",
            ),
        ] = None,
        cis_window_bp: Annotated[
            int | None,
            typer.Option(
                "--cis-window-bp",
                help="A variant within this distance of a gene's start on the same chromosome is cis "
                "(default 1000000); every other pair is trans.",
            ),
        ] = None,
        min_maf: Annotated[
            float | None, typer.Option("--min-maf", help="Minimum minor allele frequency.")
        ] = None,
        hotspot_bin_bp: Annotated[
            int | None,
            typer.Option(
                "--hotspot-bin-bp",
                min=1,
                help="Bin width that thins trans signals for hotspots (default 1000000).",
            ),
        ] = None,
        hotspot_min_traits: Annotated[
            int | None,
            typer.Option(
                "--hotspot-min-traits",
                min=1,
                help="Genes with a trans signal that make a bin a hotspot (default 5; the report "
                "lets you change it).",
            ),
        ] = None,
    ) -> None:
        # Without a config file or --pairs, --variants defines a trans reference; otherwise
        # it replaces the reference's variants for the trans rectangle.
        variants_as_reference = variants is not None and pairs is None and config is None
        ov = flag_overrides(
            geno_list=geno_list,
            pfile=pfile,
            pheno=pheno,
            pheno_format=pheno_format,
            cov=cov,
            cov_format=cov_format,
            cat_cov_names=cat_cov_names,
            trait_annot=trait_annot,
            trait_id_rule=trait_id_rule,
            pairs=pairs,
            ref_kind="trans" if variants_as_reference else "cis",
            variants=variants if variants_as_reference else None,
            dataset_name=dataset_name,
            reference_name=reference_name,
        )
        t = "pipelines.trans."
        if variants is not None and not variants_as_reference:
            set_nested(ov, t + "variants", str(variants.resolve()))
        if traits is not None:
            # 'reference', 'all', 'reference:<name>', 'file:<path>' or a path; relative paths
            # are resolved against the current directory by the config loader.
            set_nested(ov, t + "traits", traits)
        set_nested(ov, t + "rint", rint)
        set_nested(ov, t + "cis_window_bp", cis_window_bp)
        set_nested(ov, t + "min_maf", min_maf)
        set_nested(ov, t + "hotspot_bin_bp", hotspot_bin_bp)
        set_nested(ov, t + "hotspot_min_traits", hotspot_min_traits)
        if cov_include is not None:
            names = [c.strip() for c in cov_include.split(",") if c.strip()]
            if not names:
                raise typer.BadParameter(
                    "give at least one covariate name", param_hint="--cov-include"
                )
            set_nested(ov, t + "strategies", [{"name": "default", "cov": {"include": names}}])
        execute(
            _opts(
                config=config,
                out=out,
                jobs=jobs,
                threads=threads,
                dataset=dataset,
                strategy=strategy,
                dry_run=dry_run,
                until=until,
                only=only,
                force=force,
                keep_going=keep_going,
                plotlyjs=plotlyjs,
                overrides=ov,
                reports=reports,
            ),
            ["trans"],
        )

    return trans_cmd


trans_app.command(
    "run",
    help="Evaluate one set of inputs (genotypes, phenotypes, covariates, phenotype transform, "
    "reference) over a variant x gene rectangle, up to the summary tables: replicated reference "
    "pairs (cis and trans), hotspots and calibration. No HTML: use 'trans compare' or 'trans report'.",
)(_trans_command(reports=False))


@trans_app.command("compare")
def trans_compare(
    runs: Annotated[
        list[Path],
        typer.Argument(
            help="Run directories written by 'anqchor trans run'.", exists=True, file_okay=False
        ),
    ],
    out: Annotated[Path, typer.Option("--out", "-o", help="Comparison directory to write.")],
    label: Annotated[
        list[str] | None,
        typer.Option(
            "--label",
            help="Label for each run, in order (repeat once per run). Default: the run directory names (<parent>_<name> when two are the same).",
        ),
    ] = None,
    report: Annotated[
        bool, typer.Option("--report/--no-report", help="Also write report/trans.html.")
    ] = True,
    plotlyjs: O.PlotlyJs = None,
) -> None:
    """Compare finished trans runs (any inputs) and write comparison tables and the HTML report."""
    import sys

    from anqchor.cli.execute import EXIT_CONFIG, fail, render_trans_comparison
    from anqchor.steps import StepError
    from anqchor.steps.trans_compare import compare_runs

    try:
        summary = compare_runs(
            list(runs), out, labels=list(label) if label else None, argv=sys.argv
        )
    except StepError as e:
        raise fail(str(e), EXIT_CONFIG) from e
    for w in summary.get("warnings", []):
        typer.secho(f"warning: {w}", fg=typer.colors.YELLOW, err=True)
    typer.echo(f"compared {len(summary['runs'])} run(s): {', '.join(summary['labels'])}")
    typer.echo(f"tables: {out}")
    if report:
        typer.echo(f"trans report: {render_trans_comparison(out, plotlyjs)}")


@trans_app.command("report")
def trans_report(
    out: Annotated[
        Path,
        typer.Option("--out", "-o", help="A trans run directory or a 'trans compare' directory."),
    ],
    plotlyjs: O.PlotlyJs = None,
) -> None:
    """(Re-)render report/trans.html from existing tables. Runs no analysis."""
    from anqchor.cli.execute import EXIT_CONFIG, fail, render_report, render_trans_comparison
    from anqchor.steps.compare import comparison_kind

    kind = comparison_kind(out)
    if kind == "cis":
        raise fail(f"{out} is a cis comparison; use 'anqchor cis report --out {out}'", EXIT_CONFIG)
    path = render_trans_comparison(out, plotlyjs) if kind else render_report(out, "trans", plotlyjs)
    typer.echo(f"trans report: {path}")


trans_app.command(
    "all",
    help="'trans run' followed by 'trans report' on the same run directory (one run, no comparison).",
)(_trans_command(reports=True))


# ---------------------------------------------------------------- run / steps / report
def _enabled(config: Path | None) -> list[str]:
    from anqchor.cli.execute import EXIT_CONFIG, fail
    from anqchor.config.load import ConfigError, load_config

    if config is None:
        raise fail("--config is required", EXIT_CONFIG)
    try:
        cfg = load_config(config)
    except ConfigError as e:
        raise fail(str(e), EXIT_CONFIG) from e
    pl = cfg.enabled_pipelines()
    if not pl:
        raise fail("no pipeline is enabled in the configuration", EXIT_CONFIG)
    return pl


def _pick(config: Path | None, pipeline: list[str] | None) -> list[str]:
    from anqchor.cli.execute import EXIT_CONFIG, fail
    from anqchor.pipelines import PIPELINES

    if pipeline:
        bad = [p for p in pipeline if p not in PIPELINES]
        if bad:
            raise fail(f"unknown pipeline(s): {', '.join(bad)}", EXIT_CONFIG)
        return [p for p in PIPELINES if p in pipeline]
    return _enabled(config)


PipelineOpt = Annotated[
    list[str] | None,
    typer.Option(
        "--pipeline", "-p", help="match, cis or trans (repeatable; default: enabled pipelines)."
    ),
]


@app.command("run")
def run_all(
    config: O.Config = None,
    out: O.Out = None,
    jobs: O.Jobs = 1,
    threads: O.Threads = None,
    pipeline: PipelineOpt = None,
    dataset: O.Dataset = None,
    strategy: O.Strategy = None,
    dry_run: O.DryRun = False,
    until: O.Until = None,
    only: O.Only = None,
    force: O.Force = False,
    keep_going: O.KeepGoing = False,
    plotlyjs: O.PlotlyJs = None,
) -> None:
    """Run every pipeline enabled in the configuration, sharing common steps."""
    execute(
        _opts(
            config=config,
            out=out,
            jobs=jobs,
            threads=threads,
            dataset=dataset,
            strategy=strategy,
            dry_run=dry_run,
            until=until,
            only=only,
            force=force,
            keep_going=keep_going,
            plotlyjs=plotlyjs,
        ),
        _pick(config, pipeline),
    )


@app.command("steps")
def steps_cmd(
    config: O.Config = None,
    out: O.Out = None,
    pipeline: PipelineOpt = None,
    dataset: O.Dataset = None,
    strategy: O.Strategy = None,
) -> None:
    """List the steps of a run with their status (done, stale, missing)."""
    show_steps(
        _opts(config=config, out=out, dataset=dataset, strategy=strategy), _pick(config, pipeline)
    )


@app.command("report")
def report_cmd(
    config: O.Config = None,
    out: O.Out = None,
    pipeline: PipelineOpt = None,
    plotlyjs: O.PlotlyJs = None,
    force: Annotated[bool, typer.Option("--force", help="Re-render even if up to date.")] = False,
) -> None:
    """Re-render HTML reports from existing summaries (runs missing upstream steps if needed)."""
    from anqchor.cli.execute import prepare_run
    from anqchor.engine.minimake import stamp_path

    pipes = _pick(config, pipeline)
    opts = _opts(config=config, out=out, plotlyjs=plotlyjs)
    # report/index.html exists only when several pipelines run together
    targets = [f"report.{p}" for p in pipes] + (["report.index"] if len(pipes) > 1 else [])
    if force:
        _, run_dir, _, _, _ = prepare_run(opts, pipes)
        for name in targets:
            st = run_dir / stamp_path(name)
            if st.exists():
                st.unlink()
    execute(opts, pipes, targets=targets)


@app.command("_step", hidden=True)
def step_cmd(name: str, run_dir: Annotated[Path, typer.Option("--run-dir")] = Path(".")) -> None:
    """Execute one Python step (used by the generated Makefile)."""
    import sys

    from anqchor.steps import StepError, run_step

    try:
        run_step(name, run_dir)
    except StepError as e:
        typer.secho(f"error in step {name}: {e}", fg=typer.colors.RED, err=True)
        sys.exit(1)


def _register_optional() -> None:
    """Register command groups implemented in separate modules."""
    from anqchor.cli import prepare as prepare_cli
    from anqchor.cli import qpgen as qpgen_cli

    app.add_typer(qpgen_cli.app, name="qpgen")
    app.add_typer(prepare_cli.app, name="prepare")


_register_optional()


def main() -> None:
    app()


if __name__ == "__main__":
    main()
