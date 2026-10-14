"""``anqchor qpgen build|check``."""

from __future__ import annotations

from pathlib import Path
from typing import Annotated

import typer

from anqchor.qpgen import build as qbuild
from anqchor.qpgen.locate import QpgenNotFound, check, locate, missing_flags, submodule_dir

app = typer.Typer(help="Build and check the qpgen toolkit.", no_args_is_help=True)


@app.command("build")
def build_cmd(
    src: Annotated[
        Path | None,
        typer.Option("--src", help="qpgen checkout (default: the extern/qpgen submodule)."),
    ] = None,
    jobs: Annotated[int, typer.Option("--jobs", "-j", min=1, help="Parallel compile jobs.")] = 4,
    force: Annotated[bool, typer.Option("--force", help="Rebuild htslib and qgenlib too.")] = False,
    cmake_arg: Annotated[
        list[str] | None, typer.Option("--cmake-arg", help="Extra argument for cmake (repeatable).")
    ] = None,
) -> None:
    """Build qpgentools from the submodule (htslib, qgenlib, then qpgen)."""
    src = src or submodule_dir()
    if src is None:
        typer.secho(
            "error: no qpgen checkout found; pass --src or run from a source checkout",
            fg="red",
            err=True,
        )
        raise typer.Exit(4)
    for p in qbuild.preflight():
        typer.secho(f"warning: {p}", fg="yellow", err=True)
    try:
        binary = qbuild.build(src, jobs=jobs, force=force, cmake_args=cmake_arg or [])
    except qbuild.BuildError as e:
        typer.secho(f"error: {e}", fg="red", err=True)
        raise typer.Exit(4) from e
    typer.echo(str(binary))


@app.command("check")
def check_cmd(
    qpgen: Annotated[
        Path | None, typer.Option("--qpgen", help="Binary to check (default: discovery order).")
    ] = None,
) -> None:
    """Locate qpgentools, verify its subcommands, and print its provenance."""
    from anqchor.engine.provenance import qpgen_provenance

    try:
        loc = locate(qpgen)
    except QpgenNotFound as e:
        typer.secho(f"error: {e}", fg="red", err=True)
        raise typer.Exit(4) from e
    try:
        missing = check(loc.path)
        flags = [] if missing else missing_flags(loc.path)
    except QpgenNotFound as e:
        typer.secho(f"error: {e}", fg="red", err=True)
        raise typer.Exit(4) from e
    prov = qpgen_provenance(loc.path)
    typer.echo(f"binary:   {loc.path}  (found via {loc.source})")
    typer.echo(f"sha256:   {prov.get('sha256')}")
    if prov.get("git_commit"):
        typer.echo(f"commit:   {prov['git_commit']}")
    if missing:
        typer.secho(f"missing subcommands: {', '.join(missing)}", fg="red", err=True)
        raise typer.Exit(4)
    if flags:
        typer.secho(
            f"missing options: {', '.join(flags)} (categorical covariates need qpgen 5b811e2 or later; "
            "update extern/qpgen and run 'anqchor qpgen build')",
            fg="red",
            err=True,
        )
        raise typer.Exit(4)
    typer.echo("status:   ok (pair-prs, match-prs-pheno, pair-assoc, rect-assoc available)")
