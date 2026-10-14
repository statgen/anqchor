"""``anqchor prepare ...``: standalone input preparation utilities.

* ``index-pvar``: build the indexed pvar (``.pvar.idx.gz`` + ``.tbi``) qpgen needs.
* ``pgen-list``: write a qpgen pgen list from a per-chromosome prefix template.
* ``pheno``: convert a phenotype matrix in any supported layout to regenie.
* ``cov``: convert and encode a covariate matrix to numeric regenie layout.
"""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Annotated

import typer

from anqchor.io.cov import MISSING_POLICIES, read_and_encode, write_cov, write_cov_raw
from anqchor.io.genotypes import FULL_END, first_pvar_chrom, index_pvar
from anqchor.io.pheno import FORMATS, convert_to_regenie
from anqchor.io.tables import InputError, write_lines

app = typer.Typer(help="Standalone input preparation utilities.", no_args_is_help=True)

EXIT_INPUT = 2


def _fail(msg: str) -> typer.Exit:
    typer.secho(f"error: {msg}", fg=typer.colors.RED, err=True)
    return typer.Exit(EXIT_INPUT)


def _check_format(fmt: str) -> None:
    if fmt not in FORMATS:
        raise _fail(f"--format must be one of {', '.join(FORMATS)}")


@app.command("index-pvar")
def index_pvar_cmd(
    pfile: Annotated[
        Path, typer.Option("--pfile", help="PLINK2 prefix; reads PREFIX.pvar or PREFIX.pvar.zst.")
    ],
    force: Annotated[
        bool, typer.Option("--force", help="Overwrite an existing PREFIX.pvar.idx.gz.")
    ] = False,
) -> None:
    """Build PREFIX.pvar.idx.gz (9 columns with a 1-based INDEX) and its tabix index."""
    try:
        out, n = index_pvar(pfile, force=force)
    except InputError as e:
        raise _fail(str(e)) from e
    typer.echo(f"wrote {out} ({n} variants) and {out}.tbi")


def expand_chroms(spec: str) -> list[str]:
    """``1-22,X`` → ``['1', …, '22', 'X']`` (numeric ranges expand; other tokens kept)."""
    out: list[str] = []
    for tok in (t.strip() for t in spec.split(",")):
        if not tok:
            continue
        m = re.fullmatch(r"(\d+)-(\d+)", tok)
        if m:
            a, b = int(m.group(1)), int(m.group(2))
            if b < a:
                raise InputError(f"chromosome range {tok!r} runs backwards")
            out += [str(i) for i in range(a, b + 1)]
        else:
            out.append(tok)
    return out


@app.command("pgen-list")
def pgen_list_cmd(
    prefix_template: Annotated[
        str,
        typer.Option("--prefix-template", help="Prefix with '{chrom}', e.g. 'geno/chr{chrom}'."),
    ],
    chroms: Annotated[str, typer.Option("--chroms", help="Chromosomes, e.g. '1-22,X'.")],
    out: Annotated[Path, typer.Option("--out", help="Output pgen list.")],
    no_check: Annotated[
        bool, typer.Option("--no-check", help="Do not check that the files exist.")
    ] = False,
) -> None:
    """Write a qpgen pgen list (CHROM BEG END PREFIX) with absolute prefixes.

    With checking (the default), each prefix must have .pgen, .psam, .pvar.idx.gz and
    .tbi, and CHROM is taken from the first pvar record so it matches the pvar spelling.
    """
    if "{chrom}" not in prefix_template:
        raise _fail("--prefix-template must contain '{chrom}'")
    try:
        names = expand_chroms(chroms)
    except InputError as e:
        raise _fail(str(e)) from e
    lines: list[str] = []
    problems: list[str] = []
    for c in names:
        prefix = Path(prefix_template.replace("{chrom}", c)).expanduser().resolve()
        chrom = c
        if not no_check:
            miss = [
                s
                for s in (".pgen", ".psam", ".pvar.idx.gz", ".pvar.idx.gz.tbi")
                if not Path(f"{prefix}{s}").exists()
            ]
            if miss:
                problems.append(f"{prefix}: missing {', '.join(miss)}")
                continue
            chrom = first_pvar_chrom(f"{prefix}.pvar.idx.gz")
        lines.append(f"{chrom}\t1\t{FULL_END}\t{prefix}")
    if problems:
        raise _fail("genotype files not found:\n  " + "\n  ".join(problems))
    write_lines(out, lines)
    typer.echo(f"wrote {out} ({len(lines)} entries)")


Missing = Annotated[
    list[str] | None,
    typer.Option("--missing-str", help="String meaning 'missing' (repeatable; default NA)."),
]


@app.command("pheno")
def pheno_cmd(
    inp: Annotated[Path, typer.Option("--in", exists=True, dir_okay=False, help="Input matrix.")],
    fmt: Annotated[
        str, typer.Option("--format", help="regenie | tensorqtl | tsv-sample-row | tsv-sample-col")
    ],
    out: Annotated[Path, typer.Option("--out", help="Output regenie matrix (.tsv.gz).")],
    missing_str: Missing = None,
) -> None:
    """Convert a phenotype matrix to regenie layout (FID IID trait...)."""
    _check_format(fmt)
    try:
        res = convert_to_regenie(inp, fmt, out, missing_str or ["NA"])
    except InputError as e:
        raise _fail(str(e)) from e
    m = res.matrix
    how = "copied unchanged" if res.fast_path else "converted"
    typer.echo(
        f"wrote {out}: {len(m.samples)} samples x {len(m.features)} traits ({how}; {res.n_missing} missing values)"
    )


@app.command("cov")
def cov_cmd(
    inp: Annotated[
        Path, typer.Option("--in", exists=True, dir_okay=False, help="Input covariate matrix.")
    ],
    fmt: Annotated[
        str, typer.Option("--format", help="regenie | tensorqtl | tsv-sample-row | tsv-sample-col")
    ],
    out: Annotated[Path, typer.Option("--out", help="Output regenie covariates (.tsv.gz).")],
    categorical: Annotated[
        list[str] | None,
        typer.Option("--categorical", help="Covariate to treat as categorical (repeatable)."),
    ] = None,
    missing: Annotated[
        str, typer.Option("--missing", help="error | drop_samples | impute_mean")
    ] = "error",
    missing_str: Missing = None,
    json_out: Annotated[
        Path | None, typer.Option("--json", help="Also write the encoding summary as JSON.")
    ] = None,
    encode: Annotated[
        bool,
        typer.Option(
            "--encode/--keep-categorical",
            help="Write categorical covariates as numeric indicator columns (--encode) instead of their text levels (default; qpgen encodes them with --cat-cov-names).",
        ),
    ] = False,
) -> None:
    """Convert covariates to regenie layout: validate values, drop constant columns, handle
    missing values. Categorical covariates keep their text levels by default."""
    _check_format(fmt)
    if missing not in MISSING_POLICIES:
        raise _fail(f"--missing must be one of {', '.join(MISSING_POLICIES)}")
    try:
        enc = read_and_encode(inp, fmt, categorical or [], missing, missing_str or ["NA"])
        if enc.values.shape[1] == 0:
            raise InputError("no usable covariate column left (all constant or missing)")
    except InputError as e:
        raise _fail(str(e)) from e
    if encode:
        write_cov(out, enc)
        n_cols = enc.values.shape[1]
    else:
        assert enc.raw is not None
        write_cov_raw(out, enc.raw)
        n_cols = enc.raw.shape[1]
    if json_out is not None:
        json_out.write_text(json.dumps(enc.as_json(), indent=2) + "\n", encoding="utf-8")
    cats = ", ".join(f"{k} {list(v)}" for k, v in enc.categorical.items()) or "none"
    typer.echo(
        f"wrote {out}: {enc.values.shape[0]} samples x {n_cols} covariates; categorical: {cats}"
    )
    if enc.categorical and not encode:
        typer.echo(
            f"categorical covariates keep their text levels; use --cat-cov-names {','.join(enc.categorical_names)}"
        )
    if enc.dropped_constant:
        typer.echo(f"dropped constant: {', '.join(enc.dropped_constant)}")
