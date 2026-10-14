"""Builders of tiny synthetic inputs for the io/prepare/strategy unit tests.

Imported by ``tests/unit/test_io_*.py``, ``test_prepare_*.py`` and ``test_strategy_*.py``
after adding this directory to ``sys.path``.
"""

from __future__ import annotations

import gzip
import subprocess
from pathlib import Path
from typing import Any

import pytest

from anqchor.io.genotypes import find_tool
from anqchor.steps import StepContext
from anqchor.util.log import get_logger

SAMPLES = ["S1", "S2", "S3", "S4", "S5", "S6"]

# (chrom number, pos, ref, alt); written with the requested chromosome style
RECORDS = [
    ("1", 100, "A", "G"),
    ("1", 200, "C", "T"),
    ("1", 200, "C", "CA"),
    ("1", 300, "G", "A"),
    ("2", 50, "T", "C"),
    ("2", 150, "AT", "A"),
]

needs_htslib = pytest.mark.skipif(
    not (find_tool("bgzip") and find_tool("tabix")), reason="bgzip/tabix not available"
)


def write_text(path: Path, text: str) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.suffix == ".gz":
        with gzip.open(path, "wt", encoding="utf-8") as fh:
            fh.write(text)
    else:
        path.write_text(text, encoding="utf-8")
    return path


def tsv(rows: list[list[Any]]) -> str:
    return "".join("\t".join(str(x) for x in r) + "\n" for r in rows)


def make_pivar(prefix: Path, records: list[tuple[str, int, str, str]]) -> Path:
    """Write ``prefix.pvar.idx.gz`` (bgzipped, tabix-indexed, 9 columns)."""
    bgzip, tabix = find_tool("bgzip"), find_tool("tabix")
    assert bgzip and tabix
    lines = ["#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tINDEX"]
    for i, (c, p, r, a) in enumerate(records, 1):
        lines.append(f"{c}\t{p}\t{c}:{p}:{r}:{a}\t{r}\t{a}\t.\t.\t.\t{i}")
    out = Path(str(prefix) + ".pvar.idx.gz")
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "wb") as fo:
        subprocess.run(
            [bgzip, "-c"], input=("\n".join(lines) + "\n").encode(), stdout=fo, check=True
        )
    subprocess.run([tabix, "-f", "-p", "vcf", str(out)], check=True)
    return out


def make_prefix(
    prefix: Path,
    records: list[tuple[str, int, str, str]],
    samples: list[str] = SAMPLES,
    psam_header: str = "#IID\tSEX",
) -> Path:
    """A fake PLINK2 file set: empty .pgen, a psam and an indexed pvar."""
    prefix.parent.mkdir(parents=True, exist_ok=True)
    Path(str(prefix) + ".pgen").write_bytes(b"")
    rows = [psam_header] if psam_header else []
    rows += [f"{s}\tNA" if psam_header.startswith("#IID") else f"F{s}\t{s}\tNA" for s in samples]
    Path(str(prefix) + ".psam").write_text("\n".join(rows) + "\n")
    make_pivar(prefix, records)
    return prefix


def make_geno(
    root: Path, style: str = "plain", relative: bool = True, samples: list[str] = SAMPLES
) -> Path:
    """Per-chromosome fake genotypes and a 4-column pgen list; returns the list path."""
    lines = []
    for chrom in ("1", "2"):
        name = chrom if style == "plain" else "chr" + chrom
        recs = [(name, p, r, a) for c, p, r, a in RECORDS if c == chrom]
        prefix = root / "geno" / f"set.chr{chrom}"
        make_prefix(prefix, recs, samples)
        lines.append(f"{name}\t1\t1000000000\t{'geno/set.chr' + chrom if relative else prefix}")
    lst = root / "geno.list.tsv"
    lst.write_text("\n".join(lines) + "\n")
    return lst


def step_ctx(
    run_dir: Path, name: str, params: dict[str, Any], outputs: list[Path] | None = None
) -> StepContext:
    return StepContext(
        name=name,
        run_dir=run_dir,
        params=params,
        inputs=[],
        outputs=outputs or [],
        log=get_logger(f"test.{name}"),
    )
