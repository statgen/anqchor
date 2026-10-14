"""Genotype-set inputs: qpgen pgen lists, psam sample lists, indexed pvar files.

Also the chromosome-style helpers (``chr21`` versus ``21``) used by every reader that
has to agree with the genotype pvar, and discovery of the htslib command-line tools
(``tabix``, ``bgzip``).
"""

from __future__ import annotations

import gzip
import os
import shutil
import subprocess
from collections.abc import Iterator
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Literal

from anqchor.io.tables import InputError, examples, is_int
from anqchor.util.files import is_gzip, open_text, resolve_near

ChromStyle = Literal["chr", "plain"]

FULL_END = 2_000_000_000
"""END used for whole-chromosome pgen-list entries (fits a signed 32-bit integer)."""

PVAR_IDX_HEADER = ("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "INDEX")
"""Columns of the 9-column indexed pvar that qpgen reads (``.pvar.idx.gz``)."""


# ------------------------------------------------------------------ chromosomes
def chrom_style(name: str) -> ChromStyle:
    """``chr`` when the chromosome name starts with ``chr`` (any case), else ``plain``."""
    return "chr" if name[:3].lower() == "chr" else "plain"


def normalize_chrom(name: str, style: ChromStyle) -> str:
    """Rewrite a chromosome name to ``style``: add or strip a leading ``chr``.

    ``23`` and ``X`` are not translated into each other.
    """
    if style == "chr":
        return name if name[:3].lower() == "chr" else "chr" + name
    return name[3:] if name[:3].lower() == "chr" and len(name) > 3 else name


# ---------------------------------------------------------------- htslib tools
def find_tool(name: str) -> str | None:
    """Absolute path of an htslib tool (``tabix``, ``bgzip``): ``$ANQCHOR_<NAME>``,
    then ``PATH``, then the htslib build inside the qpgen submodule. Setting
    ``ANQCHOR_NO_<NAME>=1`` hides the tool (used to exercise pure-Python fallbacks)."""
    if os.environ.get(f"ANQCHOR_NO_{name.upper()}"):
        return None
    env = os.environ.get(f"ANQCHOR_{name.upper()}")
    if env and Path(env).is_file() and os.access(env, os.X_OK):
        return env
    found = shutil.which(name)
    if found:
        return found
    from anqchor.qpgen.locate import submodule_dir

    sub = submodule_dir()
    if sub is not None:
        cand = sub / "submodules" / "htslib" / name
        if cand.is_file() and os.access(cand, os.X_OK):
            return str(cand)
    return None


# ------------------------------------------------------------------ pgen lists
@dataclass(frozen=True)
class PgenEntry:
    """One row of a qpgen pgen list with absolute paths."""

    chrom: str
    beg: int
    end: int
    pgen: Path
    psam: Path
    pivar: Path

    def as_json(self) -> dict[str, object]:
        d = asdict(self)
        return {k: str(v) if isinstance(v, Path) else v for k, v in d.items()}

    def standard_prefix(self) -> str | None:
        """``P`` when the files are ``P.pgen``, ``P.psam`` and ``P.pvar.idx.gz``, else ``None``."""
        p = str(self.pgen)
        if not p.endswith(".pgen"):
            return None
        prefix = p[: -len(".pgen")]
        if str(self.psam) == prefix + ".psam" and str(self.pivar) == prefix + ".pvar.idx.gz":
            return prefix
        return None

    def line(self, prefix: str | None = None) -> str:
        """4-column pgen-list line ``CHROM BEG END PREFIX`` (tab separated).

        The 4-column form is used because qpgen's 6-column order is ambiguous: its
        documentation says ``PGEN PSAM PIVAR`` but its parser reads ``PGEN PIVAR PSAM``.
        """
        prefix = prefix or self.standard_prefix()
        if prefix is None:
            raise ValueError(f"entry {self.pgen} has no standard prefix; link the files first")
        return "\t".join([self.chrom, str(self.beg), str(self.end), prefix])


def parse_pgen_list(path: str | os.PathLike[str]) -> list[PgenEntry]:
    """Parse a qpgen pgen list (4 fields ``CHROM BEG END PREFIX`` or 6 fields
    ``CHROM BEG END PGEN PSAM PIVAR``, whitespace separated; blank and ``#`` lines
    skipped). Relative paths are resolved near the list (:func:`resolve_near`)."""
    path = Path(path)
    entries: list[PgenEntry] = []
    with open_text(path) as fh:
        for lineno, line in enumerate(fh, 1):
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            f = s.split()
            if len(f) not in (4, 6):
                raise InputError(
                    f"{path}:{lineno}: expected 4 fields (CHROM BEG END PREFIX) or 6 fields "
                    f"(CHROM BEG END PGEN PSAM PIVAR), found {len(f)}"
                )
            if not (is_int(f[1]) and is_int(f[2])):
                raise InputError(
                    f"{path}:{lineno}: BEG and END must be integers, got {f[1]!r} {f[2]!r}"
                )
            if len(f) == 4:
                prefix = str(resolve_near(f[3], path))
                pgen, psam, pivar = (Path(prefix + s) for s in (".pgen", ".psam", ".pvar.idx.gz"))
            else:
                pgen, a, b = (resolve_near(x, path) for x in f[3:6])
                psam, pivar = _six_column_order(a, b)
            entries.append(PgenEntry(f[0], int(f[1]), int(f[2]), pgen, psam, pivar))
    if not entries:
        raise InputError(f"pgen list {path} has no entries")
    return entries


def _six_column_order(a: Path, b: Path) -> tuple[Path, Path]:
    """(psam, pivar) from fields 5 and 6 of a 6-column pgen list.

    qpgen's documentation gives ``PGEN PSAM PIVAR`` while its parser reads ``PGEN PIVAR
    PSAM``; the order is recognised from the file names (``.psam``, ``.gz``) and, failing
    that, from the content (the indexed pvar is gzip-compressed). Unrecognisable fields
    follow qpgen's parser.
    """
    sa, sb = str(a).lower(), str(b).lower()
    if sa.endswith(".psam") or sb.endswith(".gz"):
        return a, b
    if sb.endswith(".psam") or sa.endswith(".gz"):
        return b, a
    if a.exists() and b.exists():
        return (b, a) if is_gzip(a) else (a, b)
    return b, a


def entries_from_pfile(prefix: str | os.PathLike[str]) -> list[PgenEntry]:
    """One whole-chromosome entry per chromosome present in ``prefix.pvar.idx.gz``."""
    p = str(Path(prefix).expanduser().resolve())
    pgen, psam, pivar = (Path(p + s) for s in (".pgen", ".psam", ".pvar.idx.gz"))
    if not pivar.exists():
        raise InputError(
            f"indexed pvar not found: {pivar}. Create it with 'anqchor prepare index-pvar --pfile {prefix}'"
        )
    return [PgenEntry(c, 1, FULL_END, pgen, psam, pivar) for c in pvar_chroms(pivar)]


def missing_files(entry: PgenEntry) -> list[Path]:
    """Files of ``entry`` that do not exist (pgen, psam, pvar.idx.gz and its .tbi)."""
    need = [entry.pgen, entry.psam, entry.pivar, Path(str(entry.pivar) + ".tbi")]
    return [p for p in need if not p.exists()]


def read_prepared_pgen_list(path: str | os.PathLike[str]) -> list[PgenEntry]:
    """Read the absolute pgen list written by ``prepare.genotypes`` (4 or 6 columns)."""
    return parse_pgen_list(path)


def link_standard_prefix(entry: PgenEntry, prefix: Path) -> str:
    """Symlink the entry's files as ``prefix.pgen/.psam/.pvar.idx.gz(.tbi)``; returns
    the prefix. Used when the files do not share a standard prefix."""
    prefix.parent.mkdir(parents=True, exist_ok=True)
    pairs = [
        (entry.pgen, ".pgen"),
        (entry.psam, ".psam"),
        (entry.pivar, ".pvar.idx.gz"),
        (Path(str(entry.pivar) + ".tbi"), ".pvar.idx.gz.tbi"),
    ]
    for target, suffix in pairs:
        link = Path(str(prefix) + suffix)
        if link.is_symlink() or link.exists():
            link.unlink()
        link.symlink_to(target)
    return str(prefix)


# ------------------------------------------------------------------------- psam
def read_psam_iids(path: str | os.PathLike[str]) -> list[str]:
    """Sample IIDs of a psam file in file order.

    Accepts a ``#IID ...`` or ``#FID IID ...`` header and, without a header, the
    .fam-like layout (``FID IID ...``; a single column is taken as the IID).
    """
    iids: list[str] = []
    col: int | None = None
    with open_text(path) as fh:
        for line in fh:
            s = line.rstrip("\n").rstrip("\r")
            if not s.strip() or s.startswith("##"):
                continue
            f = s.split("\t") if "\t" in s else s.split()
            if col is None:
                if s.startswith("#"):
                    names = [x.lstrip("#").upper() for x in f]
                    if "IID" not in names:
                        raise InputError(f"{path}: psam header has no IID column: {s!r}")
                    col = names.index("IID")
                    continue
                col = 1 if len(f) >= 2 else 0
            if len(f) <= col:
                raise InputError(f"{path}: line with too few fields: {s!r}")
            iids.append(f[col])
    if not iids:
        raise InputError(f"{path}: no samples in psam")
    return iids


# ------------------------------------------------------------------------- pvar
def iter_pvar_records(pivar: str | os.PathLike[str]) -> Iterator[list[str]]:
    """Records of an indexed pvar as ``[CHROM, POS, ID, REF, ALT]`` (header skipped)."""
    with gzip.open(pivar, "rt", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            yield line.split("\t", 5)[:5]


def first_pvar_chrom(pivar: str | os.PathLike[str]) -> str:
    """Chromosome of the first record of an indexed pvar."""
    for rec in iter_pvar_records(pivar):
        return rec[0]
    raise InputError(f"{pivar}: no variant records")


def pvar_chroms(pivar: str | os.PathLike[str]) -> list[str]:
    """Chromosomes present in an indexed pvar, in file order (``tabix -l`` when available)."""
    tabix = find_tool("tabix")
    if tabix and Path(str(pivar) + ".tbi").exists():
        proc = subprocess.run(
            [tabix, "-l", str(pivar)], capture_output=True, text=True, check=False
        )
        if proc.returncode == 0:
            chroms = [c for c in proc.stdout.split() if c]
            if chroms:
                return chroms
    seen: list[str] = []
    for rec in iter_pvar_records(pivar):
        if not seen or seen[-1] != rec[0]:
            if rec[0] in seen:
                raise InputError(f"{pivar}: chromosome {rec[0]} is not contiguous (unsorted pvar)")
            seen.append(rec[0])
    if not seen:
        raise InputError(f"{pivar}: no variant records")
    return seen


def check_list_style(entries: list[PgenEntry], style: ChromStyle) -> None:
    """Raise when a pgen-list CHROM does not use the pvar chromosome style."""
    bad = [e.chrom for e in entries if chrom_style(e.chrom) != style]
    if bad:
        want = "with a 'chr' prefix" if style == "chr" else "without a 'chr' prefix"
        raise InputError(
            f"the pvar names chromosomes {want} (e.g. {first_pvar_chrom(entries[0].pivar)!r}) "
            f"but the pgen list CHROM column has {examples(bad)}; rewrite the list's CHROM "
            "column to match the pvar spelling (qpgen requires identical names)"
        )


# ------------------------------------------------------------------ index-pvar
def _decompressed_lines(src: Path) -> Iterator[str]:
    """Lines of a plain, gzip or zstd-compressed pvar."""
    if src.suffix == ".zst":
        zstd = shutil.which("zstd")
        if not zstd:
            raise InputError(f"{src} is zstd-compressed but the 'zstd' command is not on PATH")
        proc = subprocess.Popen([zstd, "-dc", str(src)], stdout=subprocess.PIPE, text=True)
        assert proc.stdout is not None
        yield from proc.stdout
        if proc.wait() != 0:
            raise InputError(f"zstd failed to decompress {src}")
        return
    with open_text(src) as fh:
        yield from fh


def index_pvar(prefix: str | os.PathLike[str], force: bool = False) -> tuple[Path, int]:
    """Build ``PREFIX.pvar.idx.gz`` (+ ``.tbi``) from ``PREFIX.pvar`` or ``PREFIX.pvar.zst``.

    The output has the 9-column header of :data:`PVAR_IDX_HEADER`; missing QUAL,
    FILTER and INFO are written as ``.`` and INDEX is the 1-based record number in pvar
    order (the pgen row). A headerless pvar is read in .bim column order. Returns the
    output path and the number of records.
    """
    prefix = str(Path(prefix).expanduser().resolve())
    out = Path(prefix + ".pvar.idx.gz")
    if out.exists() and not force:
        raise InputError(f"{out} exists; pass --force to overwrite")
    src = next(
        (Path(prefix + s) for s in (".pvar", ".pvar.zst", ".pvar.gz") if Path(prefix + s).exists()),
        None,
    )
    if src is None:
        raise InputError(f"no {prefix}.pvar or {prefix}.pvar.zst found")
    bgzip, tabix = find_tool("bgzip"), find_tool("tabix")
    if not bgzip or not tabix:
        raise InputError("index-pvar needs the 'bgzip' and 'tabix' commands (htslib) on PATH")
    tmp = Path(str(out) + ".tmp")
    try:
        n = _write_pvar_idx(src, tmp, bgzip)
    except BaseException:
        tmp.unlink(missing_ok=True)
        raise
    os.replace(tmp, out)
    tbi = Path(str(out) + ".tbi")
    tbi.unlink(missing_ok=True)
    proc2 = subprocess.run(
        [tabix, "-p", "vcf", str(out)], capture_output=True, text=True, check=False
    )
    if proc2.returncode != 0:
        raise InputError(f"tabix failed on {out} (is the pvar sorted?): {proc2.stderr.strip()}")
    count = sum(1 for _ in iter_pvar_records(out))
    if count != n:
        raise InputError(f"{out}: wrote {n} records but read back {count}")
    return out, n


def _write_pvar_idx(src: Path, tmp: Path, bgzip: str) -> int:
    """Stream ``src`` into the 9-column indexed layout, bgzip-compressed into ``tmp``."""
    n = 0
    with open(tmp, "wb") as fo:
        proc = subprocess.Popen([bgzip, "-c"], stdin=subprocess.PIPE, stdout=fo)
        assert proc.stdin is not None
        try:
            cols: dict[str, int] | None = None
            buf: list[str] = ["\t".join(PVAR_IDX_HEADER)]
            for line in _decompressed_lines(src):
                s = line.rstrip("\n").rstrip("\r")
                if not s or s.startswith("##"):
                    continue
                f = s.split("\t") if "\t" in s else s.split()
                if cols is None:
                    if s.startswith("#"):
                        names = [x.lstrip("#").upper() for x in f]
                        missing = [
                            c for c in ("CHROM", "POS", "ID", "REF", "ALT") if c not in names
                        ]
                        if missing:
                            raise InputError(f"{src}: header lacks {', '.join(missing)}")
                        cols = {c: names.index(c) for c in names}
                        continue
                    # headerless pvar = .bim order: CHROM ID CM POS ALT REF
                    cols = {"CHROM": 0, "ID": 1, "POS": 3, "ALT": 4, "REF": 5}
                n += 1
                rec = [
                    f[cols["CHROM"]],
                    f[cols["POS"]],
                    f[cols["ID"]],
                    f[cols["REF"]],
                    f[cols["ALT"]],
                ]
                for c in ("QUAL", "FILTER", "INFO"):
                    rec.append(f[cols[c]] if c in cols and cols[c] < len(f) else ".")
                rec.append(str(n))
                buf.append("\t".join(rec))
                if len(buf) >= 10000:
                    proc.stdin.write(("\n".join(buf) + "\n").encode())
                    buf = []
            if buf:
                proc.stdin.write(("\n".join(buf) + "\n").encode())
        finally:
            proc.stdin.close()
            rc = proc.wait()
    if rc != 0:
        raise InputError("bgzip failed")
    if n == 0:
        raise InputError(f"{src}: no variant records")
    return n
