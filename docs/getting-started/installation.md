# Installation

anqchor is a Python package (Python ≥ 3.10) plus the qpgen C++ toolkit, which is vendored as
a git submodule (`extern/qpgen`, branch `dev`) and built from the checkout. anqchor is not on
PyPI yet; install it from a clone.

## 1. Clone with submodules

```bash
git clone --recursive https://github.com/statgen/anqchor.git
cd anqchor
# already cloned without --recursive?
git submodule update --init --recursive
```

The second command also fetches qpgen's own submodules (htslib, qgenlib, Eigen).

## 2. Install the Python package

=== "uv (recommended)"

    ```bash
    uv sync                 # creates .venv with anqchor and the dev tools
    uv run anqchor --version
    ```

    Run every command as `uv run anqchor ...`, or activate `.venv`.

=== "pip"

    ```bash
    python -m venv .venv && source .venv/bin/activate
    pip install -e .
    anqchor --version
    ```

Runtime dependencies: `typer`, `pydantic` (v2), `pyyaml`, `pandas`, `numpy`, `scipy`, `plotly`,
`jinja2`. anqchor also needs **GNU make** (3.81 or later; the macOS default works) on `PATH`, or
in `$ANQCHOR_MAKE`.

## 3. Build qpgen

```bash
anqchor qpgen build -j 8
```

This runs qpgen's documented three-stage build inside `extern/qpgen`: htslib (`autoreconf -i`,
`./configure`, `make`), qgenlib (`cmake`, `make`), then qpgen itself (`cmake`, `make`). Stages
whose artefacts already exist are skipped; `--force` rebuilds htslib and qgenlib too. The output
is `extern/qpgen/bin/qpgentools`, and the full log is `extern/qpgen/anqchor-build.log`.

Build requirements: a C++14 compiler with OpenMP, `cmake`, `make`, `autoreconf`, and the
compression libraries htslib uses (`zlib`, `bzip2`, `lzma`; optionally `libdeflate`, `libcurl`,
`libcrypto`). LAPACK/BLAS are optional.

!!! note "macOS"
    Apple Clang has no OpenMP. Install it and Homebrew LLVM before building:

    ```bash
    brew install libomp llvm cmake autoconf automake
    ```

    qpgen's CMake configuration selects Homebrew LLVM automatically. `anqchor qpgen build`
    warns before building when `libomp` or a Homebrew LLVM is missing. Extra CMake settings can
    be passed with `--cmake-arg` (repeatable), for example
    `--cmake-arg=-DCMAKE_CXX_COMPILER=/opt/homebrew/opt/llvm/bin/clang++`.

!!! tip "Using a qpgen binary built elsewhere"
    You do not have to build from the submodule. anqchor looks for `qpgentools` in this order:

    1. `--qpgen PATH` (root option) or `qpgen.binary` in the config;
    2. the `ANQCHOR_QPGEN` environment variable;
    3. `qpgentools` on `PATH`;
    4. `extern/qpgen/bin/qpgentools` of the source checkout.

    An explicit path or `ANQCHOR_QPGEN` that does not name an executable is an error (exit
    code 4); anqchor does not silently fall back to another binary.

## 4. Check the installation

```bash
anqchor qpgen check
```

```text
binary:   /path/to/anqchor/extern/qpgen/bin/qpgentools  (found via submodule)
sha256:   <sha256 of your build>
commit:   5b811e250591931cd80f9a7dfae80d47eb8381e2
status:   ok (pair-prs, match-prs-pheno, pair-assoc, rect-assoc available)
```

The command checks that the four subcommands anqchor uses exist. It prints the binary's
sha256 and the git commit of the checkout the binary sits in. The same information is recorded
in every run's `run.json`.

Then check the qpgen version:

```bash
extern/qpgen/bin/qpgentools --version
# qpgentools 0.2.0 (git 5b811e2)
```

!!! warning "anqchor needs qpgen 0.2.0 or later"
    qpgen 0.2.0 fixed several problems in the association and matching tools that change
    anqchor's numbers. The most important: `pair-assoc` now adjusts the genotypes for the
    covariates as well as the phenotypes, all association tools use `df = N - 2 - k` with k
    covariates, and non-numeric covariate values such as `XX`/`XY` are an error (unless the
    covariate is named with `--cat-cov-names`) instead of being read as 0. Older builds give
    different, and with many covariates optimistic, results.
    anqchor also passes `--cat-cov-names` for categorical covariates such as `sex`, which
    needs qpgen commit `5b811e2` or later (the pinned submodule); an older 0.2.0 build does not
    accept that option, and `anqchor qpgen check` reports it as missing options.
    `anqchor qpgen check` does not read the version, so run `qpgentools --version` yourself.
    Builds before 0.2.0 have no `--version` option at all. Details:
    [qpgen integration](../reference/qpgen.md).

## 5. Run the example

The repository ships a small example dataset in `examples/mage/data/`. The
[quickstart](quickstart.md) runs the quickest analysis on it in about 10 seconds; if that
works, the installation is complete.

## Other tools

- `bgzip` and `tabix` (htslib) are needed by `anqchor prepare index-pvar` to create indexed
  pvar files, and make the reference presence check faster (without `tabix` anqchor reads the
  whole pvar instead). anqchor looks for them in `$ANQCHOR_BGZIP` / `$ANQCHOR_TABIX`, then on `PATH`,
  then in the htslib built inside `extern/qpgen/submodules/htslib`.
- `zstd` on `PATH` is needed only to index a `.pvar.zst`.

## Building this documentation

```bash
uv sync --extra docs
uv run mkdocs serve          # http://127.0.0.1:8000
uv run mkdocs build --strict # static site in site/
```
