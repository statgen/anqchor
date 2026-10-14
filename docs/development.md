# Development

## Repository layout

```text
anqchor/
├── pyproject.toml              # hatchling; console script 'anqchor'; extras: docs; dev group: pytest ruff mypy
├── uv.lock
├── README.md  CHANGELOG.md  DESIGN.md  CLAUDE.md  QPGEN_KNOWN_ISSUES.md  LICENSE
├── mkdocs.yml  docs/           # this site
├── extern/qpgen/               # git submodule, branch dev
├── examples/mage/              # data/ (32 MB example), match.yaml, cis.yaml, trans.yaml, make_compact_data.sh, README
├── src/anqchor/
│   ├── cli/                    # main.py (root app; match run|report|all, cis run|compare|report|all,
│   │                           # trans run|compare|report|all, run/steps/report, hidden _step),
│   │                           # options.py (shared options, flag -> config overrides), execute.py,
│   │                           # qpgen.py (build/check), prepare.py (index-pvar, pgen-list, pheno, cov)
│   ├── config/                 # model.py (pydantic models), load.py (YAML, merge, path resolution)
│   ├── engine/                 # step.py (Step, StepGraph), minimake.py (Makefile), runner.py (make), provenance.py
│   ├── qpgen/                  # locate.py, build.py, commands.py (argv builders), outputs.py (readers)
│   ├── io/                     # pheno.py, cov.py, pairs.py, genotypes.py, maps.py, tables.py
│   ├── pipelines/              # layout.py (file names), common.py, match.py, cis.py, trans.py, report.py, context.py
│   ├── steps/                  # Python step kinds: prepare.py, strategy.py, match.py, cis.py, trans.py;
│   │                           # cis_compare.py, trans_compare.py ('cis compare', 'trans compare',
│   │                           # not steps), compare.py (labels, settings copies, compare.json)
│   ├── report/                 # report step kinds (match, cis, trans, index), render.py, theme.py, templates/, static/
│   ├── stats.py                # significance predicate shared by cis and trans (p < α; BH q ≤ α)
│   └── util/                   # files.py (atomic writes, gz helpers), log.py
└── tests/
    ├── unit/                   # no qpgen needed
    ├── integration/            # test_mage.py + expected_mage.json
    └── fixtures/               # small real qpgen outputs (qpgen/), io fixtures
```

## Environment

```bash
uv sync                    # Python >= 3.10, dev tools included
uv run anqchor --help
```

!!! tip
    If another virtual environment is active, uv warns that `VIRTUAL_ENV` does not match the
    project environment. Run commands as `env -u VIRTUAL_ENV uv run ...` to use the project's
    `.venv`.

## Tests

```bash
uv run pytest -q                       # everything available
uv run pytest -q -m "not integration"  # unit tests only, no qpgen needed
uv run pytest -q -m integration        # the example data end to end
```

The integration tests run on the example data shipped in `examples/mage/data/`, so they need
no download, only a built `qpgentools` (qpgen 0.2.0). They look for the binary in
`$ANQCHOR_QPGEN`, then `PATH`, then `extern/qpgen/bin/`, and skip when none is found:

```bash
ANQCHOR_QPGEN=/opt/qpgen/bin/qpgentools uv run pytest -q -m integration
```

`test_mage.py` runs every analysis once per session, all at the same time: flag-only
`match run` (also with `--mahalanobis --lambda auto`), one flag-only `cis run` per covariate
set (`nocov`, `basic`, `peer10`, `peer60`; `nocov` and `peer60` again with the `rect` engine,
and `peer60` on a subset of the reference), one flag-only `trans run` per setup
(`nocov`, `basic`, `peer10`, `peer60`, `peer60_rint`), `cis compare` over the cis runs and `trans compare` over the trans runs. It asserts the numbers in
`tests/integration/expected_mage.json` and the engine behaviour (re-runs, `--dry-run`,
`--only`, `--until`, reports), including `cis run` → `cis compare` → `cis report` and `trans run` → `trans compare` →
`trans report`. **When qpgen or the example data change, update `expected_mage.json`, not the
test.** The example data is rebuilt from the full MAGE files with
`examples/mage/make_compact_data.sh`.

## Lint and type checks

```bash
uv run ruff check .
uv run ruff format --check .
uv run mypy                  # strict on anqchor.config, anqchor.engine, anqchor.qpgen
```

## Documentation

```bash
uv sync --extra docs
uv run mkdocs serve
uv run mkdocs build --strict
```

The configuration reference is written from `config/model.py`, the CLI reference from
`--help`, and the output reference from fresh runs on the example data (`pipelines/layout.py`
and the step modules); update them together with the code. The reports are being redesigned,
so the pages carry no report screenshots for now.

## Foundation contracts

These rules keep the pieces independent; reviews check them.

1. **qpgen does the numerics.** No statistical re-implementation of what qpgen computes;
   Python does orchestration, normalisation, summaries and reports.
2. **qpgen flags appear only in `qpgen/commands.py`; qpgen output layouts only in
   `qpgen/outputs.py`.**
3. **File names come only from `pipelines/layout.py`.**
4. **A step reads only its declared inputs and params** and writes every declared output
   atomically. Every parameter that changes a result is in `params` (so it is in the signature
   and re-runs the step).
5. **Reports are rendered from summaries** (TSV/JSON with `schema_version`), never from raw
   qpgen output.
6. **Config is the single source of truth.** Flags build a config; the resolved config is
   written to the run directory. Unknown keys are errors.
7. **anqchor always hands qpgen regenie-layout matrices** and a 4-column pgen list with absolute
   prefixes.
8. The generated Makefile stays GNU make 3.81 compatible (no grouped targets, `.SHELLFLAGS`
   or `$(file …)`).

## Adding a step or a pipeline

- **A step**: add its file names to `Layout`, implement the kind with `@register` (module
  listed in `anqchor.steps.STEP_MODULES`), emit the `Step` from the pipeline's `build()`, add
  unit tests. See [Python API](reference/python-api.md#adding-a-step-kind).
- **A pipeline**: add `pipelines/<name>.py` with `build(ctx) -> list[Step]` and
  `REPORT_STEP`, add it to `PIPELINES`, add a config model under `PipelinesConfig`, `run`,
  `report` and `all` commands in `cli/main.py`, a report kind and the index tile, and
  integration expectations.
- **A config key**: add it to the pydantic model with a docstring, pass it through the step
  `params`, document it in `docs/user-guide/configuration.md`, add a unit test.
- **A qpgen flag**: add it to the builder in `qpgen/commands.py` and its test.

## Release checklist

- `CHANGELOG.md` updated; version in `pyproject.toml` (`__version__` is read from the installed metadata).
- Full test suite green with the pinned `extern/qpgen` commit; `expected_mage.json` current.
- `mkdocs build --strict` clean.
