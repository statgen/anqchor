# Python API

anqchor is a command-line tool first; the Python modules are documented here for people who
extend it. Everything below lives under `src/anqchor/`.

## Loading a configuration

```python
from anqchor.config import load_config

cfg = load_config("examples/mage/match.yaml", overrides={"pipelines": {"match": {"z_diff": 2.5}}})
cfg.enabled_pipelines()  # ['match']
cfg.reference_for("match")  # 'gtex_lcl'
cfg.datasets_for("match")  # ['rnaseq']
cfg.pipelines.match.z_diff  # 2.5
```

`load_config` merges overrides over the file, validates with the pydantic models in
`anqchor.config.model`, and returns a `Config` with absolute paths. It raises
`anqchor.config.load.ConfigError` (exit code 2 in the CLI).

## `Layout`: file names

`anqchor.pipelines.layout.Layout` is the single source of truth for every path in a run
directory (all relative to it). Producers and consumers must use it instead of formatting
paths:

```python
from anqchor.pipelines import Layout

Layout.pgen_list  # prepare/genotypes/pgen_list.tsv
Layout.ds_pheno("rnaseq")  # prepare/datasets/rnaseq/pheno.regenie.tsv.gz
Layout.dsref_pairs("rnaseq", "gtex_lcl")  # prepare/datasets/rnaseq/ref.gtex_lcl/pairs.tsv.gz
Layout.cis_summary("rnaseq", "strategies.tsv")
Layout.rect_assoc("trans", "rnaseq", "peer60", "raw")
Layout.report("cis")  # report/cis.html
```

## `Step` and `StepGraph`

```python
from anqchor.engine.step import Step

Step(
    name="cis.rnaseq.summarize",  # unique, dotted; also a make target
    kind="summarize_cis",  # registered Python implementation (or None)
    inputs=(...),  # files read (relative to the run directory)
    outputs=(...),  # files written; all must exist after the step
    commands=(),  # shell lines (for qpgen steps)
    params={...},  # JSON-serialisable; part of the signature
    pipeline="cis",
    description="replication of gtex_lcl pairs on rnaseq",
)
```

- Steps are linked by **paths**: when an input of one step is an output of another, the
  Makefile makes the consumer depend on the producer.
- `params` and `commands` form the step **signature** (`steps/<name>.json`). Anything that
  changes the result must be in them, so that a change re-runs the step and everything
  downstream.
- `StepGraph` de-duplicates identical steps (same name and signature), which is how pipelines
  share `prepare.*` and `prs.*`; a conflicting redefinition or an output produced by two steps
  is an error.

A pipeline module (`anqchor.pipelines.match|cis|trans`) exposes
`build(ctx: BuildContext) -> list[Step]` returning its steps plus the shared ones it needs;
`anqchor.pipelines.build_graph` assembles the graph and the make groups; it drops the
`report.*` steps for `run` commands (`reports=False`) and adds `report.index` only when several
pipelines run together.

## Python step implementations

A Python step's command is generated: `python -m anqchor _step <name> --run-dir .`. The runner
reads `steps/<name>.json` and calls the function registered for its `kind`:

```python
from anqchor.steps import StepContext, StepError, register
from anqchor.util.files import write_json


@register("my_summary")
def my_summary(ctx: StepContext) -> None:
    ds = ctx.params["dataset"]
    table = ctx.path(ctx.inputs[0])            # absolute path of a declared input
    ...
    if bad:
        raise StepError("actionable message")  # reported without a traceback, exit 3
    write_json(ctx.path(ctx.outputs[0]), {"schema_version": 1, ...})
```

Rules:

- Read only declared `inputs` and `params`; never read the config file.
- Write every declared output, atomically (`anqchor.util.files.atomic_write`, `write_json`,
  the `write_*` helpers in `anqchor.io`).
- Turn data problems into `StepError` (the `input_errors` decorator in `anqchor.steps.prepare`
  converts `anqchor.io.InputError`).

## Adding a step kind

1. Implement the function with `@register("<kind>")` in a module under `anqchor/steps/` or
   `anqchor/report/`, and make sure the module is listed in `anqchor.steps.STEP_MODULES`.
2. Add the output file names to `Layout`.
3. Emit the `Step` from the pipeline's `build()` with all parameters in `params`.
4. Add unit tests with small fixtures (see [Development](../development.md)).

## qpgen command builders and readers

```python
from anqchor.qpgen import commands as qc
from anqchor.qpgen.outputs import read_pair_assoc, iter_rect_assoc

argv = qc.rect_assoc(
    "qpgentools",
    genotypes=qc.Genotypes(pgen_list=Layout.pgen_list),
    pheno="pheno.tsv.gz",
    var_list="variants.txt",
    out="out.tsv.gz",
    min_maf=0.01,
)
qc.to_shell(argv)
for chunk in iter_rect_assoc("out.tsv.gz", ("BETA", "LOG10P")):
    chunk.meta, chunk.stats["LOG10P"]  # rows × traits arrays
```

`anqchor.qpgen.commands` is the only module that knows qpgen flag names, and
`anqchor.qpgen.outputs` the only one that knows qpgen output layouts.

## Comparing cis runs

`anqchor cis compare` is a thin wrapper around `anqchor.steps.cis_compare.compare_runs`:

```python
from pathlib import Path
from anqchor.steps.cis_compare import compare_runs, is_comparison_dir

summary = compare_runs(
    [Path("results/cis/nocov"), Path("results/cis/peer10")],
    Path("results/cis/compare"),
    labels=["nocov", "peer10"],
)
summary["recommended"]["label"]  # e.g. 'peer10/rint'
is_comparison_dir(Path("results/cis/compare"))  # True: compare.json exists
```

It reads each run's `anqchor.resolved.yaml` and `cis/<ds>/` summaries and writes the tables
described in [Output files](outputs.md#cis-compare); it does not render HTML (the CLI calls
`anqchor.cli.execute.render_cis_comparison` for that). Data problems raise
`anqchor.steps.StepError`.

## Comparing trans runs

`anqchor trans compare` wraps `anqchor.steps.trans_compare.compare_runs` the same way:

```python
from pathlib import Path
from anqchor.steps.compare import comparison_kind
from anqchor.steps.trans_compare import compare_runs

summary = compare_runs(
    [Path("results/trans/nocov"), Path("results/trans/peer10")], Path("results/trans/compare")
)
summary["comparability"]["same_cells"]  # True when the runs tested the same cells
comparison_kind(Path("results/trans/compare"))  # 'trans' (compare.json says so)
```

It writes the tables described in [Output files](outputs.md#trans-compare) and no HTML
(the CLI calls `anqchor.cli.execute.render_trans_comparison`).
