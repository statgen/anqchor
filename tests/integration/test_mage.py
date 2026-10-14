"""Acceptance tests on the in-repo MAGE chr21/22 example (``examples/mage/data``).

Every analysis runs once per session, all at the same time, into ``tmp_path_factory``
(about 25 s in total with qpgentools 0.2.0):

* ``match``: flag-only ``anqchor match run`` (independence scoring, rint, covariates);
* ``match_mah``: the same with ``--mahalanobis --lambda auto``;
* ``cis/<label>``: one flag-only ``anqchor cis run`` per covariate set (nocov, basic,
  peer10, peer60) with the ``pair`` engine, peer60 with ``--rint`` (peer60_rint), nocov,
  peer60 and peer60_rint again with ``rect``, and peer60 on the first 300 reference pairs
  (``subset``);
* ``trans/<label>``: one flag-only ``anqchor trans run`` per setup (nocov, basic, peer10,
  peer60, peer60 with ``--rint``) over the reference's variants x genes.

``anqchor cis compare`` then combines the cis runs and ``anqchor trans compare`` the trans
runs. Expected values live in
``expected_mage.json``. The tests are skipped when no qpgentools binary is found
(see ``tests/conftest.py``).
"""

from __future__ import annotations

import gzip
import json
import os
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import pandas as pd
import pytest
import yaml

from anqchor.config.load import load_config
from anqchor.engine.minimake import stamp_path
from anqchor.pipelines.layout import Layout as L

REPO = Path(__file__).resolve().parents[2]
EXAMPLES = REPO / "examples" / "mage"
EXPECTED: dict[str, Any] = json.loads((Path(__file__).parent / "expected_mage.json").read_text())

DS = "dataset"  # dataset and reference names of flag-only runs
REF = "reference"
CIS_LABELS = ("nocov", "basic", "peer10", "peer60", "peer60_rint")
TRANS_LABELS = CIS_LABELS
LEVELS = ("sign", "nominal", "fdr", "bonferroni", "strict")
integration = pytest.mark.integration


# --------------------------------------------------------------------------- helpers
@dataclass
class Cli:
    """Result of one ``anqchor`` invocation."""

    args: tuple[str, ...]
    returncode: int
    stdout: str
    stderr: str

    @property
    def ran(self) -> list[str]:
        """Steps listed on the ``ran:`` line (empty when nothing ran)."""
        for line in self.stdout.splitlines():
            if line.startswith("ran: "):
                rest = line[len("ran: ") :]
                return [] if rest.startswith("nothing") else rest.split(", ")
        raise AssertionError(f"no 'ran:' line in output:\n{self.stdout}\n{self.stderr}")

    def check(self) -> Cli:
        assert self.returncode == 0, f"anqchor {self.args} failed:\n{self.stdout}\n{self.stderr}"
        return self


def _env(qpgen: Path) -> dict[str, str]:
    env = {**os.environ, "ANQCHOR_QPGEN": str(qpgen)}
    env.pop("VIRTUAL_ENV", None)
    return env


def start(*args: str | Path, qpgen: Path) -> tuple[tuple[str, ...], subprocess.Popen[str]]:
    """Start ``python -m anqchor ARGS`` from the repository root (make calls back into it)."""
    a = tuple(map(str, args))
    proc = subprocess.Popen(
        [sys.executable, "-m", "anqchor", *a],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        env=_env(qpgen),
        cwd=REPO,
    )
    return a, proc


def finish(started: tuple[tuple[str, ...], subprocess.Popen[str]]) -> Cli:
    a, proc = started
    out, err = proc.communicate(timeout=600)
    return Cli(a, proc.returncode, out, err)


def anqchor(*args: str | Path, qpgen: Path, check: bool = True) -> Cli:
    res = finish(start(*args, qpgen=qpgen))
    return res.check() if check else res


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t")


def load_json(path: Path) -> Any:
    return json.loads(path.read_text())


def within(observed: float, expected: float, tol: float) -> bool:
    return abs(observed - expected) <= tol


def inputs(data: Path, *, cov: bool = True, pairs: Path | None = None) -> list[str | Path]:
    """Flag-only inputs of the example: genotypes, expression, reference, gene loci, covariates."""
    args: list[str | Path] = [
        "--geno-list",
        data / "geno" / "mage.pgen_list.tsv",
        "--pheno",
        data / "pheno" / "mage_expr.chr21_22.regenie.tsv.gz",
        "--pairs",
        pairs or data / "reference" / "gtex_lcl_eqtl.chr21_22.tsv.gz",
    ]
    if cov:
        args += ["--cov", data / "pheno" / "mage_covariates.regenie.tsv.gz"]
    return args


@dataclass
class Runs:
    """The session's finished runs, their first invocations and snapshots of their summaries."""

    base: Path
    first: dict[str, Cli] = field(default_factory=dict)
    json: dict[str, Any] = field(default_factory=dict)
    tables: dict[str, pd.DataFrame] = field(default_factory=dict)

    @property
    def match(self) -> Path:
        return self.base / "match"

    def trans(self, label: str = "peer60") -> Path:
        return self.base / "trans" / label

    def cis(self, label: str, engine: str = "pair") -> Path:
        return self.base / ("cis" if engine == "pair" else f"cis_{engine}") / label

    def match_args(self, data: Path, *extra: str | Path) -> list[str | Path]:
        return [
            "match",
            "run",
            *inputs(data),
            "--out",
            self.match,
            "-j",
            "4",
            "--threads",
            "4",
            *extra,
        ]


@pytest.fixture(scope="session")
def runs(mage_data: Path, qpgen_bin: Path, tmp_path_factory: pytest.TempPathFactory) -> Runs:
    r = Runs(tmp_path_factory.mktemp("mage"))
    data = mage_data
    sub = r.base / "reference.first300.tsv"
    with gzip.open(data / "reference" / "gtex_lcl_eqtl.chr21_22.tsv.gz", "rt") as fh:
        sub.write_text("".join(line for _, line in zip(range(301), fh, strict=False)))
    cov_sets = {
        "nocov": None,
        "basic": "sex,PC1-5",
        "peer10": "sex,PC1-5,PEER1-10",
        "peer60": "all",
    }

    def cis(
        label: str,
        engine: str,
        pairs: Path | None = None,
        cov: str | None = "all",
        rint: bool = False,
    ) -> list[str | Path]:
        args = [
            "cis",
            "run",
            *inputs(data, cov=cov is not None, pairs=pairs),
            "--trait-annot",
            data / "pheno" / "mage_genes.chr21_22.bed",
            "--assoc-engine",
            engine,
        ]
        if cov not in (None, "all"):
            args += ["--cov-include", cov]
        if rint:
            args.append("--rint")
        return [*args, "--out", r.cis(label, engine)]

    jobs: dict[str, list[str | Path]] = {
        "match": r.match_args(data),
        "match_mah": [
            "match",
            "run",
            *inputs(data),
            "--out",
            r.base / "match_mah",
            "--threads",
            "4",
            "--mahalanobis",
            "--lambda",
            "auto",
        ],
    }
    for label, cov in cov_sets.items():
        jobs[f"cis/{label}"] = cis(label, "pair", cov=cov)
    jobs["cis/peer60_rint"] = cis("peer60_rint", "pair", rint=True)
    for label in ("nocov", "peer60"):
        jobs[f"cis_rect/{label}"] = cis(label, "rect", cov=cov_sets[label])
    jobs["cis_rect/peer60_rint"] = cis("peer60_rint", "rect", rint=True)
    jobs["cis/subset"] = cis("subset", "pair", pairs=sub)
    for label in TRANS_LABELS:
        cov = cov_sets[label.removesuffix("_rint")]
        args = ["trans", "run", *inputs(data, cov=cov is not None)]
        args += ["--trait-annot", data / "pheno" / "mage_genes.chr21_22.bed"]
        if cov not in (None, "all"):
            args += ["--cov-include", cov]
        if label.endswith("_rint"):
            args.append("--rint")
        jobs[f"trans/{label}"] = [*args, "--out", r.trans(label)]
    started = {k: start(*v, qpgen=qpgen_bin) for k, v in jobs.items()}
    r.first = {k: finish(p).check() for k, p in started.items()}

    # snapshots, before the engine tests below touch the match run
    d = r.match
    for key, rel in {
        "geno": L.geno_json,
        "ref": L.ref_json(REF),
        "ds": L.ds_json(DS),
        "mapping": L.dsref_json(DS, REF),
        "match": L.match_summary(DS, "summary.json"),
        "run": L.run_json,
    }.items():
        r.json[key] = load_json(d / rel)
    r.json["match_mah"] = load_json(r.base / "match_mah" / L.match_summary(DS, "summary.json"))
    r.json["trans"] = {
        label: load_json(r.trans(label) / L.trans_summary(DS, "summary.json"))
        for label in TRANS_LABELS
    }
    r.tables["samples"] = read_tsv(d / L.match_summary(DS, "samples.tsv"))
    return r


@pytest.fixture(scope="session")
def compare(runs: Runs, qpgen_bin: Path) -> Path:
    """``anqchor cis compare`` of the four pair-engine runs (labels = directory names)."""
    out = runs.base / "cmp"
    res = anqchor(
        "cis", "compare", *(runs.cis(lab) for lab in CIS_LABELS), "--out", out, qpgen=qpgen_bin
    )
    assert "warning" not in res.stdout + res.stderr  # same reference, same pairs
    return out


# ------------------------------------------------------------------ shipped configs
@pytest.mark.parametrize("name", ["match.yaml", "trans.yaml"])
def test_example_config_is_valid(name: str) -> None:
    """The shipped example configuration validates and its input files exist (no qpgen needed)."""
    path = EXAMPLES / name
    if not path.exists():
        pytest.skip(f"{path} not shipped")
    cfg = load_config(path)
    pipeline = name.split(".")[0]
    assert pipeline in cfg.enabled_pipelines()
    files = [cfg.genotypes.pgen_list] + [r.pairs for r in cfg.references.values()]
    for ds in cfg.datasets.values():
        files += [ds.pheno, ds.cov, ds.trait_annot]
    for f in files:
        assert f is None or Path(f).exists(), f


# ------------------------------------------------------------------------- prepare
@integration
def test_prepare(runs: Runs) -> None:
    e = EXPECTED["prepare"]
    j = runs.json
    assert j["geno"]["n_samples"] == e["n_geno_samples"]
    assert j["geno"]["chrom_style"] == e["chrom_style"]
    ds = j["ds"]
    assert ds["n_samples_pheno"] == e["n_pheno_samples"]
    assert ds["n_samples_mapped"] == e["n_samples_mapped"]
    assert ds["n_pheno_unmapped"] == e["n_pheno_unmapped"]
    assert ds["n_traits"] == e["n_traits"]
    assert ds["n_cov_columns"] == e["n_cov_columns"]
    assert ds["cov_categorical"] == e["categorical"]  # sex XX/XY recoded without being asked
    for k, v in e["reference"].items():
        assert j["ref"][k] == v, k
    for k, v in e["mapping"].items():
        assert j["mapping"][k] == v, k
    # the covariate matrix handed to qpgen: regenie layout, 66 numeric columns
    with gzip.open(runs.match / L.ds_cov(DS), "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
    assert header[:2] == ["FID", "IID"] and len(header) == 2 + e["n_cov_columns"]


# --------------------------------------------------------------------------- match
@integration
def test_match_independent(runs: Runs) -> None:
    e = EXPECTED["match"]
    s = runs.json["match"]
    assert s["mahalanobis"] is False and s["rint"] is True and s["use_cov"] is True
    assert s["group_counts"] == e["group_counts"]
    assert s["verdict"] == e["verdict"]
    counts = runs.tables["samples"]["qpgen_status"].value_counts().to_dict()
    assert counts == e["status_counts"]
    assert {k: v for k, v in s["status_counts"].items() if v} == e["status_counts"]
    assert s["failing"] == [] and s["n_fail"] == 0
    # NA11919: own genotype at background level, best candidate without a clear lead
    u = e["unclear"]
    smp = runs.tables["samples"].set_index("pheno_id")
    row = smp.loc[u["pheno_id"]]
    assert row["qpgen_status"] == "UNCLEAR" and row["group"] == "WARN"
    assert row["best_id"] == u["best_id"] and row["rank_self"] >= u["min_rank_self"]
    assert row["z_gap"] < s["z_diff"] and pd.isna(row["proposed_geno_id"])  # no new ID
    assert s["n_duplicates"] == 0 and s["warnings"] == []
    assert s["n_status_disagreements"] == 0
    with gzip.open(runs.match / L.prs_matrix(REF), "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        n_rows = sum(1 for _ in fh)
    assert [n_rows, len(header) - 2] == e["prs_shape"]
    with gzip.open(runs.match / L.match_qpgen_outputs(DS)[2], "rt") as fh:
        assert sum(1 for _ in fh) - 1 == e["n_match_all_rows"]


@integration
def test_match_mahalanobis_auto_lambda(runs: Runs) -> None:
    e = EXPECTED["match"]["mahalanobis_auto_lambda"]
    s = runs.json["match_mah"]
    assert s["mahalanobis"] is True and s["lambda_mode"] == "auto"
    assert within(s["lambda_used"], e["lambda"], e["lambda_tolerance"])
    assert {k: v for k, v in s["status_counts"].items() if v} == e["status_counts"]
    assert s["failing"] == []
    samples = pd.read_csv(runs.base / "match_mah" / L.match_summary(DS, "samples.tsv"), sep="\t")
    assert sorted(samples.loc[samples["qpgen_status"] == "UNCLEAR", "pheno_id"]) == e["unclear"]


@integration
def test_match_run_writes_no_html(runs: Runs) -> None:
    assert not (runs.match / "report").exists()
    assert "match.dataset.summarize" in runs.first["match"].ran
    assert not any(n.startswith("report.") for n in runs.first["match"].ran)


# ----------------------------------------------------------------------------- cis
def check_counts(what: str, row: pd.Series, levels: dict[str, int], suffix: str = "") -> None:
    tol = EXPECTED["cis"]["tolerance"]
    for level, v in levels.items():
        got = int(row[f"n_rep_{level}{suffix}"])
        assert within(got, v, tol), f"{what} {level}{suffix}: {got} vs {v} ± {tol}"


@integration
@pytest.mark.parametrize("label", CIS_LABELS)
def test_cis_run(runs: Runs, label: str) -> None:
    """Each run evaluates one input set and one phenotype transform; rows carry the run
    directory name as label."""
    e = EXPECTED["cis"]
    exp = e["runs"][label]
    d = runs.cis(label)
    assert not (d / "report").exists()  # 'cis run' stops at the summaries
    assert not any(n.startswith("report.") for n in runs.first[f"cis/{label}"].ran)
    st = read_tsv(d / L.cis_summary(DS, "strategies.tsv")).set_index("label")
    assert list(st.index) == [label]
    assert "transform" not in st.columns
    row = st.loc[label]
    assert row["n_cov"] == exp["n_cov"] and bool(row["rint"]) == exp["rint"]
    for k, v in e["status"].items():
        assert row[k] == v, k
    check_counts(label, row, exp["counts"])
    s = load_json(d / L.cis_summary(DS, "summary.json"))
    assert s["engine"] == "pair" and s["n_variant_absent"] == 0 and s["rint"] is exp["rint"]
    long = read_tsv(d / L.cis_summary(DS, "pairs.long.tsv.gz"))
    assert len(long) == e["status"]["n_pairs"] and set(long["run"]) == {label}
    tested = long[long["status"] == "tested"]
    # z is aligned with beta and with the reference z (beta_ref / se_ref)
    assert ((tested["z"] > 0) == (tested["beta"] > 0)).all()
    assert tested["z_ref"].notna().all() and tested["log10p_ref"].notna().all()
    assert ((tested["z"] * tested["z_ref"] > 0) == (tested["sign_concordant"] == 1)).all()


@integration
def test_cis_engines_agree(runs: Runs) -> None:
    """pair-assoc and rect-assoc residualise phenotype and genotypes on the same covariates
    (and apply RINT at the same point), so they give identical replication counts; rect adds
    lambda_GC of its local controls."""
    e = EXPECTED["cis"]["rect_lambda_gc_control"]
    for label in ("nocov", "peer60", "peer60_rint"):
        a = read_tsv(runs.cis(label) / L.cis_summary(DS, "strategies.tsv")).set_index("label")
        b = read_tsv(runs.cis(label, "rect") / L.cis_summary(DS, "strategies.tsv")).set_index(
            "label"
        )
        cols = [f"n_rep_{lv}" for lv in LEVELS]
        assert a[cols].equals(b[cols]), label
        lam = float(b.loc[label, "lambda_gc_control"])
        if label == "nocov":
            assert lam >= e["nocov_min"]
        else:
            lo, hi = e["covariates_bounds"]
            assert lo <= lam <= hi


@integration
def test_cis_compare_tables(runs: Runs, compare: Path) -> None:
    e = EXPECTED["cis"]
    for name in (
        "runs.tsv",
        "strategies.tsv",
        "pairs.long.tsv.gz",
        "traits.tsv",
        "power_curve.tsv",
        "summary.json",
        "compare.json",
        "report/cis.html",
    ):
        assert (compare / name).exists(), name
    rt = read_tsv(compare / "runs.tsv").set_index("label")
    assert list(rt.index) == list(CIS_LABELS)  # labels default to the run directory names
    for label in CIS_LABELS:
        exp = e["runs"][label]
        assert rt.loc[label, "n_cov"] == exp["n_cov"]
        assert bool(rt.loc[label, "rint"]) == exp["rint"]
        got = rt.loc[label, "cov_include"]
        assert pd.isna(got) if exp["cov_include"] is None else got == exp["cov_include"], label
        assert Path(rt.loc[label, "run_dir"]) == runs.cis(label)
        settings = compare / "settings" / f"{label}.anqchor.yaml"
        assert yaml.safe_load(settings.read_text()) == yaml.safe_load(
            (runs.cis(label) / L.resolved_config).read_text()
        )
    st = read_tsv(compare / "strategies.tsv").set_index("label")
    assert list(st.index) == list(CIS_LABELS)
    assert (st["n_common_tested"] == e["status"]["n_tested"]).all()
    for label in CIS_LABELS:
        row = st.loc[label]
        check_counts(label, row, e["runs"][label]["counts"])
        check_counts(
            label, row, e["runs"][label]["counts"], "_common"
        )  # every pair tested in every run
    # effective sample size relative to the baseline (nocov): covariates add power
    rn = e["rel_neff"]
    for label in CIS_LABELS:
        assert within(st.loc[label, "rel_neff"], rn[label], rn["tolerance"]), label
    assert st.loc["nocov", "rel_neff"] == pytest.approx(1.0)
    s = load_json(compare / "summary.json")
    assert s["comparability"]["same_reference"] is True
    assert s["comparability"]["rint"] == {lab: e["runs"][lab]["rint"] for lab in CIS_LABELS}
    row = s["per_run"]["peer60"]
    assert type(row["n_common_tested"]) is int and type(row["n_rep_strict_common"]) is int
    assert s["baseline"] == "nocov"
    long = pd.read_csv(compare / "pairs.long.tsv.gz", sep="\t")
    assert len(long) == len(CIS_LABELS) * e["status"]["n_pairs"]
    assert list(dict.fromkeys(long["run"])) == list(CIS_LABELS)
    page = (compare / "report" / "cis.html").read_text(encoding="utf-8")
    assert (
        f"Comparing {len(CIS_LABELS)} runs" in page
        and "Effective sample size relative to nocov" in page
    )


@integration
def test_cis_compare_common_pairs(runs: Runs, qpgen_bin: Path, tmp_path: Path) -> None:
    """Runs on different reference pairs: the *_common columns recount on the shared pairs."""
    e = EXPECTED["cis"]["subset"]
    out = tmp_path / "cmp"
    res = anqchor(
        "cis",
        "compare",
        runs.cis("peer60"),
        runs.cis("subset"),
        runs.cis("nocov"),
        "--out",
        out,
        "--no-report",
        qpgen=qpgen_bin,
    )
    assert "different reference pairs" in res.stderr and "*_common columns" in res.stderr
    assert not (out / "report").exists()
    st = read_tsv(out / "strategies.tsv").set_index("label")
    assert (st["n_common_tested"] == e["n_pairs"]).all()
    assert st.loc["subset", "n_tested"] == e["n_pairs"]
    assert st.loc["peer60", "n_tested"] == EXPECTED["cis"]["status"]["n_tested"]
    # same samples, covariates and pairs: the full run restricted to the subset's pairs
    # gives the subset run's counts exactly
    for lv in LEVELS:
        assert st.loc["peer60", f"n_rep_{lv}_common"] == st.loc["subset", f"n_rep_{lv}"], lv
    assert st.loc["subset", "n_rep_strict"] == e["n_rep_strict_raw"]
    assert within(
        st.loc["nocov", "n_rep_strict_common"],
        e["nocov_strict_common_raw"],
        EXPECTED["cis"]["tolerance"],
    )
    assert not load_json(out / "summary.json")["comparability"]["same_reference"]


@integration
def test_cis_compare_labels(runs: Runs, qpgen_bin: Path, tmp_path: Path) -> None:
    out = tmp_path / "cmp"
    anqchor(
        "cis",
        "compare",
        runs.cis("nocov"),
        runs.cis("peer60"),
        "--out",
        out,
        "--label",
        "A",
        "--label",
        "B",
        "--no-report",
        qpgen=qpgen_bin,
    )
    assert list(read_tsv(out / "runs.tsv")["label"]) == ["A", "B"]
    assert list(read_tsv(out / "strategies.tsv")["label"]) == ["A", "B"]
    assert sorted(p.name for p in (out / "settings").iterdir()) == [
        "A.anqchor.yaml",
        "B.anqchor.yaml",
    ]
    bad = anqchor(
        "cis",
        "compare",
        runs.cis("nocov"),
        "--out",
        tmp_path / "x",
        "--label",
        "A",
        "--label",
        "B",
        qpgen=qpgen_bin,
        check=False,
    )
    assert bad.returncode == 2


@integration
def test_cis_report_rerenders(runs: Runs, compare: Path, qpgen_bin: Path) -> None:
    """'cis report' re-renders a comparison or a single run from its tables only."""
    html = compare / "report" / "cis.html"
    before = html.stat().st_mtime_ns
    res = anqchor("cis", "report", "--out", compare, qpgen=qpgen_bin)
    assert str(html) in res.stdout and html.stat().st_mtime_ns > before
    single = runs.cis("peer10")
    res = anqchor("cis", "report", "--out", single, qpgen=qpgen_bin)
    page = single / L.report("cis")
    assert page.read_text(encoding="utf-8").lstrip().lower().startswith("<!doctype html>")
    assert not (single / L.report("index")).exists()


@integration
def test_cis_settings_reproduce_run(
    runs: Runs, compare: Path, qpgen_bin: Path, tmp_path: Path
) -> None:
    """A run's saved settings (settings/<label>.anqchor.yaml = its anqchor.resolved.yaml)
    reproduce it exactly in a new directory."""
    restored = tmp_path / "basic"  # same directory name: the rows are labelled by it
    anqchor(
        "cis",
        "run",
        "--config",
        compare / "settings" / "basic.anqchor.yaml",
        "--out",
        restored,
        qpgen=qpgen_bin,
    )
    rel = L.cis_summary(DS, "pairs.long.tsv.gz")
    a = gzip.open(runs.cis("basic") / rel, "rt").read()
    b = gzip.open(restored / rel, "rt").read()
    assert a == b


# --------------------------------------------------------------------------- trans
@integration
def test_trans(runs: Runs) -> None:
    """One run per setup: the rectangle, replication by class, hotspots and calibration."""
    e = EXPECTED["trans"]
    lo, hi = e["covariate_lambda_gc_bounds"]
    for label, exp in e["runs"].items():
        t = runs.json["trans"][label]
        assert t["run"] == label and t["rint"] is label.endswith("_rint"), label
        assert (t["n_variants_tested"], t["n_traits_tested"]) == (e["n_variants"], e["n_traits"])
        assert t["cells"] == e["cells"], label
        rep = t["replication"]
        assert rep["trans"]["n_pairs"] == 0  # the GTEx pairs are all cis
        got = [rep["cis"][k] for k in ("n_tested", "n_rep_sign", "n_rep_nominal", "n_rep_strict")]
        assert all(
            within(g, w, e["known_tolerance"]) for g, w in zip(got, exp["cis"], strict=True)
        ), (
            label,
            got,
        )
        inf = t["inflation"]
        assert inf["n_cells"] == e["cells"]["cross_chrom"] and inf["n_cells_excluded"] == 0
        assert within(inf["lambda_gc"], exp["lambda_gc"], e["lambda_gc_tolerance"]), (
            label,
            inf["lambda_gc"],
        )
        h = t["hotspots"]
        assert h["default"]["n_hotspots"] == exp["hotspots_p1e-5_k5"], label
        keys = [x["key"] for x in h["thresholds"]]
        assert h["counts"][keys.index("bonferroni")] == exp["hotspots_bonferroni"], label
        if label != "nocov":
            assert lo <= inf["lambda_gc"] <= hi, label
    assert any("inflated" in w for w in runs.json["trans"]["nocov"]["warnings"])


@integration
def test_trans_run_writes_no_html_and_report_renders_it(
    runs: Runs, mage_data: Path, qpgen_bin: Path
) -> None:
    """'trans run' stops at the summaries; 'trans report' renders the HTML without qpgen;
    'trans all' on a finished run only adds the report."""
    d = runs.trans("basic")
    assert "report.trans" not in runs.first["trans/basic"].ran
    assert not (d / L.report("trans")).exists()
    res = anqchor("trans", "report", "--out", d, qpgen=qpgen_bin)
    html = d / L.report("trans")
    assert str(html) in res.stdout and "trans-QTL evaluation with anQChor" in html.read_text()
    html.unlink()
    data = mage_data
    args = [
        "trans",
        "all",
        *inputs(data),
        "--trait-annot",
        data / "pheno" / "mage_genes.chr21_22.bed",
    ]
    res = anqchor(*args, "--cov-include", "sex,PC1-5", "--out", d, qpgen=qpgen_bin)
    assert res.ran == ["report.trans"] and html.exists()
    assert not (d / L.report("index")).exists()


@integration
def test_trans_compare(runs: Runs, qpgen_bin: Path, tmp_path: Path) -> None:
    """'trans compare' of the five trans runs: one row per run with the runs' own numbers,
    no comparability warning (same rectangle), and 'trans report' re-renders the page."""
    e = EXPECTED["trans"]
    out = tmp_path / "trans_cmp"
    res = anqchor(
        "trans",
        "compare",
        *(runs.trans(lab) for lab in TRANS_LABELS),
        "--out",
        out,
        qpgen=qpgen_bin,
    )
    assert f"compared 5 run(s): {', '.join(TRANS_LABELS)}" in res.stdout
    assert "different" not in res.stderr  # the only warning: nocov's inflation
    t = read_tsv(out / L.compare_table("runs.tsv")).set_index("label")
    assert list(t.index) == list(TRANS_LABELS)
    for label, exp in e["runs"].items():
        r = t.loc[label]
        got = [r[f"cis_n_{k}"] for k in ("tested", "rep_sign", "rep_nominal", "rep_strict")]
        assert all(within(g, w, e["known_tolerance"]) for g, w in zip(got, exp["cis"], strict=True))
        assert r["cis_n_rep_strict_common"] == r["cis_n_rep_strict"] and r["trans_n_pairs"] == 0
        assert within(r["lambda_gc"], exp["lambda_gc"], e["lambda_gc_tolerance"]), label
        assert r["observed_p1e-5"] == exp["observed_p_lt_1e-5"], label
        assert r["n_hotspots"] == exp["hotspots_p1e-5_k5"], label
        assert r["n_cells"] == sum(e["cells"].values()) and bool(r["rint"]) is label.endswith(
            "_rint"
        )
    grid = read_tsv(out / L.compare_table("hotspot_grid.tsv"))
    bonf = grid[(grid["run"] == "nocov") & (grid["threshold"] == "bonferroni")]
    assert list(bonf["n_hotspots"]) == e["runs"]["nocov"]["hotspots_bonferroni"]
    s = load_json(out / L.compare_table("summary.json"))
    assert s["comparability"]["same_rectangle"] and s["comparability"]["same_reference"]
    html = out / L.report("trans")
    assert "Comparing 5 runs: trans-QTL evaluation with anQChor" in html.read_text()
    for lab in TRANS_LABELS:
        assert (out / L.compare_settings(lab)).exists()
    html.unlink()
    res = anqchor("trans", "report", "--out", out, qpgen=qpgen_bin)
    assert str(html) in res.stdout and html.exists()


# ----------------------------------------------------------------- engine behaviour
@integration
def test_run_json_and_resolved_config(runs: Runs) -> None:
    for d in (runs.match, runs.cis("nocov"), runs.trans()):
        assert (d / L.resolved_config).exists() and (d / L.run_json).exists()
    secs = runs.json["run"]["step_seconds"]
    assert set(runs.first["match"].ran) <= set(secs)
    assert all(v >= 0 for v in secs.values())


@integration
def test_no_index_html(runs: Runs) -> None:
    """A single-analysis run never writes report/index.html."""
    for d in [
        runs.match,
        *(runs.trans(lab) for lab in TRANS_LABELS),
        *(runs.cis(lab) for lab in CIS_LABELS),
    ]:
        assert not (d / L.report("index")).exists(), d
    assert not any("report.index" in c.ran for c in runs.first.values())


# These run in file order on the match run directory; each leaves it up to date.
@integration
def test_second_invocation_is_noop(runs: Runs, mage_data: Path, qpgen_bin: Path) -> None:
    assert anqchor(*runs.match_args(mage_data), qpgen=qpgen_bin).ran == []
    # the run's resolved settings describe the same run
    again = anqchor(
        "match",
        "run",
        "--config",
        runs.match / L.resolved_config,
        "--out",
        runs.match,
        qpgen=qpgen_bin,
    )
    assert again.ran == []


@integration
def test_steps_all_done(runs: Runs, qpgen_bin: Path) -> None:
    res = anqchor(
        "steps", "--config", runs.match / L.resolved_config, "-p", "match", qpgen=qpgen_bin
    )
    rows = [
        ln.split() for ln in res.stdout.splitlines() if ln and not ln.startswith("run directory")
    ]
    status = {r[0]: r[1] for r in rows if len(r) >= 2}
    assert set(runs.first["match"].ran) <= set(status)
    assert {status[n] for n in runs.first["match"].ran} == {"done"}


@integration
def test_parameter_change_reruns_only_dependent_steps(
    runs: Runs, mage_data: Path, qpgen_bin: Path
) -> None:
    expected = {"match.dataset.qpgen", "match.dataset.summarize"}
    changed = anqchor(*runs.match_args(mage_data, "--z-threshold", "3"), qpgen=qpgen_bin)
    assert set(changed.ran) == expected  # PRS and the prepared inputs are reused
    assert load_json(runs.match / L.match_summary(DS, "summary.json"))["z_threshold"] == 3.0
    assert set(anqchor(*runs.match_args(mage_data), qpgen=qpgen_bin).ran) == expected  # and back
    assert anqchor(*runs.match_args(mage_data), qpgen=qpgen_bin).ran == []
    # changing only the thread count re-runs nothing
    assert anqchor(*runs.match_args(mage_data)[:-2], "--threads", "2", qpgen=qpgen_bin).ran == []


@integration
def test_dry_run_changes_nothing(runs: Runs, mage_data: Path, qpgen_bin: Path) -> None:
    stamp = runs.match / stamp_path("match.dataset.summarize")
    stamp.unlink()
    summary = runs.match / L.match_summary(DS, "summary.json")
    before = summary.stat().st_mtime_ns
    res = anqchor(*runs.match_args(mage_data, "--dry-run"), qpgen=qpgen_bin)
    assert "dry run: would run: match.dataset.summarize" in res.stdout
    assert not stamp.exists() and summary.stat().st_mtime_ns == before
    assert anqchor(*runs.match_args(mage_data), qpgen=qpgen_bin).ran == ["match.dataset.summarize"]


@integration
def test_match_report_never_runs_qpgen(runs: Runs, tmp_path: Path) -> None:
    """'match report' renders report/match.html from the run's tables and touches nothing else,
    even when no qpgentools binary can be found."""

    def snapshot() -> dict[Path, int]:
        return {
            p: p.stat().st_mtime_ns
            for p in runs.match.rglob("*")
            if p.is_file() and p.parts[-2] != "report" and p.name != ".anqchor.lock"
        }

    before = snapshot()
    env = {
        **os.environ,
        "ANQCHOR_QPGEN": str(tmp_path / "no-such-qpgentools"),
        "PATH": "/usr/bin:/bin",
    }
    env.pop("VIRTUAL_ENV", None)
    proc = subprocess.run(
        [sys.executable, "-m", "anqchor", "match", "report", "--out", str(runs.match)],
        capture_output=True,
        text=True,
        env=env,
        cwd=REPO,
        check=False,
        timeout=120,
    )
    assert proc.returncode == 0, proc.stdout + proc.stderr
    html = runs.match / L.report("match")
    assert str(html) in proc.stdout
    assert html.read_text(encoding="utf-8").lstrip().lower().startswith("<!doctype html>")
    assert "NA11919" in html.read_text(encoding="utf-8")
    assert snapshot() == before
    assert not (runs.match / L.report("index")).exists()


@integration
def test_match_all_adds_only_the_report(runs: Runs, mage_data: Path, qpgen_bin: Path) -> None:
    (runs.match / L.report("match")).unlink(missing_ok=True)
    args = runs.match_args(mage_data)
    res = anqchor("match", "all", *args[2:], qpgen=qpgen_bin)
    assert res.ran == ["report.match"]
    assert (runs.match / L.report("match")).exists()
    assert not (runs.match / L.report("index")).exists()


@integration
def test_only_and_until(runs: Runs, mage_data: Path, qpgen_bin: Path) -> None:
    def run(*extra: str) -> Cli:
        return anqchor(*runs.match_args(mage_data, *extra), qpgen=qpgen_bin, check=False)

    assert run("--only", "match.dataset.summarize").check().ran == ["match.dataset.summarize"]
    assert run("--until", "match.dataset.summarize").check().ran == []
    assert run("--only", "no.such.step").returncode == 2
