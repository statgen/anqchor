"""Pydantic models for the anqchor configuration file.

The YAML layout is documented in ``docs/configuration.md``. Every model forbids
unknown keys so that typos fail at load time rather than being silently ignored.
Paths are stored as given; :func:`anqchor.config.load.resolve_paths` makes them
absolute relative to the configuration file before anything uses them.
"""

from __future__ import annotations

import math
import re
from collections.abc import Sequence
from pathlib import Path
from typing import Annotated, Any, Literal

from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator

NAME_PATTERN = r"^[A-Za-z0-9][A-Za-z0-9_.-]*$"
NAME_RE = re.compile(NAME_PATTERN)

PhenoFormat = Literal["regenie", "tensorqtl", "tsv-sample-row", "tsv-sample-col"]
"""Matrix layouts understood by anqchor (and by qpgen)."""

TraitIdRule = Literal["as_is", "before_colon", "strip_version", "before_colon_strip_version"]
"""How reference trait IDs are matched to dataset trait IDs.

The rule is applied to both sides and IDs match when the transformed keys agree:

* ``as_is``: identical strings.
* ``before_colon``: the part before the first ``:`` (``ENSG…5:WASH7P`` → ``ENSG…5``).
* ``strip_version``: drop a trailing ``.<digits>`` (``ENSG…5`` → ``ENSG…``).
* ``before_colon_strip_version``: both of the above.
"""


class StrictModel(BaseModel):
    """Base model: unknown keys are errors and assignment is validated."""

    model_config = ConfigDict(extra="forbid", validate_assignment=True, populate_by_name=True)


FILE_PREFIX = "file:"
"""Optional prefix marking a value as a file path (``traits: file:my_traits.txt``)."""
REFERENCE_PREFIX = "reference:"
"""Prefix naming a reference (``variants: reference:gtex_lcl``)."""


def _check_name(kind: str, name: str) -> str:
    if not NAME_RE.match(name):
        raise ValueError(
            f"invalid {kind} name {name!r}: use letters, digits, '_', '.', '-' "
            "and start with a letter or digit"
        )
    return name


def strip_file_prefix(value: Any) -> Any:
    """``"file:x"`` → ``"x"``; anything else is returned unchanged."""
    if isinstance(value, str) and value.startswith(FILE_PREFIX):
        return value[len(FILE_PREFIX) :]
    return value


class QpgenConfig(StrictModel):
    """Where to find ``qpgentools`` and how many threads it may use."""

    binary: Path | None = None
    """Explicit path to ``qpgentools``; ``None`` means auto-discovery."""
    threads: Annotated[int, Field(ge=1)] = 1
    """Threads for ``match-prs-pheno --threads`` (pair-prs, pair-assoc and rect-assoc are
    single-threaded). Passed to the step through ``$ANQCHOR_THREADS``, not as a step
    parameter, so changing it re-runs nothing."""


class GenotypeConfig(StrictModel):
    """One genotype set: a qpgen pgen list or a single PLINK2 prefix."""

    pgen_list: Path | None = None
    """qpgen ``--pgen-list`` file: ``CHROM BEG END PREFIX`` or ``CHROM BEG END PGEN PSAM PIVAR``."""
    pfile: Path | None = None
    """A single PLINK2 prefix with ``.pgen``, ``.psam`` and ``.pvar.idx.gz`` (+ ``.tbi``)."""

    @model_validator(mode="after")
    def _exactly_one(self) -> GenotypeConfig:
        if (self.pgen_list is None) == (self.pfile is None):
            raise ValueError("genotypes: give exactly one of 'pgen_list' or 'pfile'")
        return self


class ColumnMap(StrictModel):
    """Column names of a reference pair file."""

    trait: str = "trait"
    variant: str = "variant"
    beta: str = "beta"
    se: str = "se"
    pval: str | None = None
    """Optional p-value column, carried into reports. ``None`` reads a column named
    ``pval`` when the file has one."""


class ReferenceConfig(StrictModel):
    """A set of previously known xQTL effects."""

    kind: Literal["cis", "trans"] = "cis"
    pairs: Path | None = None
    """TSV with a header containing trait, variant, beta and se columns (see ``columns``)."""
    variants: Path | None = None
    """Trans references only: a variant list (one ``CHROM:POS:REF:ALT`` per line)."""
    columns: ColumnMap = Field(default_factory=ColumnMap)

    @model_validator(mode="after")
    def _inputs(self) -> ReferenceConfig:
        if self.kind == "cis" and self.pairs is None:
            raise ValueError("a cis reference needs 'pairs'")
        if self.kind == "trans" and self.pairs is None and self.variants is None:
            raise ValueError("a trans reference needs 'pairs' or 'variants'")
        return self


class DatasetConfig(StrictModel):
    """One omics dataset: a phenotype matrix plus optional covariates and annotation."""

    omics: str = "omics"
    """Free-text label shown in reports (``rnaseq``, ``proteomics``, ...)."""
    pheno: Path
    pheno_format: PhenoFormat = "regenie"
    trait_annot: Path | None = None
    """Trait loci: a BED-like file (``chrom start end id`` + anything) or the tensorQTL matrix itself."""
    cov: Path | None = None
    cov_format: PhenoFormat = "regenie"
    cov_categorical: list[str] = Field(default_factory=list)
    """Covariates to treat as categorical even if their values look numeric.

    Non-numeric covariates are detected automatically. Categorical covariates are recoded
    to integer codes because qpgen accepts numeric covariates only: two levels become 1/2,
    and k > 2 levels become k-1 indicator columns coded 1/2."""
    cov_missing: Literal["error", "drop_samples", "impute_mean"] = "error"
    """What to do with missing covariate values; qpgen's association tools need complete matrices."""
    missing_str: list[str] = Field(default_factory=lambda: ["NA"])
    """Strings that mean 'missing' in the phenotype and covariate files."""
    sample_map: Path | None = None
    """Identity matching only: a TSV with header ``geno_id pheno_id``. ``None`` means identical IDs."""
    trait_id_rule: TraitIdRule = "as_is"
    trait_map: Path | None = None
    """Explicit two-column TSV ``ref_trait pheno_trait`` (no header); overrides ``trait_id_rule``.
    ``trait_id_rule: file:<path>`` is accepted as a shorthand for ``trait_map: <path>``."""

    @model_validator(mode="before")
    @classmethod
    def _trait_rule_file(cls, data: Any) -> Any:
        if isinstance(data, dict):
            rule = data.get("trait_id_rule")
            if isinstance(rule, str) and rule.startswith(FILE_PREFIX):
                if data.get("trait_map") not in (None, rule[len(FILE_PREFIX) :]):
                    raise ValueError(
                        "give either 'trait_id_rule: file:<path>' or 'trait_map', not both"
                    )
                data = {**data, "trait_id_rule": "as_is", "trait_map": rule[len(FILE_PREFIX) :]}
        return data


class CovSpec(StrictModel):
    """Strategy covariates: a named subset of the dataset covariates or a separate file."""

    include: list[str] | None = None
    """Covariate names; ranges such as ``PEER1-15`` expand to ``PEER1`` … ``PEER15``.
    The names are honoured exactly: an unknown name is an error."""
    file: Path | None = None
    """A separate covariate file (same format rules as the dataset's ``cov``)."""
    format: PhenoFormat = "regenie"

    @model_validator(mode="after")
    def _exactly_one(self) -> CovSpec:
        if (self.include is None) == (self.file is None):
            raise ValueError("strategy cov: give exactly one of 'include' or 'file'")
        return self


class InputSetConfig(StrictModel):
    """The inputs of one analysis: covariates (and optionally another phenotype matrix).

    Base of :class:`CisStrategyConfig` (cis) and :class:`StrategyConfig` (trans)."""

    name: str
    cov: Literal["default"] | CovSpec | None = "default"
    """``default`` = all dataset covariates, ``null`` = none, or a :class:`CovSpec`."""
    pheno: Path | None = None
    """Optional alternative phenotype matrix (e.g. another normalization) for this strategy."""
    pheno_format: PhenoFormat = "regenie"

    @field_validator("cov", mode="before")
    @classmethod
    def _cov_forms(cls, v: Any) -> Any:
        """``default``, ``null``, a mapping, or a file path (``file:`` prefix optional)."""
        if isinstance(v, Path):
            return {"file": v}
        if isinstance(v, str) and v != "default":
            return {"file": strip_file_prefix(v)}
        return v

    @field_validator("name")
    @classmethod
    def _name(cls, v: str) -> str:
        return _check_name("strategy", v)


class CisStrategyConfig(InputSetConfig):
    """The inputs of a cis run. The phenotype transform is not part of it: a cis run tests
    one transform, set by ``pipelines.cis.rint`` (flag ``--rint``)."""

    @model_validator(mode="before")
    @classmethod
    def _no_transforms(cls, data: Any) -> Any:
        return _reject_transforms("cis", data)


def _reject_transforms(pipeline: str, data: Any) -> Any:
    """A strategy listing ``transforms`` (configurations written before one transform per
    run): an error that says what to write instead."""
    if isinstance(data, dict) and "transforms" in data:
        t = data.get("transforms")
        raise ValueError(
            f"pipelines.{pipeline}.strategies: 'transforms' ({t!r}) is no longer accepted; a {pipeline} run "
            "tests one phenotype transform. Remove the key (phenotypes as given), or set "
            f"'pipelines.{pipeline}.rint: true' (flag --rint) for RINT after covariate adjustment. To compare "
            "raw and RINT, make two runs"
            + (" and compare them with 'anqchor cis compare'" if pipeline == "cis" else "")
        )
    return data


class StrategyConfig(InputSetConfig):
    """The inputs of a trans run (covariates, optionally another phenotype matrix). The
    phenotype transform is not part of it: a trans run tests one transform, set by
    ``pipelines.trans.rint`` (flag ``--rint``)."""

    @model_validator(mode="before")
    @classmethod
    def _no_transforms(cls, data: Any) -> Any:
        return _reject_transforms("trans", data)


class LevelsConfig(StrictModel):
    """Replication levels. Sign concordance is always reported as the weakest level and is a
    prerequisite of every p-value level."""

    nominal: Annotated[float, Field(gt=0, lt=1)] | None = 0.05
    strict: Annotated[float, Field(gt=0, lt=1)] | None = 1e-5
    bonferroni: bool = True
    """Bonferroni at alpha 0.05 over the tested pairs (cis) or cells (trans)."""
    fdr: Annotated[float, Field(gt=0, lt=1)] | None = 0.05
    """Benjamini-Hochberg q-value cutoff."""


def _default_strategies() -> list[StrategyConfig]:
    return [StrategyConfig(name="default")]


def _default_cis_strategies() -> list[CisStrategyConfig]:
    return [CisStrategyConfig(name="default")]


class MatchConfig(StrictModel):
    """Sample identity matching (``qpgentools pair-prs`` + ``match-prs-pheno``)."""

    enabled: bool = True
    datasets: list[str] | None = None
    """Datasets to check; ``None`` = all."""
    reference: str | None = None
    """Reference for the PRS; ``None`` = the only reference (error if several)."""
    prs_samples: Path | None = None
    """One genotype sample ID per line to compute PRS for a subset; ``None`` = all genotyped samples."""
    use_cov: bool = True
    rint: bool = True
    mahalanobis: bool = False
    """Whiten correlated traits (slower). The default scores traits as independent."""
    lambda_: float | Literal["auto"] = Field(default="auto", alias="lambda")
    """With ``mahalanobis``: ``auto`` tunes λ when enough mapped samples exist, else 0.5."""
    min_self_for_auto: Annotated[int, Field(ge=1)] = 100
    z_threshold: float = 1.96
    z_diff: Annotated[float, Field(ge=0)] = 2.0
    min_weight: float = 0.0

    @field_validator("prs_samples", mode="before")
    @classmethod
    def _prs_samples(cls, v: Any) -> Any:
        """``all`` (or null) = every genotyped sample; ``file:<path>`` or a path = a subset."""
        if v == "all":
            return None
        return strip_file_prefix(v)

    @field_validator("lambda_", mode="before")
    @classmethod
    def _lambda_type(cls, v: Any) -> Any:
        if isinstance(v, bool):
            raise ValueError("lambda must be 'auto' or a number between 0 and 1")
        return v

    @field_validator("lambda_")
    @classmethod
    def _lambda(cls, v: float | str) -> float | str:
        if isinstance(v, float) and not (math.isfinite(v) and 0.0 <= v <= 1.0):
            raise ValueError("lambda must be 'auto' or between 0 and 1")
        return v


class CisConfig(StrictModel):
    """cis-xQTL power evaluation by replication of known pairs."""

    enabled: bool = True
    datasets: list[str] | None = None
    reference: str | None = None
    assoc_engine: Literal["pair", "rect"] = "pair"
    strategies: list[CisStrategyConfig] = Field(default_factory=_default_cis_strategies)
    """The run's inputs (one entry; several only in multi-pipeline configuration files)."""
    rint: bool = False
    """Apply a rank-based inverse normal transform (RINT) to the phenotypes after covariate
    adjustment (flag ``--rint``). Default: test the phenotypes as given. A run tests one
    transform; compare raw and RINT as two runs with ``anqchor cis compare``."""
    levels: LevelsConfig = Field(default_factory=LevelsConfig)
    max_unresolved_fraction: Annotated[float, Field(ge=0, le=1)] = 0.2
    """Fail when more than this fraction of reference pairs is absent or untested."""

    @model_validator(mode="after")
    def _unique(self) -> CisConfig:
        _check_unique_strategies(self.strategies)
        return self


class TransLevelsConfig(StrictModel):
    """Significance of a trans cell (a variant and a gene far apart). Each level is a
    threshold of the hotspot grid; ``nominal`` is also the default hotspot threshold of the
    report and the cut-off of ``pairs.sig.tsv.gz``. Reference-pair replication uses the cis
    levels (same sign, p < 0.05, p < 1e-5) instead."""

    nominal: Annotated[float, Field(gt=0, lt=1)] = 1e-5
    """A trans cell is significant at ``p < nominal``."""
    strict: Annotated[float, Field(gt=0, lt=1)] | None = None
    """An optional stricter threshold added to the hotspot grid (e.g. ``1e-8``)."""
    bonferroni: bool = True
    """Add ``p < 0.05 / trans cells tested`` to the hotspot grid."""

    @model_validator(mode="before")
    @classmethod
    def _no_fdr(cls, data: Any) -> Any:
        if isinstance(data, dict) and "fdr" in data:
            raise ValueError(
                "pipelines.trans.levels.fdr is no longer used: trans signals are counted at fixed "
                "p-value thresholds (nominal, strict, bonferroni); remove the key"
            )
        return data


class TransConfig(StrictModel):
    """trans-xQTL calibration and hotspot evaluation over a variant × trait rectangle."""

    enabled: bool = True
    datasets: list[str] | None = None
    reference: str | None = None
    variants: Literal["reference"] | Path = "reference"
    """``reference`` = the reference's variants, or a file with one ``CHROM:POS:REF:ALT`` per line."""
    traits: Literal["reference", "all"] | Path = "reference"
    """``reference`` = the reference's traits (the initial test rectangle), ``all``, or a file."""
    trait_chroms: Literal["genotyped", "all"] | list[str] = "genotyped"
    """Restrict tested traits by the chromosome of their locus (``traits.tsv``):
    ``genotyped`` (default) keeps traits on chromosomes present in the genotype set;
    ``all`` applies no restriction; a list (e.g. ``["21", "22"]``, either chromosome
    style) keeps only those. Traits without a locus are kept only with ``all``. With
    genome-wide genotypes ``genotyped`` keeps every annotated trait."""
    known_pairs: bool = True
    """Tag the reference pairs as known positives."""
    strategies: list[StrategyConfig] = Field(default_factory=_default_strategies)
    """The run's inputs: exactly one entry (one set of inputs per trans run)."""
    rint: bool = False
    """Apply a rank-based inverse normal transform (RINT) to the phenotypes after covariate
    adjustment (``rect-assoc --rint-after-adj``; flag ``--rint``). Default: test the
    phenotypes as given. A run tests one transform; raw and RINT are two runs."""
    cis_window_bp: Annotated[int, Field(ge=0)] = 1_000_000
    """A variant and a gene are cis when on the same chromosome with the variant within this
    distance of the gene start (TSS); every other pair with a known gene locus is trans."""
    min_maf: Annotated[float, Field(ge=0, lt=0.5)] = 0.01
    min_mac: Annotated[float, Field(ge=0)] = 1.0
    levels: TransLevelsConfig = Field(default_factory=TransLevelsConfig)
    hotspot_bin_bp: Annotated[int, Field(ge=1)] = 1_000_000
    """Width of the genomic bins that thin trans signals: a bin counts the distinct genes
    with a significant trans signal from any of its variants (genes in cis with the bin
    excluded)."""
    hotspot_min_traits: Annotated[int, Field(ge=1)] = 5
    """Genes a bin needs to be a hotspot (the report's default; the grid also shows 1, 2,
    3, 5, 10 and 20)."""
    max_chunk_vars: Annotated[int, Field(ge=1)] = 100
    max_sig_rows: Annotated[int, Field(ge=1)] = 1_000_000
    """Upper bound on rows written to ``pairs.sig.tsv.gz``."""

    @field_validator("trait_chroms", mode="before")
    @classmethod
    def _chroms_as_str(cls, v: Any) -> Any:
        """YAML ``[21, 22]`` arrives as integers; chromosome names are strings."""
        if isinstance(v, (list, tuple)):
            if not v:
                raise ValueError("trait_chroms: an empty list would select no traits")
            return [str(c) for c in v]
        return v

    @model_validator(mode="before")
    @classmethod
    def _reference_forms(cls, data: Any) -> Any:
        """Accept the documented spellings ``reference:<name>``, ``file:<path>`` and
        ``known_pairs: reference:<name> | null``; a named reference also sets ``reference``."""
        if not isinstance(data, dict):
            return data
        data = dict(data)
        named: list[str] = []
        for key in ("variants", "traits", "known_pairs"):
            if key not in data:
                continue
            v = data[key]
            if isinstance(v, str) and v.startswith(REFERENCE_PREFIX):
                named.append(v[len(REFERENCE_PREFIX) :])
                data[key] = True if key == "known_pairs" else "reference"
            elif key == "known_pairs":
                if v is None:
                    data[key] = False
                elif v == "reference":
                    data[key] = True
            else:
                data[key] = strip_file_prefix(v)
        names = set(named)
        ref = data.get("reference")
        if ref is not None:
            names.add(ref)
        if len(names) > 1:
            raise ValueError(
                "trans: variants/traits/known_pairs/reference name different references "
                f"({', '.join(sorted(names))}); one reference per trans run"
            )
        if named and ref is None:
            data["reference"] = named[0]
        return data

    @model_validator(mode="after")
    def _unique(self) -> TransConfig:
        _check_unique_strategies(self.strategies)
        if len(self.strategies) > 1:
            raise ValueError(
                f"pipelines.trans.strategies: a trans run evaluates one set of inputs, got "
                f"{len(self.strategies)} ({', '.join(s.name for s in self.strategies)}); make one run "
                "per covariate set (for example with --cov-include), each into its own --out directory"
            )
        return self


def _check_unique_strategies(strategies: Sequence[InputSetConfig]) -> None:
    names = [s.name for s in strategies]
    dup = sorted({n for n in names if names.count(n) > 1})
    if dup:
        raise ValueError(f"duplicate strategy names: {', '.join(dup)}")
    if not strategies:
        raise ValueError("at least one strategy is required")


class PipelinesConfig(StrictModel):
    match: MatchConfig = Field(default_factory=MatchConfig)
    cis: CisConfig = Field(default_factory=CisConfig)
    trans: TransConfig = Field(default_factory=TransConfig)


class ReportConfig(StrictModel):
    plotlyjs: Literal["embed", "cdn"] = "embed"
    """``embed`` writes plotly.js into each page (works offline); ``cdn`` keeps pages small."""


class Config(StrictModel):
    """Top-level anqchor configuration."""

    version: Literal[1] = 1
    project: str = "anqchor"
    out_dir: Path = Path("anqchor_out")
    qpgen: QpgenConfig = Field(default_factory=QpgenConfig)
    genotypes: GenotypeConfig
    references: dict[str, ReferenceConfig]
    datasets: dict[str, DatasetConfig]
    pipelines: PipelinesConfig = Field(default_factory=PipelinesConfig)
    report: ReportConfig = Field(default_factory=ReportConfig)

    @field_validator("references", "datasets")
    @classmethod
    def _names(cls, v: dict[str, object], info: object) -> dict[str, object]:
        kind = getattr(info, "field_name", "name")
        if not v:
            raise ValueError(f"at least one entry is required in '{kind}'")
        for key in v:
            _check_name(str(kind).rstrip("s"), key)
        return v

    @model_validator(mode="after")
    def _cross_refs(self) -> Config:
        for pname in ("match", "cis", "trans"):
            p = getattr(self.pipelines, pname)
            for ds in p.datasets or []:
                if ds not in self.datasets:
                    raise ValueError(f"pipelines.{pname}.datasets: unknown dataset {ds!r}")
            if p.reference is not None and p.reference not in self.references:
                raise ValueError(f"pipelines.{pname}.reference: unknown reference {p.reference!r}")
            # Only an enabled pipeline has to say which of several references it uses.
            if p.enabled and p.reference is None and len(self.references) > 1:
                raise ValueError(
                    f"pipelines.{pname}.reference must be set when several references are defined"
                )
        if self.pipelines.cis.enabled:
            ref = self.references[self.reference_for("cis")]
            if ref.pairs is None:
                raise ValueError("pipelines.cis needs a reference with 'pairs'")
        if self.pipelines.match.enabled:
            ref = self.references[self.reference_for("match")]
            if ref.pairs is None:
                raise ValueError(
                    "pipelines.match needs a reference with 'pairs' (effect sizes for the PRS)"
                )
        t = self.pipelines.trans
        if t.enabled:
            name = self.reference_for("trans")
            ref = self.references[name]
            if t.traits == "reference" and ref.pairs is None:
                raise ValueError(
                    f"pipelines.trans.traits: 'reference' needs pairs in reference {name!r}; "
                    "use 'all' or a trait file"
                )
        return self

    # ------------------------------------------------------------------ helpers
    def reference_for(self, pipeline: Literal["match", "cis", "trans"]) -> str:
        """Name of the reference used by ``pipeline``."""
        p = getattr(self.pipelines, pipeline)
        if p.reference is not None:
            return str(p.reference)
        return next(iter(self.references))

    def datasets_for(self, pipeline: Literal["match", "cis", "trans"]) -> list[str]:
        """Dataset names used by ``pipeline`` (all datasets when not restricted)."""
        p = getattr(self.pipelines, pipeline)
        return list(p.datasets) if p.datasets else list(self.datasets)

    def enabled_pipelines(self) -> list[str]:
        return [n for n in ("match", "cis", "trans") if getattr(self.pipelines, n).enabled]
