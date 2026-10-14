"""The run-directory layout: the single source of truth for file names.

Every path is relative to the run directory (make runs there). Pipelines and step
implementations must obtain paths from :class:`Layout` instead of formatting them,
so that producers and consumers cannot drift apart. The file formats are documented
in ``docs/outputs.md``.
"""

from __future__ import annotations

from pathlib import Path


class Layout:
    """Relative paths of every file anqchor writes in a run directory."""

    # ------------------------------------------------------------------ general
    resolved_config = Path("anqchor.resolved.yaml")
    run_json = Path("run.json")
    makefile = Path("Makefile")

    # --------------------------------------------------------------- genotypes
    geno_dir = Path("prepare/genotypes")
    pgen_list = geno_dir / "pgen_list.tsv"
    """qpgen pgen list, 4 columns ``CHROM BEG END PREFIX`` with absolute prefixes, no header.

    The 6-column form is avoided: qpgen's parser reads it as ``PGEN PIVAR PSAM`` while its
    documentation says ``PGEN PSAM PIVAR``. Prefixes without standard file names are
    symlinked under ``prepare/genotypes/links/``."""
    geno_samples = geno_dir / "samples.txt"
    """Genotype sample IIDs, one per line, in psam order."""
    geno_json = geno_dir / "genotypes.json"
    """Counts and properties: n_samples, chrom_style ('chr' or 'plain'), chromosomes, files."""

    # -------------------------------------------------------------- references
    @staticmethod
    def ref_dir(ref: str) -> Path:
        return Path("prepare/references") / ref

    @classmethod
    def ref_pairs(cls, ref: str) -> Path:
        """Normalised pairs: header ``trait variant beta se pval`` (pval may be NA);
        variant IDs in the genotype chromosome style; duplicates removed."""
        return cls.ref_dir(ref) / "pairs.tsv.gz"

    @classmethod
    def ref_variants(cls, ref: str) -> Path:
        """Header ``variant chrom pos ref alt is_indel present`` (one row per unique variant)."""
        return cls.ref_dir(ref) / "variants.tsv"

    @classmethod
    def ref_variants_present(cls, ref: str) -> Path:
        """Variants present in the genotype pvar, one ``CHROM:POS:REF:ALT`` per line, no header,
        sorted by chromosome and position."""
        return cls.ref_dir(ref) / "variants.present.txt"

    @classmethod
    def ref_json(cls, ref: str) -> Path:
        return cls.ref_dir(ref) / "reference.json"

    # ---------------------------------------------------------------- datasets
    @staticmethod
    def ds_dir(ds: str) -> Path:
        return Path("prepare/datasets") / ds

    @classmethod
    def ds_pheno(cls, ds: str) -> Path:
        """Phenotype matrix in regenie layout (``FID IID trait…``), the dataset's sample IDs."""
        return cls.ds_dir(ds) / "pheno.regenie.tsv.gz"

    @classmethod
    def ds_cov(cls, ds: str) -> Path:
        """Covariates in regenie layout, written only when the dataset has ``cov``. Categorical
        covariates keep their text levels (e.g. ``XX``/``XY``); qpgen encodes them
        (``--cat-cov-names``). Constant numeric columns are dropped; missing values are
        handled by ``cov_missing``."""
        return cls.ds_dir(ds) / "cov.regenie.tsv.gz"

    @classmethod
    def ds_cov_json(cls, ds: str) -> Path:
        """Covariate description: ``raw_columns`` (columns of ``cov.regenie.tsv.gz``),
        ``categorical_names`` (passed to qpgen as ``--cat-cov-names``), ``columns`` and
        ``groups`` (the equivalent numeric encoding, used to count covariates for the
        degrees of freedom and to resolve ``--cov-include`` names), ``dropped_*``."""
        return cls.ds_dir(ds) / "cov.json"

    @classmethod
    def ds_traits(cls, ds: str) -> Path:
        """Header ``trait chrom start end``; chrom/start/end are NA when no annotation exists.
        Chromosome names follow the genotype chromosome style."""
        return cls.ds_dir(ds) / "traits.tsv"

    @classmethod
    def ds_samples(cls, ds: str) -> Path:
        """Phenotype sample IDs (IID), one per line."""
        return cls.ds_dir(ds) / "samples.txt"

    @classmethod
    def ds_sample_map(cls, ds: str) -> Path:
        """Identity matching: ``geno_id pheno_id`` without header, restricted to IDs present
        on both sides. Identity pairs when the dataset has no ``sample_map``."""
        return cls.ds_dir(ds) / "sample_map.tsv"

    @classmethod
    def ds_json(cls, ds: str) -> Path:
        return cls.ds_dir(ds) / "dataset.json"

    # ------------------------------------------------- dataset × reference pairing
    @classmethod
    def dsref_dir(cls, ds: str, ref: str) -> Path:
        return cls.ds_dir(ds) / f"ref.{ref}"

    @classmethod
    def dsref_trait_map(cls, ds: str, ref: str) -> Path:
        """``ref_trait pheno_trait`` without header, for traits present on both sides."""
        return cls.dsref_dir(ds, ref) / "trait_map.tsv"

    @classmethod
    def dsref_pairs(cls, ds: str, ref: str) -> Path:
        """Reference pairs rewritten to dataset trait IDs: header
        ``trait variant beta se pval ref_trait present`` (all mapped pairs; ``present`` is 0/1)."""
        return cls.dsref_dir(ds, ref) / "pairs.tsv.gz"

    @classmethod
    def dsref_variants(cls, ds: str, ref: str) -> Path:
        """Present variants of the mapped pairs, one CPRA per line, sorted, no header."""
        return cls.dsref_dir(ds, ref) / "variants.txt"

    @classmethod
    def dsref_traits(cls, ds: str, ref: str) -> Path:
        """Dataset trait IDs of the mapped reference traits, one per line, no header."""
        return cls.dsref_dir(ds, ref) / "traits.txt"

    @classmethod
    def dsref_json(cls, ds: str, ref: str) -> Path:
        return cls.dsref_dir(ds, ref) / "mapping.json"

    # --------------------------------------------------------------------- PRS
    @staticmethod
    def prs_prefix(ref: str) -> Path:
        return Path("prs") / ref / "prs"

    @classmethod
    def prs_matrix(cls, ref: str) -> Path:
        return cls.prs_prefix(ref).with_name("prs.prs.tsv.gz")

    @classmethod
    def prs_se(cls, ref: str) -> Path:
        return cls.prs_prefix(ref).with_name("prs.se.tsv.gz")

    # ------------------------------------------------------------------- match
    @staticmethod
    def match_dir(ds: str) -> Path:
        return Path("match") / ds

    @classmethod
    def match_prefix(cls, ds: str) -> Path:
        return cls.match_dir(ds) / "qpgen"

    @classmethod
    def match_qpgen_outputs(cls, ds: str) -> tuple[Path, Path, Path]:
        """(weights, match.assigned, match.all) written by ``match-prs-pheno``."""
        p = cls.match_dir(ds)
        return (
            p / "qpgen.weights.tsv.gz",
            p / "qpgen.match.assigned.tsv.gz",
            p / "qpgen.match.all.tsv.gz",
        )

    @classmethod
    def match_trait_map(cls, ds: str) -> Path:
        """Trait map actually passed to qpgen (restricted to traits in the PRS matrix)."""
        return cls.match_dir(ds) / "trait_map.used.tsv"

    @classmethod
    def match_params(cls, ds: str) -> Path:
        """Parameters resolved at run time (e.g. auto-λ decision, λ chosen)."""
        return cls.match_dir(ds) / "qpgen.params.json"

    @classmethod
    def match_summary(cls, ds: str, name: str) -> Path:
        """Summaries: ``samples.tsv``, ``swaps.tsv``, ``candidates.tsv``, ``related.tsv``,
        ``weights.tsv``, ``summary.json``."""
        return cls.match_dir(ds) / name

    # --------------------------------------------------------------------- cis
    @staticmethod
    def cis_dir(ds: str) -> Path:
        return Path("cis") / ds

    @classmethod
    def strategy_dir(cls, pipeline: str, ds: str, strategy: str) -> Path:
        return Path(pipeline) / ds / "strategies" / strategy

    @classmethod
    def strategy_cov(cls, pipeline: str, ds: str, strategy: str) -> Path:
        """Strategy covariates (regenie layout); absent when the strategy uses none."""
        return cls.strategy_dir(pipeline, ds, strategy) / "cov.regenie.tsv.gz"

    @classmethod
    def strategy_pheno(cls, pipeline: str, ds: str, strategy: str) -> Path:
        """Strategy phenotype matrix (regenie); written only when the strategy overrides ``pheno``."""
        return cls.strategy_dir(pipeline, ds, strategy) / "pheno.regenie.tsv.gz"

    @classmethod
    def strategy_json(cls, pipeline: str, ds: str, strategy: str) -> Path:
        return cls.strategy_dir(pipeline, ds, strategy) / "strategy.json"

    @classmethod
    def cis_pair_assoc(cls, ds: str, strategy: str) -> Path:
        return cls.strategy_dir("cis", ds, strategy) / "pair_assoc.tsv.gz"

    @classmethod
    def rect_assoc(cls, pipeline: str, ds: str, strategy: str, transform: str) -> Path:
        return cls.strategy_dir(pipeline, ds, strategy) / f"rect_assoc.{transform}.tsv.gz"

    @classmethod
    def cis_summary(cls, ds: str, name: str) -> Path:
        """Summaries: ``pairs.long.tsv.gz``, ``traits.tsv``, ``strategies.tsv``,
        ``power_curve.tsv``, ``summary.json``."""
        return cls.cis_dir(ds) / name

    # ------------------------------------------------------------------- trans
    @staticmethod
    def trans_dir(ds: str) -> Path:
        return Path("trans") / ds

    @classmethod
    def trans_variants(cls, ds: str) -> Path:
        return cls.trans_dir(ds) / "variants.txt"

    @classmethod
    def trans_traits(cls, ds: str) -> Path:
        """Traits tested (dataset IDs). qpgen gets no ``--pheno-list`` only when
        ``traits: all`` and ``trait_chroms: all``."""
        return cls.trans_dir(ds) / "traits.txt"

    @classmethod
    def trans_lists_json(cls, ds: str) -> Path:
        """Variant/trait selection counts of ``trans.<ds>.lists``."""
        return cls.trans_dir(ds) / "lists.json"

    @classmethod
    def trans_summary(cls, ds: str, name: str) -> Path:
        """Summaries: ``summary.json``, ``reference_pairs.tsv``, ``hotspots.tsv``,
        ``hotspot_genes.tsv``, ``qq.tsv``, ``variants.tsv``, ``pairs.sig.tsv.gz``."""
        return cls.trans_dir(ds) / name

    # ------------------------------------------------------- comparison directory
    compare_manifest = Path("compare.json")
    """Manifest of a ``cis compare`` or ``trans compare`` directory: ``kind`` (``cis`` or
    ``trans``), ``runs`` (``label run_dir dataset settings``), ``labels``, ``argv``, ``cwd``.
    ``cis report`` / ``trans report`` use it to recognise a comparison directory."""

    @staticmethod
    def compare_table(name: str) -> Path:
        """A table at the top of a comparison directory. ``cis compare``: ``runs.tsv``,
        ``strategies.tsv``, ``pairs.long.tsv.gz``, ``traits.tsv``, ``power_curve.tsv``,
        ``summary.json``. ``trans compare``: ``runs.tsv``, ``hotspot_grid.tsv``,
        ``hotspots.tsv``, ``hotspot_genes.tsv``, ``qq.tsv``, ``reference_pairs.tsv``,
        ``summary.json``."""
        return Path(name)

    @staticmethod
    def compare_settings(label: str) -> Path:
        """``settings/<label>.anqchor.yaml``: a copy of a compared run's resolved settings."""
        return Path("settings") / f"{label}.anqchor.yaml"

    # ------------------------------------------------------------------ report
    @staticmethod
    def report(name: str) -> Path:
        """``report/<name>.html`` for ``match``, ``cis``, ``trans`` and ``index``."""
        return Path("report") / f"{name}.html"
