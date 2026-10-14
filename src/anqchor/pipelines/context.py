"""Shared state for pipeline step builders."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from anqchor.config.model import Config
from anqchor.pipelines.layout import Layout
from anqchor.qpgen.commands import Genotypes


@dataclass(frozen=True)
class BuildContext:
    """What a pipeline needs to emit its steps.

    Attributes:
        cfg: the resolved configuration (absolute paths).
        qpgen: absolute path of ``qpgentools``.
        layout: run-directory file names.
        qpgen_sha256: content hash of ``qpgentools``. :func:`anqchor.pipelines.build_graph`
            adds it to the parameters of every step that runs qpgen, so rebuilding the
            binary (even in place, at the same path) makes those steps and everything
            downstream stale. ``None`` (tests) leaves the steps untagged.
    """

    cfg: Config
    qpgen: Path
    layout: type[Layout] = Layout
    qpgen_sha256: str | None = None

    @property
    def threads(self) -> int:
        return self.cfg.qpgen.threads

    def genotypes(self) -> Genotypes:
        """Genotype arguments for qpgen: always the prepared, absolute pgen list."""
        return Genotypes(pgen_list=self.layout.pgen_list)
