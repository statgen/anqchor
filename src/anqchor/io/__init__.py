"""Input readers, converters and validators (the ``prepare`` layer's toolbox).

* :mod:`anqchor.io.tables`: separators, atomic writers, :class:`InputError`.
* :mod:`anqchor.io.genotypes`: pgen lists, psam IDs, indexed pvar, chromosome style.
* :mod:`anqchor.io.pheno`: phenotype/covariate matrices in the four layouts, trait loci.
* :mod:`anqchor.io.cov`: covariate encoding and selection.
* :mod:`anqchor.io.pairs`: reference pairs, CPRA parsing, variant presence checks.
* :mod:`anqchor.io.maps`: trait-ID rules, trait maps and sample maps.

Functions raise :class:`InputError` for bad data; step implementations turn it into
:class:`anqchor.steps.StepError`.
"""

from anqchor.io.tables import InputError

__all__ = ["InputError"]
