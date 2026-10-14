"""Significance predicates shared by the cis and trans pipelines.

One convention everywhere, matching the documented level rules:

* fixed p-value levels (``nominal``, ``strict``) and Bonferroni: ``p < alpha``, i.e.
  ``-log10 p > -log10 alpha``;
* ``fdr`` (Benjamini-Hochberg): the cells at or above the BH cut-off on ``-log10 p``
  (``p_(i) <= i q / m``), which in q-value terms is ``q <= alpha``.

qpgen prints ``LOG10P`` with 6 significant digits, so values exactly at a threshold
(``5`` for ``p = 1e-5``) do occur in large rectangles; both pipelines must classify them
the same way.
"""

from __future__ import annotations

import math

import numpy as np


def p_level_threshold(alpha: float) -> float:
    """``-log10 alpha``."""
    return -math.log10(float(alpha))


def passes(log10p: np.ndarray | float, threshold: float, level: str = "p") -> np.ndarray:
    """Whether ``-log10 p`` passes a level's ``-log10`` threshold.

    ``level == "fdr"`` uses the inclusive Benjamini-Hochberg rule (``>=`` the cut-off);
    every other level uses ``p < alpha`` (``>``). NaN never passes."""
    x = np.asarray(log10p)  # keep float32 chunks as they are (no copy)
    with np.errstate(invalid="ignore"):
        return x >= threshold if level == "fdr" else x > threshold
