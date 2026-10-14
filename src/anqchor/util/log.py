"""Logging setup shared by the CLI and by step implementations."""

from __future__ import annotations

import logging
import sys

LOGGER_NAME = "anqchor"


def get_logger(name: str | None = None) -> logging.Logger:
    """Return the package logger or one of its children (``anqchor.<name>``)."""
    return logging.getLogger(LOGGER_NAME if not name else f"{LOGGER_NAME}.{name}")


def setup_logging(verbose: bool = False) -> logging.Logger:
    """Configure the ``anqchor`` logger once; repeated calls only adjust the level."""
    logger = logging.getLogger(LOGGER_NAME)
    if not logger.handlers:
        handler = logging.StreamHandler(sys.stderr)
        handler.setFormatter(
            logging.Formatter(
                "[%(asctime)s] %(levelname)s %(name)s: %(message)s", "%Y-%m-%d %H:%M:%S"
            )
        )
        logger.addHandler(handler)
        logger.propagate = False
    logger.setLevel(logging.DEBUG if verbose else logging.INFO)
    return logger
