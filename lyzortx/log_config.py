"""Shared logging configuration for the lyzortx package."""

from __future__ import annotations

import logging
import time

LOG_FORMAT = "%(asctime)s %(levelname)s %(name)s: %(message)s"


def setup_logging(level: int = logging.INFO) -> None:
    """Configure the root logger with a timezone-aware timestamped format.

    Call once at the entry point of each pipeline runner (e.g., in ``main()``
    of ``run_track_*.py``).  Library modules should only create their own
    logger via ``logging.getLogger(__name__)`` — never call this function.

    Timestamps render in the computer's local timezone per AGENTS.md. The
    previous configuration used ``time.gmtime`` which produced UTC hour values
    paired with a ``%z`` offset that fell back to the local non-DST offset
    (e.g. CET +0100 during CEST), producing timestamps that matched neither
    wall clock nor UTC.
    """
    formatter = logging.Formatter(LOG_FORMAT, datefmt="%Y-%m-%d %H:%M:%S%z")
    formatter.converter = time.localtime
    handler = logging.StreamHandler()
    handler.setFormatter(formatter)
    logging.basicConfig(handlers=[handler], level=level)
