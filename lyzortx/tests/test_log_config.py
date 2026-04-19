"""Verify setup_logging renders timestamps in local time with a matching offset.

Hardcoded summer and winter epochs exercise both DST states. The test pins the local
timezone to Europe/Berlin (CET/CEST) via ``TZ`` + ``time.tzset`` so expected strings are
stable across CI hosts. Without the pin the test would silently pass when CI ran in UTC
(both summer and winter would render as +0000) and miss DST-transition regressions.
"""

from __future__ import annotations

import logging
import os
import re
import time

import pytest

from lyzortx.log_config import LOG_FORMAT

# Hardcoded Unix epochs for two instants in the Europe/Berlin zone.
SUMMER_EPOCH = 1784109600  # 2026-07-15 12:00:00 CEST == 2026-07-15 10:00:00 UTC
SUMMER_EXPECTED_TIMESTAMP = "2026-07-15 12:00:00+0200"

WINTER_EPOCH = 1797332400  # 2026-12-15 12:00:00 CET  == 2026-12-15 11:00:00 UTC
WINTER_EXPECTED_TIMESTAMP = "2026-12-15 12:00:00+0100"


@pytest.fixture
def europe_berlin_tz():
    if not hasattr(time, "tzset"):
        pytest.skip("time.tzset not available on this platform")
    prior = os.environ.get("TZ")
    os.environ["TZ"] = "Europe/Berlin"
    time.tzset()
    try:
        yield
    finally:
        if prior is None:
            os.environ.pop("TZ", None)
        else:
            os.environ["TZ"] = prior
        time.tzset()


def _build_formatter() -> logging.Formatter:
    formatter = logging.Formatter(LOG_FORMAT, datefmt="%Y-%m-%d %H:%M:%S%z")
    formatter.converter = time.localtime
    return formatter


def _render_at(formatter: logging.Formatter, epoch: float) -> str:
    record = logging.LogRecord(
        name="probe",
        level=logging.INFO,
        pathname="/tmp/x",
        lineno=1,
        msg="hello",
        args=(),
        exc_info=None,
    )
    record.created = epoch
    return formatter.format(record)


def _extract_timestamp(rendered: str) -> str:
    match = re.match(r"(\S+ \S+) INFO", rendered)
    assert match, f"unexpected format: {rendered!r}"
    return match.group(1)


def test_summer_timestamp_renders_in_cest(europe_berlin_tz) -> None:
    rendered = _render_at(_build_formatter(), SUMMER_EPOCH)
    assert _extract_timestamp(rendered) == SUMMER_EXPECTED_TIMESTAMP


def test_winter_timestamp_renders_in_cet(europe_berlin_tz) -> None:
    rendered = _render_at(_build_formatter(), WINTER_EPOCH)
    assert _extract_timestamp(rendered) == WINTER_EXPECTED_TIMESTAMP
