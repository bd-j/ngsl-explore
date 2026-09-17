"""Reading `data/sample.csv`, and coercing its fields to numbers.

Every figure and fitting script starts by looking a star up in the sample table
and pulling a few catalog columns out of it as floats. That had grown four
copies of the row lookup and four of the float coercion, and the copies had
drifted -- in both cases into two genuinely different behaviours, which is why
this module exposes two functions rather than one of each.

WHY TWO FLOAT COERCIONS. A missing or unparseable CSV field has to become
something, and the right something depends on what happens next:

  as_float   -> np.nan   for values that go into numpy. `np.array([as_float(...)
               for ...])` stays a float array, arithmetic propagates, and
               `lo <= as_float(x) <= hi` is False for a missing value without a
               guard. `None` here would silently produce an object-dtype array
               and break np.isfinite.
  opt_float  -> None     for optional scalars that are tested before use, as
               `if v is None: continue`. `is None` is unambiguous and cannot be
               confused with a real 0.0 the way a truthiness test can. nan here
               would need np.isnan at every call site and is TRUTHY, so it slips
               straight past `if v:`.

Neither is a better default; substituting one for the other breaks the other's
call sites, so they are named for what they return.
"""
import csv
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
SAMPLE_CSV = ROOT / 'data' / 'sample.csv'
_ROWS = {}


def _load():
    if not _ROWS:
        with open(SAMPLE_CSV) as fh:
            for r in csv.DictReader(fh):
                _ROWS[r['star']] = r
    return _ROWS


def sample_row(star, required=True):
    """-> the star's row of data/sample.csv as a dict.

    required=True raises KeyError for a star that is not in the sample, which is
    what a single-star entry point wants: a mistyped name should fail loudly
    rather than draw an empty figure. required=False returns {} instead, for
    callers sweeping the whole sample that would rather fall back on `.get`
    defaults than abort the run over one missing row.
    """
    row = _load().get(star)
    if row is None:
        if required:
            raise KeyError(f'{star} is not in {SAMPLE_CSV}')
        return {}
    return row


def all_stars(tier=None):
    """-> sorted star names, optionally restricted to one tier."""
    return sorted(s for s, r in _load().items()
                  if tier is None or r['tier'] == tier)


def as_float(x):
    """-> float, or np.nan if missing or unparseable. For numpy consumers."""
    try:
        v = float(x)
    except (TypeError, ValueError):
        return np.nan
    return v if np.isfinite(v) else np.nan


def opt_float(x):
    """-> float, or None if missing or unparseable. For `if v is None` guards."""
    try:
        v = float(x)
    except (TypeError, ValueError):
        return None
    return v if np.isfinite(v) else None
