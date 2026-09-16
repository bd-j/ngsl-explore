"""Where a per-star figure goes.

Stars whose catalog [M/H] is below the grid's -0.5 floor cannot be represented
by any node, so their figures are not comparable with the rest: the fit has to
spend some other parameter to make up the difference, and on this sample that
shows up as v sin i running to the 300 km/s ceiling and as no node satisfying
both legs at once (PLAN.md). Keeping them in their own directory stops them
being read as ordinary results.

This lives in code rather than being a one-off `mv` for the usual reason: the
three figure scripts all take --all, and the next regeneration would put them
straight back into figures/.

WHICH CATALOG. The cut uses `mh_ngsl`, because that is the column the rest of
the pipeline already uses to place a star on the grid
(fitting/observations.py, explore/plot_ebv_teff.py:nodes_for). The two
catalogues disagree for two of the six -- XSL puts HD143459 at -0.19 and
HD164967 at -0.29 against NGSL's -0.60 -- so those two are here on the NGSL
value alone. The other four are below -0.5 in both. No star sits exactly on
-0.5, so the boundary is unambiguous.
"""
import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
MH_FLOOR = -0.5                  # the model grid's lowest [M/H] node
METAL_POOR_DIR = 'below_grid_mh'
_CACHE = {}


def catalog_mh(star, column='mh_ngsl'):
    """-> the star's catalog [M/H], or None if it is not in the sample."""
    if not _CACHE:
        for r in csv.DictReader(open(ROOT / 'data' / 'sample.csv')):
            _CACHE[r['star']] = r
    row = _CACHE.get(star)
    if row is None:
        return None
    try:
        return float(row[column])
    except (KeyError, TypeError, ValueError):
        return None


def below_grid(star, column='mh_ngsl'):
    """True if the catalog puts this star below the grid's [M/H] floor."""
    z = catalog_mh(star, column)
    return z is not None and z < MH_FLOOR


def figure_dir(star):
    """-> the directory this star's figures belong in (created if needed)."""
    d = ROOT / 'figures'
    if below_grid(star):
        d = d / METAL_POOR_DIR
    d.mkdir(parents=True, exist_ok=True)
    return d


def figure_path(prefix, star, suffix='.png'):
    """-> figures/[below_grid_mh/]<prefix>_<star><suffix>."""
    return figure_dir(star) / f'{prefix}_{star}{suffix}'


def metal_poor_stars(column='mh_ngsl'):
    """-> sorted list of the sample stars below the grid floor."""
    if not _CACHE:
        catalog_mh('')
    return sorted(s for s, r in _CACHE.items()
                  if r['tier'] != 'rejected' and below_grid(s, column))
