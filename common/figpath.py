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
(fitting/observations.py, explore/plot_ebv_teff.py:nodes_for). No star sits
exactly on -0.5, so the boundary itself is unambiguous.

WHY THERE IS AN EXCEPTION LIST, AND WHY IT IS NOT A FORMULA. The catalog value
is a proxy for the question that actually matters -- can the grid reproduce
this star's metal lines? -- and for two stars the proxy is wrong. An empirical
rule was tried and rejected. Counting features whose observed depth is
SHALLOWER than any grid [M/H] can reach is the right physical signal, but it is
confounded by v sin i: HD128801 and HD106304 sit at the 300 km/s ceiling, which
smears the MODEL lines flat, so their features register as "deeper than the
grid" and they would score as WELL fit -- exactly backwards for the two stars
most in need of segregating. So the exception is a short list with a reason per
star, which is reviewable in a way a mis-tuned threshold is not.
"""
import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
MH_FLOOR = -0.5                  # the model grid's lowest [M/H] node
METAL_POOR_DIR = 'below_grid_mh'

# Stars the catalog cut sends below the grid but that the models in fact fit.
# Judged on the metal-line panels (figures/metal_lines_<star>.png), which show
# directly whether the -0.5 node reproduces the observed line depths.
IN_GRID_ANYWAY = {
    'HD164967': 'catalogs disagree (NGSL -0.60, XSL -0.29) and the metal lines '
                'are well fit at the floor node: 1 of 8 features outside the '
                'grid range, NONE shallower than the grid can reach',
    'HD143459': 'catalogs disagree (NGSL -0.60, XSL -0.19). Marginal but kept: '
                '2 of 8 features are shallower than the grid can reach, so it '
                'probably does sit near -0.5 -- which is a node the grid HAS, '
                'not a metallicity it cannot represent. XSL -0.19 is the '
                'outlier here, not NGSL -0.60',
}
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
    """True if this star is out of the grid's [M/H] reach.

    The catalog cut, minus the reviewed exceptions in IN_GRID_ANYWAY.
    """
    if star in IN_GRID_ANYWAY:
        return False
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
