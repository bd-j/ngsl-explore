"""Where a figure goes.

Stars whose catalog [M/H] is below the grid's -0.5 floor cannot be represented
by any node, so their figures are not comparable with the rest: the fit has to
spend some other parameter to make up the difference, and on this sample that
shows up as v sin i running to the 300 km/s ceiling and as no node satisfying
both legs at once (PLAN.md). Each group therefore gets its own directory --
`figures/fits_mh_in_grid/` and `figures/fits_mh_below_grid/` -- so that neither
is read as the other, and neither sits loose in `figures/` alongside the
survey and instrument figures that are not per-star fits.

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
from pathlib import Path

from common.sample import sample_row, all_stars, opt_float

ROOT = Path(__file__).resolve().parent.parent
MH_FLOOR = -0.5                  # the model grid's lowest [M/H] node
BELOW_GRID_DIR = 'fits_mh_below_grid'   # [M/H] outside the grid's reach
IN_GRID_DIR = 'fits_mh_in_grid'         # representable by some node
LIBRARY_DIR = 'explore_libraries'       # survey figures, not per-star fits

# Stars the catalog cut sends below the grid but that the models in fact fit.
# Judged on the metal-line panels (`metal_lines_<star>.png`, in whichever
# directory this module sends them), which show directly whether the -0.5 node
# reproduces the observed line depths.
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


def catalog_mh(star, column='mh_ngsl'):
    """-> the star's catalog [M/H], or None if it is not in the sample."""
    return opt_float(sample_row(star, required=False).get(column))


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
    d = ROOT / 'figures' / (BELOW_GRID_DIR if below_grid(star) else IN_GRID_DIR)
    d.mkdir(parents=True, exist_ok=True)
    return d


def figure_path(prefix, star, suffix='.png'):
    """-> figures/fits_mh_{in,below}_grid/<prefix>_<star><suffix>."""
    return figure_dir(star) / f'{prefix}_{star}{suffix}'


def library_figure_path(name):
    """-> figures/explore_libraries/<name>, for the library survey figures.

    These characterise a LIBRARY -- NGSL or UVES-POP against the models, the
    Pickles atlas, the S/N survey -- rather than fitting one star, so they are
    not results in the sense the two fit directories are, and mixing them in
    with per-star figures is what made `figures/` unreadable in the first
    place. Same reason as figure_dir: the scripts regenerate, so the choice
    belongs here and not in four separate string literals.
    """
    d = ROOT / 'figures' / LIBRARY_DIR
    d.mkdir(parents=True, exist_ok=True)
    return d / name


def metal_poor_stars(column='mh_ngsl'):
    """-> sorted list of the sample stars below the grid floor."""
    return [s for s in all_stars()
            if sample_row(s)['tier'] != 'rejected' and below_grid(s, column)]
