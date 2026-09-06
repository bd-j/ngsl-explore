"""Shared Balmer-discontinuity estimator, so NGSL and Pickles are comparable.

D = 2.5 log10( Fc_red(3646) / Fc_blue(3646) )

Both continua are extrapolated to the Balmer limit itself before taking the
ratio, so a sloping SED cancels. Each side uses an upper-envelope percentile
fit, which tracks the continuum peaks between absorption lines rather than
being dragged down by line cores.

LIMITATION: meaningful only for hot stars (Teff above ~6000 K), where the
continuum is close to linear on both sides of the limit. For cool stars the
continuum is strongly curved and heavily line-blanketed, and D is dominated
by SED slope rather than by the Balmer discontinuity - it should not be read
as a break amplitude there.
"""
import numpy as np

BALMER = 3646.0
BLUE = (3350.0, 3630.0)
RED = (3700.0, 4150.0)


def _envelope(w, f, lo, hi, pct, deg=1):
    m = (w > lo) & (w < hi) & np.isfinite(f)
    if m.sum() < 8:
        return None
    x, y = w[m], f[m]
    keep = y > np.percentile(y, pct)
    if keep.sum() < 3:
        keep = np.ones_like(y, dtype=bool)
    return float(np.polyval(np.polyfit(x[keep], y[keep], deg), BALMER))


def balmer_discontinuity(w, f):
    """-> (D in mag, blue continuum at 3646, red continuum at 3646)."""
    fc_blue = _envelope(w, f, *BLUE, pct=80)
    fc_red = _envelope(w, f, *RED, pct=88)
    if not fc_blue or not fc_red or fc_blue <= 0 or fc_red <= 0:
        return None, fc_blue, fc_red
    return 2.5 * np.log10(fc_red / fc_blue), fc_blue, fc_red
