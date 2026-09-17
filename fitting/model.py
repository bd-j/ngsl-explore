"""Forward model: grid interpolation -> reddening -> broadening -> data grid.

Deliberately NO free continuum polynomial. A multiplicative polynomial would
absorb exactly the information that constrains Teff and E(B-V) -- the continuum
shape -- leaving only line profiles to carry the temperature. The normalization
is therefore a single scalar, and the fit relies on the spectrophotometry being
good (NGSL is space-based and calibrated to ~3%; UVES-POP quotes 1.5-4%).

The consequence is that Teff and E(B-V) are covariant: both tilt the continuum.
They are separated by (a) the Balmer break amplitude, which responds to Teff but
barely to a smooth reddening law, (b) a wide wavelength baseline, and (c) a
one-sided prior on E(B-V) from a dust map. See fit.py.
"""
import sys
from pathlib import Path

import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

ROOT = Path(__file__).resolve().parent.parent


class Grid:
    """Packed model grid with trilinear interpolation in log flux."""

    def __init__(self, path=None):
        d = np.load(path or ROOT / 'models' / 'grid.npz')
        self.teff, self.logg, self.mh = d['teff'], d['logg'], d['mh']
        self.wave, self.logflux = d['wave'], d['logflux']
        self.filled = d['filled']
        self.resolution = float(d['resolution'])

    def bounds(self):
        return dict(teff=(self.teff.min(), self.teff.max()),
                    logg=(self.logg.min(), self.logg.max()),
                    mh=(self.mh.min(), self.mh.max()))

    def _axis(self, arr, x):
        """Bracketing indices and weight, clipped to the grid edge."""
        if len(arr) == 1:
            return 0, 0, 0.0
        i = int(np.clip(np.searchsorted(arr, x) - 1, 0, len(arr) - 2))
        w = (x - arr[i]) / (arr[i + 1] - arr[i])
        return i, i + 1, float(np.clip(w, 0.0, 1.0))

    def interp(self, teff, logg, mh):
        """-> f_lambda on the grid's own wavelength array."""
        i0, i1, a = self._axis(self.teff, teff)
        j0, j1, b = self._axis(self.logg, logg)
        k0, k1, c = self._axis(self.mh, mh)
        out = np.zeros_like(self.logflux[0, 0, 0], dtype=np.float64)
        for i, wi in ((i0, 1 - a), (i1, a)):
            if wi == 0:
                continue
            for j, wj in ((j0, 1 - b), (j1, b)):
                if wj == 0:
                    continue
                for k, wk in ((k0, 1 - c), (k1, c)):
                    if wk == 0:
                        continue
                    if not self.filled[i, j, k]:
                        raise ValueError(
                            f'grid cell ({self.teff[i]:.0f}, {self.logg[j]:.2f}, '
                            f'{self.mh[k]:+.2f}) is missing; the cube is incomplete')
                    out += wi * wj * wk * self.logflux[i, j, k]
        return 10.0 ** out


# REMOVED: `instrument()` and `forward()` used to live here, a second copy of
# the forward model that nothing called. They documented the NGSL profile as
# "constant R = 600", which stopped being true and was never updated, so the
# file sat there contradicting common/lsf.py. fitting/predict.py is the forward
# model; common.lsf.broaden_ngsl is the instrument profile. One of each.
