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
from scipy.ndimage import gaussian_filter1d

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common.extinction_ccm import redden
from common.lsf import (broaden_rot, broaden_R, broaden_ngsl,
                        broaden, NGSL_R_MEASURED)

ROOT = Path(__file__).resolve().parent.parent
C_KMS = 2.99792458e5


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


def instrument(w, f, kind, value):
    """Instrumental broadening. Delegates to common.lsf so the fitter and the
    comparison figures cannot drift apart.

    kind='R'      constant resolving power (value = R)
    kind='fwhm'   constant FWHM in Angstroms (value = FWHM)
    kind='ngsl'   the MEASURED NGSL profile: constant R = 600 (see docs/DATA.md).
                  Not the STIS-table profile, which is 1.7-1.8x too narrow for
                  the delivered spectra.
    kind='ngsl_tab'  the STIS-table profile, for comparison only.

    Both 'R' and 'ngsl' were previously computed here by hand and were wrong:
    the constant-R sigma ignored the grid spacing (R=600 came out as R=83), and
    'ngsl' still applied the superseded 3.85 A tabulated width.
    """
    if kind == 'ngsl':
        return broaden_ngsl(w, f)
    if kind == 'ngsl_tab':
        return broaden_ngsl(w, f, tabulated=True)
    if not value or value <= 0:
        return f
    if kind == 'R':
        return broaden_R(w, f, value)
    if kind == 'fwhm':
        return broaden(w, f, value)
    raise ValueError(f'unknown instrument kind: {kind}')


def forward(grid, wave_obs, teff, logg, mh, ebv=0.0, vsini=0.0, rv=0.0,
            inst_kind='R', inst_value=0.0, r_v=3.1):
    """Model flux on the observed wavelength grid, up to a scalar normalization.

    Order: interpolate -> rotate (star frame) -> Doppler shift -> redden
    (observer frame) -> instrument -> resample. Rotation before the instrument
    profile because they are physically sequential, and the two do not commute
    once the rotation profile is not Gaussian.
    """
    w = grid.wave
    f = grid.interp(teff, logg, mh)
    f = broaden_rot(w, f, vsini)
    ws = w * (1.0 + rv / C_KMS)
    if ebv:
        f = redden(ws, f, ebv, r_v)
    f = instrument(ws, f, inst_kind, inst_value)
    return np.interp(wave_obs, ws, f)
