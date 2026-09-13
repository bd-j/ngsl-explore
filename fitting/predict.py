"""predict(theta, observations) -> one prediction per observation.

One physical model, N projections. Everything that is a property of the star
lives in `theta`; everything that is a property of an instrument lives on the
Observation. That split is what lets the node scan and an MCMC/nested sampler
share the identical forward model rather than each keeping its own copy.

Order of operations, and why:

    interpolate (or look up a node) -> rotate -> Doppler shift -> redden
    -> instrument -> project onto the dataset

Rotation precedes the instrument profile because they are physically sequential
and do not commute once the rotation profile is not Gaussian. Reddening sits
between them: it happens in the ISM, after the star and before the telescope.

NODE-EXACT FAST PATH. When Teff, log g and [M/H] land on grid nodes the cube is
read directly with no interpolation, which is what the 1705-node scan wants. The
same function interpolates when they do not, so switching to a sampler changes
the caller and not the physics. The node test is a tolerance, never floating
point equality.
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common.extinction_ccm import redden
from common.lsf import broaden_rot, broaden_R

C_KMS = 2.99792458e5

# ABSOLUTE tolerances for "is this parameter on a grid node". Deliberately far
# below any axis step (100 K, 0.2 dex, 0.2 dex): the fast path must trigger only
# for a value the caller meant to place on a node. A tolerance scaled to the step
# would quietly accept 10241 K as the 10200 K node and return the wrong spectrum
# while reporting an exact lookup.
NODE_ATOL = dict(teff=1e-3, logg=1e-6, mh=1e-6)


class Prediction:
    """A model projected onto one observation, before calibration."""

    def __init__(self, obs, value, wavelength=None):
        self.name, self.star = obs.name, obs.star
        self.value = value                  # same units/shape as obs.flux
        self.wavelength = wavelength

    def __repr__(self):
        return f'<Prediction {self.star}/{self.name} n={np.size(self.value)}>'


def _on_node(axis, x, atol):
    """Index of the node at x, or None if x is not on a node within atol."""
    i = int(np.argmin(np.abs(np.asarray(axis) - x)))
    return i if abs(axis[i] - x) <= atol else None


def spectrum_at(grid, teff, logg, mh):
    """f_lambda on the grid's own wavelength array.

    Reads the node directly when all three parameters sit on nodes; otherwise
    trilinear in log flux (the Planck function is exponential in 1/T, so log
    flux is much closer to linear in the grid parameters than flux is).
    """
    i = _on_node(grid.teff, teff, NODE_ATOL['teff'])
    j = _on_node(grid.logg, logg, NODE_ATOL['logg'])
    k = _on_node(grid.mh, mh, NODE_ATOL['mh'])
    if i is not None and j is not None and k is not None:
        if not grid.filled[i, j, k]:
            raise ValueError(f'node ({grid.teff[i]:.0f}, {grid.logg[j]:.2f}, '
                             f'{grid.mh[k]:+.2f}) is missing from the cube')
        return 10.0 ** grid.logflux[i, j, k]
    return grid.interp(teff, logg, mh)


def instrument(w, f, resolution):
    """Smooth to a dataset's resolution."""
    kind = resolution[0]
    if kind is None:
        return f
    if kind == 'R':
        return broaden_R(w, f, float(resolution[1]))
    if kind == 'R_segments':
        # Each arm has its own constant-R kernel. Segments are convolved
        # separately and stitched; the model grid is log-sampled so a
        # constant-R kernel is a constant number of pixels within a segment.
        out = np.array(f, float)
        for lo, hi, R in resolution[1]:
            seg = (w >= lo) & (w < hi)
            if seg.any():
                out[seg] = broaden_R(w, f, float(R))[seg]
        return out
    raise ValueError(f'unknown resolution kind: {kind!r}')


def project(obs, w, f):
    """Model on the grid's wavelengths -> the observation's sampling."""
    if obs.filters is not None:
        from common.photometry import project as phot_project
        return phot_project(w, f, obs.filters)
    return np.interp(obs.wavelength, w, f)


def predict(theta, observations, grid):
    """-> [Prediction] in the same order as `observations`.

    theta keys: teff, logg, mh, and optionally ebv, vsini, rv, r_v.
    `rv` is overridden per observation by `obs.rv_fixed` when that is not None,
    which is how XSL (delivered rest-frame) and photometry stay at zero while
    NGSL keeps a free velocity.
    """
    teff, logg, mh = theta['teff'], theta['logg'], theta['mh']
    ebv = float(theta.get('ebv', 0.0))
    vsini = float(theta.get('vsini', 0.0))
    r_v = float(theta.get('r_v', 3.1))

    w0 = grid.wave
    f0 = spectrum_at(grid, teff, logg, mh)
    f0 = broaden_rot(w0, f0, vsini)          # at the star, before everything
    if ebv:
        f0 = redden(w0, f0, ebv, r_v)        # in the ISM

    out = []
    for obs in observations:
        rv = obs.rv_fixed if obs.rv_fixed is not None else float(theta.get('rv', 0.0))
        w = w0 * (1.0 + rv / C_KMS) if rv else w0
        f = instrument(w, f0, obs.resolution)
        out.append(Prediction(obs, project(obs, w, f), wavelength=obs.wavelength))
    return out
