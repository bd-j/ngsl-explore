"""Line spread functions and broadening kernels, shared by the explore
scripts and the fitter.

The NGSL LSF is set by a fixed dispersion per grating, so it is constant in
ANGSTROMS within a grating and jumps at the splices. SYNTHE's output grid is
logarithmic, so a fixed sigma in pixels would instead be a constant-R kernel --
the wrong thing here.
"""
import numpy as np
from scipy.ndimage import gaussian_filter1d

C_KMS = 2.99792458e5

# FWHM in Angstroms per grating (FWHM_px x dispersion, data/stis_lsf_resolution.csv).
# G230LB has no published LSF: the value is the 2-px sampling lower bound.
NGSL_LSF = [(1675., 3058., 2.75), (3058., 5647., 3.85), (5647., 10198., 8.09)]


def broaden(w, f, fwhm_A, step=0.01):
    """Gaussian of constant FWHM in ANGSTROMS (resample to linear grid first)."""
    wl = np.arange(w[0], w[-1], step)
    sm = gaussian_filter1d(np.interp(wl, w, f), fwhm_A / 2.3548 / step,
                           mode='nearest')
    return np.interp(w, wl, sm)


def broaden_ngsl(w, f):
    """The NGSL instrument profile: piecewise constant in Angstroms per grating."""
    out = np.array(f, dtype=float)
    for lo, hi, fwhm in NGSL_LSF:
        seg = (w >= lo) & (w < hi)
        if seg.any():
            out[seg] = broaden(w, f, fwhm)[seg]
    return out


def rot_kernel(dl, lam0, vsini, eps=0.6):
    """Gray (2005) rotational broadening profile, linear limb darkening."""
    dlL = lam0 * vsini / C_KMS
    x = dl / dlL
    k = np.zeros_like(x)
    ok = np.abs(x) < 1
    k[ok] = (2 * (1 - eps) * np.sqrt(1 - x[ok] ** 2)
             + 0.5 * np.pi * eps * (1 - x[ok] ** 2)) / (np.pi * dlL * (1 - eps / 3))
    return k / k.sum()


def broaden_rot(w, f, vsini, eps=0.6):
    """Apply rotational broadening on a log-lambda (constant-R) grid."""
    if not vsini or vsini <= 0:
        return f
    dl = float(np.median(np.diff(w)))
    n = int(np.ceil(float(np.median(w)) * vsini / C_KMS / dl)) * 2 + 1
    g = (np.arange(n) - n // 2) * dl
    return np.convolve(f, rot_kernel(g, float(np.median(w)), vsini), mode='same')
