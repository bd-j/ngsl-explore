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

# --- NGSL line spread function -------------------------------------------
#
# TABULATED: FWHM in Angstroms per grating, from the STIS LSF tables
# (FWHM_px x dispersion; data/stis_lsf_resolution.csv). Constant in Angstroms
# within a grating, jumping at the splices. G230LB has no published LSF, so its
# entry is the 2-px sampling lower bound.
NGSL_LSF_TABULATED = [(1675., 3058., 2.75), (3058., 5647., 3.85),
                      (5647., 10198., 8.09)]

# MEASURED: the delivered NGSL v2 spectra are substantially broader than that,
# and broader in a different FUNCTIONAL FORM. Measured against XSL -- which
# observes the same stars at ~10x the resolution, so no model is involved --
# the effective profile is constant in VELOCITY at R = 600 +/- 40 over
# 3900-8700 A, with no jump at the G430L/G750L splice:
#
#     lambda   FWHM     implied R          tabulated
#      3900    6.60 A      591               3.85 A
#      4400    6.99        629               3.85
#      4900    7.79        629               3.85
#      6600   12.68        521               8.09
#      8700   14.48        601               8.09
#
# Constant-R describes this with 7% scatter; constant-Angstrom needs 41%.
# The tabulated values are the single-exposure optical LSF; the delivered
# spectra are co-adds of two dithered exposures resampled onto a common grid,
# which broadens the profile beyond it. Use MEASURED for anything comparing to
# the delivered spectra. See explore/ngsl_lsf_from_xsl.py and docs/DATA.md.
NGSL_R_MEASURED = 600.0

# Back-compatible name; now the measured profile.
NGSL_LSF = NGSL_LSF_TABULATED


def broaden(w, f, fwhm_A, step=0.01):
    """Gaussian of constant FWHM in ANGSTROMS (resample to linear grid first)."""
    wl = np.arange(w[0], w[-1], step)
    sm = gaussian_filter1d(np.interp(wl, w, f), fwhm_A / 2.3548 / step,
                           mode='nearest')
    return np.interp(w, wl, sm)


def broaden_R(w, f, R):
    """Convolve to constant resolving power on a log-lambda grid."""
    step = float(np.median(np.diff(np.log(w))))
    return gaussian_filter1d(f, (1.0 / R) / step / 2.3548, mode='nearest')


def broaden_ngsl(w, f, tabulated=False):
    """The NGSL instrument profile.

    Default is the MEASURED profile: constant R = 600, from matching XSL to
    NGSL for three stars in common (no model involved). Pass tabulated=True for
    the STIS-table profile (constant in Angstroms per grating), which describes
    the single-exposure optics but is ~1.7-1.8x too narrow for the delivered
    co-added spectra.
    """
    if not tabulated:
        return broaden_R(w, f, NGSL_R_MEASURED)
    out = np.array(f, dtype=float)
    for lo, hi, fwhm in NGSL_LSF_TABULATED:
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


def degrade_to(w, f, fwhm_from_A, fwhm_to_A, step=0.01):
    """Convolve a spectrum of resolution `fwhm_from_A` to `fwhm_to_A`.

    The kernel is the QUADRATURE DIFFERENCE, not the target width: convolving
    an already-resolved spectrum with the full target FWHM over-broadens it.
    At the Balmer break XSL is 0.37 A FWHM and NGSL is 3.85 A, so the kernel is
    sqrt(3.85^2 - 0.37^2) = 3.83 A -- a small correction here, but not in the
    VIS where the two are closer.
    """
    if fwhm_to_A <= fwhm_from_A:
        return np.array(f, float)
    k = np.sqrt(fwhm_to_A ** 2 - fwhm_from_A ** 2)
    wl = np.arange(w[0], w[-1], step)
    sm = gaussian_filter1d(np.interp(wl, w, f), k / 2.3548 / step, mode='nearest')
    return np.interp(w, wl, sm)


def rebin_to_pixels(w_in, f_in, w_out):
    """Integrate onto the output pixel grid, conserving flux.

    Not the same as interpolating: a coarse detector pixel averages the flux
    falling across its width, and sampling a high-resolution spectrum at the
    pixel centres instead would keep narrow features a real detector would
    smear out.
    """
    edges = np.empty(len(w_out) + 1)
    edges[1:-1] = 0.5 * (w_out[1:] + w_out[:-1])
    edges[0] = w_out[0] - 0.5 * (w_out[1] - w_out[0])
    edges[-1] = w_out[-1] + 0.5 * (w_out[-1] - w_out[-2])
    csum = np.concatenate([[0.0], np.cumsum(np.diff(w_in) * 0.5 *
                                            (f_in[1:] + f_in[:-1]))])
    lo = np.interp(edges[:-1], w_in, csum)
    hi = np.interp(edges[1:], w_in, csum)
    width = np.diff(edges)
    out = np.where(width > 0, (hi - lo) / width, np.nan)
    out[(edges[:-1] < w_in[0]) | (edges[1:] > w_in[-1])] = np.nan
    return out
