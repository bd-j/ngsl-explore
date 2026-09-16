"""Line spread functions and broadening kernels, shared by the explore
scripts and the fitter.

The NGSL LSF is set by a fixed dispersion per grating, so it is constant in
ANGSTROMS within a grating and jumps at the splices. SYNTHE's output grid is
logarithmic, so a fixed sigma in pixels would instead be a constant-R kernel --
the wrong thing here.
"""
import numpy as np
from scipy.ndimage import gaussian_filter1d
from scipy.signal import fftconvolve

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
NGSL_R_SIGMA = 40.0      # star-to-star scatter of the measurement

# MEASURED SHAPE. The single Gaussian above is the wrong FUNCTIONAL FORM, not
# just the wrong width. Fitting a family of profiles against XSL -- no model
# involved -- on a control window with no Balmer line in it, and scoring each on
# the Balmer cores it was NOT fitted to (explore/ngsl_lsf_shape.py):
#
#   profile          n par   control rms   leftover Balmer core excess
#   Gaussian             1        0.0068        +2.33%
#   Gaussian * tophat    2        0.0067        +2.32%
#   Gaussian + Gaussian  3        0.0059        +0.75%
#   Gaussian + Lorentz   3        0.0058        +0.27%
#   MOFFAT               2        0.0058        +0.04%
#
# The tophat is the physically obvious candidate -- NGSL v2 spectra are co-adds
# of two DITHERED exposures resampled onto a common grid, and both the dither
# and the pixel are boxes. It fits a sensible 1.84 A box (~1.3 pixels) and
# changes nothing, because a box convolved with a Gaussian still has
# Gaussian-fast wings. The pedestal is not resampling; it is a heavy-tailed
# halo, of the kind grating scatter produces. Power beyond +/-10 A: Gaussian
# 0.01%, Gaussian*tophat 0.01%, Moffat 2.9%.
#
# AND THE CORE IS THE TABULATED ONE. Fitted per grating, the Moffat core comes
# out 4.02 +/- 0.59 A for G430L and 8.34 +/- 0.91 A for G750L, against the STIS
# tabulated 3.85 and 8.09 -- agreement to 3-5%. That resolves the disagreement
# explore/ngsl_lsf_from_xsl.py records as unexplained: a single Gaussian needed
# 1.7-1.9x the tabulated width because it was absorbing a tail it had no way to
# represent. So the profile is the PUBLISHED STIS core plus a halo, and the
# width is constant in ANGSTROMS per grating, as the tables say -- NOT constant
# in R. Fitting R = 1074 from one window and applying it as constant-R made the
# kernel sharper at 3800 A than anything that had been tested, and made the
# Balmer core excess worse rather than better.
# WHAT THE BALMER CORES CAN AND CANNOT SETTLE. An earlier version of this note
# claimed a Moffat removed 86% of the Balmer core excess. That was wrong: the
# core excess is degenerate with the EFFECTIVE WIDTH, not the shape. Holding
# everything else fixed and varying only the width, a plain Gaussian runs from
# +16.8% at 3.85 A to -1.2% at 7.0 A. Any profile can be tuned to zero it. So
# the cores are not evidence for this profile and are not used as such; the
# evidence is the control-window rms at free width, and the agreement of the
# fitted core with the independently published STIS value.
NGSL_MOFFAT_BETA = 1.6            # beta->1 Lorentzian, beta->inf Gaussian
# Core FWHM fitted against XSL per grating, NOT tuned on the Balmer region:
# 4.02 +/- 0.59 A (G430L) and 8.34 +/- 0.91 A (G750L). G230LB is not fitted --
# no sample star has XSL coverage below 3500 A -- so it keeps its tabulated
# 2-pixel lower bound, and nothing in this project uses it (the grid starts at
# 3200 A and the bluest band at 3220 A).
NGSL_MOFFAT_LSF = [(1675., 3058., 2.75), (3058., 5647., 4.02),
                   (5647., 10198., 8.34)]

# Back-compatible name; now the measured profile.
NGSL_LSF = NGSL_LSF_TABULATED


def broaden(w, f, fwhm_A, step=0.01):
    """Gaussian of constant FWHM in ANGSTROMS (resample to linear grid first)."""
    wl = np.arange(w[0], w[-1], step)
    sm = gaussian_filter1d(np.interp(wl, w, f), fwhm_A / 2.3548 / step,
                           mode='nearest')
    return np.interp(w, wl, sm)


def broaden_R(w, f, R):
    """Convolve to constant resolving power on a log-lambda grid (Gaussian)."""
    step = float(np.median(np.diff(np.log(w))))
    return gaussian_filter1d(f, (1.0 / R) / step / 2.3548, mode='nearest')


def moffat_kernel(n_pix, fwhm_pix, beta):
    """Normalised Moffat on a pixel grid, truncated at n_pix half-width."""
    a = fwhm_pix / (2.0 * np.sqrt(2.0 ** (1.0 / beta) - 1.0))
    x = np.arange(-n_pix, n_pix + 1, dtype=float)
    k = (1.0 + (x / a) ** 2) ** (-beta)
    return k / k.sum()


def broaden_moffat(w, f, fwhm_A, beta=NGSL_MOFFAT_BETA, step=0.05, ntrunc=40.0):
    """Moffat of constant FWHM in ANGSTROMS (resample to a linear grid first).

    `ntrunc` is the half-width in FWHM. A beta ~ 1.6 Moffat has heavy tails, so
    truncating early discards the very power that distinguishes it from a
    Gaussian; the kernel is renormalised after truncation so no flux is lost.

    `step` trades accuracy against speed, and the trade is steeper than it
    looks. Measured against a 0.02 A reference, the MAX relative error (which
    sits at the sharpest line cores) runs 1.6e-3 at 0.05 A, 7.4e-3 at 0.1 and
    6.3e-2 at 0.4 -- so a step chosen to "oversample the 4 A kernel" is not good
    enough, because what has to be resolved is the model's own line cores, not
    the kernel. 0.05 A it is, at 3.7 ms per call.

    The convolution is FFT-based for the same reason: at 0.05 A a direct
    convolution with this kernel cost 4.4 s per model evaluation, which is 25
    hours for the node scan.

    Edge handling replicates the end values, matching
    gaussian_filter1d(mode='nearest'). It matters more here than for a Gaussian:
    this kernel is hundreds of pixels wide, and zero-padding would darken the
    blue end of the grid -- exactly where the 3220-3385 A band sits, the longest
    lever the fit has on reddening.
    """
    wl = np.arange(w[0], w[-1], step)
    y = np.interp(wl, w, f)
    n = int(np.ceil(ntrunc * fwhm_A / step))
    k = moffat_kernel(n, fwhm_A / step, beta)
    pad = np.concatenate([np.full(n, y[0]), y, np.full(n, y[-1])])
    sm = fftconvolve(pad, k, mode='same')[n:n + y.size]
    return np.interp(w, wl, sm)


def broaden_ngsl_moffat(w, f, beta=NGSL_MOFFAT_BETA, segments=None, ntrunc=40.0):
    """The measured NGSL profile: STIS core per grating plus a Moffat tail.

    Each segment is convolved over its own range plus a kernel half-width of
    margin, rather than over the whole spectrum and then masked -- the latter
    did the full convolution once per grating for every model evaluation, which
    the node scan does 61 times per node.
    """
    w = np.asarray(w, float)
    out = np.array(f, dtype=float)
    for lo, hi, fwhm in (segments or NGSL_MOFFAT_LSF):
        seg = (w >= lo) & (w < hi)
        if not seg.any():
            continue
        pad = ntrunc * fwhm
        sub = (w >= lo - pad) & (w < hi + pad)
        sm = broaden_moffat(w[sub], np.asarray(f, float)[sub], fwhm, beta,
                            ntrunc=ntrunc)
        out[seg] = sm[seg[sub]]
    return out


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
