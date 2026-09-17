"""Line spread functions and broadening kernels, shared by the explore
scripts and the fitter.

The NGSL profile is measured, not assumed, and there is ONE function that
applies it: `broaden_ngsl`. See `to_ngsl_pixels` for why the kernel and the
pixel belong together, and docs/LSF.md for the measurement.
"""
import numpy as np
from scipy.ndimage import gaussian_filter1d
from scipy.signal import fftconvolve

C_KMS = 2.99792458e5

# --- NGSL -----------------------------------------------------------------
#
# (lo, hi, dispersion A/px) per grating, measured from the delivered v2
# wavelength grid. One definition; everything else keys off it.
NGSL_SEGMENTS = [('G230LB', 1675., 3058., 1.373),
                 ('G430L', 3058., 5647., 2.747),
                 ('G750L', 5647., 10198., 4.879)]

# The TABULATED STIS LSF, FWHM in Angstroms per grating (FWHM_px x dispersion,
# data/stis_lsf_resolution.csv). Kept for comparison only -- it is NOT what the
# delivered spectra have. Fitted as a fixed profile against XSL it leaves an
# rms of 1.89% and a +14% spike in every Balmer core, against 0.93% and +5.9%
# for the profile adopted below. G230LB's entry is a 2-pixel lower bound: the
# only tabulated LSFs at those wavelengths are for the MAMA G230L, a different
# detector.
NGSL_LSF_TABULATED = [(1675., 3058., 2.75), (3058., 5647., 3.85),
                      (5647., 10198., 8.09)]

# THE MEASURED PROFILE: a Moffat, core constant in ANGSTROMS per grating.
#
# Measured by `explore/ngsl_lsf.py` against XSL -- two observations of the same
# nine stars, no stellar model anywhere in it -- fitting seven profile families
# with free widths and a free wavelength shift, scored on the rms percent
# residual over the whole grating. Median over the nine primary stars, G430L:
#
#   profile        npar   rms %   core mean %   FWHM A   power >10 A
#   stis (fixed)      0   1.886       +4.74       4.04       0.00
#   gauss             1   1.085       +0.67       6.21       0.01
#   gauss (x) tophat  2   1.085       +0.67       6.21       0.01
#   stis (x) tophat   1   1.057       +0.78       6.44       0.01
#   stis (x) gauss    1   1.047       +0.65       6.03       0.03
#   MOFFAT            2   0.930       -0.02       3.54       2.44
#   gauss + gauss     3   0.908       +0.08       5.33       3.19
#
# Why the Moffat and not the marginally better two-Gaussian: its second
# parameter is a MEASUREMENT. Over the nine stars beta runs 1.40 to 2.06, while
# the two-Gaussian's broad component scatters from 25 A to 13574 A -- that
# profile has three parameters and only two of them mean anything.
#
# Why constant in ANGSTROMS. Fitting each sub-window separately and regressing
# the width on wavelength as FWHM ~ lambda^alpha (alpha = 0 constant in A,
# alpha = 1 constant in R), over four Balmer-anchored G430L windows:
#
#   Moffat core     alpha = +0.14      tabulated STIS   alpha = +0.19
#   single Gaussian alpha = +0.67
#
# The Moffat core is flat and tracks the tables' own mild wavelength
# dependence. A single Gaussian looks nearly constant in R -- which is exactly
# how this project once arrived at "R = 600, constant in velocity". That was an
# artefact of the functional form: a one-parameter profile forced to represent
# a core plus a halo drifts with wavelength as the halo's relative weight
# changes. NGSL_R_MEASURED = 600 has been REMOVED rather than kept for
# reference, because it is neither a width nor a resolution.
NGSL_MOFFAT_BETA = 1.52           # beta->1 Lorentzian, beta->inf Gaussian
NGSL_BETA_SIGMA = 0.16            # star-to-star NMAD, 18 star-grating fits

# Core FWHM in Angstroms per grating, fitted UNDER PIXEL INTEGRATION. A width
# is meaningless without that convention: the same data fitted with centre
# sampling returns 4.06 A for G430L, the 3.54 A below with the 2.747 A pixel's
# equivalent Gaussian added in quadrature.
#
#   G430L  3.54 +/- 0.16 A   (9 stars)
#   G750L  8.38 +/- 0.45 A   (9 stars)
#
# G430L is measured over 3700-5647 A -- XSL starts at 3501 A, so the bluest
# band this project uses (3220-3385 A) is an extrapolation. At alpha = +0.14
# that extrapolation is worth -2.4% in the core width, which is inside the
# star-to-star scatter. G230LB is not measured at all and keeps a 2-pixel
# placeholder; nothing in this project uses it (the model grid starts at
# 3200 A).
NGSL_LSF_CORE = [(1675., 3058., 1.37), (3058., 5647., 3.54),
                 (5647., 10198., 8.38)]


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

    `ntrunc` is the half-width in FWHM. A beta ~ 1.5 Moffat has heavy tails, so
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


def broaden_ngsl(w, f, beta=NGSL_MOFFAT_BETA, segments=None, ntrunc=40.0):
    """THE NGSL line spread function. Apply this and nothing else.

    A Moffat of core FWHM constant in Angstroms within each grating, jumping at
    the splices -- the measured profile, see NGSL_LSF_CORE above and
    docs/LSF.md. This is the only NGSL kernel in the project: there is no
    Gaussian alternative and no `tabulated=True` switch, because having two
    live broadening paths is what previously let the fitter and the comparison
    figures disagree about the instrument. To compare against the tabulated
    profile, pass NGSL_LSF_TABULATED as `segments` explicitly and say so.

    IT MUST BE FOLLOWED BY PIXEL INTEGRATION. The widths were fitted that way;
    applying this kernel and then sampling at pixel centres under-smooths the
    model by exactly the pixel, which is a real error with no symptom except a
    residual. `to_ngsl_pixels` does both together and is what callers should
    normally use.

    Each segment is convolved over its own range plus a kernel half-width of
    margin, rather than over the whole spectrum and then masked -- the latter
    did the full convolution once per grating for every model evaluation, which
    the node scan does 61 times per node.
    """
    w = np.asarray(w, float)
    f = np.asarray(f, float)
    out = np.array(f, dtype=float)
    for lo, hi, fwhm in (segments or NGSL_LSF_CORE):
        seg = (w >= lo) & (w < hi)
        if not seg.any():
            continue
        pad = ntrunc * fwhm
        sub = (w >= lo - pad) & (w < hi + pad)
        sm = broaden_moffat(w[sub], f[sub], fwhm, beta, ntrunc=ntrunc)
        out[seg] = sm[seg[sub]]
    return out


def to_ngsl_pixels(w, f, w_out, beta=NGSL_MOFFAT_BETA, segments=None,
                   ntrunc=40.0):
    """Model spectrum -> NGSL pixels: the measured LSF and the pixel, together.

    The two are a matched pair. `broaden_ngsl`'s widths were fitted against XSL
    with the comparison spectrum INTEGRATED onto NGSL pixels, so a caller that
    applies the kernel and then interpolates has silently changed the profile.
    Keeping both inside one function is the point: it is the one call that
    cannot be got half right.
    """
    return rebin_to_pixels(w, broaden_ngsl(w, f, beta, segments, ntrunc), w_out)


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
