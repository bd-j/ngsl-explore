"""Line spread functions and broadening kernels, shared by the explore
scripts and the fitter.

The NGSL profile is measured, not assumed, and there is ONE function that
applies it: `broaden_ngsl`. See `to_ngsl_pixels` for why the kernel and the
pixel belong together, and docs/LSF.md for the measurement.
"""
from pathlib import Path

import numpy as np
from scipy.ndimage import gaussian_filter1d
from scipy.signal import fftconvolve

ROOT = Path(__file__).resolve().parent.parent
C_KMS = 2.99792458e5

# --- NGSL -----------------------------------------------------------------
#
# (lo, hi, dispersion A/px) per grating, measured from the delivered v2
# wavelength grid. One definition; everything else keys off it.
NGSL_SEGMENTS = [('G230LB', 1675., 3058., 1.373),
                 ('G430L', 3058., 5647., 2.747),
                 ('G750L', 5647., 10198., 4.879)]

# THE TABULATED STIS LSF. Not the applied profile -- it does not describe the
# delivered v2 spectra -- but the loader lives here because three places need
# these files and three parsers of the same columns is how the fitter and the
# figures came to disagree about the instrument in the first place.
#
# The tables give a model LSF for four slit widths. ALL FOUR SHARE A CORE --
# 4.04 A at G430L's midpoint -- and differ only in their wings, because a wide
# slit admits scattered light a narrow one cuts off. That makes them a ruler for
# how much halo NGSL's DELIVERED profile carries, and it is the strongest
# evidence that the halo is real instrumental scattered light rather than a
# fitting artefact. Measured against XSL over nine stars (explore/ngsl_lsf.py,
# G430L, pixel integration):
#
#   profile           npar   rms %   core mean %   FWHM A   power >10 A
#   stis2 52x2.0 fix     0   1.958       -4.74       4.05       8.75
#   stis  52x0.2 fix     0   1.886       +4.74       4.04       0.00
#   gauss                1   1.085       +0.67       6.21       0.01
#   gauss (x) tophat     2   1.085       +0.67       6.21       0.01
#   stis (x) tophat      1   1.057       +0.78       6.44       0.01
#   stis (x) gauss       1   1.047       +0.65       6.03       0.03
#   stis05 52x0.5 fix    0   1.045       -0.11       4.04       2.56
#   MOFFAT               2   0.930       -0.02       3.54       2.44   <- adopted
#   gauss + gauss        3   0.908       +0.08       5.33       3.19
#
# The core residual changes SIGN across the aperture bracket (+4.74% at 52x0.2,
# -4.74% at 52x2.0), so the answer lies strictly between them; 52x2.0 is already
# too broad, since a free Gaussian or box added to it fits ZERO extra
# broadening; and the halo the free Moffat fits (2.44% beyond +/-10 A) is the
# halo 52x0.5 has (2.56%). NGSL observed through 52x0.2, whose tabulated profile
# has no halo at all, so the delivered spectra behave like a slit 2.5x wider
# than the one used -- see docs/LSF.md, which does not claim to explain that.
NGSL_APERTURE = '52x0.2'          # NGSL's slit; the default for stis_table()
NGSL_STIS_ANCHORS = {
    'G430L': {3200.: 'LSF_G430L_3200.txt', 5500.: 'LSF_G430L_5500.txt'},
    'G750L': {7000.: 'LSF_G750L_7000.txt'},
}

# FWHM in Angstroms per grating from those tables (FWHM_px x dispersion,
# data/stis_lsf_resolution.csv). G230LB's entry is a 2-pixel lower bound: the
# only tabulated LSFs at those wavelengths are for the MAMA G230L, a different
# detector.
NGSL_LSF_TABULATED = [(1675., 3058., 2.75), (3058., 5647., 3.85),
                      (5647., 10198., 8.09)]

# THE ADOPTED PROFILE: a Moffat, core constant in ANGSTROMS per grating.
#
# Lowest rms of anything with two parameters or fewer, and the only profile
# whose second parameter is a MEASUREMENT: over the nine stars beta runs 1.40 to
# 2.06, while the two-Gaussian's marginally better fit has a broad component
# scattering from 25 A to 13574 A -- three parameters, two of which mean
# anything.
#
# Constant in ANGSTROMS, not in R. Fitting each sub-window separately and
# regressing width on wavelength as FWHM ~ lambda^alpha over four
# Balmer-anchored G430L windows:
#
#   Moffat core     alpha = +0.14      tabulated STIS   alpha = +0.19
#   single Gaussian alpha = +0.67
#
# A single Gaussian looks nearly constant in R, which is how this project once
# arrived at "R = 600, constant in velocity". That was an artefact of the
# functional form: a one-parameter profile standing in for a core plus a halo
# drifts with wavelength as the halo's relative weight changes.
# NGSL_R_MEASURED = 600 has been REMOVED rather than kept for reference,
# because it is neither a width nor a resolution.
NGSL_MOFFAT_BETA = 1.52           # beta->1 Lorentzian, beta->inf Gaussian
NGSL_BETA_SIGMA = 0.16            # star-to-star NMAD, 18 star-grating fits

# (lo, hi, core FWHM in A, dispersion in A/px). Core fitted UNDER PIXEL
# INTEGRATION over nine stars -- a width is meaningless without that convention,
# and the same data fitted with centre sampling returns 4.06 A for G430L rather
# than 3.54.
#
#   G430L  3.54 +/- 0.16 A   (9 stars)
#   G750L  8.38 +/- 0.45 A   (9 stars)
#
# G430L is measured over 3700-5647 A -- XSL starts at 3501 A -- so the bluest
# band this project uses (3220-3385 A) is an extrapolation, worth -2.4% in the
# core width at alpha = +0.14, inside the star-to-star scatter. G230LB is not
# measured and keeps a 2-pixel placeholder; nothing uses it.
NGSL_LSF_CORE = [(1675., 3058., 1.37, 1.373),
                 (3058., 5647., 3.54, 2.747),
                 (5647., 10198., 8.38, 4.879)]

# Kernel truncation, as a half-width in DETECTOR PIXELS rather than in FWHM.
#
# Pixels because that is the unit the instrument works in, and because it keeps
# the kernel's reach comparable to the tabulated STIS profiles, which stop at
# +/-20 px. A Moffat has no natural edge, so where it is cut is a choice, and
# cutting it in FWHM made the G750L kernel reach 335 A while G430L's reached
# 142 A -- the same profile behaving differently in the two gratings for no
# instrumental reason.
#
# 15 px is 41 A on G430L and 73 A on G750L. It keeps the profile's effect LOCAL,
# which is the point: the empirical STIS LSFs are compact, and a kernel that
# reaches 300 A is making a claim about scattered light at a distance that
# nothing here measures. The cost is small and was measured rather than assumed.
# Convolving a line-rich spectrum on one uniform grid, no segment edges
# involved, and comparing against a +/-300 A (109 px) kernel over an interior
# window:
#
#   truncation        max rel error   rms rel error
#   +/-10 px  27 A       2.4e-04         5.8e-05
#   +/-15 px  41 A       9.8e-05         2.6e-05      <- adopted
#   +/-20 px  55 A       5.8e-05         1.4e-05
#   +/-40 px 110 A       1.2e-05         3.2e-06
#
# At +/-15 px the profile is 1.6e-4 of peak on G430L and 3.8e-4 on G750L, and
# the power discarded is 0.14% and 0.25% -- renormalised away by moffat_kernel,
# so no flux is lost, only reach. An error of 1e-4 is two orders below the 1%
# calibration floor the NGSL bands carry, so this is a free choice made on
# physical grounds rather than a trade.
NGSL_TRUNC_PX = 15.0


# --- the tabulated STIS profile ----# --- the tabulated STIS profile -------------------------------------------

_STIS_CACHE = {}


def stis_table(fname, aperture=NGSL_APERTURE):
    """-> (rel_pixel, response) for one aperture column, area-normalised.

    ONE reader for these files. `common.lsf`, `explore/ngsl_lsf.py` and
    `explore/lsf_resolution.py` all need them, and three parsers of the same
    columns is how the fitter and the figures came to disagree about the
    instrument in the first place.
    """
    key = (fname, aperture)
    if key not in _STIS_CACHE:
        lines = (ROOT / 'data' / 'stis_lsf' / fname).read_text().splitlines()
        # header is 'Rel pixel  52x0.1  52x0.2 ...', so the aperture columns
        # start at data column 1.
        cols = lines[1].split()[2:]
        j = cols.index(aperture) + 1
        d = np.array([[float(v) for v in l.split()] for l in lines[2:]
                      if l.strip()])
        x, y = d[:, 0], np.clip(d[:, j], 0.0, None)
        _STIS_CACHE[key] = (x, y / np.trapezoid(y, x))
    return _STIS_CACHE[key]


def stis_kernel(anchors, lam, disp, dl, aperture=NGSL_APERTURE):
    """The tabulated LSF at `lam`, on the offset grid `dl` (Angstroms).

    `anchors` is {anchor wavelength: filename}. Interpolated linearly in
    wavelength between them and held constant outside, and scaled from the
    table's relative-pixel abscissa to Angstroms by the grating dispersion.

    """
    ws = sorted(anchors)
    ys = [np.interp(dl / disp, *stis_table(anchors[w0], aperture),
                    left=0.0, right=0.0) for w0 in ws]
    if len(ws) == 1:
        k = ys[0]
    else:
        i = int(np.clip(np.searchsorted(ws, lam) - 1, 0, len(ws) - 2))
        t = float(np.clip((lam - ws[i]) / (ws[i + 1] - ws[i]), 0.0, 1.0))
        k = (1 - t) * ys[i] + t * ys[i + 1]
    tot = k.sum()
    return k / tot if tot > 0 else k


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


def broaden_moffat(w, f, fwhm_A, half_A, beta=NGSL_MOFFAT_BETA, step=0.05):
    """Moffat of constant FWHM in ANGSTROMS (resample to a linear grid first).

    `half_A` is the kernel half-width in ANGSTROMS, and it is REQUIRED rather
    than defaulted. A Moffat has no natural edge, so where it is cut is a
    choice with consequences -- truncating early discards the very power that
    distinguishes it from a Gaussian -- and a caller that has not thought about
    it should be made to. `broaden_ngsl` derives it from NGSL_TRUNC_PX and the
    grating dispersion. The kernel is renormalised after truncation, so no flux
    is lost, only reach.

    `step` trades accuracy against speed, and the trade is steeper than it
    looks. Measured against a 0.02 A reference, the MAX relative error (which
    sits at the sharpest line cores) runs 1.6e-3 at 0.05 A, 7.4e-3 at 0.1 and
    6.3e-2 at 0.4 -- so a step chosen to "oversample the 4 A kernel" is not good
    enough, because what has to be resolved is the model's own line cores, not
    the kernel. 0.05 A it is.

    The convolution is FFT-based for the same reason: at 0.05 A a direct
    convolution with this kernel cost 4.4 s per model evaluation, which is 25
    hours for the node scan.

    Edge handling replicates the end values, matching
    gaussian_filter1d(mode='nearest'). It matters more here than for a Gaussian:
    this kernel is hundreds of points wide, and zero-padding would darken the
    blue end of the grid -- exactly where the 3220-3385 A band sits, the longest
    lever the fit has on reddening.
    """
    wl = np.arange(w[0], w[-1], step)
    y = np.interp(wl, w, f)
    n = int(np.ceil(half_A / step))
    k = moffat_kernel(n, fwhm_A / step, beta)
    pad = np.concatenate([np.full(n, y[0]), y, np.full(n, y[-1])])
    sm = fftconvolve(pad, k, mode='same')[n:n + y.size]
    return np.interp(w, wl, sm)


def broaden_ngsl(w, f, beta=NGSL_MOFFAT_BETA, segments=None,
                 trunc_px=NGSL_TRUNC_PX):
    """THE NGSL line spread function. Apply this and nothing else.

    A Moffat of core FWHM constant in Angstroms within each grating, jumping at
    the splices -- the measured profile, see NGSL_LSF_CORE above and
    docs/LSF.md. This is the only NGSL kernel in the project: there is no
    Gaussian alternative and no `tabulated=True` switch, because having two
    live broadening paths is what previously let the fitter and the comparison
    figures disagree about the instrument. To compare against a tabulated STIS
    profile, use `stis_kernel` explicitly and say so.

    IT MUST BE FOLLOWED BY PIXEL INTEGRATION. The widths were fitted that way;
    applying this kernel and then sampling at pixel centres under-smooths the
    model by exactly the pixel, which is a real error with no symptom except a
    residual. `to_ngsl_pixels` does both together and is what callers should
    normally use.

    The kernel is truncated at `trunc_px` DETECTOR PIXELS either side, so its
    reach scales with the grating's dispersion rather than with its own width.
    Each segment is convolved over its own range plus that margin, rather than
    over the whole spectrum and then masked -- the latter did the full
    convolution once per grating for every model evaluation, which the node scan
    does 61 times per node.
    """
    w = np.asarray(w, float)
    f = np.asarray(f, float)
    out = np.array(f, dtype=float)
    for lo, hi, fwhm, disp in (segments or NGSL_LSF_CORE):
        seg = (w >= lo) & (w < hi)
        if not seg.any():
            continue
        half_A = trunc_px * disp
        sub = (w >= lo - half_A) & (w < hi + half_A)
        sm = broaden_moffat(w[sub], f[sub], fwhm, half_A, beta)
        out[seg] = sm[seg[sub]]
    return out


def to_ngsl_pixels(w, f, w_out, **kw):
    """Model spectrum -> NGSL pixels: the measured LSF and the pixel, together.

    The two are a matched pair. `broaden_ngsl`'s widths were fitted against XSL
    with the comparison spectrum INTEGRATED onto NGSL pixels, so a caller that
    applies the kernel and then interpolates has silently changed the profile.
    Keeping both inside one function is the point: it is the one call that
    cannot be got half right.
    """
    return rebin_to_pixels(w, broaden_ngsl(w, f, **kw), w_out)


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
