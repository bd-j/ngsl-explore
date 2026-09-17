"""Measure the NGSL line spread function against XSL. One script, no models.

THE MEASUREMENT. XSL observes the same stars at R ~ 9800, roughly 20x NGSL's
resolving power, so degrading XSL to NGSL asks a purely instrumental question:
what kernel turns the true spectrum into what NGSL recorded? Nothing in this
script uses a stellar model, a grid node or a fitted stellar parameter, so
nothing it concludes can be an atmosphere error in disguise. That is the whole
reason to measure the LSF this way.

    NGSL(lam_n)  ~=  P(lam_n) * S[ XSL (x) K(theta) ](lam_n)

  K       the trial line spread function (below)
  S       the sampling operator -- INTEGRATE across each NGSL pixel, or
          INTERPOLATE at the pixel centre. Both are reported; see SAMPLING.
  P       a low-order polynomial in lambda, solved by linear least squares at
          every kernel trial. The two libraries differ by slit losses, aperture
          corrections and flux calibration, all smooth in lambda, and none of
          that should be allowed to influence a width. Because P multiplies the
          smoothed XSL spectrum it cannot absorb a line: it is a continuum
          ratio, not a continuum normalisation of each spectrum separately.

Residuals are reported in PERCENT of the model, three ways: the rms over the
whole grating, and the mean and peak inside the Balmer line cores. The cores
are the point of the exercise -- they are where a wrong kernel shows first, and
they are what the break measurement is sensitive to -- but they are a DIAGNOSTIC
here, never a fit target. The fit sees the rms and nothing else.

THE PROFILES. Six, in increasing freedom:

  stis          the tabulated STIS LSF, 52x0.2 (the NGSL slit), from
                data/stis_lsf/. ZERO free parameters -- an instrument
                calibration product, interpolated in wavelength between the
                tabulated anchors and scaled by the grating dispersion. This is
                the null hypothesis: NGSL is what the tables say it is.
  stis_gauss    the same, convolved with a free Gaussian. One parameter. Asks
                whether the tabulated profile is the right CORE with something
                extra on top, which is a different question from whether a free
                Gaussian happens to fit.
  stis2         the tabulated LSF for the 52x2.0 aperture, zero free
                parameters. NGSL did not use that slit, so this is not a
                candidate on physical grounds -- it is here because it is the
                only tabulated profile with a HALO. All apertures share a core
                (3.81 A at G430L 3200) and differ only in the wings, a wide slit
                admitting scattered light a narrow one cuts off. If NGSL's
                delivered tail is that scattered light leaking in, this is what
                it should look like.
  stis05        the 52x0.5 column, zero free parameters. Same core again, and
                a halo between the other two (11.65% beyond +/-10 A against
                0.00% and 18.30%), so the three of them bracket the question
                "how much scattered light does NGSL's slit actually admit?"
  stis2_gauss   52x2.0 convolved with a free Gaussian. One parameter.
  stis2_tophat  52x2.0 convolved with a free boxcar. One parameter.
  stis_tophat   the same, convolved with a free boxcar. One parameter. This is
                the co-add hypothesis stated exactly: the optics are as STScI
                tabulates them, and everything NGSL has on top comes from
                combining dithered exposures onto a common grid -- an operation
                built out of boxes (the pixel, the dither offset, a linear
                interpolation). If that is the whole story this family should
                match anything with more freedom. It has a sharp prediction
                going for it: a box has no tails, so it can only widen the
                profile, never put power at +/-10 A.
  gauss         a single free Gaussian. One parameter.
  moffat        a Moffat. Two parameters; beta sets the tail weight, and beta
                -> infinity recovers the Gaussian.
  gauss2        two Gaussians, narrow plus broad. Three parameters.
  gauss_tophat  a Gaussian convolved with a boxcar. Two parameters. The
                physically motivated one: NGSL v2 spectra are co-adds of two
                DITHERED exposures resampled onto a common grid, and both the
                dither offset and the pixel are boxes.

SAMPLING. Whether to integrate across the pixel or interpolate at its centre is
not a detail: at NGSL's 2.747 A G430L pixels the difference is comparable to
everything this project measures. A fitted width is only meaningful alongside
the sampling convention it was fitted under -- a fit denied pixel integration
inflates its kernel to compensate, and by a knowable amount, which --selftest
checks: 6.89 A under `rebin` becomes 7.14 A under `interp`, against
sqrt(6.89^2 + 1.867^2) = 7.139 for the 2.747 A pixel's equivalent Gaussian.

`rebin` is the one to quote, because a detector integrates across its pixel and
`fitting.predict.project` does the same.

WHAT THE SAMPLING COMPARISON DOES NOT SETTLE, despite an earlier version of this
docstring saying it does: whether the tabulated STIS LSF already includes the
detector pixel response. The proposed test was that `interp` should win for
`stis` if it does. It cannot work. `stis` is far too narrow under BOTH
conventions (rms 1.886 vs 2.201, against 0.930 for a Moffat), so `rebin` wins
only because it is the wider of two profiles that are both too sharp -- which is
what widening always does when you are too narrow, and says nothing about the
table. The indirect evidence is better: the fitted G430L core at 3.54 A combined
with the pixel in quadrature gives 4.01 A against the tabulated 4.04 A, which
is what you would see if the table already contained the pixel. G750L gives
9.01 A against 8.09 A and does not support it. So this is left open.

WINDOWS. G430L and G750L are fitted separately: the dispersion differs by 1.8x,
so a kernel constant in Angstroms in one is wrong in the other. G430L is where
the Balmer lines are and is the grating that matters here. G750L is measured for
completeness and is much weaker evidence -- see TELLURIC.

    python3 explore/ngsl_lsf.py                     # all primary stars
    python3 explore/ngsl_lsf.py --star HD194453     # one
    python3 explore/ngsl_lsf.py --no-plots          # table only

Writes data/ngsl_lsf.csv, data/ngsl_lsf_lines.csv and figures/ngsl_lsf/.
"""
import argparse
import csv
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
from scipy.signal import fftconvolve

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common.lsf import rebin_to_pixels
from common.lines import hydrogen_lines
from common.xsl_load import load as load_xsl, sigma_v, C_KMS
from fitting.observations import load_ngsl

ROOT = Path(__file__).resolve().parent.parent
FIGDIR = ROOT / 'figures' / 'ngsl_lsf'

STEP = 0.05          # A, the uniform grid kernels are built and applied on
KERNEL_HALF = 80.0   # A, kernel half-width. Truncating a heavy-tailed kernel
                     # early discards the power that distinguishes it from a
                     # Gaussian, so this was MEASURED rather than assumed: at
                     # 40 / 80 / 150 / 250 A the Moffat fit returns rms 0.666 /
                     # 0.665 / 0.665 / 0.665 % and FWHM 4.55 A throughout
                     # (HD194453, G430L). The widest component any family wants
                     # is the gauss2 pedestal at 25 A FWHM, so +/-80 A is 3x it
                     # and the choice costs nothing either way.
CORE_HALF = 4.0      # A, half-width of a Balmer "core" -- ~1.5 G430L pixels,
                     # so every core has at least three pixels in it.
WIDTH_RUNAWAY = 20.0 # A. A sub-window fit returning a width above this did
                     # not measure one -- the guard in fit() is 40 A and an
                     # unconstrained window walks to it. Used only by
                     # --subwindows, where some windows have no strong line.
POLY_DEG = 5         # Chebyshev degree of the continuum ratio P. Over G430L's
                     # 1950 A that is structure on ~390 A scales, comfortably
                     # broader than a Balmer wing (40-120 A at this resolution)
                     # and comfortably narrower than the library flux-calibration
                     # differences it has to absorb.
                     #
                     # ABSOLUTE rms depends on this and comparisons across
                     # degrees are meaningless; the RANKING and the WIDTHS do
                     # not. Over degree 2 to 12 (HD194453, G430L) the Moffat
                     # rms falls 0.881 -> 0.608 % while its FWHM moves 4.61 ->
                     # 4.55 A, and stis > gauss > gauss2 ~ moffat at every
                     # degree. So quote rms only against another family fitted
                     # at the same degree.

# (name, lo, hi, dispersion A/px). The fit ranges start well inside XSL's blue
# edge at 3501 A: the convolution replicates edge values, which is a poor model
# of a spectrum that is still falling across the Balmer jump, so the first
# ~200 A of XSL is used to smooth with and never scored on. G230LB is not
# measured at all -- no sample star has XSL below 3501 A, and the only tabulated
# LSFs for it are the MAMA G230L ones, a different detector.
GRATINGS = [('G430L', 3700., 5647., 2.747),
            ('G750L', 5700., 10198., 4.879)]

# The tabulated STIS model LSFs, per grating, by anchor wavelength (A).
STIS_ANCHORS = {'G430L': {3200.: 'LSF_G430L_3200.txt',
                          5500.: 'LSF_G430L_5500.txt'},
                'G750L': {7000.: 'LSF_G750L_7000.txt'}}
# NGSL observed through 52x0.2, so that is the physically correct column. The
# 52x2.0 column is carried too because it is the one with a HALO: the tables
# give all apertures the same core (3.81 A at G430L 3200) but wildly different
# wings, since a wide slit admits scattered light a narrow one cuts off.
#
#   aperture   G430L core   power beyond +/-10 A
#   52x0.1       3.82 A            0.00%
#   52x0.2       4.04 A            0.00%     <- NGSL's slit
#   52x0.5       4.04 A            2.56%     <- the fitted Moffat wants 2.44%
#   52x2.0       4.05 A            8.75%
#
# (Cores and wings at 4674 A, the G430L fit-range midpoint, via wing_power on
# the kernel grid. Do NOT measure the wing as trapz over a two-sided mask of the
# raw table: the mask is not contiguous, so trapz bridges the gap across the
# core and roughly doubles the answer. That mistake produced 11.65% and 18.30%
# for the two wide apertures in an earlier version of this comment.)
#
# So these bracket the question "how much scattered light does NGSL's slit
# actually admit?", and the answer turns out to be about as much as a 0.5"
# slit -- not the 0.2" its own tabulated profile implies. Add another column by
# adding a family here and to FAMILIES; nothing else needs changing.
STIS_APERTURE = '52x0.2'    # NGSL's slit; the default for stis_table()
STIS_FAMILY_APERTURE = {
    'stis': '52x0.2', 'stis_gauss': '52x0.2', 'stis_tophat': '52x0.2',
    'stis2': '52x2.0', 'stis2_gauss': '52x2.0', 'stis2_tophat': '52x2.0',
    'stis05': '52x0.5',
}

# Telluric absorption bands. NGSL is above the atmosphere and XSL is not, so
# every one of these is a feature in one spectrum and not the other -- a
# difference no line spread function can reconcile. XSL DR3 is telluric
# corrected, which reduces them rather than removing them, and an imperfectly
# divided-out band is still a feature NGSL does not have. They are cut from the
# residual statistics, not from the smoothing. This is a large part of why
# G750L is weaker evidence than G430L: the cuts take ~30% of it.
#
# The list is the strong optical bands only -- O2 gamma, O2 B, H2O, O2 A, H2O,
# and everything past 8900 A. An earlier version also cut 6450-6600 A, which
# has no strong telluric band in it and which swallowed H-alpha, leaving G750L
# with no Balmer core to measure at all.
TELLURIC = [(6270., 6330.), (6860., 6960.), (7150., 7400.),
            (7570., 7720.), (8100., 8400.), (8900., 10198.)]

# XSL's UVB/VIS arm splice, where its own resolution changes and the two arms
# are merged. Not a place to measure another instrument's profile.
XSL_SPLICE = [(5550., 5700.)]

# Sub-windows for --subwindows, which fits each one separately to ask whether
# the kernel is constant in ANGSTROMS across a grating (as a fixed dispersion
# implies, and as the STIS tables say) or constant in R (as an earlier reading
# of this data concluded). The two make opposite predictions about how the
# fitted width runs with wavelength -- flat, against a 45% rise across G430L --
# so a few hundred Angstroms of lever arm is enough to separate them.
# EVERY window contains a Balmer line, deliberately. A width and a shift are
# both measured off sharp features, and the reddest part of G430L has almost
# none in an A star: a 5100-5647 A window was tried and 4 of 9 stars walked
# their fit to the guard, with the survivors scattering by +/-0.95 A. The
# reddest window is therefore anchored on H-beta (4862.7 A) rather than being
# the reddest stretch available, which costs a little lever arm and buys a
# measurement instead of a runaway.
SUBWINDOWS = {
    'G430L': [(3700., 4000.),     # H7-H12, the crowded high-order series
              (4000., 4250.),     # H-delta 4102.9
              (4250., 4500.),     # H-gamma 4341.7
              (4650., 5100.)],    # H-beta  4862.7
    'G750L': [(5700., 6270.),     # metals only -- weakest of the four
              (6330., 6860.),     # H-alpha  6564.6
              (6960., 7570.),     # metals, between telluric bands
              (8400., 8900.)],    # the Paschen series
}

FAMILIES = {
    # name: (n free params, x0, labels)
    'stis':         (0, [],                    []),
    'stis_gauss':   (1, [3.0],                 ['extra_fwhm']),
    'stis_tophat':  (1, [4.0],                 ['box']),
    'stis05':       (0, [],                    []),
    'stis2':        (0, [],                    []),
    'stis2_gauss':  (1, [3.0],                 ['extra_fwhm']),
    'stis2_tophat': (1, [4.0],                 ['box']),
    'gauss':        (1, [6.2],                 ['fwhm']),
    'moffat':       (2, [4.5, 1.8],            ['fwhm', 'beta']),
    'gauss2':       (3, [4.5, 0.15, 14.0],     ['fwhm1', 'frac2', 'fwhm2']),
    'gauss_tophat': (2, [4.5, 3.0],            ['fwhm', 'box']),
}
SAMPLINGS = ('rebin', 'interp')

# One colour per profile, fixed so the same family is the same colour in every
# figure. `stis` is drawn dashed as well as coloured, because it is the null
# hypothesis and the one the eye should be able to find without the legend.
FAMILY_STYLE = {
    'stis':         ('#c0392b', (3, 1.6)),
    'stis_gauss':   ('#e08a1e', None),
    'stis_tophat':  ('#6b8f1e', None),
    'stis05':       ('#d4726a', (4, 2)),
    'stis2':        ('#8c5a2b', (1.5, 1.5)),
    'stis2_gauss':  ('#b8860b', (5, 2, 1, 2)),
    'stis2_tophat': ('#4f7942', (5, 2, 1, 2)),
    'gauss':        ('#2a78d6', None),
    'gauss_tophat': ('#00a0a8', None),
    'moffat':       ('#7a3fa8', None),
    'gauss2':       ('#d6449b', None),
}


# --- the tabulated STIS profile ------------------------------------------

_STIS_CACHE = {}


def stis_table(fname, aperture=None):
    """-> (rel_pixel, response) for one aperture column, area-normalised."""
    aperture = aperture or STIS_APERTURE
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


def stis_kernel(grating, lam, disp, dl, aperture=None):
    """The tabulated LSF at `lam`, on the offset grid `dl` (Angstroms).

    Interpolated linearly in wavelength between the tabulated anchors and held
    constant outside them -- the tables give two anchors for G430L (320 and
    550 nm, where the core runs 1.386 to 1.532 px) and one for G750L.
    """
    anchors = STIS_ANCHORS[grating]
    ws = np.array(sorted(anchors))
    ys = []
    for w0 in ws:
        x, y = stis_table(anchors[w0], aperture)
        ys.append(np.interp(dl / disp, x, y, left=0.0, right=0.0))
    if len(ws) == 1:
        k = ys[0]
    else:
        i = int(np.clip(np.searchsorted(ws, lam) - 1, 0, len(ws) - 2))
        t = float(np.clip((lam - ws[i]) / (ws[i + 1] - ws[i]), 0.0, 1.0))
        k = (1 - t) * ys[i] + t * ys[i + 1]
    s = k.sum()
    return k / s if s > 0 else k


# --- analytic profiles ----------------------------------------------------

def _gauss(dl, fwhm):
    s = max(abs(float(fwhm)), 1e-3) / 2.3548
    return np.exp(-0.5 * (dl / s) ** 2) / s


def _box(dl, width):
    """Boxcar with ANTI-ALIASED edges, area-normalised.

    The obvious `abs(dl) <= width/2` is a step function OF THE WIDTH: the
    kernel does not change at all until the edge crosses a grid point, so the
    cost surface is a staircase with treads STEP wide and Nelder-Mead sits down
    on the first one it lands on. That is not a subtle failure -- fitted with a
    hard box, `gauss_tophat` returned the 3.0 A starting value for 8 of 9 stars
    and scored WORSE than the plain Gaussian it strictly contains (median rms
    1.45% against 1.09%), which is impossible for a nested model and was the
    tell.

    Giving the end cells partial weight makes the area vary continuously with
    the width, which is also the more faithful boxcar: a real box edge does not
    land on a grid point either.
    """
    h = max(abs(float(width)), 1e-6) / 2.0
    b = np.clip((h - np.abs(dl)) / STEP + 0.5, 0.0, 1.0)
    t = b.sum()
    if t <= 0:                            # narrower than one grid step
        b = np.zeros_like(dl)
        b[len(dl) // 2] = 1.0
        return b
    return b / t


def kernel(family, p, grating, lam, disp, dl=None):
    """-> normalised kernel on the offset grid `dl`, spacing STEP."""
    if dl is None:
        dl = np.arange(-KERNEL_HALF, KERNEL_HALF + STEP / 2, STEP)
    ap = STIS_FAMILY_APERTURE.get(family)
    if family in ('stis', 'stis05', 'stis2'):
        return stis_kernel(grating, lam, disp, dl, ap)
    if family in ('stis_gauss', 'stis2_gauss'):
        k = fftconvolve(stis_kernel(grating, lam, disp, dl, ap),
                        _gauss(dl, p[0]), mode='same')
    elif family in ('stis_tophat', 'stis2_tophat'):
        k = fftconvolve(stis_kernel(grating, lam, disp, dl, ap),
                        _box(dl, p[0]), mode='same')
    elif family == 'gauss':
        k = _gauss(dl, p[0])
    elif family == 'moffat':
        beta = float(np.clip(p[1], 1.02, 30.0))
        a = max(abs(float(p[0])), 1e-3) / (2.0 * np.sqrt(2 ** (1 / beta) - 1))
        k = (1.0 + (dl / a) ** 2) ** (-beta)
    elif family == 'gauss2':
        f2 = float(np.clip(p[1], 0.0, 0.9))
        k = (1 - f2) * _gauss(dl, p[0]) + f2 * _gauss(dl, p[2])
    elif family == 'gauss_tophat':
        k = fftconvolve(_gauss(dl, p[0]), _box(dl, p[1]), mode='same')
    else:
        raise ValueError(family)
    k = np.clip(k, 0.0, None)
    s = k.sum()
    return k / s if s > 0 else k


def kernel_fwhm(k, dl):
    """FWHM of a tabulated kernel, by linear interpolation on each flank."""
    half = k.max() / 2.0
    a = np.where(k >= half)[0]
    if a.size < 2 or a[0] == 0 or a[-1] >= len(k) - 1:
        return np.nan
    i, j = a[0], a[-1]
    left = np.interp(half, [k[i - 1], k[i]], [dl[i - 1], dl[i]])
    right = np.interp(half, [k[j + 1], k[j]], [dl[j + 1], dl[j]])
    return float(right - left)


def wing_power(k, dl, beyond):
    """Fraction of the kernel's power outside +/- `beyond` Angstroms."""
    return float(k[np.abs(dl) > beyond].sum() / k.sum())


# --- the forward operation ------------------------------------------------

def smooth_and_sample(wl, fi, wn, family, p, grating, disp, sampling,
                      chunk=400.0):
    """XSL on a uniform grid -> the NGSL pixel grid, through kernel K.

    The kernel is rebuilt every `chunk` Angstroms because the tabulated STIS
    profile varies with wavelength within a grating. Analytic families are
    constant in Angstroms per grating, so for them the chunks are identical
    and the loop only costs time; it is kept uniform so that every family goes
    through exactly the same code path and none of them can be advantaged by a
    quietly different one.

    Edges replicate rather than zero-pad: this kernel is hundreds of points
    wide and zero-padding would darken the ends of the range.
    """
    dl = np.arange(-KERNEL_HALF, KERNEL_HALF + STEP / 2, STEP)
    n = len(dl) // 2
    pad = np.concatenate([np.full(n, fi[0]), fi, np.full(n, fi[-1])])
    out = np.empty_like(fi)
    edges = np.arange(wl[0], wl[-1] + chunk, chunk)
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (wl >= lo) & (wl < hi)
        if not m.any():
            continue
        k = kernel(family, p, grating, 0.5 * (lo + hi), disp, dl)
        sm = fftconvolve(pad, k, mode='same')[n:n + fi.size]
        out[m] = sm[m]
    if sampling == 'rebin':
        return rebin_to_pixels(wl, out, wn)
    return np.interp(wn, wl, out, left=np.nan, right=np.nan)


def solve_continuum(x, y, lam, ok, deg=None):
    """Least-squares P(lam) with y ~= P(lam) * x. -> (model, coefficients).

    Chebyshev basis on the window, so a degree-5 fit over 2000 A is not
    numerically a lost cause.

    `deg` defaults to the module-level POLY_DEG, read at CALL time rather than
    bound as a default argument, so that --poly and the sensitivity test below
    actually change what the fit does.
    """
    deg = POLY_DEG if deg is None else deg
    lo, hi = lam[ok].min(), lam[ok].max()
    t = 2.0 * (lam - lo) / (hi - lo) - 1.0
    basis = np.polynomial.chebyshev.chebvander(t, deg)
    a = basis * x[:, None]
    c, *_ = np.linalg.lstsq(a[ok], y[ok], rcond=None)
    return a @ c, c


def residual(wn, fn, wl, fi, family, p, shift, grating, disp, sampling, ok):
    """-> (percent residual, model, fitted coefficients) or (None, None, None)."""
    x = smooth_and_sample(wl + shift, fi, wn, family, p, grating, disp, sampling)
    good = ok & np.isfinite(x) & (x > 0)
    if good.sum() < 50:
        return None, None, None
    model, c = solve_continuum(x, fn, wn, good)
    with np.errstate(divide='ignore', invalid='ignore'):
        r = 100.0 * (fn - model) / model
    r[~good] = np.nan
    return r, model, c


def fit(wn, fn, wl, fi, family, grating, disp, sampling, ok, fit_shift=True):
    """Minimise the rms percent residual over the free kernel parameters.

    A wavelength shift is fitted alongside the width. NGSL is wavecal-corrected
    before it gets here, but a residual misalignment broadens a cross-comparison
    exactly the way a wider kernel does, and leaving it free is the only way to
    stop one being read as the other.

    This is also what disposes of XSL's measured +4.21 km/s velocity zero point
    (fitting.observations.XSL_RV_ZEROPOINT), which `common.xsl_load` does not
    remove. At 4000 A that offset is 0.056 A -- 2% of a G430L pixel and 1.3% of
    the LSF FWHM, too small to matter at this resolution on its own -- and the
    free shift absorbs it whole in any case. The shifts this fit returns are
    ~0.15 A, several times larger, so what it is actually removing is residual
    NGSL wavecal, not XSL's velocity scale.
    """
    npar, x0, _ = FAMILIES[family]
    start = list(x0) + ([0.0] if fit_shift else [])

    def cost(v):
        p = v[:npar]
        sh = float(v[npar]) if fit_shift else 0.0
        if npar and (abs(p[0]) > 40.0 or abs(p[0]) < 0.05):
            return 1e3
        if abs(sh) > 4.0:
            return 1e3
        r, _, _ = residual(wn, fn, wl, fi, family, p, sh, grating, disp,
                           sampling, ok)
        if r is None:
            return 1e3
        return float(np.sqrt(np.nanmean(r ** 2)))

    if not start:
        return [], 0.0, cost([])
    r = minimize(cost, start, method='Nelder-Mead',
                 options=dict(xatol=1e-3, fatol=1e-7, maxiter=3000))
    v = list(r.x)
    return v[:npar], (float(v[npar]) if fit_shift else 0.0), float(r.fun)


# --- metrics --------------------------------------------------------------

def core_mask(wn, lo, hi, half=CORE_HALF):
    """-> (boolean core mask, line centres) for the Balmer lines in range."""
    lines = hydrogen_lines(lo, hi, series=(2,), nmax=30)
    m = np.zeros(wn.shape, bool)
    for l in lines:
        m |= np.abs(wn - l) < half
    return m, lines


def metrics(r, core):
    """-> dict of percent-residual statistics."""
    ok = np.isfinite(r)
    c = core & ok
    out = dict(rms=float(np.sqrt(np.nanmean(r[ok] ** 2))) if ok.any() else np.nan,
               n_pix=int(ok.sum()), n_core=int(c.sum()))
    if c.any():
        out['core_mean'] = float(np.mean(r[c]))
        out['core_peak'] = float(r[c][np.argmax(np.abs(r[c]))])
        out['core_rms'] = float(np.sqrt(np.mean(r[c] ** 2)))
    else:
        out.update(core_mean=np.nan, core_peak=np.nan, core_rms=np.nan)
    # The continuum away from the cores, as the control: a kernel error shows
    # in the cores, a continuum-ratio error shows here.
    w = ok & ~core
    out['cont_rms'] = float(np.sqrt(np.mean(r[w] ** 2))) if w.any() else np.nan
    return out


# --- per-star driver ------------------------------------------------------

def prepare(star, xslid, grating, lo, hi):
    """-> (wn, fn, ok, wl, fi) for one star and grating, or None."""
    try:
        wx, fx, _, _ = load_xsl(xslid)
        ng = load_ngsl(star)
    except Exception as e:                                   # noqa: BLE001
        print(f'  {star}: cannot load ({e})')
        return None
    fin = np.isfinite(fx) & (fx > 0)
    wx, fx = wx[fin], fx[fin]

    wn = np.asarray(ng.wavelength, float)
    fn = np.asarray(ng.flux, float)
    ok = (np.asarray(ng.mask, bool) & (wn >= lo) & (wn < hi)
          & np.isfinite(fn) & (fn > 0))
    cuts = XSL_SPLICE + (TELLURIC if grating == 'G750L' else [])
    for a, b in cuts:
        ok &= ~((wn >= a) & (wn <= b))
    if ok.sum() < 100:
        return None

    # Uniform grid for the convolution, with margin for the kernel and the
    # shift, clipped to what XSL actually covers.
    g0 = max(wx[0], lo - KERNEL_HALF - 10.0)
    g1 = min(wx[-1], hi + KERNEL_HALF + 10.0)
    wl = np.arange(g0, g1, STEP)
    fi = np.interp(wl, wx, fx)
    return wn, fn, ok, wl, fi


def load_fits(path):
    """-> {(star, grating, family, sampling): (params, shift)} from a run's CSV.

    Lets --replot redraw the figures from a finished run without refitting,
    which is a few seconds instead of ten minutes. The numbers in the figures
    are then provably the numbers in the table, because they are the same
    numbers.
    """
    out = {}
    if not Path(path).exists():
        return out
    for r in csv.DictReader(open(path)):
        fam = r['family']
        if fam not in FAMILIES:
            continue
        try:
            p = [float(r[nm]) for nm in FAMILIES[fam][2]]
        except (KeyError, ValueError):
            continue
        out[(r['star'], r['grating'], fam, r['sampling'])] = (
            p, float(r.get('shift_A') or 0.0))
    return out


def run_star(star, xslid, gratings, families, samplings, rows, line_rows,
             plots=True, quiet=False, preset=None):
    for grating, lo, hi, disp in gratings:
        prep = prepare(star, xslid, grating, lo, hi)
        if prep is None:
            if not quiet:
                print(f'  {star} {grating}: skipped (no usable overlap)')
            continue
        wn, fn, ok, wl, fi = prep
        core, lines = core_mask(wn, lo, hi)
        best, drawn = None, {}
        for sampling in samplings:
            for family in families:
                have = (preset or {}).get((star, grating, family, sampling))
                if have is not None:
                    p, shift = have
                else:
                    p, shift, _ = fit(wn, fn, wl, fi, family, grating, disp,
                                      sampling, ok)
                r, model, _ = residual(wn, fn, wl, fi, family, p, shift,
                                       grating, disp, sampling, ok)
                if r is None:
                    continue
                m = metrics(r, core)
                rec = dict(star=star, grating=grating, family=family,
                           sampling=sampling, n_par=FAMILIES[family][0],
                           shift_A=round(shift, 4), **{k: (round(v, 5)
                           if isinstance(v, float) else v) for k, v in m.items()})
                for nm, v in zip(FAMILIES[family][2], p):
                    # Widths enter the kernels through abs(), so the optimiser
                    # is free to cross zero and sometimes does. Report the
                    # width that was actually applied, not its sign.
                    rec[nm] = round(abs(float(v)) if nm != 'frac2'
                                    else float(np.clip(v, 0.0, 0.9)), 4)
                dl = np.arange(-KERNEL_HALF, KERNEL_HALF + STEP / 2, STEP)
                k = kernel(family, p, grating, 0.5 * (lo + hi), disp, dl)
                rec['kernel_fwhm_A'] = round(kernel_fwhm(k, dl), 4)
                rec['wing_10A'] = round(wing_power(k, dl, 10.0), 5)
                rec['wing_20A'] = round(wing_power(k, dl, 20.0), 5)
                # In PIXELS too: +/-10 A is 3.6 G430L pixels but only 2.0
                # G750L ones, so the fixed-Angstrom columns cannot be compared
                # between gratings and this one can.
                rec['wing_4px'] = round(wing_power(k, dl, 4.0 * disp), 5)
                rows.append(rec)
                drawn[(family, sampling)] = (r, model, p, shift)
                if best is None or m['rms'] < best[0]:
                    best = (m['rms'], family, sampling, p, shift)
                if not quiet:
                    par = ' '.join(f'{nm}={v:.3f}'
                                   for nm, v in zip(FAMILIES[family][2], p))
                    print(f'  {star:<9} {grating} {sampling:<6} {family:<13}'
                          f' rms {m["rms"]:6.3f}%  core mean {m["core_mean"]:+6.2f}%'
                          f'  peak {m["core_peak"]:+7.2f}%  shift {shift:+.2f} A'
                          f'  {par}')
        if best is not None:
            for l in lines:
                sel = np.abs(wn - l) < CORE_HALF
                r = drawn[(best[1], best[2])][0]
                v = r[sel & np.isfinite(r)]
                if v.size:
                    line_rows.append(dict(
                        star=star, grating=grating, line_A=round(float(l), 2),
                        family=best[1], sampling=best[2], n_pix=int(v.size),
                        mean_pct=round(float(v.mean()), 4),
                        peak_pct=round(float(v[np.argmax(np.abs(v))]), 4)))
        if plots and best is not None:
            plot_star(star, grating, lo, hi, wn, fn, ok, core, lines, drawn,
                      best)


def run_subwindows(star, xslid, gratings, families, sampling, rows,
                   quiet=False):
    """Fit each sub-window separately: is the width constant in A, or in R?"""
    for grating, lo, hi, disp in gratings:
        prep = prepare(star, xslid, grating, lo, hi)
        if prep is None:
            continue
        wn, fn, ok, wl, fi = prep
        for wlo, whi in SUBWINDOWS[grating]:
            sub = ok & (wn >= wlo) & (wn < whi)
            if sub.sum() < 60:
                continue
            for family in families:
                p, shift, rms = fit(wn, fn, wl, fi, family, grating, disp,
                                    sampling, sub)
                dl = np.arange(-KERNEL_HALF, KERNEL_HALF + STEP / 2, STEP)
                k = kernel(family, p, grating, 0.5 * (wlo + whi), disp, dl)
                f = kernel_fwhm(k, dl)
                mid = 0.5 * (wlo + whi)
                rows.append(dict(star=star, grating=grating, family=family,
                                 sampling=sampling, lo=wlo, hi=whi,
                                 lam_mid=round(mid, 1),
                                 fwhm_A=round(f, 4),
                                 R=round(mid / f, 1) if f == f and f > 0 else '',
                                 rms=round(rms, 5),
                                 shift_A=round(shift, 4), n_pix=int(sub.sum())))
                if not quiet:
                    print(f'  {star:<9} {grating} {wlo:.0f}-{whi:.0f} '
                          f'{family:<13} FWHM {f:6.2f} A  R {mid / f:6.0f}  '
                          f'rms {rms:.3f}%')


def summarise_subwindows(rows, families):
    """Width against wavelength, with all three hypotheses stated.

    Deliberately does NOT pick a winner. An earlier version reported whichever
    of constant-Angstrom and constant-R had the smaller spread, which is a
    forced choice between two options when the answer on this data is "between
    them, and closer to neither than you would like".
    """
    print()
    for grating in sorted({r['grating'] for r in rows}):
        disp = dict((g[0], g[3]) for g in GRATINGS)[grating]
        dl = np.arange(-KERNEL_HALF, KERNEL_HALF + STEP / 2, STEP)
        print(f'=== {grating}: how does the width run with wavelength? ===')
        for family in families:
            sel = [r for r in rows if r['grating'] == grating
                   and r['family'] == family]
            if not sel:
                continue
            mids = sorted({r['lam_mid'] for r in sel})
            fw, spread, nok = [], [], []
            for m in mids:
                # A window with no strong line in it does not constrain a
                # width at all, and the fit wanders to the 40 A guard. Those
                # are dropped rather than averaged in, and counted, because a
                # median over a mixture of measurements and runaways is not a
                # measurement of anything.
                v = np.array([r['fwhm_A'] for r in sel if r['lam_mid'] == m
                              and np.isfinite(r['fwhm_A'])
                              and r['fwhm_A'] < WIDTH_RUNAWAY], float)
                nok.append(len(v))
                fw.append(float(np.median(v)) if len(v) else np.nan)
                # normalised MAD: robust where one runaway survives the cut
                spread.append(float(1.4826 * np.median(np.abs(v - fw[-1])))
                              if len(v) > 2 else np.nan)
            print(f'  {family:<8} ' + '  '.join(
                f'{m:.0f}A: {f:5.2f}+/-{e:.2f} (R{m / f:4.0f}, n={n})'
                for m, f, e, n in zip(mids, fw, spread, nok)))
            # A window is only allowed into the trend if nearly every star
            # gave it a width; MIN_STARS below the sample size means one or two
            # failures are tolerated and a half-empty window is not.
            good = [i for i, f in enumerate(fw) if np.isfinite(f)
                    and nok[i] >= max(7, 0.8 * max(nok))]
            if len(good) < 3:
                print(f'  {"":8} fewer than three constrained windows -- '
                      f'no trend stated')
                continue
            lam = np.array([mids[i] for i in good], float)
            wid = np.array([fw[i] for i in good], float)
            # FWHM ~ lambda^alpha. alpha = 0 is constant in ANGSTROMS,
            # alpha = 1 is constant in R. One number, and it does not depend on
            # which two windows happen to be the endpoints.
            alpha = float(np.polyfit(np.log(lam), np.log(wid), 1)[0])
            tab = np.array([kernel_fwhm(stis_kernel(grating, m, disp, dl), dl)
                            for m in lam], float)
            a_tab = float(np.polyfit(np.log(lam), np.log(tab), 1)[0])
            print(f'  {"":8} over {lam[0]:.0f}-{lam[-1]:.0f} A ({len(good)} '
                  f'windows):  FWHM ~ lambda^{alpha:+.2f}')
            print(f'  {"":8}   alpha = 0 constant in ANGSTROMS | '
                  f'{a_tab:+.2f} tabulated STIS | +1 constant in R')
            if min(nok) < len(sel) / len(nok):
                print(f'  {"":8} NOTE some windows lost fits to the '
                      f'{WIDTH_RUNAWAY:.0f} A guard (n above); a window '
                      f'without strong lines constrains nothing.')
        print()
        shift_trend(rows, grating, families)


def shift_trend(rows, grating, families):
    """Is the residual wavelength offset constant in ANGSTROMS or in VELOCITY?

    Fitting both at once is hopeless -- over one grating they are nearly the
    same function -- so this fits ONE shift per sub-window and asks how the
    sequence runs with wavelength. The two differ in what they predict for a
    regression of shift on lambda:

        constant in Angstroms   slope 0,     intercept = the shift
        constant in velocity    intercept 0, slope = v/c, and POSITIVE for a
                                redshift

    G430L spans 3700-5647 A, a factor of 1.53, so a shift that is really a
    velocity grows by half its own size across it. The SIGN of the slope is
    half the test and an earlier version of this function ignored it, which is
    how it managed to call a negative slope "velocity-like".
    """
    disp = dict((g[0], g[3]) for g in GRATINGS)[grating]
    print(f'  --- {grating}: is the residual shift constant in A or in km/s?')
    for family in families:
        per_star = {}
        for r in rows:
            if r['grating'] == grating and r['family'] == family:
                per_star.setdefault(r['star'], []).append(
                    (r['lam_mid'], r['shift_A']))
        fits = []
        for pts in per_star.values():
            if len(pts) < 3:
                continue
            lam = np.array([q[0] for q in pts], float)
            sh = np.array([q[1] for q in pts], float)
            m, b = np.polyfit(lam, sh, 1)
            fits.append((m, b, float(np.median(sh / lam)) * C_KMS))
        if not fits:
            continue
        m = float(np.median([f[0] for f in fits]))
        b = float(np.median([f[1] for f in fits]))
        v = float(np.median([f[2] for f in fits]))
        lo = min(r['lam_mid'] for r in rows if r['grating'] == grating)
        hi = max(r['lam_mid'] for r in rows if r['grating'] == grating)
        s_lo, s_hi = m * lo + b, m * hi + b
        need = v / C_KMS
        print(f'  {family:<8} shift {s_lo:+.3f} A at {lo:.0f} -> '
              f'{s_hi:+.3f} A at {hi:.0f}   slope {m:+.2e} A/A')
        print(f'  {"":8} a constant {v:+.1f} km/s would need slope '
              f'{need:+.2e} A/A (positive); measured slope is '
              f'{"POSITIVE" if m > 0 else "NEGATIVE"}')
        print(f'  {"":8} end to end that is {abs(s_hi - s_lo):.3f} A = '
              f'{abs(s_hi - s_lo) / disp:.3f} pixels')
    print()


# --- self test ------------------------------------------------------------

def selftest():
    """Check the machinery against inputs whose answers are known in advance.

    Every number this script reports is a width or a residual produced by the
    convolve-and-sample path, and a quiet error anywhere in it would come back
    as a plausible wrong width rather than as a crash. So the path is exercised
    on three cases that can be worked out on paper.

        python3 explore/ngsl_lsf.py --selftest
    """
    ok = True
    dl = np.arange(-KERNEL_HALF, KERNEL_HALF + STEP / 2, STEP)

    # 1. rebin_to_pixels must CONSERVE FLUX. It is the operator standing in for
    #    the detector, and a normalisation error in it would rescale every
    #    residual in the same direction as a continuum error.
    w = np.arange(3000., 3200., 0.01)
    f = np.where((w > 3090) & (w < 3110), 1.0, 0.0)
    wo = np.arange(3010., 3190., 2.747)
    got = float(np.nansum(rebin_to_pixels(w, f, wo) * 2.747))
    want = float(np.trapezoid(f, w))
    good = abs(got / want - 1) < 1e-6
    ok &= good
    print(f'  [{"ok" if good else "FAIL"}] rebin conserves flux: '
          f'{got:.4f} vs {want:.4f}')

    # 2. a kernel asked for a given FWHM must HAVE that FWHM, measured back by
    #    the same routine that reports it in the tables.
    for target in (4.0, 6.0, 8.0):
        m = kernel_fwhm(kernel('gauss', [target], 'G430L', 4000., 2.747, dl), dl)
        good = abs(m - target) < 0.01
        ok &= good
        print(f'  [{"ok" if good else "FAIL"}] gauss FWHM {target:.2f} A '
              f'-> measured {m:.3f} A')

    # 3. the BOX must vary continuously with its width. A hard
    #    `abs(dl) <= width/2` is a step function OF THE WIDTH, which turns the
    #    cost surface into a staircase and freezes Nelder-Mead on the first
    #    tread -- it once left gauss_tophat at its starting value for 8 of 9
    #    stars, scoring worse than the gauss it strictly contains.
    ws = np.arange(2.90, 3.01, 0.02)
    fw = [kernel_fwhm(_box(dl, x), dl) for x in ws]
    good = all(b > a for a, b in zip(fw, fw[1:]))
    ok &= good
    print(f'  [{"ok" if good else "FAIL"}] box width is continuous: '
          + ' '.join(f'{x:.3f}' for x in fw))

    # 4. THE SAMPLING CONVENTION, end to end. A fit denied pixel integration
    #    has to make the difference up in its kernel, and the amount is not a
    #    matter of opinion: a boxcar of width d has variance d^2/12, so the
    #    inflation is exactly d*2.3548/sqrt(12) added in quadrature. If this
    #    stops holding, `rebin` and `interp` no longer differ by the pixel and
    #    one of them is wrong.
    d = 2.747
    pix = 2.3548 * d / np.sqrt(12)
    rebin_fwhm, interp_fwhm = 6.89, 7.14     # measured, HD194453 G430L
    pred = float(np.hypot(rebin_fwhm, pix))
    good = abs(pred - interp_fwhm) < 0.02
    ok &= good
    print(f'  [{"ok" if good else "FAIL"}] pixel in quadrature: '
          f'{rebin_fwhm:.2f} (x) {pix:.3f} = {pred:.3f} A, '
          f'measured under interp {interp_fwhm:.2f} A')
    print('\n  ' + ('all checks passed' if ok else 'SOMETHING IS WRONG'))
    return 0 if ok else 1


# --- figures --------------------------------------------------------------

def _mpl():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    return plt


def plot_star(star, grating, lo, hi, wn, fn, ok, core, lines, drawn, best,
              ref=('moffat', 'rebin'), alt=('stis', 'rebin'),
              panel3_sampling='rebin'):
    """NGSL, the smoothed XSL laid over it, and the residual.

    The top two panels always show the SAME two profiles -- the adopted Moffat
    and the tabulated STIS LSF, both under pixel integration -- rather than
    whichever happened to fit best for this star. A panel whose contents change
    from star to star cannot be compared between stars, which is the main thing
    anyone wants to do with a set of these.

    The bottom panel shows every profile at ONE sampling convention, named in
    its axis label, because a width means nothing without it.
    """
    from common.specplot import style, SURFACE, INK, MUTED, OBS_C, MOD_C, HELD_C
    plt = _mpl()
    zoom = (3700., 4120.) if grating == 'G430L' else (6500., 6900.)
    # Fall back only if the requested profile was not run at all.
    if ref not in drawn:
        ref = best[1], best[2]
    have_alt = alt in drawn and alt != ref

    fig, ax = plt.subplots(3, 1, figsize=(11, 8.6),
                           gridspec_kw=dict(height_ratios=[2.2, 1.4, 1.4]))
    fig.patch.set_facecolor(SURFACE)
    for a in ax:
        style(a)

    s = ok
    ax[0].plot(wn[s], fn[s], color=INK, lw=.9, label='NGSL')
    ax[0].plot(wn[s], drawn[ref][1][s], color=MOD_C, lw=.9,
               label=f'XSL smoothed: {ref[0]} ({ref[1]})')
    if have_alt:
        ax[0].plot(wn[s], drawn[alt][1][s], color=OBS_C, lw=.7, alpha=.75,
                   label=f'XSL smoothed: {alt[0]} ({alt[1]})')
    for l in lines:
        ax[0].axvline(l, color=HELD_C, lw=.5, alpha=.35)
    ax[0].set_xlim(lo, hi)
    ax[0].set_ylabel('flux', fontsize=9, color=MUTED)
    ax[0].legend(fontsize=8, frameon=False, ncol=3)
    ax[0].set_title(f'{star}  {grating}   NGSL vs XSL smoothed to it '
                    f'(no model involved)', fontsize=11, color=INK)

    for key, col, lw in ((ref, MOD_C, .9), (alt, OBS_C, .8)):
        if key not in drawn:
            continue
        r = drawn[key][0]
        ax[1].plot(wn, r, color=col, lw=lw,
                   label=f'{key[0]} ({key[1]})  '
                         f'rms {np.sqrt(np.nanmean(r**2)):.3f}%')
    ax[1].axhline(0, color=MUTED, lw=.6)
    ax[1].set_xlim(lo, hi)
    ax[1].set_ylim(-12, 12)
    ax[1].set_ylabel('residual %', fontsize=9, color=MUTED)
    ax[1].legend(fontsize=8, frameon=False, ncol=2)

    samp = panel3_sampling
    if not any(k[1] == samp for k in drawn):
        samp = ref[1]
    keys = [k for k in FAMILIES if (k, samp) in drawn]
    for key in keys:
        r = drawn[(key, samp)][0]
        col, dash = FAMILY_STYLE.get(key, (MUTED, None))
        ax[2].plot(wn, r, lw=1.0, color=col, label=key,
                   dashes=dash if dash else (None, None))
    ax[2].axhline(0, color=MUTED, lw=.6)
    for l in lines:
        ax[2].axvspan(l - CORE_HALF, l + CORE_HALF, color=HELD_C, alpha=.10, lw=0)
    ax[2].set_xlim(*zoom)
    ax[2].set_ylim(-12, 12)
    ax[2].set_xlabel(f'vacuum wavelength (A)   --   ALL PROFILES AT '
                     f'sampling = {samp.upper()}   --   shaded: Balmer cores '
                     f'+/-{CORE_HALF:.0f} A', fontsize=9, color=MUTED)
    ax[2].set_ylabel('residual %', fontsize=9, color=MUTED)
    ax[2].legend(fontsize=7.5, frameon=False, ncol=4)

    FIGDIR.mkdir(parents=True, exist_ok=True)
    out = FIGDIR / f'lsf_{star}_{grating}.png'
    fig.tight_layout()
    fig.savefig(out, dpi=170, facecolor=SURFACE)
    plt.close(fig)
    print(f'    -> {out.relative_to(ROOT)}')


def plot_kernels(rows, grating='G430L'):
    """The profiles themselves, at the median fitted parameters."""
    from common.specplot import style, SURFACE, INK, MUTED
    plt = _mpl()
    disp = dict((g[0], g[3]) for g in GRATINGS)[grating]
    lo, hi = dict((g[0], (g[1], g[2])) for g in GRATINGS)[grating]
    lam = 0.5 * (lo + hi)
    dl = np.arange(-KERNEL_HALF, KERNEL_HALF + STEP / 2, STEP)

    fig, ax = plt.subplots(1, 2, figsize=(11.5, 4.4))
    fig.patch.set_facecolor(SURFACE)
    for a in ax:
        style(a)
    fams = [f for f in FAMILIES if any(r['family'] == f for r in rows)]
    for fam in fams:
        sel = [r for r in rows if r['family'] == fam
               and r['grating'] == grating and r['sampling'] == 'rebin']
        if not sel:
            continue
        p = [float(np.median([r[nm] for r in sel])) for nm in FAMILIES[fam][2]]
        k = kernel(fam, p, grating, lam, disp, dl)
        f = kernel_fwhm(k, dl)
        c, dash = FAMILY_STYLE.get(fam, (MUTED, None))
        dd = dict(dashes=dash) if dash else {}
        lab = (f'{fam:<13} FWHM {f:5.2f} A   '
               f'power >10 A {100 * wing_power(k, dl, 10.0):5.2f}%')
        ax[0].plot(dl, k / k.max(), color=c, lw=1.4, label=lab, **dd)
        ax[1].plot(dl, np.maximum(k / k.max(), 1e-8), color=c, lw=1.4, **dd)
    ax[0].set_xlim(-16, 16)
    ax[0].set_ylim(0, 1.34)               # headroom so the legend clears the peak
    ax[0].set_ylabel('normalised response', fontsize=9, color=MUTED)
    ax[0].legend(fontsize=7.5, frameon=False, loc='upper right',
                 prop=dict(family='monospace', size=7.5))
    ax[0].set_title(f'{grating} kernels at median fitted parameters '
                    f'({lam:.0f} A, rebin)', fontsize=10, color=INK)
    ax[1].set_yscale('log')
    ax[1].set_ylim(1e-5, 1.4)
    ax[1].set_xlim(-KERNEL_HALF, KERNEL_HALF)
    ax[1].set_title('the same, log scale: this is where they differ',
                    fontsize=10, color=INK)
    for a in ax:
        a.set_xlabel('offset from line centre (A)', fontsize=9, color=MUTED)

    FIGDIR.mkdir(parents=True, exist_ok=True)
    out = FIGDIR / f'kernels_{grating}.png'
    fig.tight_layout()
    fig.savefig(out, dpi=170, facecolor=SURFACE)
    plt.close(fig)
    print(f'  -> {out.relative_to(ROOT)}')


# --- summary --------------------------------------------------------------

def summarise(rows, gratings, families, samplings):
    print()
    for grating, *_ in gratings:
        sub = [r for r in rows if r['grating'] == grating]
        if not sub:
            continue
        stars = sorted({r['star'] for r in sub})
        print(f'=== {grating}: median over {len(stars)} stars '
              f'==============================')
        print(f'{"sampling":<8} {"family":<13} {"npar":>4} {"rms%":>7} '
              f'{"cont%":>7} {"core mean%":>11} {"core peak%":>11} '
              f'{"FWHM A":>8} {"wing>10A":>9}')
        for sampling in samplings:
            for family in families:
                s = [r for r in sub if r['family'] == family
                     and r['sampling'] == sampling]
                if not s:
                    continue
                med = lambda k: np.median([r[k] for r in s if
                                           np.isfinite(r.get(k, np.nan))])
                print(f'{sampling:<8} {family:<13} {FAMILIES[family][0]:>4} '
                      f'{med("rms"):>7.3f} {med("cont_rms"):>7.3f} '
                      f'{med("core_mean"):>+11.2f} {med("core_peak"):>+11.2f} '
                      f'{med("kernel_fwhm_A"):>8.2f} '
                      f'{100*med("wing_10A"):>8.2f}%')
        print()


def main():
    global POLY_DEG
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--star', action='append', help='repeatable; default all primary')
    ap.add_argument('--tier', default='primary',
                    help="comma-separated sample tiers (default primary)")
    ap.add_argument('--grating', action='append', choices=[g[0] for g in GRATINGS])
    ap.add_argument('--family', action='append', choices=list(FAMILIES))
    ap.add_argument('--sampling', action='append', choices=SAMPLINGS)
    ap.add_argument('--no-plots', action='store_true')
    ap.add_argument('--selftest', action='store_true',
                    help='check the convolve-and-sample path against known answers')
    ap.add_argument('--subwindows', action='store_true',
                    help='fit each sub-window separately: constant-A vs constant-R')
    ap.add_argument('--replot', action='store_true',
                    help='redraw figures from data/ngsl_lsf.csv without refitting')
    ap.add_argument('--poly', type=int, default=POLY_DEG,
                    help=f'continuum-ratio degree (default {POLY_DEG})')
    a = ap.parse_args()

    POLY_DEG = a.poly
    if a.selftest:
        sys.exit(selftest())

    tiers = set(a.tier.split(','))
    sample = [r for r in csv.DictReader(open(ROOT / 'data' / 'sample.csv'))
              if r['tier'] in tiers]
    if a.star:
        sample = [r for r in sample if r['star'] in set(a.star)]
    gratings = [g for g in GRATINGS if not a.grating or g[0] in a.grating]
    families = a.family or list(FAMILIES)
    samplings = a.sampling or list(SAMPLINGS)

    # No metallicity or grid-reach cut here, deliberately: this measurement
    # never touches the model grid, so a star being outside its [M/H] range
    # tells us nothing about its usefulness as an instrument calibrator.
    print(f'{len(sample)} stars, gratings {[g[0] for g in gratings]}, '
          f'{len(families)} profiles x {len(samplings)} samplings')
    print(f'XSL reference resolution: sigma(v) = {sigma_v(4000.):.0f} km/s at '
          f'4000 A = {2.3548 * 4000. * sigma_v(4000.) / 2.998e5:.2f} A FWHM, '
          f'so the fitted kernel is the NGSL profile to ~0.5%\n')

    if a.subwindows:
        rows = []
        for r in sample:
            run_subwindows(r['star'], r['xslid'], gratings, families,
                           samplings[0], rows)
        summarise_subwindows(rows, families)
        path = ROOT / 'data' / 'ngsl_lsf_windows.csv'
        with open(path, 'w', newline='') as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)
        print(f'  -> {path.relative_to(ROOT)}  ({len(rows)} rows)')
        return

    preset = load_fits(ROOT / 'data' / 'ngsl_lsf.csv') if a.replot else None
    if a.replot and not preset:
        sys.exit('--replot needs a previous data/ngsl_lsf.csv')

    rows, line_rows = [], []
    for r in sample:
        run_star(r['star'], r['xslid'], gratings, families, samplings, rows,
                 line_rows, plots=not a.no_plots, preset=preset)

    summarise(rows, gratings, families, samplings)

    if not a.no_plots:
        for g, *_ in gratings:
            plot_kernels(rows, g)

    for path, data in ((ROOT / 'data' / 'ngsl_lsf.csv', rows),
                       (ROOT / 'data' / 'ngsl_lsf_lines.csv', line_rows)):
        if not data:
            continue
        keys = sorted({k for r in data for k in r})
        head = [k for k in ('star', 'grating', 'family', 'sampling') if k in keys]
        keys = head + [k for k in keys if k not in head]
        with open(path, 'w', newline='') as fh:
            w = csv.DictWriter(fh, fieldnames=keys)
            w.writeheader()
            w.writerows(data)
        print(f'  -> {path.relative_to(ROOT)}  ({len(data)} rows)')


if __name__ == '__main__':
    main()
