"""UVES-POP smoothed to NGSL pixels, for the 13 stars in both libraries.

The 13 come from data/uves_ngsl_overlap.csv (explore/uves_ngsl_overlap.py).
This is a data-to-data comparison: same star, same epoch-independent quantity,
two instruments -- HST/STIS at a 3.54 A Moffat core and VLT/UVES delivered on a
0.1 A grid. Smoothing the high-resolution spectrum to the low-resolution one is
what makes them comparable, and it goes in ONE direction only: UVES -> NGSL,
through `common.lsf.to_ngsl_pixels`. There is no way to go the other way.

The UVES-POP spectra are fetched on first run (~2.9 MB each, 13 stars) into
data/uves_pop/, which .gitignore keeps out of the repository.

Writes figures/explore_libraries/uves_ngsl_break.png     (Balmer break, 3646 A)
       figures/explore_libraries/uves_ngsl_hepsilon.png  (H-epsilon, 3971 A)
       data/uves_ngsl_compare.csv                        (the fitted numbers)

SIGN CONVENTION: the residual is (NGSL - UVES@NGSL) / UVES@NGSL, in percent.
POSITIVE MEANS NGSL IS BRIGHTER than UVES-POP smoothed to NGSL's resolution.
Neither is a model, so neither is "truth" -- this is a disagreement between two
flux calibrations, not an error in one of them.

THREE THINGS HAVE TO BE RIGHT BEFORE THE RESIDUAL MEANS ANYTHING

1. Frames. UVES-POP's delivered spectra are in the OBSERVED frame, which this
   project had not established before. Verified on the two fastest stars: the
   Ca II H and K minima sit at rest*(1+v/c) for the catalogued v, recovering
   +123.0 and +122.6 km/s for HD076932 (catalog +119.76) and +107.7 and +107.5
   for HD063077 (catalog +106.93). So the catalog RV is removed here.

2. NGSL's wavelength zero point. data/ngsl_wavecal.csv covers the 14 sample
   stars and NONE of these 13, so `apply_wavecal` would give air->vacuum only
   and leave the ~0.5-1.2 A G430L residual it was written to remove. Instead
   the shift is FITTED here, per star, against the smoothed UVES spectrum --
   which absorbs NGSL's own RV and its wavecal error together, since nothing in
   this comparison can separate them. It is reported, not hidden.

3. The grey factor. The two libraries differ by a flat scale factor, as NGSL
   and XSL do (docs/DATA.md): ground-based UVES must correct slit losses and
   NGSL need not. UVES is scaled to the NGSL flux scale over the normalization
   window so the panels compare SHAPE, and the factor is printed per star.

THE ALIGNMENT IS GOOD ENOUGH FOR THE H-EPSILON PANELS, and this was checked
rather than assumed. Refitting the shift locally over 3800-4300 A, on top of
the global one, moves it by -9.6 to +12.0 km/s across the sample -- 0.16 A at
H-epsilon at worst, about 6% of an NGSL pixel and 5% of the 3.54 A core. So the
structure the residual panels show around H-epsilon, which reaches 10% on some
stars, is NOT a wavelength artefact: an offset that small cannot put that much
flux error into a feature that wide. What it IS this figure does not say.

WHAT IS NOT CORRECTED: UVES-POP's own resolution is not deconvolved before the
NGSL kernel is applied. The delivered 0.1 A grid resolves ~0.2 A FWHM against
NGSL's 3.54 A core, so the quadrature error is 0.16% of the core width -- an
order below the +/-0.16 A star-to-star scatter of the core itself. Same
argument as `xsl_to_ngsl` in explore/plot_ngsl_xsl.py, where it is 0.4%.

SEVEN OF THE THIRTEEN ARE CATALOGUED VARIABLES, and the two libraries observed
them years apart -- so for those stars there is a floor on how well any
comparison can agree, and it is not an instrumental floor. The GCVS type comes
out of the UVES-POP file itself and is printed on every panel. It is not a
small effect here: the three largest Balmer-break residuals in the sample are
FW CMa (BE, a Be star), nu Vir (SRB) and S Mon (IA:), which are also the three
most strongly variable types present. With n = 13 that ordering is an
observation, not a trend worth fitting -- but a panel showing a Be star should
not be read as a statement about flux calibration.

GAIA XP IS AVAILABLE BUT OFF BY DEFAULT (--xp turns it on, 5 of the 13 stars:
HD058343, HD063077, HD076932, HD111786, HD115617). It was tried and it did not
earn its place on this figure.

Why it does not help. XP is R ~ 52 at the break against NGSL's 3.54 A core, so
the break is carried by about a dozen 2 nm samples; its range starts at 3360 A,
only 286 A blueward of the break, where its calibration is least tested; and
XP/NGSL carries a known INSTRUMENTAL pattern of up to 4.8%, V-shaped about the
BP/RP join, already measured in explore/xp_vs_ngsl.py (docs/DATA.md). Measured
here over 3400-4000 A, three of the five stars give the SAME -1.5% blue-to-red
tilt across the break -- that is the instrumental shape, not the stars. The
other two are the Be star (-11.6%) and HD115617 (+7.2%, unexplained).

What the code still provides, and why it is kept: the retrieval by source_id
(NOT by cone search -- these stars move up to 27" between J2000 and DR3's epoch
2016, so a 5" cone misses five of them outright), and the effective XP width
measured against NGSL. See XP_FWHM_A.

Which 5 stars have XP is not a brightness cut in any simple sense: HD115617 at
G = 4.53 has it while HD047839 at G = 4.54 and HD099648 at G = 4.68 do not, and
HD142703 at G = 6.05 does not while HD111786 at G = 6.09 does. HD206778 is not
in Gaia DR3 at all -- at V = 2.40 it is beyond the bright limit, confirmed by
finding no DR3 source within 10" of its epoch-2016 position.

H-EPSILON IS BLENDED WITH INTERSTELLAR Ca II H. They are 1.64 A apart (3971.24
and 3969.60 vacuum) and NGSL's core is 3.54 A, so in the NGSL panels they are
ONE feature and cannot be separated. The UVES spectrum resolves them, which is
the whole reason the zoom is worth drawing -- but no H-epsilon depth read off
an NGSL pixel here is a hydrogen measurement.
"""
import csv
import sys
import argparse
import warnings
import urllib.request
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.io.fits.verify import VerifyWarning
from scipy.ndimage import median_filter

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common.lsf import (broaden, broaden_ngsl, rebin_to_pixels,
                        to_ngsl_pixels, NGSL_TRUNC_PX)
from common.lines import air_to_vac, balmer_member, ISM_LINES
from common.uves_pop_load import load as load_uves
from common.figpath import library_figure_path
from explore.plot_ngsl_vs_model import (BALMER, OBS_C, SURFACE, INK, MUTED,
                                        GRID)

ROOT = Path(__file__).resolve().parent.parent
UVES_C = '#1baf7a'          # raw UVES-POP, as delivered
UVESD_C = '#4a3aa7'         # UVES degraded to NGSL: dashed, to read as derived
JOIN_C = '#b07a2a'          # UVES setting joins: instrumental, not stellar
XP_C = '#c23b6e'            # Gaia XP: a third instrument, not a smoothed version
XP_DIR = ROOT / 'data' / 'gaia_xp'

# XP's effective resolution over the break panel, MEASURED against NGSL on the
# 5 stars that have both: convolve NGSL with a Gaussian, integrate onto the XP
# 2 nm grid, fit one scale factor, and scan the width. Best fits come out 50,
# 68, 70, 98 and 50 A -- median 68, so 70 is adopted. R ~ 52 at 3646 A.
#
# THE CHOICE BARELY MATTERS, which is the useful part. The residual rms changes
# by less than 0.2 PERCENTAGE POINTS between each star's own best width and a
# common 70 A, while the residual itself runs 5-12%. So what the ratio panel
# shows is a flux difference, not a resolution mismatch, and it would look the
# same for any width in the measured range. Structure on the ~70 A scale is
# still not interpretable; the broad level and the size of the break are.
XP_FWHM_A = 70.0
C_KMS = 2.99792458e5

HEPS = balmer_member(7)                  # 3971.24 A vacuum
CAH = ISM_LINES['Ca II H']               # 3969.60 A vacuum, 1.64 A away

# Balmer-break panel. Wider than plot_ngsl_vs_model's ZOOM_B = (3600, 4000) so
# that both sides of the break are on screen: the point here is the size of the
# discontinuity, and a panel that starts at 3600 shows only its red side.
BREAK_XLIM = (3400.0, 4000.0)
# Normalised BLUEWARD of the break only, the project's convention (common/
# lines.py WINS_B): with windows on both sides a break-amplitude difference
# splits into a half-size offset on each side instead of showing up whole.
BREAK_WINS = [(3400., 3620.)]

# H-epsilon panel. +/-30 A clears H8 (3890.2) and H-delta (4102.9); the
# normalisation wings stop short of the Ca II H + H-epsilon blend and of
# interstellar Ca II K at 3934.8.
HEPS_XLIM = (HEPS - 30.0, HEPS + 30.0)
HEPS_WINS = [(3944., 3962.), (3982., 4000.)]

# UVES SETTING JOINS -- where two instrument settings are spliced together in
# the merged product. MEASURED HERE, not taken from the instrument manual,
# because the delivered file keeps only one exposure's header.
#
# Three independent handles agree:
#
#  1. HD138716 is MISSING one exposure entirely. Its coverage stops at 3859.2 A
#     and resumes at 4784.1 A, which brackets the missing setting directly --
#     no modelling, just where the star has no data and the other twelve do.
#  2. Bad-pixel dropouts cluster at the same wavelengths across stars: four have
#     one in 3733.6-3753.5 A, and HD111786 has one at 4781.9-4782.7 A, which is
#     HD138716's resume point to within 2 A.
#  3. The stacked error/sqrt(flux) of all 13 steps by a factor 1.66 -- x2.8 in
#     effective exposure -- between 3750 and 3850 A, and nowhere else in the
#     blue.
#
# And the same comparison says the REDDER setting contributes nothing below
# ~3800 A: HD138716, which lacks it, has the same relative noise as the other
# twelve from 3300 to 3800 A (ratio 0.88-1.11).
#
# THERE IS NO JOIN BLUEWARD OF THIS ONE, and that was tested rather than
# assumed. Over 3200-3740 A the stacked log-noise sits on a smooth 4th-order
# trend to 1.5% rms, and the largest step over any 30 A window anywhere in that
# range is +0.019 in ln -- against +0.51 for the real join at 3750-3850 A, so
# the join is 28x bigger than the largest thing that could hide below it. The
# steady fall from 1.44 at 3200 A to 0.59 at 3750 A is the instrument's blue
# throughput rising, not a splice. UVES-POP's delivered product simply starts
# at 3200 A inside the bluest setting, so the join at 3734-3859 A is the
# bluest one the library has.
#
# (lo, hi, label): lo is the bluest common dropout, hi the coverage edge.
SETTING_JOINS = [(3733.6, 3859.2, 'setting join'),
                 (4781.9, 4784.1, 'setting join')]

# Shift search. +/-150 km/s is +/-2.0 A at the break, which covers both NGSL's
# fitted wavecal residuals (0.07 to 1.22 A on the sample stars) and the range
# of stellar RV in this set (-9 to +120 km/s).
RV_GRID = np.arange(-150.0, 150.01, 5.0)
SHIFT_WIN = (3700.0, 4700.0)     # line-rich, inside G430L, avoids the break
SHIFT_WIN_WIDE = (3300.0, 5600.0)   # fallback: the whole of G430L
MIN_SHIFT_PX = 120               # usable NGSL pixels for a shift fit to count
CONT_PX = 101                    # median-filter width for continuum removal

# UVES-POP spectra carry NaN gaps -- mostly single bad pixels, but HD138716 has
# a 925 A hole at 3859-4784 A that swallows H-epsilon entirely. They cannot be
# left in: broaden_moffat convolves with an FFT, so ONE NaN anywhere poisons the
# whole segment and every panel would come out blank.
#
# Runs up to GAP_FILL_MAX are interpolated across, which is what they are --
# isolated dropouts, 0.6-0.9 A at worst in 12 of the 13 stars. Anything wider is
# a real hole and is blanked back out AFTER smoothing, out to the kernel's
# reach, so no pixel whose value was built partly from invented flux is drawn.
# Bridging a hole silently is exactly the failure common/specplot.py was fixed
# for: a flat segment at continuum level reads as an observation with no line.
GAP_FILL_MAX = 2.0                              # A
HOLE_MARGIN = NGSL_TRUNC_PX * 2.747             # 41 A, the G430L kernel reach


def ngsl_spectrum(fname):
    """-> (vacuum wavelength, flux) for one NGSL v2 file.

    air->vacuum only. `apply_wavecal` is deliberately NOT used: it has no entry
    for any of these 13 stars and would warn once per star while changing
    nothing. The zero point is fitted in `fit_shift` instead.
    """
    d = fits.getdata(ROOT / 'data' / 'spectra' / fname)
    w = air_to_vac(d['WAVELENGTH'].astype(float))
    f = d['FLUX'].astype(float)
    # Bad pixels are NaN'd rather than deleted. Dropping them would close the
    # grid over the gap, and every later step -- the plot, the rebin, the
    # residual -- would interpolate straight across it without saying so.
    return w, np.where(np.isfinite(f) & (f > 0), f, np.nan)


def load_xp(star):
    """-> (wavelength A, F_lambda cgs) for Gaia DR3 XP sampled, or None.

    Written by explore/fetch_gaia.py: 336-1020 nm on a 2 nm grid, already in
    VACUUM (Gaia has no air path) and already converted to erg/s/cm^2/A, so
    nothing is done to it here.

    NOT put on the NGSL wavelength scale, and NGSL is NOT degraded to it. XP is
    R ~ 30-100 -- far COARSER than NGSL's 3.54 A core, so the smoothing in this
    module runs the wrong way for it and there is no honest pixel-wise residual
    against NGSL without convolving NGSL down, which would make a different
    figure. It is drawn as delivered, on top.

    The rest-frame shift applied to NGSL is ignored for XP on purpose: the
    largest shift in this sample is 71 km/s, which is 0.9 A at the break,
    against an XP resolution element of order 100 A.
    """
    f = XP_DIR / f'{star}_xp.csv'
    if not f.exists():
        return None
    d = np.genfromtxt(f, delimiter=',', names=True)
    return (np.asarray(d['wavelength_A_vacuum'], float),
            np.asarray(d['flam'], float))


def fetch_uves(name, url):
    """Download one UVES-POP spectrum if it is not already on disk.

    The 13 are ~2.9 MB each and .gitignore keeps them out of the repository, so
    a clean checkout has none of them. Written to a temporary name and moved
    into place only on success, so an interrupted fetch cannot leave a
    zero-length file that later reads as a corrupt spectrum.
    """
    p = ROOT / 'data' / 'uves_pop' / f'{name}.fits.gz'
    if p.exists() and p.stat().st_size > 0:
        return p
    p.parent.mkdir(parents=True, exist_ok=True)
    tmp = p.with_suffix('.part')
    print(f'  fetching {name} ...', flush=True)
    req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
    with urllib.request.urlopen(req, timeout=300) as r, open(tmp, 'wb') as fh:
        fh.write(r.read())
    tmp.replace(p)
    print(f'    {p.stat().st_size / 1e6:.1f} MB -> {p.name}', flush=True)
    return p


def uves_meta(name):
    """-> (GCVS name, GCVS variability type, DATE-OBS) from the UVES-POP file.

    These headers are not FITS-standard -- cards out of order, and NaN written
    unquoted in a couple of dozen photometry keywords -- so the verify is needed
    and its ~80 warnings per file are silenced rather than printed.
    """
    with fits.open(ROOT / 'data' / 'uves_pop' / f'{name}.fits.gz') as f:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', VerifyWarning)
            f[0].verify('fix')
        date = str(f[0].header.get('DATE-OBS', ''))[:10]
        t, cn = f[2].data, f[2].columns.names
        g = lambda c: (str(t[c][0]).strip() if c in cn else '')
        return g('GCVS'), g('GCVS_TYPE'), date


def fill_small_gaps(w, f, max_A=GAP_FILL_MAX):
    """-> (flux with short NaN runs interpolated, [(lo, hi)] of the long ones)."""
    bad = ~np.isfinite(f)
    if not bad.any():
        return np.array(f, float), []
    out = np.array(f, float)
    idx = np.where(bad)[0]
    runs = np.split(idx, np.where(np.diff(idx) > 1)[0] + 1)
    good = ~bad
    out[bad] = np.interp(w[bad], w[good], f[good])
    holes = [(float(w[r[0]]), float(w[r[-1]])) for r in runs
             if float(w[r[-1]] - w[r[0]]) > max_A]
    return out, holes


def blank_near_holes(w, f, holes, margin=HOLE_MARGIN):
    """NaN out every pixel the kernel could have pulled invented flux into."""
    out = np.array(f, float)
    for lo, hi in holes:
        out[(w > lo - margin) & (w < hi + margin)] = np.nan
    return out


def normalised(w, f, win_px=CONT_PX):
    """Continuum-removed spectrum, for cross-correlation only.

    NaNs are bridged for the median filter and restored afterwards, so a gap
    neither poisons the continuum nor quietly acquires a value.
    """
    g = np.isfinite(f)
    if g.sum() < win_px:
        return np.full_like(f, np.nan)
    filled = np.where(g, f, np.interp(w, w[g], f[g]))
    c = median_filter(filled, size=win_px, mode='nearest')
    return np.where(g & (c > 0), f / c, np.nan)


def fit_shift(wn, fn, wu, fu_b, lo=SHIFT_WIN[0], hi=SHIFT_WIN[1]):
    """-> (shift km/s, sigma km/s) aligning NGSL to the smoothed UVES spectrum.

    `fu_b` is UVES ALREADY CONVOLVED with the NGSL kernel but not yet rebinned,
    so each trial only costs a rebin. Convolving once and shifting afterwards is
    exact here: the kernel is constant in Angstroms across G430L, so it commutes
    with a shift over this range. The pair the widths were fitted under --
    kernel then pixel integration -- is still what every trial applies.

    RETURNS THE CORRECTION TO APPLY TO NGSL, not the offset between them: the
    value v such that wn * (1 + v/c) puts NGSL onto the UVES rest frame. So a
    NEGATIVE shift means NGSL's features sit REDWARD of where they belong.

    The search itself runs the other way round -- it redshifts the UVES template
    until it matches NGSL -- and the result is negated on the way out. That
    negation is the whole content of this convention, and `selftest` injects a
    known shift specifically to pin its sign.
    """
    k = (wn > lo) & (wn < hi) & np.isfinite(fn)
    if k.sum() < MIN_SHIFT_PX:
        return np.nan, np.nan
    yn = normalised(wn[k], fn[k])
    chi, npix = [], []
    for v in RV_GRID:
        yu = rebin_to_pixels(wu * (1 + v / C_KMS), fu_b, wn[k])
        g = np.isfinite(yu) & np.isfinite(yn)
        if g.sum() < MIN_SHIFT_PX:
            chi.append(np.nan)
            npix.append(0)
            continue
        chi.append(float(np.nanmean((yn[g] - normalised(wn[k], yu)[g]) ** 2)))
        npix.append(int(g.sum()))
    chi = np.asarray(chi)
    if not np.isfinite(chi).any():
        return np.nan, np.nan
    j = int(np.nanargmin(chi))
    if j in (0, len(RV_GRID) - 1):          # minimum at the edge is a failure
        return np.nan, np.nan
    step = RV_GRID[1] - RV_GRID[0]
    y0, y1, y2 = chi[j - 1], chi[j], chi[j + 1]
    d = y0 - 2 * y1 + y2
    if not (d > 0):
        return np.nan, np.nan
    v_best = float(RV_GRID[j] + 0.5 * (y0 - y2) / d * step)

    # The metric is a MEAN SQUARED RESIDUAL, not a chi2: the two spectra carry
    # no common error model, so its absolute scale is arbitrary and the raw
    # curvature would give a nonsense error bar -- the first version of this
    # returned +/-1700 to +/-8200 km/s. Rescaling so the minimum equals the
    # number of pixels makes it a chi2 with the residual scatter standing in for
    # the noise, which is the usual convention when the errors are unknown.
    # It assumes the fit is a good one, so it is a PRECISION, not an accuracy:
    # it says nothing about the NGSL wavecal error it is measuring.
    # y1 == 0 is an exact match, which only happens in `selftest`, where the
    # NGSL side was built from the template itself. The shift is still defined;
    # the precision is simply not limited by scatter there.
    n = float(npix[j])
    sig = float(np.sqrt(2.0 / (d / y1 * n)) * step) if y1 > 0 else 0.0
    return -v_best, sig


def selftest():
    """Inject a known shift and recover it. Guards the sign, not just the size.

    A cross-correlation that returns the right magnitude with the wrong sign
    looks entirely plausible on a figure, and would put every H-epsilon core in
    this comparison on the wrong side of the line.
    """
    w = np.arange(3600.0, 4800.0, 0.05)
    f = np.ones_like(w)
    rng = np.random.default_rng(0)
    for c in rng.uniform(3650, 4750, 400):          # a fake line forest
        f -= 0.35 * np.exp(-0.5 * ((w - c) / 0.4) ** 2)
    wn = np.arange(3600.0, 4800.0, 2.747)           # NGSL G430L pixels
    fb = broaden_ngsl(w, f)
    for truth in (-60.0, -12.0, 0.0, 25.0, 90.0):
        # The fake NGSL spectrum is the template REDSHIFTED by `truth`, so the
        # correction that undoes it is -truth. Checking the returned number
        # would only test the magnitude, so the correction is then APPLIED and
        # the residual offset refitted: that is what the figures rely on, and
        # a sign error survives the first check but not the second.
        fn = rebin_to_pixels(w * (1 + truth / C_KMS), fb, wn)
        g = np.isfinite(fn)
        got, _ = fit_shift(wn[g], fn[g], w, fb)
        assert abs(got + truth) < 3.0, f'injected {truth}, recovered {got}'
        corrected = wn[g] * (1 + got / C_KMS)
        left, _ = fit_shift(corrected, fn[g], w, fb)
        assert abs(left) < 3.0, (f'injected {truth}, corrected by {got}, '
                                 f'{left} km/s left over -- sign is wrong')
    print('  selftest: injected shifts recovered and removed to < 3 km/s')


def prepare(star, row):
    """-> dict of everything one panel needs, or None if the star is unusable."""
    wn, fn = ngsl_spectrum(row['ngsl_file'])
    wu, fu, _ = load_uves(row['uves_name'])          # air->vacuum applied

    # 1. UVES to the rest frame with the catalogued RV, verified in the module
    #    docstring. HD058343 is the only star without one; SIMBAD gives it
    #    -4.5 km/s, which is 0.06 A at H-epsilon against a 3.54 A NGSL core, so
    #    leaving it at zero is below anything this figure can show.
    rv = float(row['rv_uves_kms']) if row['rv_uves_kms'] else 0.0
    wu_rest = wu / (1 + rv / C_KMS)

    # 2. NGSL's zero point, fitted against the smoothed UVES spectrum.
    fu_fill, holes = fill_small_gaps(wu_rest, fu)
    fu_b = broaden_ngsl(wu_rest, fu_fill)
    shift, shift_err = fit_shift(wn, fn, wu_rest, fu_b)
    if not np.isfinite(shift):
        # HD138716's 925 A hole leaves too little of the default window. Falling
        # back to the whole of G430L is worth doing before giving up, and which
        # window was used is reported rather than absorbed.
        shift, shift_err = fit_shift(wn, fn, wu_rest, fu_b, *SHIFT_WIN_WIDE)
    # fit_shift returns the CORRECTION, so this application is a plain multiply.
    wn_rest = wn * (1 + (shift if np.isfinite(shift) else 0.0) / C_KMS)

    # The final product goes through to_ngsl_pixels -- the one call that does
    # the kernel and the pixel together -- rather than reusing fu_b, so what is
    # plotted is what a caller elsewhere would get.
    ud = blank_near_holes(wn_rest, to_ngsl_pixels(wu_rest, fu_fill, wn_rest),
                          holes)
    gcvs, vtype, udate = uves_meta(row['uves_name'])
    xp = load_xp(star)
    return dict(star=star, row=row, wn=wn_rest, fn=fn, wu=wu_rest, fu=fu,
                ud=ud, rv=rv, shift=shift, shift_err=shift_err, holes=holes,
                xp=xp, gcvs=gcvs, vtype=vtype, uves_date=udate,
                ngsl_date=row.get('ngsl_obsdate', ''))


def epoch_note(d):
    """-> 'NGSL 2004 / UVES 2001' -- the two libraries are years apart."""
    n, u = (d['ngsl_date'] or '?')[:4], (d['uves_date'] or '?')[:4]
    return f'NGSL {n} / UVES {u}'


def covered(d, xlim, frac=0.5):
    """True if enough of the panel survives the holes to be worth drawing."""
    k = (d['wn'] > xlim[0]) & (d['wn'] < xlim[1])
    return k.sum() > 0 and np.isfinite(d['ud'][k]).mean() >= frac


def scale_to_ngsl(d, wins):
    """-> (UVES@NGSL on the NGSL flux scale, the grey factor, the mask used)."""
    k = np.zeros_like(d['wn'], bool)
    for lo, hi in wins:
        k |= (d['wn'] >= lo) & (d['wn'] <= hi)
    k &= np.isfinite(d['ud']) & (d['ud'] > 0) & (d['fn'] > 0)
    grey = float(np.median(d['ud'][k] / d['fn'][k])) if k.sum() > 10 else 1.0
    return d['ud'] / grey, grey, k


def panel(ax, rax, d, xlim, wins, marks, raw_alpha=0.85, raw_lw=0.5,
          raw_floor=False, show_xp=False):
    """One star: flux above, fractional residual below. -> (grey, rms%)."""
    uds, grey, kn = scale_to_ngsl(d, wins)
    resid = np.where(np.isfinite(uds) & (uds > 0),
                     (d['fn'] - uds) / uds * 100.0, np.nan)
    inw = (d['wn'] > xlim[0]) & (d['wn'] < xlim[1])
    iwu = (d['wu'] > xlim[0]) & (d['wu'] < xlim[1])

    for a in (ax, rax):
        a.set_facecolor(SURFACE)
        a.fill_between(d['wn'], 0, 1, where=kn & inw,
                       transform=a.get_xaxis_transform(), color=MUTED,
                       alpha=.13, lw=0, zorder=0)
        for m_, ls in marks:
            a.axvline(m_, color=MUTED, ls=ls, lw=1, zorder=1)
        # Setting joins, drawn only where they fall inside the panel.
        for lo, hi, _ in SETTING_JOINS:
            if hi < xlim[0] or lo > xlim[1]:
                continue
            a.axvspan(lo, hi, color=JOIN_C, alpha=.12, lw=0, zorder=0)
            for e in (lo, hi):
                a.axvline(e, color=JOIN_C, lw=0.9, alpha=.85, zorder=1)
        a.set_xlim(*xlim)
        a.grid(alpha=.25, color=GRID, lw=.7)
        a.tick_params(labelsize=7, colors=MUTED)
        for sp in a.spines.values():
            sp.set_color(GRID)

    if raw_alpha > 0:
        ax.plot(d['wu'][iwu], d['fu'][iwu] / grey, color=UVES_C, lw=raw_lw,
                alpha=raw_alpha, zorder=3, label='UVES-POP, as delivered')
    ax.plot(d['wn'][inw], d['fn'][inw], color=OBS_C, lw=1.3, zorder=4,
            label='NGSL v2')
    ax.plot(d['wn'][inw], uds[inw], color=UVESD_C, lw=1.4, ls='--', zorder=5,
            label='UVES-POP @ NGSL LSF + pixels')

    # Gaia XP, drawn AS DELIVERED on top. Only the flux panel: there is no XP
    # residual here, because XP is coarser than NGSL and a pixel-wise residual
    # would require degrading NGSL to it. The markers are the real 2 nm samples,
    # so the reader can see how few points carry the break.
    grey_xp = ''
    xp_resid = None
    if show_xp and d.get('xp') is not None:
        wx, fx = d['xp']
        onto = np.interp(d['wn'], wx, fx, left=np.nan, right=np.nan)
        kx = kn & np.isfinite(onto) & (onto > 0)
        if kx.sum() > 5:
            grey_xp = float(np.median(onto[kx] / d['fn'][kx]))
            mx = (wx > xlim[0]) & (wx < xlim[1])
            xps = fx[mx] / grey_xp                 # XP on the NGSL flux scale
            ax.plot(wx[mx], xps, color=XP_C, lw=1.5, marker='o',
                    ms=2.6, mew=0, zorder=6, alpha=.95,
                    label=f'Gaia DR3 XP, as delivered  (/{grey_xp:.3f})')
            # NGSL smoothed DOWN to XP and integrated onto XP's own samples.
            # This is the one place the smoothing runs the other way, and it is
            # a different operation from to_ngsl_pixels: a plain Gaussian of
            # measured width, not the NGSL Moffat, because the target is Gaia's
            # profile and not NGSL's.
            gn = np.isfinite(d['fn'])
            ng_xp = rebin_to_pixels(d['wn'][gn],
                                    broaden(d['wn'][gn], d['fn'][gn],
                                            XP_FWHM_A, step=2.0), wx[mx])
            xp_resid = np.where(np.isfinite(ng_xp) & (xps > 0),
                                (ng_xp - xps) / xps * 100.0, np.nan)

    # The top comes from the SMOOTHED curves only. Letting the raw spectrum set
    # it would hand the axis to the line forest between the Balmer members,
    # which peaks well above any continuum and would flatten everything else.
    vals = np.concatenate([d['fn'][inw], uds[inw][np.isfinite(uds[inw])]])
    ylo, yhi = np.nanpercentile(vals, [0.5, 99.5])
    if raw_floor and raw_alpha > 0:
        # ...but the bottom comes from the RESOLVED core, which is the whole
        # point of the zoom: at NGSL's 3.54 A the core is filled in by the
        # profile, and at UVES's 0.1 A grid it is not, so the depth difference
        # between the solid and the green line IS the measurement. Median
        # filtered first so a single cosmic-ray pixel cannot set the axis.
        raw = d['fu'][iwu] / grey
        raw = raw[np.isfinite(raw)]
        if raw.size > 20:
            ylo = min(ylo, float(np.min(median_filter(raw, size=5,
                                                      mode='nearest'))))
    pad = .10 * (yhi - ylo)
    ax.set_ylim(ylo - pad, yhi + pad)

    rax.axhline(0, color=MUTED, lw=1, zorder=2)
    rax.fill_between(d['wn'][inw], resid[inw], 0, color=OBS_C, alpha=.28,
                     lw=0, zorder=3)
    rax.plot(d['wn'][inw], resid[inw], color=OBS_C, lw=1.0, zorder=4,
             label='NGSL − UVES@NGSL')
    if xp_resid is not None:
        rax.plot(wx[mx], xp_resid, color=XP_C, lw=1.3, marker='o', ms=2.2,
                 mew=0, zorder=6, label=f'NGSL@XP − XP  ({XP_FWHM_A:.0f} Å)')
    fin = resid[inw][np.isfinite(resid[inw])]
    rms = float(np.sqrt(np.nanmean(fin ** 2))) if fin.size else np.nan
    # Range set by the structure, not by a fixed number, but never tighter than
    # +/-5% -- a panel autoscaled to a flat residual magnifies its own noise.
    allr = fin if xp_resid is None else np.concatenate(
        [fin, xp_resid[np.isfinite(xp_resid)]])
    if allr.size:
        m = max(5.0, float(np.nanpercentile(np.abs(allr), 99)) * 1.25)
        rax.set_ylim(-m, m)
    return grey, rms, grey_xp


def make_figure(data, xlim, wins, marks, title, sub, outname,
                raw_alpha=0.85, raw_lw=0.5, raw_floor=False, show_xp=False):
    """-> {star: (grey, rms%)} for the stars actually drawn."""
    ncol = 3
    nrow = int(np.ceil(len(data) / ncol))
    # The header grew to six lines; reserving a fixed FRACTION of the figure
    # let it run into the first row of panel titles. Reserve INCHES instead,
    # from the actual line count, and add them to the figure height so the
    # panels keep the same size however long the caption gets.
    nlines = 1 + f'{title}\n{sub}'.count('\n')
    head_in = 0.45 + nlines * 11.5 * 1.7 / 72.0
    fig_h = 4.1 * nrow + head_in + 0.62
    fig = plt.figure(figsize=(16.5, fig_h))
    fig.patch.set_facecolor(SURFACE)
    gs = fig.add_gridspec(nrow * 2, ncol, height_ratios=[2.2, 1] * nrow,
                          hspace=.42, wspace=.20)
    stats, xp_grey = {}, {}
    legend_done = False
    for i, d in enumerate(data):
        r, c = divmod(i, ncol)
        ax = fig.add_subplot(gs[2 * r, c])
        rax = fig.add_subplot(gs[2 * r + 1, c], sharex=ax)
        grey, rms, grey_xp = panel(ax, rax, d, xlim, wins, marks, raw_alpha,
                                   raw_lw, raw_floor, show_xp)
        stats[d['star']] = (grey, rms)
        if grey_xp != '':
            xp_grey[d['star']] = grey_xp
        sh = (f"{d['shift']:+.0f} km/s" if np.isfinite(d['shift'])
              else 'shift FAILED')
        var = f"   {d['gcvs']} ({d['vtype']})" if d['vtype'] else ''
        ax.set_title(f"{d['star']}   {d['row']['sptype'] or '—'}{var}\n"
                     f"grey /{grey:.3f}   NGSL {sh}   {epoch_note(d)}",
                     fontsize=8, color=(INK if not d['vtype'] else '#8a3324'))
        ax.tick_params(labelbottom=False)
        if c == 0:
            ax.set_ylabel(r'F$_\lambda$ (NGSL scale)', fontsize=8, color=INK)
            rax.set_ylabel('(NGSL − UVES)/UVES [%]', fontsize=8, color=INK)
        # The legend goes on the first panel that actually shows every series,
        # otherwise the XP entry would be missing from it.
        if not legend_done and (not show_xp or d.get('xp') is not None):
            ax.legend(fontsize=6.5, loc='best', framealpha=.92)
            if show_xp:
                rax.legend(fontsize=6, loc='best', framealpha=.92, ncol=2)
            legend_done = True
        if i >= len(data) - ncol:
            rax.set_xlabel(r'Wavelength [$\AA$, vacuum, rest]', fontsize=8,
                           color=INK)
    fig.suptitle(f'{title}\n{sub}', fontsize=11.5, color=INK,
                 linespacing=1.7, y=1 - 0.18 / fig_h, va='top')
    # Bottom reserve in INCHES too, for the same reason as the header: the
    # x-label needs a fixed amount of space, not a fraction of a figure
    # whose height changes with the number of rows and caption lines.
    fig.subplots_adjust(left=.055, right=.985, top=1 - head_in / fig_h,
                        bottom=0.62 / fig_h)
    out = library_figure_path(outname)
    fig.savefig(out, dpi=150, facecolor=SURFACE)
    plt.close(fig)
    n_xp = len(xp_grey)
    print(f'  {len(data)} panels'
          f'{f" ({n_xp} with Gaia XP)" if n_xp else ""} -> {out}')
    return stats


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--skip-selftest', action='store_true')
    # Gaia XP is OFF by default. It was tried on the break figure and it did not
    # earn its place: XP is R ~ 52 there against NGSL's 3.54 A core, its range
    # starts only 286 A blueward of the break, and what the ratio panel showed
    # was mostly XP's own calibration pattern -- three of the five stars give
    # the same -1.5% blue-to-red tilt, which is the instrumental XP/NGSL shape
    # already measured in explore/xp_vs_ngsl.py, not anything about the stars.
    # The code is kept because the retrieval and the measured 70 A width are
    # reusable, not because the figure was better with it.
    ap.add_argument('--xp', action='store_true',
                    help='overlay Gaia DR3 XP on the break figure (5 stars)')
    a = ap.parse_args()
    if not a.skip_selftest:
        selftest()

    rows = list(csv.DictReader(open(ROOT / 'data' / 'uves_ngsl_overlap.csv')))
    # The overlap table does not carry observation dates; the NGSL catalog does,
    # and the gap between the two epochs is what a variable star's panel has to
    # be read against.
    obsdate = {r['target']: r['obsdate'] for r in
               csv.DictReader(open(ROOT / 'data' / 'ngsl_catalog.csv'))}
    for r in rows:
        r['ngsl_obsdate'] = obsdate.get(r['ngsl_target'], '')
    data = []
    for r in rows:
        star = r['ngsl_target']
        try:
            fetch_uves(r['uves_name'], r['uves_spec_url'])
        except Exception as e:                       # noqa: BLE001 - report, skip
            print(f'  {star}: UVES spectrum unavailable ({e}), skipped')
            continue
        d = prepare(star, r)
        sh = (f"{d['shift']:+6.1f} +/- {d['shift_err']:.1f}"
              if np.isfinite(d['shift']) else '        FAILED')
        # Only holes that touch a panel are worth printing. Every UVES-POP
        # spectrum has a dozen inter-order gaps beyond ~6800 A, and reporting
        # those made all 13 stars look equally affected when only one is.
        near = [(lo, hi) for lo, hi in d['holes']
                if hi > BREAK_XLIM[0] - HOLE_MARGIN
                and lo < HEPS_XLIM[1] + HOLE_MARGIN]
        hole = (f"   {len(near)} hole(s) in range, widest "
                f"{max(hi - lo for lo, hi in near):.0f} A" if near else '')
        print(f"  {star:<10} RV {d['rv']:+7.2f} km/s   "
              f"NGSL shift {sh} km/s{hole}")
        data.append(d)
    print(f'{len(data)} stars prepared\n')

    # A star whose UVES coverage is mostly hole in a panel is left OUT of that
    # figure rather than drawn as a blank cell: HD138716's 925 A gap takes out
    # H-epsilon completely while leaving its Balmer break untouched.
    brk_set = [d for d in data if covered(d, BREAK_XLIM)]
    hep_set = [d for d in data if covered(d, HEPS_XLIM)]
    brk_names = {d['star'] for d in brk_set}
    hep_names = {d['star'] for d in hep_set}
    for d in data:
        for names, where, xl in ((hep_names, 'H-epsilon', HEPS_XLIM),
                                 (brk_names, 'the Balmer break', BREAK_XLIM)):
            if d['star'] in names:
                continue
            k = (d['wn'] > xl[0]) & (d['wn'] < xl[1])
            print(f"  {d['star']}: dropped from {where} — UVES covers only "
                  f"{np.isfinite(d['ud'][k]).mean():.0%} of that window "
                  f"({len(d['holes'])} hole(s), widest "
                  f"{max((hi - lo for lo, hi in d['holes']), default=0):.0f} A)")

    brk = make_figure(
        brk_set, BREAK_XLIM, BREAK_WINS, [(BALMER, '--')],
        f'Balmer break ({BALMER:.0f} A), UVES-POP against NGSL',
        'Shaded: the normalisation window, blueward of the break only.\n'
        'Residual positive = NGSL is BRIGHTER than UVES-POP smoothed to the '
        'NGSL LSF and integrated onto NGSL pixels.\n'
        'Panel titles in red are catalogued variables (GCVS type in brackets), '
        'observed years apart by the two libraries.\n'
        'Vertical gold lines at 3733.6 and 3859.2 A bracket the UVES setting '
        'join — measured here, and the bluest one the library has.'
        + ('\nPink: Gaia DR3 XP where it exists (5 of 13), as delivered; the '
           'ratio panel adds NGSL smoothed DOWN to XP (70 Å Gaussian).'
           if a.xp else ''),
        'uves_ngsl_break.png', raw_alpha=0.45, raw_lw=0.4, show_xp=a.xp)
    hep = make_figure(
        hep_set, HEPS_XLIM, HEPS_WINS, [(HEPS, '--'), (CAH, ':')],
        'H-epsilon (3971.2 A), UVES-POP against NGSL',
        "dashed: H-epsilon (3971.2 A).   dotted: interstellar Ca II H "
        "(3969.6 A), 1.64 A away and UNRESOLVED by NGSL's 3.54 A core —\n"
        "the NGSL feature is the BLEND of the two, so no depth read off an "
        "NGSL pixel here is a hydrogen measurement.\n"
        'Residual positive = NGSL is BRIGHTER than UVES-POP smoothed to NGSL.  '
        'Panel titles in red are catalogued variables.\n'
        'The y-axis reaches down to the UVES-POP core at native resolution; '
        'its top is set by the smoothed curves.',
        'uves_ngsl_hepsilon.png', raw_floor=True)

    with open(ROOT / 'data' / 'uves_ngsl_compare.csv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=[
            'star', 'sptype', 'rv_uves_kms', 'ngsl_shift_kms',
            'ngsl_shift_err_kms', 'ngsl_shift_A_at_3971', 'grey_break',
            'rms_break_pct', 'grey_hepsilon', 'rms_hepsilon_pct',
            'gcvs_name', 'gcvs_type', 'ngsl_obsdate', 'uves_obsdate',
            'notes'])
        w.writeheader()
        for d in data:
            sh = d['shift']
            g_b, r_b = brk.get(d['star'], ('', ''))
            g_h, r_h = hep.get(d['star'], ('', ''))
            w.writerow(dict(
                star=d['star'], sptype=d['row']['sptype'],
                rv_uves_kms=round(d['rv'], 2),
                ngsl_shift_kms=round(sh, 1) if np.isfinite(sh) else '',
                ngsl_shift_err_kms=(round(d['shift_err'], 1)
                                    if np.isfinite(d['shift_err']) else ''),
                ngsl_shift_A_at_3971=(round(HEPS * sh / C_KMS, 3)
                                      if np.isfinite(sh) else ''),
                grey_break=round(g_b, 4) if g_b != '' else '',
                rms_break_pct=round(r_b, 2) if r_b != '' else '',
                grey_hepsilon=round(g_h, 4) if g_h != '' else '',
                rms_hepsilon_pct=round(r_h, 2) if r_h != '' else '',
                gcvs_name=d['gcvs'], gcvs_type=d['vtype'],
                ngsl_obsdate=d['ngsl_date'], uves_obsdate=d['uves_date'],
                notes=';'.join(filter(None, [
                    '' if d['star'] in hep_names else 'uves_gap_at_hepsilon',
                    f"variable_{d['vtype']}" if d['vtype'] else '']))))
    print('  -> data/uves_ngsl_compare.csv')

    print(f'\n{"star":<10}{"shift km/s":>12}{"A@3971":>9}{"grey":>8}'
          f'{"rms break %":>13}{"rms Heps %":>12}')
    for d in data:
        sh = d['shift']
        g_b, r_b = brk.get(d['star'], ('', ''))
        _, r_h = hep.get(d['star'], ('', ''))
        sa = f'{HEPS * sh / C_KMS:+.2f}' if np.isfinite(sh) else '—'
        ss = f'{sh:+.1f}' if np.isfinite(sh) else 'FAILED'
        print(f'{d["star"]:<10}{ss:>12}{sa:>9}'
              f'{(f"{g_b:.3f}" if g_b != "" else "—"):>8}'
              f'{(f"{r_b:.2f}" if r_b != "" else "—"):>13}'
              f'{(f"{r_h:.2f}" if r_h != "" else "no cover"):>12}')


if __name__ == '__main__':
    main()
