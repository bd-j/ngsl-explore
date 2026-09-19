"""The Balmer series as a diagnostic: one window per member, model against XSL.

Shared by `explore/check_predict.py` (H10 and H11, which sit inside the
held-out NGSL break window) and `explore/plot_balmer_lines.py` (the whole
series). It exists for the same reason `common/specplot.py` does: the two
figures must not grow separate ideas of where a line's window is or what its
core residual means.

WHY THE HIGH ORDERS ARE WORTH A PANEL AT ALL

The held-out NGSL Balmer window carries a line-core excess of order +5 to +10%,
and that number is DEGENERATE WITH THE INSTRUMENT PROFILE: at HD194453 it runs
+2.65% at a 3.54 A Moffat core to -4.20% at 7.00 A (docs/LSF.md). Measured at
NGSL's 1.4 A pixels it therefore cannot separate hydrogen physics from an error
in the LSF. XSL resolves the same lines with a 0.39 A FWHM at H10 -- sigma = 13
km/s in the UVB, so ~10x narrower than the NGSL core and a tenth as sensitive
to getting it wrong -- which is what makes it an independent check rather than
a prettier version of the same measurement.

WHAT IS FITTED AND WHAT IS NOT

  H-alpha .. H-delta   FITTED. XSL conditions on these four, +/-50 A with the
                       core masked (fitting/observations.py). Their panels show
                       the fit's OWN calibrated model, never a second solve
                       here, so what is drawn is what the likelihood saw.
  H-epsilon and up     NOT FITTED, by anything. They blend into one another, so
                       no local continuum is defined for them, and they sit
                       inside the held-out NGSL window. Every pixel of their
                       panels is a prediction from both datasets.

For the unfitted members the model still has to be put on a local scale before
a residual means anything -- XSL is ground-based and its continuum is
marginalised away precisely because slit losses make it untrustworthy. That
scale is solved HERE, on the wings, with the core masked out, and it is not a
continuum: there is no continuum to find between blended lines. It is an
order-1 correction to the RATIO of observation to model, and the model already
carries the same blended neighbours. It costs two parameters out of the window
and feeds nothing back into any fit.
"""
import sys
from dataclasses import replace
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common.lines import BALMER, balmer_member, balmer_neighbourhood

# The members XSL actually conditions on. Written as indices and CHECKED
# against the named wavelengths below, so the two cannot drift apart -- the
# same failure mode that once let the fitter and the figures disagree about the
# NGSL LSF.
FITTED_N = (3, 4, 5, 6)                  # H-alpha, H-beta, H-gamma, H-delta
for _n, _name in zip(FITTED_N, ('Halpha', 'Hbeta', 'Hgamma', 'Hdelta')):
    assert abs(balmer_member(_n) - BALMER[_name]) < 0.2, _name

# How far up the series to go. H13 leaves 17 usable wing pixels and H14 leaves
# 2, against ~95 at H11, so the local scale stops being measurable rather than
# merely noisy. The limit is enforced per star anyway by MIN_WING_PX below --
# this is only where to stop looking.
NMAX = 12

# Core mask, as a multiple of the MEASURED half-depth half-width of the
# observed core. A fixed width cannot work across the series: the measured
# half-width runs 6.2 A at H8 to 3.4 A at H12 on this sample, so the 6.0 A the
# fit uses at H-delta is too narrow at H8 and nearly twice too wide at H12.
#
# Measured from the OBSERVATION, never scaled from the fitted v sin i, and this
# sample shows why. HD106304, HD128801 and HD117880 are all fitted at 150-200
# km/s, which is 3.3-4.4 A of broadening at H-alpha -- but their observed
# H-alpha cores are 2.4-2.6 A half-width and 0.23-0.25 deep, as sharp as
# HD194453's at v sin i = 0. The data says these are not fast rotators; the
# fitted value is the scan spending v sin i on a metallicity below the grid
# floor (see common/figpath.py). A mask taken from the fit would have put a
# 10 A window around a 2.5 A core on exactly the three stars whose fits are
# least trustworthy, and the core residual is measured inside it.
CORE_MASK_FACTOR = 1.3
# ... capped at half the NARROWER window half-width, so the mask can never eat
# the wings it is there to preserve.
CORE_MASK_MAX_FRAC = 0.5
MIN_WING_PX = 40                # below this the order-1 scale is not measurable


def window(n):
    """(lo, hi) for Balmer member n.

    The fitted members keep the fit's own +/-50 A so the panel and the
    likelihood cannot disagree; the rest get the midpoint window, which is the
    only boundary blended lines define. See common.lines.balmer_neighbourhood.
    """
    if n in FITTED_N:
        from fitting.observations import XSL_BALMER_HALFWIDTH
        lam = balmer_member(n)
        return lam - XSL_BALMER_HALFWIDTH, lam + XSL_BALMER_HALFWIDTH
    return balmer_neighbourhood(n)


def good_pixels(obs):
    """Pixels with finite, positive flux and uncertainty.

    NOT obs.mask: the whole point of these panels is the wavelengths the mask
    excludes, so the only cut here is 'is there data'.
    """
    w = np.asarray(obs.wavelength, float)
    f = np.asarray(obs.flux, float)
    e = np.asarray(obs.uncertainty, float)
    return (np.isfinite(w) & np.isfinite(f) & np.isfinite(e)
            & (f > 0) & (e > 0))


SMOOTH_PX = 5           # running median before the width is measured
DEPTH_NEAR_A = 2.0      # the depth is read this close to the line centre


def _running_median(y, k=SMOOTH_PX):
    """Median filter, edges left alone. ~0.6 A at XSL's 0.127 A pixels."""
    y = np.asarray(y, float)
    if y.size < k or k < 3:
        return y
    out = y.copy()
    win = np.lib.stride_tricks.sliding_window_view(y, k)
    out[k // 2:k // 2 + win.shape[0]] = np.median(win, axis=1)
    return out


def core_halfwidth(obs, n, good=None):
    """Half-depth half-width of the OBSERVED core, in Angstroms, or None.

    Measured against the 95th percentile of the window rather than a fitted
    continuum, because between blended lines there is no continuum to fit. That
    makes this a robust width, not a physical one -- which is all a mask needs.

    Three things here are not the obvious implementation, and each replaced a
    version that misbehaved on this sample:

      * the profile is median-filtered first. Without it the depth came from a
        single pixel, and on the low-S/N stars that pixel is a noise spike.
      * the depth is read from the pixels NEAR THE CENTRE, not from the window
        minimum, which on a blended window can belong to a neighbour.
      * the width is found by walking OUTWARD from the centre to the first
        crossing back above half depth, not from the extent of every pixel
        below it. The extent version let any deep metal line elsewhere in the
        window set the answer: it returned 21.8 A for one star's H-delta
        against 7.6-9.6 A for the other eleven, and 2.1 to 10.4 A across the
        sample at H-alpha. A mask that varies 5x between stars for no physical
        reason is not a mask, and it feeds the core residual directly.
    """
    lam, (lo, hi) = balmer_member(n), window(n)
    w = np.asarray(obs.wavelength, float)
    good = good_pixels(obs) if good is None else good
    sel = good & (w >= lo) & (w <= hi)
    if sel.sum() < 10:
        return None
    ww = w[sel]
    prof = _running_median(np.asarray(obs.flux, float)[sel])
    prof = prof / np.percentile(prof, 95)

    near = np.abs(ww - lam) < DEPTH_NEAR_A
    if near.sum() < 3:
        return None
    i0 = int(np.flatnonzero(near)[np.argmin(prof[near])])
    half = 1.0 - 0.5 * (1.0 - float(prof[i0]))

    left = np.flatnonzero(prof[:i0 + 1] > half)
    right = np.flatnonzero(prof[i0:] > half)
    i_lo = left[-1] + 1 if left.size else 0
    i_hi = i0 + right[0] - 1 if right.size else prof.size - 1
    if i_hi <= i_lo:
        return None
    return 0.5 * float(ww[i_hi] - ww[i_lo])


def core_mask_width(obs, n, good=None):
    """The half-width to mask around member n, in Angstroms."""
    lam, (lo, hi) = balmer_member(n), window(n)
    hw = core_halfwidth(obs, n, good)
    if hw is None:
        return None
    cap = CORE_MASK_MAX_FRAC * min(lam - lo, hi - lam)
    return float(min(CORE_MASK_FACTOR * hw, cap))


def member_panel(obs, model_value, n, fit_cal=None):
    """XSL around Balmer member n, model put on a local scale -> dict or None.

    Returns
      n, lam, window, core       the geometry, all derived
      view                       an Observation whose mask is the pixels that
                                 SET the scale, so common.specplot's
                                 solid-fitted / faded-masked convention draws
                                 the right thing with no extra argument
      cal                        the calibrated model, full length
      resid                      (obs - cal) / cal, NaN outside the window
      used, pred, incore         boolean masks: scale-setting pixels, predicted
                                 pixels, and the masked core within them
      source                     'fit' (the likelihood's own calibration) or
                                 'local' (solved here, on the wings)

    `fit_cal` is the calibrated model from the fit itself. It is used for the
    four fitted members and ignored for the rest, which nothing fits.
    """
    from fitting.calibration import solve, segment_pixels
    from fitting.observations import XSL_BALMER_ORDER

    lam, (lo, hi) = balmer_member(n), window(n)
    w = np.asarray(obs.wavelength, float)
    good = good_pixels(obs)
    inwin = good & (w >= lo) & (w <= hi)
    if inwin.sum() < 20:
        return None
    core = core_mask_width(obs, n, good)
    if core is None:
        return None
    # Two widths, deliberately different.
    #
    #   core     what the local scale is NOT allowed to see. 1.3x the measured
    #            half-depth width, so the solve has margin against a core that
    #            is broader than it looks.
    #   core_hw  what the quoted residual is measured over: the core itself.
    #
    # They were one number at first, and the statistic was wrong for it. The
    # median over the full 1.3x mask returned EXACTLY ZERO for a +5% excess
    # injected into the inner +/-3 A of H9 and H10, because that excess covered
    # only 39% and 47% of the masked region and a median ignores a minority --
    # while the same injection at H11 and H12, where the mask is narrower, read
    # back +5.000%. A statistic that finds an injected signal at one member and
    # not at the next is useless for the one thing this figure is for, which is
    # comparing the members against each other.
    #
    # A median, not a mean: these windows are peppered with metal lines the
    # models get wrong by 10% over a pixel or two, and one of them should not
    # be allowed to set a core residual.
    core_hw = min(core_halfwidth(obs, n, good) or core, core)

    if n in FITTED_N and fit_cal is not None:
        cal = np.asarray(fit_cal, float)
        used = inwin & np.asarray(obs.mask, bool)
        view, source = obs, 'fit'
    else:
        seg = dict(name=f'H{n}', order=XSL_BALMER_ORDER, domain=(lo, hi),
                   ranges=[(lo, lam - core), (lam + core, hi)])
        used = good & segment_pixels(w, seg)
        if used.sum() < MIN_WING_PX:
            return None
        # One segment per call, so two adjacent members -- whose windows share
        # an edge by construction -- are solved independently and can never
        # claim the same pixel. check_segments would reject them in one solve.
        view = replace(obs, mask=used, calibration=('segments', [seg]))
        # fill_domains puts the same polynomial across the masked core, which
        # is the only reason there is a model to draw through the line at all.
        cal, _ = solve(view, model_value, fill_domains=True)
        source = 'local'

    ok = inwin & np.isfinite(cal) & (cal > 0)
    resid = np.full(w.shape, np.nan)
    resid[ok] = (np.asarray(obs.flux, float)[ok] - cal[ok]) / cal[ok]
    used = used & ok
    pred = ok & ~used
    # `incore` is intersected with `pred`, not taken from the measured width
    # alone. The measured mask is wider than the fit's fixed 6 A at H-beta
    # through H-delta, so without this the 'core' statistic would quietly
    # include pixels the likelihood had already fitted -- a prediction that was
    # partly a fit, which is the one thing these panels exist to keep apart.
    return dict(n=n, lam=lam, window=(lo, hi), core=core, core_hw=core_hw,
                view=view, cal=cal, resid=resid, used=used, pred=pred,
                incore=pred & (np.abs(w - lam) < core_hw), source=source)


def core_stat(d):
    """-> (median %, rms %, n) of the residual over the core, +/- `core_hw`.

    That region sits strictly inside the masked `core`, so no scale was solved
    on any of it and this is a prediction -- and because the wings were
    levelled by the order-1 term, it is already a core-MINUS-wing excess, the
    same quantity explore/check_predict.py's core_excess() reports for NGSL.
    The two are still not the same NUMBER: the NGSL one is measured under a
    single scalar from the bands and carries a continuum error this one has
    already absorbed.
    """
    sel = d['incore'] & np.isfinite(d['resid'])
    if not sel.any():
        return None
    r = 100 * d['resid'][sel]
    return float(np.median(r)), float(np.std(r)), int(sel.sum())


def wing_stat(d):
    """-> median % of the residual on the scale-setting wings.

    Should be ~0 by construction for a 'local' panel. It is printed anyway: a
    value that is not near zero means the order-1 term could not level the
    window, which is the signature of a model wing shape that is wrong.
    """
    sel = d['used'] & np.isfinite(d['resid'])
    if not sel.any():
        return None
    return float(np.median(100 * d['resid'][sel]))


def merge_members(ds):
    """Adjacent members drawn as ONE panel -> a dict of the same shape.

    H10 and H11 share a window edge by construction (the midpoint of the gap
    between them), so they tile a single stretch of spectrum and are much
    easier to read side by side than in two panels. They are still SOLVED
    separately -- one local scale each, in `member_panel` -- because a single
    order-1 term across both would be asked to follow the blend structure of
    three lines at once. This only stitches the results for drawing.
    """
    ds = [d for d in ds if d is not None]
    if not ds:
        return None
    ds = sorted(ds, key=lambda d: d['lam'])
    w = np.asarray(ds[0]['view'].wavelength, float)
    cal = np.full(w.shape, np.nan)
    resid = np.full(w.shape, np.nan)
    used = np.zeros(w.shape, bool)
    pred = np.zeros(w.shape, bool)
    incore = np.zeros(w.shape, bool)
    for d in ds:
        lo, hi = d['window']
        inwin = (w >= lo) & (w <= hi) & np.isnan(cal)   # first claim wins
        cal[inwin] = np.asarray(d['cal'], float)[inwin]
        resid[inwin] = d['resid'][inwin]
        used |= d['used'] & inwin
        pred |= d['pred'] & inwin
        incore |= d['incore'] & inwin
    return dict(members=ds, lam=[d['lam'] for d in ds],
                window=(ds[0]['window'][0], ds[-1]['window'][1]),
                view=replace(ds[0]['view'], mask=used), cal=cal, resid=resid,
                used=used, pred=pred, incore=incore,
                source=ds[0]['source'] if len({d['source'] for d in ds}) == 1
                else 'mixed')
