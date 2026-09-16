"""Is the NGSL Balmer-core excess NLTE, or is it NGSL's line spread function?

The held-out Balmer window shows the observed high-order Balmer lines sitting
+3 to +10% ABOVE the models -- shallower cores than the models predict. That was
recorded as "almost certainly NLTE in hydrogen, which SYNTHE does not treat for
H" (CAVEATS.md). The evidence against that reading is simple: at XSL's R ~ 9800
the SAME models fit the FULL H-gamma profile, core included, for the metal-rich
stars. A physical NLTE core deficit would show up at both resolutions.

So the excess is probably instrumental, and this asks the question with NO model
in it at all:

    take XSL, degrade it to NGSL's resolution, rebin to NGSL pixels,
    and compare directly against NGSL.

Both are observations of the same star. If NGSL's cores come out shallower than
XSL's do after degradation, the cores are being filled in by the instrument, and
the model is exonerated.

The kernel width is fitted in a CONTROL window (4200-4600 A) that excludes the
Balmer region entirely, then applied to it. So the Balmer comparison is a
prediction of the control fit, not a fit to the thing being tested -- otherwise
a free width would simply absorb whatever it found there.

HISTORICAL NOTE, and a correction. This script was written to test whether a
winged profile explains the core excess, and an earlier version of it concluded
that one "removes 86%" of it. THAT CONCLUSION WAS WRONG. The core excess is
degenerate with the profile's effective WIDTH, not diagnostic of its SHAPE:
holding all else fixed, a plain Gaussian spans +16.8% at 3.85 A to -1.2% at
7.0 A. Any profile can be tuned to zero it, so it is evidence for none of them.

The profile actually adopted (common.lsf.broaden_ngsl_moffat) was chosen on
different grounds -- control-window rms at free width, and the fitted core
matching the independently published STIS value -- see explore/ngsl_lsf_shape.py
and docs/CAVEATS.md. What this script still establishes, and what it was worth
writing for, is the MODEL-FREE half: NGSL's cores are only ~0.8% shallower than
XSL degraded to the same resolution, so whatever fills them is instrumental and
not NLTE.

    python3 explore/ngsl_core_excess.py
"""
import csv
import sys
from pathlib import Path

import numpy as np
from scipy.ndimage import gaussian_filter1d
from scipy.optimize import minimize
from scipy.special import erfc

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from fitting.model import Grid
from fitting.observations import load_ngsl, BREAK_WINDOW
from fitting.predict import predict, spectrum_at
from common.extinction_ccm import redden
from common.lsf import broaden_rot
from fitting.scan import best_node
from common.xsl_load import load as load_xsl, resolving_power
from common.lsf import rebin_to_pixels
from common.lines import hydrogen_lines
from common.figpath import below_grid

ROOT = Path(__file__).resolve().parent.parent
CONTROL = (4200., 4600.)          # line-rich, and outside the held-out window
CORE_HALF = 3.0                   # A; "core" = within this of a Balmer line
# Only the RESOLVED high-order lines. Below ~3700 A the series crowds together
# faster than NGSL's ~6 A profile can separate, so every pixel there is "core"
# and the core/continuum contrast stops meaning anything. H7-H12 it is.
TEST = (3700., 4000.)
NMAX = 12
DEG = 3                           # continuum polynomial order per window


def cnorm(w, f, deg=DEG):
    ok = np.isfinite(f) & (f > 0)
    if ok.sum() < deg + 2:
        return np.full_like(np.asarray(f, float), np.nan)
    return f / np.polyval(np.polyfit(w[ok], f[ok], deg), w)


def kernel(w, f, fwhm_A, wing=0.0, wing_scale=4.0, step=0.02):
    """Gaussian of `fwhm_A`, optionally plus a broad Lorentzian-ish wing.

    `wing` is the fraction of the kernel's area in the broad component, which
    is a Gaussian `wing_scale` times wider -- a cheap stand-in for extended
    wings that keeps everything analytic and positive.
    """
    wl = np.arange(w[0], w[-1], step)
    fi = np.interp(wl, w, f)
    s = fwhm_A / 2.3548 / step
    core = gaussian_filter1d(fi, s, mode='nearest')
    if wing <= 0:
        out = core
    else:
        broad = gaussian_filter1d(fi, s * wing_scale, mode='nearest')
        out = (1.0 - wing) * core + wing * broad
    return np.interp(w, wl, out)


def degrade(wx, fx, wn, fwhm_A, wing=0.0):
    """XSL -> NGSL sampling at a given kernel."""
    sm = kernel(wx, fx, fwhm_A, wing)
    return rebin_to_pixels(wx, sm, wn)


def fit_control(wx, fx, wn, fn, fit_wing=False):
    """-> (fwhm, wing, rms) fitted on CONTROL only."""
    lo, hi = CONTROL
    sn = (wn > lo) & (wn < hi) & np.isfinite(fn) & (fn > 0)
    sx = (wx > lo - 60) & (wx < hi + 60) & np.isfinite(fx) & (fx > 0)
    if sn.sum() < 20 or sx.sum() < 200:
        return np.nan, np.nan, np.nan
    yn = cnorm(wn[sn], fn[sn])

    def rms(p):
        fw = abs(p[0])
        wg = float(np.clip(p[1], 0.0, 0.6)) if fit_wing else 0.0
        d = degrade(wx[sx], fx[sx], wn[sn], fw, wg)
        yd = cnorm(wn[sn], d)
        r = yn - yd
        return float(np.sqrt(np.nanmean(r ** 2)))

    p0 = [7.0, 0.15 if fit_wing else 0.0]
    res = minimize(rms, p0, method='Nelder-Mead',
                   options=dict(xatol=1e-3, fatol=1e-6, maxiter=400))
    fw = abs(res.x[0])
    wg = float(np.clip(res.x[1], 0.0, 0.6)) if fit_wing else 0.0
    return fw, wg, float(res.fun)


def power_beyond(fwhm_A, wing, radius_A, wing_scale=4.0):
    """Fraction of the two-Gaussian kernel's area outside +/-radius.

    This, not the wing's WIDTH, is what the data constrain. Fitting the width
    freely sends three of eight stars to the bound and improves the rms by 2%:
    a small fraction in a very broad wing and a larger fraction in a moderate
    one look alike once the spectrum is continuum-normalised, because a broad
    enough wing is just a pedestal. The power at +/-10 A is stable across both
    parameterisations (2.9% fixed-width vs 3.2% free); the power beyond +/-20 A
    is not (0.2% vs 1.5%) and should not be quoted.
    """
    s1 = fwhm_A / 2.3548
    s2 = s1 * wing_scale
    return ((1 - wing) * erfc(radius_A / (s1 * np.sqrt(2)))
            + wing * erfc(radius_A / (s2 * np.sqrt(2))))


def core_mask(w, lo, hi, half=CORE_HALF):
    """True at the high-order Balmer line cores inside [lo, hi]."""
    lines = [l for l in hydrogen_lines(lo, hi, series=(2,), nmax=NMAX)]
    m = np.zeros_like(w, bool)
    for l in lines:
        m |= np.abs(w - l) < half
    return m, lines


def main():
    grid = Grid()
    rows = [r for r in csv.DictReader(open(ROOT / 'data' / 'sample.csv'))
            if r['tier'] in ('primary', 'secondary')]
    lo, hi = TEST
    print(f'Test window {lo:.0f}-{hi:.0f} A (resolved H7-H12); control {CONTROL[0]:.0f}'
          f'-{CONTROL[1]:.0f} A; cores = +/-{CORE_HALF:.0f} A of a Balmer line\n')
    print(f'{"star":<10}{"FWHM_G":>8}{"wing":>6}{"FWHM_w":>8}'
          f'{"NGSL-XSL":>10}{"NGSL-mod":>10}{"mod Gauss":>10}{"mod+wing":>10}')
    print(f'{"":<10}{"[A]":>8}{"frac":>6}{"[A]":>8}'
          + ' '.join(f'{x:>9}' for x in ('core-cont',) * 4))
    print('-' * 80)
    out = []
    for r in rows:
        star = r['star']
        try:
            wx, fx, ex, _ = load_xsl(r['xslid'])
            ng = load_ngsl(star)
        except Exception as exc:
            print(f'{star:<10}  FAILED {type(exc).__name__}')
            continue
        wn, fn = np.asarray(ng.wavelength, float), np.asarray(ng.flux, float)
        okx = np.isfinite(fx) & (fx > 0)
        wx, fx = wx[okx], fx[okx]

        fwg, _, rg = fit_control(wx, fx, wn, fn, fit_wing=False)
        fww, wing, rw = fit_control(wx, fx, wn, fn, fit_wing=True)
        if not np.isfinite(fwg):
            print(f'{star:<10}  no control fit')
            continue

        sn = (wn > lo) & (wn < hi) & np.isfinite(fn) & (fn > 0)
        sx = (wx > lo - 60) & (wx < hi + 60)
        yn = cnorm(wn[sn], fn[sn])
        cm, _ = core_mask(wn[sn], lo, hi)

        cells = []
        for fw, wg in ((fwg, 0.0), (fww, wing)):
            yd = cnorm(wn[sn], degrade(wx[sx], fx[sx], wn[sn], fw, wg))
            res = 100 * (yn - yd) / yd
            cells.append((np.nanmedian(res[cm]) - np.nanmedian(res[~cm]),
                          np.nanmedian(res[~cm])))

        b = best_node(star)
        th = dict(teff=b['teff'], logg=b['logg'], mh=b['mh'], ebv=b['ebv'],
                  vsini=b['vsini'])
        mv = predict(th, [ng], grid)[0].value[sn]
        ym = cnorm(wn[sn], mv)
        resm = 100 * (yn - ym) / ym
        mk = np.nanmedian(resm[~cm])
        mc = np.nanmedian(resm[cm]) - mk

        # THE TEST. Re-convolve the MODEL with the winged profile measured
        # against XSL in the control window -- a region with no Balmer line in
        # it -- and ask whether the core excess survives. Nothing here is fitted
        # to the cores, so this is a prediction.
        w0 = grid.wave
        f0 = broaden_rot(w0, spectrum_at(grid, b['teff'], b['logg'], b['mh']),
                         b['vsini'])
        if b['ebv']:
            f0 = redden(w0, f0, b['ebv'], 3.1)
        sw = (w0 > lo - 80) & (w0 < hi + 80)
        mwing = []
        for fw, wg in ((fwg, 0.0), (fww, wing)):
            mm = rebin_to_pixels(w0[sw], kernel(w0[sw], f0[sw], fw, wg), wn[sn])
            ymm = cnorm(wn[sn], mm)
            rr = 100 * (yn - ymm) / ymm
            mwing.append(np.nanmedian(rr[cm]) - np.nanmedian(rr[~cm]))

        print(f'{star:<10}{fwg:>8.2f}{wing:>6.2f}{fww:>8.2f}'
              f'{cells[0][0]:>10.2f}{mc:>10.2f}'
              f'{mwing[0]:>10.2f}{mwing[1]:>10.2f}'
              + ('   (below grid)' if below_grid(star) else ''))
        out.append(dict(star=star, fwhm_gauss=fwg, wing=wing, fwhm_wing=fww,
                        rms_gauss=rg, rms_wing=rw,
                        xsl_core=cells[0][0], xsl_cont=cells[0][1],
                        model_core=mc, model_cont=mk,
                        xslwing_core=cells[1][0],
                        model_gauss_core=mwing[0], model_wing_core=mwing[1],
                        below_grid=below_grid(star)))

    if not out:
        return
    print('-' * 80)
    ing = [o for o in out if not o['below_grid']]
    g = ing if ing else out
    mg = np.median([o['model_gauss_core'] for o in g])
    mw = np.median([o['model_wing_core'] for o in g])
    xg = np.median([o['xsl_core'] for o in g])
    print(f'\n  Inside the grid (n={len(g)}), median core-minus-continuum:\n')
    print(f'    NGSL vs the MODEL, Gaussian R=600 profile   {mg:+6.2f}%   '
          f'<- what was called NLTE')
    print(f'    NGSL vs XSL degraded, same Gaussian         {xg:+6.2f}%   '
          f'<- no model involved')
    print(f'    NGSL vs the MODEL, Gaussian + {np.median([o["wing"] for o in g]):.2f} wing  '
          f'      {mw:+6.2f}%   <- profile measured on 4200-4600 A')
    pw = [power_beyond(o['fwhm_wing'], o['wing'], 10.0) for o in g]
    p5 = [power_beyond(o['fwhm_wing'], o['wing'], 5.0) for o in g]
    gauss10 = power_beyond(float(np.median([o['fwhm_wing'] for o in g])), 0.0, 10.0)
    print(f'\n    Two-Gaussian shape: core FWHM '
          f'{np.median([o["fwhm_wing"] for o in g]):.2f} A, '
          f'{100 * np.median(p5):.0f}% of the power beyond +/-5 A, '
          f'{100 * np.median(pw):.1f}% beyond +/-10 A')
    print(f'    (a pure Gaussian of that core width: {100 * gauss10:.4f}% '
          f'beyond +/-10 A)')
    print('\n    WHAT THIS DOES AND DOES NOT SHOW.')
    print('      DOES: the NGSL cores sit only '
          f'{xg:+.2f}% away from XSL degraded to the same resolution, with no')
    print('        model involved. Whatever fills them is instrumental, not NLTE.')
    print('      DOES NOT: identify the profile. The core excess is degenerate')
    print('        with effective WIDTH -- a plain Gaussian spans +16.8% to -1.2%')
    print('        across plausible widths -- so it is evidence for no particular')
    print('        shape. An earlier version of this script claimed a winged')
    print('        profile "removes 86%"; that was a width result, and is')
    print('        retracted. The adopted profile was chosen on the control-window')
    print('        rms and the published STIS core (explore/ngsl_lsf_shape.py).')
    p = ROOT / 'data' / 'ngsl_core_excess.csv'
    with open(p, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0]))
        w.writeheader()
        w.writerows(out)
    print(f'  -> {p.relative_to(ROOT)}')


if __name__ == '__main__':
    main()
