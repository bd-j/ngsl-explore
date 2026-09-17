"""SUPERSEDED -- kept for provenance, not for use. See explore/superseded/README.md

Which line-spread PROFILE does NGSL actually have?

`common.lsf` applies a single Gaussian at R = 600. That is known to be wrong in
shape, not just width: the observed high-order Balmer cores are filled relative
to a Gaussian-convolved model by ~2.3%, and degrading XSL to NGSL's resolution
reproduces most of it with no model involved (docs/CAVEATS.md,
explore/superseded/ngsl_core_excess.py). A single Gaussian of the same core width puts
0.0006% of its power beyond +/-10 A; the data want ~3%.

This fits a family of profiles, so the choice is made on evidence rather than on
whichever one was tried first. XSL is the reference throughout -- no model is
involved in the fit.

  gauss             the current model. One parameter.
  gauss_tophat      Gaussian convolved with a boxcar. The physically motivated
                    one: NGSL v2 spectra are co-adds of two DITHERED exposures
                    resampled onto a common grid, and both the dither offset and
                    the pixel are boxes. Predicts a flat-topped core, and wings
                    that still fall off as a Gaussian.
  gauss_gauss       Gaussian plus a broader Gaussian. No mechanism behind it;
                    included because it was tried first.
  gauss_lorentz     Gaussian plus a Lorentzian. The classic scattered-light /
                    grating-halo shape, and the only member here with genuinely
                    heavy tails.
  moffat            Moffat. One shape parameter controls the tail weight.

Fitted on a CONTROL window (4200-4600 A) with no Balmer line in it, then scored
on a held-out TEST window (3700-4000 A) by how much Balmer core excess each one
leaves behind. Fitting and scoring on different windows is the whole point: any
profile with enough freedom can flatten the residual it was fitted to.

    python3 explore/superseded/ngsl_lsf_shape.py
"""
import csv
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
from scipy.signal import fftconvolve

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))
from fitting.model import Grid
from fitting.observations import load_ngsl
from fitting.predict import spectrum_at
from fitting.scan import best_node
from common.extinction_ccm import redden
from common.lsf import broaden_rot, rebin_to_pixels
from common.lines import hydrogen_lines
from common.xsl_load import load as load_xsl
from common.figpath import below_grid

ROOT = Path(__file__).resolve().parent.parent.parent
CONTROL = (4200., 4600.)
TEST = (3700., 4000.)
STEP = 0.05                       # A, the grid kernels are built on
HALF = 60.0                       # A, kernel half-width
CORE_HALF = 3.0
NMAX = 12

# name -> (n_params, x0, bounds-ish scaling). Widths are FWHM in Angstroms.
FAMILY = {
    'gauss':         (1, [6.5]),
    'gauss_tophat':  (2, [4.5, 4.0]),
    'gauss_gauss':   (3, [5.2, 0.12, 4.0]),
    'gauss_lorentz': (3, [5.0, 0.10, 6.0]),
    'moffat':        (2, [6.0, 2.5]),
}


def profile(kind, p, dl=None):
    """-> normalised kernel on a uniform grid of spacing STEP."""
    if dl is None:
        dl = np.arange(-HALF, HALF + STEP / 2, STEP)
    if kind == 'gauss':
        s = abs(p[0]) / 2.3548
        k = np.exp(-0.5 * (dl / s) ** 2)
    elif kind == 'gauss_tophat':
        s = abs(p[0]) / 2.3548
        g = np.exp(-0.5 * (dl / s) ** 2)
        box = (np.abs(dl) <= abs(p[1]) / 2.0).astype(float)
        if box.sum() == 0:
            box[len(dl) // 2] = 1.0
        k = fftconvolve(g, box / box.sum(), mode='same')
    elif kind == 'gauss_gauss':
        s = abs(p[0]) / 2.3548
        f = float(np.clip(p[1], 0.0, 0.8))
        s2 = s * max(abs(p[2]), 1.2)
        k = ((1 - f) * np.exp(-0.5 * (dl / s) ** 2) / s
             + f * np.exp(-0.5 * (dl / s2) ** 2) / s2)
    elif kind == 'gauss_lorentz':
        s = abs(p[0]) / 2.3548
        f = float(np.clip(p[1], 0.0, 0.8))
        gam = max(abs(p[2]), 0.3)
        k = ((1 - f) * np.exp(-0.5 * (dl / s) ** 2) / (s * np.sqrt(2 * np.pi))
             + f * (gam / np.pi) / (dl ** 2 + gam ** 2))
    elif kind == 'moffat':
        beta = float(np.clip(p[1], 1.05, 12.0))
        a = abs(p[0]) / (2.0 * np.sqrt(2 ** (1 / beta) - 1))
        k = (1 + (dl / a) ** 2) ** (-beta)
    else:
        raise ValueError(kind)
    k = np.clip(k, 0, None)
    tot = k.sum()
    return k / tot if tot > 0 else k


def apply_profile(w, f, kind, p):
    wl = np.arange(w[0], w[-1], STEP)
    fi = np.interp(wl, w, f)
    sm = fftconvolve(fi, profile(kind, p), mode='same')
    return np.interp(w, wl, sm)


def cnorm(w, f, deg=3):
    ok = np.isfinite(f) & (f > 0)
    if ok.sum() < deg + 2:
        return np.full_like(np.asarray(f, float), np.nan)
    return f / np.polyval(np.polyfit(w[ok], f[ok], deg), w)


def fit_one(kind, wx, fx, wn, fn, lo, hi):
    sn = (wn > lo) & (wn < hi) & np.isfinite(fn) & (fn > 0)
    sx = (wx > lo - 2 * HALF) & (wx < hi + 2 * HALF) & np.isfinite(fx) & (fx > 0)
    if sn.sum() < 20 or sx.sum() < 200:
        return None, np.nan
    yn = cnorm(wn[sn], fn[sn])

    def rms(p):
        d = rebin_to_pixels(wx[sx], apply_profile(wx[sx], fx[sx], kind, p), wn[sn])
        return float(np.sqrt(np.nanmean((yn - cnorm(wn[sn], d)) ** 2)))

    r = minimize(rms, FAMILY[kind][1], method='Nelder-Mead',
                 options=dict(xatol=1e-3, fatol=1e-8, maxiter=2000))
    return r.x, float(r.fun)


def core_excess(kind, p, star, grid, ng, b):
    """Balmer core-minus-continuum left by this profile, on the TEST window."""
    lo, hi = TEST
    wn = np.asarray(ng.wavelength, float)
    sn = (wn > lo) & (wn < hi) & np.asarray(ng.mask, bool)
    yn = cnorm(wn[sn], np.asarray(ng.flux, float)[sn])
    core = np.zeros(sn.sum(), bool)
    for l in hydrogen_lines(lo, hi, series=(2,), nmax=NMAX):
        core |= np.abs(wn[sn] - l) < CORE_HALF
    w0 = grid.wave
    f0 = broaden_rot(w0, spectrum_at(grid, b['teff'], b['logg'], b['mh']),
                     b['vsini'])
    if b['ebv']:
        f0 = redden(w0, f0, b['ebv'], 3.1)
    sx = (w0 > lo - 2 * HALF) & (w0 < hi + 2 * HALF)
    m = rebin_to_pixels(w0[sx], apply_profile(w0[sx], f0[sx], kind, p), wn[sn])
    ym = cnorm(wn[sn], m)
    r = 100 * (yn - ym) / ym
    return float(np.nanmedian(r[core]) - np.nanmedian(r[~core]))


def main():
    grid = Grid()
    rows = [r for r in csv.DictReader(open(ROOT / 'data' / 'sample.csv'))
            if r['tier'] in ('primary', 'secondary') and not below_grid(r['star'])]
    kinds = list(FAMILY)
    fits = {k: [] for k in kinds}
    print(f'Fitted on {CONTROL[0]:.0f}-{CONTROL[1]:.0f} A (no Balmer line), '
          f'scored on {TEST[0]:.0f}-{TEST[1]:.0f} A\n')
    print(f'{"star":<10}' + ''.join(f'{k:>16}' for k in kinds))
    print(f'{"":<10}' + ''.join(f'{"rms   core%":>16}' for _ in kinds))
    print('-' * (10 + 16 * len(kinds)))
    out = []
    for r in rows:
        star = r['star']
        try:
            wx, fx, _, _ = load_xsl(r['xslid'])
            ng = load_ngsl(star)
        except Exception:
            continue
        wn, fn = np.asarray(ng.wavelength, float), np.asarray(ng.flux, float)
        ok = np.isfinite(fx) & (fx > 0)
        wx, fx = wx[ok], fx[ok]
        b = best_node(star)
        line = f'{star:<10}'
        rec = dict(star=star)
        for k in kinds:
            p, rms = fit_one(k, wx, fx, wn, fn, *CONTROL)
            if p is None:
                line += f'{"--":>16}'
                continue
            ce = core_excess(k, p, star, grid, ng, b)
            fits[k].append((rms, ce, p))
            rec[f'{k}_rms'] = rms
            rec[f'{k}_core'] = ce
            rec[f'{k}_par'] = ' '.join(f'{v:.3f}' for v in p)
            line += f'{rms:>9.4f}{ce:>7.2f}'
        print(line)
        out.append(rec)
    print('-' * (10 + 16 * len(kinds)))
    print(f'{"median":<10}' + ''.join(
        f'{np.median([x[0] for x in fits[k]]):>9.4f}'
        f'{np.median([x[1] for x in fits[k]]):>7.2f}' if fits[k] else f'{"--":>16}'
        for k in kinds))
    print(f'{"n params":<10}' + ''.join(f'{FAMILY[k][0]:>16d}' for k in kinds))
    print('\n  Lower rms = better description of the CONTROL window.')
    print('  |core%| near zero = the Balmer cores it was NOT fitted to come out right.')
    best = min((k for k in kinds if fits[k]),
               key=lambda k: abs(np.median([x[1] for x in fits[k]])))
    print(f'\n  Smallest leftover core excess: {best} '
          f'({np.median([x[1] for x in fits[best]]):+.2f}%), '
          f'median parameters '
          + ' '.join(f'{v:.3f}' for v in np.median([x[2] for x in fits[best]], axis=0)))
    p = ROOT / 'data' / 'superseded' / 'ngsl_lsf_shape.csv'
    with open(p, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=sorted({k for r in out for k in r}))
        w.writeheader()
        w.writerows(out)
    print(f'  -> {p.relative_to(ROOT)}')


if __name__ == '__main__':
    main()
