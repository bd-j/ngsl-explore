"""Measure the NGSL line spread function directly, using XSL as the reference.

No model is involved. XSL resolves the same star ~10x better than NGSL, so
convolving XSL with a free-width kernel and rebinning onto NGSL pixels asks a
purely instrumental question: what profile turns the true spectrum into what
NGSL recorded?

This matters because two earlier results disagreed. The STIS LSF tables give
3.85 A FWHM for G430L, but fitting an effective width against the ATLAS12
models preferred ~7.0-7.5 A. That fit was discounted at the time because the
Balmer cores carry a genuine flux excess the LTE models cannot produce, so a
broader kernel could be absorbing physics. Using XSL instead removes the models
from the question entirely.

Both spectra are continuum-normalized in each window before comparison, so the
grey flux-calibration offset between the libraries (-11% to +4%, star
dependent) cannot influence the width. A wavelength shift is fitted jointly
with the width, so a residual wavecal error cannot masquerade as broadening.

Writes data/ngsl_lsf_measured.csv
"""
import csv
import sys
from pathlib import Path

import numpy as np
from astropy.io import fits
from scipy.optimize import minimize

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common.lsf import degrade_to, rebin_to_pixels
from common.ngsl_wavecal import apply_wavecal, load_table
from common.xsl_load import load as load_xsl, resolving_power

ROOT = Path(__file__).resolve().parent.parent

# (label, lo, hi, tabulated FWHM in A). Windows are chosen line-rich: a
# featureless continuum constrains no width at all.
WINDOWS = [('G430L 3700-4100', 3700., 4100., 3.85),
           ('G430L 4200-4600', 4200., 4600., 3.85),
           ('G430L 4700-5100', 4700., 5100., 3.85),
           ('G750L 6400-6800', 6400., 6800., 8.09),
           ('G750L 8400-9000', 8400., 9000., 8.09)]

PAIRS = [('HD194453', 'X0196', 'h_stis_ngsl_hd194453_v2.fits'),
         ('HD143459', 'X0386', 'h_stis_ngsl_hd143459_v2.fits'),
         ('HD128801', 'X0243', 'h_stis_ngsl_hd128801_v2.fits')]


def cnorm(w, f, deg=3):
    """Continuum-normalize, so only line shape drives the fit."""
    ok = np.isfinite(f) & (f > 0)
    if ok.sum() < deg + 2:
        return f
    return f / np.polyval(np.polyfit(w[ok], f[ok], deg), w)


def match(wx, fx, wn, fn, lo, hi, fwhm_xsl):
    """-> (best FWHM in A, best shift in A, rms). Width and shift fitted jointly."""
    seg = (wn > lo) & (wn < hi) & np.isfinite(fn) & (fn > 0)
    if seg.sum() < 30:
        return np.nan, np.nan, np.nan
    wns, fns = wn[seg], cnorm(wn[seg], fn[seg])
    # Crop XSL to the window plus a margin wide enough for the widest kernel.
    # degrade_to resamples onto a 0.01 A grid; over the full 3500-24770 A range
    # that is 2.1M points per function evaluation and the fit never finishes.
    pad = 120.0
    kx = (wx > lo - pad) & (wx < hi + pad)
    wx, fx = wx[kx], fx[kx]
    if len(wx) < 100:
        return np.nan, np.nan, np.nan

    def resid(p):
        fw, dv = p
        if not (0.5 < fw < 25.0) or abs(dv) > 3.0:
            return 1e3
        d = degrade_to(wx, fx, fwhm_xsl, fw)
        r = rebin_to_pixels(wx + dv, d, wns)
        ok = np.isfinite(r) & (r > 0)
        if ok.sum() < 20:
            return 1e3
        return float(np.std(cnorm(wns[ok], r[ok]) - fns[ok]))

    best = None
    for fw0 in (3.0, 5.0, 7.0, 9.0):
        r = minimize(resid, [fw0, 0.0], method='Nelder-Mead',
                     options=dict(xatol=1e-2, fatol=1e-5, maxiter=400))
        if best is None or r.fun < best.fun:
            best = r
    return float(best.x[0]), float(best.x[1]), float(best.fun)


rows = []
print('NGSL LSF measured against XSL (no model involved)\n')
print(f'{"window":<18}{"tabulated":>10}' +
      ''.join(f'{s:>18}' for s, _, _ in PAIRS))
data = {}
for star, xid, nf in PAIRS:
    d = fits.getdata('data/spectra/' + nf)
    wn = apply_wavecal(d['WAVELENGTH'].astype(float), star, load_table())
    wx, fx, ex, h = load_xsl(xid)
    data[star] = (wn, d['FLUX'].astype(float), wx, fx)

for lbl, lo, hi, tab in WINDOWS:
    line = f'{lbl:<18}{tab:>8.2f} A'
    for star, _, _ in PAIRS:
        wn, fn, wx, fx = data[star]
        mid = 0.5 * (lo + hi)
        fw, dv, rms = match(wx, fx, wn, fn, lo, hi, mid / resolving_power(mid))
        rows.append(dict(window=lbl, star=star, lo=lo, hi=hi,
                         tabulated_A=tab,
                         measured_A=round(fw, 2) if np.isfinite(fw) else '',
                         shift_A=round(dv, 2) if np.isfinite(dv) else '',
                         ratio=round(fw / tab, 2) if np.isfinite(fw) else '',
                         rms=round(rms, 4) if np.isfinite(rms) else ''))
        line += (f'{fw:>10.2f} A{dv:>+7.2f}' if np.isfinite(fw)
                 else f'{"--":>18}')
    print(line)

with open(ROOT / 'data' / 'ngsl_lsf_measured.csv', 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0]))
    w.writeheader(); w.writerows(rows)

g430 = [r['measured_A'] for r in rows if r['window'].startswith('G430L') and r['measured_A'] != '']
g750 = [r['measured_A'] for r in rows if r['window'].startswith('G750L') and r['measured_A'] != '']
print(f'\nG430L: measured {np.median(g430):.2f} A  vs tabulated 3.85  '
      f'-> {np.median(g430)/3.85:.2f}x   (R = {3646/np.median(g430):.0f} at the break)')
print(f'G750L: measured {np.median(g750):.2f} A  vs tabulated 8.09  '
      f'-> {np.median(g750)/8.09:.2f}x')
print('\n-> data/ngsl_lsf_measured.csv')
