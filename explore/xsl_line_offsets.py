"""Velocity offset of each line in each star: XSL observed minus model.

The point is to decide how the offset should be treated, so the table is built
to separate the two possibilities:

  a per-STAR offset   residual radial velocity left by XSL's rest-frame
                      reduction. Would move every line in a star together.
  a per-LINE offset   a wavelength error in the Kurucz list, or a blend whose
                      component ratio the model gets wrong. Would repeat in
                      every star.

So the table carries both margins and a two-way variance decomposition. If the
star margin dominates, fit an RV; if the line margin dominates, the line list is
at fault and an RV would be fitting the wrong thing.

Measured by CROSS-CORRELATION over each panel's window, not by a single line
centroid. The centroid version left most of the table empty: a fast rotator
washes its lines below any sensible depth cut, and where two lines sit closer
than the 0.41 A resolution element the minimum is not attributable to either.
Cross-correlation uses every line in the window at once, which is also what the
eye is judging when a panel "looks shifted".

Both spectra are divided by a low-order polynomial first, so no flux calibration
is needed and a continuum error cannot masquerade as a shift.

Blended lines are flagged. They are not measuring wavelength at all -- when the
model gets the relative strengths of a close pair wrong, the composite minimum
moves, which looks identical to a shift.

Writes data/xsl_line_offsets.csv

    python3 explore/xsl_line_offsets.py
"""
import argparse
import csv
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from fitting.model import Grid
from fitting.observations import XSL_METAL_WINDOW, MODEL_BAD_REGIONS
from common.xsl_load import load as xsl_load
from common.lsf import broaden_R
from common.species import (read_lines, atmosphere_point, abundances,
                            ion_fraction, ELEMENT, ROMAN)

ROOT = Path(__file__).resolve().parent.parent
C_KMS = 2.99792458e5
MIN_DEPTH = 0.03          # skip anything shallower; a noise minimum is not a line
MIN_CORR = 0.60           # reject a window whose peak correlation is this poor
# A cross-correlation needs structure. Several feature windows are 1-5 A wide,
# which gave peaks anywhere in the search range (-70 to +68 km/s) purely from
# noise; the window is therefore widened about the panel centre to at least
# +/-CC_HALFWIDTH so the correlation has several lines to work with.
CC_HALFWIDTH = 10.0
ISOLATION = 0.35          # A; another comparable line nearer than this = blended
SEARCH = 0.45             # A; half-width of the centroid search


def xcorr_shift(x, y, xm, ym, lo, hi, maxshift=1.2, step=0.01, deg=3):
    """Shift in A that the MODEL must move to match the data, over [lo, hi].

    Positive = the data sits redward of the model.
    """
    s = (x > lo) & (x < hi) & np.isfinite(y) & (y > 0)
    for blo, bhi, _why in MODEL_BAD_REGIONS:     # never correlate on these
        s &= ~((x >= blo) & (x <= bhi))
    if s.sum() < 25:
        return np.nan, 0.0
    xs, ys = x[s], y[s]
    yn = ys / np.polyval(np.polyfit(xs, ys, deg), xs)
    mi = np.interp(xs, xm, ym)
    mn = mi / np.polyval(np.polyfit(xs, mi, deg), xs)
    if np.std(yn) < 1e-6 or np.std(mn) < 1e-6:
        return np.nan, 0.0
    shifts = np.arange(-maxshift, maxshift + step / 2, step)
    cc = np.array([np.corrcoef(yn, np.interp(xs, xs + d, mn))[0, 1]
                   for d in shifts])
    j = int(np.nanargmax(cc))
    if j in (0, len(shifts) - 1) or not np.isfinite(cc[j]):
        return np.nan, float(np.nanmax(cc))
    return float(shifts[j]), float(cc[j])


def centroid(x, y, lam, half=SEARCH):
    """Parabolic vertex through the three lowest points near lam.

    Validated by injection: a model shifted by +0.100 A is recovered at +0.104.
    An earlier version had the sign of the correction reversed, which reported
    every offset backwards.
    """
    s = (x > lam - half) & (x < lam + half)
    if s.sum() < 5:
        return np.nan, 0.0
    xs, ys = x[s], y[s]
    i = int(np.argmin(ys))
    if i == 0 or i == len(xs) - 1:
        return np.nan, 0.0
    y0, y1, y2 = ys[i - 1], ys[i], ys[i + 1]
    h = 0.5 * (xs[i + 1] - xs[i - 1])
    den = y0 - 2 * y1 + y2
    cen = xs[i] + (0.5 * h * (y0 - y2) / den if den != 0 else 0.0)
    depth = 1.0 - y1 / np.percentile(ys, 90)
    return cen, float(depth)


def panel_lines(T, ne, eps, nfeat=8):
    """One entry per PANEL of figures/metal_lines_*.png.

    Columns must correspond one-to-one with the panels, so the features are the
    top-N of data/metal_sensitivity.csv in the same order the figure uses, and
    the measured line is the deepest model line inside each feature window --
    which is the feature the panel is showing.
    """
    feats = sorted(csv.DictReader(open(ROOT / 'data' / 'metal_sensitivity.csv')),
                   key=lambda r: float(r['depth_change']))[:nfeat]
    out = []
    theta = 5040.0 / T
    for f in feats:
        lo, hi = float(f['lam_lo']), float(f['lam_hi'])
        rows = []
        for lam, gf, z, st, el in read_lines(lo, hi):
            if z not in ELEMENT or st > 2:
                continue
            fi = ion_fraction(z, st, T, ne)
            if fi <= 0:
                continue
            rows.append((gf + eps.get(z, -12.0) + np.log10(fi) - theta * el,
                         lam, f'{ELEMENT[z]} {ROMAN[st]}'))
        if not rows:
            continue
        rows.sort(reverse=True)
        s0, lam0, sp0 = rows[0]
        near = min((abs(l - lam0) for s, l, _ in rows[1:]
                    if s > s0 - 1.0 and abs(l - lam0) > 0.03), default=99.0)
        out.append(dict(lam=lam0, species=sp0, sep=near,
                        blended=near < ISOLATION,
                        # An explicit fit window wins: widening the 4134 panel
                        # to +/-10 A pushed it straight into the 4120-4127 and
                        # 4137.5-4143 ranges that MODEL_BAD_REGIONS excludes
                        # precisely because the models get them wrong.
                        **(dict(zip(('lo', 'hi'),
                                    XSL_METAL_WINDOW[float(f['lam_center'])]))
                           if float(f['lam_center']) in XSL_METAL_WINDOW else
                           dict(lo=min(lo, float(f['lam_center']) - CC_HALFWIDTH),
                                hi=max(hi, float(f['lam_center']) + CC_HALFWIDTH))),
                        panel=float(f['lam_center']),
                        dd=float(f['depth_change'])))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--ebv', type=float, default=0.0)
    ap.add_argument('--n', type=int, default=8, help='number of panels')
    a = ap.parse_args()

    grid = Grid()
    atm = ROOT / 'models' / 'work' / 'HD194453.atm'
    T, ne = atmosphere_point(atm)
    eps = abundances(atm)
    lines = panel_lines(T, ne, eps, a.n)
    print(f'{len(lines)} panels of figures/metal_lines_*.png '
          f'({sum(1 for l in lines if l["blended"])} blended)\n')

    stars = [r for r in csv.DictReader(open(ROOT / 'data' / 'sample.csv'))
             if r['tier'] != 'rejected']
    table, kept = {}, []
    for r in stars:
        star = r['star']
        try:
            w, f, e, _ = xsl_load(r['xslid'])
        except Exception as exc:
            print(f'  {star}: no XSL ({type(exc).__name__})')
            continue
        t = float(grid.teff[np.argmin(np.abs(grid.teff - float(r['teff_ngsl'])))])
        lg = float(grid.logg[np.argmin(np.abs(grid.logg - float(r['logg_ngsl'])))])
        z = float(grid.mh[np.argmin(np.abs(grid.mh - float(r['mh_ngsl'])))])
        mod = broaden_R(grid.wave, grid.interp(t, lg, z), 9800.)
        row = {}
        for L in lines:
            d, r = xcorr_shift(w, f, grid.wave, mod, L['lo'], L['hi'])
            if np.isfinite(d) and r > MIN_CORR:
                row[L['lam']] = C_KMS * d / L['panel']
        table[star] = row
        kept.append(star)

    print('  ' + f'{"panel":<11}' + ''.join(f'{L["panel"]:>8.1f}' for L in lines))
    print('  ' + f'{"line meas.":<11}' + ''.join(f'{L["lam"]:>8.1f}' for L in lines))
    print('  ' + f'{"species":<11}'
          + ''.join(f'{L["species"].replace(" ", ""):>8}' for L in lines))
    print('  ' + f'{"blended?":<11}'
          + ''.join(f'{("BLEND" if L["blended"] else "-"):>8}' for L in lines))
    print('  ' + '-' * (11 + 8 * len(lines)))
    for s in kept:
        cells = ''.join(f'{table[s][L["lam"]]:>8.1f}' if L['lam'] in table[s]
                        else f'{"--":>8}' for L in lines)
        iso = [table[s][L['lam']] for L in lines
               if L['lam'] in table[s] and not L['blended']]
        print('  ' + f'{s:<11}' + cells
              + (f'   | star mean (isolated) {np.mean(iso):+6.1f}' if iso else ''))
    print('  ' + '-' * (11 + 8 * len(lines)))
    means, sds = [], []
    for L in lines:
        v = [table[s][L['lam']] for s in kept if L['lam'] in table[s]]
        means.append(np.mean(v) if v else np.nan)
        sds.append(np.std(v) if len(v) > 1 else np.nan)
    print('  ' + f'{"LINE MEAN":<11}' + ''.join(f'{m:>8.1f}' for m in means))
    print('  ' + f'{"LINE SD":<11}' + ''.join(f'{d:>8.1f}' for d in sds))

    # two-way decomposition on the isolated lines only
    iso_lines = [L for L in lines if not L['blended']]
    M = np.array([[table[s].get(L['lam'], np.nan) for L in iso_lines]
                  for s in kept])
    ok = ~np.isnan(M)
    if ok.sum() > 6:
        grand = np.nanmean(M)
        star_eff = np.nanmean(M, axis=1) - grand
        line_eff = np.nanmean(M, axis=0) - grand
        resid = M - grand - star_eff[:, None] - line_eff[None, :]
        print(f'\n  Two-way decomposition, ISOLATED lines only '
              f'({M.shape[0]} stars x {M.shape[1]} lines):')
        print(f'    grand mean          {grand:+6.2f} km/s')
        print(f'    star-to-star sd     {np.nanstd(star_eff):6.2f} km/s   '
              f'<- a residual RV would live here')
        print(f'    line-to-line sd     {np.nanstd(line_eff):6.2f} km/s   '
              f'<- a line-list error would live here')
        print(f'    unexplained sd      {np.nanstd(resid):6.2f} km/s')
        print(f'\n    XSL pixel ~10 km/s, resolution FWHM ~31 km/s (UVB)')

    out = ROOT / 'data' / 'xsl_line_offsets.csv'
    with open(out, 'w', newline='') as fh:
        wr = csv.writer(fh)
        cols = ['star', 'panel_center', 'lam_vac', 'species', 'blended',
                'sep_A', 'offset_kms']
        wr.writerow(cols)
        for s in kept:
            for L in lines:
                if L['lam'] in table[s]:
                    rec = dict(star=s, panel_center=f'{L["panel"]:.1f}',
                               lam_vac=f'{L["lam"]:.3f}', species=L['species'],
                               blended='yes' if L['blended'] else 'no',
                               sep_A=f'{L["sep"]:.2f}',
                               offset_kms=f'{table[s][L["lam"]]:.2f}')
                    wr.writerow([rec[c] for c in cols])
    print(f'\n  -> {out.relative_to(ROOT)}')


if __name__ == '__main__':
    main()
