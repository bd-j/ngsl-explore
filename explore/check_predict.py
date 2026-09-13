"""Validate observations.py + predict.py end to end on one star.

Predicts the same star two ways -- at its nominal catalog parameters
(interpolated) and at the nearest grid node (the node-exact fast path the 1705
node scan will use) -- projects both onto NGSL, XSL and Gaia photometry, and
plots prediction against observation for each.

This is the smoke test before any long scan: one work unit, end to end, with the
artifact checked rather than the exit code.

    python3 explore/check_predict.py --star HD194453
"""
import argparse
import csv
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from fitting.model import Grid
from fitting.observations import (load_ngsl, load_xsl, load_photometry,
                                  load_xp, BREAK_WINDOW)
from fitting.predict import predict, spectrum_at, _on_node, NODE_ATOL
from fitting.calibration import solve, residual, chi2
from common.photometry import overlaps, coverage

ROOT = Path(__file__).resolve().parent.parent
BALMER, PASCHEN = 3646.0, 8205.9
DUST_WINDOW = (3200.0, 3400.0)       # line-free, NGSL-only dust constraint

OBS_C, MOD_C, MOD2_C = '#2a78d6', '#eb6834', '#7a3fa8'
SURFACE, INK, MUTED, GRIDC = '#fcfcfb', '#22262b', '#6b7280', '#dfe3e8'


def style(ax):
    ax.set_facecolor(SURFACE)
    ax.grid(alpha=.25, color=GRIDC, lw=.7)
    ax.tick_params(labelsize=8, colors=MUTED)
    for s in ax.spines.values():
        s.set_color(GRIDC)


def nearest_node(grid, teff, logg, mh):
    return (float(grid.teff[np.argmin(np.abs(grid.teff - teff))]),
            float(grid.logg[np.argmin(np.abs(grid.logg - logg))]),
            float(grid.mh[np.argmin(np.abs(grid.mh - mh))]))


def sample_row(star):
    for r in csv.DictReader(open(ROOT / 'data' / 'sample.csv')):
        if r['star'] == star:
            return r
    raise KeyError(star)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--star', default='HD194453')
    ap.add_argument('--vsini', type=float, default=0.0)
    ap.add_argument('--ebv', type=float, default=0.0)
    a = ap.parse_args()

    row = sample_row(a.star)
    grid = Grid()
    print(f'{a.star}: grid {len(grid.teff)}x{len(grid.logg)}x{len(grid.mh)}, '
          f'{grid.wave[0]:.0f}-{grid.wave[-1]:.0f} A at R={grid.resolution:.0f}')

    obs = []
    for fn, label in ((load_ngsl, 'ngsl'), (load_xsl, 'xsl'),
                      (load_xp, 'xp'), (load_photometry, 'phot')):
        try:
            o = fn(a.star)
            obs.append(o)
            print(f'  {o!r}')
        except Exception as exc:
            print(f'  no {label}: {type(exc).__name__}: {exc}')

    xp = next((o for o in obs if o.name == 'xp'), None)
    if xp is not None:
        print(f'\n  Gaia XP bands (held-out window {BREAK_WINDOW[0]:.0f}-'
              f'{BREAK_WINDOW[1]:.0f} A excluded):')
        print(f'    kept    : {", ".join(xp.meta["names"])}')
        print(f'    dropped : {", ".join(xp.meta["excluded"]) or "none"}')

    phot = next((o for o in obs if o.name == 'phot'), None)
    if phot is not None:
        print('\n  filter coverage against the model grid '
              f'({grid.wave[0]:.0f}-{grid.wave[-1]:.0f} A):')
        for nm, cov, lo, hi in coverage(grid.wave, phot.filters):
            print(f'    {nm:<10} {lo:>6.0f}-{hi:<6.0f} '
                  f'{"covered" if cov else "NOT COVERED -> NaN"}')
        ov = overlaps(phot.filters, BREAK_WINDOW)
        print(f'    filters overlapping the held-out break '
              f'{BREAK_WINDOW}: {", ".join(ov) if ov else "none"}')

    # nominal (interpolated) and nearest node (fast path)
    t0 = (float(row['teff_ngsl']), float(row['logg_ngsl']), float(row['mh_ngsl']))
    tn = nearest_node(grid, *t0)
    thetas = [('nominal (interp)', dict(teff=t0[0], logg=t0[1], mh=t0[2],
                                        ebv=a.ebv, vsini=a.vsini)),
              ('nearest node', dict(teff=tn[0], logg=tn[1], mh=tn[2],
                                    ebv=a.ebv, vsini=a.vsini))]
    print(f'\n  nominal  Teff={t0[0]:.0f} logg={t0[1]:.2f} [M/H]={t0[2]:+.2f}'
          f'   -> on-node? '
          + str(all(_on_node(getattr(grid, k), v, NODE_ATOL[k]) is not None
                    for k, v in zip(('teff', 'logg', 'mh'), t0))))
    print(f'  node     Teff={tn[0]:.0f} logg={tn[1]:.2f} [M/H]={tn[2]:+.2f}'
          f'   -> on-node? '
          + str(all(_on_node(getattr(grid, k), v, NODE_ATOL[k]) is not None
                    for k, v in zip(('teff', 'logg', 'mh'), tn))))

    results = {}
    for label, th in thetas:
        preds = predict(th, obs, grid)
        results[label] = {}
        print(f'\n  {label}:')
        for o, p in zip(obs, preds):
            from fitting.calibration import usable
            u = usable(o, p.value)
            if u.sum() < o.ndata:
                dropped = ([n for n, ok in zip(o.meta.get('names', []), u) if not ok]
                           if o.filters is not None else None)
                print(f'    {o.name:<5} {o.ndata - int(u.sum())} of {o.ndata} points '
                      f'dropped: model not finite'
                      + (f' ({", ".join(dropped)} run past the grid)' if dropped else ''))
            cal, c = solve(o, p.value)
            r = residual(o, cal)
            x2 = chi2(o, cal)
            results[label][o.name] = (cal, r, c, u)
            rms = np.nanstd(r[u]) * 100
            print(f'    {o.name:<5} chi2/N = {x2 / max(int(u.sum()), 1):8.2f}   '
                  f'residual rms = {rms:6.2f}%   ncoeff={len(c)}   n={int(u.sum())}')
            if o.name == 'ngsl':
                for wlo, whi, nm in ((DUST_WINDOW[0], DUST_WINDOW[1], 'dust 3200-3400'),
                                     (BREAK_WINDOW[0], BREAK_WINDOW[1], 'break (held out)'),
                                     (4000., 9400., 'red 4000-9400')):
                    s = o.mask & (o.wavelength > wlo) & (o.wavelength < whi)
                    if s.sum() > 5:
                        print(f'            {nm:<16} rms {np.nanstd(r[s])*100:5.2f}%  '
                              f'median {np.nanmedian(r[s])*100:+6.2f}%')

    figure(a.star, obs, results, row)


def figure(star, obs, results, row):
    ng = next(o for o in obs if o.name == 'ngsl')
    xs = next((o for o in obs if o.name == 'xsl'), None)
    ph = next((o for o in obs if o.name == 'phot'), None)

    fig = plt.figure(figsize=(12.5, 13.5))
    gs = fig.add_gridspec(5, 2, height_ratios=[2.0, 1.0, 1.6, 1.6, 1.4],
                          hspace=.42, wspace=.22)
    fig.patch.set_facecolor(SURFACE)
    labels = list(results)

    # --- NGSL full range + residual --------------------------------------
    ax, rax = fig.add_subplot(gs[0, :]), fig.add_subplot(gs[1, :])
    for a_ in (ax, rax):
        style(a_)
        for lam in (BALMER, PASCHEN):
            a_.axvline(lam, color=MUTED, ls='--', lw=1)
        a_.axvspan(*BREAK_WINDOW, color='#c0392b', alpha=.09, lw=0)
        a_.axvspan(*DUST_WINDOW, color='#2a78d6', alpha=.09, lw=0)
    m = ng.mask
    ax.plot(ng.wavelength[m], ng.flux[m], color=OBS_C, lw=1.3, label='NGSL observed')
    for lab, c in zip(labels, (MOD_C, MOD2_C)):
        cal = results[lab]['ngsl'][0]
        ax.plot(ng.wavelength[m], cal[m], color=c, lw=1.1, label=f'model, {lab}')
    ax.set_yscale('log')
    ax.set_ylabel(r'F$_\lambda$ [erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$]', fontsize=9,
                  color=INK)
    ax.set_title(f'{star} — NGSL, one free scalar (blue band: dust window; '
                 f'red band: held-out break)', fontsize=10, color=INK)
    ax.legend(fontsize=8, loc='upper right', framealpha=.92)
    for lab, c in zip(labels, (MOD_C, MOD2_C)):
        rax.plot(ng.wavelength[m], results[lab]['ngsl'][1][m] * 100, color=c, lw=.9,
                 label=lab)
    rax.axhline(0, color=MUTED, lw=1)
    rax.set_ylim(-12, 12)
    rax.set_ylabel('(obs−model)/model [%]', fontsize=9, color=INK)
    rax.set_xlabel(r'Wavelength [$\AA$, vacuum]', fontsize=9, color=INK)
    rax.legend(fontsize=8, loc='upper right', ncol=2, framealpha=.92)

    # --- NGSL break zoom, and the dust window ----------------------------
    for col, (lo, hi, ttl) in enumerate((
            (3300., 3900., 'Balmer break region (the prediction)'),
            (7900., 8600., 'Paschen break region'))):
        axz = fig.add_subplot(gs[2, col])
        style(axz)
        axz.axvspan(*BREAK_WINDOW, color='#c0392b', alpha=.09, lw=0)
        axz.axvspan(*DUST_WINDOW, color='#2a78d6', alpha=.09, lw=0)
        s = m & (ng.wavelength > lo) & (ng.wavelength < hi)
        axz.plot(ng.wavelength[s], ng.flux[s], color=OBS_C, lw=1.4)
        for lab, c in zip(labels, (MOD_C, MOD2_C)):
            axz.plot(ng.wavelength[s], results[lab]['ngsl'][0][s], color=c, lw=1.1)
        for lam in (BALMER, PASCHEN):
            if lo < lam < hi:
                axz.axvline(lam, color=MUTED, ls='--', lw=1)
        axz.set_xlim(lo, hi)
        axz.set_title(ttl, fontsize=9, color=INK)
        axz.set_xlabel(r'$\lambda$ [$\AA$]', fontsize=8, color=INK)
        axz.set_ylabel(r'F$_\lambda$', fontsize=8, color=INK)

    # --- XSL, continuum marginalised -------------------------------------
    if xs is not None:
        for col, (lo, hi, ttl) in enumerate((
                (4290., 4400., r'XSL: H$\gamma$ 4341 (continuum marginalised)'),
                (5150., 5220., 'XSL: Mg I b 5167-5183'))):
            axx = fig.add_subplot(gs[3, col])
            style(axx)
            s = xs.mask & (xs.wavelength > lo) & (xs.wavelength < hi)
            axx.plot(xs.wavelength[s], xs.flux[s], color=OBS_C, lw=1.1,
                     label='XSL observed')
            for lab, c in zip(labels, (MOD_C, MOD2_C)):
                axx.plot(xs.wavelength[s], results[lab]['xsl'][0][s], color=c,
                         lw=1.0, label=f'model, {lab}')
            axx.set_xlim(lo, hi)
            axx.set_title(ttl, fontsize=9, color=INK)
            axx.set_xlabel(r'$\lambda$ [$\AA$]', fontsize=8, color=INK)
            axx.set_ylabel(r'F$_\lambda$', fontsize=8, color=INK)
            if col == 0:
                axx.legend(fontsize=7.5, loc='lower left', framealpha=.92)

    # --- Gaia XP bands: the dust lever ------------------------------------
    xp = next((o for o in obs if o.name == 'xp'), None)
    if xp is not None:
        axp = fig.add_subplot(gs[4, 0])
        axp2 = fig.add_subplot(gs[4, 1])
        for a_ in (axp, axp2):
            style(a_)
        lam = np.array([f.wave_effective for f in xp.filters])
        u0 = results[labels[0]]['xp'][3]
        axp.errorbar(lam, xp.flux, yerr=xp.uncertainty, fmt='o', color=OBS_C,
                     ms=6, capsize=3, label='Gaia XP, banded')
        for lab, c in zip(labels, (MOD_C, MOD2_C)):
            cal = results[lab]['xp'][0]
            axp.plot(lam[u0], cal[u0], 's', color=c, ms=6, label=f'model, {lab}')
        axp.set_yscale('log')
        axp.set_ylabel('band flux', fontsize=9, color=INK)
        axp.set_xlabel(r'effective $\lambda$ [$\AA$]', fontsize=9, color=INK)
        axp.set_title('Gaia XP in bands (LSF-free), one free scalar',
                      fontsize=9, color=INK)
        axp.legend(fontsize=7.5, framealpha=.92)

        for lab, c in zip(labels, (MOD_C, MOD2_C)):
            r = results[lab]['xp'][1]
            axp2.errorbar(lam[u0], r[u0] * 100,
                          yerr=100 * xp.uncertainty[u0] / xp.flux[u0],
                          fmt='o-', color=c, lw=1.1, ms=5, capsize=3, label=lab)
        axp2.axhline(0, color=MUTED, lw=1)
        axp2.axvspan(*BREAK_WINDOW, color='#c0392b', alpha=.09, lw=0)
        axp2.set_xlabel(r'effective $\lambda$ [$\AA$]', fontsize=9, color=INK)
        axp2.set_ylabel('(obs-model)/model [%]', fontsize=9, color=INK)
        axp2.set_title('XP residual - this is the dust lever (err bars = 1% floor)',
                       fontsize=9, color=INK)
        axp2.legend(fontsize=7.5, framealpha=.92)

    fig.suptitle(
        f'{star}   predict() check: nominal Teff={float(row["teff_ngsl"]):.0f} / '
        f'log g={row["logg_ngsl"]} / [M/H]={row["mh_ngsl"]}   '
        f'(XSL: {row["teff_xsl"]} / {row["logg_xsl"]} / {row["mh_xsl"]})',
        fontsize=11, color=INK)
    out = ROOT / 'figures' / f'predict_check_{star}.png'
    fig.savefig(out, dpi=170, facecolor=SURFACE, bbox_inches='tight')
    plt.close(fig)
    print(f'\n  -> {out.relative_to(ROOT)}')


if __name__ == '__main__':
    main()
