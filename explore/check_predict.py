"""Validate the conditioning / held-out split end to end on one star.

The design this checks:

  CONDITION ON   NGSL collapsed into synthetic bands (continuum -> dust and
                 continuum shape) + XSL with its continuum marginalised away
                 (line profiles -> Teff, log g, v sin i, immune to reddening).
  PREDICT        the NGSL spectrum inside the Balmer break (3550-4000 A) and
                 inside the Paschen region (8180-9500 A), using the scalar
                 solved on the BANDS. Neither window is ever fitted.

Every NGSL pixel is therefore used at most once, and both breaks are
predictions rather than fits. The hydrogen lines are used from neither -- XSL
resolves them ~16x better for these same stars.

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
from fitting.observations import (load_ngsl, conditioning_set, heldout,
                                  ngsl_band_edges, BREAK_WINDOW, PASCHEN_WINDOW)
from fitting.predict import predict
from fitting.calibration import solve, residual, chi2, usable

ROOT = Path(__file__).resolve().parent.parent
BALMER, PASCHEN = 3646.0, 8205.9
OBS_C, MOD_C, MOD2_C = '#2a78d6', '#eb6834', '#7a3fa8'
SURFACE, INK, MUTED, GRIDC = '#fcfcfb', '#22262b', '#6b7280', '#dfe3e8'
BAND_C, HELD_C = '#2a78d6', '#c0392b'


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

    row, grid = sample_row(a.star), Grid()
    bands = ngsl_band_edges()
    print(f'{a.star}: {len(bands)} NGSL bands, '
          f'{bands[0][1]:.0f}-{bands[-1][2]:.0f} A')
    print(f'  held out: Balmer {BREAK_WINDOW}, Paschen {PASCHEN_WINDOW}')

    cond = conditioning_set(a.star)
    nb = cond[0]
    xs = next((o for o in cond if o.name == 'xsl'), None)
    spec = load_ngsl(a.star)
    held = {'Balmer': heldout(a.star, window=BREAK_WINDOW),
            'Paschen': heldout(a.star, window=PASCHEN_WINDOW)}
    print('  conditioning: ' + ', '.join(f'{o.name}(n={o.ndata})' for o in cond))
    print('  held out    : ' + ', '.join(f'{k}(n={v.ndata})'
                                         for k, v in held.items()))

    t0 = (float(row['teff_ngsl']), float(row['logg_ngsl']), float(row['mh_ngsl']))
    tn = nearest_node(grid, *t0)
    thetas = [('nominal (interp)', dict(teff=t0[0], logg=t0[1], mh=t0[2],
                                        ebv=a.ebv, vsini=a.vsini)),
              ('nearest node', dict(teff=tn[0], logg=tn[1], mh=tn[2],
                                    ebv=a.ebv, vsini=a.vsini))]
    print(f'\n  nominal Teff={t0[0]:.0f} logg={t0[1]:.2f} [M/H]={t0[2]:+.2f}  '
          f'-> node Teff={tn[0]:.0f} logg={tn[1]:.2f} [M/H]={tn[2]:+.2f}')

    results, xsl_cal = {}, {}
    for label, th in thetas:
        print(f'\n  {label}:')
        pb = predict(th, [nb], grid)[0]
        cal_b, c = solve(nb, pb.value)
        rb, ub = residual(nb, cal_b), usable(nb, pb.value)
        print(f'    bands   chi2/N = {chi2(nb, cal_b) / max(int(ub.sum()), 1):7.2f}'
              f'   rms = {100 * np.nanstd(rb[ub]):5.2f}%   n={int(ub.sum())}')

        # the held-out windows use the scalar solved on the BANDS, never refit
        ps = predict(th, [spec], grid)[0]
        model_spec = ps.value * c[0]
        hres = {}
        for k, h in held.items():
            hm = h.mask & np.isfinite(model_spec) & (model_spec > 0)
            r = np.full_like(model_spec, np.nan)
            r[hm] = (h.flux[hm] - model_spec[hm]) / model_spec[hm]
            hres[k] = r
            print(f'    {k:<7} PREDICTED  median {100 * np.nanmedian(r[hm]):+6.2f}%'
                  f'   rms {100 * np.nanstd(r[hm]):5.2f}%   n={int(hm.sum())}')

        rx = None
        if xs is not None:
            px = predict(th, [xs], grid)[0]
            cal_x, _ = solve(xs, px.value)
            rx, ux = residual(xs, cal_x), usable(xs, px.value)
            xsl_cal[label] = cal_x
            print(f'    xsl     chi2/N = '
                  f'{chi2(xs, cal_x) / max(int(ux.sum()), 1):7.2f}'
                  f'   rms = {100 * np.nanstd(rx[ux]):5.2f}%   n={int(ux.sum())}')
        results[label] = dict(cal_b=cal_b, rb=rb, ub=ub, scalar=c[0],
                              model_spec=model_spec, hres=hres)

    figure(a.star, nb, spec, xs, held, results, xsl_cal, row, bands)


def figure(star, nb, spec, xs, held, results, xsl_cal, row, bands):
    labels = list(results)
    fig = plt.figure(figsize=(12.5, 14))
    gs = fig.add_gridspec(5, 2, height_ratios=[1.9, 1.1, 1.5, 1.5, 1.4],
                          hspace=.45, wspace=.22)
    fig.patch.set_facecolor(SURFACE)
    m = spec.mask

    def shade(ax):
        for _, lo, hi in bands:
            ax.axvspan(lo, hi, color=BAND_C, alpha=.10, lw=0)
        for win in (BREAK_WINDOW, PASCHEN_WINDOW):
            ax.axvspan(*win, color=HELD_C, alpha=.13, lw=0)

    ax, rax = fig.add_subplot(gs[0, :]), fig.add_subplot(gs[1, :])
    for a_ in (ax, rax):
        style(a_)
        shade(a_)
        for lam in (BALMER, PASCHEN):
            a_.axvline(lam, color=MUTED, ls='--', lw=1)
    ax.plot(spec.wavelength[m], spec.flux[m], color=OBS_C, lw=1.2,
            label='NGSL observed')
    for lab, c in zip(labels, (MOD_C, MOD2_C)):
        ax.plot(spec.wavelength[m], results[lab]['model_spec'][m], color=c,
                lw=1.0, label=f'model, {lab} (scalar from bands)')
    ax.set_yscale('log')
    ax.set_ylabel(r'F$_\lambda$', fontsize=9, color=INK)
    ax.set_title(f'{star} — blue = bands used for conditioning;  '
                 'red = held out and PREDICTED', fontsize=10, color=INK)
    ax.legend(fontsize=8, loc='upper right', framealpha=.92)
    for lab, c in zip(labels, (MOD_C, MOD2_C)):
        r = (spec.flux - results[lab]['model_spec']) / results[lab]['model_spec']
        rax.plot(spec.wavelength[m], r[m] * 100, color=c, lw=.8, label=lab)
    rax.axhline(0, color=MUTED, lw=1)
    rax.set_ylim(-12, 12)
    rax.set_ylabel('(obs−model)/model [%]', fontsize=9, color=INK)
    rax.set_xlabel(r'Wavelength [$\AA$, vacuum]', fontsize=9, color=INK)
    rax.legend(fontsize=8, loc='upper left', ncol=2, framealpha=.92)

    for col, (nm, win, lam0) in enumerate((('Balmer', BREAK_WINDOW, BALMER),
                                           ('Paschen', PASCHEN_WINDOW, PASCHEN))):
        axz = fig.add_subplot(gs[2, col])
        style(axz)
        shade(axz)
        lo, hi = win[0] - 350, win[1] + 150
        s = m & (spec.wavelength > lo) & (spec.wavelength < hi)
        axz.plot(spec.wavelength[s], spec.flux[s], color=OBS_C, lw=1.3)
        for lab, c in zip(labels, (MOD_C, MOD2_C)):
            axz.plot(spec.wavelength[s], results[lab]['model_spec'][s],
                     color=c, lw=1.0)
        axz.axvline(lam0, color=MUTED, ls='--', lw=1)
        axz.set_xlim(lo, hi)
        med = 100 * np.nanmedian(results[labels[-1]]['hres'][nm])
        axz.set_title(f'{nm} — HELD OUT (node residual {med:+.1f}%)',
                      fontsize=9, color=INK)
        axz.set_xlabel(r'$\lambda$ [$\AA$]', fontsize=8, color=INK)
        axz.set_ylabel(r'F$_\lambda$', fontsize=8, color=INK)

    if xs is not None:
        for col, (lo, hi, ttl) in enumerate((
                (4290., 4400., r'XSL: H$\gamma$ (continuum marginalised)'),
                (5150., 5220., 'XSL: Mg I b'))):
            axx = fig.add_subplot(gs[3, col])
            style(axx)
            s = xs.mask & (xs.wavelength > lo) & (xs.wavelength < hi)
            axx.plot(xs.wavelength[s], xs.flux[s], color=OBS_C, lw=1.0,
                     label='XSL')
            for lab, c in zip(labels, (MOD_C, MOD2_C)):
                axx.plot(xs.wavelength[s], xsl_cal[lab][s], color=c, lw=.9,
                         label=f'model, {lab}')
            axx.set_xlim(lo, hi)
            axx.set_title(ttl, fontsize=9, color=INK)
            axx.set_xlabel(r'$\lambda$ [$\AA$]', fontsize=8, color=INK)
            if col == 0:
                axx.legend(fontsize=7, loc='lower left', framealpha=.92)

    axb = fig.add_subplot(gs[4, :])
    style(axb)
    lam = np.array([f.wave_effective for f in nb.filters])
    for lab, c in zip(labels, (MOD_C, MOD2_C)):
        u = results[lab]['ub']
        axb.errorbar(lam[u], results[lab]['rb'][u] * 100,
                     yerr=100 * nb.uncertainty[u] / nb.flux[u],
                     fmt='o-', color=c, lw=1.1, ms=5, capsize=3, label=lab)
    axb.axhline(0, color=MUTED, lw=1)
    for win in (BREAK_WINDOW, PASCHEN_WINDOW):
        axb.axvspan(*win, color=HELD_C, alpha=.13, lw=0)
    axb.set_xlabel(r'band effective $\lambda$ [$\AA$]', fontsize=9, color=INK)
    axb.set_ylabel('(obs−model)/model [%]', fontsize=9, color=INK)
    axb.set_title('NGSL band residuals — the dust lever '
                  '(error bars = 1% calibration floor)', fontsize=9, color=INK)
    axb.legend(fontsize=8, framealpha=.92)

    fig.suptitle(f'{star}   conditioning = NGSL bands + XSL lines;   '
                 f'Balmer and Paschen held out   '
                 f'(nominal {row["teff_ngsl"]} / {row["logg_ngsl"]} / '
                 f'{row["mh_ngsl"]})', fontsize=11, color=INK)
    out = ROOT / 'figures' / f'predict_check_{star}.png'
    fig.savefig(out, dpi=170, facecolor=SURFACE, bbox_inches='tight')
    plt.close(fig)
    print(f'\n  -> {out.relative_to(ROOT)}')


if __name__ == '__main__':
    main()
