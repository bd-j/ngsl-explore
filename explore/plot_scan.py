"""Read results/<star>/scan.npz and answer the project's question.

The scan stores, for every one of the 1705 grid nodes, the whole conditional
curve in each nuisance direction. This turns that into:

  1. where the data puts the star in (Teff, log g, [M/H]) -- and in particular
     whether it piles up against the grid's [M/H] = -0.5 floor, which is the
     leading suspect for the unphysical reddening the fixed-[M/H] sweep found;
  2. what E(B-V) the bands want once Teff, log g and [M/H] are all free;
  3. THE HELD-OUT BREAK PREDICTION, marginalised rather than plugged in.

Panel 4 is the point of the whole project. Each node+E(B-V) combination
predicts a Balmer and a Paschen residual in a region the fit never saw, and
also has a Delta chi^2. Plotting prediction against Delta chi^2 shows whether
every acceptable fit predicts the same break -- which is a far stronger
statement than the best fit happening to land near zero.

A NOTE ON SHARPNESS. The XSL leg has ~2400 pixels against the bands' 13, and
independent-pixel chi^2 treats every one as independent when the residuals are
correlated over the LSF and by the continuum. The nominal Delta chi^2 = 1
contour is therefore far too tight -- single-node in the fixed-[M/H] sweep.
Contours here are drawn at wider levels and labelled nominal; the honest
uncertainty on the break prediction is the SPREAD of panel 4, not the curvature
at its minimum.

    python3 explore/plot_scan.py --star HD194453
    python3 explore/plot_scan.py --all
"""
import argparse
import csv
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from fitting.scan import posterior as scan_posterior

from common.figpath import figure_path, below_grid

ROOT = Path(__file__).resolve().parent.parent
SURFACE, INK, MUTED, GRIDC = '#fcfcfb', '#22262b', '#6b7280', '#dfe3e8'
C_BAL, C_PAS, C_CAT = '#2a78d6', '#eb6834', '#1a9e5c'
# The legs are combined as chi^2/dof (fitting.scan.posterior), so the displayed
# statistic is Delta[sum of chi^2/dof], not a chi^2. Multiplying by the band
# count puts it in BAND-EQUIVALENT chi^2 -- "how many band-pixels' worth of
# misfit" -- which is the only unit here anyone has a feel for, and which makes
# the usual 2-parameter levels readable again. They are still not calibrated
# confidence regions: the residuals are correlated.
BAND_EQUIV = 13.0
LEVELS = tuple(x / BAND_EQUIV for x in (2.30, 6.17, 11.8))
GOOD = 9.0 / BAND_EQUIV   # within which a node counts as "acceptable"


def style(ax):
    ax.set_facecolor(SURFACE)
    ax.grid(alpha=.22, color=GRIDC, lw=.7)
    ax.tick_params(labelsize=8, colors=MUTED)
    for s in ax.spines.values():
        s.set_color(GRIDC)


def sample_row(star):
    for r in csv.DictReader(open(ROOT / 'data' / 'sample.csv')):
        if r['star'] == star:
            return r
    return {}


def fnum(x):
    try:
        v = float(x)
        return v if np.isfinite(v) else None
    except (TypeError, ValueError):
        return None


def load(star, scale='profile', weight='inverse_dof'):
    """-> dict with the total Delta chi^2 over (node, E) and the axes.

    The posterior itself is built by fitting.scan.posterior, so this figure and
    explore/check_predict.py cannot disagree about which node is best.
    """
    d, total = scan_posterior(star, scale, weight)
    return dict(d=d, lnl=total, teff=d['teff'], logg=d['logg'], mh=d['mh'],
                ebv=d['ebv'], vsini=d['vsini'],
                dchi2=2.0 * (np.nanmax(total) - total),
                rbal=d['resid_balmer'], rpas=d['resid_paschen'])


def profile(dchi2, keep_axes):
    """Min over every axis except `keep_axes` (a tuple), NaN-safe."""
    drop = tuple(a for a in range(dchi2.ndim) if a not in keep_axes)
    with np.errstate(all='ignore'):
        return np.nanmin(dchi2, axis=drop)


def heat(ax, X, Y, d):
    d = d - np.nanmin(d)
    lo, hi = 0.3 / BAND_EQUIV, 300.0 / BAND_EQUIV
    im = ax.pcolormesh(X, Y, np.clip(d, lo, hi), cmap='magma',
                       shading='nearest', norm=LogNorm(vmin=lo, vmax=hi),
                       alpha=.92)
    cs = ax.contour(X, Y, d, levels=LEVELS, colors='w', linewidths=1.0)
    ax.clabel(cs, fmt=dict(zip(LEVELS, ('1σ', '2σ', '3σ'))), fontsize=7)
    i, j = np.unravel_index(np.nanargmin(d), d.shape)
    return im, i, j


def figure(star, S, out):
    row = sample_row(star)
    teff, logg, mh, ebv = S['teff'], S['logg'], S['mh'], S['ebv']
    dchi2 = S['dchi2']
    fig, axes = plt.subplots(1, 4, figsize=(19.6, 4.9), constrained_layout=True)
    fig.patch.set_facecolor(SURFACE)
    for ax in axes:
        style(ax)

    # 1. Teff x [M/H] -- does it pile up on the grid floor?
    ax = axes[0]
    X, Y = np.meshgrid(teff, mh, indexing='ij')
    im, i, j = heat(ax, X, Y, profile(dchi2, (0, 2)))
    ax.plot(teff[i], mh[j], marker='*', ms=15, color='w', markeredgecolor='k',
            markeredgewidth=.6, zorder=5,
            label=f'{teff[i]:.0f} K, [M/H]={mh[j]:+.2f}')
    zc = fnum(row.get('mh_ngsl'))
    if zc is not None:
        if mh.min() <= zc <= mh.max():
            ax.axhline(zc, color=C_CAT, ls=':', lw=1.2)
            ax.text(teff[0], zc, ' catalog', fontsize=7, color=C_CAT, va='bottom')
        else:
            ax.annotate(f'catalog [M/H]={zc:+.2f} — off the grid',
                        xy=(teff[0], mh.min()), xytext=(3, 3),
                        textcoords='offset points', fontsize=7, color=C_CAT)
    ax.set_ylabel('[M/H]', fontsize=9, color=INK)
    ax.set_title('Teff × [M/H]  (log g, E(B−V) profiled)', fontsize=10, color=INK)
    ax.legend(fontsize=7.5, loc='upper left', framealpha=.88)

    # 2. Teff x log g
    ax = axes[1]
    X, Y = np.meshgrid(teff, logg, indexing='ij')
    im, i, j = heat(ax, X, Y, profile(dchi2, (0, 1)))
    ax.plot(teff[i], logg[j], marker='*', ms=15, color='w', markeredgecolor='k',
            markeredgewidth=.6, zorder=5,
            label=f'{teff[i]:.0f} K, log g={logg[j]:.2f}')
    gc = fnum(row.get('logg_ngsl'))
    if gc is not None and logg.min() <= gc <= logg.max():
        ax.axhline(gc, color=C_CAT, ls=':', lw=1.2)
        ax.text(teff[0], gc, ' catalog', fontsize=7, color=C_CAT, va='bottom')
    ax.set_ylabel('log g', fontsize=9, color=INK)
    ax.set_title('Teff × log g  ([M/H], E(B−V) profiled)', fontsize=10, color=INK)
    ax.legend(fontsize=7.5, loc='upper left', framealpha=.88)

    # 3. Teff x E(B-V)
    ax = axes[2]
    X, Y = np.meshgrid(teff, ebv, indexing='ij')
    im, i, j = heat(ax, X, Y, profile(dchi2, (0, 3)))
    ax.plot(teff[i], ebv[j], marker='*', ms=15, color='w', markeredgecolor='k',
            markeredgewidth=.6, zorder=5,
            label=f'{teff[i]:.0f} K, E(B−V)={ebv[j]:.3f}')
    for key, lab in (('ebv_phot', 'photometric'), ('ebv_sf11', 'SF11 column')):
        v = fnum(row.get(key))
        if v is None or (key == 'ebv_sf11' and row.get('ebv_map_useful') == 'no'):
            continue
        if ebv.min() <= v <= ebv.max():
            ax.axhline(v, color=C_CAT if key == 'ebv_phot' else '#c0392b',
                       ls=':', lw=1.2)
            ax.text(teff[0], v, f' {lab}', fontsize=7, va='bottom',
                    color=C_CAT if key == 'ebv_phot' else '#c0392b')
    ax.set_ylabel('E(B−V)', fontsize=9, color=INK)
    ax.set_title('Teff × E(B−V)  (log g, [M/H] profiled)', fontsize=10, color=INK)
    ax.legend(fontsize=7.5, loc='upper left', framealpha=.88)
    fig.colorbar(im, ax=axes[:3].tolist(), shrink=.85,
                 label='Δ(Σ χ²/dof)  (log)')

    # 4. THE POINT: held-out prediction vs how well that model fits
    ax = axes[3]
    ok = np.isfinite(dchi2)
    stats = {}
    for arr, c, lab in ((S['rbal'], C_BAL, 'Balmer'), (S['rpas'], C_PAS, 'Paschen')):
        m = ok & np.isfinite(arr)
        x, y = 100 * arr[m], dchi2[m]
        ax.scatter(x, y, s=4, c=c, alpha=.25, lw=0, label=f'{lab} (held out)')
        g = y <= GOOD
        if g.any():
            lo, hi = np.percentile(x[g], [2.5, 97.5])
            best = x[g][np.argmin(y[g])]
            stats[lab] = (best, lo, hi, int(g.sum()))
            ax.axvspan(lo, hi, color=c, alpha=.10, lw=0)
    ax.axvline(0, color=INK, lw=1.2)
    ax.axhline(GOOD, color=MUTED, ls='--', lw=1)
    ax.text(ax.get_xlim()[0], GOOD, f' {GOOD * BAND_EQUIV:.0f} band-equiv. χ²',
            fontsize=7, color=MUTED, va='bottom')
    ax.set_yscale('log')
    ax.set_ylim(0.3 / BAND_EQUIV, 3e3 / BAND_EQUIV)
    ax.set_xlim(-12, 12)
    ax.set_xlabel('predicted held-out residual  (obs−model)/model [%]',
                  fontsize=9, color=INK)
    ax.set_ylabel('Δ(Σ χ²/dof) of that model', fontsize=9, color=INK)
    ax.set_title('Held-out break prediction vs fit quality', fontsize=10,
                 color=INK)
    ax.legend(fontsize=7.5, loc='upper right', framealpha=.88)
    txt = '\n'.join(f'{k}: {v[0]:+.2f}%  [{v[1]:+.2f}, {v[2]:+.2f}] '
                    f'over {v[3]} models' for k, v in stats.items())
    if txt:
        ax.text(.03, .03, txt, transform=ax.transAxes, fontsize=7.5, color=INK,
                va='bottom', bbox=dict(fc='white', ec=GRIDC, alpha=.9))

    fig.suptitle(
        f'{star} ({row.get("tier", "?")}) — all {S["lnl"][..., 0].size} grid nodes, '
        'no interpolation;  E(B−V) and v sin i fitted at every node\n'
        'legs combined as χ²/dof, so 13 NGSL bands carry the same total weight '
        'as 2342 XSL pixels;  contours are nominal and NOT calibrated — '
        'the spread in panel 4, not the curvature, is the honest uncertainty',
        fontsize=10.5, color=INK, linespacing=1.5)
    fig.savefig(out, dpi=160, facecolor=SURFACE)
    plt.close(fig)
    return stats


def leg_tension(d, thresh=1.5):
    """Can ANY node fit both legs at once? -> dict of the weighting-free facts.

    The leg weighting decides how a conflict between the continuum (bands) and
    the line profiles (XSL) is resolved. It cannot tell you whether there is a
    conflict. This can: it asks whether the grid contains any node at which both
    legs reach chi^2/N < thresh, with E(B-V) and v sin i free.

    Where the answer is no, the weighting is choosing WHICH failure you see
    rather than removing it -- and the tell is v s in i running to the ceiling,
    which is broadening being spent to wash out model lines that are too strong
    because [M/H] is pinned at the grid floor.
    """
    nb, nx = int(d['n_bands']), int(d['n_xsl'])
    cb = np.nanmin(d['chi2_bands'] / nb, axis=-1)        # best over E, per node
    cx = np.nanmin(d['chi2_xsl'] / nx, axis=-1)          # best over v, per node
    both = np.isfinite(cb) & np.isfinite(cx) & (cb < thresh) & (cx < thresh)
    worst = np.maximum(cb, cx)
    i, j, k = np.unravel_index(np.nanargmin(worst), worst.shape)
    iv = int(np.nanargmin(np.where(np.isfinite(d['chi2_xsl'][i, j, k]),
                                   d['chi2_xsl'][i, j, k], np.inf)))
    return dict(n_joint=int(both.sum()), floor_bands=float(np.nanmin(cb)),
                floor_xsl=float(np.nanmin(cx)),
                joint_bands=float(cb[i, j, k]), joint_xsl=float(cx[i, j, k]),
                joint_vsini=float(d['vsini'][iv]),
                feasible=bool(both.sum() > 0))


def sample_summary(rows, out):
    """The sample-level result: is a failed break prediction PREDICTABLE?

    Panel 1 is the finding. The conditioning fit quality (NGSL bands, 13 points
    the fit DID see) tracks the held-out break error (a region it never saw) at
    r = +0.96. So the break prediction can be trusted exactly when the
    conditioning fit is acceptable -- which is checkable without ever looking at
    the held-out data, and is what makes this a prediction rather than a hope.

    Panel 3 says why the failures fail: every one of them is a star the grid
    cannot reach in [M/H], so the continuum blanketing is wrong. XSL cannot
    detect that (its continuum is marginalised away, and it fits those stars at
    chi2/n < 1.2), which is precisely why the bands leg has to be there.
    """
    fig, axes = plt.subplots(1, 3, figsize=(15.6, 5.0), constrained_layout=True)
    fig.patch.set_facecolor(SURFACE)
    for ax in axes:
        style(ax)
    mk = {'primary': 'o', 'secondary': 's'}

    ax = axes[0]
    x = np.array([r['chi2n_bands'] for r in rows])
    y = np.array([abs(r['balmer']) for r in rows])
    for r in rows:
        c = '#c0392b' if r['below_grid'] else C_BAL
        ax.plot(r['chi2n_bands'], abs(r['balmer']), marker=mk[r['tier']], ms=7,
                color=c, mfc='none' if r['below_grid'] else c, mew=1.5, lw=0)
        ax.annotate(r['star'].replace('HD', ''),
                    (r['chi2n_bands'], abs(r['balmer'])), xytext=(5, 2),
                    textcoords='offset points', fontsize=6.5, color=MUTED)
    ax.axvline(1.0, color=MUTED, ls=':', lw=1.1)
    ax.text(1.0, ax.get_ylim()[1], ' acceptable fit ', rotation=90, fontsize=7,
            color=MUTED, va='top', ha='right')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('NGSL band χ²/N  —  the CONDITIONING fit', fontsize=9, color=INK)
    ax.set_ylabel('|held-out Balmer residual| [%]', fontsize=9, color=INK)
    ax.set_title(f'A bad break prediction is visible in the conditioning data\n'
                 f'r = {np.corrcoef(np.log10(x), np.log10(y))[0, 1]:+.2f} '
                 f'(log–log, n={len(rows)})', fontsize=9.5, color=INK)
    ax.plot([], [], marker='o', lw=0, color='#c0392b', mfc='none', mew=1.5,
            label='catalog [M/H] < −0.5 — below the grid')
    ax.plot([], [], marker='o', lw=0, color=C_BAL, label='inside the grid')
    ax.legend(fontsize=7.5, loc='upper left', framealpha=.9)

    ax = axes[1]
    order = sorted(rows, key=lambda r: r['balmer'])
    for n, r in enumerate(order):
        for key, c, off in (('balmer', C_BAL, -.16), ('paschen', C_PAS, .16)):
            lo, hi = r[key + '_lo'], r[key + '_hi']
            ax.plot([lo, hi], [n + off] * 2, color=c, lw=3, alpha=.35,
                    solid_capstyle='round')
            ax.plot(r[key], n + off, marker=mk[r['tier']], ms=5.5, color=c)
    ax.axvline(0, color=INK, lw=1.2)
    ax.axvspan(-1, 1, color=MUTED, alpha=.08, lw=0)
    ax.set_yticks(np.arange(len(order)))
    ax.set_yticklabels([r['star'] + ('*' if r['below_grid'] else '')
                        for r in order], fontsize=7.5, color=MUTED)
    ax.set_ylim(-.7, len(order) - .3)
    ax.set_xlim(-20, 32)
    ax.plot([], [], marker='o', lw=0, color=C_BAL, label='Balmer')
    ax.plot([], [], marker='o', lw=0, color=C_PAS, label='Paschen')
    ax.legend(fontsize=7.5, loc='lower right', framealpha=.9)
    ax.set_xlabel('held-out residual [%]  (bar = 95% of acceptable models)',
                  fontsize=9, color=INK)
    ax.set_title('The held-out prediction, per star', fontsize=10, color=INK)

    ax = axes[2]
    for n, r in enumerate(order):
        ax.plot(r['mh'], n, marker=mk[r['tier']], ms=7, lw=0,
                color='#c0392b' if r['below_grid'] else C_BAL,
                mfc='none' if r['below_grid'] else None, mew=1.5)
        if r['mh_cat'] is not None:
            ax.plot(max(r['mh_cat'], -2.1), n, marker='|', ms=9, color=C_CAT)
    ax.axvline(-0.5, color='#c0392b', ls='--', lw=1.2)
    ax.text(-0.5, len(order) - .5, ' grid floor ', rotation=90, fontsize=7,
            color='#c0392b', va='top', ha='right')
    ax.plot([], [], marker='|', lw=0, color=C_CAT, label='catalog [M/H]')
    ax.legend(fontsize=7.5, loc='lower right', framealpha=.9)
    ax.set_yticks(np.arange(len(order)))
    ax.set_yticklabels([], fontsize=7)
    ax.set_ylim(-.7, len(order) - .3)
    ax.set_xlim(-2.2, .5)
    ax.set_xlabel('[M/H] of the best node', fontsize=9, color=INK)
    ax.set_title('Why they fail: the grid cannot reach them', fontsize=10,
                 color=INK)

    ing = [r for r in rows if not r['below_grid']]
    lowz = [r for r in rows if r['below_grid']]
    mb = np.median([abs(r['balmer']) for r in ing]) if ing else float('nan')
    fig.suptitle(
        'Node scan — all 1705 grid nodes per star, no interpolation;  '
        'legs combined as χ²/dof;  E(B−V) and v sin i fitted at every node\n'
        f'* and open symbols = catalog [M/H] < −0.5, BELOW the grid floor '
        f'({len(lowz)} of {len(rows)}) — kept separate, not pooled.  '
        f'Inside the grid ({len(ing)}): held-out |Balmer| median {mb:.2f}%',
        fontsize=10.5, color=INK, linespacing=1.5)
    fig.savefig(out, dpi=165, facecolor=SURFACE)
    plt.close(fig)
    print(f'-> {out.relative_to(ROOT)}')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--star', default='HD194453')
    ap.add_argument('--all', action='store_true')
    a = ap.parse_args()
    stars = ([p.parent.name for p in sorted((ROOT / 'results').glob('*/scan.npz'))]
             if a.all else [a.star])
    rows = []
    for s in stars:
        try:
            S = load(s)
        except FileNotFoundError:
            print(f'{s}: no scan.npz')
            continue
        out = figure_path('scan', s)
        st = figure(s, S, out)
        d2 = S['dchi2']
        i, j, k, e = np.unravel_index(np.nanargmin(d2), d2.shape)
        edge = []
        if k in (0, len(S['mh']) - 1):
            edge.append('[M/H] edge')
        if i in (0, len(S['teff']) - 1):
            edge.append('Teff edge')
        if j in (0, len(S['logg']) - 1):
            edge.append('log g edge')
        if e in (0, len(S['ebv']) - 1):
            edge.append('E(B−V) edge')
        print(f'{s:<10} best node Teff={S["teff"][i]:.0f} log g={S["logg"][j]:.2f} '
              f'[M/H]={S["mh"][k]:+.2f}  E(B-V)={S["ebv"][e]:.3f}'
              + (f'   AT: {", ".join(edge)}' if edge else ''))
        for nm, v in st.items():
            print(f'             {nm:<8} held out {v[0]:+.2f}%  '
                  f'95% of acceptable models in [{v[1]:+.2f}, {v[2]:+.2f}]%')
        d = S['d']
        v = int(np.nanargmin(d['chi2_xsl'][i, j, k]))
        pz = profile(d2, (2,))
        rec = dict(star=s, tier=sample_row(s).get('tier', 'primary'),
                   teff=float(S['teff'][i]), logg=float(S['logg'][j]),
                   mh=float(S['mh'][k]), ebv=float(S['ebv'][e]),
                   vsini=float(S['vsini'][v]),
                   mh_cat=fnum(sample_row(s).get('mh_ngsl')),
                   pressed=bool(k == 0 and (pz[1] - pz[0]) > 2.30),
                   below_grid=bool(below_grid(s)),
                   chi2n_bands=float(d['chi2_bands'][i, j, k, e]) / int(d['n_bands']),
                   chi2n_xsl=float(d['chi2_xsl'][i, j, k, v]) / int(d['n_xsl']))
        for key, nm in (('balmer', 'Balmer'), ('paschen', 'Paschen')):
            b = st.get(nm)
            rec[key], rec[key + '_lo'], rec[key + '_hi'] = (
                b[0], b[1], b[2]) if b else (np.nan,) * 3
        lt = leg_tension(d)
        rec.update(lt)
        vmax = float(S['vsini'][-1])
        if rec['vsini'] >= vmax:
            print(f'             v sin i AT THE CEILING ({vmax:.0f} km/s) — '
                  f'broadening is absorbing a model error, not measuring rotation')
        if not lt['feasible']:
            print(f'             NO NODE fits both legs at chi2/N < 1.5 '
                  f'(floors: bands {lt["floor_bands"]:.2f}, xsl '
                  f'{lt["floor_xsl"]:.2f}; best joint {lt["joint_bands"]:.2f}/'
                  f'{lt["joint_xsl"]:.2f}) — the leg weighting is choosing which '
                  f'failure you see, not removing it')
        rows.append(rec)
        print(f'             -> {out.relative_to(ROOT)}')
    if len(rows) > 2:
        sample_summary(rows, ROOT / 'figures' / 'scan_sample.png')
        # Report the two groups SEPARATELY and never pool them: a star whose
        # catalog [M/H] is below the grid floor cannot be represented by any
        # node, so its residuals are not measuring the same thing.
        for lab, g in (('inside the grid', [r for r in rows if not r['below_grid']]),
                       ('BELOW the grid [M/H] floor',
                        [r for r in rows if r['below_grid']])):
            if not g:
                continue
            b = np.abs([r['balmer'] for r in g])
            pa = np.abs([r['paschen'] for r in g])
            cb = [r['chi2n_bands'] for r in g]
            nf = [r['star'] for r in g if not r['feasible']]
            vc = [r['star'] for r in g if r['vsini'] >= 300]
            print(f'\n  {lab} (n={len(g)}): {", ".join(r["star"] for r in g)}')
            print(f'    held out |Balmer| median {np.median(b):.2f}%  '
                  f'max {b.max():.2f}%   |Paschen| median {np.median(pa):.2f}%')
            print(f'    band chi2/N median {np.median(cb):.2f}  max {max(cb):.2f}')
            print(f'    no joint leg solution: {len(nf)}/{len(g)}'
                  + (f' ({", ".join(nf)})' if nf else ''))
            if vc:
                print(f'    v sin i at the ceiling: {", ".join(vc)}')
    return rows


if __name__ == '__main__':
    main()
