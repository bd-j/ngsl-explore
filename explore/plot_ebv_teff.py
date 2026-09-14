"""The chi^2 surface in Teff and E(B-V), for each leg of the conditioning set.

Three panels, because the whole design rests on the two legs constraining
different directions:

  NGSL bands   the continuum. Teff and E(B-V) are covariant here -- both tilt
               it -- so this panel should show a long diagonal valley, not a
               closed minimum. Its slope is the degeneracy.
  XSL lines    locally normalised line profiles, continuum marginalised away.
               A smooth reddening law cannot change them, so this panel should
               show near-VERTICAL contours: a Teff constraint that does not care
               about dust. That is what makes the break prediction possible.
  combined     where the diagonal valley meets the vertical wall.

Teff runs over the grid's own nodes, no interpolation, exactly as the node scan
will.

COMBINING THE TWO. Adding raw chi^2 would let XSL dominate by sheer pixel count
(2476 against 13 bands), so each leg is rescaled by its own chi^2_min/dof --
the standard "inflate the errors until the fit is acceptable" convention.

The inflation is floored at 1: it may widen an error bar, never shrink one. The
NGSL bands come out at chi^2/n = 0.14, and rescaling that to 1 would DEFLATE
their errors by a factor 2.7, claiming a precision the 1% calibration floor was
deliberately set not to claim. A chi^2/dof below 1 means the error bars are
conservative, which is not licence to shrink them.

The Delta-chi^2 contours are drawn at the nominal 2-parameter levels (2.30,
6.17, 11.8) but they are NOT calibrated confidence regions: the residuals are
correlated, so the true regions are larger. Read them as a relative map.

    python3 explore/plot_ebv_teff.py --star HD194453
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
from fitting.observations import conditioning_set
from fitting.predict import predict
from fitting.calibration import solve, usable, chi2

ROOT = Path(__file__).resolve().parent.parent
SURFACE, INK, MUTED, GRIDC = '#fcfcfb', '#22262b', '#6b7280', '#dfe3e8'
CAT_C = {'teff_ngsl': '#2a78d6', 'teff_xsl': '#7a3fa8',
         'ebv_phot': '#1a9e5c', 'ebv_miles': '#d68910',
         'ebv_sf11': '#c0392b', 'ebv_sfd98': '#c0392b'}
LEVELS = (2.30, 6.17, 11.8)          # nominal 1/2/3 sigma, two parameters


def style(ax):
    ax.set_facecolor(SURFACE)
    ax.tick_params(labelsize=8, colors=MUTED)
    for s in ax.spines.values():
        s.set_color(GRIDC)


def sample_row(star):
    for r in csv.DictReader(open(ROOT / 'data' / 'sample.csv')):
        if r['star'] == star:
            return r
    raise KeyError(star)


def fnum(x):
    try:
        v = float(x)
        return v if np.isfinite(v) else None
    except (TypeError, ValueError):
        return None


def surface(obs, grid, teffs, ebvs, logg, mh, vsini):
    """-> chi2[nT, nE] and the per-point count, for one observation."""
    out = np.full((len(teffs), len(ebvs)), np.nan)
    ndata = 0
    for i, t in enumerate(teffs):
        for j, e in enumerate(ebvs):
            th = dict(teff=float(t), logg=logg, mh=mh, ebv=float(e), vsini=vsini)
            try:
                p = predict(th, [obs], grid)[0]
                cal, _ = solve(obs, p.value)
            except (ValueError, np.linalg.LinAlgError):
                continue
            u = usable(obs, cal)
            if u.sum() < 3:
                continue
            out[i, j] = chi2(obs, cal)
            ndata = max(ndata, int(u.sum()))
    return out, ndata


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--star', default='HD194453')
    ap.add_argument('--vsini', type=float, default=0.0)
    ap.add_argument('--ebv-max', type=float, default=0.12)
    ap.add_argument('--nebv', type=int, default=41)
    a = ap.parse_args()

    row, grid = sample_row(a.star), Grid()
    logg = float(grid.logg[np.argmin(np.abs(grid.logg - float(row['logg_ngsl'])))])
    mh = float(grid.mh[np.argmin(np.abs(grid.mh - float(row['mh_ngsl'])))])
    teffs = grid.teff
    ebvs = np.linspace(0.0, a.ebv_max, a.nebv)
    print(f'{a.star}: {len(teffs)} Teff nodes x {len(ebvs)} E(B-V), '
          f'at log g={logg:.2f}, [M/H]={mh:+.2f}, v sin i={a.vsini:.0f}')

    cond = conditioning_set(a.star)
    surf, ndat = {}, {}
    for o in cond:
        surf[o.name], ndat[o.name] = surface(o, grid, teffs, ebvs, logg, mh,
                                             a.vsini)
        k = np.nanmin(surf[o.name])
        print(f'  {o.name:<12} n={ndat[o.name]:>5}  chi2_min={k:.1f}  '
              f'chi2_min/n={k / max(ndat[o.name], 1):.2f}')

    # inflate each leg by its own chi2_min/dof before adding; never deflate
    scaled = {}
    for k, c in surf.items():
        cmin = np.nanmin(c)
        infl = max(1.0, cmin / max(ndat[k], 1))
        scaled[k] = (c - cmin) / infl
        print(f'  {k:<12} error inflation x{np.sqrt(infl):.2f}'
              + ('  (floored at 1)' if infl == 1.0 else ''))
    # Measure the degeneracy slope from the bands surface: the Teff that
    # minimises chi^2 at each E(B-V). Independently predicted as +126 K per
    # 0.01 mag from a synthetic test, so this is a check on that.
    if 'ngsl_bands' in scaled:
        c = scaled['ngsl_bands']
        ridge = np.array([teffs[np.nanargmin(c[:, j])] for j in range(len(ebvs))])
        ok = np.isfinite(ridge)
        slope = np.polyfit(ebvs[ok], ridge[ok], 1)[0] * 0.01
        print(f'  degeneracy ridge: {slope:+.0f} K per 0.01 mag of E(B-V) '
              f'(synthetic prediction was +126)')
    total = np.nansum(np.dstack(list(scaled.values())), axis=2)
    panels = [(o.name, scaled[o.name]) for o in cond] + [('combined', total)]

    fig, axes = plt.subplots(1, len(panels), figsize=(4.9 * len(panels), 5.0),
                             sharey=True, constrained_layout=True)
    fig.patch.set_facecolor(SURFACE)
    axes = np.atleast_1d(axes)
    X, Y = np.meshgrid(teffs, ebvs, indexing='ij')

    for ax, (name, d) in zip(axes, panels):
        style(ax)
        d = d - np.nanmin(d)
        from matplotlib.colors import LogNorm
        im = ax.pcolormesh(X, Y, np.clip(d, 0.3, 300), cmap='magma',
                           shading='nearest', alpha=.9,
                           norm=LogNorm(vmin=0.3, vmax=300))
        cs = ax.contour(X, Y, d, levels=LEVELS, colors='w', linewidths=1.1)
        ax.clabel(cs, fmt={l: s for l, s in zip(LEVELS, ('1σ', '2σ', '3σ'))},
                  fontsize=7)
        i, j = np.unravel_index(np.nanargmin(d), d.shape)
        # For XSL the E(B-V) direction is flat by construction, so the minimum's
        # E(B-V) is not a measurement and must not be presented as one.
        flat = name == 'xsl'
        lab = (f'min: {teffs[i]:.0f} K, E(B−V) unconstrained' if flat
               else f'min: {teffs[i]:.0f} K, {ebvs[j]:.3f}')
        ax.plot(teffs[i], ebvs[j] if not flat else np.mean(ebvs), marker='*',
                ms=15, color='w', markeredgecolor='k', markeredgewidth=.6,
                zorder=5, label=lab)
        title = {'ngsl_bands': 'NGSL bands — continuum (dust)',
                 'xsl': 'XSL lines — dust-immune',
                 'combined': 'combined (each leg rescaled)'}.get(name, name)
        ax.set_title(title, fontsize=10, color=INK)
        ax.set_xlabel('Teff [K]', fontsize=9, color=INK)
        ax.legend(fontsize=7.5, loc='upper left', framealpha=.85)

        # catalog markers
        for key, lab in (('teff_ngsl', 'NGSL Teff'), ('teff_xsl', 'XSL Teff')):
            v = fnum(row.get(key))
            if v and teffs.min() <= v <= teffs.max():
                ax.axvline(v, color=CAT_C[key], ls='--', lw=1.2, alpha=.9)
                ax.text(v, ebvs[-1], f' {lab}', rotation=90, va='top',
                        ha='right', fontsize=7, color=CAT_C[key])
        for key, lab in (('ebv_phot', 'photometric'), ('ebv_miles', 'MILES'),
                         ('ebv_sf11', 'SF11 (upper bound)')):
            v = fnum(row.get(key))
            if v is None:
                continue
            if ebvs.min() <= v <= ebvs.max():
                ax.axhline(v, color=CAT_C[key], ls=':', lw=1.2, alpha=.9)
                ax.text(teffs[0], v, f' {lab}', va='bottom', fontsize=7,
                        color=CAT_C[key])
            elif v < ebvs.min():
                # a negative photometric E(B-V) is unphysical and means
                # "consistent with zero"; say so rather than silently dropping it
                ax.annotate(f'{lab} {v:+.3f} (≈0, off scale)',
                            xy=(teffs[0], ebvs.min()), xytext=(3, 3),
                            textcoords='offset points', fontsize=7,
                            color=CAT_C[key], va='bottom')
    axes[0].set_ylabel('E(B−V) [mag]', fontsize=9, color=INK)
    fig.colorbar(im, ax=axes.tolist(), shrink=.85, label='Δχ² (log scale)')

    fig.suptitle(
        f'{a.star} — χ² surface in Teff and E(B−V)   '
        f'(log g={logg:.2f}, [M/H]={mh:+.2f}, v sin i={a.vsini:.0f} km/s held)\n'
        'contours are nominal 2-parameter levels, NOT calibrated confidence '
        'regions — residuals are correlated',
        fontsize=11, color=INK, linespacing=1.5)
    out = ROOT / 'figures' / f'ebv_teff_{a.star}.png'
    fig.savefig(out, dpi=170, facecolor=SURFACE)
    plt.close(fig)
    print(f'  -> {out.relative_to(ROOT)}')


if __name__ == '__main__':
    main()
