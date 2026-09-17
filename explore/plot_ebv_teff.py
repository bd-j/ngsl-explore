"""The chi^2 surface in Teff, E(B-V) and v sin i, for each leg of the
conditioning set.

Three panels, because the whole design rests on the two legs constraining
different directions:

  NGSL bands   the continuum, over (Teff, E(B-V)). Both parameters tilt it, so
               this panel should show a long diagonal valley, not a closed
               minimum. Its slope IS the degeneracy.
  XSL lines    locally normalised line profiles, continuum marginalised away,
               over (Teff, v sin i) -- a Teff constraint that does not care
               about dust. That is what makes the break prediction possible.
  combined     where the diagonal valley meets XSL's vertical wall.

Teff runs over the grid's own nodes, no interpolation, exactly as the node scan
will.

WHY THE TWO LEGS USE DIFFERENT AXES. A 3-D scan is unnecessary because the
dependences factorise, measured on HD194453 at 10250 K rather than assumed:

  XSL chi^2 over E(B-V) = 0 -> 0.12       Delta chi^2 = +2.0   (of 6362, 2342 px)
  band chi^2 over v sin i = 0 -> 300      Delta chi^2 = -0.04  (of 3.41, 13 bands)
  XSL chi^2 over v sin i = 0 -> 300       Delta chi^2 = +17000

So each leg is scanned only over the axes it responds to. The +2.0 discarded
from XSL is below the nominal 1-sigma level and is in any case a constraint on
XSL's CONTINUUM, which slit losses make untrustworthy and which the segmented
polynomial exists to marginalise away. It is re-measured per star and printed,
so the approximation is checked rather than trusted.

v sin i must be FITTED, not held at 0: it moves XSL's chi^2 by four orders of
magnitude more than the dust does, and a rotator forced to 0 would bias the one
Teff measurement the design depends on. It is profiled out of the XSL leg, and
the profile is the panel.

COMBINING THE TWO. Adding raw chi^2 would let XSL dominate by sheer pixel count
(~2400 against 13 bands), so each leg is rescaled by its own chi^2_min/dof --
the standard "inflate the errors until the fit is acceptable" convention.

The inflation is floored at 1: it may widen an error bar, never shrink one. The
NGSL bands come out at chi^2/n ~ 0.14, and rescaling that to 1 would DEFLATE
their errors by a factor 2.7, claiming a precision the 1% calibration floor was
deliberately set not to claim. A chi^2/dof below 1 means the error bars are
conservative, which is not licence to shrink them.

The Delta-chi^2 contours are drawn at the nominal 2-parameter levels (2.30,
6.17, 11.8) but they are NOT calibrated confidence regions: the residuals are
correlated, so the true regions are larger. Read them as a relative map.

    python3 explore/plot_ebv_teff.py --star HD194453
    python3 explore/plot_ebv_teff.py --all          # sample + summary figure
    python3 explore/plot_ebv_teff.py --summary-only # redraw from the CSV
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
from fitting.model import Grid
from fitting.observations import conditioning_set
from fitting.predict import predict
from fitting.calibration import solve, usable, chi2

from common.figpath import figure_path, below_grid
from common.sample import sample_row, opt_float

ROOT = Path(__file__).resolve().parent.parent
SURFACE, INK, MUTED, GRIDC = '#fcfcfb', '#22262b', '#6b7280', '#dfe3e8'
CAT_C = {'teff_ngsl': '#2a78d6', 'teff_xsl': '#7a3fa8',
         'ebv_phot': '#1a9e5c', 'ebv_miles': '#d68910',
         'ebv_sf11': '#c0392b', 'ebv_sfd98': '#c0392b'}
LEVELS = (2.30, 6.17, 11.8)          # nominal 1/2/3 sigma, two parameters
VMIN, VMAX = 0.3, 300.0              # shared colour scale across panels

# Coarse near 0 because XSL's FWHM (~31 km/s in the UVB) cannot resolve a
# rotation slower than ~15 km/s. Imported rather than copied so this and the
# node scan cannot drift apart.
from fitting.scan import VSINI_GRID
SCAN_CSV = ROOT / 'data' / 'ebv_teff_scan.csv'


def style(ax):
    ax.set_facecolor(SURFACE)
    ax.tick_params(labelsize=8, colors=MUTED)
    for s in ax.spines.values():
        s.set_color(GRIDC)


def sample_rows(tiers=('primary', 'secondary')):
    return [r for r in csv.DictReader(open(ROOT / 'data' / 'sample.csv'))
            if r['tier'] in tiers]


def one_chi2(obs, grid, theta):
    """chi^2 of one observation at one parameter vector, or (nan, 0)."""
    try:
        p = predict(theta, [obs], grid)[0]
        cal, _ = solve(obs, p.value)
    except (ValueError, np.linalg.LinAlgError):
        return np.nan, 0
    u = usable(obs, cal)
    if u.sum() < 3:
        return np.nan, 0
    return chi2(obs, cal), int(u.sum())


def surface(obs, grid, teffs, logg, mh, xkey, xvals, **fixed):
    """chi^2[nT, nX] over Teff and one other axis, plus the pixel count."""
    out = np.full((len(teffs), len(xvals)), np.nan)
    ndata = 0
    for i, t in enumerate(teffs):
        for j, x in enumerate(xvals):
            th = dict(teff=float(t), logg=logg, mh=mh, **{xkey: float(x)},
                      **fixed)
            out[i, j], n = one_chi2(obs, grid, th)
            ndata = max(ndata, n)
    return out, ndata


def nodes_for(row, grid):
    """Nearest grid node in log g and [M/H], with the catalog value kept.

    8 of 13 stars fall outside the grid, almost all in [M/H] (the grid stops at
    -0.5, the sample reaches -1.92). Clamping is the agreed first pass, but the
    clamped value is returned alongside the catalog one so the figure can say it
    is fitting at -0.5 a star catalogued at -1.9 rather than implying a fit.
    """
    lg_cat, mh_cat = opt_float(row['logg_ngsl']), opt_float(row['mh_ngsl'])
    lg = float(grid.logg[np.argmin(np.abs(grid.logg - lg_cat))])
    mh = float(grid.mh[np.argmin(np.abs(grid.mh - mh_cat))])
    # Clamped = the catalog value is OUTSIDE the grid, which is a real caveat.
    # Snapping to the nearest node is not: with a 0.2 dex step every value moves
    # by up to 0.1, and reporting that as clamping cried wolf on 13 of 13 stars.
    clamped = dict(
        logg=lg_cat if not grid.logg.min() <= lg_cat <= grid.logg.max() else None,
        mh=mh_cat if not grid.mh.min() <= mh_cat <= grid.mh.max() else None)
    return lg, mh, lg_cat, mh_cat, clamped


def interval(x, curve, level=1.0):
    """[lo, hi] where a 1-parameter profile sits within `level` of its minimum.

    The minimum is subtracted here rather than by the caller: a profile of the
    COMBINED surface does not reach zero, because the two legs minimise at
    different Teff, and comparing it to `level` raw returned an empty interval
    for every star.

    Returns node edges, so a single-node interval reports that node twice.
    Nominal: the residuals are correlated, so the true interval is wider.
    """
    c = np.asarray(curve, float)
    if not np.isfinite(c).any():
        return np.nan, np.nan
    ok = np.isfinite(c) & (c - np.nanmin(c) <= level)
    w = np.asarray(x, float)[ok]
    return float(w.min()), float(w.max())


def run_star(star, grid, a):
    """Sweep one star, draw its figure, return the row for the summary."""
    row = sample_row(star)
    logg, mh, lg_cat, mh_cat, clamped = nodes_for(row, grid)
    teffs = grid.teff
    ebvs = np.linspace(0.0, a.ebv_max, a.nebv)
    vsinis = VSINI_GRID
    clamp = ''.join(f'  {k} CLAMPED from {v:+.2f} (outside grid)'
                    for k, v in clamped.items() if v is not None)
    print(f'{star} ({row["tier"]}): {len(teffs)} Teff x {len(ebvs)} E(B-V) '
          f'x {len(vsinis)} v sin i, at log g={logg:.2f}, [M/H]={mh:+.2f}'
          + clamp)

    cond = {o.name: o for o in conditioning_set(star)}
    if 'ngsl_bands' not in cond or 'xsl' not in cond:
        print(f'  skipped: legs present = {sorted(cond)}')
        return None

    # Each leg over the axes it responds to (see the module docstring).
    c_b, n_b = surface(cond['ngsl_bands'], grid, teffs, logg, mh,
                       'ebv', ebvs, vsini=a.vsini_bands)
    c_x, n_x = surface(cond['xsl'], grid, teffs, logg, mh,
                       'vsini', vsinis, ebv=0.0)
    if not np.isfinite(c_b).any() or not np.isfinite(c_x).any():
        print('  skipped: no finite chi^2')
        return None

    kb, kx = np.nanmin(c_b), np.nanmin(c_x)
    print(f'  ngsl_bands   n={n_b:>5}  chi2_min={kb:9.1f}  '
          f'chi2_min/n={kb / max(n_b, 1):.2f}')
    print(f'  xsl          n={n_x:>5}  chi2_min={kx:9.1f}  '
          f'chi2_min/n={kx / max(n_x, 1):.2f}')

    # error inflation, floored at 1 -- may widen an error bar, never shrink one
    infl_b = max(1.0, kb / max(n_b, 1))
    infl_x = max(1.0, kx / max(n_x, 1))
    d_b = (c_b - kb) / infl_b
    d_x = (c_x - kx) / infl_x
    for nm, infl in (('ngsl_bands', infl_b), ('xsl', infl_x)):
        print(f'  {nm:<12} error inflation x{np.sqrt(infl):.2f}'
              + ('  (floored at 1)' if infl == 1.0 else ''))

    # Profile v sin i out of the XSL leg; the profile is a pure Teff constraint
    # and broadcasts across E(B-V).
    prof_x = np.nanmin(d_x, axis=1)
    vhat = np.array([vsinis[np.nanargmin(d_x[i])] if np.isfinite(d_x[i]).any()
                     else np.nan for i in range(len(teffs))])
    total = d_b + prof_x[:, None]

    # The degeneracy slope, from the bands: the Teff minimising chi^2 at each
    # E(B-V). Predicted +126 K per 0.01 mag from a synthetic test, so this is a
    # check on that rather than a free parameter.
    ridge = np.array([teffs[np.nanargmin(d_b[:, j])] if np.isfinite(d_b[:, j]).any()
                      else np.nan for j in range(len(ebvs))])
    ok = np.isfinite(ridge)
    slope = (np.polyfit(ebvs[ok], ridge[ok], 1)[0] * 0.01 if ok.sum() > 2
             else np.nan)
    print(f'  degeneracy ridge: {slope:+.0f} K per 0.01 mag of E(B-V)')

    it, ie = np.unravel_index(np.nanargmin(total), total.shape)
    iv = int(np.nanargmin(d_x[it])) if np.isfinite(d_x[it]).any() else 0
    t_hat, e_hat, v_hat = teffs[it], ebvs[ie], vsinis[iv]

    # Re-measure the factorisation cost AT THIS STAR's solution: how much XSL
    # chi^2 really moves across the whole E(B-V) range. Printed, not assumed.
    lo_e, _ = one_chi2(cond['xsl'], grid, dict(teff=float(t_hat), logg=logg,
                                               mh=mh, ebv=0.0,
                                               vsini=float(v_hat)))
    hi_e, _ = one_chi2(cond['xsl'], grid, dict(teff=float(t_hat), logg=logg,
                                               mh=mh, ebv=float(ebvs[-1]),
                                               vsini=float(v_hat)))
    dust_cost = (hi_e - lo_e) / infl_x
    print(f'  XSL dust sensitivity discarded: dchi2={dust_cost:+.2f} over '
          f'E(B-V)=0-{ebvs[-1]:.2f} (nominal 1sigma is 2.30)')

    t_lo, t_hi = interval(teffs, np.nanmin(total, axis=1))
    e_lo, e_hi = interval(ebvs, np.nanmin(total, axis=0))
    # v sin i CONDITIONAL on the best Teff, matching v_hat. Profiling it over
    # Teff instead gave a point estimate and an interval that were different
    # quantities: HD128801 reported 250 km/s with an interval of [300, 300].
    v_lo, v_hi = interval(vsinis, d_x[it])
    at = []
    if e_hat <= ebvs[0] + 1e-9:
        at.append('ebv_floor')
    if e_hat >= ebvs[-1] - 1e-9:
        at.append('ebv_ceiling')
    if t_hat <= teffs[0] or t_hat >= teffs[-1]:
        at.append('teff_edge')
    if v_hat >= vsinis[-1] - 1e-9:
        at.append('vsini_ceiling')
    print(f'  combined min: Teff={t_hat:.0f} [{t_lo:.0f},{t_hi:.0f}] K, '
          f'E(B-V)={e_hat:.3f} [{e_lo:.3f},{e_hi:.3f}], '
          f'v sin i={v_hat:.0f} [{v_lo:.0f},{v_hi:.0f}] km/s'
          + (f'   AT BOUNDARY: {"+".join(at)}' if at else ''))

    draw_star(star, row, teffs, ebvs, vsinis, d_b, d_x, total, logg, mh,
              clamped, t_hat, e_hat, v_hat, slope, dust_cost, vhat)

    return dict(
        star=star, tier=row['tier'],
        teff=f'{t_hat:.0f}', teff_lo=f'{t_lo:.0f}', teff_hi=f'{t_hi:.0f}',
        ebv=f'{e_hat:.4f}', ebv_lo=f'{e_lo:.4f}', ebv_hi=f'{e_hi:.4f}',
        vsini=f'{v_hat:.0f}', vsini_lo=f'{v_lo:.0f}', vsini_hi=f'{v_hi:.0f}',
        logg_node=f'{logg:.2f}', mh_node=f'{mh:+.2f}',
        logg_cat=f'{lg_cat:.2f}', mh_cat=f'{mh_cat:+.2f}',
        chi2_bands=f'{kb:.2f}', n_bands=n_b,
        chi2_xsl=f'{kx:.1f}', n_xsl=n_x,
        ridge_K_per_001=f'{slope:.0f}', xsl_dust_dchi2=f'{dust_cost:+.2f}',
        at_boundary='+'.join(at), clamped='+'.join(
            k for k, v in clamped.items() if v is not None),
        teff_ngsl=row['teff_ngsl'], teff_xsl=row['teff_xsl'],
        ebv_phot=row['ebv_phot'], ebv_miles=row['ebv_miles'],
        ebv_sf11=row['ebv_sf11'], ebv_map_useful=row['ebv_map_useful'])


def heat(ax, X, Y, d, levels=LEVELS):
    d = d - np.nanmin(d)
    im = ax.pcolormesh(X, Y, np.clip(d, VMIN, VMAX), cmap='magma',
                       shading='nearest', alpha=.9,
                       norm=LogNorm(vmin=VMIN, vmax=VMAX))
    cs = ax.contour(X, Y, d, levels=levels, colors='w', linewidths=1.1)
    ax.clabel(cs, fmt=dict(zip(levels, ('1σ', '2σ', '3σ'))), fontsize=7)
    return im


def draw_star(star, row, teffs, ebvs, vsinis, d_b, d_x, total, logg, mh,
              clamped, t_hat, e_hat, v_hat, slope, dust_cost, vhat):
    fig, axes = plt.subplots(1, 3, figsize=(15.0, 5.1), constrained_layout=True)
    fig.patch.set_facecolor(SURFACE)
    XE, YE = np.meshgrid(teffs, ebvs, indexing='ij')
    XV, YV = np.meshgrid(teffs, vsinis, indexing='ij')

    for ax, (name, d, x, y) in zip(axes, (
            ('ngsl_bands', d_b, XE, YE), ('xsl', d_x, XV, YV),
            ('combined', total, XE, YE))):
        style(ax)
        im = heat(ax, x, y, d)
        i, j = np.unravel_index(np.nanargmin(d), d.shape)
        yv = (vsinis if name == 'xsl' else ebvs)[j]
        lab = (f'min: {teffs[i]:.0f} K, {yv:.0f} km/s' if name == 'xsl'
               else f'min: {teffs[i]:.0f} K, E(B−V)={yv:.3f}')
        ax.plot(teffs[i], yv, marker='*', ms=15, color='w',
                markeredgecolor='k', markeredgewidth=.6, zorder=5, label=lab)
        ax.set_title({'ngsl_bands': 'NGSL bands — continuum (dust)',
                      'xsl': 'XSL lines — dust-immune Teff + v sin i',
                      'combined': 'combined (each leg rescaled)'}[name],
                     fontsize=10, color=INK)
        ax.set_xlabel('Teff [K]', fontsize=9, color=INK)

        if name == 'xsl':
            # the ridge of best-fit v sin i, i.e. what is profiled away
            ok = np.isfinite(vhat)
            ax.plot(teffs[ok], vhat[ok], color='#7fd4ff', lw=1.0, ls='--',
                    alpha=.8, label='v sin i profile')
            ax.set_ylabel('v sin i [km s$^{-1}$]', fontsize=9, color=INK)
        ax.legend(fontsize=7.5, loc='upper left', framealpha=.85)

        for key, lab in (('teff_ngsl', 'NGSL Teff'), ('teff_xsl', 'XSL Teff')):
            v = opt_float(row.get(key))
            if v and teffs.min() <= v <= teffs.max():
                ax.axvline(v, color=CAT_C[key], ls='--', lw=1.2, alpha=.9)
                ax.text(v, ax.get_ylim()[1], f' {lab}', rotation=90, va='top',
                        ha='right', fontsize=7, color=CAT_C[key])
        if name == 'xsl':
            continue
        for key, lab in (('ebv_phot', 'photometric'), ('ebv_miles', 'MILES'),
                         ('ebv_sf11', 'SF11 (upper bound)')):
            v = opt_float(row.get(key))
            if v is None:
                continue
            if key in ('ebv_sf11', 'ebv_sfd98') and row.get('ebv_map_useful') == 'no':
                ax.annotate(f'{lab} {v:.2f} — map unusable at b={row["gal_b"]}°',
                            xy=(teffs[0], ebvs[-1]), xytext=(3, -3),
                            textcoords='offset points', fontsize=7,
                            color=CAT_C[key], va='top')
                continue
            if ebvs.min() <= v <= ebvs.max():
                ax.axhline(v, color=CAT_C[key], ls=':', lw=1.2, alpha=.9)
                ax.text(teffs[0], v, f' {lab}', va='bottom', fontsize=7,
                        color=CAT_C[key])
            else:
                # a value off the scale is stated, not silently dropped: a
                # negative photometric E(B-V) means "consistent with zero"
                side = 'below' if v < ebvs.min() else 'above'
                yy = ebvs.min() if v < ebvs.min() else ebvs.max()
                ax.annotate(f'{lab} {v:+.3f} (off scale {side})',
                            xy=(teffs[0], yy), xytext=(3, 3 if side == 'below' else -12),
                            textcoords='offset points', fontsize=7,
                            color=CAT_C[key], va='bottom')

    axes[0].set_ylabel('E(B−V) [mag]', fontsize=9, color=INK)
    axes[2].set_ylabel('E(B−V) [mag]', fontsize=9, color=INK)
    fig.colorbar(im, ax=axes.tolist(), shrink=.85, label='Δχ² (log scale)')

    clamp = ''.join(
        f'  —  {"[M/H]" if k == "mh" else "log g"} held at the grid edge '
        f'{mh if k == "mh" else logg:+.2f}, catalog {v:+.2f}'
        for k, v in clamped.items() if v is not None)
    fig.suptitle(
        f'{star} ({row["tier"]}) — χ² surfaces   '
        f'(log g={logg:.2f}, [M/H]={mh:+.2f} held){clamp}\n'
        f'ridge {slope:+.0f} K per 0.01 mag;  XSL dust sensitivity discarded '
        f'Δχ²={dust_cost:+.2f} over the full E(B−V) range;  '
        'contours are nominal 2-parameter levels, NOT calibrated — residuals '
        'are correlated',
        fontsize=10.5, color=INK, linespacing=1.5)
    out = figure_path('ebv_teff', star)
    fig.savefig(out, dpi=165, facecolor=SURFACE)
    plt.close(fig)
    print(f'  -> {out.relative_to(ROOT)}')


def draw_summary(rows):
    """The sample question: is the E(B-V) offset common to all the stars?

    A common offset points at NGSL's calibration or at the extinction law; a
    scatter points at the individual stars. That is the whole reason for running
    the sweep on more than one star, so the left panel is deliberately
    fitted-against-catalog with the 1:1 line, not two separate histograms.

    Stars whose solution sits on a boundary are drawn as LIMITS (arrows) and
    excluded from the mean offset. Averaging a value that the scan was not free
    to move would bias the one number the panel exists to produce.
    """
    fig, axes = plt.subplots(1, 3, figsize=(15.4, 5.2), constrained_layout=True)
    fig.patch.set_facecolor(SURFACE)
    for ax in axes:
        style(ax)
        ax.grid(alpha=.25, color=GRIDC, lw=.7)
    mark = {'primary': 'o', 'secondary': 's'}

    def col(r, key='ebv_phot'):
        return CAT_C[key] if r['tier'] == 'primary' else MUTED

    def hit(r, what):
        return what in (r.get('at_boundary') or '')

    # 1. fitted vs photometric E(B-V) --------------------------------------
    ax = axes[0]
    lim = [-0.03, 0.45]
    ax.plot(lim, lim, color=MUTED, lw=1.0, ls='--', alpha=.7, label='1:1')
    clean = []
    for r in rows:
        pv = opt_float(r['ebv_phot'])
        if pv is None:
            continue
        e, lo, hi = opt_float(r['ebv']), opt_float(r['ebv_lo']), opt_float(r['ebv_hi'])
        lim_hi, lim_lo = hit(r, 'ebv_ceiling'), hit(r, 'ebv_floor')
        if below_grid(r['star']):
            # below the grid floor: plotted, but never in the mean offset
            ax.plot(pv, e, marker='x', ms=7, color='#c0392b', mew=1.6)
            continue
        if lim_hi or lim_lo:
            ax.plot(pv, e, marker='^' if lim_hi else 'v', ms=7, color=col(r),
                    mfc='none', mew=1.4)
        else:
            ax.errorbar(pv, e, yerr=[[max(e - lo, 0)], [max(hi - e, 0)]],
                        marker=mark[r['tier']], ms=6, lw=0, elinewidth=1.1,
                        color=col(r), capsize=2)
            if r['tier'] == 'primary':
                clean.append(e - pv)
        ax.annotate(r['star'].replace('HD', '')
                    + ('*' if below_grid(r['star']) else ''),
                    (pv, e), xytext=(4, 3), textcoords='offset points',
                    fontsize=6.5, color=MUTED)
    if clean:
        m, sd = np.mean(clean), np.std(clean)
        ax.axline((0, m), slope=1, color=CAT_C['ebv_sf11'], lw=1.3, alpha=.85,
                  label=(f'primary, off-boundary (n={len(clean)}):\n'
                         f'mean {m:+.3f} ± {sd / np.sqrt(len(clean)):.3f}, '
                         f'scatter {sd:.3f}'))
    ax.plot([], [], marker='^', lw=0, mfc='none', color=MUTED, mew=1.3,
            label='on a boundary — a limit, not a value')
    ax.plot([], [], marker='x', lw=0, color='#c0392b', mew=1.6,
            label='catalog [M/H] < −0.5 — excluded from the mean')
    ax.set_xlim(lim)
    ax.set_ylim(-0.01, 0.32)
    ax.set_xlabel('catalog photometric E(B−V)', fontsize=9, color=INK)
    ax.set_ylabel('fitted E(B−V)', fontsize=9, color=INK)
    ax.set_title('Is the dust offset common to the sample?', fontsize=10,
                 color=INK)
    ax.legend(fontsize=7, loc='upper left', framealpha=.9)

    # 2. fitted vs catalog Teff --------------------------------------------
    ax = axes[1]
    tl = [7800, 11800]
    ax.plot(tl, tl, color=MUTED, lw=1.0, ls='--', alpha=.7, label='1:1')
    for key, lab in (('teff_ngsl', 'NGSL catalog'), ('teff_xsl', 'XSL catalog')):
        xs, ys = [], []
        for r in rows:
            c, t = opt_float(r[key]), opt_float(r['teff'])
            if c is None:
                continue
            xs.append(c)
            ys.append(t)
            ax.errorbar(c, t, yerr=[[max(t - opt_float(r['teff_lo']), 0)],
                                    [max(opt_float(r['teff_hi']) - t, 0)]],
                        marker=mark[r['tier']], ms=5.5, lw=0, elinewidth=1.0,
                        color=CAT_C[key], capsize=2, alpha=.9,
                        mfc='none' if hit(r, 'teff_edge') else CAT_C[key])
        if xs:
            d = np.array(ys) - np.array(xs)
            ax.plot([], [], marker='o', lw=0, color=CAT_C[key],
                    label=f'{lab}: median Δ {np.median(d):+.0f} K, '
                          f'scatter {np.std(d):.0f} K')
    for r in rows:
        ax.annotate(r['star'].replace('HD', ''),
                    (opt_float(r['teff_xsl']), opt_float(r['teff'])), xytext=(4, 3),
                    textcoords='offset points', fontsize=6.5, color=MUTED)
    for e in (8500, 11500):
        ax.axhline(e, color=MUTED, lw=.8, ls=':', alpha=.7)
    ax.text(tl[0] + 60, 8500, ' grid edge', fontsize=6.5, color=MUTED,
            va='bottom')
    ax.set_xlim(tl)
    ax.set_ylim(8300, 11700)
    ax.set_xlabel('catalog Teff [K]', fontsize=9, color=INK)
    ax.set_ylabel('fitted Teff [K]  (grid nodes)', fontsize=9, color=INK)
    ax.set_title('Teff: the fit against the two catalogs', fontsize=10,
                 color=INK)
    ax.legend(fontsize=7, loc='upper left', framealpha=.9)

    # 3. v sin i, which nobody has published for these stars ---------------
    # explore/vsini_mask_test.py refits each star at 4 core-mask widths. A real
    # rotation cannot depend on where the mask edge is, so the ones that move
    # are drawn open: they are measuring the NLTE core mismatch, not the star.
    ax = axes[2]
    stab = {}
    mt = ROOT / 'data' / 'vsini_mask_test.csv'
    if mt.exists():
        stab = {r['star']: r for r in csv.DictReader(open(mt))}
    order = sorted(rows, key=lambda r: opt_float(r['vsini']))
    for k, r in enumerate(order):
        v, lo, hi = opt_float(r['vsini']), opt_float(r['vsini_lo']), opt_float(r['vsini_hi'])
        st = stab.get(r['star'])
        bad = st is not None and st['vsini_stable'] == 'no'
        if bad:      # the full range the mask width moves it over
            vv = [opt_float(st[c]) for c in st if c.startswith('vsini_mask')]
            ax.plot([min(vv), max(vv)], [k, k], color=CAT_C['ebv_sf11'], lw=3,
                    alpha=.35, solid_capstyle='round', zorder=1)
        else:
            ax.plot([lo, hi], [k, k], color=GRIDC, lw=3, solid_capstyle='round',
                    zorder=1)
        ax.plot(v, k, marker=mark[r['tier']], ms=6.5, zorder=3,
                color=CAT_C['ebv_sf11'] if bad else
                (CAT_C['teff_xsl'] if r['tier'] == 'primary' else MUTED),
                mfc='none' if (bad or hit(r, 'vsini_ceiling')) else None,
                mew=1.5)
        note = ' '.join(x for x in (
            'mask-driven' if bad else '',
            '[M/H] clamped' if r['clamped'] else '') if x)
        if note:
            ax.annotate(note, (VSINI_GRID[-1], k), xytext=(-4, 0),
                        textcoords='offset points', ha='right', va='center',
                        fontsize=6.5, color=CAT_C['ebv_sf11'])
    if stab:
        nb = sum(1 for r in rows if stab.get(r['star'], {}).get('vsini_stable') == 'no')
        ax.plot([], [], marker='o', lw=0, mfc='none', mew=1.5,
                color=CAT_C['ebv_sf11'],
                label=f'{nb} move with the core mask — not measurements')
        ax.legend(fontsize=7, loc='lower right', framealpha=.9)
    ax.axvline(15, color=CAT_C['ebv_sf11'], ls=':', lw=1.2)
    ax.text(15, len(order) - .5, ' XSL resolution limit ', fontsize=7,
            color=CAT_C['ebv_sf11'], va='top', rotation=90, ha='right')
    ax.set_yticks(np.arange(len(order)))
    ax.set_yticklabels([r['star'] for r in order], fontsize=7.5, color=MUTED)
    ax.set_ylim(-.7, len(order) - .3)
    ax.set_xlim(-12, VSINI_GRID[-1] + 12)
    ax.set_xlabel('v sin i [km s$^{-1}$]  (nominal Δχ² ≤ 1 at the best Teff)',
                  fontsize=9, color=INK)
    ax.set_title('v sin i from the XSL lines', fontsize=10, color=INK)
    

    n_p = sum(1 for r in rows if r['tier'] == 'primary')
    fig.suptitle(
        f'E(B−V)–Teff sweep over the sample — {len(rows)} stars '
        f'({n_p} primary, circles; {len(rows) - n_p} secondary, squares); '
        '* = catalog [M/H] < −0.5, below the grid floor\n'
        'log g and [M/H] held at the nearest node to the NGSL catalog value; '
        'error bars are nominal Δχ² ≤ 1 profiles, NOT calibrated — the '
        'star-to-star scatter is the honest uncertainty',
        fontsize=10.5, color=INK, linespacing=1.5)
    out = ROOT / 'figures' / 'ebv_teff_sample.png'
    fig.savefig(out, dpi=165, facecolor=SURFACE)
    plt.close(fig)
    print(f'-> {out.relative_to(ROOT)}')


def report(rows):
    """Print the sample table and the question the sweep exists to answer."""
    hdr = (f'{"star":<10}{"tier":<6}{"Teff [lo,hi]":>20}'
           f'{"E(B-V) [lo,hi]":>23}{"vsini [lo,hi]":>17}'
           f'{"E_phot":>8}{"dE":>8}{"x2/n_b":>8}{"x2/n_x":>8}  flags')
    print('\n' + hdr)
    print('-' * len(hdr))
    for r in rows:
        pv = opt_float(r['ebv_phot'])
        d = opt_float(r['ebv']) - pv if pv is not None else np.nan
        flags = ' '.join(x for x in (r['at_boundary'].replace('+', ' '),
                                     ('[M/H]@edge' if r['clamped'] else '')) if x)
        print(f'{r["star"]:<10}{r["tier"][:4]:<6}'
              + f'{float(r["teff"]):>6.0f} [{float(r["teff_lo"]):.0f},{float(r["teff_hi"]):.0f}]'.rjust(20)
              + f'{float(r["ebv"]):>6.3f} [{float(r["ebv_lo"]):.3f},{float(r["ebv_hi"]):.3f}]'.rjust(23)
              + f'{float(r["vsini"]):>4.0f} [{float(r["vsini_lo"]):.0f},{float(r["vsini_hi"]):.0f}]'.rjust(17)
              + (f'{pv:>8.3f}' if pv is not None else f'{"-":>8}')
              + (f'{d:>+8.3f}' if np.isfinite(d) else f'{"-":>8}')
              + f'{float(r["chi2_bands"]) / max(int(r["n_bands"]), 1):>8.2f}'
              + f'{float(r["chi2_xsl"]) / max(int(r["n_xsl"]), 1):>8.2f}'
              + (f'  {flags}' if flags else ''))
    print('-' * len(hdr))
    # Stars whose catalog [M/H] is below the grid's -0.5 floor are reported
    # separately and never pooled: no node can represent them, so the fit pays
    # for the mismatch somewhere else and their E(B-V) is not measuring the
    # same thing. See common/figpath.py.
    for grp, keep in (('inside the grid', lambda r: not below_grid(r['star'])),
                      ('BELOW the grid [M/H] floor', lambda r: below_grid(r['star']))):
        for tier in ('primary', 'secondary'):
            sub = [r for r in rows if r['tier'] == tier and keep(r)
                   and not (r['at_boundary'] or '')
                   and opt_float(r['ebv_phot']) is not None]
            d = [opt_float(r['ebv']) - opt_float(r['ebv_phot']) for r in sub]
            if not d:
                continue
            print(f'  {grp}, {tier}, off-boundary: fitted - photometric E(B-V) '
                  f'= {np.mean(d):+.3f} +/- {np.std(d) / np.sqrt(len(d)):.3f} '
                  f'(sem), scatter {np.std(d):.3f}, n={len(d)}')
    for what, msg in (('ebv_ceiling', 'E(B-V) ceiling'),
                      ('ebv_floor', 'E(B-V) floor (= consistent with zero)'),
                      ('teff_edge', 'Teff grid edge'),
                      ('vsini_ceiling', 'v sin i ceiling')):
        bad = [r['star'] for r in rows if what in (r['at_boundary'] or '')]
        if bad:
            print(f'  AT THE {msg}: {", ".join(bad)} -- limits, not measurements')
    cl = [r['star'] for r in rows if r['clamped']]
    if cl:
        print(f'  [M/H] HELD AT THE GRID EDGE: {", ".join(cl)}')
    # SF11 is the integrated Galactic column to infinity. A star inside the
    # Galaxy cannot be reddened by more than that, so exceeding it is not a
    # tension to weigh against other evidence -- it is unphysical, and marks a
    # solution that has gone wrong somewhere else.
    over = [(r['star'], opt_float(r['ebv']), opt_float(r['ebv_sf11']), bool(r['clamped']))
            for r in rows
            if r.get('ebv_map_useful') != 'no' and opt_float(r['ebv_sf11'])
            and opt_float(r['ebv']) > opt_float(r['ebv_sf11'])]
    if over:
        print('  ABOVE THE SF11 TOTAL GALACTIC COLUMN (unphysical):')
        for st, e, sf, cla in over:
            print(f'    {st:<10} fitted {e:.3f} vs column {sf:.3f} '
                  f'({e / sf:.1f}x)' + ('   [M/H] clamped' if cla else ''))
        ncl = sum(1 for *_, c in over if c)
        print(f'    {ncl} of {len(over)} are [M/H]-clamped -- the grid floor is '
              f'the first thing to suspect')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--star', default='HD194453')
    ap.add_argument('--all', action='store_true',
                    help='every primary + secondary star, then the summary')
    ap.add_argument('--summary-only', action='store_true',
                    help=f'redraw the summary from {SCAN_CSV.name}')
    ap.add_argument('--vsini-bands', type=float, default=0.0,
                    help='v sin i for the bands leg; they are insensitive to it')
    ap.add_argument('--ebv-max', type=float, default=0.30,
                    help='must clear the SF11 upper bounds (to 0.28 here) or '
                         'reddened stars pile up on the ceiling')
    ap.add_argument('--nebv', type=int, default=121)
    a = ap.parse_args()

    if a.summary_only:
        rows = list(csv.DictReader(open(SCAN_CSV)))
        report(rows)
        draw_summary(rows)
        return

    grid = Grid()
    if not a.all:
        r = run_star(a.star, grid, a)
        if r:
            report([r])
        return

    rows = []
    for s in sample_rows():
        try:
            r = run_star(s['star'], grid, a)
        except Exception as exc:                 # a missing file must not stop
            print(f'{s["star"]}: FAILED {type(exc).__name__}: {exc}')
            continue
        if r:
            rows.append(r)
        print()
    if not rows:
        print('no stars succeeded')
        return
    with open(SCAN_CSV, 'w', newline='') as fh:
        wr = csv.DictWriter(fh, fieldnames=list(rows[0]))
        wr.writeheader()
        wr.writerows(rows)
    print(f'-> {SCAN_CSV.relative_to(ROOT)}  ({len(rows)} stars)')
    report(rows)
    draw_summary(rows)


if __name__ == '__main__':
    main()
