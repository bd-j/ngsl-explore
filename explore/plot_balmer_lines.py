"""One panel per Balmer member, model against XSL, with its residual.

The companion to `explore/plot_metal_lines.py`: that figure asks whether the
grid can reach the observed metal-line depths, this one asks whether the models
reproduce the hydrogen PROFILES -- which is a different question with a
different answer, and the one the whole project turns on.

Each member gets two panels:

  * the window, XSL against the model, drawn with the house convention (solid
    where the pixels set the model's scale, faded where they did not, the
    masked core tinted);
  * the fractional residual underneath, with the SAME split, so a residual on
    a pixel that helped fit the scale cannot be read as a prediction.

and the figure ends with the series summary -- core residual against upper
level n -- which is the thing no single panel can show.

WHAT IS A PREDICTION HERE AND WHAT IS NOT

  H-alpha .. H-delta   FITTED. XSL conditions on these four. Their panels show
                       the fit's own calibrated model, so what is drawn is what
                       the likelihood saw. The masked CORE is still predicted.
  H-epsilon and up     Fitted by nothing, in either dataset: XSL excludes them
                       (they blend, so no local continuum is defined) and NGSL
                       holds the whole 3550-4000 A window out. Only an order-1
                       scale is solved on each one's wings, here, for display.

So every core in this figure is a prediction, and from H-epsilon up so is
every wing. See common/balmer.py for the window and core-mask definitions.

WHY IT IS WORTH A FIGURE. The held-out NGSL Balmer residual carries a line-core
excess, and at NGSL's 1.4 A pixels that number is degenerate with the
instrument profile -- it runs +2.65% to -4.20% at HD194453 as the Moffat core
goes 3.54 -> 7.00 A (docs/LSF.md). XSL resolves the same lines ~16x better,
where that sensitivity is ~10x smaller, and it resolves the whole series at
once rather than one line at a time.

    python3 explore/plot_balmer_lines.py --star HD194453
    python3 explore/plot_balmer_lines.py --all
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
from fitting.observations import load_xsl
from fitting.predict import predict
from fitting.calibration import solve
from common.balmer import (member_panel, core_stat, wing_stat, FITTED_N, NMAX,
                           MIN_WING_PX)
from common.lines import BALMER, ISM_LINES, balmer_member
from common.specplot import (spectrum_panel, residual_panel, style, MOD_C,
                             HELD_C, SURFACE, INK, MUTED)
from common.figpath import figure_path
from common.sample import sample_row

ROOT = Path(__file__).resolve().parent.parent

# Fixed residual range, as the break panels use, so that the panels are
# comparable BETWEEN members and between stars. Per-panel autoscaling would
# make a 1% core residual and a 10% one look identical, which is the one
# comparison this figure exists to support.
RLIM = (-10., 10.)

# Members drawn, blue-ward first so the figure reads down the series the way
# the break does. Anything unusable for a given star drops out in `run`.
MEMBERS = tuple(range(3, NMAX + 1))

# Greek names for the four the literature names; the rest are Hn.
NAMES = {3: 'Hα', 4: 'Hβ', 5: 'Hγ', 6: 'Hδ', 7: 'Hε'}


def member_name(n):
    return NAMES.get(n, f'H{n}')


def ism_marks(lo, hi):
    """Interstellar lines inside a window, as (lambda, label) markers.

    H-epsilon is 1.6 A from Ca II H and its window also holds Ca II K, so its
    'core' is a blend of a stellar hydrogen line with two interstellar ones.
    That is not a detail -- it is why H-epsilon is excluded from the fit -- and
    a panel that does not say so invites reading the blend as a model error.
    """
    return [(lam, name) for name, lam in ISM_LINES.items() if lo < lam < hi]


def run(star, a):
    """One star: solve, report the series, and draw its figure."""
    row, grid = sample_row(star), Grid()
    node = None
    if not a.no_scan and (ROOT / 'results' / star / 'scan.npz').exists():
        from fitting.scan import best_node
        node = best_node(star)
    ebv = a.ebv if a.ebv is not None else (node['ebv'] if node else 0.0)
    vsini = a.vsini if a.vsini is not None else (node['vsini'] if node else 0.0)

    if node is not None:
        teff, logg, mh = node['teff'], node['logg'], node['mh']
        src = 'node-scan ML'
    else:
        teff = float(grid.teff[np.argmin(np.abs(grid.teff - float(row['teff_ngsl'])))])
        logg = float(grid.logg[np.argmin(np.abs(grid.logg - float(row['logg_ngsl'])))])
        mh = float(grid.mh[np.argmin(np.abs(grid.mh - float(row['mh_ngsl'])))])
        src = 'nearest node to catalog'

    # the fit's OWN view of XSL -- default windows, bad regions dropped -- so
    # the four fitted members are drawn exactly as the likelihood saw them
    obs = load_xsl(star)
    pv = predict(dict(teff=teff, logg=logg, mh=mh, ebv=ebv, vsini=vsini),
                 [obs], grid)[0].value
    fit_cal, _ = solve(obs, pv, fill_domains=True)

    print(f'{star}: {src} Teff={teff:.0f} log g={logg:.2f} [M/H]={mh:+.2f}, '
          f'E(B-V)={ebv:.3f}, v sin i={vsini:.0f}')
    print(f'  XSL n={obs.ndata}; fitted members {[member_name(n) for n in FITTED_N]}')

    panels = []
    for n in MEMBERS:
        d = member_panel(obs, pv, n, fit_cal=fit_cal)
        if d is not None:
            panels.append(d)
    dropped = [n for n in MEMBERS if n not in {d['n'] for d in panels}]
    if dropped:
        # Said out loud rather than silently skipped: at the top of the series
        # the wings run out (H13 leaves ~17 usable pixels against ~95 at H11),
        # so there is nothing left to solve an order-1 scale on.
        print(f'  not drawn (fewer than {MIN_WING_PX} wing px, or no data): '
              + ', '.join(member_name(n) for n in dropped))

    # 'mask' is what the local scale was kept away from; 'over' is the
    # narrower region the core residual is actually measured on. See
    # common/balmer.py for why they are not the same number.
    print(f'\n  {"line":>6}{"lambda":>10}{"role":>12}{"mask":>7}{"over":>7}'
          f'{"core":>9}{"rms":>8}{"n":>6}{"wings":>8}')
    stats = {}
    for d in panels:
        cs, ws = core_stat(d), wing_stat(d)
        if cs is None:
            continue
        stats[d['n']] = cs
        role = 'FITTED' if d['source'] == 'fit' else 'predicted'
        print(f'  {member_name(d["n"]):>6}{d["lam"]:>10.2f}{role:>12}'
              f'{d["core"]:>7.1f}{d["core_hw"]:>7.1f}{cs[0]:>+9.2f}'
              f'{cs[1]:>8.2f}{cs[2]:>6d}{ws:>+8.2f}')
    # One star is one star. The trend across the series is what the summary
    # panel is for; whether it is a property of THIS sample needs --all and the
    # twelve of them side by side, not this line of output.
    print(f'\n  {len(stats)} members measured; core residual runs '
          f'{min(v[0] for v in stats.values()):+.2f}% to '
          f'{max(v[0] for v in stats.values()):+.2f}%')

    figure(star, panels, stats, dict(teff=teff, logg=logg, mh=mh, ebv=ebv,
                                     vsini=vsini, src=src, node=node))
    return dict(star=star, stats=stats)


def figure(star, panels, stats, p):
    """Two columns of (window, residual) pairs, then the series summary."""
    ncol = 2
    nrow = int(np.ceil(len(panels) / ncol))
    fig = plt.figure(figsize=(13.0, 2.9 * nrow + 3.4))
    fig.patch.set_facecolor(SURFACE)
    # Nested, not one flat grid of 2*nrow rows. A residual belongs hard against
    # the window it describes while one member must stand clear of the next, and
    # a single hspace cannot do both: flat, every panel's title landed on the
    # wavelength axis of the pair above it.
    outer = fig.add_gridspec(nrow + 1, ncol, height_ratios=[2.3] * nrow + [1.6],
                             hspace=.40, wspace=.19)

    last_in_col = {}                 # so only the bottom pair carries an axis
    for i, d in enumerate(panels):
        last_in_col[i % ncol] = i
    for i, d in enumerate(panels):
        r, c = divmod(i, ncol)
        inner = outer[r, c].subgridspec(2, 1, height_ratios=[1.6, 0.8],
                                        hspace=.07)
        lo, hi = d['window']
        cs = stats.get(d['n'])
        fitted = d['source'] == 'fit'
        role = ('FITTED (±50 Å, core masked)' if fitted
                else 'PREDICTED — nothing fits this line')
        marks = ism_marks(lo, hi)
        head = f'{member_name(d["n"])}  {d["lam"]:.1f} $\\AA$ — {role}'
        if cs is not None:
            head += f';  core {cs[0]:+.2f}%'
        if marks:
            head += '   ⚠ ' + ' + '.join(m[1] for m in marks) + ' in window'

        ax = fig.add_subplot(inner[0])
        spectrum_panel(ax, d['view'], [('model', d['cal'])], lo, hi,
                       obs_label=('XSL (solid = fitted)' if fitted
                                  else 'XSL (solid = sets the local scale)'),
                       legend=(i == 0), fontsize=9)
        ax.set_title(head, fontsize=9, color=(HELD_C if marks else INK))
        ax.set_xlabel('')
        ax.tick_params(labelbottom=False)
        ax.axvline(d['lam'], color=MUTED, ls='--', lw=1)

        axr = fig.add_subplot(inner[1], sharex=ax)
        residual_panel(axr, d['view'].wavelength,
                       [(d['resid'], d['used'], d['pred'], [d['incore']])],
                       lo, hi,
                       ylim=RLIM, marks=marks, ylabel='data−model [%]',
                       xlabel=(i == last_in_col[c]),
                       bands=[(d['lam'] - d['core'], d['lam'] + d['core'])])
        axr.axvline(d['lam'], color=MUTED, ls='--', lw=1)

    # --- the series summary --------------------------------------------
    #
    # The panels above show one line each; the question the figure is really
    # asking -- does the core residual run with the order of the line, as an
    # NLTE or a Stark-broadening error would, or is it flat, as a calibration
    # error would -- only exists across the series.
    axs = fig.add_subplot(outer[-1, :])
    style(axs)
    ns = sorted(stats)
    med = [stats[n][0] for n in ns]
    # The bar is the SCATTER inside the core, not an error on the median: XSL
    # pixels are correlated over its LSF (~4 px) and by the local scale solved
    # on the wings, so n_px is not a count of independent measurements and
    # rms/sqrt(n) would be an uncertainty this figure has not earned.
    rms = [stats[n][1] for n in ns]
    fitted = [n in FITTED_N for n in ns]
    axs.errorbar(ns, med, yerr=rms, fmt='none', ecolor=MOD_C, alpha=.45,
                 capsize=3, lw=1)
    axs.plot(ns, med, '-', color=MOD_C, lw=1.2, zorder=2)
    axs.plot([n for n, f in zip(ns, fitted) if f],
             [m for m, f in zip(med, fitted) if f], 'o', color=MOD_C, ms=7,
             mfc=SURFACE, mew=1.6, zorder=3, label='fitted line (core still predicted)')
    axs.plot([n for n, f in zip(ns, fitted) if not f],
             [m for m, f in zip(med, fitted) if not f], 'o', color=MOD_C, ms=7,
             zorder=3, label='predicted line (nothing fits it)')
    axs.axhline(0, color=MUTED, lw=1)
    axs.set_xticks(ns)
    axs.set_xticklabels([member_name(n) for n in ns])
    axs.set_xlabel('Balmer member (upper level n)', fontsize=9, color=INK)
    axs.set_ylabel('core residual, median [%]', fontsize=9, color=INK)
    axs.set_title('The series — core (data − model)/model against line order;  '
                  'bars are the SCATTER within the core, not an error on the '
                  'median', fontsize=9, color=INK)
    axs.legend(fontsize=8, framealpha=.92)

    node = p['node']
    line2 = (f'{p["src"]}: Teff={p["teff"]:.0f} / log g={p["logg"]:.2f} / '
             f'[M/H]={p["mh"]:+.2f} / E(B−V)={p["ebv"]:.3f} / '
             f'v sin i={p["vsini"]:.0f} km/s'
             + ('   ⚠ [M/H] AT THE GRID FLOOR'
                if node and node['at_mh_floor'] else ''))
    fig.suptitle(
        f'{star} — the Balmer series in XSL\n{line2}\n'
        'every core is a prediction; from Hε up so is every wing — '
        'shaded = masked core, faded = pixels that set the scale\n'
        'residuals are (DATA − MODEL)/model, so positive means the '
        'observation is brighter than the model',
        fontsize=11, color=INK, linespacing=1.5)
    out = figure_path('balmer_lines', star)
    fig.savefig(out, dpi=170, facecolor=SURFACE, bbox_inches='tight')
    plt.close(fig)
    print(f'  -> {out.relative_to(ROOT)}')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--star', default='HD194453')
    ap.add_argument('--all', action='store_true',
                    help='every primary + secondary star in the sample')
    # None, not 0: the default is the node scan's own solution when it exists.
    ap.add_argument('--vsini', type=float, default=None)
    ap.add_argument('--ebv', type=float, default=None)
    ap.add_argument('--no-scan', action='store_true')
    a = ap.parse_args()
    stars = ([r['star'] for r in csv.DictReader(open(ROOT / 'data' / 'sample.csv'))
              if r['tier'] in ('primary', 'secondary')] if a.all else [a.star])
    out = []
    for st in stars:
        try:
            out.append(run(st, a))
        except Exception as exc:
            print(f'{st}: FAILED {type(exc).__name__}: {exc}')
        print()
    return out


if __name__ == '__main__':
    main()
