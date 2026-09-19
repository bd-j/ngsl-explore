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

The figure also carries an H11 + H10 row at XSL resolution. Those two members
lie inside the held-out NGSL window and are fitted by nothing in either
dataset, so they are predictions on the same footing as the breaks -- and they
are the one view of the Balmer core excess that is not dominated by the NGSL
instrument profile. The whole series is in explore/plot_balmer_lines.py.

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
from common.specplot import (spectrum_panel, residual_panel, offscale_note,
                             style, OBS_C, MOD_C, MOD2_C, BAND_C, HELD_C,
                             SURFACE, INK, MUTED)
from common.lines import hydrogen_lines
from common.balmer import member_panel, merge_members, core_stat, wing_stat

from common.figpath import figure_path
from common.sample import sample_row

ROOT = Path(__file__).resolve().parent.parent
CORE_HALF_A = 4.0        # Balmer 'core' half-width, ~1.5 G430L pixels
                         # -- the same definition explore/ngsl_lsf.py uses
BALMER, PASCHEN = 3646.0, 8205.9

# The two high-order members shown at XSL resolution. Both sit INSIDE the
# held-out NGSL window (3550-4000 A) and neither is fitted by anything -- XSL
# conditions on H-alpha through H-delta and the metal windows only -- so every
# pixel of that panel is a prediction from both datasets, exactly like the
# break panels above it.
#
# They are here because the NGSL core excess cannot be interpreted on its own.
# It is degenerate with the instrument profile: at HD194453 it runs +2.65% at a
# 3.54 A Moffat core to -4.20% at 7.00 A (docs/LSF.md). XSL resolves the same
# lines with a 0.39 A FWHM at H10, ~10x narrower than the NGSL core, so the
# same excess measured there is ~10x less sensitive to the LSF -- which is what
# makes it an independent check rather than a restatement.
#
# NOTE the two numbers are not the same quantity and must not be read as one.
# core_excess() below is measured against the NGSL window's own continuum under
# a single scalar from the bands; the XSL core residual is measured against an
# order-1 scale solved on that line's own wings, which absorbs a continuum
# error the NGSL number still carries. They agree on the SIGN and shape of a
# core excess, not on its size.
H10_H11_N = (11, 10)
H10_H11_RLIM = (-15., 15.)      # the Balmer break panel's range, deliberately


def nearest_node(grid, teff, logg, mh):
    return (float(grid.teff[np.argmin(np.abs(grid.teff - teff))]),
            float(grid.logg[np.argmin(np.abs(grid.logg - logg))]),
            float(grid.mh[np.argmin(np.abs(grid.mh - mh))]))


def scan_row(star):
    """The star's row from the E(B-V)-Teff sweep, if it has been run."""
    p = ROOT / 'data' / 'ebv_teff_scan.csv'
    if not p.exists():
        return None
    for r in csv.DictReader(open(p)):
        if r['star'] == star:
            return r
    return None


def scan_params(star):
    """Maximum-likelihood node from the full node scan, or None.

    This is the preferred source: unlike data/ebv_teff_scan.csv it has Teff,
    log g, [M/H], E(B-V) and v sin i all determined TOGETHER, so the figure
    shows a single self-consistent model rather than one parameter borrowed
    from a fit that held the others fixed.
    """
    if not (ROOT / 'results' / star / 'scan.npz').exists():
        return None
    from fitting.scan import best_node
    return best_node(star)


def run(star, a):
    """One star: fit, report, and draw its figure."""
    ebv, vsini = a.ebv, a.vsini
    node = None if a.no_scan else scan_params(star)
    sc = None if (a.no_scan or node is not None) else scan_row(star)
    if node is not None:
        src = (f'node scan ML: Teff={node["teff"]:.0f} log g={node["logg"]:.2f} '
               f'[M/H]={node["mh"]:+.2f}')
        ebv = node['ebv'] if ebv is None else ebv
        vsini = node['vsini'] if vsini is None else vsini
    elif sc is not None:
        src = f'ebv_teff_scan.csv (fitted at Teff={float(sc["teff"]):.0f} K)'
        ebv = float(sc['ebv']) if ebv is None else ebv
        vsini = float(sc['vsini']) if vsini is None else vsini
    else:
        src = 'command line'
    ebv = 0.0 if ebv is None else ebv
    vsini = 0.0 if vsini is None else vsini

    row, grid = sample_row(star), Grid()
    bands = ngsl_band_edges()
    print(f'{star}: {len(bands)} NGSL bands, '
          f'{bands[0][1]:.0f}-{bands[-1][2]:.0f} A')
    print(f'  E(B-V)={ebv:.3f}, v sin i={vsini:.0f} km/s  <- {src}')
    print(f'  held out: Balmer {BREAK_WINDOW}, Paschen {PASCHEN_WINDOW}')

    cond = conditioning_set(star)
    nb = cond[0]
    xs = next((o for o in cond if o.name == 'xsl'), None)
    spec = load_ngsl(star)
    held = {'Balmer': heldout(star, window=BREAK_WINDOW),
            'Paschen': heldout(star, window=PASCHEN_WINDOW)}
    print('  conditioning: ' + ', '.join(f'{o.name}(n={o.ndata})' for o in cond))
    print('  held out    : ' + ', '.join(f'{k}(n={v.ndata})'
                                         for k, v in held.items()))

    t0 = (float(row['teff_ngsl']), float(row['logg_ngsl']), float(row['mh_ngsl']))
    # The stellar parameters must come from the SAME fit as the dust. They sit
    # on a common degeneracy ridge (+93 K per 0.01 mag), so pairing one fit's
    # E(B-V) with another's Teff double-counts the reddening -- it moved the
    # predicted Balmer residual from +1.2% to +2.5% on HD194453 for no physical
    # reason.
    if node is not None:
        tn = (node['teff'], node['logg'], node['mh'])
        node_label = 'ML node (scan)'
    elif sc is not None and ebv == float(sc['ebv']):
        tn = (float(sc['teff']), float(sc['logg_node']), float(sc['mh_node']))
        node_label = 'scan solution'
    else:
        tn = nearest_node(grid, *t0)
        node_label = 'nearest node'
    # ONE model: the fitted solution. A second curve at the catalog parameters
    # used to be drawn alongside it, but it had to borrow this fit's E(B-V) and
    # v sin i -- nothing refit them for the catalog node -- so its residuals
    # were not those of any fit that was actually performed, and printing them
    # next to the real ones invited reading a handicapped comparison as a
    # result. The catalog values are in data/sample.csv; the node-scan figures
    # (scan_<star>.png) mark them against the likelihood surface, which is the
    # honest place for that comparison.
    thetas = [(node_label, dict(teff=tn[0], logg=tn[1], mh=tn[2],
                                ebv=ebv, vsini=vsini))]
    print(f'\n  catalog Teff={t0[0]:.0f} logg={t0[1]:.2f} [M/H]={t0[2]:+.2f}  '
          f'-> used Teff={tn[0]:.0f} logg={tn[1]:.2f} [M/H]={tn[2]:+.2f}')
    return _fit_and_draw(star, grid, row, bands, nb, xs, spec, held, thetas,
                         ebv, vsini, node)


def core_excess(h, r, hm):
    """Median residual in the Balmer cores minus the median outside them."""
    w = np.asarray(h.wavelength, float)
    core = np.zeros(w.shape, bool)
    for l in hydrogen_lines(float(w.min()), float(w.max()), series=(2,),
                            nmax=20):
        core |= np.abs(w - l) < CORE_HALF_A
    cc, ww = hm & core, hm & ~core
    if not cc.any() or not ww.any():
        return None
    return float(100 * (np.nanmedian(r[cc]) - np.nanmedian(r[ww])))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--star', default='HD194453')
    ap.add_argument('--all', action='store_true',
                    help='every primary + secondary star in the sample')
    # Defaults of None, not 0. E(B-V) = 0 produced a figure whose 13
    # conditioning bands ran -5.5% to +5.0% across the spectrum -- the
    # unreddened tilt, which a single scalar cannot absorb and is not meant to
    # -- and the held-out break panels inherited it, which reads as a broken
    # normalisation rather than as the null dust hypothesis it is. So the
    # default is now the sweep's own solution when it exists, and the figure
    # says which it used.
    ap.add_argument('--vsini', type=float, default=None)
    ap.add_argument('--ebv', type=float, default=None)
    ap.add_argument('--no-scan', action='store_true',
                    help='ignore data/ebv_teff_scan.csv; E(B-V)=0 unless given')
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


def _fit_and_draw(star, grid, row, bands, nb, xs, spec, held, thetas, ebv,
                  vsini, node=None):
    results, xsl_cal, hiord = {}, {}, {}
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
            if k == 'Balmer':
                # CORE minus CONTINUUM, reported separately because a window
                # median mixes two different things: the continuum shape across
                # the break, which is the question, and the line-core residual,
                # which is sensitive to the instrument profile and to hydrogen
                # physics instead. Quoting one number for both is what let the
                # core excess be argued about for months.
                #
                # It is also DEGENERATE WITH THE LSF WIDTH -- at HD194453 it
                # runs +2.65% at a 3.54 A Moffat core to -4.20% at 7.00 A -- so
                # it is only interpretable alongside the profile in use. See
                # docs/LSF.md.
                ce = core_excess(h, r, hm)
                if ce is not None:
                    print(f'    {"":<7} core - continuum  {ce:+6.2f}%')

        rx = None
        if xs is not None:
            px = predict(th, [xs], grid)[0]
            # fill_domains so the model is drawn through the masked Balmer core
            # too; it cannot affect chi2, which selects on obs.mask
            cal_x, _ = solve(xs, px.value, fill_domains=True)
            rx, ux = residual(xs, cal_x), usable(xs, px.value)
            xsl_cal[label] = cal_x
            print(f'    xsl     chi2/N = '
                  f'{chi2(xs, cal_x) / max(int(ux.sum()), 1):7.2f}'
                  f'   rms = {100 * np.nanstd(rx[ux]):5.2f}%   n={int(ux.sum())}')

            # H11 and H10 at XSL resolution. Nothing fits them, so these are
            # predictions on the same footing as the break windows -- and they
            # are the LSF-independent version of the core excess the Balmer
            # panel reports. See H10_H11_N.
            merged = merge_members([member_panel(xs, px.value, n)
                                    for n in H10_H11_N])
            if merged is not None:
                hiord[label] = merged
                for d in merged['members']:
                    cs, ws = core_stat(d), wing_stat(d)
                    if cs is None:
                        continue
                    print(f'    H{d["n"]:<6d} PREDICTED  core   {cs[0]:+6.2f}%'
                          f'   rms {cs[1]:5.2f}%   n={cs[2]}'
                          f'   (over ±{d["core_hw"]:.1f} A, masked ±'
                          f'{d["core"]:.1f} A; wings {ws:+.2f}%)')
        results[label] = dict(cal_b=cal_b, rb=rb, ub=ub, scalar=c[0],
                              model_spec=model_spec, hres=hres)

    figure(star, nb, spec, xs, held, results, xsl_cal, row, bands,
           fixed=dict(ebv=ebv, vsini=vsini), node=node, hiord=hiord)
    labels = list(results)
    return dict(star=star, ebv=ebv, vsini=vsini,
                balmer=float(np.nanmedian(results[labels[0]]['hres']['Balmer'])),
                paschen=float(np.nanmedian(results[labels[0]]['hres']['Paschen'])))


def figure(star, nb, spec, xs, held, results, xsl_cal, row, bands,
           fixed=None, node=None, hiord=None):
    labels = list(results)
    # The high-order row is conditional: a star with no XSL, or whose XSL does
    # not reach 3760 A, gets the figure it always got rather than two blank
    # panels claiming a check that was not made.
    hiord = hiord or {}
    heights = [1.9, 1.5, 0.95, 1.4] + ([1.35, 0.85] if hiord else []) + [1.3]
    fig = plt.figure(figsize=(12.5, 15 + (3.0 if hiord else 0)))
    # hspace was .48, which already let every row's title sit on the xlabel of
    # the row above it; a seventh row and a two-line title made that collision
    # unreadable rather than merely untidy.
    gs = fig.add_gridspec(len(heights), 2, height_ratios=heights,
                          hspace=.62, wspace=.22)
    fig.patch.set_facecolor(SURFACE)
    m = spec.mask

    def shade(ax):
        for _, lo, hi in bands:
            ax.axvspan(lo, hi, color=BAND_C, alpha=.10, lw=0)
        for win in (BREAK_WINDOW, PASCHEN_WINDOW):
            ax.axvspan(*win, color=HELD_C, alpha=.13, lw=0)

    ax = fig.add_subplot(gs[0, :])
    style(ax)
    shade(ax)
    for lam in (BALMER, PASCHEN):
        ax.axvline(lam, color=MUTED, ls='--', lw=1)
    ax.plot(spec.wavelength[m], spec.flux[m], color=OBS_C, lw=1.2,
            label='NGSL observed')
    for lab, c in zip(labels, (MOD_C, MOD2_C)):
        ax.plot(spec.wavelength[m], results[lab]['model_spec'][m], color=c,
                lw=1.0, label=f'model, {lab} (scalar from bands)')
    ax.set_yscale('log')
    ax.set_ylabel(r'F$_\lambda$', fontsize=9, color=INK)
    ax.set_title(f'{star} — blue = bands used for conditioning;  '
                 'red = held out and PREDICTED', fontsize=10, color=INK)
    ax.set_xlabel(r'Wavelength [$\AA$, vacuum]', fontsize=9, color=INK)
    ax.legend(fontsize=8, loc='upper right', framealpha=.92)

    # Panel limits are per-break, not a uniform offset: the Balmer panel starts
    # at 3450 A so the sub-limit continuum band that anchors the prediction is
    # visible without the rest of the blue squeezing the break itself.
    # Fixed residual ranges, so the break panels are comparable between stars.
    # A per-panel autoscale made a 0.5% residual and a 15% one look identical,
    # which is the opposite of what these panels are for. Balmer gets the wider
    # range because the high-order line cores sit at +5 to +10% there (the known
    # NLTE excess), while Paschen's continuum rarely leaves +/-3%.
    panels = (('Balmer', BREAK_WINDOW, BALMER, (3450., 4150.), (-15., 15.)),
              ('Paschen', PASCHEN_WINDOW, PASCHEN, (7830., 9650.), (-8., 8.)))
    for col, (nm, win, lam0, xlim, rlim) in enumerate(panels):
        axz = fig.add_subplot(gs[1, col])
        axr = fig.add_subplot(gs[2, col], sharex=axz)
        lo, hi = xlim
        sel = m & (spec.wavelength > lo) & (spec.wavelength < hi)
        for a_ in (axz, axr):
            style(a_)
            shade(a_)
            a_.axvline(lam0, color=MUTED, ls='--', lw=1)
            a_.set_xlim(lo, hi)
        axz.plot(spec.wavelength[sel], spec.flux[sel], color=OBS_C, lw=1.3,
                 label='NGSL observed')
        for lab, c in zip(labels, (MOD_C, MOD2_C)):
            axz.plot(spec.wavelength[sel], results[lab]['model_spec'][sel],
                     color=c, lw=1.0, label=f'model, {lab}')
        med = 100 * np.nanmedian(results[labels[0]]['hres'][nm])
        axz.set_title(f'{nm} — HELD OUT and predicted '
                      f'({labels[0]}: {med:+.2f}%)', fontsize=9, color=INK)
        axz.set_ylabel(r'F$_\lambda$', fontsize=8, color=INK)
        axz.tick_params(labelbottom=False)
        if col == 0:
            axz.legend(fontsize=7, loc='upper left', framealpha=.92)

        # residual panel: ONLY inside the held-out window is a prediction --
        # everything outside it was used to set the scalar, so the two must be
        # drawn differently or the figure invites reading a fit as a prediction.
        inwin = (spec.wavelength >= win[0]) & (spec.wavelength <= win[1])
        for lab, c in zip(labels, (MOD_C, MOD2_C)):
            full = ((spec.flux - results[lab]['model_spec'])
                    / results[lab]['model_spec'])
            out = sel & ~inwin
            axr.plot(spec.wavelength[out], full[out] * 100, color=c, lw=.7,
                     alpha=.35)
            hin = sel & inwin
            axr.plot(spec.wavelength[hin], full[hin] * 100, color=c, lw=1.2)
            axr.axhline(100 * np.nanmedian(full[hin]), color=c, ls=':', lw=1)
        axr.axhline(0, color=MUTED, lw=1)
        r_all = 100 * np.concatenate([
            (((spec.flux - results[l]['model_spec']) / results[l]['model_spec'])
             [sel & inwin]) for l in labels])
        axr.set_ylim(*rlim)
        # A fixed range can hide a point off the top. Say so rather than let the
        # panel imply the residual stayed inside it.
        n_out = int(np.sum(np.isfinite(r_all) & ((r_all < rlim[0]) | (r_all > rlim[1]))))
        if n_out:
            axr.text(0.99, 0.04, f'{n_out} px outside ±{rlim[1]:.0f}% '
                                 f'(to {np.nanmax(np.abs(r_all)):.0f}%)',
                     transform=axr.transAxes, ha='right', va='bottom',
                     fontsize=6.5, color=HELD_C,
                     bbox=dict(fc='white', ec='none', alpha=.8, pad=1.5))
        axr.set_ylabel('(obs−model)/model [%]', fontsize=8, color=INK)
        axr.set_xlabel(r'$\lambda$ [$\AA$]', fontsize=8, color=INK)

    if xs is not None:
        for col, (lo, hi, ttl) in enumerate((
                (4292., 4392., r'XSL: H$\gamma$ wings — core masked (NLTE)'),
                (4125., 4139., 'XSL: Si II 4129/4132 + Fe II — fitted window'))):
            spectrum_panel(fig.add_subplot(gs[3, col]), xs,
                           [(f'model, {lab}', xsl_cal[lab]) for lab in labels],
                           lo, hi, title=ttl, legend=(col == 0))

    if hiord:
        first = hiord[labels[0]]
        lo, hi = first['window']
        axh = fig.add_subplot(gs[4, :])
        spectrum_panel(
            axh, first['view'],
            [(f'model, {lab}', hiord[lab]['cal']) for lab in labels], lo, hi,
            obs_label='XSL (solid = sets the local scale)')
        axh.set_xlabel('')
        axh.tick_params(labelbottom=False)
        # Legend bottom RIGHT and the line names at the TOP: the bottom left of
        # this panel is H11's core and the tops of the line centres are the one
        # part of an absorption panel guaranteed to be empty.
        axh.legend(fontsize=7, loc='lower right', framealpha=.92)
        bits = []
        for d in first['members']:
            axh.axvline(d['lam'], color=MUTED, ls='--', lw=1)
            axh.text(d['lam'], 0.97, f'H{d["n"]}', fontsize=8, color=INK,
                     transform=axh.get_xaxis_transform(), ha='center',
                     va='top',
                     bbox=dict(fc=SURFACE, ec='none', alpha=.85, pad=1.5))
            cs = core_stat(d)
            if cs is not None:
                bits.append(f'H{d["n"]} {cs[0]:+.2f}%')
        axh.set_title(
            'XSL, ~16× the NGSL resolution — H11 and H10 are inside the '
            'held-out window and fitted by NOTHING\n'
            'core residual ' + ', '.join(bits)
            + '   (shaded = masked core, predicted; the local scale is solved '
              'on the wings only)', fontsize=9, color=INK, linespacing=1.4)

        axhr = fig.add_subplot(gs[5, :], sharex=axh)
        residual_panel(
            axhr, first['view'].wavelength,
            [(hiord[lab]['resid'], hiord[lab]['used'], hiord[lab]['pred'],
              [m['incore'] for m in hiord[lab]['members']])
             for lab in labels], lo, hi, ylim=H10_H11_RLIM,
            bands=[(d['lam'] - d['core'], d['lam'] + d['core'])
                   for d in first['members']])
        for d in first['members']:
            axhr.axvline(d['lam'], color=MUTED, ls='--', lw=1)

    axb = fig.add_subplot(gs[len(heights) - 1, :])
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

    fixed = fixed or {}
    # E(B-V) and v sin i are HELD, not fitted, in this check -- and the held-out
    # residuals move a lot with both, so the title must say what they were or
    # the numbers in it cannot be compared between runs.
    if node is not None:
        line2 = (f'ML node: Teff={node["teff"]:.0f} / log g={node["logg"]:.2f} / '
                 f'[M/H]={node["mh"]:+.2f} / E(B−V)={node["ebv"]:.3f} / '
                 f'v sin i={node["vsini"]:.0f} km/s   — all fitted together '
                 f'over 1705 nodes'
                 + ('   ⚠ [M/H] AT THE GRID FLOOR' if node['at_mh_floor'] else ''))
        line3 = ('E(B−V) and v sin i are fitted at every node; the held-out '
                 'Balmer and Paschen windows are never fitted')
    else:
        line2 = (f'nominal Teff={row["teff_ngsl"]} / log g={row["logg_ngsl"]} / '
                 f'[M/H]={row["mh_ngsl"]}')
        line3 = (f'held fixed: E(B−V)={fixed.get("ebv", 0.0):.3f}, '
                 f'v sin i={fixed.get("vsini", 0.0):.0f} km/s')
    fig.suptitle(
        f'{star}   conditioning = NGSL bands + XSL lines;   '
        f'Balmer and Paschen held out\n{line2}\n{line3}',
        fontsize=10.5, color=INK, linespacing=1.5)
    out = figure_path('predict_check', star)
    fig.savefig(out, dpi=170, facecolor=SURFACE, bbox_inches='tight')
    plt.close(fig)
    print(f'\n  -> {out.relative_to(ROOT)}')


if __name__ == '__main__':
    main()
