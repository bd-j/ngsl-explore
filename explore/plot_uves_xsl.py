"""UVES-POP smoothed to XSL resolution, for the 9 stars in both libraries.

The 9 come from data/uves_xsl_overlap.csv (explore/uves_xsl_overlap.py). Both
libraries are absolutely flux calibrated, so this asks how far apart their
flux calibrations are on the same star -- a data-to-data comparison with no
model anywhere in it.

THE SMOOTHING GOES UVES-POP -> XSL, WHICH IS THE ONLY POSSIBLE DIRECTION.
At the Balmer break XSL is R ~ 9800 (sigma_v = 13 km/s, and XSL quotes sigma,
not FWHM) while UVES-POP delivered is R ~ 18,000 on its 0.1 A grid and
R = 80,000 native. UVES-POP is the sharper of the two by about 1.9x, so
smoothing XSL to UVES-POP would be a deconvolution. Same rule as
`to_ngsl_pixels` in plot_uves_ngsl.py: the high-resolution spectrum goes to
the low-resolution one, never the reverse.

THE KERNEL IS CONSTANT IN VELOCITY, because XSL's LSF is -- unlike NGSL's,
which is constant in Angstroms. So the convolution runs on a log-lambda grid
(common.lsf.broaden_R) and the spectrum is then integrated onto XSL's own
pixels with `rebin_to_pixels`, not interpolated.

The kernel is the QUADRATURE DIFFERENCE, sigma_v = sqrt(13^2 - sigma_uves^2),
where sigma_uves combines UVES-POP's native R = 80,000 (1.59 km/s) with its
0.1 A delivered pixel (a boxcar, 8.2 km/s wide at the break, so 2.37 km/s rms).
That gives 12.68 km/s against the 13.0 target: UVES-POP's own resolution is a
2.5% correction to the kernel width and the comparison does not rest on it.
Taking the 2-pixel figure (R ~ 18,000, 16.4 km/s) instead would give 25.8 km/s
-- but that number is a SAMPLING limit, not an LSF, and using it as one would
over-broaden by a factor 2.

SIGN CONVENTION: the residual is (XSL - UVES@XSL) / UVES@XSL, in percent.
POSITIVE MEANS XSL IS BRIGHTER than UVES-POP smoothed to XSL's resolution.
Neither is a model, so neither is "truth" -- this is a disagreement between two
flux calibrations, not an error in one of them.

THREE THINGS HAVE TO BE RIGHT BEFORE THE RESIDUAL MEANS ANYTHING

1. Frames. XSL is REST-frame (the RV is already removed) and UVES-POP's
   delivered spectra are in the OBSERVED frame -- established in
   plot_uves_ngsl.py against the Ca II minima. So the catalog RV is removed
   from UVES-POP here and nothing is done to XSL. The RVs run -171 to +101
   km/s across these 9, which is up to 2.1 A at the break: not optional.

2. Air vs vacuum. Both libraries are in AIR as delivered and both loaders
   convert to vacuum by default, so this is handled but worth naming, since
   the two conversions have to be the SAME one (common.lines.air_to_vac).

3. Slit-loss correction -- and THIS ONE BITES, because it is exactly the
   quantity being compared. XSL ships variants: only the plain `_merged.fits`
   has LOSS_COR = True. THREE OF THESE NINE DO NOT --

       HD102212   _ncl_ncge   LOSS_COR = False
       HD099648   _scl        LOSS_COR = False
       Betelgeuse _ncl_ncge   LOSS_COR = False

   For those three the XSL flux carries an uncorrected slit loss, so their
   grey factor is NOT a statement about either library's calibration. They are
   drawn, because the SHAPE is still comparable after normalisation, but their
   grey factor is flagged in the table and excluded from the summary. The flag
   is read from the header per star, not hard-coded from this list.

WHAT THIS FIGURE CANNOT SETTLE: the two libraries observed these stars years
apart and several are catalogued variables -- Betelgeuse (SRC) and nu Vir (SRB)
above all. For those there is a floor on how well any comparison can agree and
it is not an instrumental floor. The GCVS type comes out of the UVES-POP file
and is printed on every panel; a variable's panel is not a calibration
measurement.

THE PANEL STARTS AT 3502 A, not the 3400 A of the NGSL version, because XSL's
blue limit is 3501.0 A on all nine. That leaves 144 A blueward of the break to
normalise on, against 246 A there.

Writes figures/explore_libraries/uves_xsl_break.png     (Balmer break, 3646 A)
       figures/explore_libraries/uves_xsl_hepsilon.png  (H-epsilon, 3971 A)
       data/uves_xsl_compare.csv                        (the numbers)
"""
import csv
import sys
import argparse
import warnings
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.ndimage import median_filter

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common.lsf import broaden_R, rebin_to_pixels
from common.lines import balmer_member, ISM_LINES
from common.xsl_load import load as load_xsl, loss_corrected
from common.uves_pop_load import load as load_uves, NATIVE_R, GRID_STEP
from common.figpath import library_figure_path
from explore.plot_uves_ngsl import (fetch_uves, uves_meta, fill_small_gaps,
                                    SETTING_JOINS)
from explore.plot_ngsl_vs_model import (BALMER, OBS_C, SURFACE, INK, MUTED,
                                        GRID)

ROOT = Path(__file__).resolve().parent.parent
UVES_C = '#1baf7a'          # raw UVES-POP, as delivered
UVESD_C = '#4a3aa7'         # UVES degraded to XSL: dashed, to read as derived
JOIN_C = '#b07a2a'          # UVES setting joins: instrumental, not stellar
C_KMS = 2.99792458e5

XSL_SIGMA_V = 13.0          # km/s, UVB arm; XSL quotes sigma, NOT FWHM
HEPS = balmer_member(7)                  # 3971.24 A vacuum
CAH = ISM_LINES['Ca II H']               # 3969.60 A vacuum, 1.64 A away

# XSL's blue limit is 3501.0 A on every one of these 9, so the break panel
# cannot start at the NGSL version's 3400 A; 3502 keeps it inside the data on
# all of them. The normalisation window stops 26 A short of the break.
BREAK_XLIM = (3502.0, 4000.0)
BREAK_WINS = [(3510., 3620.)]
HEPS_XLIM = (HEPS - 30.0, HEPS + 30.0)
HEPS_WINS = [(3944., 3962.), (3982., 4000.)]
GAP_FILL_MAX = 2.0          # A, as in plot_uves_ngsl


def uves_sigma_v(lam):
    """UVES-POP's own delivered sigma(v) in km/s: native LSF + the 0.1 A pixel.

    The native R = 80,000 contributes c/R/2.3548 = 1.59 km/s. The delivered
    0.1 A grid contributes a boxcar of width c*0.1/lam, whose rms is that over
    sqrt(12). At 3646 A the two give 1.59 and 2.37 -> 2.86 km/s combined.
    """
    s_native = C_KMS / NATIVE_R / 2.3548
    s_pixel = (C_KMS * GRID_STEP / lam) / np.sqrt(12.0)
    return np.hypot(s_native, s_pixel)


def uves_to_xsl(w_u, f_u, w_x, lam0):
    """UVES-POP -> XSL resolution and XSL pixels. ONE direction only.

    lam0 is where the kernel is evaluated. It is wavelength dependent only
    through the pixel term, which over 3500-4000 A moves the kernel from 12.76
    to 12.82 km/s -- 0.5%, so one kernel per panel is enough and the panel
    centre is used.
    """
    s_k = np.sqrt(XSL_SIGMA_V ** 2 - uves_sigma_v(lam0) ** 2)
    # Resample onto a uniform log-lambda grid, where a velocity-constant
    # kernel is a constant number of pixels. R = 200,000 oversamples the
    # 30 km/s FWHM being applied by a factor ~20.
    n = int(np.log(w_u[-1] / w_u[0]) * 200000.0) + 1
    wl = w_u[0] * np.exp(np.arange(n) / 200000.0)
    fl = np.interp(wl, w_u, f_u)
    sm = broaden_R(wl, fl, C_KMS / (2.3548 * s_k))
    return rebin_to_pixels(wl, sm, w_x), s_k


def prepare(row):
    """-> dict of everything one panel needs, or None if the star cannot load."""
    star = row['uves_name']
    xid = row['xslid'].split(';')[0]          # first epoch; repeats noted below
    wx, fx, ex, hx = load_xsl(xid)            # rest-frame, vacuum, FLUX
    wu, fu, eu = load_uves(star)              # observed frame, vacuum

    # UVES-POP is in the OBSERVED frame; XSL is already at rest.
    rv = float(row['rv_uves_kms']) if row['rv_uves_kms'] else 0.0
    wu = wu / (1.0 + rv / C_KMS)

    # NaN dropouts poison an FFT convolution, so short runs are interpolated
    # and long ones are blanked back out after smoothing.
    fu, holes = fill_small_gaps(wu, fu, GAP_FILL_MAX)
    gcvs, vtype, date = uves_meta(star)
    return dict(star=star, row=row, xslid=xid, n_ep=int(row['n_xsl_epochs']),
                wx=wx, fx=fx, wu=wu, fu=fu, rv=rv, holes=holes,
                gcvs=gcvs, vtype=vtype, date=date,
                loss_cor=loss_corrected(xid))


# A panel is only drawn if its GREY FACTOR can be measured, which needs real
# flux from both libraries inside the normalisation windows. Checking the
# wavelength endpoints is not enough: UVES-POP's Betelgeuse spectrum runs
# 3200-10250 A like the others but has a 552 A HOLE at 3200.9-3753.3 A, so it
# spans the break panel while carrying no data across it. Drawn without this
# test it came out as an empty panel with grey = nan.
MIN_NORM_PTS = 50            # UVES samples at 10/A, so this is a low bar


def usable(d, xlim, wins):
    """-> (ok, reason). Real data from BOTH libraries in the windows."""
    for w, f, lab in ((d['wx'], d['fx'], 'XSL'), (d['wu'], d['fu'], 'UVES')):
        if w[0] > xlim[0] or w[-1] < xlim[1]:
            return False, f'{lab} does not span {xlim[0]:.0f}-{xlim[1]:.0f} A'
    # The holes were interpolated across by fill_small_gaps, so test the
    # ORIGINAL hole list, not the filled flux.
    n_u = 0
    for lo, hi in wins:
        m = (d['wu'] >= lo) & (d['wu'] <= hi)
        keep = np.ones(int(m.sum()), bool)
        for a, b in d['holes']:
            keep &= ~((d['wu'][m] > a) & (d['wu'][m] < b))
        n_u += int(keep.sum())
    if n_u < MIN_NORM_PTS:
        gap = ', '.join(f'{a:.0f}-{b:.0f}' for a, b in d['holes']
                        if b > xlim[0] and a < xlim[1])
        return False, (f'UVES has {n_u} real points in the normalisation '
                       f'windows (hole at {gap} A)')
    n_x = 0
    for lo, hi in wins:
        m = (d['wx'] >= lo) & (d['wx'] <= hi)
        n_x += int(np.sum(np.isfinite(d['fx'][m]) & (d['fx'][m] > 0)))
    if n_x < 10:
        return False, f'XSL has {n_x} positive points in the windows'
    return True, str()


def panel(ax, rax, d, xlim, wins, marks):
    """-> (grey, rms%). Draws XSL, raw UVES, UVES@XSL and the residual."""
    lam0 = 0.5 * (xlim[0] + xlim[1])
    inx = (d['wx'] >= xlim[0] - 60) & (d['wx'] <= xlim[1] + 60)
    wx, fx = d['wx'][inx], d['fx'][inx]
    ud, s_k = uves_to_xsl(d['wu'], d['fu'], wx, lam0)

    # Blank the long UVES holes back out, so no pixel built from invented flux
    # is drawn. Margin is the kernel's reach: 4 sigma at 12.7 km/s.
    for lo, hi in d['holes']:
        m = 4.0 * s_k / C_KMS * lam0
        ud[(wx > lo - m) & (wx < hi + m)] = np.nan

    # Normalise on the windows: the grey factor is measured there and divided
    # out, so the panels compare SHAPE and the factor is reported separately.
    sel = np.zeros(wx.shape, bool)
    for lo, hi in wins:
        sel |= (wx >= lo) & (wx <= hi)
    ok = sel & np.isfinite(ud) & (ud > 0) & np.isfinite(fx) & (fx > 0)
    grey = float(np.median(fx[ok] / ud[ok])) if ok.sum() > 10 else np.nan
    uds = ud * grey                                  # UVES on XSL's flux scale

    inw = (wx >= xlim[0]) & (wx <= xlim[1])
    resid = np.where(np.isfinite(uds) & (uds > 0),
                     (fx - uds) / uds * 100.0, np.nan)

    for lo, hi in wins:
        ax.axvspan(lo, hi, color=GRID, alpha=.55, lw=0, zorder=1)
    for x, ls in marks:
        ax.axvline(x, color=MUTED, lw=1, ls=ls, zorder=2)
        rax.axvline(x, color=MUTED, lw=1, ls=ls, zorder=2)
    # UVES setting joins are instrumental, and they land inside the break panel.
    for lo, hi, _ in SETTING_JOINS:
        if hi > xlim[0] and lo < xlim[1]:
            ax.axvspan(lo, hi, color=JOIN_C, alpha=.13, lw=0, zorder=1)
            rax.axvspan(lo, hi, color=JOIN_C, alpha=.13, lw=0, zorder=1)

    iwu = (d['wu'] >= xlim[0]) & (d['wu'] <= xlim[1])
    ax.plot(d['wu'][iwu], d['fu'][iwu] * grey, color=UVES_C, lw=.5, alpha=.75,
            zorder=3, label='UVES-POP, as delivered')
    ax.plot(wx[inw], fx[inw], color=OBS_C, lw=1.1, zorder=5, label='XSL DR3')
    ax.plot(wx[inw], uds[inw], color=UVESD_C, lw=1.2, ls='--', zorder=6,
            label=f'UVES-POP @ XSL ({s_k:.1f} km/s)')

    vals = np.concatenate([fx[inw], uds[inw][np.isfinite(uds[inw])]])
    ylo, yhi = np.nanpercentile(vals, [0.5, 99.5])
    raw = (d['fu'][iwu] * grey)
    raw = raw[np.isfinite(raw)]
    if raw.size > 20:
        ylo = min(ylo, float(np.min(median_filter(raw, size=5, mode='nearest'))))
    pad = .10 * (yhi - ylo)
    ax.set_ylim(ylo - pad, yhi + pad)

    rax.axhline(0, color=MUTED, lw=1, zorder=2)
    rax.fill_between(wx[inw], resid[inw], 0, color=OBS_C, alpha=.28, lw=0,
                     zorder=3)
    rax.plot(wx[inw], resid[inw], color=OBS_C, lw=1.0, zorder=4)
    fin = resid[inw][np.isfinite(resid[inw])]
    rms = float(np.sqrt(np.nanmean(fin ** 2))) if fin.size else np.nan
    if fin.size:
        m = max(5.0, float(np.nanpercentile(np.abs(fin), 99)) * 1.25)
        rax.set_ylim(-m, m)
    return grey, rms


def make_figure(data, xlim, wins, marks, title, sub, outname):
    """-> {star: (grey, rms%)} for the stars actually drawn."""
    ncol = 3
    nrow = int(np.ceil(len(data) / ncol))
    nlines = 1 + f'{title}\n{sub}'.count('\n')
    head_in = 0.45 + nlines * 11.5 * 1.7 / 72.0
    fig_h = 4.1 * nrow + head_in + 0.62
    fig = plt.figure(figsize=(16.5, fig_h))
    fig.patch.set_facecolor(SURFACE)
    gs = fig.add_gridspec(nrow * 2, ncol, height_ratios=[2.2, 1] * nrow,
                          hspace=.42, wspace=.20)
    stats = {}
    for i, d in enumerate(data):
        r, c = divmod(i, ncol)
        ax = fig.add_subplot(gs[2 * r, c])
        rax = fig.add_subplot(gs[2 * r + 1, c], sharex=ax)
        grey, rms = panel(ax, rax, d, xlim, wins, marks)
        stats[d['star']] = (grey, rms)
        var = f"   {d['gcvs']} ({d['vtype']})" if d['vtype'] else ''
        # A star whose XSL flux is not slit-loss corrected gets it in the title,
        # because its grey factor must not be read as a calibration difference.
        nolc = '' if d['loss_cor'] else '   NO SLIT-LOSS CORR'
        ep = f"   {d['n_ep']} XSL epochs" if d['n_ep'] > 1 else ''
        ax.set_title(f"{d['star']} / {d['row']['xsl_name']}{var}{nolc}\n"
                     f"grey x{grey:.3f}   UVES RV {d['rv']:+.0f} km/s   "
                     f"{d['xslid']}{ep}",
                     fontsize=8,
                     color=(INK if (d['loss_cor'] and not d['vtype'])
                            else '#8a3324'))
        ax.tick_params(labelbottom=False)
        if c == 0:
            ax.set_ylabel(r'F$_\lambda$ (XSL scale)', fontsize=8, color=INK)
            rax.set_ylabel('(XSL − UVES)/UVES [%]', fontsize=8, color=INK)
        if i == 0:
            ax.legend(fontsize=6.5, loc='best', framealpha=.92)
        if i >= len(data) - ncol:
            rax.set_xlabel(r'Wavelength [$\AA$, vacuum, rest]', fontsize=8,
                           color=INK)
    fig.suptitle(f'{title}\n{sub}', fontsize=11.5, color=INK, linespacing=1.7,
                 y=1 - 0.18 / fig_h, va='top')
    fig.subplots_adjust(left=.055, right=.985, top=1 - head_in / fig_h,
                        bottom=0.62 / fig_h)
    out = library_figure_path(outname)
    fig.savefig(out, dpi=150, facecolor=SURFACE)
    plt.close(fig)
    print(f'  {len(data)} panels -> {out}')
    return stats


def selftest():
    """The kernel and the sign, against injected values.

    1. Degrading a spectrum that is ALREADY at XSL resolution must be a no-op
       in width: feed a Gaussian line of known sigma and check the recovered
       width grows by the quadrature amount and no more.
    2. The residual sign: make XSL 5% brighter by construction and require
       (XSL - UVES)/UVES = +5.
    """
    lam0, s_in = 3700.0, 8.0                       # km/s
    w = lam0 * np.exp(np.arange(-40000, 40000) / 400000.0)
    prof = np.exp(-0.5 * ((np.log(w / lam0) * C_KMS) / s_in) ** 2)
    out, s_k = uves_to_xsl(w, prof, w, lam0)
    v = np.log(w / lam0) * C_KMS
    g = np.isfinite(out)
    s_out = np.sqrt(np.sum(out[g] * v[g] ** 2) / np.sum(out[g]))
    want = np.hypot(s_in, s_k)
    assert abs(s_out - want) < 0.25 * want * 0.05, (s_out, want)
    print(f'  kernel: {s_in:.1f} km/s in, kernel {s_k:.2f}, out {s_out:.2f} '
          f'vs expected {want:.2f} km/s  OK')

    uves = np.ones_like(w)
    xsl = 1.05 * np.ones_like(w)
    r = (xsl - uves) / uves * 100.0
    assert abs(r.mean() - 5.0) < 1e-9
    print('  sign: XSL 5% brighter -> residual +5.0%  OK '
          '(positive means XSL is brighter)')


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--skip-selftest', action='store_true')
    a = ap.parse_args()
    if not a.skip_selftest:
        print('Self-test:')
        selftest()

    rows = list(csv.DictReader(open(ROOT / 'data' / 'uves_xsl_overlap.csv')))
    print(f'\n{len(rows)} stars in both libraries')
    data = []
    for r in rows:
        try:
            fetch_uves(r['uves_name'], r['uves_spec_url'])
            d = prepare(r)
        except Exception as e:                    # noqa: BLE001 - report, skip
            print(f"  {r['uves_name']}: unavailable ({e}), skipped")
            continue
        lc = '' if d['loss_cor'] else '  NO SLIT-LOSS CORRECTION'
        vt = f"  {d['vtype']}" if d['vtype'] else ''
        print(f"  {d['star']:<12} {d['xslid']}  RV {d['rv']:+7.1f} km/s{lc}{vt}")
        data.append(d)

    brk, hep = [], []
    for nm, keep, xl, wn in (('break', brk, BREAK_XLIM, BREAK_WINS),
                             ('H-epsilon', hep, HEPS_XLIM, HEPS_WINS)):
        for d in data:
            ok, why = usable(d, xl, wn)
            if ok:
                keep.append(d)
            else:
                print(f'  {nm}: {d["star"]} dropped -- {why}')

    sub_common = ('Positive residual = XSL brighter.  '
                  'Grey factor divided out per star, on the shaded windows.  '
                  'Red title = catalogued variable or no slit-loss correction.')
    print('\nBalmer break:')
    s_brk = make_figure(
        brk, BREAK_XLIM, BREAK_WINS, [(BALMER, '--')],
        'UVES-POP vs XSL DR3 at the Balmer break — two flux calibrations, '
        'same star, no model',
        'UVES-POP smoothed to XSL (sigma_v = 13 km/s) and integrated onto XSL '
        'pixels; the reverse is impossible, UVES-POP is the sharper.\n'
        + sub_common, 'uves_xsl_break.png')
    print('H-epsilon:')
    s_hep = make_figure(
        hep, HEPS_XLIM, HEPS_WINS, [(HEPS, '--'), (CAH, ':')],
        'UVES-POP vs XSL DR3 at H-epsilon (3971.24 Å, blended with '
        'interstellar Ca II H at 3969.60 Å)',
        'UVES-POP smoothed to XSL (sigma_v = 13 km/s) and integrated onto XSL '
        'pixels.  At 30.6 km/s FWHM the two lines are 1.64 Å apart and only '
        'marginally separated.\n' + sub_common, 'uves_xsl_hepsilon.png')

    out = []
    for d in data:
        gb, rb = s_brk.get(d['star'], (np.nan, np.nan))
        gh, rh = s_hep.get(d['star'], (np.nan, np.nan))
        f = lambda v, n: ('' if not np.isfinite(v) else round(float(v), n))
        out.append(dict(
            star=d['star'], xsl_name=d['row']['xsl_name'], xslid=d['xslid'],
            n_xsl_epochs=d['n_ep'], loss_corrected=d['loss_cor'],
            gcvs_type=d['vtype'], uves_date=d['date'],
            rv_kms=round(d['rv'], 2),
            grey_break=f(gb, 4), rms_break_pct=f(rb, 2),
            grey_heps=f(gh, 4), rms_heps_pct=f(rh, 2),
            # Ratio of the two grey factors = how much the XSL/UVES flux ratio
            # changes between the two normalisation windows, 3565 -> 3972 A.
            # Those windows straddle the break, so this is the libraries'
            # disagreement ACROSS THE BREAK, not a broadband colour term.
            d_blue_red_pct=f(100.0 * (gh / gb - 1.0), 2)
            if np.isfinite(gb) and np.isfinite(gh) and gb != 0 else '',
            teff_xsl=d['row']['teff_xsl'], teff_uves=d['row']['teff_uves']))
    out.sort(key=lambda r: r['star'])
    p = ROOT / 'data' / 'uves_xsl_compare.csv'
    with open(p, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0]))
        w.writeheader()
        w.writerows(out)
    print(f'-> {p}')

    print(f'\n{"star":<12}{"lossc":>6}{"var":>6}{"grey_brk":>10}{"rms%":>7}'
          f'{"grey_Heps":>11}{"rms%":>7}{"blue-red%":>11}')
    for r in out:
        print(f'{r["star"]:<12}{str(r["loss_corrected"]):>6}'
              f'{(r["gcvs_type"] or "-"):>6}{str(r["grey_break"]):>10}'
              f'{str(r["rms_break_pct"]):>7}{str(r["grey_heps"]):>11}'
              f'{str(r["rms_heps_pct"]):>7}{str(r["d_blue_red_pct"]):>11}')

    # The summary excludes the stars whose XSL flux is not slit-loss corrected:
    # their grey factor is an uncorrected aperture loss, not a calibration
    # difference, and pooling them would corrupt the one number worth quoting.
    g = [r['grey_break'] for r in out
         if r['loss_corrected'] and r['grey_break'] != '']
    if g:
        print(f'\nGrey factor XSL/UVES-POP at the break, '
              f'{len(g)} slit-loss-corrected stars: '
              f'median {np.median(g):.3f}, range {min(g):.3f} to {max(g):.3f}')
    bad = [r['star'] for r in out if not r['loss_corrected']]
    if bad:
        print(f'EXCLUDED from that summary (XSL LOSS_COR = False): '
              f'{", ".join(bad)}')

    # The grey factor is NOT the whole story: it is measured blueward of the
    # break, and the same stars do not keep it redward. This is the number the
    # residual panels show as a rise from 0 at 3500 A.
    t = [r['d_blue_red_pct'] for r in out
         if r['loss_corrected'] and r['d_blue_red_pct'] != '']
    if t:
        print(f'\nXSL/UVES flux ratio, 3565 A window -> 3972 A window '
              f'(i.e. across the break), {len(t)} slit-loss-corrected stars:\n'
              f'  median {np.median(t):+.1f}%, range {min(t):+.1f} to '
              f'{max(t):+.1f}%  -- positive = XSL relatively brighter to the RED'
              f'\n  So the two calibrations are NOT related by a grey factor '
              f'over this range.')


if __name__ == '__main__':
    sys.exit(main())
