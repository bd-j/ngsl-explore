"""NGSL + XSL + model, for the stars that appear in both libraries.

Same layout as plot_ngsl_vs_model.py, with the XSL spectrum overlaid: an
independent, ground-based, ~10x higher-resolution observation of the same star.

Each observation is put on the same footing before comparison:
  NGSL  air -> vacuum plus the fitted per-grating residual (ngsl_wavecal)
  XSL   air -> vacuum; already rest-frame, so no RV to remove
  model degraded separately to each instrument's LSF -- these differ in KIND,
        not just width: NGSL is constant in Angstroms per grating (fixed
        dispersion), XSL is constant in velocity (sigma(v) = 13 km/s in UVB).

Writes figures/ngsl_xsl_<star>.png
"""
import csv
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.ndimage import gaussian_filter1d

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common.lsf import (broaden_ngsl, degrade_to, rebin_to_pixels,
                        NGSL_LSF_TABULATED, NGSL_R_MEASURED)
from common.ngsl_wavecal import apply_wavecal, load_table
from common.xsl_load import load as load_xsl, resolving_power
from grid.make_model import hnu_to_flam
from explore.plot_ngsl_vs_model import (BALMER, PASCHEN, ZOOM_B, ZOOM_P,
                                        WINS_B, WINS_P, WINS_FULL, norm_mask,
                                        OBS_C, MOD_C, SURFACE, INK, MUTED, GRID)

ROOT = Path(__file__).resolve().parent.parent
XSL_C = '#1baf7a'          # categorical slot 3; slots 1/2 are NGSL and model
XSLD_C = '#4a3aa7'         # XSL degraded to NGSL: dashed, to read as derived


def xsl_to_ngsl(wx, fx, wn):
    """XSL convolved to the NGSL LSF and rebinned onto NGSL pixels.

    This is the apples-to-apples line: if a feature in the NGSL residual is a
    resolution artifact, it must appear here too, since this is what NGSL would
    have recorded had it observed the XSL spectrum. If it does NOT appear, the
    feature is in the NGSL data, not in the smoothing.

    XSL's own resolution is removed in quadrature, and the rebinning integrates
    across each NGSL pixel rather than sampling at its centre.
    """
    # Use the MEASURED NGSL profile (constant R = 600), not the STIS-table
    # one: matching XSL to NGSL for stars in common shows the delivered spectra
    # are 1.7-1.8x broader than the tables and constant in velocity, not in
    # Angstroms. See explore/ngsl_lsf_from_xsl.py.
    out = np.full_like(wn, np.nan, dtype=float)
    for lo, hi, _ in NGSL_LSF_TABULATED:
        seg_n = (wn >= lo) & (wn < hi)
        if not seg_n.any():
            continue
        mid = 0.5 * (max(lo, wx[0]) + min(hi, wx[-1]))
        fwhm_xsl = mid / resolving_power(mid)        # XSL FWHM in A here
        fwhm_ngsl = mid / NGSL_R_MEASURED
        deg = degrade_to(wx, fx, fwhm_xsl, fwhm_ngsl)
        out[seg_n] = rebin_to_pixels(wx, deg, wn[seg_n])
    return out
FULL = (3300.0, 9400.0)
WAVECAL = load_table()


def xsl_broaden(w, f):
    """Degrade a model to the XSL LSF: constant in velocity, so a fixed pixel
    width on a log grid. sigma(v) = 13 km/s (UVB) / 11 (VIS)."""
    out = np.array(f, float)
    step = np.median(np.diff(np.log(w)))          # model grid is log-lambda
    for lo, hi, sig in [(3000., 5600., 13.0), (5600., 10200., 11.0)]:
        seg = (w >= lo) & (w < hi)
        if seg.any():
            out[seg] = gaussian_filter1d(f, (sig / 2.99792458e5) / step,
                                         mode='nearest')[seg]
    return out


def make_figure(star, xslid, cat):
    d = fits.getdata('data/spectra/' + cat[star]['file'])
    wn = apply_wavecal(d['WAVELENGTH'].astype(float), star, WAVECAL)
    fn = d['FLUX'].astype(float)
    wx, fx, ex, hdr = load_xsl(xslid)

    wm, hnu, _ = np.loadtxt(f'models/work/{star}.spec', unpack=True)
    fm = hnu_to_flam(wm, hnu)
    m_ngsl = broaden_ngsl(wm, fm)
    m_xsl = xsl_broaden(wm, fm)

    lo, hi = max(FULL[0], wm[0], wx[0]), min(FULL[1], wm[-1], wx[-1])
    kn = (wn > lo) & (wn < hi) & np.isfinite(fn) & (fn > 0)
    wn, fn = wn[kn], fn[kn]
    kx = (wx > lo) & (wx < hi)
    wx, fx = wx[kx], fx[kx]
    xd = xsl_to_ngsl(wx, fx, wn)          # XSL as NGSL would have seen it

    cols = [((lo, hi), 'full range', None, WINS_FULL, 'upper right'),
            (ZOOM_B, f'Balmer break ({BALMER:.0f} $\\AA$)', BALMER, WINS_B, 'lower right'),
            (ZOOM_P, f'Paschen break ({PASCHEN:.0f} $\\AA$)', PASCHEN, WINS_P, 'lower left')]
    fig, axes = plt.subplots(6, 1, figsize=(11.5, 17.5),
                             gridspec_kw={'height_ratios': [2.1, 1] * 3})
    fig.patch.set_facecolor(SURFACE)

    for i, (xlim, ttl, mark, wins, legloc) in enumerate(cols):
        ax, rax = axes[2 * i], axes[2 * i + 1]
        kmn = norm_mask(wn, wins)
        kmx = norm_mask(wx, wins)
        if kmn.sum() < 20:
            kmn = np.ones_like(wn, bool)
        if kmx.sum() < 20:
            kmx = np.ones_like(wx, bool)
        mn = np.interp(wn, wm, m_ngsl)
        mx = np.interp(wx, wm, m_xsl)
        mn = mn * np.median(fn[kmn] / mn[kmn])
        good = np.isfinite(xd)
        kd = kmn & good

        # The two libraries differ by a GREY scale factor -- flat in wavelength,
        # star-dependent, -11% to +4% here. It is an absolute flux-calibration
        # difference (ground-based XSL must correct slit losses; NGSL need not),
        # and smoothing cannot cause it since convolution conserves flux.
        # Scaling each series to its OWN model made the flux panels disagree by
        # that factor and look like a shape difference. Everything is now tied
        # to the NGSL scale so the panels compare shape, and the grey factor is
        # reported in the legend instead of hidden.
        grey = float(np.median(xd[kd] / fn[kd])) if kd.sum() > 10 else 1.0
        xds = xd / grey
        mx_s = mx * np.median(fx[kmx] / mx[kmx]) / grey
        fx_s = fx / grey
        rn = (fn - mn) / mn
        rx = (fx - mx * np.median(fx[kmx] / mx[kmx])) / (mx * np.median(fx[kmx] / mx[kmx]))
        rd = np.where(good, (xds - mn) / mn, np.nan)

        for a in (ax, rax):
            a.set_facecolor(SURFACE)
            a.fill_between(wn, 0, 1, where=kmn, transform=a.get_xaxis_transform(),
                           color=MUTED, alpha=.13, lw=0, zorder=0)
            for m_ in ([BALMER, PASCHEN] if mark is None else [mark]):
                a.axvline(m_, color=MUTED, ls='--', lw=1, zorder=1)
            a.set_xlim(*xlim)
            a.grid(alpha=.25, color=GRID, lw=.7)
            a.tick_params(labelsize=8, colors=MUTED)
            for sp in a.spines.values():
                sp.set_color(GRID)
        ax.plot(wx, fx_s, color=XSL_C, lw=0.8, zorder=5,
                label=f'XSL DR3  (R~{resolving_power(np.mean(xlim)):.0f}), '
                      f'/{grey:.3f} to the NGSL flux scale')
        ax.plot(wn, fn, color=OBS_C, lw=1.3, zorder=4, label='NGSL (wavecal applied)')
        ax.plot(wn, xds, color=XSLD_C, lw=1.4, ls='--', zorder=6,
                label='XSL degraded to NGSL LSF + pixels')
        ax.plot(wn, mn, color=MOD_C, lw=1.3, zorder=3, label='ATLAS12, NGSL LSF')
        inw = (wn > xlim[0]) & (wn < xlim[1])
        iwx = (wx > xlim[0]) & (wx < xlim[1])
        allv = np.concatenate([fn[inw], mn[inw], fx_s[iwx]])
        ylo, yhi = np.nanpercentile(allv, [0.2, 99.8])
        ax.set_ylim(ylo - .08 * (yhi - ylo), yhi + .08 * (yhi - ylo))
        ax.set_title(f'{star} — {ttl}', fontsize=10, color=INK)
        ax.set_ylabel(r'F$_\lambda$  (locally scaled)', fontsize=9, color=INK)
        ax.legend(fontsize=7.5, loc=legloc, framealpha=.92)

        rax.axhline(0, color=MUTED, lw=1, zorder=2)
        rax.plot(wx, rx * 100, color=XSL_C, lw=0.7, zorder=3, label='XSL - model')
        rax.plot(wn, rd * 100, color=XSLD_C, lw=1.2, ls='--', zorder=6,
                 label='XSL@NGSL - model')
        rax.fill_between(wn, rn * 100, 0, color=OBS_C, alpha=.28, lw=0, zorder=4)
        rax.plot(wn, rn * 100, color=OBS_C, lw=.9, zorder=5, label='NGSL - model')
        rax.set_xlabel(r'Wavelength [$\AA$, vacuum]', fontsize=9, color=INK)
        rax.set_ylabel('(obs - model) / model  [%]', fontsize=9, color=INK)
        rax.set_ylim(*np.nanpercentile(np.concatenate([rn[inw], rx[iwx],
                                       rd[inw][np.isfinite(rd[inw])]]) * 100,
                                       [1, 99]) + np.array([-5, 5]))
        if i == 1:
            rax.legend(fontsize=7, loc='upper right', framealpha=.9)
            stats = (np.nanstd(rn[inw]) * 100, np.nanstd(rx[iwx]) * 100,
                     np.nanstd(rd[inw]) * 100)

    r = cat[star]
    fig.suptitle(
        f'{star}   NGSL Teff={float(r["teff_ngsl_K"]):.0f} K  log g={r["logg_ngsl"]}  '
        f'[M/H]={float(r["m_h_ngsl"]):+.1f}\n'
        f'Balmer-region residual RMS:  NGSL {stats[0]:.1f}%   '
        f'XSL {stats[1]:.1f}%   XSL@NGSL {stats[2]:.1f}%',
        fontsize=11, color=INK, linespacing=1.5)
    fig.tight_layout(rect=[0, 0, 1, 0.962])
    out = ROOT / 'figures' / f'ngsl_xsl_{star}.png'
    fig.savefig(out, dpi=200, facecolor=SURFACE)
    plt.close(fig)
    print(f'  {star}: RMS  NGSL {stats[0]:.1f}%  XSL {stats[1]:.1f}%  '
          f'XSL@NGSL {stats[2]:.1f}%  -> {out.name}')


def main():
    cat = {r['star']: r for r in
           csv.DictReader(open(ROOT / 'data' / 'balmer_candidates.csv'))}
    xsl = {}
    for r in csv.DictReader(open(ROOT / 'data' / 'xsl_all.csv')):
        key = r['star'].replace(' ', '')
        if key in cat and (ROOT / 'data' / 'xsl' / 'XSL_DR3_release' /
                           f'xsl_spectrum_{r["xslid"]}_merged.fits').exists():
            xsl.setdefault(key, r['xslid'])
    have = {s: i for s, i in xsl.items()
            if (ROOT / 'models' / 'work' / f'{s}.spec').exists()}
    print(f'{len(have)} stars in both NGSL and XSL with a model:')
    for s, i in sorted(have.items()):
        make_figure(s, i, cat)


if __name__ == '__main__':
    main()
