"""Compare ATLAS12 models against UVES-POP spectra.

The processing chain differs from the NGSL comparison in three ways, all forced
by the data:

  dereddening   UVES-POP publishes a FITTED E(B-V) per star (not a map column),
                so the observation is dereddened with CCM89, R_V = 3.1.

  rotation      At this resolution rotation dominates the line profile. Models
                are broadened with each star's published v sin i using a Gray
                (2005) rotation profile with linear limb darkening. This is the
                library's own spectroscopic fit; refitting stays open.

  instrument    UVES-POP is R = 80,000 natively but delivered rebinned onto a
                0.1 A linear grid, so the model gets the native Gaussian AND a
                0.1 A boxcar. The boxcar dominates: 0.1 A vs a 0.046 A
                resolution element at 3646 A.

Writes figures/uves_vs_model_<star>.png
"""
import csv
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d, uniform_filter1d

sys.path.insert(0, str(Path(__file__).resolve().parent))
from make_atlas_model import hnu_to_flam, SYNTHE_R
from uves_pop_load import load as load_uves, params as uves_params
from extinction_ccm import deredden
from plot_model_vs_obs import (hydrogen_lines, BALMER, PASCHEN, H_MASK_A,
                               OBS_C, MOD_C, SURFACE, INK, MUTED, GRID)

ROOT = Path(__file__).resolve().parent.parent
UVES_R = 80000.0            # native resolving power
GRID_STEP = 0.1             # A, delivered sampling
C_KMS = 2.99792458e5

FULL = (3300.0, 9200.0)
ZOOM_B = (3600.0, 4000.0)
ZOOM_P = (7600.0, 9050.0)
WINS_FULL = [(4400., 4800.), (5000., 5600.), (6100., 6400.), (6800., 7400.)]
# Break panels are normalized on the BLUE side of the break only. With bands on
# both sides the fit splits the difference and a break-amplitude error shows up
# as a half-size offset on each side; anchoring blueward makes the red side of
# the residual read directly as the break mismatch.
WINS_B = [(3400., 3620.)]
WINS_P = [(7700., 8150.)]


def rot_kernel(dl, lam0, vsini, eps=0.6):
    """Gray (2005) rotational profile, linear limb darkening."""
    dlL = lam0 * vsini / C_KMS
    x = dl / dlL
    k = np.zeros_like(x)
    ok = np.abs(x) < 1
    k[ok] = (2 * (1 - eps) * np.sqrt(1 - x[ok] ** 2)
             + 0.5 * np.pi * eps * (1 - x[ok] ** 2)) / (np.pi * dlL * (1 - eps / 3))
    return k / k.sum()


def process_model(wm, fm, vsini, wout):
    """Rotation + native LSF + 0.1 A rebin, onto the delivered grid."""
    f = fm
    if vsini and vsini > 0:                    # rotation: constant in velocity
        dl = np.median(np.diff(wm))
        n = int(np.ceil(np.median(wm) * vsini / C_KMS / dl)) * 2 + 1
        g = (np.arange(n) - n // 2) * dl
        f = np.convolve(f, rot_kernel(g, float(np.median(wm)), vsini), mode='same')
    # native instrument profile: constant R, so fixed pixels on SYNTHE's log grid
    f = gaussian_filter1d(f, (SYNTHE_R / UVES_R) / 2.3548, mode='nearest')
    # delivered product is rebinned to 0.1 A: a boxcar, applied on a linear grid
    wl = np.arange(wm[0], wm[-1], GRID_STEP / 5.0)
    fl = uniform_filter1d(np.interp(wl, wm, f), int(round(5)), mode='nearest')
    return np.interp(wout, wl, fl)


def norm_mask(w, wins):
    keep = np.zeros_like(w, dtype=bool)
    for lo, hi in wins:
        keep |= (w > lo) & (w < hi)
    for lam in hydrogen_lines(w.min(), w.max()):
        keep &= np.abs(w - lam) > H_MASK_A
    return keep


def gap_free(f):
    return np.isfinite(f)


def make_figure(star, p):
    wo, fo, eo = load_uves(star)                       # air -> vacuum inside
    ebv = float(p['ebv']) if p['ebv'] != '' else 0.0
    vsini = float(p['vsini']) if p['vsini'] != '' else 0.0
    # UVES-POP has real coverage gaps: a dichroic gap near 5750-5844 A, a
    # 8515-8690 A gap, and growing echelle-order gaps redward of 9100 A. One
    # star (HD162678) is missing 3859-4779 A outright. Blank them with NaN
    # rather than dropping the points, so lines BREAK at the gaps instead of
    # being drawn straight across them -- which previously looked like a
    # spurious flux feature at 8500 A.
    fo = np.where(np.isfinite(fo) & (fo > 0), fo, np.nan)
    fo = deredden(wo, fo, ebv)                         # CCM89, R_V = 3.1

    wm, hnu, _ = np.loadtxt(f'models/work/{star}.spec', unpack=True)
    fm = hnu_to_flam(wm, hnu)
    lo, hi = max(FULL[0], wm[0]), min(FULL[1], wm[-1])
    m = (wo > lo) & (wo < hi)
    wo, fo = wo[m], fo[m]
    gapfrac = float(np.mean(~np.isfinite(fo)))
    mi0 = process_model(wm, fm, vsini, wo)

    cols = [((lo, hi), 'full range', None, WINS_FULL, 'upper right'),
            (ZOOM_B, f'Balmer break ({BALMER:.0f} $\\AA$)', BALMER, WINS_B, 'lower right'),
            (ZOOM_P, f'Paschen break ({PASCHEN:.0f} $\\AA$)', PASCHEN, WINS_P, 'lower left')]
    fig, axes = plt.subplots(6, 1, figsize=(11.5, 17.5),
                             gridspec_kw={'height_ratios': [2.1, 1] * 3})
    fig.patch.set_facecolor(SURFACE)
    rms = np.nan
    cover = []

    for i, (xlim, ttl, mark, wins, legloc) in enumerate(cols):
        keep = norm_mask(wo, wins) & gap_free(fo)
        if keep.sum() < 20:
            keep = gap_free(fo)
        mi = mi0 * np.median(fo[keep] / mi0[keep])
        resid = (fo - mi) / mi
        ax, rax = axes[2 * i], axes[2 * i + 1]
        for a in (ax, rax):
            a.set_facecolor(SURFACE)
            a.fill_between(wo, 0, 1, where=keep, transform=a.get_xaxis_transform(),
                           color=MUTED, alpha=.13, lw=0, zorder=0)
            for m_ in ([BALMER, PASCHEN] if mark is None else [mark]):
                a.axvline(m_, color=MUTED, ls='--', lw=1, zorder=1)
            a.set_xlim(*xlim)
            a.grid(alpha=.25, color=GRID, lw=.7)
            a.tick_params(labelsize=8, colors=MUTED)
            for sp in a.spines.values():
                sp.set_color(GRID)
        ax.fill_between(wo, 0, 1, where=keep, transform=ax.get_xaxis_transform(),
                        color=MUTED, alpha=.13, lw=0, zorder=0,
                        label='normalization windows (H masked)')
        ax.plot(wo, fo, color=OBS_C, lw=1.0, zorder=4,
                label=f'UVES-POP, dereddened E(B-V)={ebv:.3f}')
        ax.plot(wo, mi, color=MOD_C, lw=1.0, zorder=3,
                label=f'ATLAS12 + SYNTHE, v sin i={vsini:.0f} km/s')
        inw = (wo > xlim[0]) & (wo < xlim[1]) & gap_free(fo)
        ylo, yhi = np.nanpercentile(np.concatenate([fo[inw], mi[inw]]), [0.2, 99.8])
        pad = .08 * (yhi - ylo)
        ax.set_ylim(ylo - pad, yhi + pad)
        ax.set_title(f'{star} — {ttl}', fontsize=10, color=INK)
        ax.set_ylabel(r'F$_\lambda$  (locally scaled)', fontsize=9, color=INK)
        ax.legend(fontsize=7.5, loc=legloc, framealpha=.92)
        rax.axhline(0, color=MUTED, lw=1, zorder=2)
        rax.fill_between(wo, resid * 100, 0, color=OBS_C, alpha=.28, lw=0, zorder=3)
        rax.plot(wo, resid * 100, color=OBS_C, lw=.8, zorder=4)
        rax.set_xlabel(r'Wavelength [$\AA$, vacuum]', fontsize=9, color=INK)
        rax.set_ylabel('(obs - model) / model  [%]', fontsize=9, color=INK)
        rax.set_ylim(*np.nanpercentile(resid[inw] * 100, [1, 99]) + np.array([-4, 4]))
        if i == 1:
            rms = np.nanstd(resid[inw]) * 100
            cover.append(float(np.mean(gap_free(fo)[(wo > xlim[0]) & (wo < xlim[1])])))

    fig.suptitle(
        f'{star}:  Teff={float(p["teff"]):.0f} K   log g={float(p["logg"]):.2f}   '
        f'[Fe/H]={float(p["fe_h"]):+.2f}   v sin i={vsini:.0f} km/s   '
        f'E(B-V)={ebv:.3f}      Balmer-region residual RMS {rms:.1f}%'
        f'      Balmer coverage {100*cover[0]:.0f}%'
        + ('   [GAPS: {:.0f}% of range missing]'.format(100 * gapfrac)
           if gapfrac > 0.02 else ''),
        fontsize=11, color=INK)
    fig.tight_layout(rect=[0, 0, 1, 0.975])
    out = ROOT / 'figures' / f'uves_vs_model_{star}.png'
    fig.savefig(out, dpi=200, facecolor=SURFACE)
    plt.close(fig)
    print(f'  {star}: Balmer RMS={rms:.1f}%, Balmer coverage {100*cover[0]:.0f}%, '
          f'gaps {100*gapfrac:.0f}% of full range  -> {out.name}')


def main():
    rows = {r['star']: r for r in
            csv.DictReader(open(ROOT / 'data' / 'uves_pop_selected.csv'))}
    stars = [s for s in sorted(rows)
             if (ROOT / 'models' / 'work' / f'{s}.spec').exists()]
    print(f'{len(stars)} of {len(rows)} UVES-POP models ready:')
    for s in stars:
        make_figure(s, rows[s])


if __name__ == '__main__':
    main()
