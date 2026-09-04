"""Compare ATLAS12 models against observed NGSL spectra at the hydrogen breaks.

Three things this figure gets right, each of which was wrong at some point:

1. The smoothing kernel is wavelength dependent. The NGSL LSF is set by a fixed
   dispersion per grating, so it is constant in ANGSTROMS within a grating and
   jumps at the splices (G430L 3.85 A -> G750L 8.09 A). SYNTHE's grid is
   logarithmic, so a fixed sigma in pixels would instead be a constant-R kernel.

2. Flux scaling excludes hydrogen. Model and observation are matched by one
   multiplicative factor measured in continuum windows with every Balmer and
   Paschen line masked, so the lines under test never set the normalization.

3. Each break panel is renormalized locally, in windows bracketing that break,
   so a small continuum tilt far away cannot masquerade as a break mismatch.
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

sys.path.insert(0, str(Path(__file__).resolve().parent))
from balmer_metric import balmer_discontinuity
from make_atlas_model import hnu_to_flam
from ngsl_wavecal import apply_wavecal, load_table

WAVECAL = load_table()

ROOT = Path(__file__).resolve().parent.parent
BALMER = 3646.0                      # 911.7635 * 4, vacuum
PASCHEN = 8205.9                     # 911.7635 * 9, vacuum
FULL = (3200.0, 9400.0)
ZOOM_B = (3600.0, 4000.0)
ZOOM_P = (7600.0, 9100.0)
H_MASK_A = 20.0

# Normalization windows. The full-range panel uses broad continuum; each break
# panel renormalizes in windows bracketing that break.
WINS_FULL = [(4400., 4800.), (5000., 5600.), (6100., 6400.), (6800., 7400.)]
# Break panels are normalized on the BLUE side of the break only. With bands on
# both sides the fit splits the difference and a break-amplitude error shows up
# as a half-size offset on each side; anchoring blueward makes the red side of
# the residual read directly as the break mismatch.
WINS_B = [(3400., 3620.)]
WINS_P = [(7700., 8150.)]

# NGSL LSF FWHM in Angstroms per grating (FWHM_px x dispersion, from
# data/stis_lsf_resolution.csv). G230LB has no published LSF: 2-px lower bound.
NGSL_LSF = [(1675., 3058., 2.75), (3058., 5647., 3.85), (5647., 10198., 8.09)]

OBS_C, MOD_C, SURFACE = '#2a78d6', '#eb6834', '#fcfcfb'
INK, MUTED, GRID = '#22262b', '#6b7280', '#dfe3e8'


def hydrogen_lines(wmin, wmax, series=(2, 3), nmax=40):
    out = []
    for m in series:
        for n in range(m + 1, nmax):
            lam = 911.7635 / (1.0 / m ** 2 - 1.0 / n ** 2)
            if wmin <= lam <= wmax:
                out.append(lam)
    return np.array(out)


def broaden(w, f, fwhm_A, step=0.01):
    """Gaussian of constant FWHM in ANGSTROMS (resample to linear grid first)."""
    wl = np.arange(w[0], w[-1], step)
    sm = gaussian_filter1d(np.interp(wl, w, f), fwhm_A / 2.3548 / step,
                           mode='nearest')
    return np.interp(w, wl, sm)


def broaden_ngsl(w, f):
    """The NGSL instrument profile: piecewise constant in Angstroms per grating."""
    out = np.array(f, dtype=float)
    for lo, hi, fwhm in NGSL_LSF:
        seg = (w >= lo) & (w < hi)
        if seg.any():
            out[seg] = broaden(w, f, fwhm)[seg]
    return out


def norm_mask(wo, wins):
    """Points inside the given windows and clear of every hydrogen line."""
    keep = np.zeros_like(wo, dtype=bool)
    for lo, hi in wins:
        keep |= (wo > lo) & (wo < hi)
    for lam in hydrogen_lines(wo.min(), wo.max()):
        keep &= np.abs(wo - lam) > H_MASK_A
    return keep


def scaled(wo, fo, mi, wins):
    """Model scaled to the observation in the given windows -> (model, resid, mask)."""
    keep = norm_mask(wo, wins)
    if keep.sum() < 20:
        keep = np.ones_like(wo, dtype=bool)
    m = mi * np.median(fo[keep] / mi[keep])
    return m, (fo - m) / m, keep


def shade(ax, wo, keep, label=None):
    """Shade the normalization windows across the full height of an axis."""
    ax.fill_between(wo, 0, 1, where=keep, transform=ax.get_xaxis_transform(),
                    color=MUTED, alpha=.13, lw=0, zorder=0, label=label)


def style(ax):
    ax.set_facecolor(SURFACE)
    ax.grid(alpha=.25, color=GRID, lw=.7)
    ax.tick_params(labelsize=8, colors=MUTED)
    for s in ax.spines.values():
        s.set_color(GRID)


def make_figure(star, cat):
    d = fits.getdata('data/spectra/' + cat[star]['file'])
    # NGSL is delivered in AIR with a per-grating zero point (no wavecals were
    # taken with the stellar exposures). Convert and correct before comparing.
    wo_all = apply_wavecal(d['WAVELENGTH'].astype(float), star, WAVECAL)
    fo_all = d['FLUX'].astype(float)
    wm, hnu, _ = np.loadtxt(f'models/work/{star}.spec', unpack=True)
    fm = hnu_to_flam(wm, hnu)

    lo, hi = max(FULL[0], wm[0]), min(FULL[1], wm[-1])
    narrow = (lo > FULL[0] + 1) or (hi < FULL[1] - 1)
    m = (wo_all > lo) & (wo_all < hi) & np.isfinite(fo_all) & (fo_all > 0)
    wo, fo = wo_all[m], fo_all[m]
    mi0 = np.interp(wo, wm, broaden_ngsl(wm, fm))

    cols = [((lo, hi), 'full range', None, WINS_FULL, 'upper right'),
            (ZOOM_B, f'Balmer break ({BALMER:.0f} $\\AA$)', BALMER, WINS_B, 'lower right'),
            (ZOOM_P, f'Paschen break ({PASCHEN:.0f} $\\AA$)', PASCHEN, WINS_P, 'lower left')]

    fig, axes = plt.subplots(6, 1, figsize=(11.5, 17.5),
                             gridspec_kw={'height_ratios': [2.1, 1] * 3})
    fig.patch.set_facecolor(SURFACE)
    d_obs = d_mod = None

    for i, (xlim, ttl, mark, wins, legloc) in enumerate(cols):
        mi, resid, keep = scaled(wo, fo, mi0, wins)
        ax, rax = axes[2 * i], axes[2 * i + 1]
        for a in (ax, rax):
            style(a)
            shade(a, wo, keep)
            for m_ in ([BALMER, PASCHEN] if mark is None else [mark]):
                a.axvline(m_, color=MUTED, ls='--', lw=1, zorder=1)
            a.set_xlim(*xlim)
        shade(ax, wo, keep, label='normalization windows (H masked)')
        ax.plot(wo, fo, color=OBS_C, lw=1.3, zorder=4, label='NGSL observed')
        ax.plot(wo, mi, color=MOD_C, lw=1.3, zorder=3,
                label='ATLAS12 + SYNTHE, NGSL LSF')
        ax.set_title(f'{star} — {ttl}', fontsize=10, color=INK)
        ax.set_ylabel(r'F$_\lambda$  (locally scaled)', fontsize=9, color=INK)
        # Scale y to the data actually inside this panel's wavelength range,
        # not to the whole spectrum.
        inwin = (wo > xlim[0]) & (wo < xlim[1])
        ylo, yhi = np.percentile(np.concatenate([fo[inwin], mi[inwin]]), [0.2, 99.8])
        pad = 0.08 * (yhi - ylo)
        ax.set_ylim(ylo - pad, yhi + pad)
        ax.legend(fontsize=7.5, loc=legloc, framealpha=.92)
        rax.axhline(0, color=MUTED, lw=1, zorder=2)
        rax.fill_between(wo, resid * 100, 0, color=OBS_C, alpha=.28, lw=0, zorder=3)
        rax.plot(wo, resid * 100, color=OBS_C, lw=.9, zorder=4)
        rax.set_xlabel(r'Wavelength [$\AA$, vacuum]', fontsize=9, color=INK)
        rax.set_ylabel('(obs - model) / model  [%]', fontsize=9, color=INK)
        sel = (wo > xlim[0]) & (wo < xlim[1])
        rax.set_ylim(*np.percentile(resid[sel] * 100, [1, 99]) + np.array([-4, 4]))
        if i == 1:
            d_obs, _, _ = balmer_discontinuity(wo, fo)
            d_mod, _, _ = balmer_discontinuity(wo, mi)
            rms = np.std(resid[sel]) * 100

    r = cat[star]
    def _wc(g):
        c = WAVECAL[(star, g)]
        return (f'{g} {c["a"]:+.2f}' +
                (f'{c["b"]*1000:+.2f}/kA' if c['b'] else ''))
    wc = ', '.join(_wc(g) for g in ('G430L', 'G750L') if (star, g) in WAVECAL)
    fig.suptitle(
        f'{star}:  Teff={float(r["teff_ngsl_K"]):.0f} K   log g={r["logg_ngsl"]}   '
        f'[M/H]={float(r["m_h_ngsl"]):+.1f}      '
        f'D$_{{Balmer}}$  obs {d_obs:.3f} / model {d_mod:.3f} ({d_mod-d_obs:+.3f})'
        f'      Balmer-region residual RMS {rms:.1f}%'
        f'      wavecal: air$\\to$vac, {wc} $\\AA$'
        + (f'   [model only {lo:.0f}-{hi:.0f} $\\AA$]' if narrow else ''),
        fontsize=11, color=INK)
    fig.tight_layout(rect=[0, 0, 1, 0.975])
    out = ROOT / 'figures' / f'model_vs_obs_{star}.png'
    fig.savefig(out, dpi=200, facecolor=SURFACE)
    plt.close(fig)
    print(f'  {star}: D obs={d_obs:.3f} model={d_mod:.3f} ({d_mod-d_obs:+.3f}), '
          f'Balmer RMS={rms:.1f}%  -> {out.name}')


def main():
    all_rows = list(csv.DictReader(open(ROOT / 'data' / 'balmer_candidates.csv')))
    cat = {r['star']: r for r in all_rows}
    dropped = [(r['star'], r['reject_reason']) for r in all_rows
               if r.get('selected') == 'no']
    for star, why in dropped:
        print(f'  skipping {star}: {why}')
    stars = [r['star'] for r in all_rows if r.get('selected') != 'no'
             and (ROOT / 'models' / 'work' / f'{r["star"]}.spec').exists()]
    print(f'{len(stars)} model(s) available:')
    for s in stars:
        make_figure(s, cat)


if __name__ == '__main__':
    main()
