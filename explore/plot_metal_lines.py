"""Panels for the most [M/H]-sensitive features, model against XSL.

One panel per feature from `data/metal_sensitivity.csv` (written by
explore/metal_sensitivity.py). Each shows:

  * the XSL observation, solid where it is actually fitted and faded where the
    mask excludes it, so an exclusion cannot be misread as missing data;
  * the model at the star's nearest grid node;
  * a shaded band spanning the grid's [M/H] range at that Teff / log g.

That last one is the point. It shows how much of the feature's variation the
grid can actually reach, so a star sitting outside the band is visibly outside
the grid rather than merely fitting badly -- which matters here, because 8 of
the 13 sample stars fall outside the grid in [M/H] and the scan can only pile
up on the boundary for them.

Each model variant gets its own calibration solve, as the fit does, so what is
being compared is line depth and not continuum placement.

    python3 explore/plot_metal_lines.py --star HD194453
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
from common.specplot import spectrum_panel, style, BAND_C, HELD_C
from common.species import (atmosphere_point, abundances, species_label,
                            dominant_species)

ROOT = Path(__file__).resolve().parent.parent
OBS_C, MOD_C, BAND_C = '#2a78d6', '#eb6834', '#7a3fa8'
SURFACE, INK, MUTED, GRIDC, HELD_C = '#fcfcfb', '#22262b', '#6b7280', '#dfe3e8', '#c0392b'


def style(ax):
    ax.set_facecolor(SURFACE)
    ax.grid(alpha=.25, color=GRIDC, lw=.7)
    ax.tick_params(labelsize=8, colors=MUTED)
    for s in ax.spines.values():
        s.set_color(GRIDC)


def sample_row(star):
    for r in csv.DictReader(open(ROOT / 'data' / 'sample.csv')):
        if r['star'] == star:
            return r
    raise KeyError(star)


def top_features(n):
    rows = list(csv.DictReader(open(ROOT / 'data' / 'metal_sensitivity.csv')))
    rows.sort(key=lambda r: float(r['depth_change']))
    return rows[:n]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--star', default='HD194453')
    ap.add_argument('--n', type=int, default=8)
    ap.add_argument('--vsini', type=float, default=0.0)
    ap.add_argument('--ebv', type=float, default=0.0)
    a = ap.parse_args()

    row, grid = sample_row(a.star), Grid()
    # 'all' windows AND drop_bad=False: features excluded from the fit because
    # the models get them wrong are exactly the ones worth LOOKING at, so the
    # prediction panels keep every one of them.
    obs = load_xsl(a.star, metals='all', drop_bad=False)
    feats = top_features(a.n)

    teff = float(grid.teff[np.argmin(np.abs(grid.teff - float(row['teff_ngsl'])))])
    logg = float(grid.logg[np.argmin(np.abs(grid.logg - float(row['logg_ngsl'])))])
    mh = float(grid.mh[np.argmin(np.abs(grid.mh - float(row['mh_ngsl'])))])
    mh_lo, mh_hi = float(grid.mh.min()), float(grid.mh.max())
    print(f'{a.star}: node Teff={teff:.0f} log g={logg:.2f} [M/H]={mh:+.2f}; '
          f'grid [M/H] spans {mh_lo:+.1f} to {mh_hi:+.1f}')
    print(f'  XSL n={obs.ndata}; {len(feats)} features')

    def model_for(z):
        th = dict(teff=teff, logg=logg, mh=z, ebv=a.ebv, vsini=a.vsini)
        cal, _ = solve(obs, predict(th, [obs], grid)[0].value)
        return cal

    cal_node, cal_lo, cal_hi = model_for(mh), model_for(mh_lo), model_for(mh_hi)

    # Dominant species per window, DERIVED from the Kurucz list with a full
    # Saha-Boltzmann weight at the model's own line-forming conditions -- not
    # assigned from memory. Ranking on log gf alone would return Co I and Nb I,
    # which have the most transitions in this range and are entirely ionised
    # away at 11,600 K.
    atm = ROOT / 'models' / 'work' / f'{a.star}.atm'
    species = {}
    if atm.exists():
        T_line, ne_line = atmosphere_point(atm)
        eps = abundances(atm)
        print(f'  line-forming point: T={T_line:.0f} K, Ne={ne_line:.2e} cm^-3')
        for f in feats:
            species[float(f['lam_center'])] = species_label(
                float(f['lam_lo']), float(f['lam_hi']), T_line, ne_line, eps)
    else:
        print(f'  no {atm.name}: species not identified')

    # Where does the observed depth sit relative to what the grid can reach?
    # Eyeballing the band is not enough: a feature the model gets wrong at EVERY
    # grid metallicity is a line-list problem, not a metallicity measurement,
    # and the two look similar in a plot.
    print(f'\n  {"species":>12}{"lambda":>9}{"obs":>9}{"[M/H]=-0.5":>12}'
          f'{"[M/H]=+0.3":>12}  verdict')
    verdicts = {}
    for f in feats:
        lo, hi = float(f['lam_lo']), float(f['lam_hi'])
        sel = obs.mask & (obs.wavelength >= lo) & (obs.wavelength <= hi)
        if sel.sum() < 3:
            continue
        # depth relative to each model's own local level, so the continuum
        # solve cannot masquerade as a depth difference
        def dep(y):
            return 1.0 - np.mean(y[sel]) / np.percentile(y[sel], 90)
        d_obs, d_lo, d_hi = dep(obs.flux), dep(cal_lo), dep(cal_hi)
        band_lo, band_hi = min(d_lo, d_hi), max(d_lo, d_hi)
        if d_obs < band_lo - 0.005:
            v = 'SHALLOWER than the grid can reach'
        elif d_obs > band_hi + 0.005:
            v = 'DEEPER than the grid can reach'
        else:
            v = 'within grid range'
        verdicts[float(f['lam_center'])] = v
        print(f'  {species.get(float(f["lam_center"]), "?"):>12}'
              f'{float(f["lam_center"]):>9.1f}{d_obs:>9.3f}{d_lo:>12.3f}'
              f'{d_hi:>12.3f}  {v}')
    n_out = sum(1 for v in verdicts.values() if 'grid can reach' in v)
    print(f'\n  {n_out}/{len(verdicts)} features fall outside the grid\'s '
          f'[M/H] range entirely')

    ncol = 2
    nrow = int(np.ceil(len(feats) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(12.0, 2.35 * nrow))
    fig.patch.set_facecolor(SURFACE)
    axes = np.atleast_1d(axes).ravel()

    for ax, f in zip(axes, feats):
        lo, hi = float(f['lam_lo']), float(f['lam_hi'])
        cen, dd = float(f['lam_center']), float(f['depth_change'])
        pad = max(6.0, 0.6 * (hi - lo))
        v = verdicts.get(cen, '')
        flag = '  \u2014 outside grid' if 'grid can reach' in v else ''
        sp = species.get(cen, '')
        drawn = spectrum_panel(
            ax, obs, [(f'model [M/H]={mh:+.1f}', cal_node)], lo - pad, hi + pad,
            title=f'{sp}   {cen:.1f} $\\AA$   ($\\Delta$depth {dd:+.3f}, '
                  f'{hi - lo:.0f} $\\AA$ wide){flag}',
            band=(lo, hi))
        if not drawn:
            continue
        ax.set_title(ax.get_title(), fontsize=9,
                     color=(HELD_C if flag else '#22262b'))
        # the band the grid can reach, which is what this figure is for
        inr = (obs.wavelength > lo - pad) & (obs.wavelength < hi + pad)
        ax.fill_between(obs.wavelength[inr], cal_lo[inr], cal_hi[inr],
                        color='#7a3fa8', alpha=.22, lw=0, zorder=0,
                        label=f'grid [M/H] {mh_lo:+.1f} to {mh_hi:+.1f}')

    for ax in axes[len(feats):]:
        ax.set_visible(False)
    axes[0].legend(fontsize=7, loc='lower left', framealpha=.92)

    fig.suptitle(
        f'{a.star} — the {len(feats)} most [M/H]-sensitive features in XSL\n'
        f'node Teff={teff:.0f} / log g={logg:.2f} / [M/H]={mh:+.2f}    '
        f'held fixed: E(B-V)={a.ebv:.3f}, v sin i={a.vsini:.0f} km/s    '
        f'shaded band = what the grid can reach',
        fontsize=11, color=INK, linespacing=1.5)
    fig.tight_layout(rect=[0, 0, 1, 1 - 0.055 * (6.0 / nrow)])
    out = ROOT / 'figures' / f'metal_lines_{a.star}.png'
    fig.savefig(out, dpi=170, facecolor=SURFACE, bbox_inches='tight')
    plt.close(fig)
    print(f'  -> {out.relative_to(ROOT)}')


if __name__ == '__main__':
    main()
