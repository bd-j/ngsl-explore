"""Rank spectral features by their sensitivity to [M/H], from the models alone.

Which lines should the XSL fit actually look at? Reputation is a poor guide at
10,000 K -- Mg I b is a standard metallicity diagnostic in cool stars and is
close to useless here, because magnesium is largely ionised. So the regions are
chosen by measurement: take two grid models differing only in [M/H], broaden
both to the instrument, normalise each by its own continuum, and rank features
by how much their DEPTH changes.

Depth rather than flux, and continuum-normalised, because the fit marginalises
XSL's continuum away: a feature that only shifts the continuum level carries no
information once that is done.

Excluded up front, each for a stated reason:

  hydrogen +/- H_EXCLUDE     the Balmer lines get their own windows; they
                             constrain Teff and log g, not [M/H]
  the held-out break window  never conditioned on (see fitting.observations)
  Ca II H and K              inside the held-out window AND carrying an
                             interstellar component on these sightlines. K is
                             the single most [M/H]-sensitive feature in the
                             whole spectrum, which is exactly what makes it
                             dangerous: an ISM line read as stellar metallicity
                             biases [M/H] in the same direction as the reddening,
                             so the error would look self-consistent.
  Na I D                     interstellar for the same reason

Each feature is also labelled with its dominant species, DERIVED from the
Kurucz line list with a full Saha-Boltzmann weight at the model atmosphere's own
line-forming conditions (common/species.py). Identifications must not be
assigned from memory: an earlier version of this analysis carried labels written
from expectation and half of them were wrong.

Writes data/metal_sensitivity.csv

    python3 explore/metal_sensitivity.py
    python3 explore/metal_sensitivity.py --check-stability
"""
import argparse
import csv
import sys
from pathlib import Path

import numpy as np
from scipy.ndimage import maximum_filter1d, uniform_filter1d

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from fitting.model import Grid
from fitting.observations import BREAK_WINDOW, xsl_resolution_segments
from common.lines import hydrogen_lines, ISM_LINES
from common.species import atmosphere_point, abundances, species_label

ROOT = Path(__file__).resolve().parent.parent
XSL_RANGE = (3501.0, 9500.0)
H_EXCLUDE = 25.0            # A, half-width around every hydrogen line
ISM_EXCLUDE = 6.0           # A, half-width around interstellar lines
MIN_DEPTH_CHANGE = 0.02     # keep features above this
MERGE_GAP = 3.0             # A; nearer than this counts as one feature
CONT_WINDOW = 0.004         # d(ln lambda) for the running-max continuum


def normalised(w, f):
    """Continuum-normalise with a running upper envelope.

    A running max over a window wider than any line, then smoothed. Not exact,
    but both metallicities get the identical treatment, so what survives in the
    difference is the line depth change rather than the envelope's own error.
    """
    step = float(np.median(np.diff(np.log(w))))
    n = max(3, int(CONT_WINDOW / step))
    return f / uniform_filter1d(maximum_filter1d(f, n), n)


def sensitivity(grid, teff, logg, mh_lo, mh_hi, resolution):
    """-> (wave, depth_change) with excluded regions zeroed."""
    from common.lsf import broaden_R
    w = grid.wave
    s = (w >= XSL_RANGE[0]) & (w <= XSL_RANGE[1])
    w = w[s]

    def prep(mh):
        f = grid.interp(teff, logg, mh)[s]
        # one constant-R kernel per arm, as the XSL observation is treated
        out = np.array(f, float)
        for lo, hi, R in resolution:
            seg = (w >= lo) & (w < hi)
            if seg.any():
                out[seg] = broaden_R(w, f, float(R))[seg]
        return normalised(w, out)

    d = prep(mh_hi) - prep(mh_lo)        # negative where the metal-rich is deeper

    keep = np.ones_like(w, bool)
    for lam in hydrogen_lines(XSL_RANGE[0] - 200, XSL_RANGE[1] + 200, series=(2, 3)):
        keep &= np.abs(w - lam) > H_EXCLUDE
    keep &= ~((w >= BREAK_WINDOW[0]) & (w <= BREAK_WINDOW[1]))
    for lam in ISM_LINES.values():
        keep &= np.abs(w - lam) > ISM_EXCLUDE
    return w, np.where(keep, d, 0.0)


def features(w, d, nmax=4000):
    """Group the most sensitive pixels into contiguous features."""
    idx = np.argsort(d)[:nmax]
    out = []
    for i in sorted(idx):
        if out and w[i] - out[-1][1] < MERGE_GAP:
            out[-1][1] = w[i]
            out[-1][2] = min(out[-1][2], d[i])
        else:
            out.append([w[i], w[i], d[i]])
    out = [f for f in out if abs(f[2]) >= MIN_DEPTH_CHANGE]
    out.sort(key=lambda f: f[2])
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--teff', type=float, default=10200.)
    ap.add_argument('--logg', type=float, default=3.8)
    ap.add_argument('--mh-lo', type=float, default=-0.5)
    ap.add_argument('--mh-hi', type=float, default=0.3)
    ap.add_argument('--top', type=int, default=15)
    ap.add_argument('--check-stability', action='store_true',
                    help='repeat over the Teff/log g grid and report rank stability')
    ap.add_argument('--union', action='store_true',
                    help='write the UNION of the top features across Teff, which '
                         'is what the XSL windows should use')
    ap.add_argument('--union-top', type=int, default=40)
    ap.add_argument('--pad', type=float, default=2.0,
                    help='A added to each side of a feature to make a window')
    ap.add_argument('--atm', default='models/work/HD194453.atm',
                    help='model atmosphere supplying T, N_e and abundances for '
                         'the species identification')
    a = ap.parse_args()

    grid = Grid()
    res = xsl_resolution_segments()
    w, d = sensitivity(grid, a.teff, a.logg, a.mh_lo, a.mh_hi, res)
    feat = features(w, d)

    atm = ROOT / a.atm
    ident = {}
    if atm.exists():
        T_line, ne_line = atmosphere_point(atm)
        eps = abundances(atm)
        print(f'species from {atm.name}: T={T_line:.0f} K, '
              f'Ne={ne_line:.2e} cm^-3 at tau_5000 = 2/3')
        for f in feat:
            ident[f[0]] = species_label(f[0], f[1], T_line, ne_line, eps)
    else:
        print(f'{atm} not found: species not identified')

    print(f'[M/H] sensitivity at {a.teff:.0f} K / log g {a.logg}, '
          f'{a.mh_lo:+.1f} vs {a.mh_hi:+.1f}, XSL resolution')
    print(f'excluded: H +/-{H_EXCLUDE:.0f} A, break window '
          f'{BREAK_WINDOW[0]:.0f}-{BREAK_WINDOW[1]:.0f} A, '
          f'{", ".join(ISM_LINES)} (ISM)\n')
    print(f'  {"lambda":>9} {"width":>7} {"d(depth)":>9}  species')
    for f in feat[:a.top]:
        print(f'  {0.5 * (f[0] + f[1]):>9.1f} {f[1] - f[0]:>7.1f} {f[2]:>9.3f}'
              f'  {ident.get(f[0], "?")}')
    blue = sum(1 for f in feat if f[0] < 4600)
    print(f'\n  {len(feat)} features above {MIN_DEPTH_CHANGE}; '
          f'{blue} blueward of 4600 A')

    with open(ROOT / 'data' / 'metal_sensitivity.csv', 'w', newline='') as fh:
        wr = csv.writer(fh)
        wr.writerow(['lam_center', 'lam_lo', 'lam_hi', 'width', 'depth_change',
                     'species'])
        for f in feat:
            wr.writerow([f'{0.5 * (f[0] + f[1]):.2f}', f'{f[0]:.2f}',
                         f'{f[1]:.2f}', f'{f[1] - f[0]:.2f}', f'{f[2]:.4f}',
                         ident.get(f[0], '')])
    print(f'  -> data/metal_sensitivity.csv ({len(feat)} rows)')

    if a.check_stability:
        check_stability(grid, res, a, feat)
    if a.union:
        union_windows(grid, res, a)


def check_stability(grid, res, a, ref):
    """Does the ranking hold across the grid's Teff / log g range?

    Worth testing rather than assuming: if the ranking moved a lot, the XSL
    metal windows would have to be chosen per star instead of once.
    """
    from scipy.stats import spearmanr
    ref_top = [0.5 * (f[0] + f[1]) for f in ref[:30]]
    print(f'\nRank stability of the top 30 features across the grid:')
    print(f'  {"Teff":>7}{"logg":>7}{"shared":>9}{"spearman":>10}')
    for teff in (9000., 10200., 11000.):
        for logg in (3.2, 3.8, 4.6):
            w2, d2 = sensitivity(grid, teff, logg, a.mh_lo, a.mh_hi, res)
            f2 = features(w2, d2)
            cen2 = [0.5 * (f[0] + f[1]) for f in f2[:30]]
            # match by wavelength proximity, then correlate the two orderings
            pairs = []
            for i, c1 in enumerate(ref_top):
                j = int(np.argmin([abs(c1 - c2) for c2 in cen2])) if cen2 else -1
                if j >= 0 and abs(c1 - cen2[j]) < MERGE_GAP:
                    pairs.append((i, j))
            rho = (spearmanr([p[0] for p in pairs], [p[1] for p in pairs]).statistic
                   if len(pairs) > 3 else np.nan)
            print(f'  {teff:>7.0f}{logg:>7.1f}{len(pairs):>6}/30{rho:>10.2f}')
    print('\n  shared = how many of the reference top 30 are still in the top 30;')
    print('  spearman = whether those keep the same relative order.')




def union_windows(grid, res, a):
    """Merge the top features from several Teff into one window list.

    The ranking is NOT stable in Teff: at 9000 K only 11-17 of the reference
    top 30 survive and their order is uncorrelated (Spearman -0.04 to 0.22),
    because the ionisation balance shifts. The sample spans 8759-10885 K, so a
    window set optimised at one temperature is wrong for the ends of it.

    Taking the union costs almost nothing -- a window where the line happens to
    be weak simply contributes little to the likelihood -- whereas omitting a
    window that matters at 9000 K loses that star's metallicity constraint
    outright. Asymmetric risk, so take the union.
    """
    spans = []
    for teff in (9000., 10000., 11000.):
        w, d = sensitivity(grid, teff, a.logg, a.mh_lo, a.mh_hi, res)
        for f in features(w, d)[:a.union_top]:
            spans.append((f[0] - a.pad, f[1] + a.pad))
    spans.sort()
    merged = []
    for lo, hi in spans:
        if merged and lo <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], hi)
        else:
            merged.append([lo, hi])

    total = sum(hi - lo for lo, hi in merged)
    print(f'\nUnion of the top {a.union_top} features at 9000 / 10000 / 11000 K,'
          f' padded {a.pad:.0f} A:')
    print(f'  {len(merged)} windows, {total:.0f} A total')
    for lo, hi in merged:
        print(f'    {lo:8.1f} - {hi:8.1f}   ({hi - lo:5.1f} A)')
    out = ROOT / 'data' / 'xsl_metal_windows.csv'
    with open(out, 'w', newline='') as fh:
        wr = csv.writer(fh)
        wr.writerow(['lo', 'hi'])
        for lo, hi in merged:
            wr.writerow([f'{lo:.2f}', f'{hi:.2f}'])
    print(f'  -> {out.relative_to(ROOT)}')

if __name__ == '__main__':
    main()
