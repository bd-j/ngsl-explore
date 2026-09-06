"""Build ATLAS12 model atmospheres and R=10000 spectra for the Balmer-break sample.

Pipeline per star:
  1. ATLAS12  - converge a model atmosphere at (Teff, log g, [Fe/H])
  2. SYNTHE   - synthesize at R = 300000 (the code's enforced floor)
  3. smooth   - Gaussian-convolve to R = 10000 and convert to f_lambda

Why synthesis runs at R=300000 and not 10000: SYNTHE's `resolu` sets the
COMPUTATION grid, not an instrumental profile, and the code hard-refuses
anything below 300000 (RESOLU_MIN in synthe.f90). A coarse grid undersamples
line cores and then integrates them, biasing the smoothed spectrum ~10% too
absorbed. The observed resolution is applied afterwards, here.

Metallicity: ATLAS12's `zscale` is a LINEAR multiplicative factor on all
Z >= 3 (it sets XRELATIVE(3:99) = log10(zscale)), so a solar-scaled [Fe/H]
is zscale = 10**[Fe/H].

Usage:
  export ATLAS12=/path/to/atlas12
  python3 explore/make_atlas_model.py --all
  python3 explore/make_atlas_model.py --star HD194453 --teff 10241 --logg 3.9 --feh 0.0
"""
import argparse
import csv
import re
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

C_ANG = 2.99792458e18                     # c in Angstrom/s
SYNTHE_R = 300000.0                       # SYNTHE floor; also its native grid
ROOT = Path(__file__).resolve().parent.parent
MODELS = ROOT / 'models'
# Default starting atmosphere: the solar model shipped in atlas12/workdir.
# It is a poor start for A stars -- SCALE_MODEL rescales T linearly by the Teff
# ratio, and from 5777 K the deep layers of a ~10000 K model start so far off
# that the temperature correction diverges (layer 80 ran 12000 -> 120000 K over
# 8 iterations, then spun). Point --start at a model near the target Teff.
START_ATM = 'ap00t5777g4.44at12.dat'

# C3K v2.3 starting grid: 493 converged ATLAS12 models, Teff 2500-50000 K,
# log g -1.0 to 5.0, at [Fe/H]=0. Near 10000 K it is spaced 250 K by 0.5 dex,
# so every Balmer-break target starts within ~125 K and 0.2 dex of a converged
# model of the right structure. Filenames encode the parameters as t#####g#.##.
GRID_DIR = ('grids/c3k_v2.3/at12_feh+0.00_afe+0.0/atm')
GRID_RE = re.compile(r'_t(\d{5})g(-?\d\.\d{2})\.atm$')


def nearest_grid_model(griddir, teff, logg):
    """Closest converged grid model to (teff, logg).

    Distance is normalized by the grid spacing near the A stars (250 K, 0.5 dex)
    so neither axis dominates; a 250 K error and a 0.5 dex error count the same.
    """
    best, bestd = None, None
    for f in sorted(Path(griddir).glob('*.atm')):
        m = GRID_RE.search(f.name)
        if not m:
            continue
        t, g = float(m.group(1)), float(m.group(2))
        d = ((t - teff) / 250.0) ** 2 + ((g - logg) / 0.5) ** 2
        if bestd is None or d < bestd:
            best, bestd, bt, bg = f, d, t, g
    if best is None:
        sys.exit(f'ERROR: no parseable .atm models in {griddir}')
    return best, bt, bg


def atlas_home():
    h = os.environ.get('ATLAS12')
    if not h:
        sys.exit('ERROR: set $ATLAS12 to the atlas12 checkout root')
    return Path(h)


def run(cmd, cwd, logfile):
    """Run a command, streaming stdout+stderr to logfile. Returns elapsed seconds."""
    t0 = time.time()
    with open(logfile, 'w') as fh:
        p = subprocess.run(cmd, cwd=cwd, stdout=fh, stderr=subprocess.STDOUT)
    if p.returncode != 0:
        print(f'    FAILED (exit {p.returncode}) - see {logfile}')
        return None
    return time.time() - t0


def teff_ladder(t_start, t_target, max_step=2500.0):
    """Optional intermediate Teff rungs (--ladder).

    Passing teff= makes ATLAS12 regrid through SCALE_MODEL, which rescales T
    linearly by the Teff ratio - a rough first guess, but the iteration
    recovers from it, so a direct jump normally works and just takes its time.
    Rungs are only a fallback for a run that will not converge. When stepping
    up to a hot target, keep rungs above 8000 K where possible: below that
    ATLAS12 reads the TiO and H2O lists (~4 GB) for no benefit here."""
    if abs(t_target - t_start) <= max_step:
        return [t_target]
    n = int(np.ceil(abs(t_target - t_start) / max_step))
    return list(np.linspace(t_start, t_target, n + 1)[1:])


def make_atmosphere(star, teff, logg, feh, vturb, numit, ladder, workdir, A12,
                    start_atm=START_ATM):
    """Converge an ATLAS12 atmosphere. Returns path to the .atm file."""
    zscale = 10.0 ** feh
    start = workdir / Path(start_atm).name
    if not start.exists():
        src = Path(start_atm)
        if not src.is_absolute():
            src = A12 / 'workdir' / start_atm
        if not src.exists():
            sys.exit(f'ERROR: starting atmosphere not found: {src}')
        shutil.copy(src, start)

    steps = teff_ladder(5777.0, teff) if ladder else [teff]
    current = start.name
    for i, t in enumerate(steps):
        last = (i == len(steps) - 1)
        base = f'{star}' if last else f'{star}_step{i}'
        # log g is walked in step with Teff; only the final step must land exactly
        g = logg if last else 4.44 + (logg - 4.44) * (i + 1) / len(steps)
        # Intermediate rungs only need a usable structure to hand upward,
        # not convergence; the final rung gets the full iteration budget.
        nit = numit if last else max(10, numit // 2)
        cmd = [str(A12 / 'bin' / 'atlas12c.exe'), current, base,
               f'numit={nit}', f'vturb={vturb}',
               f'teff={t:.0f}', f'logg={g:.2f}', f'zscale={zscale:.6f}']
        print(f'    ATLAS12 step {i+1}/{len(steps)}: Teff={t:.0f} logg={g:.2f}'
              f' zscale={zscale:.4f} numit={nit} ... ', end='', flush=True)
        dt = run(cmd, workdir, workdir / f'{base}.atlas.log')
        if dt is None:
            return None
        print(f'{dt/60:.1f} min')
        current = f'{base}.atm'
    return workdir / f'{star}.atm'


def make_spectrum(star, atm, wlbeg, wlend, workdir, A12):
    """Synthesize at the SYNTHE floor resolution. Returns path to the .spec file."""
    cmd = [str(A12 / 'bin' / 'synthe.exe'), atm.name,
           f'wlbeg={wlbeg}', f'wlend={wlend}', f'resolu={SYNTHE_R:.0f}']
    print(f'    SYNTHE {wlbeg}-{wlend} nm at R={SYNTHE_R:.0f} ... ', end='', flush=True)
    dt = run(cmd, workdir, workdir / f'{star}.synthe.log')
    if dt is None:
        return None
    print(f'{dt/60:.1f} min')
    return workdir / f'{star}.spec'


def smooth_to_R(flux, model_R, obs_R):
    """Gaussian-smooth a constant-R log-lambda spectrum to resolving power obs_R.
    Same convention as atlas12/tools/mann_lib.py."""
    from scipy.ndimage import gaussian_filter1d
    return gaussian_filter1d(flux, model_R / (obs_R * 2.3548), mode='nearest')


def hnu_to_flam(wl_A, hnu):
    """Kurucz Eddington flux H_nu -> surface f_lambda [erg/s/cm^2/A]."""
    return 4.0 * np.pi * hnu * C_ANG / wl_A ** 2


def postprocess(star, spec, obs_R, outdir):
    """Convolve to obs_R, convert to f_lambda, write a plain CSV."""
    w, hnu, hcont = np.loadtxt(spec, unpack=True)
    flam, cont = hnu_to_flam(w, hnu), hnu_to_flam(w, hcont)
    out = outdir / f'{star}_R{obs_R:.0f}.csv'
    np.savetxt(out, np.column_stack([
        w, smooth_to_R(flam, SYNTHE_R, obs_R), smooth_to_R(cont, SYNTHE_R, obs_R)]),
        delimiter=',', header='wavelength_A_vacuum,flam,flam_continuum',
        comments='', fmt=['%.4f', '%.6e', '%.6e'])
    print(f'    -> {out.name}  ({len(w)} points, {w[0]:.1f}-{w[-1]:.1f} A)')
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--all', action='store_true',
                    help='run every star in data/balmer_candidates.csv')
    ap.add_argument('--star'); ap.add_argument('--teff', type=float)
    ap.add_argument('--logg', type=float); ap.add_argument('--feh', type=float)
    ap.add_argument('--start',
                    help='explicit starting atmosphere (path, or filename in '
                         '$ATLAS12/workdir). Overrides --grid.')
    ap.add_argument('--grid', default=GRID_DIR,
                    help='directory of converged models to start from; the one '
                         'nearest (Teff, log g) is chosen per star. Relative '
                         'paths resolve under $ATLAS12. Default: C3K v2.3.')
    ap.add_argument('--vturb', type=float, default=2.0, help='km/s (default 2)')
    ap.add_argument('--numit', type=int, default=30)
    # Wide enough to carry both hydrogen breaks: Balmer 3646 A and
    # Paschen 8206 A (vacuum), with continuum either side of each.
    ap.add_argument('--wlbeg', type=float, default=320.0, help='nm (default 320)')
    ap.add_argument('--wlend', type=float, default=950.0, help='nm (default 950)')
    ap.add_argument('--obs-R', type=float, default=10000.0)
    ap.add_argument('--ladder', action='store_true',
                    help='walk up to Teff in bounded rungs. Not normally needed: '
                         'the direct teff= jump converges, it just takes 0.5-1 h. '
                         'Reach for this only if a direct run fails to converge.')
    ap.add_argument('--workdir',
                    help='directory for atmospheres and spectra '
                         '(default models/work)')
    ap.add_argument('--skip-atlas', action='store_true',
                    help='reuse an existing .atm and only redo synthesis')
    args = ap.parse_args()

    A12 = atlas_home()
    workdir = Path(args.workdir) if args.workdir else MODELS / 'work'
    workdir.mkdir(parents=True, exist_ok=True)

    if args.all:
        # honour the selection flag set by candidate_table.py (reddening cut etc.)
        targets = [(r['star'], float(r['teff_ngsl_K']), float(r['logg_ngsl']),
                    float(r['m_h_ngsl']))
                   for r in csv.DictReader(open(ROOT / 'data' / 'balmer_candidates.csv'))
                   if r.get('selected') != 'no']
    elif args.star and args.teff and args.logg is not None and args.feh is not None:
        targets = [(args.star, args.teff, args.logg, args.feh)]
    else:
        ap.error('give --all, or --star with --teff --logg --feh')

    print(f'{len(targets)} model(s); ATLAS12={A12}\n')
    for star, teff, logg, feh in targets:
        print(f'{star}: Teff={teff:.0f} logg={logg} [Fe/H]={feh:+.2f}')
        if args.start:
            start_atm, note = args.start, 'explicit'
        else:
            gd = Path(args.grid)
            if not gd.is_absolute():
                gd = A12 / gd
            gm, gt, gg = nearest_grid_model(gd, teff, logg)
            start_atm = str(gm)
            note = f'grid t{gt:.0f} g{gg:.2f} (d={teff-gt:+.0f} K, {logg-gg:+.2f} dex)'
        print(f'    start: {Path(start_atm).name}  [{note}]')
        atm = workdir / f'{star}.atm'
        if not (args.skip_atlas and atm.exists()):
            atm = make_atmosphere(star, teff, logg, feh, args.vturb, args.numit,
                                  args.ladder, workdir, A12, start_atm)
            if atm is None:
                continue
        spec = make_spectrum(star, atm, args.wlbeg, args.wlend, workdir, A12)
        if spec is None:
            continue
        postprocess(star, spec, args.obs_R, MODELS)
        print()


if __name__ == '__main__':
    main()
