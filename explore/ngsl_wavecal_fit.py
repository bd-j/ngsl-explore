"""Fit NGSL's per-grating wavelength residual against the model grid.

NGSL took no wavecal exposures with the stellar spectra -- the delivery readme
states zero points were derived per spectrum from the positions of strong
stellar features -- so every star carries its own wavelength error and it has to
be measured per star. This script measures it; `common.ngsl_wavecal.apply_wavecal`
applies it. Same split as `explore/ngsl_lsf.py` measuring what `common/lsf.py`
applies.

Writes data/ngsl_wavecal.csv

WHAT IS MEASURED. For each grating's sub-windows, the velocity that best aligns
the model to the observation, converted to Angstroms at the window's mean
wavelength, then fitted per grating as `a + b*(lambda - 4000)`. G430L gets the
linear term because its residual is a clean monotonic ramp -- roughly +0.8 A at
the Balmer break falling through zero by ~5000 A -- and a constant fitted over
3300-5600 A leaves ~1 A uncorrected exactly where the break is. G750L scatters
window to window without a trend, so it gets a robust constant.

SIGN. A positive velocity moves the MODEL redward to meet the data, so a
positive result means the observation sits redward of truth and `apply_wavecal`
SUBTRACTS the fitted s(lambda). `--selftest` checks this against an injected
shift rather than leaving it to be reasoned about.

MEASURED ON AIR->VACUUM WAVELENGTHS ONLY. The fit deliberately does NOT use
`load_ngsl`, which applies the table this script writes: measuring through it
would return the residual left after the previous correction, and refitting
would report ~0 for any star already in the table while leaving the stars
missing from it uncorrected. The observation here is built from the FITS with
`air_to_vac` and nothing else, so `a` and `b` are absolute.

WHICH MODEL. The grid at the star's node-scan ML parameters, with its E(B-V) and
v sin i, broadened with the measured NGSL profile and INTEGRATED onto the
detector pixels -- the same forward model the fitter uses, via `predict`. The
node comes from a scan run under the previous table, which is safe because the
quantity that selects the node is insensitive to this shift: the NGSL leg is 13
broad tophat bands, whose chi2/N moves by <2% across the full range of offsets
measured here, and the XSL leg never sees NGSL's wavelength solution at all.
`--node-sensitivity` re-fits at the catalog node to show the coefficients do not
depend on that choice.
"""
import argparse
import copy
import csv
import sys
from pathlib import Path

import numpy as np
from astropy.io import fits

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common.lines import air_to_vac
from common.ngsl_wavecal import SEGMENTS, LAM_REF
from common.sample import sample_row, all_stars
from fitting.model import Grid
from fitting.observations import (Observation, NGSL_SNR_CEILING,
                                  MODEL_RANGE)
from common.lsf import NGSL_MOFFAT_BETA
from fitting.predict import predict
from fitting.calibration import solve, usable
from fitting.scan import best_node

ROOT = Path(__file__).resolve().parent.parent
C_KMS = 2.99792458e5
RV_GRID = np.arange(-120.0, 120.01, 3.0)    # +-1.5 A at the break
MIN_PIX = 25                                # usable pixels for a window to count


def ngsl_file(star):
    """-> the star's NGSL filename, from the sample or the full catalog."""
    try:
        return sample_row(star)['ngsl_file']
    except KeyError:
        pass
    with open(ROOT / 'data' / 'ngsl_catalog.csv') as fh:
        for r in csv.DictReader(fh):
            if r['target'] == star:
                return r['file']
    raise KeyError(f'no NGSL spectrum listed for {star}')


def raw_observation(star):
    """The NGSL spectrum on AIR->VACUUM wavelengths, with no wavecal applied.

    Built from the FITS rather than through `load_ngsl`, for two reasons: that
    function applies the very table this script writes, and it requires a row in
    sample.csv, which the catalog-only stars do not have. Everything else -- the
    S/N ceiling, the usable-pixel mask, the instrument profile -- matches it.
    """
    d = fits.getdata(ROOT / 'data' / 'spectra' / ngsl_file(star))
    w = air_to_vac(d['WAVELENGTH'].astype(float))
    f = d['FLUX'].astype(float)
    e = np.maximum(d['STATERR'].astype(float), f / NGSL_SNR_CEILING)
    ok = (np.isfinite(w) & np.isfinite(f) & np.isfinite(e) & (f > 0) & (e > 0)
          & (w >= MODEL_RANGE[0]) & (w <= MODEL_RANGE[1]))
    return Observation(name='ngsl', star=star, wavelength=w, flux=f,
                       uncertainty=e, mask=ok,
                       resolution=('ngsl', NGSL_MOFFAT_BETA),
                       calibration=('scalar',), rv_fixed=None)


def _window(obs, lo, hi):
    o = copy.copy(obs)
    o.mask = obs.mask & (obs.wavelength >= lo) & (obs.wavelength <= hi)
    return o


def _chi2(obs, grid, theta):
    p = predict(theta, [obs], grid)[0]
    cal, _ = solve(obs, p.value)
    u = usable(obs, cal)
    if u.sum() < MIN_PIX:
        return np.nan, 0
    r = (obs.flux[u] - cal[u]) / obs.uncertainty[u]
    return float(np.sum(r * r)), int(u.sum())


def window_offset(obs, grid, theta, lo, hi, rvs=RV_GRID):
    """-> (offset in A, sigma in A, mean lambda, n) for one sub-window, or nans.

    nan when the window is not covered, has too few pixels, or puts its chi2
    minimum at the edge of the search -- the last is a failure, not a value.
    """
    w = _window(obs, lo, hi)
    if w.mask.sum() < MIN_PIX or lo < float(grid.wave[0]):
        return np.nan, np.nan, np.nan, 0
    c = np.array([_chi2(w, grid, dict(theta, rv=float(v)))[0] for v in rvs])
    if not np.isfinite(c).any():
        return np.nan, np.nan, np.nan, 0
    j = int(np.nanargmin(c))
    if j in (0, len(rvs) - 1):
        return np.nan, np.nan, np.nan, 0
    y0, y1, y2 = c[j - 1], c[j], c[j + 1]
    d = y0 - 2 * y1 + y2
    if not (d > 0):
        return np.nan, np.nan, np.nan, 0
    step = rvs[1] - rvs[0]
    rv = float(rvs[j] + 0.5 * (y0 - y2) / d * step)
    sig_rv = float(np.sqrt(2.0 / d) * step)
    lam = float(np.mean(obs.wavelength[w.mask]))
    n = int(w.mask.sum())
    return lam * rv / C_KMS, lam * sig_rv / C_KMS, lam, n


def fit_star(star, grid, theta, verbose=False):
    """-> {grating: dict(a, b, n, rms) or None} for one star."""
    obs = raw_observation(star)
    out = {}
    for name, _, _, deg, wins in SEGMENTS:
        xs, ys = [], []
        for lo, hi in wins:
            s, _, lam, n = window_offset(obs, grid, theta, lo, hi)
            if np.isfinite(s):
                xs.append(lam - LAM_REF)
                ys.append(s)
            if verbose:
                print(f'      {name:7s} {lo:.0f}-{hi:.0f}  '
                      + (f'{s:+.3f} A  (n={n})' if np.isfinite(s) else '--'))
        if not xs:
            out[name] = None
            continue
        xs, ys = np.array(xs), np.array(ys)
        if deg == 1 and len(xs) >= 3:
            b, a = np.polyfit(xs, ys, 1)
            rms = float(np.std(ys - (a + b * xs)))
        else:
            a, b = float(np.median(ys)), 0.0     # median: robust to one bad window
            rms = float(np.std(ys - a))
        out[name] = dict(a=float(a), b=float(b), n=len(xs), rms=rms)
    return out


def comparison_stars():
    """-> stars outside the fitted sample that still need a correction.

    `explore/plot_ngsl_vs_model.py` draws the NGSL-only comparison from
    data/balmer_candidates.csv, which includes HD040573 -- a star XSL never
    observed, so it has no sample row and no node scan. Leaving it out of the
    table would put it back in exactly the silently-uncorrected state this
    script exists to prevent.

    Read from that file rather than from the previous table: a refit that
    inherits its own star list can only ever lose stars, never recover one.
    """
    path = ROOT / 'data' / 'balmer_candidates.csv'
    if not path.exists():
        return []
    with open(path) as fh:
        return [r['star'] for r in csv.DictReader(fh) if r.get('selected') != 'no']


def catalog_params(star):
    """-> (teff, logg, mh) from data/sample.csv, or data/ngsl_catalog.csv.

    The catalog fallback is for stars outside the fitted sample that are still
    plotted -- HD040573 is in the NGSL-only comparison figures but XSL never
    observed it, so it has no sample row and no node scan. Dropping it from the
    table would leave it silently uncorrected, which is the failure this whole
    exercise exists to remove.
    """
    try:
        r = sample_row(star)
        return float(r['teff_ngsl']), float(r['logg_ngsl']), float(r['mh_ngsl'])
    except KeyError:
        pass
    with open(ROOT / 'data' / 'ngsl_catalog.csv') as fh:
        for r in csv.DictReader(fh):
            if r['target'] == star:
                return float(r['teff']), float(r['logg']), float(r['logz'])
    raise KeyError(f'{star} is in neither sample.csv nor ngsl_catalog.csv')


def theta_for(star, node=True):
    """Node-scan ML parameters when a scan exists, else the catalog values.

    The fallback costs almost nothing: `--node-sensitivity` shows the fitted
    coefficients move by <0.02 A between the two, against a correction of ~0.8 A.
    """
    if node and (ROOT / 'results' / star / 'scan.npz').exists():
        n = best_node(star)
        return dict(teff=n['teff'], logg=n['logg'], mh=n['mh'],
                    ebv=n['ebv'], vsini=n['vsini'])
    t, g, m = catalog_params(star)
    return dict(teff=t, logg=g, mh=m, ebv=0.0, vsini=0.0)


def selftest(grid, star='HD194453'):
    """Inject a known shift into the data and check it comes back, with sign."""
    print(f'Self-test on {star}: inject a known shift, recover it\n')
    theta = theta_for(star)
    obs = raw_observation(star)
    lo, hi = 3700.0, 4100.0
    base, _, lam, _ = window_offset(obs, grid, theta, lo, hi)
    print(f'  window {lo:.0f}-{hi:.0f} A, baseline offset {base:+.3f} A\n')
    print(f'  {"injected":>10s} {"recovered":>10s} {"error":>8s}')
    ok = True
    for inj in (-0.60, -0.30, 0.0, +0.30, +0.60):
        o = copy.copy(obs)
        o.wavelength = obs.wavelength + inj        # move the DATA redward
        o.mask = obs.mask.copy()
        got, _, _, _ = window_offset(o, grid, theta, lo, hi)
        err = got - base - inj
        ok &= abs(err) < 0.06
        print(f'  {inj:>+10.3f} {got - base:>+10.3f} {err:>+8.3f}')
    print('\n  ' + ('[ok] sign and scale recovered to better than 0.06 A'
                    if ok else '[FAIL] injected shift not recovered'))
    return ok


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--star', help='one star (default: the whole sample)')
    ap.add_argument('--selftest', action='store_true')
    ap.add_argument('--node-sensitivity', action='store_true',
                    help='also fit at the catalog node, to show a,b barely move')
    ap.add_argument('--verbose', action='store_true', help='per-window offsets')
    ap.add_argument('--dry-run', action='store_true', help='do not write the CSV')
    a = ap.parse_args()

    grid = Grid()
    if a.selftest:
        raise SystemExit(0 if selftest(grid) else 1)

    if a.star:
        stars = [a.star]
    else:
        stars = [s for s in all_stars() if sample_row(s)['tier'] != 'rejected']
        stars += [s for s in comparison_stars() if s not in set(stars)]
    print(f'Fitted wavelength correction (after air->vacuum), pivot {LAM_REF:.0f} A')
    print(f'{"star":<10}{"grating":<8}{"a (A)":>8}{"b (A/1000A)":>13}'
          f'{"n":>4}{"rms":>7}' + ('  catalog-node a, b' if a.node_sensitivity else ''))
    rows = []
    for s in stars:
        theta = theta_for(s)
        if a.verbose:
            print(f'  {s}  Teff={theta["teff"]:.0f} logg={theta["logg"]:.2f} '
                  f'[M/H]={theta["mh"]:+.2f} E(B-V)={theta["ebv"]:.3f}')
        fit = fit_star(s, grid, theta, verbose=a.verbose)
        alt = (fit_star(s, grid, theta_for(s, node=False))
               if a.node_sensitivity else {})
        for g, c in fit.items():
            if c is None:
                rows.append(dict(star=s, grating=g, a_A='', b_A_per_A='',
                                 n_windows=0, fit_rms_A=''))
                continue
            rows.append(dict(star=s, grating=g, a_A=round(c['a'], 4),
                             b_A_per_A=round(c['b'], 7), n_windows=c['n'],
                             fit_rms_A=round(c['rms'], 3)))
            extra = ''
            if a.node_sensitivity and alt.get(g):
                extra = (f"   {alt[g]['a']:+.3f} {alt[g]['b'] * 1000:+.3f}"
                         f"  (d={c['a'] - alt[g]['a']:+.3f})")
            print(f'{s:<10}{g:<8}{c["a"]:>8.2f}{c["b"] * 1000:>13.3f}'
                  f'{c["n"]:>4}{c["rms"]:>7.2f}{extra}')

    if a.dry_run:
        print(f'\n[dry run] {len(rows)} rows NOT written')
        return
    out = ROOT / 'data' / 'ngsl_wavecal.csv'
    with open(out, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=['star', 'grating', 'a_A', 'b_A_per_A',
                                           'n_windows', 'fit_rms_A'])
        w.writeheader()
        w.writerows(rows)
    print(f'\n{len(rows)} rows ({len(stars)} stars) -> data/ngsl_wavecal.csv')


if __name__ == '__main__':
    main()
