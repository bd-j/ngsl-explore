"""Brute-force scan of every grid node, storing whole conditional curves.

For each of the 1705 nodes (Teff x log g x [M/H], no interpolation) this fits
only the two nuisance directions and keeps the FULL curve in each:

    chi2_bands[nt,ng,nm,nE]   NGSL bands vs E(B-V)
    chi2_xsl  [nt,ng,nm,nV]   XSL lines  vs v sin i
    resid_balmer / resid_paschen [nt,ng,nm,nE]   held out, never fitted

Keeping curves rather than argmaxes is the point: the break prediction can then
be MARGINALISED over the nuisance parameters instead of evaluated at a plug-in
value, which is the difference between an honest envelope and a misleadingly
tight curve.

WHY chi^2 AND NOT lnL. PLAN.md specified lnl_* arrays. Storing chi^2 with its
pixel count instead is strictly more informative and costs nothing: the error
model is still an open decision (independent-pixel chi^2 overstates confidence
on correlated residuals), and baking a likelihood convention into the stored
file would mean re-running the whole scan to change it. fitting/likelihood.py
converts chi^2 + N to lnL downstream.

WHAT IS NOT SCANNED, AND WHY -- each measured, not assumed:

  v sin i in the NGSL legs.  The held-out residual medians are identical to
      3 decimal places over v sin i = 0 -> 300 km/s: NGSL's R = 600 profile is
      6-13 A FWHM and a 300 km/s rotation at 3800 A is 2.5 A, so the instrument
      swamps it. The bands move by <0.006 mmag. So the E(B-V) loop runs at
      v sin i = 0 and the two legs factorise exactly.
  E(B-V) in the XSL leg.  Delta chi^2 = +2.0 across the entire 0-0.12 range,
      against +17000 for v sin i -- and what little there is is a constraint on
      XSL's CONTINUUM, which slit losses make untrustworthy and which the
      segmented polynomial exists to marginalise away.

The forward model is NOT approximated. It would have been ~2x faster to drop
the instrument profile from the band projection (band integrals are conserved
under convolution to 9e-4, a tenth of the calibration floor), but two code
paths disagreeing about the NGSL LSF is exactly how this project once had the
fitter and the figures at R=939 and R=83. Every number here comes from the same
predict/solve path the figures use; `verify()` checks that on real nodes.

The speedups are exact ones: the node spectrum and the reddening curves are
hoisted out of the inner loops, and ONE broadened spectrum serves the bands and
both held-out windows instead of the three that separate predict() calls would
compute.

    python3 fitting/scan.py --star HD194453 --verify
    python3 fitting/scan.py --all
"""
import argparse
import csv
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from fitting.model import Grid
from fitting.observations import (conditioning_set, heldout, BREAK_WINDOW,
                                  PASCHEN_WINDOW)
from fitting.predict import predict, instrument, project, spectrum_at
from fitting.calibration import solve, usable, chi2
from fitting.likelihood import curvature
from common.extinction_ccm import redden
from common.lsf import broaden_rot

C_KMS = 2.99792458e5

ROOT = Path(__file__).resolve().parent.parent
EBV_MAX, N_EBV = 0.30, 61
# Capped at 200 km/s. The 250 and 300 nodes were doing no useful work: the only
# solutions that reached them were the below-grid-[M/H] stars, where broadening
# was being spent to wash out model metal lines the grid cannot make weak enough
# (PLAN.md), and a ceiling that high just let that go further before showing up.
VSINI_MAX = 200.0
VSINI_GRID = np.array([0., 10., 20., 30., 40., 60., 80., 110., 150., 200.])


def attenuation(wave, ebvs, r_v=3.1):
    """(nE, npix) multiplicative extinction, built through `redden` itself.

    Taking the ratio of a reddened unit spectrum to the unit spectrum is not a
    reimplementation of CCM89 -- it IS CCM89, called once per E(B-V) instead of
    once per node per E(B-V). Any future change to the extinction law is picked
    up here automatically.
    """
    ones = np.ones_like(wave)
    return np.array([redden(wave, ones, float(e), r_v) if e else ones
                     for e in ebvs])


def obs_wave(obs, w0, rv=0.0):
    """The wavelength array predict() would use for this observation.

    NOT always w0. XSL carries a measured velocity zero point in `rv_fixed`
    (+4.21 km/s plus a per-star departure), and predict() shifts the grid before
    both the instrument convolution and the projection. Reusing w0 here made the
    XSL chi^2 differ from predict() by 6.6e-3 -- small enough to look like a
    tolerance and large enough to rank nodes wrongly. verify() caught it.
    """
    v = obs.rv_fixed if obs.rv_fixed is not None else rv
    return w0 * (1.0 + v / C_KMS) if v else w0


def node_curves(f0, w0, atten, nb, hb, hp, xs, vsinis, grid):
    """All four conditional curves for ONE node.

    f0 is the node spectrum on the grid's own wavelengths, unbroadened.
    """
    w_nb, w_hb, w_hp, w_xs = (obs_wave(o, w0) for o in (nb, hb, hp, xs))
    nE = atten.shape[0]
    c_b = np.full(nE, np.nan)
    r_b = np.full(nE, np.nan)
    r_p = np.full(nE, np.nan)
    scal = np.full(nE, np.nan)
    # ln|A| is a property of the MODEL at this node, so it cannot be recovered
    # from a stored chi^2 afterwards. Storing it here is what lets
    # fitting/likelihood.py marginalise the calibration coefficients without
    # re-running the scan.
    d_b = np.full(nE, np.nan)
    d_x = np.full(len(vsinis), np.nan)

    for e in range(nE):
        f = f0 * atten[e]
        # one instrument call, three consumers -- bands, Balmer, Paschen. They
        # share it only because all three carry rv_fixed=None and the same
        # ('R', 600) profile; asserted below so a change to either is caught.
        fi = instrument(w_nb, f, nb.resolution)
        pb = project(nb, w_nb, fi)
        try:
            cal, c = solve(nb, pb)
        except (ValueError, np.linalg.LinAlgError):
            continue
        if usable(nb, pb).sum() < 3:
            continue
        c_b[e] = chi2(nb, cal)
        d_b[e] = curvature(nb, pb)[0]
        scal[e] = c[0]
        for arr, h, wh in ((r_b, hb, w_hb), (r_p, hp, w_hp)):
            m = np.interp(h.wavelength, wh, fi) * c[0]
            ok = h.mask & np.isfinite(m) & (m > 0)
            if ok.any():
                arr[e] = np.nanmedian((h.flux[ok] - m[ok]) / m[ok])

    c_x = np.full(len(vsinis), np.nan)
    for v, vs in enumerate(vsinis):
        fr = broaden_rot(w0, f0, float(vs)) if vs else f0
        px = project(xs, w_xs, instrument(w_xs, fr, xs.resolution))
        try:
            calx, _ = solve(xs, px)
        except (ValueError, np.linalg.LinAlgError):
            continue
        if usable(xs, px).sum() < 3:
            continue
        c_x[v] = chi2(xs, calx)
        d_x[v] = curvature(xs, px)[0]
    return c_b, c_x, r_b, r_p, scal, d_b, d_x


def verify(star, grid, obs, atten, ebvs, vsinis, n=4, seed=0):
    """Check the scan's fast path against predict() on real nodes.

    The scan hoists work out of loops and shares one broadened spectrum between
    three projections. Those are meant to be algebraically identical to calling
    predict() per observation -- 'meant to be' is why this exists.
    """
    nb, hb, hp, xs = obs
    rng = np.random.default_rng(seed)
    w0 = grid.wave
    worst = dict(bands=0.0, xsl=0.0, balmer=0.0, paschen=0.0,
                 lndetA_b=0.0, lndetA_x=0.0)
    for _ in range(n):
        i, j, k = (rng.integers(len(grid.teff)), rng.integers(len(grid.logg)),
                   rng.integers(len(grid.mh)))
        t, lg, mh = float(grid.teff[i]), float(grid.logg[j]), float(grid.mh[k])
        c_b, c_x, r_b, r_p, scal, d_b, d_x = node_curves(
            spectrum_at(grid, t, lg, mh), w0, atten, nb, hb, hp, xs, vsinis,
            grid)
        e = int(rng.integers(len(ebvs)))
        v = int(rng.integers(len(vsinis)))
        th = dict(teff=t, logg=lg, mh=mh, ebv=float(ebvs[e]), vsini=0.0)
        ref_b = chi2(nb, solve(nb, predict(th, [nb], grid)[0].value)[0])
        calr, cr = solve(nb, predict(th, [nb], grid)[0].value)
        thx = dict(teff=t, logg=lg, mh=mh, ebv=0.0, vsini=float(vsinis[v]))
        ref_x = chi2(xs, solve(xs, predict(thx, [xs], grid)[0].value)[0])
        refs = {}
        for nm, h in (('balmer', hb), ('paschen', hp)):
            m = predict(th, [h], grid)[0].value * cr[0]
            ok = h.mask & np.isfinite(m) & (m > 0)
            refs[nm] = np.nanmedian((h.flux[ok] - m[ok]) / m[ok])
        ref_db = curvature(nb, predict(th, [nb], grid)[0].value)[0]
        ref_dx = curvature(xs, predict(thx, [xs], grid)[0].value)[0]
        for key, got, ref in (('bands', c_b[e], ref_b), ('xsl', c_x[v], ref_x),
                              ('balmer', r_b[e], refs['balmer']),
                              ('paschen', r_p[e], refs['paschen']),
                              ('lndetA_b', d_b[e], ref_db),
                              ('lndetA_x', d_x[v], ref_dx)):
            d = abs(got - ref) / max(abs(ref), 1e-30)
            worst[key] = max(worst[key], d)
        print(f'    node T={t:.0f} g={lg:.1f} z={mh:+.1f}  '
              f'E={ebvs[e]:.3f} v={vsinis[v]:.0f}  '
              + '  '.join(f'{k2} {worst[k2]:.1e}' for k2 in worst))
    bad = {k: v for k, v in worst.items() if v > 1e-9}
    if bad:
        raise SystemExit(f'  VERIFY FAILED, fast path differs from predict(): {bad}')
    print(f'  verify OK: fast path matches predict() to '
          f'{max(worst.values()):.1e} relative on {n} nodes\n')


def scan_star(star, grid, ebvs, vsinis, verify_first=False, progress=True):
    cond = {o.name: o for o in conditioning_set(star)}
    nb, xs = cond['ngsl_bands'], cond['xsl']
    hb = heldout(star, window=BREAK_WINDOW)
    hp = heldout(star, window=PASCHEN_WINDOW)
    w0 = grid.wave
    atten = attenuation(w0, ebvs)
    # node_curves reuses one broadened spectrum for the bands and both held-out
    # windows. That is only valid while the three agree on velocity and profile.
    if not (nb.rv_fixed == hb.rv_fixed == hp.rv_fixed
            and nb.resolution == hb.resolution == hp.resolution):
        raise ValueError('bands and held-out windows no longer share a '
                         'velocity/resolution; node_curves must stop sharing '
                         'one instrument call')

    nt, ng, nm = len(grid.teff), len(grid.logg), len(grid.mh)
    nE, nV = len(ebvs), len(vsinis)
    shape_e, shape_v = (nt, ng, nm, nE), (nt, ng, nm, nV)
    out = dict(chi2_bands=np.full(shape_e, np.nan),
               chi2_xsl=np.full(shape_v, np.nan),
               resid_balmer=np.full(shape_e, np.nan),
               resid_paschen=np.full(shape_e, np.nan),
               scalar=np.full(shape_e, np.nan),
               lndetA_bands=np.full(shape_e, np.nan),
               lndetA_xsl=np.full(shape_v, np.nan))

    probe = spectrum_at(grid, float(grid.teff[len(grid.teff) // 2]),
                        float(grid.logg[len(grid.logg) // 2]),
                        float(grid.mh[len(grid.mh) // 2]))
    k_b = curvature(nb, project(nb, obs_wave(nb, w0),
                                instrument(obs_wave(nb, w0), probe,
                                           nb.resolution)))[1]
    k_x = curvature(xs, project(xs, obs_wave(xs, w0),
                                instrument(obs_wave(xs, w0), probe,
                                           xs.resolution)))[1]
    print(f'{star}: {int(grid.filled.sum())} nodes x ({nE} E(B-V) + {nV} v sin i)'
          f'   n_bands={nb.ndata} n_xsl={xs.ndata} '
          f'n_balmer={hb.ndata} n_paschen={hp.ndata}  '
          f'k_bands={k_b} k_xsl={k_x}')
    if verify_first:
        verify(star, grid, (nb, hb, hp, xs), atten, ebvs, vsinis)

    t0 = time.perf_counter()
    done, shown = 0, 0
    total = int(grid.filled.sum())
    for i in range(nt):
        for j in range(ng):
            for k in range(nm):
                if not grid.filled[i, j, k]:
                    continue
                c_b, c_x, r_b, r_p, sc, d_b, d_x = node_curves(
                    spectrum_at(grid, float(grid.teff[i]), float(grid.logg[j]),
                                float(grid.mh[k])),
                    w0, atten, nb, hb, hp, xs, vsinis, grid)
                out['chi2_bands'][i, j, k] = c_b
                out['chi2_xsl'][i, j, k] = c_x
                out['resid_balmer'][i, j, k] = r_b
                out['resid_paschen'][i, j, k] = r_p
                out['scalar'][i, j, k] = sc
                out['lndetA_bands'][i, j, k] = d_b
                out['lndetA_xsl'][i, j, k] = d_x
                done += 1
        if progress and done > shown:
            shown, el = done, time.perf_counter() - t0
            print(f'  {done:>5}/{total} nodes  {el / 60:5.1f} min'
                  f'  eta {el / done * (total - done) / 60:5.1f} min', flush=True)

    d = ROOT / 'results' / star
    d.mkdir(parents=True, exist_ok=True)
    p = d / 'scan.npz'
    np.savez_compressed(
        p, teff=grid.teff, logg=grid.logg, mh=grid.mh, ebv=ebvs, vsini=vsinis,
        n_bands=nb.ndata, n_xsl=xs.ndata, n_balmer=hb.ndata,
        n_paschen=hp.ndata, k_bands=k_b, k_xsl=k_x, star=star, **out)
    print(f'  -> {p.relative_to(ROOT)}  ({p.stat().st_size / 1e6:.1f} MB, '
          f'{(time.perf_counter() - t0) / 60:.1f} min)')
    return p


def posterior(star, scale='profile', weight='inverse_dof', vsini='profile'):
    """Read results/<star>/scan.npz -> combined objective over (node, E(B-V)).

    Lives here, not in the plotting script, because two consumers now ask which
    node is best -- explore/plot_scan.py and explore/check_predict.py -- and a
    figure drawn at a different "best node" than the one the analysis reports
    is exactly the drift this project keeps having to undo.

    The two legs get DIFFERENT coefficient priors, for the reason set out in
    fitting/likelihood.py: the bands carry a single scalar, which is (R/d)^2
    times a grey factor, and a flat prior on it would impose a spurious
    ~ -4 ln(T) tilt across the grid. With k = 1 the scale-invariant marginal is
    closed form, so the bands use it. XSL's 12 segmented coefficients have no
    single amplitude to factor out, so that leg uses the flat prior.

    LEG WEIGHTING. weight='inverse_dof' divides each leg by its own pixel count,
    i.e. combines chi^2/dof rather than chi^2. Unweighted, XSL brings 2342
    pixels against the bands' 13, so node ranking was determined almost entirely
    by XSL and the bands -- the ENTIRE dust lever, and the only thing that sees
    the continuum -- barely moved it. The symptom was concrete: on HD194453 the
    unweighted maximum sat at band chi^2/N = 4.4 with a +3.4% kink in the two
    bluest bands, while a solution fitting the bands at chi^2/N = 0.52 existed
    and predicted both held-out breaks better (+0.39%/+0.15% against
    +1.11%/+2.25%). The scan was choosing the worse model on every measure that
    matters here.

    Inverse-dof is a blunt instrument, and it is a stand-in for the real
    quantity, which is the EFFECTIVE number of independent points. XSL's pixels
    are correlated over the LSF (~4 px) and by its continuum, so its effective
    dof is far below 2342; 13 banded NGSL points with a 1% calibration floor are
    much closer to independent. Weighting by 1/n asserts the two legs deserve
    equal total say, which is defensible but not derived. weight=None restores
    the raw sum.

    v sin i is PROFILED rather than marginalised whenever a weight is applied:
    once a leg's ln L is scaled by 1/n it is no longer a likelihood, and a
    logsumexp over it would be integrating something that is not a density.
    """
    from fitting.likelihood import lnlike
    d = np.load(ROOT / 'results' / star / 'scan.npz')
    lb = lnlike(d['chi2_bands'], int(d['n_bands']), d['lndetA_bands'],
                int(d['k_bands']), scale=scale, coeff_prior='scale_free',
                c0=d['scalar'])
    lx = lnlike(d['chi2_xsl'], int(d['n_xsl']), d['lndetA_xsl'],
                int(d['k_xsl']), scale=scale, coeff_prior='flat')
    # Existing scan.npz files were written with nodes up to 300 km/s. Apply the
    # cap at read time so the change takes effect without re-running the scan;
    # a scan written after this will simply have no columns to drop.
    keep = np.asarray(d['vsini'], float) <= VSINI_MAX + 1e-9
    if not keep.all():
        lx = lx[..., keep]
    if vsini == 'marginal' and weight is None:
        # a nuisance parameter, integrated out so a node that needs a finely
        # tuned rotation is not rewarded for the tuning
        m = np.nanmax(lx, axis=-1)
        lx_v = m + np.log(np.nansum(np.exp(lx - m[..., None]), axis=-1))
    else:
        lx_v = np.nanmax(lx, axis=-1)
    if weight == 'inverse_dof':
        wb, wx = 1.0 / int(d['n_bands']), 1.0 / int(d['n_xsl'])
    elif weight is None:
        wb = wx = 1.0
    else:
        raise ValueError(f'unknown weight: {weight!r}')
    return d, wb * lb + wx * lx_v[..., None]


def best_node(star, scale='profile', weight='inverse_dof'):
    """-> dict(teff, logg, mh, ebv, vsini, ...) at the maximum of the posterior.

    v sin i is reported CONDITIONAL on the winning node, which is what a figure
    drawn at those parameters needs -- the marginal used for ranking has no
    single v sin i attached to it.
    """
    d, lnl = posterior(star, scale, weight)
    i, j, k, e = np.unravel_index(np.nanargmax(lnl), lnl.shape)
    cx = d['chi2_xsl'][i, j, k]
    vs = np.asarray(d['vsini'], float)
    cx = np.where(vs <= VSINI_MAX + 1e-9, cx, np.nan)
    v = int(np.nanargmin(cx))
    return dict(star=star, teff=float(d['teff'][i]), logg=float(d['logg'][j]),
                mh=float(d['mh'][k]), ebv=float(d['ebv'][e]),
                vsini=float(d['vsini'][v]),
                chi2n_bands=float(d['chi2_bands'][i, j, k, e]) / int(d['n_bands']),
                chi2n_xsl=float(d['chi2_xsl'][i, j, k, v]) / int(d['n_xsl']),
                at_mh_floor=bool(k == 0), idx=(int(i), int(j), int(k), int(e)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--star', default='HD194453')
    ap.add_argument('--all', action='store_true')
    ap.add_argument('--verify', action='store_true',
                    help='check the fast path against predict() first')
    ap.add_argument('--ebv-max', type=float, default=EBV_MAX)
    ap.add_argument('--nebv', type=int, default=N_EBV)
    a = ap.parse_args()

    grid = Grid()
    ebvs = np.linspace(0.0, a.ebv_max, a.nebv)
    stars = ([r['star'] for r in csv.DictReader(open(ROOT / 'data' / 'sample.csv'))
              if r['tier'] in ('primary', 'secondary')] if a.all else [a.star])
    for s in stars:
        try:
            scan_star(s, grid, ebvs, VSINI_GRID, verify_first=a.verify)
        except Exception as exc:
            print(f'{s}: FAILED {type(exc).__name__}: {exc}')


if __name__ == '__main__':
    main()
