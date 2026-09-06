"""Fit an observed spectrum against the ATLAS12 grid.

Free parameters: Teff, log g, [M/H], E(B-V), v sin i, instrumental broadening,
radial velocity, and an error-scale term. The flux normalization is NOT free in
the sampler -- it is solved analytically at each call (see `_scale`), which is
exact for a single multiplicative constant and removes a dimension.

There is deliberately no continuum polynomial: it would absorb the continuum
shape that constrains Teff and E(B-V). Those two are consequently covariant,
and the fit leans on three things to separate them:

  * the Balmer break amplitude, which responds sharply to Teff and only weakly
    to a smooth reddening law,
  * as wide a wavelength baseline as the data allow,
  * a prior on E(B-V) (see `dust_prior`).

The error scale exists because NGSL's STATERR is optimistic by ~3x (it carries
propagated counting statistics only). Fixing the errors at face value would
produce absurdly tight parameter errors.
"""
import sys
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from fitting.model import Grid, forward

PARAMS = ['teff', 'logg', 'mh', 'ebv', 'vsini', 'inst', 'rv', 'lnerr']


def hydrogen_lines(wmin, wmax, series=(2, 3), nmax=40):
    out = []
    for m in series:
        for n in range(m + 1, nmax):
            lam = 911.7635 / (1.0 / m ** 2 - 1.0 / n ** 2)
            if wmin <= lam <= wmax:
                out.append(lam)
    return np.array(out)


def build_mask(wave, fit_ranges=None, exclude=None, mask_h_cores=0.0,
               h_series=(2, 3)):
    """Boolean mask of pixels to FIT.

    fit_ranges     list of (lo, hi) to include; None = everything
    exclude        list of (lo, hi) to drop (e.g. detector gaps, bad regions)
    mask_h_cores   half-width in A to drop around every hydrogen line; use this
                   to fit the break and wings while ignoring cores that LTE
                   models get wrong (they are filled by ~10% of the line EW,
                   see CAVEATS.md)
    """
    m = np.zeros_like(wave, bool) if fit_ranges else np.ones_like(wave, bool)
    for lo, hi in (fit_ranges or []):
        m |= (wave >= lo) & (wave <= hi)
    for lo, hi in (exclude or []):
        m &= ~((wave >= lo) & (wave <= hi))
    if mask_h_cores > 0:
        for lam in hydrogen_lines(wave.min(), wave.max(), h_series):
            m &= np.abs(wave - lam) > mask_h_cores
    return m


@dataclass
class DustPrior:
    """Prior on E(B-V).

    kind='upper'    U(0, value). The correct use of an SFD98/SF11 map value:
                    it is the TOTAL column through the Galactic dust layer and
                    so an upper bound for a star inside it, not an estimate.
                    For our stars the map exceeds the true reddening by up to
                    8x (HD040573: SF11 0.470 vs 0.06 photometric at 190 pc).
    kind='gaussian' N(value, sigma). Appropriate for a 3D/tomographic map
                    evaluated at the star's distance, or a measured value such
                    as the one UVES-POP publishes.
    kind='fixed'    hold E(B-V) at value.
    """
    kind: str = 'upper'
    value: float = 1.0
    sigma: float = 0.02

    def logp(self, ebv):
        if ebv < 0:
            return -np.inf
        if self.kind == 'upper':
            return 0.0 if ebv <= self.value else -np.inf
        if self.kind == 'gaussian':
            return -0.5 * ((ebv - self.value) / self.sigma) ** 2
        if self.kind == 'fixed':
            return 0.0
        raise ValueError(self.kind)


@dataclass
class FitConfig:
    inst_kind: str = 'R'                 # 'R', 'fwhm' or 'ngsl'
    bounds: dict = field(default_factory=dict)
    dust: DustPrior = field(default_factory=DustPrior)
    fixed: dict = field(default_factory=dict)
    r_v: float = 3.1


def default_bounds(grid):
    b = grid.bounds()
    return dict(teff=b['teff'], logg=b['logg'], mh=b['mh'],
                ebv=(0.0, 1.0), vsini=(0.0, 400.0), inst=(500.0, 200000.0),
                rv=(-400.0, 400.0), lnerr=(-1.0, 2.5))


def _scale(obs, err, mod):
    """Analytic maximum-likelihood multiplicative normalization."""
    iv = 1.0 / err ** 2
    den = np.sum(mod * mod * iv)
    return np.sum(obs * mod * iv) / den if den > 0 else 1.0


class SpectrumFit:
    def __init__(self, grid, wave, flux, err, mask, cfg):
        self.grid, self.cfg = grid, cfg
        self.w, self.f, self.e = wave[mask], flux[mask], err[mask]
        self.bounds = {**default_bounds(grid), **cfg.bounds}
        self.free = [p for p in PARAMS if p not in cfg.fixed]

    def unpack(self, theta):
        d = dict(zip(self.free, theta))
        d.update(self.cfg.fixed)
        return d

    def model(self, p):
        return forward(self.grid, self.w, p['teff'], p['logg'], p['mh'],
                       ebv=p.get('ebv', 0.0), vsini=p.get('vsini', 0.0),
                       rv=p.get('rv', 0.0), inst_kind=self.cfg.inst_kind,
                       inst_value=p.get('inst', 0.0), r_v=self.cfg.r_v)

    def log_prior(self, p):
        for k in self.free:
            lo, hi = self.bounds[k]
            if not (lo <= p[k] <= hi):
                return -np.inf
        return self.cfg.dust.logp(p.get('ebv', 0.0))

    def log_prob(self, theta):
        p = self.unpack(theta)
        lp = self.log_prior(p)
        if not np.isfinite(lp):
            return -np.inf
        try:
            mod = self.model(p)
        except ValueError:
            return -np.inf                 # unfilled grid cell
        if not np.all(np.isfinite(mod)) or np.all(mod <= 0):
            return -np.inf
        mod = mod * _scale(self.f, self.e, mod)
        s = np.exp(p.get('lnerr', 0.0))
        var = (self.e * s) ** 2
        return lp - 0.5 * np.sum((self.f - mod) ** 2 / var + np.log(2 * np.pi * var))

    def optimize(self, start):
        x0 = [start[k] for k in self.free]
        r = minimize(lambda t: -self.log_prob(t), x0, method='Nelder-Mead',
                     options=dict(maxiter=4000, xatol=1e-3, fatol=1e-3))
        return self.unpack(r.x), -r.fun

    def sample(self, start, nwalkers=32, nsteps=2000, burn=500, progress=False):
        import emcee
        x0 = np.array([start[k] for k in self.free])
        scatter = np.array([max(1e-3, abs(v) * 1e-3) for v in x0])
        for i, k in enumerate(self.free):
            lo, hi = self.bounds[k]
            scatter[i] = 0.01 * (hi - lo)
        pos = x0 + scatter * np.random.randn(nwalkers, len(x0))
        for i, k in enumerate(self.free):
            pos[:, i] = np.clip(pos[:, i], *self.bounds[k])
        s = emcee.EnsembleSampler(nwalkers, len(x0), self.log_prob)
        s.run_mcmc(pos, nsteps, progress=progress)
        return s.get_chain(discard=burn, flat=True), self.free
