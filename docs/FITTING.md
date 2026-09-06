# Fitting

Measure stellar parameters from the spectra rather than adopting them from a
catalog. Catalog values may be wrong, biased, or derived by a method carrying
its own systematics relative to these data — UVES-POP's, for instance, come
from fitting a PHOENIX grid, so adopting them wholesale would import a
code-to-code systematic into an ATLAS12 comparison.

Code: `fitting/model.py` (forward model), `fitting/fit.py` (likelihood,
priors, samplers). Grid: see [GRID.md](GRID.md).

## Free parameters

| parameter | role |
|---|---|
| Teff, log g, [M/H] | interpolated in the grid |
| E(B−V) | CCM89, R_V = 3.1 by default |
| v sin i | Gray (2005) rotation profile, linear limb darkening |
| instrumental | constant R, constant FWHM in Å, or the measured NGSL profile |
| RV | Doppler shift |
| error scale | multiplies the quoted uncertainties |

The **normalization is not sampled**. For a single multiplicative constant the
maximum-likelihood value is available in closed form, so it is solved at every
likelihood call and one dimension is removed at no cost.

The **error scale** is there because NGSL's `STATERR` is optimistic by ~3× — it
carries propagated counting statistics only. Fitting with the quoted errors at
face value yields parameter uncertainties that are confidently wrong.

## No continuum polynomial, and what that costs

A free multiplicative polynomial would absorb the continuum shape — which is
exactly the information that constrains Teff and E(B−V) — leaving only line
profiles to carry the temperature. It is therefore omitted, and the fit relies
on the spectrophotometry being good: NGSL is space-based and calibrated to ~3%,
UVES-POP quotes 1.5–4% absolute.

The price is that **Teff and E(B−V) are covariant**: both tilt the continuum.
Three things separate them.

1. **The Balmer break amplitude** responds sharply to Teff and only weakly to a
   smooth reddening law. This is the main lever, and the reason the break is
   worth fitting rather than masking.
2. **A wide wavelength baseline.** Over 3200–9500 Å the CCM89 curvature and a
   temperature change are distinguishable in a way they are not over 400 Å.
3. **A prior on E(B−V)** (below).

Where the data cannot break the degeneracy, sample rather than optimize: the
posterior will show the Teff–E(B−V) banana honestly, and a point estimate will
not.

## Constrain the instrumental broadening — do not fit it blind

For NGSL the profile is **measured, not assumed**: matching XSL to NGSL for
three stars in common gives **R = 600 ± 40**, constant in velocity, with no
model involved (see [DATA.md](DATA.md)). Use it.

```python
from fitting.fit import ngsl_config
cfg = ngsl_config()                      # inst_kind='R', prior N(600, 40)
cfg = ngsl_config(fixed=dict(inst=600))  # or hold it exactly
```

Leaving `inst` free re-opens its degeneracy with Teff and log g — line depth
trades against broadening — for no gain, because three stars constrain it
better than one spectrum can. This is not hypothetical: while the fitter
reimplemented its own kernels, `kind='R'` ignored the grid spacing and turned a
request for R=600 into R=83. Because `inst` was free, nothing crashed; the fit
would have absorbed the error into an absurd broadening value and dragged Teff
and log g with it. The kernels now delegate to `common.lsf`, and a prior means
such a failure shows up as a fight with the prior rather than passing silently.

Per library:

| library | `inst_kind` | value |
|---|---|---|
| NGSL | `'R'` | 600 ± 40 (or `'ngsl'`, which applies it directly) |
| XSL | `'R'` | ~9800 UVB, ~11600 VIS — constant in velocity |
| UVES-POP | — | R = 80,000 *then* a 0.1 Å rebin; neither pure kind is exact, apply the boxcar as `plot_uves_vs_model.py` does |

`FitConfig.priors` takes any `{parameter: object with .logp}`, so the same
mechanism constrains v sin i from a published value, or RV, or anything else
known independently.

## The dust prior

`DustPrior` has three modes, and the choice matters more than it looks.

| kind | form | when |
|---|---|---|
| `upper` | U(0, value) | an SFD98/SF11 **map column** |
| `gaussian` | N(value, σ) | a 3D/tomographic map at the star's distance, or a published fitted value |
| `fixed` | — | hold it |

**An SFD/SF11 value is an upper bound, not an estimate.** It is the total
column through the entire Galactic dust layer, while these stars sit inside it
at 137–317 pc. The overestimate is large and not uniform: for HD040573 SF11
gives 0.470 against 0.06 photometric, 8×. Using it as a Gaussian centre would
force the fit to absurd reddening and drag Teff with it along the degeneracy.

For a real estimate, use a **tomographic map evaluated at the star's distance**
(parallaxes are in `data/reddening.csv`). Bayestar19 is queryable through the
Argonaut API without downloading the map, but is PS1-based and covers Dec > −30
— three of the four UVES-POP stars sit near Dec −34 and fall outside it.
Lallement/Vergely and Edenhofer cover the southern sky at coarser resolution.
UVES-POP publishes its own fitted E(B−V) per star, which is a measurement along
the actual sightline and preferable to any map.

## Masking

`build_mask` takes arbitrary `fit_ranges` and `exclude` windows, plus
`mask_h_cores`: a half-width in Å dropped around every Balmer and Paschen line.

That last one exists because the observed hydrogen cores carry a genuine flux
excess of ~10% of the line equivalent width relative to these LTE models —
almost certainly NLTE in hydrogen, which the code does not treat for H. Masking
the cores lets the break and the line wings constrain Teff and log g without
the fit trying to accommodate physics the models are missing. Set it to 0 to
fit the cores deliberately, e.g. to measure that excess.

Detector gaps must also be excluded: UVES-POP has a dichroic gap at 5750–5844 Å,
inter-order gaps redward of 8515 Å, and one star (HD162678) is missing
3859–4779 Å entirely.

## Usage

```python
from fitting.model import Grid
from fitting.fit import SpectrumFit, FitConfig, DustPrior, build_mask

grid = Grid()                                    # models/grid.npz
mask = build_mask(wave, fit_ranges=[(3300, 4600)], mask_h_cores=6.0)
cfg  = FitConfig(inst_kind='ngsl',
                 dust=DustPrior('upper', 0.13),  # SF11 column for this sightline
                 fixed=dict(rv=0.0))
f    = SpectrumFit(grid, wave, flux, err, mask, cfg)

best, lnp = f.optimize(start)                    # Nelder-Mead point estimate
chain, names = f.sample(best, nwalkers=32, nsteps=2000)   # emcee posterior
```

Optimize first, then start the sampler from the optimum. Sampling is what
exposes the Teff–E(B−V) covariance; the point estimate hides it.
