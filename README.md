# Balmer break: empirical spectra vs model atmospheres

How well do model stellar atmospheres reproduce the **Balmer break** — the
continuum discontinuity at 3646 A where bound-free absorption from hydrogen
n = 2 cuts off — in real stars near 10,000 K, where the break is strongest?

The approach: take well-calibrated observed spectra of A stars, compute ATLAS12
model atmospheres and SYNTHE spectra, and compare the continuum shape across
the break. Parameters are then **fitted** rather than adopted from catalogs,
since catalog values may carry systematics of their own. The Paschen break
(8206 A) comes along for free as a consistency check.

## Original brief

> Find the Next Generation Spectral Library database of HST UV spectra of stars
> and any associated documentation. Summarize the wavelength coverage and
> resolution by each grating, the typical signal to noise, and stellar parameter
> coverage of the library. Summarize the way the spectra are flux calibrated.
> Try to find a star at around 10000 degrees Kelvin temperature that would be
> good for testing stellar spectral models of the Balmer break shape (at a
> wavelength of ~3600 AA).
>
> Data and documentation: https://archive.stsci.edu/prepds/stisngsl/

## Documentation

| | |
|---|---|
| [docs/DATA.md](docs/DATA.md) | the spectral libraries and the selected sample |
| [docs/GRID.md](docs/GRID.md) | the 1705-node ATLAS12 model grid |
| [docs/FITTING.md](docs/FITTING.md) | what conditions on what, held-out design, error budget |
| [docs/LSF.md](docs/LSF.md) | the NGSL line spread function: how it was measured and what it is |
| [docs/CAVEATS.md](docs/CAVEATS.md) | **known issues and traps — read before trusting any number** |
| [docs/STALE.md](docs/STALE.md) | retracted numbers, superseded designs, dead ends — nothing here is current |
| [PLAN.md](PLAN.md) | current state and what comes next |

Most of the work in this project turned out to be identifying ways the
comparison goes silently wrong: wavelength conventions, resolution mismatches,
reddening, peculiar stars, binaries, detector gaps. CAVEATS.md catalogues them
with symptoms and fixes.

The four reference documents — DATA, GRID, FITTING, LSF — describe the project
as it is now. Where a number or a design was retracted along the way, the record
is in STALE.md rather than in their running text.

## Layout

```
common/    shared: extinction (CCM89), break metric, LSF kernels, line lists
           and species identification, sedpy photometry, panel plotting, IO
             lsf.py           broaden_ngsl / to_ngsl_pixels -- ONE NGSL profile
grid/      model grid construction (make_model.py, build_grid.py, pack_grid.py)
fitting/   observations, forward model, calibration
             observations.py  one record per dataset + conditioning_set/heldout
             predict.py       predict(theta, observations) -> predictions
             calibration.py   the linear nuisance solve
explore/   survey, selection and figure scripts
             ngsl_lsf.py      the whole NGSL LSF measurement
             superseded/      replaced scripts, kept for provenance
data/      catalogs, derived CSVs, selected-star spectra
docs/      this documentation, plus NGSL delivery docs in ngsl_delivery/
figures/   comparison and diagnostic figures
             fits_mh_in_grid/     per-star fit figures, grid can reach the star
             fits_mh_below_grid/  the same, for stars below the [M/H] floor
                                  -- segregated, never pooled (common/figpath.py)
             explore_libraries/   library survey figures: NGSL / UVES-POP /
                                  XSL vs model, the Pickles atlas, the S/N and
                                  parameter-coverage surveys
             ngsl_lsf/            the LSF measurement
models/    ATLAS12 output and the packed grid (gitignored, ~25 GB)
```

## Pipeline

```bash
export ATLAS12=/path/to/atlas12

./explore/fetch_ngsl.sh                  # NGSL, docs, STIS LSFs, Pickles atlas
python3 explore/extract_docs.py          # PDF text (parameters live in the readme)
python3 explore/build_catalog.py         # merge headers + params + magnitudes
python3 explore/measure_snr.py           # empirical S/N, two estimators
python3 explore/crossmatch_libraries.py  # vs MILES and MaStar
python3 explore/candidate_table.py       # candidate selection + reddening cut
python3 explore/reddening.py             # E(B-V): map, photometric, fitted
python3 explore/uves_pop_astars.py       # select A stars from UVES-POP
python3 common/uves_pop_load.py          # UVES-POP parameters + reddening cut

python3 explore/build_sample.py          # the fitted sample: NGSL n XSL, 9000-11000 K
python3 explore/fetch_gaia.py            # Gaia DR3 photometry + XP spectra
python3 explore/build_sample.py          # rerun: folds Gaia RUWE into the binarity flag

python3 grid/build_grid.py --workers 9   # the model grid (~25 h, resumable)
python3 grid/pack_grid.py                # collapse it into models/grid.npz

python3 explore/plot_ngsl_vs_model.py    # NGSL comparison figures
python3 explore/plot_uves_vs_model.py    # UVES-POP comparison figures (superseded)
```

The instrument profile — measured against XSL, no model involved, see
[LSF.md](docs/LSF.md):

```bash
python3 explore/lsf_resolution.py        # the tabulated STIS LSF -> a table
python3 explore/ngsl_lsf.py --selftest   # the machinery, against known answers
python3 explore/ngsl_lsf.py              # the measurement (~15 min) + figures
python3 explore/ngsl_lsf.py --subwindows # constant-A vs constant-R, and the shift
```

Fitting — see [docs/FITTING.md](docs/FITTING.md) for what conditions on what:

```bash
python3 explore/metal_sensitivity.py --union   # XSL metal windows, by measured
                                               # [M/H] sensitivity + species
python3 explore/xsl_line_offsets.py      # XSL velocity offset per line per star
python3 explore/xp_vs_ngsl.py            # XP vs NGSL band ratios (why XP is not used)

python3 explore/check_predict.py  --star HD194453   # end-to-end smoke test
python3 explore/plot_metal_lines.py --star HD194453 # per-feature model vs data
python3 explore/plot_ebv_teff.py  --star HD194453   # chi2 surface, Teff vs E(B-V)
```

Single model for one star:

```bash
python3 grid/make_model.py --star HD194453 --teff 10241 --logg 3.9 --feh 0.0
```

`zscale = 10**[Fe/H]`, synthesis runs at R = 300,000 (SYNTHE's enforced floor)
and is convolved down afterwards, and the `.spec` flux column is Eddington
H_nu, not f_lambda.

## Findings so far

- The **continuum across the break is reproduced well**: away from hydrogen
  lines the residual is ~1%, and D_Balmer agrees to 0.02-0.08 mag. With the
  break **held out and predicted** rather than fitted (see
  [FITTING.md](docs/FITTING.md)), HD194453 comes out at **+1.2% (Balmer)** and
  **+1.1% (Paschen)** at the band-preferred reddening — with v sin i pinned and
  only one node tried, so a hint rather than a result.
- The experiment is **Teff-limited, not dust-limited**. dD/dTeff = −0.017 mag
  per 100 K against dD/dE(B−V) = +0.0027 per 0.01 mag at fixed Teff, so the
  binding requirement is σ(Teff) ≈ 100 K — about 3× better than published values
  for these stars.
- The **Balmer line-core residual is consistent with zero**, and is degenerate
  with the instrument profile. Core-minus-continuum against the models has a
  median of +0.42% over the eight in-grid stars (−1.21% to +3.41%, both signs)
  under the measured LSF — but −4.08% under a 7.00 Å core, so it can only be
  quoted alongside the profile in use. None of it is a kernel error: the same
  profile reproduces NGSL's cores from smoothed XSL to +0.01% median over 117
  line×star combinations, with no model involved. See [LSF.md](docs/LSF.md).
- **NGSL is in air**, not vacuum, with a linear-in-lambda residual per grating
  that is recalibrated against the models.
- **NGSL's `STATERR` is optimistic by ~3x**, and its delivered line spread
  function is a **Moffat**: core 3.54 A (G430L) / 8.38 A (G750L), beta = 1.52,
  constant in Angstroms per grating, with ~2.4% of its power beyond ±10 A where
  a Gaussian of the same core puts 0.01%. That halo is **real instrumental
  scattered light**: the STIS tables give all four slit widths the same core and
  different wings, and the tabulated **52x0.5** profile — zero free parameters —
  reproduces the fitted halo (2.56% against 2.44%) and nulls the Balmer core
  residual. NGSL used 52x0.2, whose tabulated profile has no halo, so the
  delivered spectra behave like a slit 2.5× wider than the one used. A single Gaussian
  forced on the same data lands at R = 600, which is neither a width nor a
  resolution: it drifts with wavelength (FWHM ∝ λ^0.67 against λ^0.14 for the
  Moffat core) and that drift is what once made the profile look constant in
  velocity. See [LSF.md](docs/LSF.md); earlier claims of R = 939 and R = 665 are
  also retracted there.
- **Model spectra must be INTEGRATED onto detector pixels**, not sampled at
  pixel centres, and a kernel width is meaningless without saying which was
  used: the same data give a 3.54 A Moffat core under integration and 4.06 A
  under centre sampling — differing by exactly the 2.747 A pixel in quadrature.
- **The models carry line-list artifacts.** Predicted (K13) O I transitions to
  n = 15–16 Rydberg levels put absorption at 4403 and 4827 Å that is absent from
  the data — the upper levels are dissolved by the plasma microfield at
  photospheric density and cannot carry a line at all.
- **Reddening is the main selection risk**, and map columns are upper bounds
  only, overshooting by up to 8× for stars this close.
