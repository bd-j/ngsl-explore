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
| [docs/FITTING.md](docs/FITTING.md) | fitting parameters, priors, degeneracies |
| [docs/CAVEATS.md](docs/CAVEATS.md) | **known issues and traps — read before trusting any number** |

Most of the work in this project turned out to be identifying ways the
comparison goes silently wrong: wavelength conventions, resolution mismatches,
reddening, peculiar stars, binaries, detector gaps. CAVEATS.md catalogues them
with symptoms and fixes.

## Layout

```
common/    shared: extinction (CCM89), break metric, LSF kernels, IO
grid/      model grid construction (make_model.py, build_grid.py, pack_grid.py)
fitting/   forward model and fitter (model.py, fit.py)
explore/   survey and comparison scripts
data/      catalogs, derived CSVs, selected-star spectra
docs/      this documentation, plus NGSL delivery docs in ngsl_delivery/
figures/   comparison and coverage figures
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

python3 grid/build_grid.py --workers 9   # the model grid (~25 h, resumable)
python3 grid/pack_grid.py                # collapse it into models/grid.npz

python3 explore/plot_model_vs_obs.py     # NGSL comparison figures
python3 explore/plot_uves_vs_model.py    # UVES-POP comparison figures
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
  lines the residual is ~1%, and D_Balmer agrees to 0.02-0.08 mag.
- The **Balmer line cores are systematically filled** relative to the LTE
  models, by ~10% of the line equivalent width, in every star. A
  flux-conservation test rules out a broadening mismatch. Most likely NLTE in
  hydrogen, which this code does not treat for H. The same signature appears in
  the Paschen lines.
- **NGSL is in air**, not vacuum, with a linear-in-lambda residual per grating
  that is recalibrated against the models.
- **NGSL's `STATERR` is optimistic by ~3x** and its resolution is R = 939 at the
  break, not the ~665 implied by pixel sampling.
- **Reddening is the main selection risk**: two of nine candidates were dropped
  on it alone.
