# Balmer break: empirical spectra vs model atmospheres

Testing how well model stellar atmospheres reproduce the **Balmer break** — the
continuum discontinuity at 3646 A where bound-free absorption from hydrogen n=2
cuts off — in real stars near 10,000 K, where the break is strongest.

The approach is to pick well-calibrated observed spectra of A stars, compute
ATLAS12 model atmospheres and SYNTHE spectra at those stars' parameters, and
compare the continuum shape across the break. The Paschen break (8206 A) comes
along for free and serves as a consistency check.

**Read [CAVEATS.md](CAVEATS.md) before trusting any number here.** Most of the
work in this project turned out to be identifying ways the comparison can go
silently wrong — wavelength conventions, resolution mismatches, reddening,
peculiar stars, binaries. They are catalogued there with symptoms and fixes.

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

## Spectral libraries

| library | what it is | resolution | coverage | role here |
|---|---|---|---|---|
| **NGSL v2** | 379 HST/STIS spectra, space-based spectrophotometry | R = 804-1343 (LSF-measured) | 1675-10198 A | primary sample; the only one reaching the UV |
| **UVES-POP** | 406 VLT/UVES spectra, re-reduced and flux-calibrated (2023) | R = 80,000 native, **~18,000 as delivered** | 3200-10250 A | high-resolution follow-up; resolves line cores |
| **MILES** | 985 ground-based spectra | FWHM 2.5 A | 3525-7500 A | independent Teff, [Fe/H] and E(B-V); 145 stars shared with NGSL |
| **Pickles** | 131 composite templates by spectral type | R ~ 500 | 1150-10620 A | reference break shape vs type and gravity; not individual stars |
| **MaStar** | MaNGA stellar library | R ~ 1800 | 3622-10354 A | checked, **zero overlap** with NGSL (disjoint in brightness) |

Model atmospheres are computed with [ATLAS12 + SYNTHE](https://github.com/cconroy20/atlas12)
(Kurucz, F90 translation), started from the C3K v2.3 grid.

## The sample

Seven stars near 10,000 K across two libraries, selected for clean model
comparison. Selection is enforced in code, not by hand: `scripts/candidate_table.py`
for NGSL and `scripts/uves_pop_load.py` for UVES-POP, both applying the same
`E(B-V) <= 0.10` cut.

### NGSL — space-based spectrophotometry, R ~ 940 at the break

| star | Teff | log g | [M/H] | E(B-V) | slit offset | notes |
|---|---|---|---|---|---|---|
| HD194453 | 10241 | 3.9 | +0.0 | -0.01 | **0.00 px** | primary target |
| HD040573 | 10200 | 4.2 | -0.4 | 0.06 | 0.36 px | classification and gravity agree |
| HD128801 | 10123 | 3.7 | -1.9 | 0.03 | 0.58 px | metal-poor comparison; in MILES |
| HD143459 | 9878 | 3.6 | -0.6 | 0.04 | 0.09 px | horizontal-branch star |
| ~~HD147550~~ | 10074 | 3.9 | -0.0 | **0.125** | 0.24 px | dropped: reddened |

HD194453 is the primary target because its **slit offset is 0.00 px** — the
wavelength-dependent slit-throughput correction, the dominant systematic on
break *shape*, is essentially null for it.

### UVES-POP — high resolution, resolves the line cores

| star | Teff | log g | [Fe/H] | v sin i | E(B-V) | S/N | notes |
|---|---|---|---|---|---|---|---|
| HD162678 | 9908 | 3.53 | +0.03 | 37 | 0.077 | 169 | slowest rotator; best for profiles |
| HD188294 | 11016 | 4.04 | +0.05 | 182 | 0.040 | 299 | best S/N and lowest reddening |
| HD162393 | 9955 | 4.05 | -0.55 | 142 | 0.065 | 132 | metal-poor comparison |
| ~~HD162817~~ | 10153 | 3.66 | +0.11 | 65 | **0.102** | 176 | dropped: reddened |
| ~~HD162630~~ | 10494 | 3.92 | -0.18 | 46 | 0.01 | 226 | dropped: spectroscopic binary |

UVES-POP spectra are **dereddened** with CCM89 (R_V = 3.1) using the library's
own fitted E(B-V), and models are rotationally broadened with each star's
published `v sin i`. Both libraries were screened against Ap/Am peculiarity and
binarity — see [CAVEATS.md](CAVEATS.md), since neither is caught by parameter
cuts alone.

The two samples are **disjoint** — no star appears in both — so they are
independent tests rather than a repeat measurement. Between them they span
9878-11016 K, log g 3.5-4.2, [M/H] -1.9 to +0.1, and rotation 37-182 km/s.

## Pipeline

```bash
export ATLAS12=/path/to/atlas12

./scripts/fetch_ngsl.sh                  # NGSL, docs, STIS LSFs, Pickles atlas
python3 scripts/extract_docs.py          # PDF text (parameters live in the readme)
python3 scripts/build_catalog.py         # merge headers + params + magnitudes
python3 scripts/measure_snr.py           # empirical S/N, two estimators
python3 scripts/crossmatch_libraries.py  # vs MILES and MaStar
python3 scripts/candidate_table.py       # candidate selection + reddening cut
python3 scripts/reddening.py             # E(B-V): map, photometric, fitted

python3 scripts/uves_pop_astars.py       # select A stars from UVES-POP
python3 scripts/uves_pop_load.py         # parameters + reddening cut

python3 scripts/make_atlas_model.py --all      # atmospheres + spectra (~8 min/star)
python3 scripts/plot_model_vs_obs.py           # NGSL comparison figures
python3 scripts/plot_uves_vs_model.py          # UVES-POP comparison figures
```

Making one model, given `(Teff, log g, [Fe/H])`:

```bash
python3 scripts/make_atlas_model.py --star HD194453 --teff 10241 --logg 3.9 --feh 0.0
```

which runs ATLAS12 from the nearest C3K grid model, synthesizes at R = 300,000
(SYNTHE's enforced floor), and convolves to R = 10,000. Note `zscale = 10**[Fe/H]`
and that the `.spec` flux column is Eddington H_nu, not f_lambda.

## Findings so far

- The **continuum across the break is reproduced well**: away from hydrogen
  lines the residual is ~1%, and D_Balmer agrees to 0.03-0.07 mag.
- The **Balmer line cores are systematically filled** relative to the LTE
  models, by ~10% of the line equivalent width, in every star. A
  flux-conservation test rules out a broadening mismatch. Most likely NLTE in
  hydrogen, which this code does not treat for H. The same signature appears in
  the Paschen lines.
- **Reddening is the main selection risk.** Unmodelled extinction depresses the
  blue side and inflates the measured break, mimicking a model failure. Two of
  the nine candidates were dropped on this alone.
- **NGSL's `STATERR` is optimistic by ~3x** and its resolution is R = 939 at the
  break, not the ~665 implied by pixel sampling — the CCD LSF is narrower than
  Nyquist.

## Layout

```
scripts/     pipeline, one step per file
data/        catalogs and derived CSVs (bulk spectra are gitignored, refetchable)
docs/        NGSL delivery documentation
figures/     comparison and coverage figures
models/      ATLAS12 atmospheres and synthesized spectra
report.html  NGSL library survey        pickles.html  Pickles atlas survey
CAVEATS.md   known issues and traps
```
