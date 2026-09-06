# Data: libraries and sample

The observed spectra this project fits, where they come from, and which stars
survived selection. Traps specific to each library — wavelength conventions,
resolution, coverage gaps, peculiar stars — are in [CAVEATS.md](CAVEATS.md).

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

## NGSL resolution: use R = 600, not the STIS tables

The STIS LSF tables give the *single-exposure optical* profile: constant in
Angstroms within a grating (3.85 A for G430L, 8.09 A for G750L), implying
R = 804-1343. **The delivered NGSL v2 spectra are 1.7-1.8x broader than that,
and broader in a different functional form.**

Measured against XSL, which observes three of these stars at ~10x the
resolution — so no model, no NLTE, nothing but instrument
(`explore/ngsl_lsf_from_xsl.py`):

| lambda | measured FWHM | implied R | tabulated |
|---|---|---|---|
| 3900 A | 6.60 A | 591 | 3.85 A |
| 4400 A | 6.99 A | 629 | 3.85 A |
| 4900 A | 7.79 A | 629 | 3.85 A |
| 6600 A | 12.68 A | 521 | 8.09 A |
| 8700 A | 14.48 A | 601 | 8.09 A |

**R = 600 +/- 40**, constant in velocity, with no jump at the G430L/G750L
splice. Constant-R describes it with 7% scatter; constant-Angstrom needs 41%.
Three stars agree independently, and the jointly fitted wavelength shifts are
below 0.2 A in G430L, so the width is not absorbing a wavecal error. Both
spectra are continuum-normalized per window, so the grey flux offset between
the libraries cannot contribute either.

The likely cause is in the delivery rather than the optics: v2 spectra are
co-adds of two dithered exposures resampled onto a common grid, which broadens
the profile beyond the single-exposure LSF the tables describe.

Using the correct profile matters a great deal. Switching the model convolution
from 3.85 A to R = 600 cut the Balmer-region residual RMS by ~2.5x across the
sample and shrank every break residual:

| star | RMS before -> after | D residual before -> after |
|---|---|---|
| HD194453 | 6.7% -> 3.1% | +0.027 -> +0.010 |
| HD040573 | 5.0% -> 1.9% | +0.084 -> +0.058 |
| HD128801 | 6.2% -> 2.4% | +0.025 -> +0.020 |
| HD143459 | 7.3% -> 2.6% | -0.060 -> -0.007 |

`common.lsf.broaden_ngsl` uses the measured profile by default;
`tabulated=True` gives the STIS-table form.

Two earlier statements are superseded. **R = 939 at the Balmer break is wrong**
— it is ~600. And the claim that the data are *undersampled* (LSF narrower than
2 pixels) is wrong: at 6.6 A FWHM over 2.744 A pixels the spectra are sampled
at ~2.4 px per resolution element, which is adequate rather than aliased.

## The two libraries differ by a grey flux offset

NGSL and XSL are both absolutely calibrated, so their fluxes can be compared
directly. For the three stars in common the ratio XSL/NGSL is flat in
wavelength and star-dependent:

| star | XSL / NGSL | spread over 3300-9100 A |
|---|---|---|
| HD194453 | 1.04 | +/-2% |
| HD143459 | 0.90 | +/-2% |
| HD128801 | 1.00 | +/-3% |

Grey, so smoothing cannot be the cause — convolution conserves flux. It is an
absolute flux-calibration difference; ground-based XSL must correct slit
losses, NGSL need not. HD143459's 10% offset exceeds both libraries' quoted
accuracies (NGSL ~3%, XSL 1.5-4%), so at least one is worse than advertised for
that star. Comparisons should normalize; do not read an absolute flux
difference between the libraries as astrophysical.

## Provenance

| dataset | source |
|---|---|
| NGSL v2 | https://archive.stsci.edu/prepds/stisngsl/ (`explore/fetch_ngsl.sh`) |
| STIS LSFs | https://www.stsci.edu/hst/instrumentation/stis/performance/spectral-resolution |
| Pickles atlas | https://archive.stsci.edu/hlsps/reference-atlases/cdbs/grid/pickles/ |
| UVES-POP | https://sl.voxastro.org/library/UVES-POP/details/ (JSON API `/api/objects/UVES-POP/`) |
| MILES | Vizier J/MNRAS/371/703 |
| MaStar | Vizier J/ApJ/883/175 |
| Dust | IRSA SFD98/SF11 service; UVES-POP fitted E(B-V) |

Delivery documentation for NGSL is in `docs/ngsl_delivery/`. The spectra of the
selected stars are tracked in the repository; the full libraries are not, and
are refetchable with the scripts above.
