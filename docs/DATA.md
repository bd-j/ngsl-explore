# Data: libraries and sample

The observed spectra this project fits, where they come from, and which stars
survived selection. Traps specific to each library — wavelength conventions,
resolution, coverage gaps, peculiar stars — are in [CAVEATS.md](CAVEATS.md).

## Spectral libraries

| library | what it is | resolution | coverage | role here |
|---|---|---|---|---|
| **NGSL v2** | 379 HST/STIS spectra, space-based spectrophotometry | **R ~ 600 as delivered** (tables say 804-1343) | 1675-10198 A | primary sample; the only one reaching the UV |
| **XSL DR3** | 830 VLT/X-shooter spectra of 683 stars (Verro+2022) | R ~ 9800 UVB, ~11600 VIS | 3500-24800 A | **23 stars shared with NGSL** — independent check on the same objects |
| **UVES-POP** | 406 VLT/UVES spectra, re-reduced and flux-calibrated (2023) | R = 80,000 native, **~18,000 as delivered** | 3200-10250 A | high-resolution follow-up; resolves line cores |
| **MILES** | 985 ground-based spectra | FWHM 2.5 A | 3525-7500 A | independent Teff, [Fe/H] and E(B-V); 145 stars shared with NGSL |
| **Pickles** | 131 composite templates by spectral type | R ~ 500 | 1150-10620 A | reference break shape vs type and gravity; not individual stars |
| **MaStar** | MaNGA stellar library | R ~ 1800 | 3622-10354 A | checked, **zero overlap** with NGSL (disjoint in brightness) |

Model atmospheres are computed with [ATLAS12 + SYNTHE](https://github.com/cconroy20/atlas12)
(Kurucz, F90 translation), started from the C3K v2.3 grid.

## The sample

**The fitted sample is the NGSL ∩ XSL intersection, 9000–11000 K**, built by
`explore/build_sample.py` → `data/sample.csv`. 13 stars: 9 primary, 3 secondary
(binaries, kept and flagged), 1 rejected.

This supersedes the earlier NGSL + UVES-POP sample described further down.
Two reasons the definition changed:

* **UVES-POP is no longer fitted.** It has no overlap with NGSL, so it can only
  ever be a separate sample rather than a cross-check on the same star, and its
  continuum normalisation is too uncertain to contribute to a break measurement.
* **The intersection is the point.** Each library supplies something the other
  cannot, and the break test needs both at once:

  | library | what only it provides |
  |---|---|
  | NGSL | space-based spectrophotometry; the 3200–3500 Å continuum slope is measured from above the atmosphere |
  | XSL | ~16× the resolving power; locally normalised line profiles are **immune to reddening**, so they carry Teff, log g and v sin i free of the dust degeneracy |

  A degree-4 polynomial absorbs CCM89 across a 1100 Å window to 3×10⁻⁵, which is
  what makes the second statement exact rather than approximate. HD040573, a
  mainstay of the earlier NGSL-only comparison, drops out: XSL never observed it.

Selection, all enforced in code with the reason recorded and rejected rows kept:

| criterion | action |
|---|---|
| Teff 9000–11000 K, satisfied by **either** catalog | select |
| Ap/Bp/Am by SIMBAD type, or [Fe/H] ≥ +0.4 with no binary flag | **reject** |
| binary — SIMBAD type **or Gaia RUWE > 1.4** | flag → secondary |
| horizontal branch (`HB*`) | flag only, keep |
| outside the model grid in Teff, log g or [M/H] | flag, recorded per axis |

Either catalog satisfying the Teff window is deliberate: XSL and NGSL disagree
by up to 2900 K for the same star (HD164257: 10885 vs 7977), which is the whole
reason parameters are fitted here. Requiring agreement would let one catalog's
systematic define the sample.

`HB*` was previously treated as peculiarity and is not. A field horizontal-branch
star at 9000–11000 K sits below the ~11500 K Grundahl jump where radiative
levitation starts, so a scaled-solar atmosphere still describes it; what actually
disqualifies these stars is low log g and low [M/H], which the grid-coverage flag
records instead. This reclassification is why HD143459, HD074721 and HD128801 are
in the sample rather than rejected.

**Gaia RUWE caught two binaries that SIMBAD did not**: HD164967 (RUWE = 8.32) and
HD147550 (1.80), both of which have clean SIMBAD object types. That is the failure
mode CAVEATS records for HD162630, caught this time.

### The binding constraint is the grid's [M/H] floor

**8 of the 13 stars fall outside the model grid**, almost all in [M/H]: the grid
stops at −0.5 while the sample reaches −1.92. A node scan cannot extrapolate — it
piles up on the boundary and returns a wall, not a measurement — so this is a
hard limit on which stars produce quotable parameters, not a bias to be corrected.

| axis | grid | sample needs | stars affected |
|---|---|---|---|
| [M/H] | −0.5 … +0.3 | down to −1.92 | HD143459, HD074721, HD164967, HD117880, HD128801, HD106304, HD164257, HD072968 |
| log g | 3.0 … 5.0 | down to 2.84 | HD128801 |
| Teff | 8500 … 11500 | down to 7977 (NGSL value) | HD166991, HD164257 |

### Superseded: the earlier NGSL + UVES-POP sample

The tables below describe the sample used for the *catalog-parameter comparison*
figures, before the fit was set up. Kept because the figures in `figures/` and
the resolution and flux-offset measurements below were made with it.

#### NGSL — space-based spectrophotometry, R ~ 600 at the break

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

#### UVES-POP — high resolution, resolves the line cores

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

### XSL — the only library that overlaps NGSL star by star

XSL DR3 (Verro et al. 2022): 830 spectra of 683 stars, ground-based
VLT/X-shooter. 74 A stars in 7000-11500 K; 11 in the Balmer window, 8 clean
after vetting (`explore/xsl_astars.py`, parameters from Arentsen et al. 2019
since the DR3 table carries only names and filenames).

Its value here is the **overlap**. UVES-POP is disjoint from NGSL, but 23 of
the 74 XSL A stars are also in NGSL, including three of the four Balmer-break
targets. Same star, two instruments, one space-based and one ground-based at
~16x the resolution — which is what made the NGSL LSF measurement below
possible without a model.

| star | XSL Teff | NGSL Teff | Δ | in our sample |
|---|---|---|---|---|
| HD 194453 | 10489 | 10241 | +248 | yes |
| HD 147550 | 10044 | 10074 | −30 | dropped (reddened) |
| HD 143459 | 10689 | 9878 | +811 | yes |
| HD 128801 | 8774 | 10123 | −1349 | yes |

Those disagreements are the case for fitting rather than adopting catalog
parameters, made concrete: up to 1450 K in Teff and 0.9 dex in log g for the
same star, from two published analyses.

**Format facts, all verified rather than assumed:**

* **Air wavelengths.** Established by cross-correlating HD194453 — in both
  libraries — against the wavecal-corrected NGSL spectrum: +1.05 A required
  against a +1.13 A air-vacuum offset at 4000 A, closing to −0.05 A after
  conversion. `common/xsl_load.py` converts to vacuum by default.
* **Rest-frame**, so the RV is already removed. Fix `rv = 0` when fitting XSL,
  unlike NGSL.
* **Wavelengths in nm**, log-sampled at ~R = 30,000 (3 px per resolution
  element). The range is spectrum-dependent because the rest-frame shift differs
  per star.
* **Resolution is quoted as sigma(v), NOT FWHM**: 13 km/s UVB, 11 VIS, 16 NIR.
  FWHM = 2.3548 sigma, so R ~ 9800 at the Balmer break — misreading this gives an
  answer 2.35x wrong. It is constant in VELOCITY, unlike NGSL's, so the two need
  different convolution kernels.
* **Dereddened flux ships alongside raw** (`FLUX_DR`), plus variants for
  slit-loss correction (`_scl`, `_ncl`, `_ncge`). We use `FLUX` with our own
  CCM89 so extinction is handled identically across all three libraries, and
  treat `FLUX_DR` as a cross-check.

**Coverage gaps are real and must be masked**, not interpolated across: a
dichroic gap at 5750-5844 A, inter-order gaps every ~150 A redward of 8515 A,
and one star (HD162678) missing 3859-4779 A outright — 43% of its Balmer
coverage. See [CAVEATS.md](CAVEATS.md), which also has the echelle order-width
measurement explaining why the blue has no gaps and the red does.

#### Getting an XSL spectrum: extract it from the tarball, do not refetch

**The whole release is already on disk** as `data/xsl/XSL_DR3_release.tar`
(772,163,584 bytes, 832 members, **606 `_merged.fits` spectra**). It is
gitignored but it is not deleted after extraction, so adding a star to the
sample needs no download at all — only an extraction. This was missed once
already, and nearly cost a second 772 MB fetch.

Spectra are keyed by XSL ID (the `xslid` column of `data/sample.csv` or
`data/xsl_all.csv`), so one star is:

```bash
tar -xf data/xsl/XSL_DR3_release.tar -C data/xsl \
    XSL_DR3_release/xsl_spectrum_X0196_merged.fits
```

Extracting several at once is the same command with more members; passing the
member list is what keeps it from unpacking all 606 (12 MB vs 736 MB on disk).
Each merged spectrum is 928 KB.

**Verify the tarball by reading it, not by its size.** `data/xsl/fetch.sh` is a
resume loop because the server drops long connections: this copy needed **three
resumes** (264 MB → 465 MB → 722 MB → complete), and a single curl exits
cleanly on a truncated file. `tar -tf` returning exit 0 with 832 members is the
check that the archive is whole; `ls -l` is not.

Only the spectra of the *selected* stars are extracted and tracked, by name, so
the tracked set cannot drift out of step with the sample. 13 are currently
extracted: the 12 non-rejected sample stars, plus **X0288, a second XSL epoch of
HD194453**. XSL observed the primary target twice (X0196 and X0288) and
`build_sample.py` takes one epoch per star, so X0288 is not in `sample.csv` — but
two independent observations of the same star with the same instrument are a
repeatability check on the XSL continuum and on v sin i, which is worth having
for the one star the whole analysis leans on.

## Gaia DR3: the independent dust lever

`explore/fetch_gaia.py` → `data/gaia_sample.csv` and `data/gaia_xp/<star>_xp.csv`.

The break test needs E(B−V) *measured*, from data statistically independent of
the break itself. Gaia supplies that from space, over almost exactly the NGSL
range, with no slit and no atmosphere. What matters is the lever arm:

| baseline | differential extinction per 0.01 mag E(B−V) |
|---|---|
| NGSL 3200–3500 Å window | 0.0035 mag |
| NGSL full 3300–9400 Å | 0.0365 mag |
| **Gaia XP 3360–10200 Å** | **0.0358 mag** |

XP sampled spectra are 343 points, 336–1020 nm at 2 nm, externally calibrated to
~1–2%. Retrieved for 11 of 13 stars; HD143459 (V = 5.53) and HD174240 have no XP.

**Verified, not assumed:** XP against NGSL for HD194453 gives a median flux ratio
of **1.011** over 4000–9000 Å — two independent space observatories agreeing to
1.1% in absolute flux, which validates both calibrations and the unit conversion.
The 10.7% scatter about that is XP's very low resolution smearing lines
differently, not a calibration disagreement; the ratio at 4000 Å drops to 0.91
because Hδ/Hγ sit there.

**Units.** XP flux is W m⁻² nm⁻¹; NGSL and the models are erg s⁻¹ cm⁻² Å⁻¹. The
factor is 10², and the wavelengths are **vacuum** — Gaia has no air path, so
unlike NGSL, XSL and UVES-POP no conversion is needed.

**Two sources, deliberately.** VizieR (`I/355/gaiadr3`) for the cone search, ESA
DataLink for the XP spectra. ESA's TAP endpoint answered a single cone search in
115 s and then degraded to failing three retries in a row; VizieR answers the
same query in 1.1 s. VizieR also returns `Source` as int64, so the 19-digit
source_id survives exactly — through ESA TAP it arrived float64-rounded, and
`fetch_gaia.py` now checks DataLink's echoed id against the requested one so a
mismatch cannot pass silently.

Caveats carried per star rather than assumed away: Gaia photometry and XP
calibration degrade below G ≈ 6 (three stars), so `XPcont`, `XPsamp`, G and the
BP/RP excess factor are all stored; and the cone is kept at 5″ with the
counterpart chosen by brightness (for an A star G ≈ V to ~0.1 mag) because the
library coordinates are epoch 2000 against DR3's 2016.0. Observed offsets came
out 0.1–2.7″.

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

## NGSL and XSL differ by a grey flux offset

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
| Gaia DR3 | VizieR I/355/gaiadr3 (photometry); ESA DataLink XP_SAMPLED (spectra) — `explore/fetch_gaia.py` |
| Dust | IRSA SFD98/SF11 service; UVES-POP fitted E(B-V) |

Delivery documentation for NGSL is in `docs/ngsl_delivery/`. The spectra of the
selected stars are tracked in the repository; the full libraries are not, and
are refetchable with the scripts above.
