# Data: libraries and sample

The observed spectra this project fits, where they come from, and which stars
survived selection. Traps specific to each library — wavelength conventions,
resolution, coverage gaps, peculiar stars — are in [CAVEATS.md](CAVEATS.md).

## Spectral libraries

| library | what it is | resolution | coverage | role here |
|---|---|---|---|---|
| **NGSL v2** | 379 HST/STIS spectra, space-based spectrophotometry | **R ~ 600 as delivered** (tables say 804-1343) | 1675-10198 A | primary sample; the only one reaching the UV |
| **XSL DR3** | 830 VLT/X-shooter spectra of 683 stars (Verro+2022) | R ~ 9800 UVB, ~11600 VIS | 3500-24800 A | **23 stars shared with NGSL** — independent check on the same objects |
| **UVES-POP** | 406 VLT/UVES spectra, re-reduced and flux-calibrated (2023) | R = 80,000 native, **~18,000 as delivered** | 3200-10250 A | high-resolution follow-up; resolves line cores. 13 stars shared with NGSL, **9 with XSL**, none in the Balmer window |
| **MILES** | 985 ground-based spectra | FWHM 2.5 A | 3525-7500 A | independent Teff, [Fe/H] and E(B-V); 145 stars shared with NGSL |
| **Pickles** | 131 composite templates by spectral type | R ~ 500 | 1150-10620 A | reference break shape vs type and gravity; not individual stars |
| **MaStar** | MaNGA stellar library | R ~ 1800 | 3622-10354 A | checked, **zero overlap** with NGSL (disjoint in brightness) |

Model atmospheres are computed with [ATLAS12 + SYNTHE](https://github.com/cconroy20/atlas12)
(Kurucz, F90 translation), started from the C3K v2.3 grid.

## The sample

**The fitted sample is the NGSL ∩ XSL intersection, 9000–11000 K**, built by
`explore/build_sample.py` → `data/sample.csv`. 13 stars: 9 primary, 3 secondary
(binaries, kept and flagged), 1 rejected.

HD194453 is the worked example throughout these documents because its NGSL
**slit offset is 0.00 px**: the wavelength-dependent slit-throughput correction
is the dominant systematic on break *shape*, and for this star it is
essentially null (see [CAVEATS.md](CAVEATS.md) for the general case).

**The intersection is the point.** Each library supplies something the other
cannot, and the break test needs both at once:

  | library | what only it provides |
  |---|---|
  | NGSL | space-based spectrophotometry; the 3200–3500 Å continuum slope is measured from above the atmosphere |
  | XSL | ~16× the resolving power; locally normalised line profiles are **immune to reddening**, so they carry Teff, log g and v sin i free of the dust degeneracy |

  A degree-4 polynomial absorbs CCM89 across a 1100 Å window to 3×10⁻⁵, which is
  what makes the second statement exact rather than approximate. Requiring both
  libraries is what sets the sample size: NGSL alone has many more A stars.

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

**`HB*` is not a peculiarity type here.** A field horizontal-branch star at
9000–11000 K sits below the ~11500 K Grundahl jump where radiative levitation
starts, so a scaled-solar atmosphere still describes it. What disqualifies these
stars is low log g and low [M/H], which the grid-coverage flag records instead
— which is why HD143459, HD074721 and HD128801 are in the sample rather than
rejected.

**Gaia RUWE caught two binaries that SIMBAD did not**: HD164967 (RUWE = 8.32) and
HD147550 (1.80), both of which have clean SIMBAD object types. That is the failure
mode CAVEATS records for HD162630, caught this time.

### The sample, star by star

Generated from `data/sample.csv`; catalog parameters as delivered, not fitted.

| star | tier | Teff NGSL / XSL | log g | [M/H] NGSL / XSL | SpType | slit px | flags | grid |
|---|---|---|---|---|---|---|---|---|
| HD167946 | primary | 10634 / 10079 | 4.30 | -0.10 / -0.47 | A0 | +0.63 | — | in |
| **HD194453** | primary | 10241 / 10489 | 3.90 | +0.00 / -0.02 | A0III | +0.00 | — | in |
| HD128801 | primary | 10123 / 8774 | 3.70 | -1.90 / -1.92 | B9 | +0.58 | HB | below |
| HD143459 | primary | 9878 / 10689 | 3.60 | -0.60 / -0.19 | A0V | +0.09 | HB | in † |
| HD117880 | primary | 9426 / 8843 | 3.70 | -0.60 / -1.65 | B9IV/V | **+0.91** ‡ | — | below |
| HD106304 | primary | 9376 / 8759 | 3.60 | -1.80 / -1.67 | B9V | -0.33 | — | below |
| HD174240 | primary | 9274 / 9262 | 3.80 | -0.20 / -0.40 | A1IV | +0.28 | map | in |
| HD074721 | primary | 8774 / 9571 | 3.30 | -0.60 / -0.57 | A0V | -0.35 | HB | below |
| HD166991 | primary | 8497 / 9008 | 4.00 | -0.30 / -0.20 | A1V | +0.45 | — | in |
| HD147550 | secondary | 10074 / 10044 | 3.90 | +0.00 / -0.26 | B9V | +0.24 | **B** (RUWE 1.80) | in |
| HD164967 | secondary | 8534 / 9351 | 4.10 | -0.60 / -0.29 | A0 | +0.37 | **B** (RUWE 8.32) | in † |
| HD164257 | secondary | 7977 / 10885 | 3.50 | -0.10 / +0.73 | A0 | -0.36 | **B** (El\*), Z+ | in |
| ~~HD072968~~ | rejected | 9253 / 9569 | 3.90 | +0.50 / +0.44 | A1VpSrCr | -0.17 | **pec**, Z+ | in |

**flags** — **B** binary, with the evidence that caught it (Gaia RUWE > 1.4, or
the SIMBAD object type); **HB** field horizontal branch, flagged but kept;
**pec** chemically peculiar, the one rejection; **Z+** [Fe/H] ≥ +0.4;
**map** the SFD/SF11 dust column is unusable at this star's galactic latitude.

**slit px** — the NGSL v2 slit offset, the star's miscentring in the 52×0.2
aperture. It is the proxy for the dominant systematic on break *shape*: the
wavelength-dependent slit-throughput correction grows with it, and the v2
correction is only reliable **below 0.9 px** ([CAVEATS.md](CAVEATS.md)). The
sign is the direction of the offset; the magnitude is what matters here.

‡ marks the one star NGSL itself flags `dataqual = suspect`. That flag is
exactly the 0.9 px cut: across all 379 catalog stars the 35 `suspect` rows span
|offset| 0.905–1.165 px and the 344 `good` rows span 0.000–0.895, with no
exception either way. HD117880 at +0.91 px is therefore only just over the line,
and it is already segregated for being below the grid's [M/H] floor.

**grid** — whether any node of the model grid can represent the star, which
decides where its figures go (`common/figpath.py`: `figures/fits_mh_in_grid/` or
`figures/fits_mh_below_grid/`). † marks the two stars the catalog cut sends
below the floor but whose metal lines are fit at the −0.5 node; the exception is
a reviewed list with a reason per star, not a threshold.

Note that `grid` and the per-axis out-of-grid flag in the table above answer
different questions: `grid` is about the [M/H] **floor** only, while HD164257 is
flagged out-of-grid because XSL puts it at +0.73, above the +0.3 **ceiling**.

### The binding constraint is the grid's [M/H] floor

**9 of the 13 stars fall outside the model grid on some axis, 8 of them in
[M/H]** — and the floor is what binds: the grid stops at −0.5 while the sample
reaches −1.92. A node scan cannot extrapolate — it piles up on the boundary and
returns a wall, not a measurement — so this is a hard limit on which stars
produce quotable parameters, not a bias to be corrected.

A star counts as outside if **either** catalog puts it outside, on the same
either-catalog principle as the Teff window.

| axis | grid | sample reaches | stars affected |
|---|---|---|---|
| [M/H] below floor | −0.5 | down to −1.92 | HD143459, HD074721, HD164967, HD117880, HD128801, HD106304 |
| [M/H] above ceiling | +0.3 | up to +0.73 | HD164257, HD072968 |
| log g | 3.0 … 5.0 | down to 2.84 | HD128801 |
| Teff | 8500 … 11500 | down to 7977 (NGSL value) | HD166991, HD164257 |

The two above the ceiling are a different problem from the six below the floor:
HD072968 is the rejected chemically peculiar star, and HD164257's +0.73 is the
XSL value for a star NGSL puts at −0.10 — the 2900 K catalog disagreement again,
not a real super-metal-rich A star.

### XSL — the library that overlaps NGSL in the Balmer window

XSL DR3 (Verro et al. 2022): 830 spectra of 683 stars, ground-based
VLT/X-shooter. 74 A stars in 7000-11500 K; 11 in the Balmer window, 8 clean
after vetting (`explore/xsl_astars.py`, parameters from Arentsen et al. 2019
since the DR3 table carries only names and filenames).

Its value here is the **overlap**. UVES-POP overlaps NGSL in only 13 stars and
none of them is an A star in the Balmer window, but 23 of the 74 XSL A stars
are also in NGSL, including three of the four Balmer-break targets. Same star,
two instruments, one space-based and one ground-based at ~16x the resolution —
which is what made the NGSL LSF measurement below possible without a model.

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

**XSL's coverage is continuous across the optical.** Checked directly: 486
points in 5750-5844 A and 610 in 8515-8690 A. Its only gaps wider than 3 A are
a handful in the NIR telluric regions near 18,400-19,150 A, star-dependent and
outside the range used here. The limit that does bite is XSL's **3501 A start**
— only 145 A blueward of the break, which is what sets the blue edge of the
comparison panels in `explore/plot_uves_xsl.py`.

The dichroic gap at 5750-5844 A, the inter-order gaps redward of 8515 A and
HD162678's 3859-4779 A hole belong to **UVES-POP**, not XSL; they are listed
under UVES-POP in [CAVEATS.md](CAVEATS.md), which is also where the
gap-handling rule lives.

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

### UVES-POP × NGSL: 13 stars in common, none in the Balmer window

`explore/uves_ngsl_overlap.py` → `data/uves_ngsl_overlap.csv`. Neither catalog
carries the other's identifiers and UVES-POP names its brightest targets by
proper name, so the match is positional: SIMBAD resolves the 406 UVES-POP names
and those coordinates are matched against the NGSL v2 header positions.

The 13 pairs all fall within 7.3", the next-nearest pair is at 36" and the one
after that at 1952", so the 10" cut is in a wide gap and the answer does not
depend on where in it the line is drawn. The 36" pair is **HD36959 (UVES-POP) /
HD36960 (NGSL)** — two components of a visual double, not one star seen twice,
and correctly excluded. 49 UVES-POP entries carry open-cluster running numbers
that SIMBAD does not resolve (46 in IC 2391, 3 in NGC 6475); the nearest NGSL
star to either field is 6.4° and 12.2° away, so none of them can be a match.

The radius is 10" rather than the 5" `build_sample.py` uses for XSL because
these are nearby bright stars and the leftover separation is **proper motion**
between the NGSL header epoch and SIMBAD's J2000 positions. The four widest
matches are exactly the four fastest movers — 171 Pup 1.72"/yr at 7.33",
61 Vir 1.51"/yr at 3.20", ε Eri 0.98"/yr at 1.75", 10 Tau 0.54"/yr at 2.11" —
which is what makes them epoch offsets rather than doubtful matches. At 5"
171 Pup would have been dropped for moving.

| star | SpT | V | Teff NGSL | Teff UVES | Δ | log g N / U | [M/H] N / [Fe/H] U | v sin i |
|---|---|---|---|---|---|---|---|---|
| HD022049 | K2V | 3.73 | 5130 | 5342 | 212 | 4.50 / 5.04 | 0.00 / -0.14 | 0.3 |
| HD022484 | F9IV-V | 4.28 | 6141 | 5866 | -275 | 4.10 / 3.60 | 0.10 / -0.45 | 1.6 |
| HD047839 | O7Ve | 4.66 | — | 45541 | — | — / 3.96 | — / 0.34 | 92.4 |
| HD058343 | B2Vne | 5.20 | — | — | — | — / — | — / — | — |
| HD063077 | G0V | 5.37 | 5926 | 5795 | -131 | 4.20 / 4.02 | -0.70 / -0.95 | 1.8 |
| HD076932 | F7-8IV-V | 5.86 | 6034 | 5852 | -182 | 4.10 / 3.59 | -0.70 / -0.96 | 4.1 |
| HD099648 | G8Iab: | 4.95 | 4811 | 4990 | 179 | 2.00 / 1.89 | -0.30 / -0.27 | 9.1 |
| HD102212 | M1III | 4.05 | 3800 | 3982 | 182 | 1.10 / 1.14 | 0.10 / -0.42 | 11.0 |
| HD111786 | A0III | 6.14 | 7598 | 7436 | -162 | 3.90 / 4.03 | -1.30 / -1.70 | 47.0 |
| HD115617 | G5V | 4.74 | 5557 | 5755 | 198 | 4.30 / 4.70 | 0.00 / -0.07 | 3.9 |
| HD138716 | K1IV | 4.61 | 4771 | 5214 | 443 | 2.90 / 3.63 | -0.10 / 0.05 | 9.2 |
| HD142703 | A2Ib/II | 6.12 | 7337 | 7285 | -52 | 3.90 / 3.85 | -1.40 / -1.77 | 92.5 |
| HD206778 | K2Ib | 2.40 | 4095 | 4483 | 388 | 1.30 / -0.05 | 0.10 / -0.32 | 10.5 |

Teff from NGSL (Castelli 2004 on Victoria-Regina isochrones, fitted to v1) and
from UVES-POP (VOXAstro PHOENIX fit) agree to within 443 K on the 11 stars both
libraries fitted, but log g and [M/H] do not: up to 1.35 dex in log g (HD206778)
and 0.55 dex in metallicity (HD022484). The two parameter scales are not
interchangeable, which is the same lesson the XSL comparison gives above.

Four rows carry a reason in the `notes` column: HD047839 (O7Ve) and HD058343
(B2Vne) are outside the range NGSL attempted a fit for, HD058343 is unfitted in
UVES-POP too, HD138716 is `dataqual = suspect` in NGSL, and HD206778's UVES-POP
log g of −0.05 is pinned to a grid edge rather than measuring a K2Ib supergiant
— take that one as a fit artefact, not a gravity.

**The library's `v` column is the radial velocity in km/s, not the V
magnitude.** The script checks this against SIMBAD on every run: across the 12
stars here it tracks the catalogued RV to within 3.7 km/s (worst case 15 Mon, a
spectroscopic binary) while missing V by up to 114 mag. It is `rv_uves_kms` in
the overlap table.

Why this overlap does nothing for the Balmer break: UVES-POP is a southern
VLT library of bright stars, and the 13 shared stars are one of nearly every
type — O7, B2, A0, A2, two F, G0, G5, G8, K1, K2 V, K2 Ib, M1. The nearest
things to A stars are HD111786 (7598 K) and HD142703 (7337 K), both well below
the 9000–11000 K window, and the only hotter star with a fit is 15 Mon at
45,500 K. For the Balmer break the useful overlap is XSL's, above.

`explore/plot_uves_ngsl.py` compares the two libraries directly for all 13 —
UVES-POP smoothed to the NGSL LSF and integrated onto NGSL pixels through
`common.lsf.to_ngsl_pixels`, one panel per star with a fractional residual,
plus an H-epsilon zoom — writing `figures/explore_libraries/uves_ngsl_break.png`,
`uves_ngsl_hepsilon.png` and `data/uves_ngsl_compare.csv`. Three results from
it are worth having here rather than only in its docstring (the `v` column
finding above came from the same script):

* **UVES-POP's delivered spectra are in the OBSERVED frame.** Verified on the
  two fastest stars: the Ca II H and K minima sit at rest×(1+v/c) for the
  catalogued v, recovering +123.0/+122.6 km/s for HD076932 (catalog +119.76)
  and +107.7/+107.5 for HD063077 (+106.93). The catalog RV must be removed.
* **The bluest UVES setting join is 3733.6–3859.2 Å**, measured three
  independent ways, and there is none blueward of it. Detail in the script.
* **NGSL's wavelength zero point is fitted per star** here, against the smoothed
  UVES spectrum, because `data/ngsl_wavecal.csv` covers none of these 13. The
  fitted shifts run −13.5 to −70.6 km/s and absorb NGSL's own RV and its
  wavecal residual together, which nothing in this comparison can separate.

Gaia XP exists for 5 of the 13 and is available behind `--xp`, but is off by
default: three of the five give the same −1.5% blue-to-red tilt across the
break, which is the instrumental XP/NGSL pattern measured in
`explore/xp_vs_ngsl.py`, not a property of the stars.

Those two would be the wrong stars anyway. SIMBAD types them `F0VkA1mA1_lB`
and `F1VskA1.5mA1.5_lB` — **λ Boo stars**, and δ Scuti pulsators besides — so
the [M/H] ≈ −1.3 to −1.8 both catalogs report for them is a surface accretion
pattern, not a halo metallicity, and their continua are not those of a normal
A star. Note also how far the NGSL table's `A0III` and `A2Ib/II` are from
SIMBAD's F0V/F1V: the `sptype` column is SIMBAD-as-of-2012 and should not be
trusted for luminosity class on these two.

### UVES-POP × XSL: 9 stars in common, and none in the Balmer window either

`explore/uves_xsl_overlap.py` → `data/uves_xsl_overlap.csv`. Same positional
method as the NGSL overlap above — SIMBAD resolves the UVES-POP names, and
those coordinates are matched against the XSL DR3 positions in
`data/xsl_all.csv`, which XSL repeats once per epoch and the script collapses
to 683 unique stars.

Both are VLT libraries, so there is no hemisphere penalty here — and the
overlap is **smaller** than NGSL's 13 anyway. The two programs simply chose
nearly disjoint targets: UVES-POP took bright nearby stars for a
high-resolution atlas, XSL a stellar-population grid weighted to cool giants,
the bulge and the Magellanic Clouds.

The 9 pairs all fall within **0.09"** and the next-nearest pair is at 430"
(HD109379 / HD 109443 — β Crv and a different star), so the cut sits in a
four-decade gap and the answer does not depend on where in it the line is
drawn. The radius is the 5" `build_sample.py` uses for XSL, not the NGSL
overlap's 10", and the sub-0.1" residuals are the reason: **both** catalogs
here carry J2000 catalog positions, so there is no proper motion to absorb,
unlike the NGSL header positions at the epoch of observation.

The 49 open-cluster entries are **resolved, not dropped**. SIMBAD knows them
under their `Cl*` designations — all 3 NGC 6475 stars as `Cl* NGC 6475 JJnn`,
and 10 of the 46 IC 2391 stars as `Cl* IC 2391 SHJM n` or `PP n` — which takes
the resolved count from 357 to 370. None of the 13 matches anything. For the
36 IC 2391 entries still unresolved the script does not assume: it asserts that
XSL has **no star within 3° of the field** (nearest is 6.54°, and the resolved
members sit 0.20–0.64° from the centre), so the assertion fails loudly if a
future XSL release lands there. This is worth flagging because
`uves_ngsl_overlap.py` makes the same argument for NGSL from a 6.4°/12.2°
margin — but for XSL the NGC 6475 field is only **1.34°** away, close enough
that the NGSL script's reasoning could not simply be carried over.

The positional match is confirmed a second, independent way, on every run:
resolving the XSL identifier as well and requiring the two SIMBAD `main_id`s
to be the same object. **9/9 agree.**

| XSL | UVES-POP | SIMBAD | ep | Teff X / U | Δ | log g X / U | Δ | [Fe/H] X / U | Δ | v sin i |
|---|---|---|---|---|---|---|---|---|---|---|
| HD 39801 | Betelgeuse | α Ori | 1 | 3654 / 3779 | +125 | 0.43 / −0.46 | −0.89 | −0.21 / −0.13 | +0.08 | 14.3 |
| BS 4517 | HD102212 | ν Vir | 1 | 3733 / 3982 | +249 | 1.43 / 1.14 | −0.29 | −0.54 / −0.42 | +0.12 | 11.0 |
| HD 99648 | HD099648 | τ Leo | 1 | 4933 / 4990 | +57 | 2.18 / 1.89 | −0.29 | −0.04 / −0.27 | −0.23 | 9.1 |
| HD 24616 | HD024616 | HD 24616 | 1 | 5020 / 5225 | +205 | 3.35 / 3.47 | +0.12 | −0.74 / −0.71 | +0.03 | 3.7 |
| HD 140283 | HD140283 | HD 140283 | 2 | 5718 / 5712 | −6 | 3.66 / 2.06 | **−1.60** | −2.43 / −3.04 | −0.61 | 0.5 |
| HD 16673 | HD016673 | HD 16673 | 1 | 6189 / 6235 | +46 | 4.27 / 4.43 | +0.16 | −0.05 / −0.11 | −0.06 | 23.3 |
| HD 84937 | HD084937 | HD 84937 | 1 | 6212 / 6390 | +178 | 4.07 / 3.58 | −0.49 | −2.03 / −2.07 | −0.04 | 8.7 |
| HD 142703 | HD142703 | HR Lib | 1 | 7219 / 7285 | +66 | 4.24 / 3.85 | −0.39 | −1.16 / −1.77 | −0.61 | 92.5 |
| HD 111786 | HD111786 | MO Hya | 1 | 7595 / 7436 | −159 | 4.28 / 4.03 | −0.25 | −1.01 / −1.70 | −0.69 | 47.0 |

Δ is **UVES-POP minus XSL** throughout; the script asserts that sign against an
injected pair on every run. `ep` is the number of XSL epochs — HD 140283 has
two (X0687, X0688), a free repeatability check.

The parameter scales behave the same way they do in the NGSL comparison:
**Teff agrees, gravity and metallicity do not.** Teff runs +66 K median with a
72 K MAD and a −159 to +249 K range, and the two widest are the two coolest
stars (ν Vir, Betelgeuse) where both grids are extrapolating. Against that,
log g disagrees by up to 1.60 dex and [Fe/H] by 0.69 dex.

Two rows should not be used without reading why:

* **HD 140283** is the outlier and the one to check before use. The two
  libraries agree on Teff to **6 K** and then disagree by **1.60 dex** in
  gravity — that is subgiant versus giant, not a small-print difference — and
  by 0.61 dex in [Fe/H]. Something is wrong on one side; the table records both
  rather than picking one.
* **Betelgeuse** carries `uves_logg=-0.46_at_grid_edge` in the `notes` column.
  Like HD206778 in the NGSL table, that value is pinned to a grid edge and is a
  fit artefact, not a gravity.

The [Fe/H] gaps on HD 142703 and HD 111786 are **not** a scale offset: both are
λ Boo stars, and their low metallicity is a surface accretion pattern rather
than a composition — see the sample notes above, where the same two stars are
discussed for the NGSL overlap.

For four stars UVES-POP's own headers carry a literature value (`LIT_*`) as a
third opinion, and only two of those have `LIT_LOGG`: on HD099648 XSL is
nearer (2.18 vs lit 2.17, UVES-POP 1.89) and on HD142703 UVES-POP is
(3.85 vs lit 3.89, XSL 4.24). **One each way, n = 2** — that settles nothing
about which scale is better and is recorded here only so the next person does
not redo it.

**Four of the nine are in all three libraries** — HD099648, HD102212, HD111786
and HD142703 also appear in `data/uves_ngsl_overlap.csv`. Those are the stars
where NGSL, XSL and UVES-POP can be put on one plot at three resolutions.

Why this overlap, too, does nothing for the Balmer break: the hottest star in
it is **7595 K**, below the 9000–11000 K window, and the rest run down to
Betelgeuse at 3654 K. As with UVES-POP × NGSL, the useful Balmer-window
overlap remains NGSL ∩ XSL.

#### Comparing the two flux calibrations

`explore/plot_uves_xsl.py` → `figures/explore_libraries/uves_xsl_break.png`,
`uves_xsl_hepsilon.png` and `data/uves_xsl_compare.csv`. Both libraries are
absolutely calibrated, so this is a data-to-data comparison with no model in it.

**The smoothing runs UVES-POP → XSL, which is the only possible direction.**
At the break XSL is R ~ 9800 (σ_v = 13 km/s — XSL quotes σ, not FWHM) while
UVES-POP delivered is R ~ 18,000 and R = 80,000 native. UVES-POP is the sharper
by ~1.9×, so smoothing XSL to UVES-POP would be a deconvolution.

The kernel is the quadrature difference and is applied **at constant velocity**,
because XSL's LSF is (NGSL's is constant in Å — the two need different kernels).
σ_kernel = √(13² − σ_UVES²) = **12.7 km/s**, where σ_UVES combines the native
R = 80,000 (1.59 km/s) with the 0.1 Å delivered pixel (2.37 km/s rms). UVES-POP's
own resolution is therefore a **2.5% correction** to the kernel and nothing here
rests on it. Taking the 2-pixel figure (16.4 km/s) instead would give a 25.8 km/s
kernel — but that number is a *sampling* limit, not an LSF, and using it as one
over-broadens by 2×. The script self-tests the kernel against an injected
Gaussian and the residual sign against an injected 5% offset.

**Three of the nine XSL spectra are not slit-loss corrected**, which matters
here more than anywhere else, because slit loss *is* absolute flux. Only the
plain `_merged.fits` files carry `LOSS_COR = True`; the `_scl` and `_ncl_ncge`
variants do not, and the loader reads the flag per star rather than assuming:

| star | XSL file | LOSS_COR | XSL/UVES-POP |
|---|---|---|---|
| HD099648 | `_scl` | False | 0.59 |
| HD102212 | `_ncl_ncge` | False | 0.016 |
| Betelgeuse | `_ncl_ncge` | False | 0.003 |

Those are factors of **1.7, 63 and 340** — not subtle, and not a calibration
disagreement. They are drawn (shape survives normalisation) but flagged in red
on the panel and excluded from every summary number.

**For the six that are slit-loss corrected, the grey factor at the break is
0.853 to 1.279, median 0.992.** Comparable to the XSL/NGSL spread (0.90–1.04)
and, as there, it should be normalised out rather than read as astrophysical.

**But unlike XSL/NGSL, the difference is not grey.** The residual panels rise
systematically from 0 at 3500 Å, and the ratio of the two normalisation windows
(3565 Å → 3972 Å, which straddle the break) is **+6.6% median, −1.5 to +11.2%**
across the six — carried as `d_blue_red_pct`. Checked over a wider baseline with
the same smoothing, XSL/UVES-POP peaks ~+7% near 3950–4330 Å, returns to ~1.00
by 5000–6000 Å and falls to ~−6% by 7600 Å. **With n = 6 and a star-to-star
spread as large as the trend, that is a characterisation, not a correction** —
the useful conclusion is only that a single grey factor does not describe these
two libraries over the break, so normalise locally.

| star | grey (break) | rms % | grey (Hε) | rms % | blue→red % |
|---|---|---|---|---|---|
| HD016673 | 1.000 | 12.5 | 1.112 | 7.4 | +11.2 |
| HD024616 | 1.020 | 6.0 | 1.005 | 4.1 | −1.5 |
| HD084937 | 0.984 | 5.0 | 1.044 | 1.8 | +6.1 |
| HD111786 | 1.279 | 3.3 | 1.314 | 2.0 | +2.7 |
| HD140283 | 0.853 | 6.5 | 0.929 | 1.9 | +9.0 |
| HD142703 | 0.906 | 5.0 | 0.970 | 1.1 | +7.1 |

**Betelgeuse has no break panel at all.** Its UVES-POP spectrum has a 552 Å
hole at 3200.9–3753.3 Å, so there is nothing to compare across the break — see
the UVES-POP gap entry in [CAVEATS.md](CAVEATS.md), which this figure is the
third project to trip over. The script's `usable` check drops it with a printed
reason rather than drawing an empty panel.

Two further limits on reading these panels: the RVs run −171 to +101 km/s and
are removed from UVES-POP (observed frame) but not from XSL (rest frame), and
four of the nine are catalogued variables observed years apart — Betelgeuse
(SRC) and ν Vir (SRB) above all, whose panels are noise. A variable's panel is
not a calibration measurement, and its title is printed in red.

### Synthetic colours against catalogue photometry

`explore/uves_ngsl_ub.py` synthesises Johnson U-B and Stromgren c1 from both
libraries and compares them with published photometry, for the 13 overlap
stars. Both U-band indices straddle the Balmer break — bessell_U runs
3050-4150 A with lambda_eff 3571 A, Stromgren u runs 3150-3775 A — so they are
broadband measures of the break, and a library that gets the break wrong must
get them wrong by a related amount.

Writes `data/uves_ngsl_ub.csv` (the comparison),
`data/uves_ngsl_photometry.csv` (every column of both catalogues, 37 of them,
kept because the uncertainties are what say whether a residual is real) and
`figures/explore_libraries/uves_ngsl_ub.png`.

| comparison | n | median | NMAD |
|---|---|---|---|
| **c1, NGSL − catalogue** (zero point removed) | 11 | **+0.000** | **0.007** |
| U−B, NGSL − catalogue | 11 | −0.029 | 0.025 |
| c1, UVES − NGSL (identical truncated band) | 12 | −0.025 | 0.058 |
| U−B, UVES − NGSL (identical truncated band) | 12 | −0.011 | 0.024 |

**Stromgren is the better index here, and by a wide margin.** NGSL reproduces
catalogue c1 to 0.007 NMAD with no offset, against 0.025 and a −0.029 offset
for U−B. Read that as a statement about the Johnson U bandpass rather than
about NGSL: U's blue edge is set by atmospheric cutoff for ground-based work
and runs to 3050 A, while Stromgren u is narrower and better defined. The U−B
offset also does not scale with c1 over −0.13 to +0.81, so it is not the break
being mis-measured — though with n = 9 and the range carried by two variables,
that is "no evidence of a trend" rather than a constraint.

**UVES-POP cannot measure Johnson U at all.** It starts at 3200 A against U's
3050 A, so 1.96% of the band's transmission-weighted integral is missing;
`common.photometry.project` returns NaN for an uncovered filter rather than a
clipped number, so this is caught. Truncated bands are therefore applied to
BOTH libraries for the like-for-like rows above, and the cost of truncating is
measured on NGSL, which covers both: 0.008 mag in U−B, 0.005 in c1. Stromgren u
loses only 0.63%, so of the two systems it is the one UVES-POP nearly reaches.

**Errors matter more than the medians on this sample.** Catalogue `e_U-B` runs
0.005 to 0.075 mag, a factor of 15. HD076932 carries the 0.075, so its −0.116
residual — the largest of any non-variable star — is 1.5 sigma and not an
outlier at all. Its c1 agrees to 0.000 with e_c1 = 0.007. A colour residual on
these stars cannot be read without opening `data/uves_ngsl_photometry.csv`.

**Both catalogues are matched by identifier, never by position.** II/215
carries B1950 coordinates and this sample moves up to 27" between epochs;
matching II/215 positionally at 15" returned c1 for 5 of 13 and widening to
120" returned 12, by picking up neighbours. The LID prefix is NOT the same in
the two: II/215 writes HD076932 as `0100076932`, II/168 as `+100076932`.
Trying only the first form made every II/168 lookup fall through to a
positional fallback silently — it still found 11 of 13, which is why it went
unnoticed for a while.

U and B come from ONE source, Mermilliod II/168, which publishes V, B−V and
U−B together. SIMBAD is queried as a cross-check and deliberately not used: it
gives HD022484 B = 5.150 against V = 4.300, so B−V = +0.85 for an F9IV-V star
that should sit near +0.57, where II/168's homogeneous value is +0.572.

The AB→Vega direction is pinned by a selftest that runs sedpy's own Vega
spectrum through the filters, where U = B = 0 by definition. It caught a real
0.023 mag error: sedpy ships two CALSPEC files and calibrates on
`alpha_lyr_stis_005`, and a glob had picked up `_011`.

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

**The absolute level agrees; the colour does not.** XP against NGSL for
HD194453 gives a median flux ratio of **1.011** over 4000–9000 Å, which
validates the mean level and the unit conversion. It does **not** validate the
colour, and the colour is the whole dust lever — see
[the BP/RP discontinuity](#gaia-xp-has-a-discontinuity-at-the-bprp-join) below.
A median flux ratio is not a colour check.

### Gaia XP has a discontinuity at the BP/RP join

Integrating XP and NGSL through **identical** tophat bands removes resolution
from the comparison entirely (a band integral is conserved under convolution),
so what is left is pure spectrophotometry. Over the 10 sample stars with both
(`explore/xp_vs_ngsl.py` → `data/xp_ngsl_bandratio.csv`, normalised at 4293 Å):

| band λ_eff | arm | XP / NGSL | star-to-star scatter |
|---|---|---|---|
| 3459 | BP | 1.033 | 0.039 |
| 4293 | BP | 1.000 | — |
| 4864 | BP | 0.979 | 0.005 |
| 5539 | BP | 0.963 | 0.005 |
| 6097 | BP | **0.952** | 0.009 |
| 6688 | RP | **0.982** | 0.012 |
| 7489 | RP | 0.985 | 0.013 |
| 8654 | RP | 0.987 | 0.012 |

The ratio falls monotonically across BP and then **jumps by +3.1% at the
BP/RP changeover near 6400 Å**, giving a V-shaped pattern with its minimum at
the join — which is what two arms with independent flux calibrations look like
when the join is imperfect.

**It is instrumental, not astrophysical.** The common pattern reaches 4.8% while
the star-to-star scatter is 1.2%, and these stars span E(B−V) from ~0 to 0.125 —
real reddening differences would appear as scatter, not as a shape shared by
every star.

**Why it matters:** the BP-side tilt alone is equivalent to **ΔE(B−V) = 0.040
mag**, against the ~0.005 mag the Balmer-break prediction needs. Left in, it
would have been read as reddening. This is why the dust constraint comes from
NGSL alone (see [FITTING.md](FITTING.md)) and XP is not used as an independent
lever. The 1.2% reproducibility does mean the pattern could be divided out if XP
were ever needed, at the cost of XP no longer being independent of NGSL.

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

## NGSL resolution

**Moved.** The line spread function has its own document:
**[LSF.md](LSF.md)**. In one line: a Moffat, core 3.54 A (G430L) / 8.38 A
(G750L), beta = 1.52, constant in Angstroms per grating, applied with pixel
integration via `common.lsf.to_ngsl_pixels`.

Earlier resolving powers quoted for NGSL -- R = 939, R = 665 and R = 600 -- are
retracted; see [STALE.md](STALE.md) for where each came from and why it is
wrong.

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
