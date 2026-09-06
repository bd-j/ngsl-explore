# Caveats and known issues

Traps found while building this analysis, each one discovered by being bitten by
it. Grouped by what they affect. Every entry states the symptom, the cause, and
what to do about it.

---

## Sample selection

### Chemically peculiar stars (Ap / Am)
Magnetic Ap stars have abundance patches and strong fields; Am stars have
anomalous metal lines. Neither is described by a standard scaled-solar model
atmosphere, so both are poor tests of model physics — and they are seductive,
because their sharp lines look ideal for profile work.

**Signature:** very low `v sin i` combined with strongly super-solar `[Fe/H]`.
In UVES-POP this flagged HD094660 (v sin i = 0, [Fe/H] = +0.99), HD125248
(9 km/s, +0.96), HD137509 (33 km/s, +0.97).

**Do:** screen on the abundance/rotation combination, then confirm against
SIMBAD spectral types — look for `Ap`, `Bp`, `pec`, or a trailing `m`
(`A0mA1Va` is Sirius, an Am star).

### Binaries
A composite spectrum cannot be fitted with a single-star model, and the flux
ratio varies with wavelength, which distorts continuum shape specifically.

**Caught:** HD162630 (`SB*`, dropped from this analysis entirely), Castor
(`**`, and `A1V+A2Vm`), Sirius / HD048915 (`SB*` and Am).

**Do:** check SIMBAD `otype` against `SB*`, `**`, `EB*`, `Al*`, `El*`, `RS*`.
Abundance and rotation cuts will **not** catch these — HD162630 looked
perfectly normal on both.

### Reddening
Unmodelled extinction depresses the blue side of a break and therefore
*inflates* the measured discontinuity. A reddened star looks exactly like a
model that under-predicts the break.

**Caught:** HD147550 at E(B-V) = 0.125 — 96% of the entire Galactic column on
its sightline, despite being only 141 pc away — and the only star in the NGSL
sample whose observed break exceeded the model (1.070 vs 1.003). Dropped;
`EBV_MAX = 0.10` in `scripts/candidate_table.py` now enforces the cut.

**Do NOT use SFD98 or SF11 map values as the reddening to a star.** They are
the *total* column through the whole Galactic dust layer, and these stars sit
inside it at 137-317 pc. The overestimate is severe and not uniform:

| star | SF11 (total column) | actual E(B-V) | overestimate |
|---|---|---|---|
| HD040573 | 0.470 | 0.06 (photometric) | 8x |
| HD162678 | 0.636 | 0.077 (fitted) | 8x |
| HD162817 | 0.547 | 0.102 (fitted) | 5x |

Use a fitted value where one exists (UVES-POP and MILES both publish them),
otherwise `(B-V)_obs - (B-V)_0` from the spectral type. Treat the map columns
as an upper bound only. See `scripts/reddening.py`.

---

## Resolution and instrument profiles

### NGSL's delivered resolution is R ~ 600, not what the tables imply
Three numbers have been quoted for this at different times; the last is right.

| source | claim |
|---|---|
| 2-px sampling limit | R ~ 665 at the break |
| STIS LSF tables (3.85 A) | R ~ 939 at the break |
| **measured against XSL** | **R = 600 +/- 40, constant in velocity** |

The tables describe the single-exposure optical LSF. The delivered v2 spectra
are co-adds of two dithered exposures resampled onto a common grid and are
1.7-1.8x broader, and constant in VELOCITY rather than in Angstroms (7% scatter
vs 41%). Measured with XSL as the reference, so no model is involved --
see docs/DATA.md and explore/ngsl_lsf_from_xsl.py.

**A dismissed result was right.** Fitting an effective LSF against the ATLAS12
models had preferred 7.0-7.5 A, and that was rejected as "laundering physics
into an instrumental parameter" because the Balmer cores carry genuine excess
flux the LTE models cannot produce. The reasoning was sound but the conclusion
was too strong: the model-free XSL measurement gives ~6.7 A in the same window.
Both things are true at once -- the cores are filled AND the profile is broader
than tabulated. When two explanations are both available, prefer the test that
removes one of them entirely rather than arguing about which dominates.

**The data are not undersampled.** An earlier claim that the LSF is narrower
than 2 pixels followed from the tabulated width; at the measured 6.6 A over
2.744 A pixels there are ~2.4 px per resolution element.

### The smoothing kernel must be wavelength dependent
The NGSL LSF is set by a fixed dispersion per grating, so it is constant in
**Angstroms** within a grating and jumps at the splices (G430L 3.85 A ->
G750L 8.09 A, a factor >2). SYNTHE's output grid is **logarithmic**, so a fixed
sigma in pixels is a constant-**R** kernel — the wrong thing.

**Bug that bit us:** a constant-R kernel anchored at 3700 A varied from 3.33 A
at 3200 to 4.37 A at 4200 (31%) against a true LSF varying 2%. Harmless inside
one grating, off by ~3x for anyone extending below 3058 A.
See `broaden_ngsl()` in `scripts/plot_ngsl_vs_model.py`.

### UVES-POP has real coverage gaps, and they are not where you expect
The delivered spectra have holes. Blanking them with NaN is mandatory: drawing
a line across a gap previously produced an apparent flux feature at 8500 A that
was very nearly investigated as physics.

Three kinds, all visible in the data:

| gap | cause | affects |
|---|---|---|
| 5750-5844 A | dichroic / arm split | all stars |
| 8515-8690 A, then every ~150 A redward | inter-order gaps | all stars, Paschen region |
| 3859-4779 A (921 A) | a missing spectral setup | HD162678 only |

The last is not an order gap -- it is 44x the local free spectral range -- and
it removes 43% of that star's Balmer coverage, which is why HD162678 cannot be
used for the break despite being the slowest rotator of the sample.

**Order width at the Balmer break.** The regular red-end gaps ARE the
inter-order gaps, so they measure the free spectral range directly. Order
spacing runs 140 A at 9331 A to 164 A at 10101 A, with the order number m
falling 67 -> 62, giving a grating invariant m*lambda = 633,000 +/- 28,000 A.
The gaps are ~28 A on a ~151 A spacing, so the detector captures ~81% of each
order there -- the missing 19% is what makes the red-end holes.

Carried to the Balmer break, m ~ 174 and the **free spectral range is ~21 A per
order**, about 210 pixels of the delivered 0.1 A grid and ~6x narrower than at
the red end. Since the detector covers a roughly fixed number of pixels per
order, orders **overlap comfortably in the blue**, which is why the Balmer
region has 100% coverage in every star while the red end is riddled with holes.
The same instrument behaves oppositely at the two ends of its range.

Caveats on that extrapolation: it assumes one echelle with m*lambda constant
across arms, and UVES blue and red arms share the echelle but use different
cross-dispersers, cameras and detectors -- so the 81% detector-coverage figure
certainly does NOT carry over. The invariant has 4.5% scatter because it is
derived by differencing gap centres rather than measuring order edges. Treat
m ~ 174 and FSR ~ 21 A as good to ~5%.

### UVES-POP is not delivered at R = 80,000
Native resolution is R = 80,000, but the archive product is resampled to a
**0.1 A linear grid**, so usable resolution at the Balmer break is
**R ~ 18,000** (2-px). Still ~19x NGSL, but quoting 80,000 for these files is
wrong.

### SYNTHE refuses R below 300,000
`RESOLU_MIN = 300000` in `synthe.f90` is a hard floor. `resolu` sets the
*computation* grid, not an instrumental profile; a coarse grid undersamples line
cores and then integrates them, over-absorbing the smoothed spectrum by ~10%
(17% in TiO bands). **Synthesize at >= 300,000 and convolve afterwards.**

---

## Wavelength conventions

| source | convention |
|---|---|
| ATLAS12 / SYNTHE output | **vacuum** (verified: Ca II K at 3934.773 vs vacuum 3934.777) |
| NGSL v2 | **air**, plus a per-grating residual — see below |
| UVES-POP | **air** — `CTYPE1 = AWAV` |

### NGSL is in air, and an early test here said otherwise
A first cross-correlation of NGSL against a model over **3300-4150 A only**
returned a best shift of -0.35 A and was read as "vacuum" (r = 0.9940 at zero
shift vs 0.9894 for air). **That conclusion was wrong.** Repeating the test in
seven windows from 3300 to 9100 A shows the required shift running -0.2 to
-3.2 A and tracking the air-vacuum curve; converting air->vacuum drops the mean
from -1.56 A to +0.15 A and the maximum from 3.18 A to 0.85 A.

**Lesson:** 3300-4150 A is the worst possible window for this test. The
air-vacuum offset there (~1 A) is comparable to the instrument's own zero-point
error, so the two are not separable. Test a wavelength convention over the
widest possible baseline, never in one narrow window.

### A linear-in-lambda residual survives the conversion
NGSL took no wavecals with the stellar exposures -- the readme states zero
points were derived per spectrum from stellar feature positions -- and the
gratings were reduced separately. After air->vacuum, G430L shows a clean
monotonic ramp, +1.11 A at 3500 A falling to +0.08 A at 5350 A.

This is a **slope, not a zero point**, and a constant offset is the wrong model
for it: fitting one at 4200-5600 A (where the residual is ~0) leaves ~1 A
uncorrected at the Balmer break, and makes the fit *worse* where it matters.
A linear-in-lambda residual of this kind is what a different air-vacuum
convention produces -- Edlen (1953/1966) vs Ciddor (1996), or different assumed
temperature/pressure for the air index -- so it is modelled as a line.

`scripts/ngsl_wavecal.py` fits air->vacuum plus a per-(star, grating)
correction: linear for G430L (fit rms 0.02-0.20 A, slope -0.5 to -1.2 A per
1000 A, consistent across stars), robust constant for G750L, which scatters
window-to-window without a clean trend. G230LB cannot be calibrated this way at
all -- the models start at 3200 A -- and is left uncorrected rather than given
a fitted number.

Residual shift, median over the sample:

| stage | mean | scatter | max abs |
|---|---|---|---|
| as delivered | -1.41 A | 0.96 | 2.95 |
| air->vacuum | +0.17 A | 0.57 | 1.11 |
| air->vacuum + fit | **+0.08 A** | **0.25** | **0.66** |

**Watch the sign.** The cross-correlation returns how far the MODEL must move
to meet the data, so the data are corrected by SUBTRACTING it. Getting this
backwards improves the break amplitude while making the line residuals worse --
a combination that should be read as a sign error, not a partial success.

## Model physics

### Balmer line cores are filled relative to LTE models
Across all NGSL stars the observed cores carry a **net flux excess of
~10% of the line equivalent width**. This is genuine added flux, not a
broadening mismatch: the excess is **unchanged (8.8% vs 9.0%)** whether the
model is smoothed with a 3.85 A or a 7.0 A kernel, and broadening conserves
flux by construction.

**Interpretation:** most likely **NLTE in hydrogen** — LTE Balmer cores in A
stars are known to come out too deep, and this code's `NLTE_MODE` covers only
Na I, Mg I, Ca I/II and Fe I, *not* H. Chromospheric emission was considered
and is disfavoured: these stars are at ~10,000 K, well above the granulation
boundary (~7500 K) where convective envelopes and chromospheric indicators
disappear.

**Not yet tested:** `USE_KP_HYDROGEN` in `synthe_module.f90:1045` (Kurucz-Peterson
vs the default Stehle-Hutcheon Stark profiles) is the obvious A/B.

**Do not** absorb this into a fitted instrumental LSF. Fitting a broader kernel
does minimize the residual (7.5 A, "R ~ 486"), but it launders physics into an
instrumental parameter and that number is not a resolution measurement.

### Rotation is not the explanation, but must still be applied
Rotational broadening was tested for the NGSL core excess and rejected: it needs
an implausible 300 km/s and still fits worse than a plain Gaussian, because a
rotation profile's flat-topped shape is wrong. **But** at UVES-POP resolution
rotation dominates the profile and must be applied. The library publishes
`v sin i` per star (0-230 km/s in the A stars); use it, and keep refitting open.

### Cross-code systematics
UVES-POP parameters come from fitting a **PHOENIX** grid (`GRID_NAME =
phx20atm`). Adopting their Teff / log g wholesale for an ATLAS12 comparison
imports a code-to-code systematic into the residuals.

### Starting atmospheres matter enormously
Running ATLAS12 from the shipped 5777 K solar model to a 10,241 K A star
**diverges**: `SCALE_MODEL` rescales T linearly by the Teff ratio, and the deep
layers run away (layer 80: 12,000 -> 119,000 K over 8 iterations, then 2.3e7 K),
after which the code spins in an inner loop indefinitely — it does not crash.

**Do:** start from a converged model near the target. With the C3K v2.3 grid
(within 123 K and 0.2 dex for every target) the same run converges in 7 minutes
with layer 80 stable at ~53,700 K.

---

## Measurement methodology

### The break metric is only meaningful above ~6000 K
`D = 2.5 log10(Fc_red / Fc_blue)` at 3646 A, with both continua extrapolated to
the limit so a sloping SED cancels. Below ~6000 K the continuum is strongly
*curved* and heavily line-blanketed, and D tracks SED slope rather than the
discontinuity — uncorrected, M0III returns 2.52, which would rank it the
strongest "break" in the Pickles atlas. It is an artifact.

An earlier version took the blue continuum as a median centred at ~3570 A
rather than extrapolating it to 3646 A, so the slope did not cancel; values
shifted by 0.02-0.09 mag when fixed. See `scripts/balmer_metric.py`.

### The break metric is resolution sensitive
The same model gives D = 1.072 at R = 300,000 and 1.001 at R = 939.
**Always degrade the model to the data's resolution before comparing.** Failing
to do this accounted for 0.071 of an apparent 0.089 mag model-data discrepancy.

### Normalize locally, and mask hydrogen
Flux scaling must exclude the features under test. An early version scaled in a
4000-4200 A window that *contains* H-delta at 4102.9 A. Each break panel is now
normalized in windows bracketing that break, with every Balmer and Paschen line
masked +/- 20 A.

### Never extrapolate a model past its synthesis range
`np.interp` returns edge values silently. A model synthesized over 3200-4200 A
and compared over 3200-9400 A produced a 15.8% residual RMS that was pure
artifact. `plot_ngsl_vs_model.py` now clips to the model's actual coverage and
says so in the figure title.

---

## Library-specific

- **NGSL v2 dropped the stellar parameters** that v1 carried as header keywords
  (`TEFF`, `LOG_G`, `LOG_Z`, `EBMV`, `DPC`) along with the `FLUX_UNRED` and
  `FLUX_10PC` columns. They survive only as a text table inside `aaareadme.pdf`;
  `scripts/build_catalog.py` parses them back out.
- **NGSL `STATERR` is optimistic by ~3x.** It holds propagated counting
  statistics only. Real pixel scatter in a line-free continuum gives S/N ~ 100,
  not the ~330 claimed. Inflate before any chi-squared.
- **NGSL slit offset is a systematic-error proxy.** The v2 throughput correction
  is only reliable below 0.9 px; 35 of 379 stars are flagged `suspect`.
- **Pickles spectra are composites**, averaged over several stars per spectral
  type. No error array, no measured gravity or abundance, and Teff is assigned
  by type rather than fitted. Fit their continuum shape; do not infer a single
  star's gravity or rotation from their line profiles.
- **MaStar does not overlap NGSL at all** (0/379). The cause is structural:
  MaStar targets r = 12-17 through SDSS fibres, NGSL runs V = 1.5-12.2.
- **Castor is a composite and its RMS is misleading.** It scores the worst
  Balmer RMS of the UVES-POP sample (10.4%) while its break continuum matches
  well. The RMS is inflated by hundreds of sharp metal-line residuals: at
  v sin i = 18 km/s the metal lines are fully resolved, so every abundance or
  gf error produces a tall narrow spike, and Castor B is Am with genuinely
  anomalous abundances. Fast rotators smear these away and score better for the
  wrong reason -- HD162393 at 142 km/s gets 3.3%. **Residual RMS rewards
  rotational smearing**; weight the break-region continuum instead. A smooth
  +2 to +5% tilt across 5600-8200 A and predominantly positive metal-line
  residuals are both consistent with dilution by a second component.
- **Near-solar and ~10,000 K are nearly exclusive in UVES-POP.** Only two stars
  of 406 satisfy |[Fe/H]| <= 0.25 with 9300-11200 K and dwarf/subgiant gravity,
  and both are fast rotators.
