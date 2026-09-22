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
`EBV_MAX = 0.10` in `explore/candidate_table.py` now enforces the cut.

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
as an upper bound only. See `explore/reddening.py`.

---

## Resolution and instrument profiles

### The NGSL line spread function lives in LSF.md

**Moved**, with the measurement rebuilt from scratch: **[LSF.md](LSF.md)**.

The trap worth carrying here is the one that produced three contradictory
numbers in this file over the project's life. **Fitting a ONE-parameter profile
to an instrument that has a core and a halo does not measure the core.** The
Gaussian inflates to split the difference, and -- because the halo's relative
weight changes across a grating -- it also DRIFTS with wavelength, which reads
as "constant in velocity, R = 600". Give the profile a tail and the core stops
moving: FWHM ~ lambda^+0.14 for a Moffat core against +0.67 for a Gaussian, on
the same data.

### The smoothing kernel must be wavelength dependent
The NGSL LSF is set by a fixed dispersion per grating, so it is constant in
**Angstroms** within a grating and jumps at the splices (G430L 3.85 A ->
G750L 8.09 A, a factor >2). SYNTHE's output grid is **logarithmic**, so a fixed
sigma in pixels is a constant-**R** kernel — the wrong thing.

**Bug that bit us:** a constant-R kernel anchored at 3700 A varied from 3.33 A
at 3200 to 4.37 A at 4200 (31%) against a true LSF varying 2%. Harmless inside
one grating, off by ~3x for anyone extending below 3058 A.
See `broaden_ngsl()` in `common/lsf.py`.

### Gaia XP has a flux discontinuity at the BP/RP join

**Symptom:** a model fitted to Gaia XP bands shows a step in the residual near
6400 A rather than the smooth tilt a reddening error would give.

**Cause:** BP and RP are separately calibrated and their join is imperfect.
Integrating XP and NGSL through identical bands over 10 stars, the ratio falls
monotonically across BP (1.033 at 3459 A to 0.952 at 6097 A) and then **jumps
+3.1% at the changeover** (0.982 at 6688 A), V-shaped with the minimum at the
join. The common pattern reaches 4.8% while star-to-star scatter is only 1.2%,
and the sample spans E(B-V) from ~0 to 0.125 -- real reddening differences would
show as scatter, not as a shape every star shares. So it is instrumental.

**Why it matters:** the BP-side tilt alone is worth **dE(B-V) = 0.040 mag**
against the ~0.005 mag the break prediction needs, so it would be read as
reddening, 8x the error budget.

**Do:** do not use Gaia XP as an independent dust lever against NGSL. The dust
constraint here comes from NGSL alone. A median flux ratio of 1.011 between the
two, which is what tempted us, validates the mean LEVEL and says nothing about
the colour -- and the colour is the entire constraint. If XP is ever needed, the
1.2% reproducibility means the pattern can be divided out
(`data/xp_ngsl_bandratio.csv`), but then XP is no longer independent of NGSL.

### UVES-POP has real coverage gaps, and they are not where you expect
**This one has cost time repeatedly — check for holes before using a star, not
after a panel comes out blank.** The delivered spectra have holes. Blanking
them with NaN is mandatory: drawing a line across a gap previously produced an
apparent flux feature at 8500 A that was very nearly investigated as physics.

**Every spectrum has holes.** Surveyed over the 22 files in `data/uves_pop/`:
13 to 15 separate holes per star, **5.9% to 40% of pixels NaN**, and every
single star has at least one hole wider than 100 A. A UVES-POP spectrum with
clean coverage does not exist in this sample.

Three kinds, all visible in the data:

| gap | cause | affects |
|---|---|---|
| 5750-5844 A | dichroic / arm split | all stars |
| 8515-8690 A, then every ~150 A redward | inter-order gaps | all stars, Paschen region |
| a whole setting, 550-930 A wide | a missing spectral setup | 3 of 22 stars — see below |

**THE BLUE IS NOT SAFE, despite what the order-width argument below predicts.**
A missing setting is not an order gap — it is tens of times the local free
spectral range — and three of the 22 stars on disk have one landing in the
Balmer region:

| star | hole | what it destroys |
|---|---|---|
| HD162678 | 3858.6-4779.4 A (921 A) | 43% of its Balmer coverage; cannot be used for the break despite being the sample's slowest rotator |
| HD138716 | 3859.2-4784.1 A (925 A) | H-epsilon entirely |
| **Betelgeuse** | **3200.9-3753.3 A (552 A)** | **the Balmer break entirely** — no data blueward of 3753 A |

HD162678 and HD138716 are missing the SAME setting; Betelgeuse is missing the
bluest one instead, plus a second at 4979-6706 A (1727 A, the widest hole in
the sample), which is how it reaches 40% NaN.

The edges line up with the instrument, not with the stars: `SETTING_JOINS` in
`explore/plot_uves_ngsl.py` measures the joins at 3733.6-3859.2 A and
4781.9-4784.1 A, and these holes begin and end exactly there. They are missing
exposures, so **which** star is affected is arbitrary and cannot be predicted
from its magnitude, colour or type — it has to be looked up per star.

**Do:** test the wavelengths you actually need, on the star you actually have.
`explore/plot_uves_xsl.py:usable` is the pattern — it checks for real (not
gap-filled) points inside the normalisation windows and drops the star from
that figure with a printed reason. Checking the array endpoints is NOT enough:
Betelgeuse runs 3200-10250 A like every other star and still has nothing across
the break. That exact mistake drew an empty panel with `grey = nan`.

Short runs (<2 A) are isolated dropouts and are interpolated across
(`fill_small_gaps`); anything wider is a real hole and must be blanked back out
AFTER smoothing, out to the kernel's reach, so no pixel built partly from
invented flux is ever drawn.

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
region has **dense order coverage** in every star while the red end is riddled
with inter-order holes. The same instrument behaves oppositely at the two ends
of its range.

**That argument covers order gaps only, and an earlier version of this entry
over-read it as "100% coverage in the blue in every star".** It is not: a
MISSING SETTING is a different failure, it owes nothing to order spacing, and
it takes out 550-930 A of the blue in three of the 22 stars (table above).

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

`explore/ngsl_wavecal_fit.py` MEASURES the correction and
`common/ngsl_wavecal.py` APPLIES it -- the same split as `explore/ngsl_lsf.py`
and `common/lsf.py`. Per (star, grating): linear for G430L (fit rms 0.04-0.46 A,
slope -0.3 to -0.8 A per 1000 A for the in-grid stars), robust constant for
G750L, which scatters window-to-window without a clean trend. G230LB cannot be
calibrated this way at all -- the models start at 3200 A -- and is left
uncorrected rather than given a fitted number.

Residual shift over the 5 G430L windows, all 13 calibrated stars:

| stage | mean | scatter | max abs |
|---|---|---|---|
| as delivered | -0.81 A | 0.54 | 1.89 |
| air->vacuum | +0.38 A | 0.50 | 1.95 |
| air->vacuum + fit | **-0.02 A** | **0.28** | **0.91** |

**THE TABLE MUST BE REFITTED WHEN THE SAMPLE CHANGES, AND A MISSING STAR USED
TO BE SILENT.** The table was fitted for the four stars of the superseded
NGSL-only sample and never regenerated for the NGSL n XSL sample, so **nine of
the thirteen stars had no row and received air->vacuum and nothing else** --
~0.8 A at the Balmer break, a quarter of a G430L pixel. Nothing failed.
`apply_wavecal` skipped the missing stars without a word, every figure still
drew, and the only symptom was an oscillating residual in the Balmer window
that looked like a modelling problem. It survived for weeks.

Two things changed as a result. `apply_wavecal` now raises a `MissingWavecal`
warning naming the star, the grating and the refit command, so the state is
visible rather than inferred. And the fit is driven from `data/sample.csv` plus
the comparison-figure list rather than from a hand-written candidate file, so
adding a star to the sample does not quietly leave it uncorrected.

What it did and did not affect, measured rather than assumed: node selection is
untouched (the NGSL leg is 13 broad tophat bands, whose chi2/N moves by <2%, and
the XSL leg never sees NGSL's wavelength solution), the adopted LSF core moved
by 0.007 A because a free shift absorbs misalignment, and the held-out break
medians moved by <0.5%. The Balmer CORE excess moved from +0.42% to +1.36%,
because a +/-4 A core mask on a 0.8 A offset samples the wings asymmetrically.

**Watch the sign.** The cross-correlation returns how far the MODEL must move
to meet the data, so the data are corrected by SUBTRACTING it. Getting this
backwards improves the break amplitude while making the line residuals worse --
a combination that should be read as a sign error, not a partial success.

## Model physics

### Balmer line cores in NGSL: the instrument, not NLTE

**This entry previously said the filled Balmer cores were "most likely NLTE in
hydrogen". That does not survive.** What forced the revision: at XSL's R ~ 9800
the same models fit the **full H-gamma profile, core included**, for the
metal-rich stars. A physical NLTE core deficit cannot be present at NGSL's
resolution and absent at R = 9800.

**A second correction, to the first version of this entry.** It claimed a winged
profile "removes 86% of the core excess". That was wrong: **the core excess is
degenerate with the effective WIDTH, not diagnostic of the SHAPE.** Changing
only the Moffat core at HD194453's ML node runs it from +2.83% at 3.54 A to
-3.83% at 7.00 A. Any profile can be tuned to zero it.

**And a third correction, about how far this entry's conclusion reaches.** The
conclusion -- "the core excess is the instrument profile, not NLTE" -- was
reached with a width that nulled the excess, so the number attached to it was
never meaningful. Measured against XSL rather than tuned, the Balmer core minus
continuum over the eight in-grid stars has a median of **+1.36%**, scattered
-1.10% to +4.22% with both signs; a bootstrap 95% interval on the median runs
-0.63% to +2.83%. Consistent with zero, so the conclusion stands; but the same
column reads -5.78% under a 7.00 A core, so no PRECISE value for a residual
excess can be quoted without naming the profile -- and that includes naming
where the profile is truncated, and which wavelength solution the data carry.
The core mask is +/-4 A, so an uncorrected ~0.8 A wavelength offset samples the
line wings asymmetrically: this column read +0.42% before the per-star wavecal
was refitted for the whole sample (see below).

The model-free half is now much better established, and is the part to rely on:
under the adopted profile, NGSL's Balmer cores are reproduced from smoothed XSL
to a median of **-0.05%** over 117 line-star combinations
(`data/ngsl_lsf_lines.csv`). Whatever is left over is not a kernel error. See
[LSF.md](LSF.md#what-the-cores-can-and-cannot-settle).

### Predicted O I Rydberg lines put spurious absorption in the models

The Kurucz list contains **predicted** (computed, not laboratory) transitions as
well as measured ones, and two of them land hard in the blue where the [M/H]
diagnostics are.

**Symptom:** the model shows a strong absorption feature the data does not have.
At 4403.4 A the model over-absorbs by **21.6% in flux**; in locally normalised
line depth it is 0.132 against an observed 0.012. A second, weaker case sits at
4827 (model 0.077, observed 0.022). For contrast, a real line in the same
spectrum behaves normally: Fe II 4550 gives model 0.085 against observed 0.108.

**Cause:** in each case a cluster of O I lines flagged `K13` whose UPPER LEVEL is
a very high Rydberg state.

| feature | transition | upper level | below the O I limit |
|---|---|---|---|
| 4403.6 / 4404.0 / 4404.7 | `(4S)3p 5P -> (4S)16s 5S` | 109334.4 cm-1 | 503 cm-1 |
| 4826.6 / 4826.8 | `(4S)3p 3P -> (4S)15d 3D` | 109348.9 cm-1 | 488 cm-1 |
| 4828.5 / 4828.6 | `(4S)3p 3P -> (4S)16s 3S` | 109341.0 cm-1 | 496 cm-1 |

Those levels cannot exist in a photosphere. At the line-forming electron density
of these models, N_e = 2.7e14 cm-3 at tau_5000 = 2/3, the Inglis-Teller estimate
puts the last surviving level at **n ~ 15**, so n = 15-16 is dissolved into the
continuum by the plasma microfield. The tabulated log gf of -0.27 to -0.83 is
also far too strong for a 3p -> 16s transition, where the oscillator strength
should fall off roughly as n^-3.

**Do:** exclude these regions from any fit window. `XSL_METAL_REJECT` in
`fitting/observations.py` records 4410.1 and 4827.3 as rejected on exactly this
ground. They are kept in the PREDICTION plots, because a feature the models get
wrong is worth looking at -- it just must not be allowed to drive a fit.

**How widespread:** scanning 3500-9500 A for lines with log gf > -2 whose upper
level is within 1500 cm-1 of the ionization limit returns 35 candidates, and the
strongest O I ones are exactly the two features already noticed. So this is a
small closed set, not a pervasive problem.

**The Si II entries that scan also returns are FALSE POSITIVES**, and the reason
matters for anyone repeating the test. Si II 5708.0 and 6701.3 (`KEP`, log gf
-0.23 and -0.25) sit only 53 cm-1 below the Si II ground-state limit, which looks
far worse than the O I cases -- but they are `3s3p4p 4D`, a **doubly excited**
configuration belonging to a series that converges on an EXCITED Si III limit,
not a high-n Rydberg state of the ground configuration. Their levels are barely
populated and they produce no measurable feature: model depth 0.004 and 0.003
against observed 0.015 and 0.009, where the real Si II 6347 gives 0.070. A
proximity-to-the-limit test alone therefore over-flags; the term designation has
to be read to tell a single-electron Rydberg series from a doubly excited one.

**Not yet checked:** whether SYNTHE has a level-dissolution or occupation
probability cutoff that should have removed the O I lines and is not being
applied.

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
shifted by 0.02-0.09 mag when fixed. See `common/balmer_metric.py`.

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
  `explore/build_catalog.py` parses them back out.
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
