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

### "R ~ 1000" is a ceiling for NGSL, not a typical value
Each grating has fixed linear dispersion, so resolving power climbs across every
segment and resets at each splice. Measured from the STIS LSF tables, actual R
runs **804 -> 1343 across G430L** and is **939 at the Balmer break** — not the
664 the 2-pixel sampling limit implies, because the CCD LSF is only
1.39-1.66 px FWHM, *narrower* than Nyquist.

**Do:** convolve models with the tabulated LSF. A 2-px Gaussian over-broadens by
~40% at 3646 A and will make a correct model look too sharp. Note also that the
data are formally undersampled, so interpolating onto a finer grid recovers
nothing real. G230LB has no published LSF — its UV resolution is unquantified.

### The smoothing kernel must be wavelength dependent
The NGSL LSF is set by a fixed dispersion per grating, so it is constant in
**Angstroms** within a grating and jumps at the splices (G430L 3.85 A ->
G750L 8.09 A, a factor >2). SYNTHE's output grid is **logarithmic**, so a fixed
sigma in pixels is a constant-**R** kernel — the wrong thing.

**Bug that bit us:** a constant-R kernel anchored at 3700 A varied from 3.33 A
at 3200 to 4.37 A at 4200 (31%) against a true LSF varying 2%. Harmless inside
one grating, off by ~3x for anyone extending below 3058 A.
See `broaden_ngsl()` in `scripts/plot_model_vs_obs.py`.

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
| NGSL v2 | **vacuum** (cross-correlation r = 0.9940 at zero shift vs 0.9894 for air) |
| UVES-POP | **air** — `CTYPE1 = AWAV` |

UVES-POP must be converted before comparison or every line sits ~1.1 A blue.
`scripts/uves_pop_load.py` converts by default. Note the NGSL determination
needed a model to settle: Balmer-centroid tests alone were inconclusive,
because the air-vac difference is 0.42 NGSL pixels.

---

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
artifact. `plot_model_vs_obs.py` now clips to the model's actual coverage and
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
- **Near-solar and ~10,000 K are nearly exclusive in UVES-POP.** Only two stars
  of 406 satisfy |[Fe/H]| <= 0.25 with 9300-11200 K and dwarf/subgiant gravity,
  and both are fast rotators.
