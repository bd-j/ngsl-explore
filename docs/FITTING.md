# Fitting

The question is not "what are this star's parameters". It is:

> **Can the model atmospheres reproduce the shape of the continuum around the
> Balmer break, while remaining consistent with everything else we know about
> the star?**

That phrasing decides the statistics. The break is not fitted — it is **held
out and predicted**, from parameters constrained by data that excludes it. The
deliverable is a prediction residual with an honest uncertainty, not a best-fit
χ².

## The design

| role | data | constrains |
|---|---|---|
| **condition** | NGSL collapsed into 13 synthetic bands, 3220–8180 Å | E(B−V), continuum shape |
| **condition** | XSL in named windows, continuum marginalised per segment | Teff, log g, [M/H], v sin i |
| **predict** | NGSL 3550–4000 Å (Balmer) | the answer |
| **predict** | NGSL 8180–9500 Å (Paschen) | a second, free test |
| *neither* | the hydrogen lines in NGSL | XSL resolves them ~16× better |

**Every NGSL pixel is used at most once.** That is what makes NGSL-derived
photometry legitimate: the band regions are removed from any spectral use, so
the continuum information is not double-counted against itself.
`fitting.observations.conditioning_set()` and `heldout()` return the legal sets
from one place, so the separation is enforced rather than remembered.

Why each instrument gets the job it does:

* **NGSL** is space-based spectrophotometry. It is the only thing here trusted
  for absolute calibration over a wide baseline, so it carries the continuum
  and therefore the dust.
* **XSL** resolves lines ~16× better, but it is ground-based and slit-loss
  corrected, so its continuum is not trusted. Marginalising a polynomial over
  it removes continuum shape and leaves line profiles — and **a locally
  normalised line profile cannot be changed by a smooth reddening law**. A
  degree-4 polynomial absorbs CCM89 across a 1100 Å window to 3×10⁻⁵. So XSL
  carries Teff, log g and v sin i *free of the dust degeneracy*, which is what
  makes the break prediction possible at all.

## The error budget — this is Teff-limited, not dust-limited

Sensitivity of the predicted Balmer discontinuity:

| parameter | dD/d(parameter) |
|---|---|
| Teff | **−0.0170 mag per 100 K** |
| log g | −0.0060 mag per 0.1 dex |
| [M/H] | −0.0011 mag per 0.1 dex |
| E(B−V) | +0.0027 mag per 0.01 mag |

At **fixed** Teff, dust barely matters. It matters through its *covariance* with
Teff: along the locus that keeps the conditioning set unchanged, dTeff/dE(B−V) =
**+126 K per 0.01 mag**, and that Teff swing is what moves the break.

| term | σ | error in predicted D |
|---|---|---|
| **Teff at published XSL precision** | 280 K | **0.048 mag** |
| Teff if the XSL line fit reaches | 100 K | 0.017 mag |
| log g | 0.1 dex | 0.006 mag |
| [M/H] | 0.2 dex | 0.004 mag |
| E(B−V) from the NGSL bands | ~0.005 | 0.001 mag |

The effect being chased is 0.02–0.08 mag. So the binding requirement is
**σ(Teff) ≈ 100 K, about 3× better than Arentsen+2019 achieve for these stars**,
and the XSL line-profile leg is the critical path. An earlier draft of this
document called the problem dust-limited; that was wrong, and the correction
matters because it changes which leg deserves the work.

## The synthetic bands

`fitting.observations.ngsl_band_edges()` — deterministic, from the analytic
Rydberg line positions, so it needs no model spectrum and cannot drift with the
grid. Hydrogen masked ±20 Å, both held-out windows removed, remaining stretches
split into ~400 Å bands (165 Å blueward of the break).

**Band edges need no line-free placement**, unlike the Gaia XP bands that
preceded them: NGSL is already at R = 600 and the model is broadened to R = 600,
so both sides carry the same LSF and there is no leakage mismatch to dodge.

**Bands stop just blueward of the Paschen break** (8206 Å). That keeps 89.5% of
the 3220–9480 lever arm — 0.0348 against 0.0389 mag of differential extinction
per 0.01 mag of E(B−V), because CCM89 is nearly flat redward of 8000 Å — and
buys the entire Paschen region back as a second untouched prediction.

**3220–3550 Å is the dust lever and is split in two.** It sits blueward of the
Balmer series limit, so it contains no hydrogen lines at all. It also sits ~2.8%
below the SYNTHE continuum from smooth metal blanketing — which is part of the
continuum shape under test, not a reason to exclude it. Two bands rather than
one add almost no leverage (0.0017 mag per 0.01 mag E(B−V) internally, against
0.0328 for blue-vs-8000 Å) but give a **shape check**: if they disagreed it
would point at the near-UV specifically — blanketing, CCM89's near-UV shape, or
G430L calibration at its blue edge — rather than at reddening.

## The XSL fit regions

XSL is fitted only inside named windows, not across its whole range.

**Balmer: Hα, Hβ, Hγ, Hδ, each ±50 Å with the core ±6 Å masked.** These are the
dust-immune Teff / log g diagnostic. Measured at XSL resolution, the wing merges
back into the continuum by 37–41 Å, and the 50%-depth core runs 0.8 Å (Hα) to
4.7 Å (Hδ), so ±50/±6 covers wings and excludes cores with margin. For a fast
rotator the core mask should grow — v sin i = 200 km/s adds 2.9 Å at Hγ.

The cores are masked because the observed Balmer cores carry a flux excess of
~10% of the line EW relative to these LTE models — almost certainly NLTE in
hydrogen, which the code does not treat for H. Fitting them would drag Teff and
log g to absorb physics the models are missing.

**Hε and higher orders are excluded.** They blend into one another so a local
continuum is not defined: the wing of H8 does not return to within 2% of the
continuum until **122 Å** from centre, against 37–41 Å for the four used. They
also sit inside the held-out break window.

**Metal windows chosen by measurement**, not reputation
(`explore/metal_sensitivity.py` → `data/xsl_metal_windows.csv`): two grid models
differing only in [M/H], broadened to XSL resolution, continuum-normalised, and
ranked by change in line *depth*. Depth and normalised, because the fit
marginalises XSL's continuum away — a feature that only shifts the continuum
level carries no information once that is done.

That measurement overturned two expectations:

* **Mg I b is the wrong magnesium diagnostic here.** At ~10,000 K magnesium is
  largely ionised: Mg I b 5167 gives −0.090 against Mg II 4481 at −0.121, and
  both trail the Fe II blends. Sensitivity is concentrated in **4000–4600 Å**,
  65% of the final window coverage.

**Species are derived, never assigned from memory** (`common/species.py`). Of
the top eight windows: five Fe II, one Ti II (4287.6), one Ti II + Fe II
(4535.3), one Cr II + O I (4827.3). Ranking on log gf alone would be worse than
useless here — the most numerous species over 4100–4900 Å are Co I (18,500
lines), V I (18,273) and Nb I (15,438), every one of them ionised below a
fraction of 10⁻⁴ at the line-forming temperature — so the weight is full
Saha–Boltzmann at the model atmosphere's own T and Nₑ (11,596 K, 2.7×10¹⁴ cm⁻³
at τ₅₀₀₀ = 2/3), with abundances read from its own table. A second species is
reported only within 0.3 dex, the partition-function uncertainty of the method:
Cr II leads O I at 4827 by 0.06 dex and is not separable, while Fe II leads by
1.15 dex at 4410 and is unambiguous.
* **The ranking is not stable in Teff.** At 9000 K only 11–17 of the reference
  top 30 survive, with Spearman −0.04 to 0.22; in log g it is stable (0.80–1.00).
  The sample spans 8759–10885 K, so the windows are a **union over 9000 / 10000 /
  11000 K** — 40 windows, 765 Å. A window where the line happens to be weak
  simply contributes little; omitting one loses a star's constraint outright, so
  the risk is asymmetric.

**Ca II H and K are excluded** despite K being the single most [M/H]-sensitive
feature in the optical (−0.209). They sit inside the held-out window *and* carry
an interstellar component on these sightlines. Being the strongest feature is
what made K dangerous rather than useful: an ISM line read as stellar
metallicity biases [M/H] in the same direction as the reddening, so the error
would look self-consistent. Na I D is out for the same reason.

## The linear calibration, solved in closed form

Every calibration policy is linear in its coefficients:

```
model_i = M_i * sum_k c_k B_k(lambda_i)
```

* `'scalar'` — one column. For a star this is (R/d)², plus any grey calibration
  error. NGSL bands get this.
* `'poly', n` — n+1 Chebyshev columns across the whole fitted range.
* `'segments'` — an independent Chebyshev per segment. **XSL gets this**: each
  Balmer window carries its own local continuum (order 1), while the metal
  windows are 5–35 Å wide, cannot each support one, and share a polynomial per
  arm (order 3). 16 coefficients over ~6400 pixels.
* `'none'` — no free columns.

A segment states the pixels it **applies to** separately from its polynomial
**domain**, because the metal windows are scattered across a whole arm while
sharing one polynomial, and must not also claim the Balmer pixels. Defining a
segment by its domain alone made those overlap and the design matrix went
rank-deficient; `check_segments` now raises on any overlap.

So χ² is quadratic in **c** and one weighted least-squares solve covers all of
them (`fitting/calibration.py`). Writing it once is the point: NGSL's
normalisation, XSL's continuum and the photometric scale are then provably the
same operation, so a convention error cannot apply to one dataset and not
another. Marginalising rather than profiling adds a −½ln|BᵀC⁻¹B| term and needs
a prior on **c**; that belongs with the likelihood and is not yet written.

**The held-out windows use the scalar solved on the bands**, never refit. That
is the whole point — the bands bracket the Balmer break on both sides, so the
prediction is asking whether the model's *jump* matches, with the continuum
level pinned next door.

## Instrumental broadening — constrained, not fitted

For NGSL the profile is **measured**: matching XSL to NGSL for three stars in
common gives **R = 600 ± 40**, constant in velocity, with no model involved (see
[DATA.md](DATA.md)). Use it.

Leaving `inst` free re-opens its degeneracy with Teff and log g for no gain.
This is not hypothetical: while the fitter reimplemented its own kernels,
`kind='R'` ignored the grid spacing and turned a request for R = 600 into R = 83.
Because `inst` was free, nothing crashed. The kernels now delegate to
`common.lsf`, and `Observation.resolution` is a property of the instrument.

**NGSL cannot measure v sin i below ~150 km/s.** Fractional model change after
grey rescaling, 3300–9400 Å at R = 600:

| change | rms | | change | rms |
|---|---|---|---|---|
| v sin i 0→40 | 0.014% | | R 600→560 (1σ prior) | 0.242% |
| v sin i 0→130 | 0.141% | | Teff +100 K | 0.720% |
| v sin i 0→180 | 0.266% | | log g +0.1 | 0.560% |

v sin i = 180 km/s is worth exactly as much as the 1σ width of the R prior. So
v sin i comes from XSL, where after continuum division the signal is 0.53% rms at
20 km/s, 1.33% at 40 and 2.22% at 80 — then **saturates** (2.77% at 150, 3.08%
at 250). XSL measures it over ~15–100 km/s and loses it above ~150.

Per library:

| library | `resolution` | value |
|---|---|---|
| NGSL | `('R', 600)` | measured, ±40 |
| XSL | `('R_segments', …)` | ~9800 UVB, ~11600 VIS — constant in velocity |

## Code

```
common/lines.py           hydrogen line positions, named Balmer members, ISM
                          lines. ONE definition -- hydrogen_lines previously
                          existed twice, in fitting/fit.py and again in
                          explore/plot_ngsl_vs_model.py.
fitting/observations.py   one record per dataset: data + resolution +
                          calibration + mask. conditioning_set() / heldout().
fitting/predict.py        predict(theta, observations) -> one prediction each.
                          Node-exact lookup when on grid nodes, trilinear
                          otherwise, same code path either way.
fitting/calibration.py    the linear solve.
common/photometry.py      sedpy filter projection, shared by model and data.
explore/check_predict.py  the end-to-end smoke test on one star.
explore/metal_sensitivity.py  ranks features by [M/H] sensitivity; --union
                          writes the XSL metal windows.
explore/plot_metal_lines.py   one panel per top [M/H]-sensitive feature,
                          with a band showing what the grid can reach.
```

`predict()` agrees with interpolation to 6×10⁻⁸ at a node and is 5× faster
(0.1 vs 0.5 ms), which is what a 1705-node scan per star needs. The node
tolerance is **absolute** (10⁻³ K): a step-scaled tolerance would have accepted
10241 K as the 10200 K node and returned the wrong spectrum while reporting an
exact lookup.

Still to write: `likelihood.py` (marginal, with the log-det term) and `scan.py`
(the 1705-node driver). See [../PLAN.md](../PLAN.md).

## Traps, all hit before being fixed

* **sedpy silently integrates a truncated bandpass.** `obj_counts_hires` has its
  "source does not span filter" assertion commented out, so a source truncated
  at 9500 Å returns `gaia_rp` **0.040 mag too faint with no warning**. The dust
  constraint is a colour, so `common.photometry.project` returns NaN instead.
* **Gaia broadband photometry is unusable with this grid.** `gaia_g`
  (3270–10500 Å) and `gaia_rp` (6160–10700 Å) run past the grid's 9500 Å limit;
  `gaia_bp` alone is exactly determined by its own free scalar, so it carries no
  information, and it spans the break besides.
* **Gaia XP and NGSL disagree in colour by up to 4.8%** with only 1.2%
  star-to-star scatter — instrumental, worth 0.040 mag in E(B−V), 8× the error
  budget (`explore/xp_vs_ngsl.py`). Their agreeing to a median flux ratio of
  1.011 validates the mean *level*, not the *colour*, and the colour is the
  whole lever. This is why the dust constraint comes from NGSL alone.
* **A band at exactly 3200 Å is silently dropped** for a spectrum starting at
  3201 Å, because `tophat` tapers ~10 Å past each edge. That cost the bluest and
  most important band once already; `NGSL_BAND_RANGE` is now inset.

## Superseded

The earlier design fitted the whole NGSL spectrum with emcee over free Teff,
log g, [M/H], E(B−V), v sin i, `inst`, RV and an error scale, with the break
*included* in the fit and a `DustPrior` doing the work of separating Teff from
reddening. It lived in `fitting/fit.py`, which has been **deleted** — it kept a
second copy of the forward model, which is exactly the drift that once turned a
request for R = 600 into R = 83.

Two things killed it. A Gaussian dust prior **does not bite**: with 1466 pixels
and a free error scale the likelihood formally measures E(B−V) to ~0.001, so
N(0.00, 0.02) came back 2σ out at 0.038. And fitting the break to then report
its residual answers a different question from the one at the top of this file.

UVES-POP is no longer fitted at all — it shares no star with NGSL, so it can
only be a separate sample rather than a cross-check on the same object, and its
continuum normalisation is too uncertain for a break measurement.
