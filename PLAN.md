# Plan

Where the work stands and what comes next. See [docs/FITTING.md](docs/FITTING.md)
for why the design is what it is.

## State at 2026-09-14

The conditioning / held-out machinery works end to end on HD194453.

```bash
python3 explore/check_predict.py --star HD194453 --ebv 0.03
```

| | E(B−V) = 0 | E(B−V) = 0.03 |
|---|---|---|
| NGSL bands χ²/N | 9.23 | **0.22** |
| Balmer, **held out** | −2.68% | **+1.22%** |
| Paschen, **held out** | +6.41% | **+1.07%** |
| XSL χ²/N | 2.50 | **2.83** |

At the band-preferred reddening both held-out breaks land near +1%, from a
scalar solved only on the bands. That is the first real hint that the answer to
the project's question may be *yes* — but v sin i is pinned at 0, Teff is pinned
at the nearest node, log g and [M/H] are fixed, and nothing has been
marginalised. It is a hint, not a result.

The χ² sweep (`explore/plot_ebv_teff.py --all`) has now been run on **all 12
non-rejected stars**. It behaves as the design requires: a diagonal degeneracy
valley from the NGSL bands, a **vertical stripe** from XSL confirming its dust
immunity, and a closed intersection. Ridge slope +88 to +97 K per 0.01 mag
across the sample.

At HD194453's **self-consistent** solution (Teff = 10400 K, E(B−V) = 0.045 — the
scan fits them jointly, so they must be quoted together) both held-out breaks
are predicted to better than half a percent:

| | catalog Teff, E(B−V)=0 | scan solution |
|---|---|---|
| NGSL bands rms | 5.6% tilt | **0.72%** |
| Balmer, **held out** | −2.68% | **+0.39%** |
| Paschen, **held out** | +6.41% | **+0.15%** |

**Done:** sample (13 stars, NGSL ∩ XSL); Gaia photometry + XP; all XSL spectra
extracted; `observations.py`, `predict.py`, `calibration.py`,
`common/{lines,species,photometry,specplot}.py`; XSL fit regions with measured
metal windows and a measured velocity zero point; `fitting/fit.py` retired.

**Figures:** `predict_check_<star>.png` and `metal_lines_<star>.png` are now
produced for **every star at its node-scan ML parameters** (`--all` on both
scripts; Teff, log g, [M/H], E(B−V) and v sin i all from the same fit, via the
shared `fitting.scan.best_node`). Species labels come from the nearest grid
atmosphere rather than a bespoke `models/work/` run, so all 12 are labelled
instead of the 3 that happened to have one.

`predict_check_<star>.png` (conditioning vs held-out),
`metal_lines_<star>.png` (per-feature, species-labelled),
`ebv_teff_<star>.png` (χ² surface).

## Next

### 1. Calibrate the leg weighting properly

`fitting.scan.posterior` now combines the two legs as **χ²/dof** rather than
χ², after the raw sum was found to be choosing badly (see Done). Inverse-dof is
a stand-in for the quantity that actually matters — the EFFECTIVE number of
independent points. XSL's pixels are correlated over the LSF (~4 px) and by its
continuum, so its effective dof is far below 2342; 13 banded NGSL points with a
1% calibration floor are much closer to independent. Weighting by 1/n asserts
the legs deserve equal total say, which is defensible but not derived.

Measure it instead: the XSL residual autocorrelation length gives an effective
dof directly, and the bands' correlation comes from the NGSL calibration
covariance. Until that is done, quote results under both weightings when they
disagree.

**And the weighting is not the whole problem.** `leg_tension()` in
`explore/plot_scan.py` asks a weighting-free question — is there ANY node at
which both legs reach χ²/N < 1.5, with E(B−V) and v sin i free? **Only 5 of 12
stars have one.** For the other 7 the weighting is choosing *which* failure you
see, not removing it.

Two things follow, and both are more important than the weighting itself:

* **XSL has an irreducible χ²/N floor at every node in the grid** — 2.34 for
  HD194453, 1.76 for HD147550, 18.3 for HD164257. No choice of Teff, log g,
  [M/H], E(B−V) or v sin i gets below it. That is model inadequacy or
  underestimated XSL uncertainties, not a node-selection problem, and summing
  raw χ² lets that systematic floor dominate which node is picked. Diagnosing
  it — NLTE cores, line-list errors, XSL error bars — is the real task.
* **v sin i running to the ceiling is the clean tell.** At the χ²/dof node for
  HD106304 and HD128801, XSL's χ² falls monotonically all the way to 300 km/s
  (HD106304: 6.33 → 2.22). That is not a rotation measurement; it is broadening
  being spent to wash out model metal lines that are too strong because [M/H]
  is pinned at the grid floor. `plot_scan.py` now flags it.

### 2. Extend the grid below [M/H] = −0.5

The earlier justification here — that stars pressed to the floor were exactly
the ones whose break prediction failed — **did not survive the weighting fix**;
under χ²/dof the six stars still at −0.5 predict both held-out breaks to better
than 1%. But the floor is still doing damage, in a different place: it is those
same stars that have no node fitting both legs, and whose v sin i runs to the
300 km/s ceiling because broadening is the only way left to weaken metal lines
the grid cannot make weak enough. So the case for extending the grid stands —
the symptom was misidentified, not the cause.

What remains: [M/H] = −0.5 is a *boundary*, not a measurement, for half the
sample, so no metallicity is actually being measured for them; and Fe II EWs
say HD117880 / HD128801 / HD106304 really are metal-poor, with Mg *not*
correspondingly weak — α-enhancement, which no scaled-solar `afe+0.0` node can
represent at any [M/H]. Needs Cannon: the raw `.spec` files exist only there,
and `pack_grid.py` is not incremental.

### 3. Split the held-out residual into continuum and line cores

Unchanged and still the right next analysis step, now with more reason: the
interior stars sit at +0.95% median Balmer, and until the core excess is
separated from the continuum it is not clear how much of even that is the known
NLTE hydrogen problem rather than a continuum error.

The residual panels show that inside 3550-4000 A the residual is dominated by
the high-order **Balmer line cores**, at +5 to +10%, while the continuum either
side of the break sits near 0-2%. The same is true of the Paschen window. A
single median over the window therefore mixes two different things:

* the **continuum shape across the break** -- the actual question, and
* the **line-core excess**, already known and understood: the observed cores
  carry ~10% of the line EW in excess of these LTE models, almost certainly
  NLTE in hydrogen, which the code does not treat for H (CAVEATS.md).

Report them separately -- continuum median, core median, and the break metric
D from `common.balmer_metric`. This also means `balmer_metric`'s blue window
(3350-3630 A, degree-1 extrapolated to 3646) needs checking against the band
edges: its slope is fitted over a range that overlaps the 3385-3550 band, so the
"independent" claim needs the windows to be disjoint.

### 4. Calibrate the error model

The node scan makes this unavoidable rather than optional. Independent-pixel χ²
with XSL's ~2400 correlated pixels gives a nominal Δχ² ≤ 1 interval of a single
grid node. `fitting/likelihood.py` now carries the error-scale profile, but the
correlation length is not in it. Until it is, quote the SPREAD of the held-out
prediction over acceptable models — which the scan stores — and not the
curvature at the minimum.

## Open decisions

* ~~**[M/H] floor**~~ — **DECIDED, by running it.** "Try the current grid first"
  was the right call. Six of 12 stars sit at −0.5, but with the legs combined as
  χ²/dof they predict the held-out breaks as well as anyone else, so the floor
  is not currently costing the experiment anything. It does mean [M/H] is a
  boundary rather than a measurement for half the sample. An earlier version of
  this entry blamed the floor for a set of failures that turned out to be the
  leg weighting — see Done.
  Fe II EW measurements independently confirm HD117880 / HD128801 / HD106304
  really are metal-poor (3–4× weaker than the −0.5 model) — and that Mg is *not*
  correspondingly weak, i.e. α-enhancement, which a scaled-solar `afe+0.0` grid
  cannot represent at any [M/H].
* **Error model.** Independent-pixel χ² overstates confidence on any smooth
  parameter. Cheapest honest treatment: a ~1% systematic floor for the bands
  (already applied) plus reporting the envelope under ±1% continuum tilt
  perturbations, which is NGSL's known calibration accuracy rather than a guess.
* **Extinction law.** CCM89 at R_V = 3.1 is assumed. 3220 Å is x = 3.1 μm⁻¹,
  where CCM89, F99 and G23 diverge most; a few-percent shape difference there is
  ~0.003 mag in E(B−V). Should be a listed systematic, with an F99 refit to bound
  it.

## Things to chase

* ~~XSL sits redward of the model~~ — **DONE**. Measured at +4.21 km/s grand
  mean over 12 stars × 6 isolated lines, and applied as a zero point plus a
  per-star departure clipped at ±2 km/s (`fitting.observations.xsl_rv`). Residual
  mean +0.31 km/s. Two stars remain offset by design: HD164967 (clipped binary)
  and HD174240 (one usable line). It is NOT air/vacuum, which would be +84 km/s.
* **Apparent shifts in the BLENDED panels are a line-ratio effect, not a
  wavelength error.** Fe II 4410 shows +15.1 km/s with another Fe II 0.61 Å
  away, and 4827 shows +3.6 km/s with O I 1.30 Å away, while every isolated line
  sits at +0.7 to +4.2. When the model gets the relative strengths of a close
  pair wrong, the composite minimum moves. Do not correct these with a
  wavelength shift.

* **Ti II 4287.6 A**: at the node-scan ML parameters this is no longer out of
  reach. Held at the CATALOG node (10241 K / 3.9 / +0.0) the observed depth
  0.035 sat above the grid's whole 0.009-0.027 span; at the ML node
  (9900 K / 3.60 / -0.10) the span is 0.011-0.032 and **0 of 8 features fall
  outside the grid**, against 1 of 8 before. So part of what looked like a line
  list or abundance problem was the fixed log g. Still worth watching -- it is
  the closest to the edge of the eight -- but it is no longer evidence of
  anything on its own. Original note follows.

  (observed depth 0.035 against 0.009 at [M/H] = -0.5 and 0.027 at +0.3). The
  other seven of the top eight sit inside the grid's reach, so this is not a
  metallicity result. Now that the species is identified as **Ti II**, the
  candidates are a gf problem in the Kurucz list, or titanium not scaling with
  [M/H] in this star -- the second would be consistent with the alpha-enhancement
  already suspected from the Fe/Mg EW split. Worth checking whether the other
  Ti II windows (4535.3, and Ti II as second at 4390.9) run the same way, which
  would point at the element rather than at one line.

* **The Paschen residual was +6.4% at E(B−V) = 0** and collapses to +1.1% at
  0.03. Worth confirming it is the dust solution and not the G750L wavecal
  (constant-only, 0.6–0.9 Å rms) or a red-end model issue.
* **HD194453's preferred E(B−V) ≈ 0.03** against a photometric value of −0.01.
  Once Teff is free this will move along the degeneracy; if it does not, the
  tension is real and needs explaining.
* **v sin i for the sample** — measured from XSL, but the value depends on what
  else is held fixed, and that is now quantified. The fixed-[M/H] sweep and the
  free node scan disagree badly for exactly the stars whose other parameters
  moved:

  | star | sweep (log g, [M/H] fixed) | node scan (free) |
  |---|---|---|
  | HD074721 | 200 | **10** |
  | HD106304 | 200 | **30** |
  | HD128801 | 250 | **60** |
  | HD117880 | 110 | **40** |

  Three of those four were the ones `explore/vsini_mask_test.py` flagged as
  mask-driven, so the earlier suspicion was right in substance — v sin i was
  absorbing parameter error — even though the mask test could not confirm the
  mechanism. Quote v s in i from the node scan, not from the sweep. None of the
  sample has a published value to check against.
* **The grid's raw `.spec` files exist only on Cannon**; `models/grid/` holds
  35 `.atm` locally. Repacking at a different resolution, or extending in
  [M/H], needs that machine. `pack_grid.py` is also not incremental — it globs
  `models/grid/*.spec` and would re-read everything.

## Done

### The node scan — all 1705 nodes, 12 stars

`fitting/scan.py` → `results/<star>/scan.npz` (gitignored, ~4 MB each, ~10 min
per star). `explore/plot_scan.py` → `scan_<star>.png`, `scan_sample.png`.

**Results are now segregated by whether the grid can reach the star.** Six of
the twelve have a catalog [M/H] below the grid's −0.5 floor, so no node can
represent them and the fit pays for the mismatch somewhere else. Their figures
live in `figures/below_grid_mh/` and their numbers are never pooled with the
rest (`common/figpath.py`; the cut uses `mh_ngsl`, the column the pipeline
already uses to place a star on the grid). Until there are lower-[M/H] nodes,
**the six inside the grid are the result** and the other six are diagnostics.

**The project's question, answered for the six stars the grid can reach:** the
held-out Balmer break is predicted to a median **|0.56%|** and Paschen to
**|0.69%|**, from a scalar solved only on bands that exclude those regions —
and Teff agrees with the published values to **91 K rms**, which is the
σ(Teff) ≈ 100 K the error budget says the break prediction needs.

**A CORRECTION.** This section first reported a split — 5 stars inside the grid
predicting to ~1% and 7 pressed to the [M/H] floor failing at +10.5% — and named
the grid's metallicity floor as the cause. **That was wrong, and the cause was
the leg combination.** Adding the two legs' raw χ² let XSL's 2342 pixels
outvote the 13 bands ~180:1, so node ranking was effectively XSL-only and the
bands — the entire dust lever, and the only leg that sees the continuum — barely
moved it. The scan was selecting models that fit the lines and wrecked the
continuum. Combining as χ²/dof instead:

Split by group, since pooling them is what this section got wrong once already:

| | \multicolumn — inside the grid (6) | | below the floor (6) | |
|---|---|---|---|---|
| | raw χ² | **χ²/dof** | raw χ² | χ²/dof |
| held-out \|Balmer\| median | 0.77% | **0.56%** | 10.65% | 0.25% |
| held-out \|Paschen\| median | 1.71% | **0.69%** | 4.85% | 0.98% |
| band χ²/N median | 0.75 | **0.23** | 49.5 | 0.33 |
| Teff vs catalogs, rms | 301 K | **91 K** | 1374 K | 432 K |
| Teff vs catalogs, bias | +56 K | **−42 K** | +1296 K | +149 K |
| E(B−V) above the SF11 column | 1 | **0** | 4 | 1 |

Teff is the independent check — neither weighting is fitted to the published
values. Inside the grid the weighting takes it from 301 K to **91 K**.

Read the below-grid columns with care: χ²/dof appears to *rescue* those stars
(|Balmer| 10.65% → 0.25%), but it buys that with v sin i running to the 300 km/s
ceiling and with no node satisfying both legs. The improvement is the fit
finding somewhere else to put the error, not the model getting it right.

**A failed prediction is visible in the conditioning data.** This survived the
weighting fix and is now stated over the whole node space rather than over 12
selected nodes, which is the stronger form: per star, the correlation between
band χ²/N and |held-out Balmer| is **median r = +0.52 (log–log)**, and pooled
over all 1705 × 61 models,

* models fitting the bands at χ²/N < 1 predict Balmer to a median **2.97%**
* models at χ²/N > 10 miss it by a median **17.2%**

That is what makes this a prediction rather than a hope: whether to trust the
break is checkable from the conditioning data, before looking at it.

**XSL cannot see a wrong continuum**, and that is why the weighting mattered so
much. Its continuum is marginalised away by construction, so it constrains line
shapes and is indifferent to the continuum shape the break prediction lives on.
Under the raw sum it fitted the worst models at χ²/N = 0.6–1.1 while their break
prediction was off by 15%. A line-based Teff alone would have looked fine on
every one of them. This is the direct argument both for keeping the bands leg
and for not letting pixel count decide how much it counts.

**HD194453's dust tension was an artifact of fixing log g and [M/H].** With them
free it lands at 9900 K, log g 3.60, [M/H] −0.10, **E(B−V) = 0.005** against a
photometric −0.01. The 0.045 from the fixed-[M/H] sweep came from pinning the
other two at catalog values. Held out: Balmer +1.11%, Paschen +2.25% — against
+0.39%/+0.15% for the fixed sweep's 10400 K solution, so the two sit at
different points along the degeneracy and neither is settled.

Design notes worth keeping. The scan stores **χ², ln|A| and the pixel count, not
lnL** — the error model is still open, and a stored likelihood would freeze a
convention that a 2 h re-run would be needed to change. ln|A| had to be stored
because it depends on the model at each node and cannot be recovered afterwards.
`verify()` checks the scan's hoisted fast path against `predict()` on random
real nodes and requires **bitwise** agreement; it caught two real bugs (below).


### E(B−V) sweep and the Teff–E(B−V) χ² surface

`explore/plot_ebv_teff.py` → `figures/ebv_teff_<star>.png`. Three panels, Teff on
the grid's own nodes × E(B−V), at fixed log g and [M/H]. It behaves as the design
requires:

* **NGSL bands**: a long diagonal degeneracy valley, as expected — both
  parameters tilt the continuum. Measured ridge slope **+93 K per 0.01 mag**,
  against **+126** from the earlier synthetic test. Same sign and order, 26%
  apart, and the difference is explicable: the synthetic number came from the
  full 3300–9400 Å spectrum, this one from 13 bands stopping at 8180 Å, so they
  weight wavelengths differently. Neither is wrong; quote the banded one for the
  banded fit.
* **XSL lines**: a **vertical stripe**. Insensitive to E(B−V), exactly as the
  dust-immunity argument requires. Its E(B−V) minimum is meaningless and the
  figure says so rather than printing a number.
* **Combined**: the intersection closes, at **10300 K, E(B−V) = 0.039**.

Error inflation is floored at 1, so it may widen an error bar and never shrink
one: the bands come out at χ²/n = 0.14 and rescaling that to 1 would deflate
them by 2.7×, claiming a precision the 1% calibration floor exists to disclaim.

**Is the dust offset common to the sample? Not cleanly — so it does not indict
NGSL.** Over the 8 primary stars whose solution is not on a boundary:

    fitted - photometric E(B-V) = +0.029 +/- 0.014 (sem), scatter 0.041

A ~2σ mean offset, but the **star-to-star scatter of 0.041 is the headline
number**: it is 8× the ~0.005 mag the break prediction needs, and it is far too
large for a shared calibration error, which would show as an offset with small
scatter. The scatter, not the mean, is what has to be explained.

Two solutions are **unphysical**: HD117880 (0.142 vs an SF11 total Galactic
column of 0.077) and HD128801 (0.055 vs 0.023). SF11 integrates to infinity, so
a star inside the Galaxy cannot exceed it. **Both are [M/H]-clamped**, which
makes the grid's −0.5 floor the first thing to suspect rather than the dust.
`report()` now checks this on every run.

**v sin i is a measurement for only 8 of the 12** (`explore/vsini_mask_test.py`).
Refitting at core-mask half-widths of 6, 10, 15 and 25 Å, four stars move with
the mask — HD117880 runs 110 → 250 km/s, HD128801 and HD106304 into the 300
ceiling. A rotation cannot depend on where the mask edge is, so those four are
measuring the NLTE core mismatch. Two checks that came out negative and are
worth not repeating: dropping the three metal windows changes v sin i in *no*
star (it never goes down, so the metal windows are not driving it), and at Hγ a
200 km/s rotation is 2.9 Å against Stark wings hundreds of Å wide — so whatever
sets v sin i is the metal lines *inside* the Balmer windows, not the wings.

**A trap that cost an hour, now closed.** `check_predict.py --ebv` defaulted to
**0.0**, so the committed figure showed the null-dust case: its 13 conditioning
bands ran −5.6% to +5.0%, and the held-out break panels inherited that tilt. It
reads as a broken normalisation, but the giveaway is that the two blue sides
have **opposite signs** (−3.68% at Balmer, +5.64% at Paschen) — a mis-solved
scalar shifts every panel the same way; only a tilt can do that. The scalar is
one number and cannot absorb a slope, which is the entire point: that slope is
the dust signal. `--ebv` now defaults to the scan solution, and takes **Teff
from the scan with it** — pairing scan dust with catalog Teff double-counts the
degeneracy and moved the predicted Balmer residual from +0.4% to +2.9%.

### The XSL fit regions

Balmer windows (Hα, Hβ, Hγ, Hδ, ±50 Å, cores ±6 Å masked) plus 40 metal windows
chosen by measured [M/H] sensitivity, with a segmented calibration: order 1 per
Balmer window, order 3 per arm across the metal windows. 16 coefficients,
~6400 pixels. `explore/metal_sensitivity.py` generates the windows.

Dominant species, derived by `common/species.py` (Saha-Boltzmann weighted
against the Kurucz list, not assigned from memory):

| species | lambda (vac) | width | d(depth) | note |
|---|---|---|---|---|
| Fe II | 4410.1 | 19 | -0.246 | |
| Cr II + O I | 4827.3 | 5 | -0.196 | separated by only 0.06 dex |
| Fe II | 4550.7 | 1 | -0.158 | |
| Ti II + Fe II | 4535.3 | 1 | -0.155 | |
| Fe II | 4390.9 | 13 | -0.145 | |
| Fe II | 4134.2 | 11 | -0.136 | |
| **Ti II** | 4287.6 | 8 | -0.133 | **deeper than any grid [M/H]** |
| Fe II | 4183.0 | 20 | -0.124 | |

**Only three metal windows are fitted** (`XSL_METAL_KEEP`), with the rejections
recorded in `XSL_METAL_REJECT`:

| fitted | why |
|---|---|
| Fe II 4550.7 | clean isolated line, model matches |
| Ti II + Fe II 4535.3 | model matches |
| Fe II + Si II 4129 (4123–4138) | shifted blueward onto the Si II 4128/4131 doublet |

Rejected: **4410.1 and 4827.3 are line-list artifacts** (see below); **4287.6**
(Ti II) is deeper than any grid [M/H] and conflicts with the Fe II windows;
4390.9 and 4183.0 are simply not needed. All of them stay in the prediction
plots — a feature the models get wrong is worth looking at and must not drive
the fit.

Still open: the core mask should scale with v sin i (200 km/s adds 2.9 Å at Hγ),
and the ranking should be re-derived once [M/H] is actually fitted rather than
assumed at the grid midpoint.

### Predicted O I Rydberg lines contaminate the models

Written up in [docs/CAVEATS.md](docs/CAVEATS.md#model-physics) -- the 4403 and
4827 features are Kurucz predicted O I transitions to n = 15/16 levels that the
plasma microfield dissolves, and the Si II entries the same scan flags are false
positives (doubly excited, not Rydberg). Both windows are rejected from fitting
in `XSL_METAL_REJECT` and kept in the prediction plots.

Done: `MODEL_BAD_REGIONS` in `fitting/observations.py` now excludes both O I
regions from XSL fitting regardless of which window contains them, gated by
`drop_bad` so the prediction plots still show them. Still to do: check whether
SYNTHE has an occupation probability cutoff that should already have removed
these lines.

