# Plan

Where the work stands and what comes next. See [docs/FITTING.md](docs/FITTING.md)
for why the design is what it is.

## State at the end of 2026-09-12

The conditioning / held-out machinery works end to end on HD194453.

```bash
python3 explore/check_predict.py --star HD194453 --ebv 0.03 --vsini 0
```

| | E(B−V) = 0 | E(B−V) = 0.03 |
|---|---|---|
| NGSL bands χ²/N | 9.23 | **0.22** |
| Balmer, **held out** | −2.68% | **+1.22%** |
| Paschen, **held out** | +6.41% | **+1.07%** |
| XSL χ²/N | 2.50 | 2.51 |

At the band-preferred reddening both held-out breaks land near +1%, from a
scalar solved only on the bands. That is the first real hint that the answer to
the project's question may be *yes* — but v sin i is still pinned at 0, Teff is
pinned at the nearest node, and nothing has been marginalised, so it is a hint
and not a result.

Done: sample (13 stars, NGSL ∩ XSL), Gaia photometry + XP, all XSL spectra
extracted, `observations.py`, `predict.py`, `calibration.py`,
`common/photometry.py`, and `explore/check_predict.py` with residual panels
under both held-out breaks.

## Tomorrow

### 1. E(B−V) sweep, minimum χ² — the immediate task

Sweep E(B−V) at fixed node and plot χ²(E(B−V)) for the bands, with the minimum
and its curvature marked. Then the same sweep in 2-D against Teff, because the
two are covariant at +126 K per 0.01 mag and a 1-D sweep hides that.

Outputs: one panel per star, plus a summary of `E(B−V)_best` and σ from the
curvature. Watch for the minimum sitting at a *grid edge* rather than a turning
point — that is a failure, not a measurement.

### 2. More grid points — the node scan

Scale from one node to all 1705. At 0.1 ms per node lookup plus the band and
XSL projections this is minutes per star, local.

Per node, store the **whole conditional likelihood curve**, not the argmax:

```
results/<star>/scan.npz
  teff, logg, mh, ebv_grid, vsini_grid
  lnl_bands[nt,ng,nm,nE]     NGSL bands vs E(B-V)
  lnl_xsl[nt,ng,nm,nV]       XSL lines vs v sin i
  resid_balmer[nt,ng,nm,nE]  held out
  resid_paschen[nt,ng,nm,nE] held out
```

~53k floats per array — trivially small and fully inspectable. Keeping the
curves rather than point estimates is what lets the break prediction be
*marginalised* over the nuisance parameters instead of evaluated at a plug-in
value, which is the difference between an honest envelope and a misleadingly
tight curve.

The separability that makes this cheap is verified: v sin i changes broadband
band fluxes by ≤0.006 mmag, and a degree-4 polynomial absorbs CCM89 over an XSL
window to 3×10⁻⁵. So the two conditionals are ~57 evaluations per node, not 806.

### 3. Split the held-out residual into continuum and line cores

The residual panels (now in the figure) show that inside 3550-4000 A the
residual is dominated by the high-order **Balmer line cores**, at +5 to +10%,
while the continuum either side of the break sits near 0-2%. The same is true
of the Paschen window. A single median over the window therefore mixes two
different things:

* the **continuum shape across the break** -- the actual question, and
* the **line-core excess**, already known and understood: the observed cores
  carry ~10% of the line EW in excess of these LTE models, almost certainly
  NLTE in hydrogen, which the code does not treat for H (CAVEATS.md).

Report them separately -- continuum median, core median, and the break metric
D from `common.balmer_metric` -- or a genuine NLTE signature will be read as a
continuum failure. This also means `balmer_metric`'s blue window (3350-3630 A,
degree-1 extrapolated to 3646) needs checking against the band edges: its slope
is fitted over a range that overlaps the 3385-3550 band, so the "independent"
claim needs the windows to be disjoint.

### 4. Define the XSL fit regions properly

Right now XSL gets one degree-4 Chebyshev across 3501-9500 A, which is a very
long baseline for five terms and eats Teff information indiscriminately.
Replace it with explicit windows at order 1-2 each, of two kinds.

**Balmer lines.** These are the dust-immune Teff / log g diagnostic and the
reason XSL is in the analysis at all -- a locally normalised profile cannot be
changed by a smooth reddening law. XSL covers H-alpha 6563, H-beta 4861,
H-gamma 4341 and H-delta 4102 with clean separation, so each gets its own window
with a low-order local continuum. Two cautions:

* The high-order members (H-epsilon 3970 down to the series limit) are blended,
  so a "local continuum" there is ill-defined and the polynomial may absorb or
  inject break-like structure. Use them only after checking that.
* Those same members sit inside the NGSL held-out window. That is *not* double
  counting -- XSL's continuum is marginalised away, so it contributes line
  SHAPE and carries no information about the jump amplitude -- but it should be
  stated rather than assumed, since it looks like a violation at first glance.
* H-alpha's red wing approaches the telluric B band at 6860 A; keep the window
  clear of it.

**Metal lines, chosen by measured [M/H] sensitivity.** Seed measurement done:
grid models at 10200 K / log g 3.8, [M/H] = -0.5 vs +0.3, broadened to R = 9800
and continuum-normalised, ranked by change in line depth (H masked +/-25 A).
237 features exceed 0.02; the strongest:

| feature | lambda (vac) | d(depth) |
|---|---|---|
| Ca II K | 3934.8 | -0.209 |
| Fe II 4549 | 4550.0 | -0.158 |
| Mg II 4481 | 4482.4 | -0.121 |
| Si II 6347 | 6348.9 | -0.113 |
| Fe II | 5057.5 | -0.119 |
| Si II 6371 | 6373.1 | -0.078 |
| Mg I b 5167 | 5169.0 | -0.090 |

Three things to carry forward from that:

* **The sensitivity is concentrated in 3900-4600 A**, dominated by Fe II, Ti II
  and Cr II blends. Any XSL metallicity constraint will come mostly from there.
* **Mg I b is NOT a good choice at these temperatures.** It was worth trying, but
  at ~10,000 K magnesium is largely ionised: Mg I b 5167 gives -0.090 against
  Mg II 4481 at -0.121, and both are well behind the Fe II blends. Use Mg II
  4481 as the magnesium diagnostic, not Mg I b.
* **Ca II K is the single most sensitive feature and must not be used naively.**
  It carries an interstellar component on these sightlines, exactly like Na I D
  (-0.076, also excluded). Either drop both or model the ISM component; do not
  let an ISM line masquerade as stellar metallicity.

Redo the ranking at more than one node before fixing the windows -- the seed is
one Teff / log g, and the answer may move across the sample.

### 5. `likelihood.py`

`marginalize_linear(design, y, ivar)` — generalising `calibration.solve` with the
−½ln|A| term and a prior on the coefficients. Plus the analytic error-scale
profile (ŝ² = χ²/N, so node ranking is by N·ln(χ²/N)), which removes `lnerr`
from the parameter vector and calibrates the confidence scaling to the actual
residual level.

## Open decisions

* **[M/H] floor.** 8 of 13 stars fall outside the grid, almost all in [M/H]
  (grid stops at −0.5, sample reaches −1.92). Agreed: **try the current grid
  first**, and see where the scan piles up on the boundary. Fe II EW measurements
  confirm HD117880 / HD128801 / HD106304 really are metal-poor (3–4× weaker than
  the −0.5 model) — and that Mg is *not* correspondingly weak, i.e.
  α-enhancement, which a scaled-solar `afe+0.0` grid cannot represent at any
  [M/H]. Only 3 primary stars sit fully inside the grid.
* **Error model.** Independent-pixel χ² overstates confidence on any smooth
  parameter. Cheapest honest treatment: a ~1% systematic floor for the bands
  (already applied) plus reporting the envelope under ±1% continuum tilt
  perturbations, which is NGSL's known calibration accuracy rather than a guess.
* **Extinction law.** CCM89 at R_V = 3.1 is assumed. 3220 Å is x = 3.1 μm⁻¹,
  where CCM89, F99 and G23 diverge most; a few-percent shape difference there is
  ~0.003 mag in E(B−V). Should be a listed systematic, with an F99 refit to bound
  it.

## Things to chase

* **The Paschen residual was +6.4% at E(B−V) = 0** and collapses to +1.1% at
  0.03. Worth confirming it is the dust solution and not the G750L wavecal
  (constant-only, 0.6–0.9 Å rms) or a red-end model issue.
* **HD194453's preferred E(B−V) ≈ 0.03** against a photometric value of −0.01.
  Once Teff is free this will move along the degeneracy; if it does not, the
  tension is real and needs explaining.
* **v sin i for the sample** — measure from XSL, where it is measurable over
  ~15–100 km/s. None of the sample has a published value in hand.
* **`fitting/fit.py`** still holds the superseded emcee path and two live bugs
  (walkers initialised outside the prior; burn-in hardcoded at 500 rather than
  τ-based). Either fix or retire it — it must not keep a second copy of the
  forward model.
* **The grid's raw `.spec` files exist only on Cannon**; `models/grid/` holds
  35 `.atm` locally. Repacking at a different resolution, or extending in
  [M/H], needs that machine. `pack_grid.py` is also not incremental — it globs
  `models/grid/*.spec` and would re-read everything.
