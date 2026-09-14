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

### 4. ~~Define the XSL fit regions~~ — DONE 2026-09-14

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

Still to do: make that scan a standing exclusion so a predicted Rydberg line
cannot enter a fit window silently, and check whether SYNTHE has an occupation
probability cutoff that should already have removed them.

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

* **Ti II 4287.6 A is deeper in HD194453 than any grid [M/H] can produce**
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
* **v sin i for the sample** — measure from XSL, where it is measurable over
  ~15–100 km/s. None of the sample has a published value in hand.
* **The grid's raw `.spec` files exist only on Cannon**; `models/grid/` holds
  35 `.atm` locally. Repacking at a different resolution, or extending in
  [M/H], needs that machine. `pack_grid.py` is also not incremental — it globs
  `models/grid/*.spec` and would re-read everything.
