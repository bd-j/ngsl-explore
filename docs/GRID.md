# Model grid

A grid of ATLAS12 + SYNTHE spectra for fitting observed spectra, so that the
stellar parameters are measured from the data rather than taken from a
catalog. The catalog values may be wrong, biased, or derived by a method with
its own systematics relative to these spectra — the point of fitting is not to
assume otherwise.

## Extent and spacing

| dimension | range | spacing | nodes |
|---|---|---|---|
| Teff | 8500 – 11500 K | 100 K | 31 |
| log g | 3.0 – 5.0 | 0.2 dex | 11 |
| [M/H] | −0.5 – +0.3 | 0.2 dex | 5 |

**1705 nodes**, ~8 min each: ~227 h serial, **~25 h at 9 workers**, ~25 GB of
raw output.

Sampling is finer than the observational precision on purpose, so that
interpolation error is negligible next to measurement error. log g is the axis
worth resolving here: the Balmer line wings are the gravity diagnostic, and
gravity is the parameter most likely to carry a catalog systematic. The +0.3
[M/H] node exists so Castor ([Fe/H] = +0.19) is interpolated, not extrapolated.

## Starting atmospheres, and why they matter

Each node is converged with the **current** ATLAS12 build, starting from the
nearest C3K v2.3 atmosphere (exact in Teff, nearest in log g). C3K is only a
starting guess — never the output. Regenerating is not cosmetic: against the
old-version C3K file for t10250g4.00, the freshly converged structure differs
by up to **97 K** (median 36 K).

The grid is finer than C3K (which is 250 K / 0.5 dex, solar only), so most
nodes start up to 125 K and 0.25 dex away. That is well inside ATLAS12's
recovery range — the observed-star runs converged in 7 min from starts 123 K
and 0.20 dex off, with final-iteration drift below 0.2 K.

What ATLAS12 does **not** recover from is a start thousands of K away. From the
shipped 5777 K solar model, a 10241 K run diverges: the deep layers run away
(layer 80: 12,000 → 119,000 K over 8 iterations, then 2.3e7 K) and the code
then spins in an inner loop rather than crashing. See CAVEATS.md.

`numit = 30` is enough: final-iteration drift is ≤ 0.2 K at every depth.

## Running it

```bash
export ATLAS12=/path/to/atlas12
python3 grid/build_grid.py --dry-run          # what would be computed
python3 grid/build_grid.py --workers 9        # run it
```

**It is resumable and incremental.** A node is skipped only if its `.spec` is
verified complete — `spec_complete()` checks that the spectrum reaches the
requested end wavelength, because existence and file size alone would accept a
truncated file from an interrupted run and admit a partial spectrum into the
grid. Anything short is deleted and recomputed. So interrupting is safe: rerun
the same command. Widening the grid later costs only the new nodes.

Memory is ~900 MB per ATLAS12 process, so worker count is bounded by RAM rather
than cores; 9 workers on a 10-core machine leaves one free.

## Output

`models/grid/<node>.{atm,spec,iter,flux}`, node named `t10000g4.00m-0.10`.

| file | each | ×1705 |
|---|---|---|
| `.spec` synthesis, R = 300,000, 3200–9500 Å | 13.7 MB | 23.4 GB |
| `.flux` emergent flux | 0.81 MB | 1.4 GB |
| `.iter` convergence diagnostics | 0.41 MB | 0.7 GB |
| `.atm` converged atmosphere | 12 KB | 20 MB |

Synthesis runs at R = 300,000 because that is SYNTHE's enforced floor
(`RESOLU_MIN`): `resolu` sets the computation grid, not an instrument profile,
and a coarse grid undersamples line cores and over-absorbs the smoothed
spectrum by ~10%. Nothing downstream needs that resolution — the observations
are R ≤ 18,000 — so the fitter loads a packed array resampled to R ≈ 50,000
(~100 MB for the whole grid) rather than the raw files.

All of `models/` is gitignored: 25 GB, and fully regenerable from this script.
