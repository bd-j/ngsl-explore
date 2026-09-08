# Running the model grid on Cannon

1705 nodes, each an independent ATLAS12 + SYNTHE run of ~7-15 min. Embarrassingly
parallel, so it maps directly onto a SLURM job array.

## Why the cluster

Locally the run does not scale. ATLAS12 holds a ~900 MB working set and streams
multi-GB line lists, so it is memory-bandwidth bound: nine local workers gave
**9.4 nodes/hour against 8.6 for one worker** — essentially no speedup, and
per-node wall time rose from 7 to 57 minutes. Spreading tasks across nodes,
rather than cores of one machine, is what actually helps.

## What to stage

**Code.** Compile ATLAS12 on Cannon (`cd src && make`). Do not copy the macOS
binaries — they are linked against Homebrew gfortran. Set `ATLAS12` to that
checkout; `grid/build_grid.py` reads `$ATLAS12/grids/...` and `$ATLAS12/bin/`.

**Starting atmospheres.** `$ATLAS12/grids/c3k_v2.3/at12_feh+0.00_afe+0.0/atm/`
(493 files, small). Required — every node starts from the nearest one.

**Line lists**, into `$ATLAS12/data/`:

| file | size | needed? |
|---|---|---|
| `gfallvac08oct17.dat` | 354 MB | yes |
| `gfpred29dec2014.bin` | 3.9 GB | yes |
| `hilines.bin` | 157 MB | yes, hot ion stages |
| `mol/` (from `mol.tar.gz`) | ~1 GB | yes — molecules are ON below 10000 K |
| `mol/tiototo2024.bin` | 3.9 GB | **no** |
| `mol/h2opokazatel.bin` | 392 MB | **no** |

TiO and H2O are skipped above `TEFF_COOL_LIMIT` = 8000 K and this grid starts at
8500 K, so ~4.3 GB need not be staged. Everything else in `data/` is small and
comes with the repo.

## Submitting

```bash
cd $PROJECT_DIR/grid/jobs
sbatch --array=0-99 cannon_grid.slurm            # 100 tasks, ~17 nodes each
```

Each task takes a strided slice (`--task i --ntasks N`), so tasks never collide
and the union covers the grid exactly once.

**It is idempotent.** A node is skipped only when `spec_complete()` confirms its
spectrum reaches the requested end wavelength — existence and file size would
accept a truncated file from a killed job. So resubmitting the same array picks
up whatever is missing:

```bash
python3 grid/build_grid.py --dry-run             # how many remain
sbatch --array=0-99 --export=ALL,ONLY_MISSING=1 cannon_grid.slurm
```

`ONLY_MISSING=1` slices the missing nodes rather than all nodes, so a
resubmission after a partial run spreads the remaining work evenly instead of
leaving most tasks with nothing to do.

## Checking the result

`build_grid.py` prints a tally of **every** status and the final grid count.
Check it. An earlier local run reported "384 ok, 0 failed" while producing 23%
of the grid: 1320 nodes returned `no-start-model` and were counted as neither.
The cause is fixed (nearest-neighbour starting models) and the tally now reports
all statuses, but confirm `grid now N/1705 complete` rather than trusting that
the job exited cleanly.

## Then

```bash
python3 grid/pack_grid.py     # -> models/grid.npz for the fitter
```
