"""Build a grid of ATLAS12+SYNTHE spectra for fitting observed spectra.

Grid (aligned to the C3K v2.3 atmosphere nodes so every point starts from an
exactly-matching converged atmosphere, which converges in ~7 min instead of
diverging as a mismatched start does):

    Teff   8500-11500 K, 100 K   (31 nodes)
    log g  3.0-5.0,      0.2 dex  (11 nodes)
    [M/H]  -0.5-+0.3,    0.2 dex  ( 5 nodes)

1705 points at ~8 min each: ~227 h serial, ~25 h at 9-way parallelism.

The sampling is finer than the C3K starting grid in every dimension (C3K is
250 K / 0.5 dex, solar only), so most nodes start from an atmosphere up to
125 K and 0.25 dex away rather than an exact match. That is well inside what
ATLAS12 recovers from: the observed-star runs converged in 7 min from starts
123 K and 0.20 dex off, with final-iteration drift below 0.2 K. What it does
NOT recover from is a start thousands of K away (see CAVEATS.md).

The +0.3 [M/H] node exists so Castor ([Fe/H] = +0.19) is interpolated rather
than extrapolated.

The runner is INCREMENTAL: a point whose .spec already exists is skipped, so
extending the grid later costs only the new nodes, and an interrupted run
resumes. Use --dry-run to see what would be computed.

Outputs one .spec per node under models/grid/, then pack_grid.py collapses them
into a single interpolatable array. The per-node smoothed CSV that make_model.py
writes for single stars is suppressed here (--no-csv): it is not downsampled, so
it is larger than the .spec it derives from, and nothing reads it -- 19 GB and
2.4 hours over the full grid for no purpose.
"""
import argparse
import itertools
import os
import re
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
GRID_DIR = ROOT / 'models' / 'grid'
C3K = 'grids/c3k_v2.3/at12_feh+0.00_afe+0.0/atm'

TEFF = np.arange(8500, 11501, 100)
LOGG = np.round(np.arange(3.0, 5.001, 0.2), 2)
MH = np.round(np.arange(-0.5, 0.301, 0.2), 2)

WLBEG, WLEND = 320.0, 950.0     # nm; wide enough for both hydrogen breaks
NUMIT = 30
VTURB = 2.0


def node_name(t, g, m):
    return f't{t:05.0f}g{g:.2f}m{m:+.2f}'


def spec_complete(path):
    """True only if the .spec runs to the requested end wavelength.

    Existence and size are not enough: a run killed mid-write leaves a
    truncated file that would otherwise be skipped on resume and silently
    enter the grid as a partial spectrum.
    """
    try:
        with open(path, 'rb') as fh:
            fh.seek(0, 2)
            if fh.tell() < 1000:
                return False
            fh.seek(-300, 2)
            last = fh.read().decode('ascii', 'replace').strip().splitlines()[-1]
        return float(last.split()[0]) >= WLEND * 10.0 - 5.0
    except Exception:
        return False


_C3K_CACHE = None


def _c3k_nodes():
    """(path, Teff, logg) for every C3K starting atmosphere."""
    global _C3K_CACHE
    if _C3K_CACHE is None:
        d = Path(os.environ['ATLAS12']) / C3K
        out = []
        for p in sorted(d.glob('at12_feh+0.00_afe+0.0_t*g*.atm')):
            m = re.search(r'_t(\d{5})g(-?\d\.\d{2})$', p.stem)
            if m:
                out.append((p, float(m.group(1)), float(m.group(2))))
        _C3K_CACHE = out
    return _C3K_CACHE


def c3k_start(t, g):
    """NEAREST C3K atmosphere in both Teff and log g.

    This used to require an EXACT Teff match, which silently destroyed a run:
    C3K is spaced 250 K and this grid is 100 K, so only nodes landing on a
    shared multiple (every 500 K, since both start at 8500) found a start. 1320
    of 1705 nodes returned 'no-start-model' and were counted as neither ok nor
    failed, so a 40-hour job reported "384 ok, 0 failed" while producing 23% of
    the grid.

    Nearest-neighbour is safe here: the offsets are at most 125 K and 0.1 dex,
    and the observed-star runs converged in 7 min from starts 123 K and 0.20 dex
    away. What ATLAS12 cannot recover from is a start thousands of K off.
    """
    nodes = _c3k_nodes()
    if not nodes:
        return None
    # normalize the two axes by their grid spacings so neither dominates
    p, _, _ = min(nodes, key=lambda n: ((n[1] - t) / 250.0) ** 2
                                       + ((n[2] - g) / 0.5) ** 2)
    return p


def run_node(args):
    t, g, m = args
    name = node_name(t, g, m)
    spec = GRID_DIR / f'{name}.spec'
    if spec_complete(spec):
        return name, 'skip', 0.0
    if spec.exists():
        spec.unlink()          # truncated leftover from an interrupted run
    start = c3k_start(t, g)
    if start is None:
        return name, 'no-start-model', 0.0
    t0 = time.time()
    cmd = [sys.executable, str(ROOT / 'grid' / 'make_model.py'),
           '--star', name, '--teff', f'{t:.0f}', '--logg', f'{g:.2f}',
           '--feh', f'{m:.2f}', '--start', str(start),
           '--wlbeg', str(WLBEG), '--wlend', str(WLEND),
           '--numit', str(NUMIT), '--vturb', str(VTURB),
           '--workdir', str(GRID_DIR), '--no-csv']
    r = subprocess.run(cmd, capture_output=True, text=True)
    dt = (time.time() - t0) / 60.0
    if r.returncode != 0 or not spec.exists():
        (GRID_DIR / f'{name}.FAILED').write_text(r.stdout + '\n' + r.stderr)
        return name, 'FAILED', dt
    return name, 'ok', dt


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--workers', type=int, default=9,
                    help='parallel ATLAS12 runs (~900 MB each; default 9 of '
                         '10 cores, leaving one free)')
    ap.add_argument('--dry-run', action='store_true')
    # Cluster job arrays: each task takes a deterministic slice of the node
    # list, so tasks never collide and the set is covered exactly once.
    ap.add_argument('--task', type=int,
                    help='0-based array task index (SLURM_ARRAY_TASK_ID)')
    ap.add_argument('--ntasks', type=int,
                    help='total number of array tasks')
    ap.add_argument('--only-missing', action='store_true',
                    help='slice the MISSING nodes rather than all nodes; use '
                         'when resubmitting so tasks share the work evenly')
    a = ap.parse_args()

    if not os.environ.get('ATLAS12'):
        sys.exit('ERROR: set $ATLAS12')
    GRID_DIR.mkdir(parents=True, exist_ok=True)
    nodes = list(itertools.product(TEFF, LOGG, MH))
    todo = [n for n in nodes
            if not spec_complete(GRID_DIR / f'{node_name(*n)}.spec')]
    if a.task is not None and a.ntasks:
        pool = todo if a.only_missing else nodes
        mine = pool[a.task::a.ntasks]          # stride, so tasks interleave
        todo = [n for n in mine
                if not spec_complete(GRID_DIR / f'{node_name(*n)}.spec')]
        print(f'array task {a.task}/{a.ntasks}: {len(mine)} assigned, '
              f'{len(todo)} to compute')
    print(f'grid: {len(TEFF)} Teff x {len(LOGG)} logg x {len(MH)} [M/H] '
          f'= {len(nodes)} nodes')
    print(f'  already present : {len(nodes) - len(todo)}')
    print(f'  to compute      : {len(todo)}')
    print(f'  estimate        : {len(todo) * 8 / 60:.1f} h serial, '
          f'{len(todo) * 8 / 60 / a.workers:.1f} h at {a.workers} workers')
    print(f'  disk            : ~{len(todo) * 14.9 / 1000:.1f} GB of raw output')
    if a.dry_run or not todo:
        return

    t0 = time.time()
    tally = {}
    with ProcessPoolExecutor(max_workers=a.workers) as ex:
        futs = {ex.submit(run_node, n): n for n in todo}
        for i, f in enumerate(as_completed(futs), 1):
            name, status, dt = f.result()
            tally[status] = tally.get(status, 0) + 1
            el = (time.time() - t0) / 60
            eta = el / i * (len(todo) - i)
            print(f'[{i:3d}/{len(todo)}] {name} {status:6s} {dt:5.1f} min '
                  f'| elapsed {el:5.1f} min, eta {eta:5.1f} min', flush=True)
    print(f'\ndone in {(time.time() - t0) / 3600:.2f} h wall')
    for k in sorted(tally):
        print(f'  {k:16s} {tally[k]:5d}')
    bad = sum(v for k, v in tally.items() if k not in ('ok', 'skip'))
    produced = sum(1 for n in nodes
                   if spec_complete(GRID_DIR / f'{node_name(*n)}.spec'))
    print(f'  grid now {produced}/{len(nodes)} complete')
    if bad:
        print(f'  WARNING: {bad} node(s) neither computed nor skipped')


if __name__ == '__main__':
    main()
