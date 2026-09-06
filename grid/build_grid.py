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
into a single interpolatable array.
"""
import argparse
import itertools
import os
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


def c3k_start(t, g):
    """Nearest C3K atmosphere: exact in Teff, nearest available log g.

    C3K is spaced 0.5 dex in log g and this grid is 0.25 dex, so half the nodes
    start 0.25 dex away. ATLAS12 recovers from that easily; what it cannot
    recover from is a start hundreds of K away in Teff (see CAVEATS.md).
    """
    a12 = Path(os.environ['ATLAS12'])
    d = a12 / C3K
    best = None
    for cand in sorted(d.glob(f'at12_feh+0.00_afe+0.0_t{t:05.0f}g*.atm')):
        gc = float(cand.stem.split('g')[-1])
        if best is None or abs(gc - g) < abs(best[1] - g):
            best = (cand, gc)
    return best[0] if best else None


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
           '--workdir', str(GRID_DIR)]
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
    a = ap.parse_args()

    if not os.environ.get('ATLAS12'):
        sys.exit('ERROR: set $ATLAS12')
    GRID_DIR.mkdir(parents=True, exist_ok=True)
    nodes = list(itertools.product(TEFF, LOGG, MH))
    todo = [n for n in nodes
            if not spec_complete(GRID_DIR / f'{node_name(*n)}.spec')]
    print(f'grid: {len(TEFF)} Teff x {len(LOGG)} logg x {len(MH)} [M/H] '
          f'= {len(nodes)} nodes')
    print(f'  already present : {len(nodes) - len(todo)}')
    print(f'  to compute      : {len(todo)}')
    print(f'  estimate        : {len(todo) * 8 / 60:.1f} h serial, '
          f'{len(todo) * 8 / 60 / a.workers:.1f} h at {a.workers} workers')
    print(f'  disk            : ~{len(todo) * 14.9 / 1000:.1f} GB of raw output')
    if a.dry_run or not todo:
        return

    done = t0 = time.time()
    ok = fail = 0
    with ProcessPoolExecutor(max_workers=a.workers) as ex:
        futs = {ex.submit(run_node, n): n for n in todo}
        for i, f in enumerate(as_completed(futs), 1):
            name, status, dt = f.result()
            if status == 'FAILED':
                fail += 1
            elif status == 'ok':
                ok += 1
            el = (time.time() - t0) / 60
            eta = el / i * (len(todo) - i)
            print(f'[{i:3d}/{len(todo)}] {name} {status:6s} {dt:5.1f} min '
                  f'| elapsed {el:5.1f} min, eta {eta:5.1f} min', flush=True)
    print(f'\ndone: {ok} ok, {fail} failed, '
          f'{(time.time() - t0) / 3600:.2f} h wall')


if __name__ == '__main__':
    main()
