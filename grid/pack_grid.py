"""Pack the node spectra into one interpolatable array for the fitter.

The raw .spec files are R = 300,000 (SYNTHE's enforced floor) and 13.7 MB each;
25 GB for the full grid. Nothing downstream needs that: the observations are
R <= 18,000. Packing resamples every node onto a common log-lambda grid at
R = 50,000 and stores float32, which is ~100 MB for 1705 nodes and loads in a
second.

Flux is stored as log10(f_lambda). Interpolation is then linear in log flux,
which is much closer to correct across a 100 K Teff step than linear in flux:
the Planck function is exponential in 1/T, so log flux is nearly linear in the
grid parameters while flux is not.

Writes models/grid.npz  (teff, logg, mh, wave, logflux[nt,ng,nm,nw], filled)
"""
import argparse
import re
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
GRID_DIR = ROOT / 'models' / 'grid'
OUT = ROOT / 'models' / 'grid.npz'
PACK_R = 50000.0
WMIN, WMAX = 3200.0, 9500.0
NODE = re.compile(r'^t(\d{5})g(\d\.\d{2})m([+-]\d\.\d{2})\.spec$')
C_ANG = 2.99792458e18


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--resolution', type=float, default=PACK_R)
    ap.add_argument('--out', default=str(OUT))
    a = ap.parse_args()

    found = {}
    for p in sorted(GRID_DIR.glob('*.spec')):
        m = NODE.match(p.name)
        if m:
            found[(int(m.group(1)), float(m.group(2)), float(m.group(3)))] = p
    if not found:
        raise SystemExit('no node spectra in models/grid/')
    teff = np.array(sorted({k[0] for k in found}), float)
    logg = np.array(sorted({k[1] for k in found}), float)
    mh = np.array(sorted({k[2] for k in found}), float)
    print(f'{len(found)} nodes: {len(teff)} Teff x {len(logg)} logg x {len(mh)} [M/H] '
          f'= {len(teff)*len(logg)*len(mh)} full cube')

    n = int(np.ceil(np.log(WMAX / WMIN) * a.resolution))
    wave = WMIN * np.exp(np.arange(n) / a.resolution)
    cube = np.full((len(teff), len(logg), len(mh), n), np.nan, np.float32)
    filled = np.zeros((len(teff), len(logg), len(mh)), bool)

    for (t, g, m), p in found.items():
        w, hnu, _ = np.loadtxt(p, unpack=True)
        flam = 4.0 * np.pi * hnu * C_ANG / w ** 2
        i, j, k = (int(np.argmin(abs(teff - t))), int(np.argmin(abs(logg - g))),
                   int(np.argmin(abs(mh - m))))
        with np.errstate(divide='ignore', invalid='ignore'):
            cube[i, j, k] = np.log10(np.interp(wave, w, flam)).astype(np.float32)
        filled[i, j, k] = True

    np.savez_compressed(a.out, teff=teff, logg=logg, mh=mh, wave=wave,
                        logflux=cube, filled=filled, resolution=a.resolution)
    mb = Path(a.out).stat().st_size / 1e6
    print(f'  wavelength: {wave[0]:.0f}-{wave[-1]:.0f} A, {n} points at R={a.resolution:.0f}')
    print(f'  filled {filled.sum()}/{filled.size} cube cells')
    print(f'  -> {a.out}  ({mb:.1f} MB)')


if __name__ == '__main__':
    main()
