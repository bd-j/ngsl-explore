"""Is each star's v sin i a measurement, or an artifact of the core mask?

The XSL Balmer windows mask the line core (+/-6 A) because the observed cores
carry ~10% more EW than these LTE models -- almost certainly NLTE in hydrogen,
which SYNTHE does not treat for H. The mask edge is therefore a boundary between
a region the model fits and one it does not, and v sin i is exactly the kind of
parameter that can be spent smearing flux across such a boundary instead of
measuring rotation.

The test: refit v sin i at 4 mask half-widths. A real rotation is a property of
the star and must not care where the mask edge is. A v sin i that tracks the
mask is measuring the mask.

This matters because at H-gamma a 200 km/s rotation is 2.9 A, against Stark
wings hundreds of A wide -- so the hydrogen wings themselves carry almost no
rotation signal, and whatever sets v sin i is the metal lines inside those
windows. For the stars clamped to the grid's [M/H] = -0.5 floor from catalog
values near -1.9, those lines are ~1.4 dex too strong in the model.

Writes data/vsini_mask_test.csv, which plot_ebv_teff.py --summary-only reads to
draw the unstable stars as open symbols.

    python3 explore/vsini_mask_test.py
"""
import csv
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from fitting.model import Grid
from fitting.observations import load_xsl
from fitting.predict import predict
from fitting.calibration import solve, chi2

ROOT = Path(__file__).resolve().parent.parent
MASKS = (6.0, 10.0, 15.0, 25.0)
VSINI = np.array([0., 10., 20., 30., 40., 60., 80., 110., 150., 200., 250., 300.])
# Two adjacent nodes of the v sin i grid. Anything inside that is the grid step,
# not a real dependence on the mask.
STABLE_SPREAD = 40.0


def best_vsini(xs, grid, T, lg, mh):
    c = [chi2(xs, solve(xs, predict(dict(teff=T, logg=lg, mh=mh, ebv=0.0,
                                         vsini=float(v)), [xs], grid)[0].value)[0])
         for v in VSINI]
    return float(VSINI[int(np.nanargmin(c))])


def main():
    grid = Grid()
    scan = ROOT / 'data' / 'ebv_teff_scan.csv'
    if not scan.exists():
        sys.exit('run plot_ebv_teff.py --all first')
    rows = list(csv.DictReader(open(scan)))
    print('v sin i vs the Balmer core mask half-width')
    print(f'{"star":<10}' + ''.join(f'{f"+/-{m:.0f}A":>9}' for m in MASKS)
          + f'{"spread":>9}   verdict')
    out = []
    for r in rows:
        T, lg, mh = (float(r['teff']), float(r['logg_node']), float(r['mh_node']))
        vs = []
        for m in MASKS:
            try:
                vs.append(best_vsini(load_xsl(r['star'], core_mask=m), grid,
                                     T, lg, mh))
            except Exception as exc:
                print(f'{r["star"]:<10}  FAILED {type(exc).__name__}: {exc}')
                vs = None
                break
        if vs is None:
            continue
        spread = max(vs) - min(vs)
        stable = spread <= STABLE_SPREAD
        print(f'{r["star"]:<10}' + ''.join(f'{v:>9.0f}' for v in vs)
              + f'{spread:>9.0f}   '
              + ('measurement' if stable else 'MASK-DRIVEN — not a measurement'))
        out.append(dict(star=r['star'],
                        **{f'vsini_mask{m:.0f}': f'{v:.0f}'
                           for m, v in zip(MASKS, vs)},
                        spread=f'{spread:.0f}',
                        vsini_stable='yes' if stable else 'no',
                        mh_clamp=f'{abs(float(r["mh_cat"]) - mh):.1f}'))

    p = ROOT / 'data' / 'vsini_mask_test.csv'
    with open(p, 'w', newline='') as fh:
        wr = csv.DictWriter(fh, fieldnames=list(out[0]))
        wr.writeheader()
        wr.writerows(out)
    bad = [o for o in out if o['vsini_stable'] == 'no']
    print(f'\n  {len(out) - len(bad)}/{len(out)} stable; mask-driven: '
          + (', '.join(o['star'] for o in bad) or 'none'))
    if bad:
        cl = [o for o in bad if float(o['mh_clamp']) > 0.5]
        print(f'  of those, {len(cl)} are clamped >0.5 dex in [M/H]'
              + (f' ({", ".join(o["star"] for o in cl)})' if cl else ''))
    print(f'  -> {p.relative_to(ROOT)}')


if __name__ == '__main__':
    main()
