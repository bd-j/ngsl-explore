"""Cross-calibrate Gaia XP against NGSL, band by band, over the whole sample.

Both are space-based and absolutely calibrated, and both are ultimately tied to
CALSPEC, so they ought to agree. They do not, in the one respect the dust
constraint depends on.

Integrating both through IDENTICAL tophat bands removes resolution from the
comparison entirely (convolution conserves a band integral), so what is left is
pure spectrophotometry. The discriminant is then simple:

  * a pattern that REPEATS star to star is instrumental,
  * a pattern that VARIES star to star can be astrophysical -- these stars span
    E(B-V) from about 0 to 0.125, so real reddening differences would show up as
    scatter, not as a common shape.

Measured over the sample the star-to-star scatter is 0.5-1.3% while the common
pattern reaches 4.8%, so the disagreement is instrumental. It is shaped like a V
with its minimum at the BP/RP join near 6400 A, which is what two arms with
independent flux calibrations look like when the join is imperfect.

Why it matters: a 5% tilt across 3459-6097 A is about 0.026 mag of E(B-V), and
the Balmer-break prediction needs E(B-V) to ~0.005. So this is 5x the error
budget, and it would have been read as reddening.

Writes data/xp_ngsl_bandratio.csv
"""
import csv
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from fitting.observations import load_xp, load_ngsl
from common.photometry import tophat, project

ROOT = Path(__file__).resolve().parent.parent
ANCHOR = 1              # band index the ratios are normalised to (4050-4550 A)
BP_RP_JOIN = 6400.0     # Gaia's BP/RP changeover


def main():
    rows = [r for r in csv.DictReader(open(ROOT / 'data' / 'sample.csv'))
            if r['tier'] != 'rejected']
    ratios, stars, lam, missing = [], [], None, []

    for r in rows:
        star = r['star']
        try:
            xp, ng = load_xp(star), load_ngsl(star)
        except (FileNotFoundError, KeyError) as exc:
            missing.append(f'{star} ({type(exc).__name__})')
            continue
        m = ng.mask
        filt = [tophat(*b) for b in xp.meta['bands']]
        nb = project(ng.wavelength[m], ng.flux[m], filt)
        ratio = xp.flux / nb
        ratio = ratio / ratio[ANCHOR]
        if lam is None:
            lam = np.array([f.wave_effective for f in filt])
        ratios.append(ratio)
        stars.append(star)

    R = np.array(ratios)
    mean, scatter = np.nanmean(R, axis=0), np.nanstd(R, axis=0)

    print('Gaia XP / NGSL through identical bands, normalised at '
          f'{lam[ANCHOR]:.0f} A\n')
    print('  ' + f'{"star":<11}' + ''.join(f'{x:>8.0f}' for x in lam))
    for s, row in zip(stars, R):
        print('  ' + f'{s:<11}' + ''.join(f'{v:>8.3f}' for v in row))
    print('  ' + '-' * (11 + 8 * len(lam)))
    print('  ' + f'{"MEAN":<11}' + ''.join(f'{v:>8.3f}' for v in mean))
    print('  ' + f'{"SCATTER":<11}' + ''.join(f'{v:>8.3f}' for v in scatter))
    if missing:
        print(f'\n  no XP: {", ".join(missing)}')

    worst = np.nanmax(np.abs(mean - 1.0))
    med_scatter = np.nanmedian(scatter[scatter > 0])
    print(f'\n  common pattern reaches {100 * worst:.1f}%; '
          f'star-to-star scatter is {100 * med_scatter:.1f}% (median)')
    print('  pattern >> scatter  =>  instrumental, not astrophysical'
          if worst > 3 * med_scatter else
          '  pattern comparable to scatter => not clearly instrumental')

    with open(ROOT / 'data' / 'xp_ngsl_bandratio.csv', 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['band_lam_eff', 'arm', 'xp_over_ngsl', 'scatter', 'n_stars'])
        for i, le in enumerate(lam):
            w.writerow([f'{le:.1f}', 'BP' if le < BP_RP_JOIN else 'RP',
                        f'{mean[i]:.5f}', f'{scatter[i]:.5f}',
                        int(np.isfinite(R[:, i]).sum())])
    print(f'\n  {len(lam)} bands, {len(stars)} stars '
          f'-> data/xp_ngsl_bandratio.csv')


if __name__ == '__main__':
    main()
