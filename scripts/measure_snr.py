"""Measure S/N per star in wavelength bands, two independent ways:
  snr_stat  = median(FLUX/STATERR)         formal, from propagated counting errors
  snr_der   = DER_SNR estimator            empirical, from local pixel-to-pixel scatter
                                           (Stoehr et al. 2008, ST-ECF NL 42)
Writes data/ngsl_snr.csv
"""
import glob, csv
import numpy as np
from astropy.io import fits

BANDS = [('2000_2500', 2000, 2500), ('2500_3000', 2500, 3000),
         ('3000_3600', 3000, 3600), ('3600_4000', 3600, 4000),
         ('4000_5000', 4000, 5000), ('5000_5500', 5000, 5500),
         ('5700_7000', 5700, 7000), ('7000_9000', 7000, 9000),
         ('9000_10000', 9000, 10000)]


def der_snr(flux):
    """Median-based S/N from pixel-to-pixel scatter; insensitive to smooth features."""
    f = flux[np.isfinite(flux)]
    if len(f) < 10:
        return np.nan
    signal = np.median(f)
    noise = 0.6052697 * np.median(np.abs(2.0 * f[2:-2] - f[:-4] - f[4:]))
    return signal / noise if noise > 0 else np.nan


rows = []
for path in sorted(glob.glob('data/spectra/*.fits')):
    with fits.open(path) as h:
        name = h[0].header['TARGNAME'].strip()
        d = h[1].data
    w, fl, er = d['WAVELENGTH'], d['FLUX'], d['STATERR']
    r = {'target': name}
    for label, lo, hi in BANDS:
        m = (w >= lo) & (w < hi) & np.isfinite(fl)
        if m.sum() < 10:
            r[f'stat_{label}'] = r[f'der_{label}'] = ''
            continue
        good = m & (er > 0) & (fl > 0)
        r[f'stat_{label}'] = round(float(np.median(fl[good] / er[good])), 1) if good.sum() > 10 else ''
        v = der_snr(fl[m])
        r[f'der_{label}'] = round(float(v), 1) if np.isfinite(v) else ''
    rows.append(r)

cols = ['target'] + [f'{p}_{b[0]}' for p in ('stat', 'der') for b in BANDS]
with open('data/ngsl_snr.csv', 'w', newline='') as fh:
    w_ = csv.DictWriter(fh, fieldnames=cols)
    w_.writeheader()
    w_.writerows(rows)
print(f'{len(rows)} stars -> data/ngsl_snr.csv')

print(f'\n{"band":<12}{"median S/N (formal)":>22}{"median S/N (DER_SNR)":>23}')
for label, lo, hi in BANDS:
    s = [r[f'stat_{label}'] for r in rows if r[f'stat_{label}'] != '']
    d = [r[f'der_{label}'] for r in rows if r[f'der_{label}'] != '']
    q = lambda a: f'{np.median(a):.0f} [{np.percentile(a,16):.0f}-{np.percentile(a,84):.0f}]' if a else 'n/a'
    print(f'{label:<12}{q(s):>22}{q(d):>23}')
