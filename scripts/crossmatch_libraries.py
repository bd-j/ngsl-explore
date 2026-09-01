"""Positional cross-match of NGSL v2 against other empirical spectral libraries.

  MILES  : J/MNRAS/371/703  (Sanchez-Blazquez+ 2006), 985 stars, 3525-7500 A, FWHM 2.5 A
  MaStar : J/ApJ/883/175    (Yan+ 2019), MaNGA Stellar Library first release, 3622-10354 A

Writes data/ngsl_crossmatch.csv
"""
import csv
import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord
from astroquery.vizier import Vizier

RADIUS = 5 * u.arcsec

cat = list(csv.DictReader(open('data/ngsl_catalog.csv')))
ngsl = SkyCoord([float(r['ra']) for r in cat] * u.deg,
                [float(r['dec']) for r in cat] * u.deg)

v = Vizier(columns=['*', '_RAJ2000', '_DEJ2000'], row_limit=-1)


def fetch(catid):
    t = v.get_catalogs(catid)[0]
    print(f'  {catid}: {len(t)} rows, columns={t.colnames[:12]}')
    return t


def match(table, label):
    c = SkyCoord(table['_RAJ2000'], table['_DEJ2000'], unit=(u.deg, u.deg))
    idx, d2d, _ = ngsl.match_to_catalog_sky(c)
    hit = d2d < RADIUS
    print(f'  -> {hit.sum()}/{len(cat)} NGSL stars matched in {label}')
    return idx, d2d, hit, table


print('Fetching catalogs...')
miles = fetch('J/MNRAS/371/703')
mastar = fetch('J/ApJ/883/175')

print('\nMatching within 5"...')
mi_idx, mi_d, mi_hit, mi_t = match(miles, 'MILES')
ma_idx, ma_d, ma_hit, ma_t = match(mastar, 'MaStar')

namecol = lambda t: next((c for c in ('Name', 'Star', 'ID', 'MaStar', 'HD', 'SimbadName')
                          if c in t.colnames), t.colnames[0])
mi_n, ma_n = namecol(mi_t), namecol(ma_t)

def val(table, row_i, col, nd=2):
    """MILES value or '' when absent/masked."""
    if col not in table.colnames:
        return ''
    x = table[col][row_i]
    if np.ma.is_masked(x) or not np.isfinite(x):
        return ''
    return round(float(x), nd)


rows = []
for i, r in enumerate(cat):
    j = mi_idx[i]
    rows.append(dict(
        target=r['target'], sptype=r['sptype'], teff=r['teff'], vmag=r['vmag'],
        in_miles='Y' if mi_hit[i] else 'N',
        miles_id=str(mi_t[mi_n][j]) if mi_hit[i] else '',
        miles_sep=round(float(mi_d[i].arcsec), 2) if mi_hit[i] else '',
        miles_teff=val(mi_t, j, 'Teff', 0) if mi_hit[i] else '',
        miles_logg=val(mi_t, j, 'logg') if mi_hit[i] else '',
        miles_feh=val(mi_t, j, '[Fe/H]') if mi_hit[i] else '',
        miles_ebv=val(mi_t, j, 'E(B-V)', 3) if mi_hit[i] else '',
        in_mastar='Y' if ma_hit[i] else 'N',
        mastar_id=str(ma_t[ma_n][ma_idx[i]]) if ma_hit[i] else '',
        mastar_sep=round(float(ma_d[i].arcsec), 2) if ma_hit[i] else ''))

with open('data/ngsl_crossmatch.csv', 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0]))
    w.writeheader()
    w.writerows(rows)

both = sum(1 for r in rows if r['in_miles'] == 'Y' and r['in_mastar'] == 'Y')
print(f'\nIn both MILES and MaStar: {both}')
print('-> data/ngsl_crossmatch.csv')

print('\n=== Balmer-break candidates ===')
print(f'{"target":<10}{"MILES":>8}{"MaStar":>9}   ids')
for r in rows:
    if r['target'] in ('HD194453', 'HD040573', 'HD147550', 'HD143459'):
        print(f'{r["target"]:<10}{r["in_miles"]:>8}{r["in_mastar"]:>9}   '
              f'{r["miles_id"]} {r["mastar_id"]}')
