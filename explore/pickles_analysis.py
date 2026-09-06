"""Build a Pickles (1998) atlas catalog and measure the Balmer break in it.

Pickles spectra are COMPOSITES: each is an average of several observed stars
of the same spectral type, resampled to 5 A/pixel at R ~ 500, normalized to
1.0 at 5556 A. They are templates by spectral type, not individual stars, and
carry no error array.

Teff comes from Table 2 of the CDBS AA_README (itself Table 2 of Pickles 1998),
which covers a subset of the 131 spectra.

Writes data/pickles_catalog.csv
"""
import csv
import pathlib
import sys
import re
import numpy as np
from astropy.io import fits

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))
from common.balmer_metric import balmer_discontinuity


def parse_teff(path='data/pickles/AA_README'):
    """Rows look like: 'pickles_uk_9    A0V     9549.93'."""
    teff = {}
    row = re.compile(r'^\s*pickles_uk_(\d+)\s+(\S+)\s+([\d.]+)\s*$')
    for line in open(path):
        m = row.match(line.rstrip())
        if m:
            teff[int(m.group(1))] = (m.group(2), float(m.group(3)))
    return teff


def break_metrics(w, f):
    d, fc_blue, fc_red = balmer_discontinuity(w, f)
    return fc_blue, fc_red, d


teff_tab = parse_teff()
idx = fits.getdata('data/pickles/pickles_index.fits')

rows = []
for i, rec in enumerate(idx, start=1):
    stem = rec['FILENAME'].strip()
    sptype = rec['SPTYPE'].strip()
    num = int(stem.split('_')[-1])
    d = fits.getdata(f'data/pickles/{stem}.fits')
    w, f = d['WAVELENGTH'].astype(float), d['FLUX'].astype(float)

    fc_blue, fc_red, d_b = break_metrics(w, f)
    t_type, t_eff = teff_tab.get(num, ('', ''))

    # luminosity class and metallicity flag from the type string
    lum = re.search(r'(VII|VI|IV|V|III|II|I)$', sptype)
    rows.append({
        'num': num,
        'file': stem + '.fits',
        'sptype': sptype,
        'lum_class': lum.group(1) if lum else '',
        'metallicity': {'w': 'weak', 'r': 'rich'}.get(sptype[0], 'solar'),
        'teff_K': round(t_eff) if t_eff != '' else '',
        'teff_sptype_check': t_type,
        'wave_min': round(float(w[0])),
        'wave_max': round(float(w[-1])),
        'd_balmer_mag': round(d_b, 3) if d_b is not None else '',
    })

with open('data/pickles_catalog.csv', 'w', newline='') as fh:
    wtr = csv.DictWriter(fh, fieldnames=list(rows[0]))
    wtr.writeheader()
    wtr.writerows(rows)

print(f'{len(rows)} Pickles spectra -> data/pickles_catalog.csv')
print(f'  with Teff: {sum(1 for r in rows if r["teff_K"] != "")}')
print(f'  metallicity: ' + ', '.join(
    f'{k}={sum(1 for r in rows if r["metallicity"]==k)}'
    for k in ('solar', 'weak', 'rich')))

print('\n=== Near 10,000 K (8500-12500 K) ===')
print(f'{"file":<16}{"sptype":<10}{"lum":>5}{"Teff":>8}{"D_Balmer":>10}')
near = [r for r in rows if r['teff_K'] != '' and 8500 <= r['teff_K'] <= 12500]
for r in sorted(near, key=lambda r: -r['teff_K']):
    print(f'{r["file"]:<16}{r["sptype"]:<10}{r["lum_class"]:>5}'
          f'{r["teff_K"]:>8}{r["d_balmer_mag"]:>10}')

print('\n=== A0 luminosity sequence: break vs gravity at fixed Teff ===')
for r in sorted(rows, key=lambda r: r['num']):
    if r['sptype'].startswith('A0') and r['teff_K'] != '':
        print(f'  {r["sptype"]:<8} Teff={r["teff_K"]:>5} K   '
              f'D_Balmer={r["d_balmer_mag"]}')
