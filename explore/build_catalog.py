"""Merge NGSL v2 FITS headers with the readme model-fit table and the SIMBAD
magnitude table into a single catalog: data/ngsl_catalog.csv

Sources:
  docs/ngsl_delivery/aaareadme.txt   section 4 model fits (Teff, log g, log Z, alpha) - Castelli(2004)
                       models on Victoria-Regina isochrones, from NGSL v1
  docs/ngsl_delivery/Table_V2.txt    B, V, spectral type as extracted from SIMBAD
  data/spectra/*.fits  v2 primary headers (offset, data quality, coords, dates)
"""
import re, csv, glob
from astropy.io import fits

# --- readme section 4: model-fit parameters -------------------------------
# "HD017081       good      1   0.029  13057.    3.6   -0.5    n"
# "HD000886       good      3"                     (no fit performed)
fit = {}
row = re.compile(
    r'^\s*(\S+)\s+(good|suspect)\s+([123])'
    r'(?:\s+([\d.]+)\s+([\d.]+)\.?\s+(-?[\d.]+)\s+(-?[\d.]+)\s+([na]))?\s*$')
for line in open('docs/ngsl_delivery/aaareadme.txt'):
    m = row.match(line.rstrip())
    if not m:
        continue
    name, dq, fq, rms, teff, logg, logz, alpha = m.groups()
    fit[name.upper()] = dict(fit_quality=fq, fit_rms=rms or '', teff=teff or '',
                             logg=logg or '', logz=logz or '', alpha=alpha or '')

# --- Table_V2: SIMBAD magnitudes and spectral types ------------------------
# "BD+112998   16 30 16.7824 +10 59 51.741   9.70    9.07    F8"
mag = {}
trow = re.compile(r'^\s*(\S+)\s+(\d\d \d\d [\d.]+)\s+([+-]\d\d \d\d [\d.]+)'
                  r'\s+(\S+)\s+(\S+)\s*(.*?)\s*$')
for line in open('docs/ngsl_delivery/Table_V2.txt'):
    m = trow.match(line.rstrip())
    if not m:
        continue
    name, ra, dec, b, v, sp = m.groups()
    num = lambda s: '' if s in ('~', '') else s
    mag[name.upper()] = dict(bmag=num(b), vmag=num(v),
                             sptype=sp.strip() if sp.strip() != '~' else '')

# --- v2 FITS headers -------------------------------------------------------
rows = []
for path in sorted(glob.glob('data/spectra/*.fits')):
    h = fits.getheader(path)
    name = h['TARGNAME'].strip()
    key = name.upper()
    r = dict(target=name, file=path.split('/')[-1],
             ra=h['RA'], dec=h['DEC'], hipparcos=h['HIPP'],
             obsdate=h['OBSDATE'], offset_px=h['OFFSETPX'],
             dataqual=h['DATAQUAL'].strip(), minwave=h['MINWAVE'],
             maxwave=h['MAXWAVE'], maxflux=h['MAXFLUX'])
    r.update(mag.get(key, dict(bmag='', vmag='', sptype='')))
    r.update(fit.get(key, dict(fit_quality='', fit_rms='', teff='',
                               logg='', logz='', alpha='')))
    rows.append(r)

cols = ['target', 'file', 'ra', 'dec', 'hipparcos', 'obsdate', 'sptype',
        'bmag', 'vmag', 'teff', 'logg', 'logz', 'alpha', 'fit_quality',
        'fit_rms', 'offset_px', 'dataqual', 'minwave', 'maxwave', 'maxflux']
with open('data/ngsl_catalog.csv', 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=cols)
    w.writeheader()
    w.writerows(rows)

n = len(rows)
print(f'{n} stars -> data/ngsl_catalog.csv')
print(f'  readme fit table parsed : {len(fit)} entries')
print(f'  Table_V2 parsed         : {len(mag)} entries')
for f in ['teff', 'vmag', 'sptype']:
    print(f'  with {f:7s}: {sum(1 for r in rows if r[f] != "")}/{n}')
print(f'  dataqual good           : {sum(1 for r in rows if r["dataqual"]=="good")}/{n}')
