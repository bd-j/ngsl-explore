"""E(B-V) for the NGSL Balmer-break sample and the UVES-POP A stars.

Three independent estimates, because none alone is sufficient:

  SFD98 / SF11   Total Galactic column through the whole dust layer, from the
                 IRSA service. These are UPPER LIMITS for our stars, which sit
                 at a few hundred pc, well inside the layer. SF11 is the
                 Schlafly & Finkbeiner (2011) recalibration of SFD98 and is
                 the one to prefer (SFD98 over-predicts by ~14%).

  photometric    E(B-V) = (B-V)_obs - (B-V)_0, with the intrinsic colour from
                 the spectral type. Direct, but only as good as the type and
                 the photometry.

  fitted         UVES-POP publishes its own E(B-V) from spectrophotometric
                 fitting; MILES publishes one for its stars. These are real
                 measurements along the actual line of sight to the star.

Writes data/reddening.csv
"""
import csv
import re
import urllib.parse
import urllib.request
from pathlib import Path

import numpy as np
from astropy.io import fits
from astroquery.simbad import Simbad

ROOT = Path(__file__).resolve().parent.parent
IRSA = 'https://irsa.ipac.caltech.edu/cgi-bin/DUST/nph-dust?locstr={ra}+{dec}+equ+j2000'

# Intrinsic (B-V)_0 for main-sequence/giant B8-A3, Pecaut & Mamajek (2013)
# and Fitzgerald (1970); adequate for a first-order reddening estimate.
BV0 = {'B8': -0.11, 'B9': -0.07, 'B9.5': -0.05, 'A0': 0.00, 'A1': 0.03,
       'A2': 0.06, 'A3': 0.09, 'A5': 0.15}


def irsa_ebv(ra, dec):
    """-> (SFD98, SF11) mean E(B-V) in the reference aperture."""
    url = IRSA.format(ra=ra, dec=dec)
    req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
    with urllib.request.urlopen(req, timeout=90) as r:
        x = r.read().decode('utf-8', 'replace')
    def grab(tag):
        m = re.search(r'<' + tag + r'>\s*([\d.]+)\s*\(mag\)\s*</' + tag + r'>', x)
        return float(m.group(1)) if m else np.nan
    return grab('refPixelValueSFD'), grab('refPixelValueSandF')


def bv0_from_type(sp):
    if not sp:
        return None
    m = re.match(r'([OBAFGKM]\d?(?:\.\d)?)', sp.strip())
    if not m:
        return None
    return BV0.get(m.group(1))


_sb = Simbad()
_sb.add_votable_fields('plx_value')


def distance_pc(name):
    """Parallax distance, to judge how much of the dust column is in front."""
    try:
        q = _sb.query_object(name)
        cn = {c.lower(): c for c in q.colnames}
        key = next((cn[k] for k in cn if 'plx' in k and 'err' not in k), None)
        plx = float(q[key][0])
        return round(1000.0 / plx, 1) if np.isfinite(plx) and plx > 0 else ''
    except Exception:
        return ''


rows = []

# --- NGSL Balmer-break sample -------------------------------------------
cat = list(csv.DictReader(open(ROOT / 'data' / 'balmer_candidates.csv')))
for r in cat:
    sfd, sf11 = irsa_ebv(r['ra_deg'], r['dec_deg'])
    bv0 = bv0_from_type(r['sptype_simbad'])
    ebv_phot = (float(r['b_minus_v']) - bv0
                if bv0 is not None and r['b_minus_v'] else np.nan)
    rows.append(dict(
        sample='NGSL', star=r['star'], sptype=r['sptype_simbad'],
        teff=r['teff_ngsl_K'], ra=r['ra_deg'], dec=r['dec_deg'],
        ebv_sfd98=round(sfd, 4), ebv_sf11=round(sf11, 4),
        ebv_phot=('' if np.isnan(ebv_phot) else round(ebv_phot, 3)),
        ebv_fitted=r['ebv_miles'] or '', fitted_src='MILES' if r['ebv_miles'] else '',
        dist_pc=distance_pc(r['star'])))

# --- UVES-POP A stars ---------------------------------------------------
for p in sorted((ROOT / 'data' / 'uves_pop').glob('*.fits.gz')):
    star = p.name.replace('.fits.gz', '')
    with fits.open(p) as f:
        t = f[2].data
        ra, dec = float(t['RA'][0]), float(t['DEC'][0])
        ebv, teff = float(t['EBV'][0]), float(t['TEFF'][0])
    sfd, sf11 = irsa_ebv(ra, dec)
    rows.append(dict(
        sample='UVES-POP', star=star, sptype='', teff=round(teff),
        ra=round(ra, 5), dec=round(dec, 5), ebv_sfd98=round(sfd, 4),
        ebv_sf11=round(sf11, 4), ebv_phot='',
        ebv_fitted=round(ebv, 4), fitted_src='UVES-POP',
        dist_pc=distance_pc(star)))

with open(ROOT / 'data' / 'reddening.csv', 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0]))
    w.writeheader()
    w.writerows(rows)

hdr = ['sample', 'star', 'teff', 'dist_pc', 'ebv_sfd98', 'ebv_sf11',
       'ebv_phot', 'ebv_fitted', 'fitted_src']
print(''.join(f'{h:>12}' for h in hdr))
for r in rows:
    print(''.join(f'{str(r[h]):>12}' for h in hdr))
print('\nSFD98/SF11 are TOTAL line-of-sight columns through the whole Galactic')
print('dust layer: UPPER LIMITS, since these stars sit inside it. Compare each')
print('against ebv_phot / ebv_fitted, and against the distance.')
