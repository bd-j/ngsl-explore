"""Build the Balmer-break candidate table as CSV.

The five NGSL stars near 10,000 K that survive screening for data quality,
peculiarity and multiplicity. Columns combine the NGSL model fits, the
measured spectra, and independent MILES parameters where the star is in MILES.

Writes data/balmer_candidates.csv
"""
import csv
import pathlib
import re
import sys
import numpy as np
from astropy.io import fits
from astroquery.simbad import Simbad

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))
from common.balmer_metric import balmer_discontinuity

CANDS = ['HD194453', 'HD040573', 'HD147550', 'HD128801', 'HD143459']

# Reddening rejection threshold. Unmodelled extinction depresses the blue side
# of the break and inflates the measured D, so a reddened star looks like a
# model failure. HD147550 sits at E(B-V) = 0.125 -- 96% of the ENTIRE Galactic
# column on its sightline despite being only 141 pc away -- and was the sole
# star whose observed break exceeded the model. See CAVEATS.md.
EBV_MAX = 0.10

# Intrinsic (B-V)_0 by spectral type, Pecaut & Mamajek (2013) / Fitzgerald (1970)
BV0 = {'B8': -0.11, 'B9': -0.07, 'B9.5': -0.05, 'A0': 0.00, 'A1': 0.03,
       'A2': 0.06, 'A3': 0.09, 'A5': 0.15}


def ebv_photometric(sptype, b_minus_v):
    """E(B-V) = (B-V)_obs - (B-V)_0 from the spectral type."""
    if b_minus_v == '' or not sptype:
        return ''
    m = re.match(r'([OBAFGKM]\d?(?:\.\d)?)', sptype.strip())
    if not m or m.group(1) not in BV0:
        return ''
    return round(float(b_minus_v) - BV0[m.group(1)], 3)

NOTES = {
    'HD194453': 'RECOMMENDED. Best-centred star in range (offset 0.00 px); '
                'unreddened; solar-scaled. SIMBAD types A0III vs fit logg 3.9 '
                '- leave logg free.',
    'HD040573': 'Best backup. Classification and gravity agree (B9.5V, logg 4.2); '
                'D_Balmer within 0.003 mag of HD194453.',
    'HD147550': 'Sound but slightly reddened (B-V +0.05 vs -0.07 intrinsic for B9V).',
    'HD128801': 'Only candidate with independent MILES parameters, but NGSL and '
                'MILES disagree by 0.7 dex on abundance; largest slit offset.',
    'HD143459': 'EXCLUDED: horizontal-branch star, E(B-V)=0.106, and NGSL/MILES '
                'temperatures disagree by ~1000 K.',
}

cat = {r['target']: r for r in csv.DictReader(open('data/ngsl_catalog.csv'))}
xm = {r['target']: r for r in csv.DictReader(open('data/ngsl_crossmatch.csv'))}

# Current SIMBAD types, which differ from the types NGSL extracted for Table_V2.
# The disagreement matters: HD194453 is A0III in SIMBAD but the NGSL fit returns
# logg 3.9, so its luminosity class is not settled.
sb = Simbad()
sb.add_votable_fields('sp_type')
simbad_sp = {}
for name in CANDS:
    try:
        t = sb.query_object(name)
        col = next((c for c in t.colnames if c.lower() == 'sp_type'), None)
        simbad_sp[name] = str(t[col][0]).strip() if t is not None and col else ''
    except Exception:
        simbad_sp[name] = ''

rows = []
for name in CANDS:
    r, x = cat[name], xm[name]
    with fits.open('data/spectra/' + r['file']) as h:
        d = h[1].data
    w, fl, er = d['WAVELENGTH'], d['FLUX'], d['STATERR']

    d_balmer, _, _ = balmer_discontinuity(w, fl)

    # S/N at the break, both scales: STATERR is optimistic by ~3x (see report)
    b = (w > 3250) & (w < 3600)
    sn_formal = np.median(fl[b] / er[b])
    resid = fl[b] - np.polyval(np.polyfit(w[b], fl[b], 3), w[b])
    sn_real = np.median(fl[b]) / np.std(resid)

    bv = (float(r['bmag']) - float(r['vmag'])) if r['bmag'] and r['vmag'] else ''
    bv_r = round(bv, 3) if bv != '' else ''
    ebv = ebv_photometric(simbad_sp.get(name, ''), bv_r)
    reject = ('E(B-V)={:.3f} > {:.2f}'.format(ebv, EBV_MAX)
              if ebv != '' and ebv > EBV_MAX else '')
    rows.append({
        'star': name,
        'sptype_ngsl_table': r['sptype'],
        'sptype_simbad': simbad_sp.get(name, ''),
        'teff_ngsl_K': round(float(r['teff'])),
        'logg_ngsl': float(r['logg']),
        'm_h_ngsl': float(r['logz']),
        'alpha_enhanced': {'a': 'yes', 'n': 'no'}.get(r['alpha'], ''),
        'fit_quality': r['fit_quality'],
        'fit_rms': r['fit_rms'],
        'in_miles': x['in_miles'],
        'teff_miles_K': x['miles_teff'],
        'feh_miles': x['miles_feh'],
        'ebv_miles': x['miles_ebv'],
        'vmag': r['vmag'],
        'b_minus_v': bv_r,
        'ebv_phot': ebv,
        'selected': 'no' if reject else 'yes',
        'reject_reason': reject,
        'offset_px': round(float(r['offset_px']), 2),
        'dataqual': r['dataqual'],
        'd_balmer_mag': round(float(d_balmer), 3),
        'snr_3250_3600_staterr': round(float(sn_formal)),
        'snr_3250_3600_measured': round(float(sn_real)),
        'hipparcos': r['hipparcos'],
        'ra_deg': round(float(r['ra']), 6),
        'dec_deg': round(float(r['dec']), 6),
        'file': r['file'],
        'notes': NOTES[name],
    })

with open('data/balmer_candidates.csv', 'w', newline='') as fh:
    w_ = csv.DictWriter(fh, fieldnames=list(rows[0]))
    w_.writeheader()
    w_.writerows(rows)

print(f'{len(rows)} candidates -> data/balmer_candidates.csv\n')
hdr = ['star', 'sptype_simbad', 'teff_ngsl_K', 'logg_ngsl', 'm_h_ngsl',
       'ebv_phot', 'offset_px', 'd_balmer_mag', 'snr_3250_3600_measured',
       'selected']
print(''.join(f'{h:>13}' for h in hdr))
for r in rows:
    print(''.join(f'{str(r[h]):>13}' for h in hdr))
