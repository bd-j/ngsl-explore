"""Overlap between NGSL, MILES and the Pickles atlas.

The comparison is necessarily of two different kinds:

  NGSL vs MILES   individual stars, so a real positional cross-match
                  (done in crossmatch_libraries.py: 145/379 shared)
  Pickles vs both Pickles spectra are COMPOSITES averaged over several stars
                  per spectral type, with no published per-star input list in
                  Vizier, so overlap can only be spectral-type coverage:
                  which Pickles templates have real-star counterparts.

Writes data/library_overlap.csv
"""
import csv
import re
from collections import Counter

import astropy.units as u
import numpy as np
from astroquery.vizier import Vizier

LUM = ['I', 'II', 'III', 'IV', 'V']
CLASSES = list('OBAFGKM')

sp_re = re.compile(r'^[wr]?([OBAFGKM])\s*(\d+(?:\.\d+)?)?')
# Ia / Iab / Ib are supergiant subclasses and all map to luminosity class I;
# they must precede the bare alternatives or 'I' matches first and truncates.
lum_re = re.compile(r'(Iab|Ia|Ib|VII|VI|IV|V|III|II|I)')


def parse(sp):
    """-> (class letter, subclass, luminosity class); '' where unparseable."""
    if not sp:
        return None
    sp = sp.strip()
    m = sp_re.match(sp)
    if not m:
        return None
    cls = m.group(1)
    sub = float(m.group(2)) if m.group(2) else np.nan
    tail = sp[m.end():]
    lm = lum_re.search(tail)
    lum = lm.group(1) if lm else ''
    return cls, sub, ('I' if lum in ('Ia', 'Iab', 'Ib') else lum)


# --- the three libraries ---------------------------------------------------
ngsl = [r['sptype'] for r in csv.DictReader(open('data/ngsl_catalog.csv'))]
pick = [r['sptype'] for r in csv.DictReader(open('data/pickles_catalog.csv'))]
miles_t = Vizier(columns=['*'], row_limit=-1).get_catalogs('J/MNRAS/371/703')[0]
miles = [str(s) for s in miles_t['SpType']]

libs = {'NGSL': ngsl, 'MILES': miles, 'Pickles': pick}
parsed = {k: [p for p in (parse(s) for s in v) if p] for k, v in libs.items()}

print('=== parsed spectral types ===')
for k, v in libs.items():
    print(f'  {k:<8} {len(parsed[k])}/{len(v)} types parsed')

# --- coverage grid: class x luminosity -------------------------------------
print('\n=== Coverage: spectral class x luminosity class ===')
print('    (NGSL stars / MILES stars / Pickles templates)')
hdr = f'{"":<6}' + ''.join(f'{l:>18}' for l in LUM) + f'{"no lum class":>18}'
print(hdr)
rows = []
for cls in CLASSES:
    line = f'{cls:<6}'
    for l in LUM + ['']:
        cnt = {k: sum(1 for c, s, lm in parsed[k] if c == cls and lm == l)
               for k in libs}
        line += f'{cnt["NGSL"]:>6}/{cnt["MILES"]:>5}/{cnt["Pickles"]:>4}'
        rows.append(dict(spectral_class=cls, lum_class=l or 'none',
                         ngsl_stars=cnt['NGSL'], miles_stars=cnt['MILES'],
                         pickles_templates=cnt['Pickles'],
                         all_three='Y' if all(cnt[k] for k in libs) else 'N'))
    print(line)

with open('data/library_overlap.csv', 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0]))
    w.writeheader()
    w.writerows(rows)

both = [r for r in rows if r['all_three'] == 'Y']
print(f'\ncells populated in all three libraries: {len(both)}/{len(rows)}')
print('  ' + ', '.join(f'{r["spectral_class"]}{r["lum_class"]}' for r in both))

# --- Pickles templates near 10,000 K and their real-star counterparts ------
print('\n=== Pickles templates near 10,000 K: how many real stars back them? ===')
pk = [r for r in csv.DictReader(open('data/pickles_catalog.csv'))
      if r['teff_K'] and 8500 <= float(r['teff_K']) <= 12500]
print(f'{"template":<10}{"Teff":>7}{"D_Balmer":>10}{"NGSL":>7}{"MILES":>7}   matching NGSL stars')
ngsl_rows = list(csv.DictReader(open('data/ngsl_catalog.csv')))
for r in sorted(pk, key=lambda r: -float(r['teff_K'])):
    p = parse(r['sptype'])
    if not p:
        continue
    cls, sub, lm = p
    def match(c, s, l):
        return c == cls and (np.isnan(sub) or np.isnan(s) or abs(s - sub) <= 1) and l == lm
    n_hits = [x['target'] for x in ngsl_rows
              if (q := parse(x['sptype'])) and match(*q)]
    m_hits = sum(1 for q in parsed['MILES'] if match(*q))
    print(f'{r["sptype"]:<10}{r["teff_K"]:>7}{r["d_balmer_mag"]:>10}'
          f'{len(n_hits):>7}{m_hits:>7}   {", ".join(n_hits[:4])}')
