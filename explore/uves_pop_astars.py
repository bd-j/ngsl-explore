"""Select A stars from the re-calibrated UVES-POP library (VOXAstro-SL).

Library: 406 stars, R = 80,000, from 3200 A, absolute flux calibration good to
1.5-4%, E(B-V) estimated per star. Re-reduction of Bagnulo et al. (2003) UVES
data described in 2023ApJS..266...11B.

Why it matters here: at R=80,000 the Balmer line cores are resolved ~85x better
than NGSL (R~950), so the core-filling excess measured against the ATLAS12
models can be seen as a profile rather than inferred from integrated flux.
The catalog also carries v sin i, which the model comparison needs.

Writes data/uves_pop_astars.csv (A stars) and data/uves_pop_all.csv (all 406).
"""
import csv
import json
import urllib.request
from pathlib import Path

API = 'https://sl.voxastro.org/api/objects/UVES-POP/'
SPEC = 'http://gal-02.sai.msu.ru/uves-pop/model_spec/fit_res/v221115/{name}.fits.gz'
ROOT = Path(__file__).resolve().parent.parent

# A-type main sequence spans roughly 7500-10000 K; widen slightly to catch
# late B and early F at the boundaries, since the library's Teff are fitted.
A_RANGE = (7000.0, 11500.0)
# The NGSL Balmer-break sample, for a direct-overlap check.
NGSL_SAMPLE = ['HD194453', 'HD040573', 'HD147550', 'HD128801', 'HD143459']

COLS = ['name', 'teff', 'e_teff', 'logg', 'e_logg', 'fe_h', 'e_fe_h', 'a_fe',
        'vsini', 'e_vsini', 'v', 'eso_starid', 'spec_url']


def fetch():
    req = urllib.request.Request(API, headers={'User-Agent': 'Mozilla/5.0'})
    with urllib.request.urlopen(req, timeout=90) as r:
        return json.load(r)


def row(o):
    d = {c: o.get(c, '') for c in COLS if c != 'spec_url'}
    d['spec_url'] = SPEC.format(name=o['name'])
    return d


def main():
    objs = fetch()
    print(f'UVES-POP: {len(objs)} stars')

    with open(ROOT / 'data' / 'uves_pop_all.csv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=COLS)
        w.writeheader()
        w.writerows(row(o) for o in objs)

    a = sorted((o for o in objs
                if o.get('teff') and A_RANGE[0] <= o['teff'] <= A_RANGE[1]),
               key=lambda o: -o['teff'])
    with open(ROOT / 'data' / 'uves_pop_astars.csv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=COLS)
        w.writeheader()
        w.writerows(row(o) for o in a)
    print(f'A stars ({A_RANGE[0]:.0f}-{A_RANGE[1]:.0f} K): {len(a)}'
          f'  -> data/uves_pop_astars.csv')

    print(f'\n=== near 10,000 K, the Balmer-break window (9000-11000 K) ===')
    print(f'{"name":<14}{"Teff":>7}{"logg":>6}{"[Fe/H]":>8}{"vsini":>7}')
    near = [o for o in a if 9000 <= o['teff'] <= 11000]
    for o in near:
        print(f'{o["name"]:<14}{o["teff"]:>7.0f}{o["logg"]:>6.2f}'
              f'{o["fe_h"]:>+8.2f}{o["vsini"]:>7.0f}')
    print(f'  {len(near)} stars')

    names = {o['name'] for o in objs}
    hit = [s for s in NGSL_SAMPLE if s in names]
    print(f'\nNGSL Balmer-break sample present in UVES-POP: '
          f'{hit if hit else "none — the two libraries are disjoint here"}')


if __name__ == '__main__':
    main()
