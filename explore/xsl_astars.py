"""Select Balmer-break candidates from the X-shooter Spectral Library (XSL) DR3.

XSL DR3: 830 spectra of 683 stars, 3500-24800 A at R ~ 10,000, arm-combined
(Verro et al. 2022, A&A 660, A34). Ground-based (VLT/X-shooter), so unlike NGSL
its flux calibration is not space-based -- but it resolves the Balmer lines
~10x better than NGSL and covers the break.

Parameters come from Arentsen et al. (2019), who derived Teff, log g and [Fe/H]
for the XSL stars; the DR3 table itself carries only names and filenames.

Every candidate is vetted for binarity and chemical peculiarity against SIMBAD
and given a reddening estimate, because neither is caught by a parameter cut:
in UVES-POP a spectroscopic binary (HD162630) looked entirely normal on
abundance and rotation, and reddening cost two of nine candidates.

Writes data/xsl_astars.csv (the A-star window) and data/xsl_all.csv.
"""
import csv
import sys
from pathlib import Path

import numpy as np
from astroquery.vizier import Vizier

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common.dust import irsa_ebv, simbad_info

ROOT = Path(__file__).resolve().parent.parent
DR3 = 'J/A+A/660/A34'            # Verro+2022, the DR3 spectrum list
PARAMS = 'J/A+A/627/A138'        # Arentsen+2019, atmospheric parameters
A_RANGE = (7000.0, 11500.0)
BALMER_WINDOW = (9000.0, 11500.0)
SPEC_URL = 'http://xsl.u-strasbg.fr/data/DR3/{fname}.fits'


def norm(n):
    return ' '.join(str(n).split()).upper().replace('  ', ' ')


def val(x):
    try:
        v = float(x)
        return v if np.isfinite(v) and abs(v) < 1e6 else np.nan
    except (TypeError, ValueError):
        return np.nan


def main():
    v = Vizier(row_limit=-1, columns=['*', '_RAJ2000', '_DEJ2000'])
    dr3 = v.get_catalogs(DR3)[0]
    par = v.get_catalogs(PARAMS)[0]
    print(f'XSL DR3: {len(dr3)} spectra;  Arentsen+2019: {len(par)} stars')

    pmap = {}
    for r in par:
        pmap[norm(r['HNAME'])] = dict(
            teff=val(r['Teff']), e_teff=val(r['e_Teff']),
            logg=val(r['logg']), e_logg=val(r['e_logg']),
            feh=val(r['[Fe/H]']), e_feh=val(r['e_[Fe/H]']))

    rows, seen = [], {}
    for r in dr3:
        name = str(r['Star']).strip()
        p = pmap.get(norm(r['SimbadName'])) or pmap.get(norm(name)) or {}
        rows.append(dict(
            star=name, xslid=str(r['XSLID']).strip(),
            simbad=str(r['SimbadName']).strip(),
            filename=str(r['FileName']).strip(),
            ra=float(r['_RAJ2000']), dec=float(r['_DEJ2000']),
            teff=p.get('teff', np.nan), e_teff=p.get('e_teff', np.nan),
            logg=p.get('logg', np.nan), e_logg=p.get('e_logg', np.nan),
            feh=p.get('feh', np.nan), e_feh=p.get('e_feh', np.nan),
            comment=str(r['Com']).strip()))
    with open(ROOT / 'data' / 'xsl_all.csv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader(); w.writerows(rows)
    have = sum(1 for r in rows if np.isfinite(r['teff']))
    print(f'  parameters matched for {have}/{len(rows)} spectra')

    # one entry per star in the A window (XSL repeats stars across epochs)
    cand = {}
    for r in rows:
        t = r['teff']
        if not np.isfinite(t) or not (A_RANGE[0] <= t <= A_RANGE[1]):
            continue
        cand.setdefault(r['star'], []).append(r)
    print(f'  A-star window {A_RANGE[0]:.0f}-{A_RANGE[1]:.0f} K: '
          f'{len(cand)} stars, {sum(len(v) for v in cand.values())} spectra')

    out = []
    print(f'\n{"star":<16}{"Teff":>7}{"logg":>6}{"[Fe/H]":>8}{"n_ep":>5}'
          f'  {"SIMBAD":<16}{"otype":<7}{"SFD98":>7}{"SF11":>7}  flags')
    for star, eps in sorted(cand.items(), key=lambda kv: -kv[1][0]['teff']):
        e = eps[0]
        info = simbad_info(e['simbad'] or star)
        try:
            sfd, sf11 = irsa_ebv(e['ra'], e['dec'])
        except Exception:
            sfd = sf11 = np.nan
        flags = []
        if info.get('binary'):
            flags.append('BINARY')
        if info.get('peculiar'):
            flags.append('PECULIAR')
        rec = dict(star=star, xslid=e['xslid'], n_epochs=len(eps),
                   teff=round(e['teff']), e_teff=round(e['e_teff']) if np.isfinite(e['e_teff']) else '',
                   logg=round(e['logg'], 2), e_logg=round(e['e_logg'], 2) if np.isfinite(e['e_logg']) else '',
                   feh=round(e['feh'], 2), e_feh=round(e['e_feh'], 2) if np.isfinite(e['e_feh']) else '',
                   ra=round(e['ra'], 5), dec=round(e['dec'], 5),
                   sp_type=info.get('sp_type', ''), otype=info.get('otype', ''),
                   dist_pc=info.get('dist_pc', ''),
                   ebv_sfd98=round(sfd, 4) if np.isfinite(sfd) else '',
                   ebv_sf11=round(sf11, 4) if np.isfinite(sf11) else '',
                   binary='yes' if info.get('binary') else 'no',
                   peculiar='yes' if info.get('peculiar') else 'no',
                   filename=e['filename'],
                   spec_url=SPEC_URL.format(fname=e['filename']),
                   comment=e['comment'])
        out.append(rec)
        if BALMER_WINDOW[0] <= e['teff'] <= BALMER_WINDOW[1]:
            print(f'{star:<16}{e["teff"]:>7.0f}{e["logg"]:>6.2f}{e["feh"]:>8.2f}'
                  f'{len(eps):>5}  {info.get("sp_type",""):<16}'
                  f'{info.get("otype",""):<7}'
                  f'{(f"{sfd:.3f}" if np.isfinite(sfd) else "--"):>7}'
                  f'{(f"{sf11:.3f}" if np.isfinite(sf11) else "--"):>7}'
                  f'  {",".join(flags)}')

    with open(ROOT / 'data' / 'xsl_astars.csv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0]))
        w.writeheader(); w.writerows(out)
    clean = [r for r in out if r['binary'] == 'no' and r['peculiar'] == 'no'
             and BALMER_WINDOW[0] <= r['teff'] <= BALMER_WINDOW[1]]
    print(f'\n{len(out)} A stars -> data/xsl_astars.csv')
    print(f'in the Balmer window and clean of binarity/peculiarity: {len(clean)}')


if __name__ == '__main__':
    main()
