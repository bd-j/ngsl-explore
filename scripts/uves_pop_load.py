"""Load UVES-POP spectra and their parameters, on the same footing as the models.

Two conversions are mandatory before any comparison:

  air -> vacuum   The FITS headers declare CTYPE1 = AWAV. ATLAS12/SYNTHE output
                  and the NGSL spectra are both vacuum, so UVES-POP must be
                  converted or every line sits ~1.1 A blue of where it belongs.

  resolution      The library is R = 80,000 natively, but the delivered product
                  is resampled to a 0.1 A linear grid, so the usable resolution
                  is ~R = 18,000 at the Balmer break (2-px), not 80,000.

Parameters come in two flavours in the FITS metadata table and both are kept:
VSINI/TEFF/LOGG/FE_H are the library's own spectral fit (against a PHOENIX grid,
GRID_NAME = phx20atm), LIT_* are literature values. Use VSINI as the starting
point for model broadening; refit if the profiles demand it.

Writes data/uves_pop_selected.csv
"""
import csv
from pathlib import Path

import numpy as np
from astropy.io import fits

ROOT = Path(__file__).resolve().parent.parent
UVES = ROOT / 'data' / 'uves_pop'
NATIVE_R = 80000.0
GRID_STEP = 0.1           # A, delivered sampling

# Same reddening cut as the NGSL sample (candidate_table.EBV_MAX). UVES-POP
# publishes a fitted E(B-V) with an error, so in principle these could be
# dereddened rather than rejected -- but the sample is held to one standard.
# HD162817 (E(B-V) = 0.102) was dropped under this cut.
EBV_MAX = 0.10


def air_to_vac(w):
    """Ciddor (1996) via the IAU standard inverse; good to <1e-4 A here."""
    s = 1e4 / np.asarray(w, dtype=float)
    n = (1 + 0.05792105 / (238.0185 - s * s) + 0.00167917 / (57.362 - s * s))
    return np.asarray(w, dtype=float) * n


def load(star, to_vacuum=True):
    """-> (wavelength, flux, error). Flux is erg/cm^2/s/A, absolute."""
    with fits.open(UVES / f'{star}.fits.gz') as f:
        h = f[0].header
        n = h['NAXIS1']
        w = h['CRVAL1'] + (np.arange(n) + 1 - h['CRPIX1']) * h['CDELT1']
        fl, er = np.asarray(f[0].data, float), np.asarray(f[1].data, float)
    if to_vacuum:
        w = air_to_vac(w)
    return w, fl, er


def params(star):
    with fits.open(UVES / f'{star}.fits.gz') as f:
        t = f[2].data
        g = lambda c: (float(t[c][0]) if c in f[2].columns.names
                       and np.ndim(t[c][0]) == 0 else '')
        s = lambda c: (str(t[c][0]).strip() if c in f[2].columns.names else '')
        return dict(
            star=star, teff=g('TEFF'), e_teff=g('E_TEFF'), logg=g('LOGG'),
            e_logg=g('E_LOGG'), fe_h=g('FE_H'), e_fe_h=g('E_FE_H'),
            vsini=g('VSINI'), e_vsini=g('E_VSINI'), ebv=g('EBV'),
            ebv_err=g('EBV_ERR'), lit_teff=g('LIT_TEFF'), lit_logg=g('LIT_LOGG'),
            lit_fe_h=g('LIT_FE_H'), lit_vsini=g('LIT_VSINI'),
            grid=s('GRID_NAME'))


def effective_R(lam):
    """Usable resolving power of the delivered grid (2-px), capped at native."""
    return np.minimum(lam / (2 * GRID_STEP), NATIVE_R)


if __name__ == '__main__':
    stars = sorted(p.name.replace('.fits.gz', '') for p in UVES.glob('*.fits.gz'))
    rows = [params(s) for s in stars]
    for r, s in zip(rows, stars):
        w, fl, er = load(s)
        r['wave_min'] = round(float(w[0]), 1)
        r['wave_max'] = round(float(w[-1]), 1)
        r['snr_median'] = round(float(np.nanmedian(fl / er)))
        b = (w > 3600) & (w < 4100)
        r['snr_balmer'] = round(float(np.nanmedian((fl / er)[b])))
        e = r['ebv']
        r['selected'] = 'no' if (e != '' and e > EBV_MAX) else 'yes'
        r['reject_reason'] = (f'E(B-V)={e:.3f} > {EBV_MAX:.2f}'
                              if r['selected'] == 'no' else '')
    cols = list(rows[0])
    with open(ROOT / 'data' / 'uves_pop_selected.csv', 'w', newline='') as fh:
        wtr = csv.DictWriter(fh, fieldnames=cols)
        wtr.writeheader()
        wtr.writerows(rows)
    print(f'{len(rows)} stars -> data/uves_pop_selected.csv\n')
    hdr = ['star', 'teff', 'logg', 'fe_h', 'vsini', 'lit_vsini', 'ebv',
           'snr_balmer']
    print(''.join(f'{h:>12}' for h in hdr))
    for r in rows:
        print(''.join(
            f'{(f"{r[h]:.2f}" if isinstance(r[h], float) else str(r[h])):>12}'
            for h in hdr))
    print(f'\neffective R at 3646 A = {effective_R(3646.0):.0f} '
          f'(NGSL is ~950; native UVES is 80000 before resampling)')
