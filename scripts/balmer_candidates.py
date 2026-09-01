"""Plot the Balmer break region for the finalist ~10,000 K stars and quantify
the break amplitude. Writes figures/balmer_break_candidates.png
"""
import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from astropy.io import fits

CANDS = ['HD194453', 'HD040573', 'HD147550', 'HD143459']
cat = {r['target']: r for r in csv.DictReader(open('data/ngsl_catalog.csv'))}

fig, axes = plt.subplots(2, 2, figsize=(13, 8), sharex=True)
print(f'{"target":<10}{"Teff":>7}{"logg":>6}{"[M/H]":>7}{"alpha":>7}{"off_px":>8}{"D_Balmer":>10}')

for ax, name in zip(axes.ravel(), CANDS):
    r = cat[name]
    with fits.open('data/spectra/' + r['file']) as h:
        d = h[1].data
    w, fl, er = d['WAVELENGTH'], d['FLUX'], d['STATERR']
    m = (w > 3200) & (w < 4200)
    ax.plot(w[m], fl[m], lw=0.9, color='#1f3b73')
    ax.fill_between(w[m], fl[m] - er[m], fl[m] + er[m], color='#1f3b73', alpha=0.25, lw=0)
    ax.axvline(3646, color='#c1440e', ls='--', lw=1, label='Balmer limit 3646 A')

    # Balmer decrement: mean flux just blue vs just red of the limit
    blue = np.median(fl[(w > 3500) & (w < 3620)])
    red = np.median(fl[(w > 3680) & (w < 3800)])
    d_b = 2.5 * np.log10(red / blue)
    sn = lambda lo, hi: np.median(fl[(w > lo) & (w < hi)] / er[(w > lo) & (w < hi)])
    s35, s37 = sn(3450, 3600), sn(3650, 3800)
    print(f'{name:<10}{float(r["teff"]):7.0f}{float(r["logg"]):6.1f}'
          f'{float(r["logz"]):+7.1f}{r["alpha"]:>7}{float(r["offset_px"]):8.2f}{d_b:10.3f}')

    alpha = {'a': r', $\alpha$-enh', 'n': ''}.get(r['alpha'], '')
    ax.set_title(f'{name}   {r["sptype"]}   Teff={float(r["teff"]):.0f} K   '
                 f'logg={r["logg"]}   [M/H]={float(r["logz"]):+.1f}{alpha}   '
                 f'offset={float(r["offset_px"]):+.2f} px', fontsize=9.5)
    ax.set_ylabel(r'F$_\lambda$ [erg cm$^{-2}$ s$^{-1}$ $\AA^{-1}$]', fontsize=8)
    ax.legend(fontsize=7, loc='lower right')
    ax.tick_params(labelsize=8)

for ax in axes[1]:
    ax.set_xlabel(r'Wavelength [$\AA$]')
fig.suptitle('NGSL v2: Balmer break region for ~10,000 K candidates  (NGSL metallicities are scaled-solar [M/H], not [Fe/H])', fontsize=11)
fig.tight_layout()
fig.savefig('figures/balmer_break_candidates.png', dpi=140)
print('\n-> figures/balmer_break_candidates.png')
