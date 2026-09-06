"""Figures for the Pickles atlas report."""
import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from astropy.io import fits

UV, TEAL, AMBER, INK3 = '#4a3fb5', '#0f6d75', '#b07d10', '#727a94'
pk = list(csv.DictReader(open('data/pickles_catalog.csv')))
by_type = {r['sptype']: r for r in pk}


def spec(rec):
    d = fits.getdata('data/pickles/' + rec['file'])
    return d['WAVELENGTH'].astype(float), d['FLUX'].astype(float)


# --- Fig 1: A0 luminosity sequence at ~constant Teff ------------------------
fig, ax = plt.subplots(1, 2, figsize=(13, 4.8))
seq = [('A0V', UV, 'V  dwarf'), ('A0IV', TEAL, 'IV  subgiant'),
       ('A0III', AMBER, 'III  giant'), ('A0I', '#b0303a', 'I  supergiant')]
for sp, col, lbl in seq:
    r = by_type[sp]
    w, f = spec(r)
    m = (w > 3200) & (w < 4300)
    norm = np.median(f[(w > 4000) & (w < 4200)])
    ax[0].plot(w[m], f[m] / norm, color=col, lw=1.4,
               label=f'{sp}  {r["teff_K"]} K   D={r["d_balmer_mag"]}')
    v = float(r['d_balmer_mag'])
    ax[1].bar(lbl, v, color=col, width=.6)
    ax[1].text(lbl, v + .02, f'{v:.3f}', ha='center', fontsize=9)
ax[0].axvline(3646, color=INK3, ls='--', lw=1)
ax[0].text(3646, 0.30, ' Balmer limit', fontsize=8, color=INK3)
ax[0].set_xlabel(r'Wavelength [$\AA$]')
ax[0].set_ylabel(r'F$_\lambda$ / F$_\lambda$(4000-4200 $\AA$)')
ax[0].set_title('Pickles A0 sequence: gravity at nearly fixed Teff (9550-9727 K)',
                fontsize=10)
ax[0].legend(fontsize=8, loc='lower right')
ax[1].set_ylabel(r'D$_{Balmer}$ [mag]')
ax[1].set_ylim(0, 1.35)
ax[1].set_title('Break amplitude vs luminosity class', fontsize=10)
ax[1].tick_params(labelsize=9)
fig.tight_layout()
fig.savefig('figures/pickles_a0_sequence.png', dpi=140)

# --- Fig 2: D_Balmer vs Teff across the atlas -------------------------------
fig, ax = plt.subplots(figsize=(9, 5.2))
marks = {'V': ('o', UV), 'IV': ('s', TEAL), 'III': ('^', AMBER),
         'II': ('D', '#7a5ea8'), 'I': ('v', '#b0303a')}
for lc, (mk, col) in marks.items():
    # Metric is only interpretable above ~6000 K; see balmer_metric.py.
    xs = [(float(r['teff_K']), float(r['d_balmer_mag'])) for r in pk
          if r['teff_K'] and r['d_balmer_mag'] and r['lum_class'] == lc
          and float(r['teff_K']) >= 6000]
    if xs:
        ax.scatter(*zip(*xs), marker=mk, c=col, s=42, label=f'class {lc}',
                   edgecolor='white', linewidth=.6, zorder=3)
ax.axvspan(9000, 11500, color=AMBER, alpha=.13, zorder=0)
ax.set_xlim(5800, 42000)
ax.text(10200, 0.12, 'Balmer-break\ntest window', ha='center', fontsize=8, color=INK3)
ax.set_xscale('log')
ax.set_xticks([6000, 8000, 10000, 15000, 20000, 40000])
ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax.set_xlabel('Teff [K]')
ax.set_ylabel(r'D$_{Balmer}$ [mag]')
ax.set_title('Balmer break across the Pickles atlas, Teff > 6000 K\n(below that the metric tracks SED slope, not the discontinuity)', fontsize=10.5)
ax.legend(fontsize=8)
ax.grid(alpha=.15)
fig.tight_layout()
fig.savefig('figures/pickles_break_vs_teff.png', dpi=140)

# --- Fig 3: Pickles A0V vs the NGSL recommendation --------------------------
fig, ax = plt.subplots(figsize=(9.5, 5))
w, f = spec(by_type['A0V'])
m = (w > 3300) & (w < 4200)
ax.plot(w[m], f[m] / np.median(f[(w > 4000) & (w < 4200)]), color=AMBER, lw=2,
        label='Pickles A0V template (R~500, 5 $\\AA$/px, composite)')
d = fits.getdata('data/spectra/h_stis_ngsl_hd194453_v2.fits')
wn, fn = d['WAVELENGTH'], d['FLUX']
mn = (wn > 3300) & (wn < 4200)
ax.plot(wn[mn], fn[mn] / np.median(fn[(wn > 4000) & (wn < 4200)]), color=UV,
        lw=1.0, alpha=.9,
        label='NGSL HD194453 (R~940 at break, single star)')
ax.axvline(3646, color=INK3, ls='--', lw=1)
ax.set_xlabel(r'Wavelength [$\AA$]')
ax.set_ylabel(r'F$_\lambda$ / F$_\lambda$(4000-4200 $\AA$)')
ax.set_title('Template vs individual star: the same break at two resolutions',
             fontsize=11)
ax.legend(fontsize=8.5, loc='upper left')
fig.tight_layout()
fig.savefig('figures/pickles_vs_ngsl.png', dpi=140)

# --- Fig 4: three-library coverage grid -------------------------------------
ov = list(csv.DictReader(open('data/library_overlap.csv')))
CL, LU = list('OBAFGKM'), ['I', 'II', 'III', 'IV', 'V']
fig, axes = plt.subplots(1, 3, figsize=(14, 3.9))
for ax_, (key, title, cmap) in zip(axes, [
        ('ngsl_stars', 'NGSL (individual stars)', 'Purples'),
        ('miles_stars', 'MILES (individual stars)', 'GnBu'),
        ('pickles_templates', 'Pickles (composite templates)', 'YlOrBr')]):
    M = np.array([[int(next(r[key] for r in ov if r['spectral_class'] == c
                            and r['lum_class'] == l)) for l in LU] for c in CL])
    im = ax_.imshow(M, cmap=cmap, aspect='auto',
                    norm=matplotlib.colors.PowerNorm(0.5, vmin=0, vmax=max(M.max(), 1)))
    ax_.set_xticks(range(len(LU)), LU)
    ax_.set_yticks(range(len(CL)), CL)
    ax_.set_title(title, fontsize=10)
    ax_.set_xlabel('luminosity class')
    for i in range(len(CL)):
        for j in range(len(LU)):
            ax_.text(j, i, M[i, j], ha='center', va='center', fontsize=8.5,
                     color='white' if M[i, j] > M.max() * .55 else '#333')
axes[0].set_ylabel('spectral class')
fig.suptitle('Where the three libraries overlap: 26 of 42 class x luminosity cells are in all three',
             fontsize=11)
fig.tight_layout()
fig.savefig('figures/library_coverage_grid.png', dpi=140)
print('-> 4 figures written')
