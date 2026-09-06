"""Stellar-parameter and S/N coverage of NGSL v2.
Writes figures/parameter_coverage.png and figures/snr_vs_wavelength.png
"""
import csv, collections
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

cat = list(csv.DictReader(open('data/ngsl_catalog.csv')))
snr = {r['target']: r for r in csv.DictReader(open('data/ngsl_snr.csv'))}
xm = {r['target']: r for r in csv.DictReader(open('data/ngsl_crossmatch.csv'))}
f = lambda r, k: float(r[k]) if r[k] not in ('', None) else np.nan

teff = np.array([f(r, 'teff') for r in cat])
logg = np.array([f(r, 'logg') for r in cat])
logz = np.array([f(r, 'logz') for r in cat])
vmag = np.array([f(r, 'vmag') for r in cat])
good = np.array([r['dataqual'] == 'good' for r in cat])
ok = np.isfinite(teff)

print('=== NGSL v2 stellar parameter coverage ===')
print(f'  stars total                 : {len(cat)}')
print(f'  data quality "good"         : {good.sum()}  ("suspect": {(~good).sum()})')
print(f'  with model-fit parameters   : {ok.sum()}')
print(f'  in MILES                    : {sum(1 for r in cat if xm[r["target"]]["in_miles"]=="Y")}')
for lbl, a, unit in [('Teff', teff, 'K'), ('log g', logg, 'dex'),
                     ('[M/H]', logz, 'dex'), ('V', vmag, 'mag')]:
    v_ = a[np.isfinite(a)]
    print(f'  {lbl:<6}: {v_.min():8.2f} to {v_.max():8.2f} {unit:4s}'
          f' median {np.median(v_):8.2f}')

print('\n  Teff distribution:')
for lo, hi in [(3000,4000),(4000,5000),(5000,6000),(6000,7000),(7000,8000),
               (8000,10000),(10000,15000),(15000,50000)]:
    n = ((teff >= lo) & (teff < hi)).sum()
    print(f'    {lo:5d}-{hi:5d} K : {n:3d}  {"#"*int(n/2)}')

sp = collections.Counter(r['sptype'][:1] for r in cat if r['sptype'][:1] in 'OBAFGKM')
print('\n  spectral types:', ', '.join(f'{k}={sp[k]}' for k in 'OBAFGKM' if k in sp))

fig, ax = plt.subplots(1, 3, figsize=(15, 4.4))
s = ax[0].scatter(teff[ok], logg[ok], c=logz[ok], cmap='viridis', s=22,
                  edgecolor='k', linewidth=0.3)
ax[0].invert_xaxis(); ax[0].invert_yaxis()
ax[0].set_xlabel('Teff [K]'); ax[0].set_ylabel('log g')
ax[0].set_title('Kiel diagram, coloured by [M/H]')
plt.colorbar(s, ax=ax[0], label='[M/H] (scaled-solar)')

ax[1].scatter(teff[ok], logz[ok], s=22, c='#1f3b73', edgecolor='k', linewidth=0.3)
ax[1].invert_xaxis()
ax[1].set_xlabel('Teff [K]'); ax[1].set_ylabel('[M/H]')
ax[1].set_title('Metallicity coverage: [M/H] (scaled-solar, not [Fe/H])')

ax[2].hist(teff[ok], bins=np.arange(3000, 15000, 500),
           color='#1f3b73', edgecolor='white')
ax[2].set_xlabel('Teff [K]'); ax[2].set_ylabel('N stars')
ax[2].set_title(f'Teff distribution (N={ok.sum()})')
ax[2].axvspan(9500, 10500, color='#c1440e', alpha=0.25)
fig.suptitle('NGSL v2 parameter coverage: Teff, log g and scaled-solar [M/H] '
             '(Castelli models on Victoria-Regina isochrones; 107/324 stars alpha-enhanced)')
fig.tight_layout(); fig.savefig('figures/parameter_coverage.png', dpi=140)

# --- S/N vs wavelength ---
BANDS = [('2000_2500',2250),('2500_3000',2750),('3000_3600',3300),('3600_4000',3800),
         ('4000_5000',4500),('5000_5500',5250),('5700_7000',6350),
         ('7000_9000',8000),('9000_10000',9500)]
fig, ax = plt.subplots(figsize=(8.5, 5))
for pre, lbl, col in [('stat', 'formal (STATERR, optimistic ~3x)', '#c1440e'),
                      ('der', 'DER_SNR (pessimistic: counts lines as noise)', '#1f3b73')]:
    med, p16, p84, xs = [], [], [], []
    for b, x in BANDS:
        a = np.array([f(snr[r['target']], f'{pre}_{b}') for r in cat])
        a = a[np.isfinite(a)]
        if len(a) < 5: continue
        xs.append(x); med.append(np.median(a))
        p16.append(np.percentile(a, 16)); p84.append(np.percentile(a, 84))
    ax.plot(xs, med, 'o-', color=col, label=lbl)
    ax.fill_between(xs, p16, p84, color=col, alpha=0.18, lw=0)
ax.axhline(100, ls=':', c='k', lw=1)
ax.text(9300, 108, 'true S/N ~100/pix (measured)', fontsize=8, ha='right')
ax.axvspan(3500, 3800, color='#888', alpha=0.2)
ax.text(3650, 3.5, 'Balmer\nbreak', fontsize=8, ha='center')
ax.set_yscale('log'); ax.set_xlabel(r'Wavelength [$\AA$]')
ax.set_ylabel('S/N per pixel'); ax.legend(fontsize=8)
ax.set_title('NGSL v2 signal-to-noise (median, 16-84th percentile band)')
fig.tight_layout(); fig.savefig('figures/snr_vs_wavelength.png', dpi=140)
print('\n-> figures/parameter_coverage.png, figures/snr_vs_wavelength.png')
