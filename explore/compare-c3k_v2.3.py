import sys
from pathlib import Path
import numpy as np
import h5py

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from fitting.model import Grid
from fitting.predict import spectrum_at
from fitting.observations import load_ngsl
from common.lsf import broaden_R, broaden_rot, to_ngsl_pixels
from common.extinction_ccm import redden

#NODE_ATOL["logt"] = 1e-6  # log10(teff) is the axis in the c3k cube, not teff itself

c3k_dir = Path("/Users/bjohnson/Projects/c3k-fsps-lib/output/c3k_v2.3/vt10_allfal")

balmer = Grid()
c3k = h5py.File(str(c3k_dir / "spec" / "c3k_v2.3_feh+0.00_afe+0.0.spec.h5"), 'r')

params = c3k["parameters"][:]

star = dict(teff=10500, logg=4.2, feh=-0.1, ebv=0.085, vsini=60, name="HD167946")
star["logt"] = np.log10(star["teff"])

ngsl = load_ngsl(star["name"])

sel = np.ones(len(params), dtype=bool)
for p in ["logt", "logg", "feh"]:
    delta = star[p] - params[p]
    diff = np.min(np.abs(delta))
    sel = sel & (np.isclose(np.abs(delta), diff))
ind = np.where(sel)[0]
print(params[ind], star)

balmer_spectrum = spectrum_at(balmer, star["teff"], star["logg"], star["feh"])
c3k_spectrum = c3k["spectra"][ind[0], :]


balmer_wave = balmer.wave
c3k_wave = c3k["wavelengths"][:]

valid = (3200 < c3k_wave) & (c3k_wave < 10000)
c3k_spectrum = c3k_spectrum[valid]
c3k_wave = c3k_wave[valid]
c3k_spectrum /= c3k_wave**2 # to f_lambda

c3k_spectrum = redden(c3k_wave, c3k_spectrum, star["ebv"], 3.1)
c3k_ngsl = to_ngsl_pixels(c3k_wave, c3k_spectrum, ngsl.wavelength)
c3k_r10k = broaden_R(c3k_wave, c3k_spectrum, 10000,)

balmer_spectrum = redden(balmer_wave, balmer_spectrum, star["ebv"], 3.1)
balmer_ngsl = to_ngsl_pixels(balmer_wave, balmer_spectrum, ngsl.wavelength)
balmer_r10k = broaden_R(balmer_wave, balmer_spectrum, 10000,)


node_spectrum = spectrum_at(balmer, 10**params[ind][0]["logt"], params[ind][0]["logg"], params[ind][0]["feh"])
node_spectrum = redden(balmer_wave, node_spectrum, star["ebv"], 3.1)
node_r10k = broaden_R(balmer_wave, node_spectrum, 10000,)
node_ngsl = to_ngsl_pixels(balmer_wave, node_spectrum, ngsl.wavelength)


# renorm
balmer_ngsl *= ngsl.flux[1100:1200].mean() / balmer_ngsl[1100:1200].mean()
c3k_ngsl *= ngsl.flux[1100:1200].mean() / c3k_ngsl[1100:1200].mean()
node_ngsl *= ngsl.flux[1100:1200].mean() / node_ngsl[1100:1200].mean()
renorm = c3k_r10k[(c3k_wave < 3500) & (c3k_wave > 3200)].mean() / node_r10k[(balmer_wave < 3500) & (balmer_wave > 3200)].mean()
node_r10k *= renorm

node_r10k_int = np.interp(c3k_wave, balmer_wave, node_r10k)

import matplotlib.pyplot as pl
pl.style.use("via")
pl.ion()
fig, axes = pl.subplots(2, 1, figsize=(14, 7), sharex=True)
ax = axes[0]
ax.plot(ngsl.wavelength, ngsl.flux, label="NGSL", color="k", alpha=0.8,)
ax.plot(ngsl.wavelength, balmer_ngsl, label="Balmer", color="dodgerblue", alpha=0.8)
ax.plot(ngsl.wavelength, c3k_ngsl, label=f"C3K_v2.3\nlogg={params[ind][0]['logg']},[Fe/H]={params[ind][0]['feh']}", color="r", alpha=0.8)
ax.plot(ngsl.wavelength, node_ngsl, label="Balmer @ C3K", color="cyan", alpha=0.8)
#ax.set_xlabel("Wavelength [Angstrom]")
ax.set_ylabel("Flux [erg/s/cm^2/Angstrom]", fontsize=13)
ax.set_title(f"Star: {star['name']} (Teff={star['teff']}, logg={star['logg']}, [Fe/H]={star['feh']})")
ax.set_xlim(3600, 4000)
ax.set_ylim(0.3e-11, 1.19e-11)
ax.legend()

ax = axes[1]
ax.plot(ngsl.wavelength,(c3k_ngsl - node_ngsl)/node_ngsl*100, label="C3K - Balmer@C3K", color="red", alpha=0.8)
ax.plot(ngsl.wavelength,(balmer_ngsl - ngsl.flux)/ngsl.flux*100, label="Balmer - NGSL", color="dodgerblue", alpha=0.8)
ax.plot(ngsl.wavelength, np.zeros_like(ngsl.wavelength), color="k", alpha=0.8, ls="--")
ax.set_xlim(3600, 4000)
ax.set_ylim(-10, 10)
ax.set_xlabel("Wavelength [Angstrom]")
ax.legend()
ax.set_ylabel("Residual [%]", fontsize=13)
fig.savefig(f"../figures/compare_c3k_v2.3_{star['name']}.png", dpi=300, bbox_inches="tight")

# --------
fig, axes = pl.subplots(2, 1, figsize=(14, 7), sharex=True)
ax = axes[0]
ax.plot(c3k_wave, c3k_r10k, label="C3K", color="k", alpha=0.8,)
ax.plot(balmer_wave, node_r10k, label="Balmer @ C3K", color="cyan", alpha=0.8)
ax.set_xlim(3600, 4000)
#ax.set_ylim(0.3e-11, 1.19e-11)
ax.set_ylabel("Flux [erg/s/cm^2/Angstrom]", fontsize=13)
ax.set_title(f'Teff={10**params[ind][0]["logt"]:.0f} K, logg={params[ind][0]["logg"]:.2f}, [Fe/H]={params[ind][0]["feh"]:.2f}')
ax.legend()

ax = axes[1]
ax.plot(c3k_wave, (c3k_r10k - node_r10k_int)/node_r10k_int*100, label="C3K - Balmer@C3K", color="red", alpha=0.8)
ax.set_xlim(3600, 4000)
ax.set_ylim(-15, 15)
ax.set_xlabel("Wavelength [Angstrom]")
ax.legend()
ax.set_ylabel("Residual [%]", fontsize=13)
fig.savefig(f"../figures/compare_c3k_v2.3_R10K.png", dpi=300, bbox_inches="tight")
