"""Load X-shooter Spectral Library (XSL) DR3 spectra.

Format (Verro et al. 2022, A&A 660, A34):
  WAVE      nm, REST-FRAME, log10-sampled at ~R = 30,000 (3 px per resolution
            element). Rest-frame means the RV is already removed, so RV should
            be fixed at 0 when fitting XSL, unlike NGSL.
  FLUX      erg/s/cm^2/A, corrected for slit losses
  FLUX_DR   the same, additionally corrected for Galactic extinction
  ERR       uncertainty

Resolution is quoted as sigma(v), NOT FWHM: 13 km/s UVB, 11 VIS, 16 NIR.
So at the Balmer break (UVB) FWHM = 2.3548 x 13 = 30.6 km/s, i.e. R ~ 9800.
This is constant in VELOCITY -- unlike NGSL, whose LSF is set by a fixed
dispersion per grating and is constant in ANGSTROMS with jumps at the splices.
The two libraries need different convolution kernels.

Overlap regions are smoothed to the worse of the two arms (13 km/s across
UVB/VIS), so the resolution is not uniform across a splice.
"""
from pathlib import Path

import numpy as np
from astropy.io import fits

ROOT = Path(__file__).resolve().parent.parent
XSL = ROOT / 'data' / 'xsl' / 'XSL_DR3_release'

# sigma(v) in km/s per arm, and the wavelength ranges they cover (Angstroms)
ARMS = [('UVB', 3000., 5600., 13.0), ('VIS', 5600., 10200., 11.0),
        ('NIR', 10200., 24800., 16.0)]
C_KMS = 2.99792458e5


def sigma_v(wave_A):
    """sigma(v) in km/s at each wavelength, from the arm it falls in."""
    out = np.full(np.shape(wave_A), np.nan, float)
    for _, lo, hi, s in ARMS:
        m = (np.asarray(wave_A) >= lo) & (np.asarray(wave_A) < hi)
        out[m] = s
    return out


def resolving_power(wave_A):
    """R = c / FWHM(v); FWHM = 2.3548 sigma. Careful: XSL quotes sigma."""
    return C_KMS / (2.3548 * sigma_v(wave_A))


def air_to_vac(w):
    s = 1e4 / np.asarray(w, float)
    n = 1 + 0.05792105 / (238.0185 - s * s) + 0.00167917 / (57.362 - s * s)
    return np.asarray(w, float) * n


def load(xslid, dereddened=False, to_vacuum=True):
    """-> (wave_A, flux, err, header). Wavelengths converted nm -> Angstrom.

    XSL is in AIR. Established by cross-correlating HD194453 -- which is in both
    XSL and NGSL -- against the wavecal-corrected NGSL spectrum: XSL needs
    +1.05 A to match, against a +1.13 A air-vacuum offset at 4000 A. Converted
    to vacuum by default so all three libraries and the models share one scale.

    dereddened=False returns FLUX (slit-loss corrected only) so extinction can
    be handled the same way as for NGSL and UVES-POP; True returns XSL's own
    FLUX_DR, useful as an independent check.
    """
    p = XSL / f'xsl_spectrum_{xslid}_merged.fits'
    with fits.open(p) as f:
        d = f[1].data
        w = d['WAVE'].astype(float) * 10.0
        if to_vacuum:
            w = air_to_vac(w)
        fl = d['FLUX_DR' if dereddened else 'FLUX'].astype(float)
        er = d['ERR'].astype(float)
        hdr = dict(f[0].header)
    ok = np.isfinite(fl) & (fl > 0)
    return w[ok], fl[ok], er[ok], hdr
