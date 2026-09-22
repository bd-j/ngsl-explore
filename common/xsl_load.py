"""Load X-shooter Spectral Library (XSL) DR3 spectra.

Format (Verro et al. 2022, A&A 660, A34):
  WAVE      nm, REST-FRAME, log10-sampled at ~R = 30,000 (3 px per resolution
            element). Rest-frame means the RV is already removed, so RV should
            be fixed at 0 when fitting XSL, unlike NGSL.
  FLUX      erg/s/cm^2/A, corrected for slit losses -- BUT ONLY WHEN THE
            HEADER SAYS SO, see the filename variants below
  FLUX_DR   the same, additionally corrected for Galactic extinction
  ERR       uncertainty

NOT EVERY SPECTRUM IS SLIT-LOSS CORRECTED, and the filename says which:

  <id>_merged.fits            LOSS_COR = True,  columns FLUX, FLUX_DR, ERR
  <id>_merged_scl.fits        LOSS_COR = False, columns FLUX, FLUX_SC, ERR
  <id>_merged_ncl.fits        LOSS_COR = False
  <id>_merged_ncge.fits       EXT_AVG  = False, no FLUX_DR
  <id>_merged_ncl_ncge.fits   both False, columns FLUX, ERR only

Only 606 of the 830 DR3 spectra are the plain `_merged.fits`, so a path built
as f'xsl_spectrum_{xslid}_merged.fits' is wrong for 224 of them -- use
`spectrum_path`, which reads the filename recorded in data/xsl_all.csv.

THIS MATTERS FOR ABSOLUTE FLUX, NOT FOR SHAPE. A spectrum with
LOSS_COR = False carries an uncorrected slit loss, so its absolute level
cannot be compared against another library's; its shape still can, after
normalisation. `load` returns the header so callers can check, and
`loss_corrected` reads the flag directly.

Resolution is quoted as sigma(v), NOT FWHM: 13 km/s UVB, 11 VIS, 16 NIR.
So at the Balmer break (UVB) FWHM = 2.3548 x 13 = 30.6 km/s, i.e. R ~ 9800.
This is constant in VELOCITY -- unlike NGSL, whose LSF is set by a fixed
dispersion per grating and is constant in ANGSTROMS with jumps at the splices.
The two libraries need different convolution kernels.

Overlap regions are smoothed to the worse of the two arms (13 km/s across
UVB/VIS), so the resolution is not uniform across a splice.
"""
import sys
from pathlib import Path

import numpy as np
from astropy.io import fits

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common.lines import air_to_vac

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


_FILENAMES = {}


def spectrum_path(xslid):
    """-> Path of this XSL ID's spectrum, whichever correction variant it is.

    The filename is read from data/xsl_all.csv rather than assumed, because
    224 of the 830 spectra carry a _scl/_ncl/_ncge suffix (see the module
    docstring). Falls back to a glob so a spectrum extracted without the
    catalog still loads.
    """
    if not _FILENAMES:
        import csv
        with open(ROOT / 'data' / 'xsl_all.csv') as fh:
            for r in csv.DictReader(fh):
                _FILENAMES[r['xslid']] = r['filename']
    name = _FILENAMES.get(xslid)
    if name and (XSL / name).exists():
        return XSL / name
    hits = sorted(XSL.glob(f'xsl_spectrum_{xslid}_merged*.fits'))
    if not hits:
        raise FileNotFoundError(
            f'no XSL spectrum for {xslid} in {XSL}. Extract it from the '
            f'tarball: tar -xf data/xsl/XSL_DR3_release.tar -C data/xsl '
            f'XSL_DR3_release/{name or "xsl_spectrum_" + xslid + "_merged.fits"}')
    return hits[0]


def loss_corrected(xslid):
    """-> True if this spectrum's FLUX carries the slit-loss correction."""
    with fits.open(spectrum_path(xslid)) as f:
        return bool(f[0].header.get('LOSS_COR', False))


def load(xslid, dereddened=False, to_vacuum=True):
    """-> (wave_A, flux, err, header). Wavelengths converted nm -> Angstrom.

    XSL is in AIR. Established by cross-correlating HD194453 -- which is in both
    XSL and NGSL -- against the wavecal-corrected NGSL spectrum: XSL needs
    +1.05 A to match, against a +1.13 A air-vacuum offset at 4000 A. Converted
    to vacuum by default so all three libraries and the models share one scale.

    dereddened=False returns FLUX so extinction can be handled the same way as
    for NGSL and UVES-POP; True returns XSL's own dereddened column, useful as
    an independent check. Check header['LOSS_COR'] before using FLUX as an
    ABSOLUTE flux: it is False for the _scl/_ncl variants.
    """
    p = spectrum_path(xslid)
    with fits.open(p) as f:
        d = f[1].data
        w = d['WAVE'].astype(float) * 10.0
        if to_vacuum:
            w = air_to_vac(w)
        if dereddened:
            # The dereddened column is FLUX_DR in the plain files and FLUX_SC
            # in the _scl ones; the _ncge variants have neither.
            col = next((c for c in ('FLUX_DR', 'FLUX_SC')
                        if c in d.columns.names), None)
            if col is None:
                raise KeyError(
                    f'{p.name} has no dereddened column (columns: '
                    f'{d.columns.names}); it is an _ncge variant, so use '
                    f'dereddened=False and redden the model instead')
        else:
            col = 'FLUX'
        fl = d[col].astype(float)
        er = d['ERR'].astype(float)
        hdr = dict(f[0].header)
    ok = np.isfinite(fl) & (fl > 0)
    return w[ok], fl[ok], er[ok], hdr
