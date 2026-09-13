"""Observations: one record per dataset, each knowing how to be compared.

An `Observation` carries the data plus three policies that are properties of the
INSTRUMENT, not of the fit:

  resolution    how to smooth a model to this dataset
  calibration   which linear nuisance terms are marginalised away
  masking       which pixels are used at all

Keeping them on the observation rather than in the likelihood is what lets one
`predict()` and one likelihood serve NGSL, XSL and photometry without a branch
per dataset -- and it is why the NGSL LSF cannot silently disagree between the
fitter and the figures, which it has done before (a request for R=600 became
R=83).

The calibration policies, and why each dataset gets the one it does:

  'none'     trust the absolute flux as delivered.
  'scalar'   one multiplicative constant. For a star this is (R/d)^2 plus any
             grey calibration error. NGSL gets this: its SHAPE is trusted -- it
             is space-based spectrophotometry, which is the only thing here
             trusted for absolute calibration over a wide baseline -- but its
             absolute level carries the unknown stellar radius and distance.
  'poly', n  a Chebyshev of order n, marginalised. XSL gets this: it resolves
             lines ~16x better than NGSL but is ground-based and must correct
             slit losses, so its continuum is not trusted. Marginalising it
             removes continuum shape from the likelihood and leaves the line
             profiles -- which are immune to reddening, and therefore carry
             Teff, log g and v sin i free of the dust degeneracy.

Every dataset gets its own calibration, so nothing depends on two instruments
agreeing in absolute flux. They happen to agree to 1.1% (Gaia XP vs NGSL for
HD194453), but the fit does not lean on that.
"""
import csv
import sys
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common.ngsl_wavecal import apply_wavecal, load_table
from common.xsl_load import load as xsl_load, ARMS, C_KMS

ROOT = Path(__file__).resolve().parent.parent

# The model grid covers 3200-9500 A (grid/pack_grid.py), and CCM89 here is
# implemented for 3030-9090 A on the optical branch plus the IR branch above it.
# Nothing outside the grid can be predicted, so it is dropped at load time
# rather than extrapolated: np.interp would return edge values in silence, which
# once produced a 15.8% residual that was pure artifact (CAVEATS.md).
MODEL_RANGE = (3200.0, 9500.0)

# NGSL's quoted STATERR holds propagated counting statistics only and is
# optimistic by ~3x; real continuum scatter gives S/N ~ 100 against a claimed
# ~330. The error model handles the amplitude, but a floor is applied here so no
# single pixel claims implausible weight.
NGSL_SNR_CEILING = 200.0


@dataclass
class Observation:
    """One dataset for one star."""
    name: str                       # 'ngsl', 'xsl', 'phot'
    star: str
    flux: np.ndarray                # f_lambda, or maggies for photometry
    uncertainty: np.ndarray
    mask: np.ndarray                # True = use in the likelihood
    wavelength: np.ndarray = None   # vacuum A; None for photometry
    filters: list = None            # sedpy filters, photometry only
    resolution: tuple = (None,)     # ('R', 600.) | ('R_segments', [(lo,hi,R)..])
    calibration: tuple = ('scalar',)
    rv_fixed: float = None          # None = free; 0.0 = already rest-frame
    meta: dict = field(default_factory=dict)

    @property
    def ndata(self):
        return int(self.mask.sum())

    def __repr__(self):
        w = ('' if self.wavelength is None else
             f' {self.wavelength[self.mask].min():.0f}-'
             f'{self.wavelength[self.mask].max():.0f}A')
        return (f'<Observation {self.star}/{self.name} n={self.ndata}{w} '
                f'res={self.resolution[0]} cal={self.calibration[0]}>')


def _sample_row(star):
    for r in csv.DictReader(open(ROOT / 'data' / 'sample.csv')):
        if r['star'] == star:
            return r
    raise KeyError(f'{star} is not in data/sample.csv')


def load_ngsl(star, exclude=None, snr_ceiling=NGSL_SNR_CEILING):
    """NGSL v2 spectrum, air->vacuum and wavecal-corrected.

    NGSL is delivered in AIR with a per-grating residual that is a SLOPE, not a
    zero point (+1.11 A at 3500 falling to +0.08 at 5350), because no wavecals
    were taken with the stellar exposures. `apply_wavecal` removes both.
    """
    from astropy.io import fits
    row = _sample_row(star)
    d = fits.getdata(ROOT / 'data' / 'spectra' / row['ngsl_file'])
    w = apply_wavecal(d['WAVELENGTH'].astype(float), star, load_table())
    f = d['FLUX'].astype(float)
    e = d['STATERR'].astype(float)

    ok = (np.isfinite(w) & np.isfinite(f) & np.isfinite(e) & (f > 0) & (e > 0)
          & (w >= MODEL_RANGE[0]) & (w <= MODEL_RANGE[1]))
    e = np.maximum(e, f / snr_ceiling)
    for lo, hi in (exclude or []):
        ok &= ~((w >= lo) & (w <= hi))

    return Observation(
        name='ngsl', star=star, wavelength=w, flux=f, uncertainty=e, mask=ok,
        resolution=('R', 600.0),        # MEASURED against XSL, not tabulated
        calibration=('scalar',),        # shape trusted, absolute level is (R/d)^2
        rv_fixed=None,
        meta=dict(offset_px=row['offset_px'], dataqual=row['dataqual'],
                  snr_ceiling=snr_ceiling, file=row['ngsl_file']))


def xsl_resolution_segments():
    """XSL resolving power per arm.

    XSL quotes sigma(v), NOT FWHM: 13 km/s UVB, 11 VIS, 16 NIR. FWHM = 2.3548
    sigma, so the UVB is R ~ 9800 -- reading the quoted number as a FWHM gives
    an answer 2.35x wrong.
    """
    return [(lo, hi, C_KMS / (2.3548 * sig)) for _, lo, hi, sig in ARMS]


def load_xsl(star, poly_order=4, exclude=None):
    """XSL DR3 spectrum: vacuum, rest-frame, continuum marginalised.

    Rest-frame means the RV is already removed, so it is FIXED at 0 -- unlike
    NGSL. The continuum gets a marginalised polynomial because XSL is
    ground-based and slit-loss corrected; what survives is the line profiles.
    """
    row = _sample_row(star)
    w, f, e, hdr = xsl_load(row['xslid'])          # nm->A, air->vacuum
    ok = (np.isfinite(w) & np.isfinite(f) & np.isfinite(e) & (f > 0) & (e > 0)
          & (w >= MODEL_RANGE[0]) & (w <= MODEL_RANGE[1]))
    for lo, hi in (exclude or []):
        ok &= ~((w >= lo) & (w <= hi))
    return Observation(
        name='xsl', star=star, wavelength=w, flux=f, uncertainty=e, mask=ok,
        resolution=('R_segments', xsl_resolution_segments()),
        calibration=('poly', poly_order),
        rv_fixed=0.0,
        meta=dict(xslid=row['xslid'], poly_order=poly_order))


# Gaia DR3 integrated photometry. Space-based, so trusted for absolute
# calibration over a wide baseline -- the same reason NGSL is. NOTE gaia_bp
# spans the Balmer break, so including it puts break information into the
# conditioning set; `common.photometry.overlaps` reports which filters do.
GAIA_FILTERS = ('gaia_g', 'gaia_bp', 'gaia_rp')

# Gaia DR3 G/BP/RP are in the VEGAMAG system as published. sedpy works in AB, so
# the published magnitudes are converted with the DR3 passband zero points
# (Gaia DR3 documentation, Table 5.2: G, BP, RP AB offsets).
GAIA_VEGA_TO_AB = dict(gaia_g=0.118, gaia_bp=0.056, gaia_rp=0.379)

# Calibration floor. Gaia's quoted photometric errors are millimag; the real
# limit for stars this bright is the passband calibration, and the whole dust
# argument turns on a colour, so an over-tight error here would be read as a
# reddening measurement it cannot support.
GAIA_MAG_FLOOR = 0.01


def load_photometry(star, names=GAIA_FILTERS, mag_floor=GAIA_MAG_FLOOR):
    """Gaia DR3 G/BP/RP as maggies -> Observation with sedpy filters."""
    from common.photometry import filter_set, mag_to_maggies
    g = None
    for r in csv.DictReader(open(ROOT / 'data' / 'gaia_sample.csv')):
        if r['star'] == star:
            g = r
            break
    if g is None:
        raise KeyError(f'{star} is not in data/gaia_sample.csv')

    key = dict(gaia_g='gmag', gaia_bp='bpmag', gaia_rp='rpmag')
    ekey = dict(gaia_g='e_gmag', gaia_bp='e_bpmag', gaia_rp='e_rpmag')
    mags, errs = [], []
    for n in names:
        m = float(g[key[n]]) + GAIA_VEGA_TO_AB[n]
        try:
            em = float(g[ekey[n]])
        except (KeyError, ValueError):
            em = 0.0
        mags.append(m)
        errs.append(max(em, mag_floor))
    maggies, munc = mag_to_maggies(np.array(mags), np.array(errs))

    return Observation(
        name='phot', star=star, flux=maggies, uncertainty=munc,
        mask=np.isfinite(maggies), filters=filter_set(names),
        resolution=(None,),             # a bandpass integral needs no LSF
        calibration=('scalar',),        # colours only; (R/d)^2 is marginalised
        rv_fixed=0.0,
        meta=dict(names=list(names), vegamag=[float(g[key[n]]) for n in names],
                  mag_floor=mag_floor, source_id=g['source_id']))


# The held-out region: predicted, never conditioned on. 3550-4000 A rather than
# a narrow window around 3646, because the high-order Balmer lines crowd
# together from H-epsilon 3970 down to the series limit, and Ca II H and K sit
# at 3934/3968. Anything in here is Balmer-break information.
BREAK_WINDOW = (3550.0, 4000.0)

# Bands for the Gaia XP spectrum. Edges are placed in line-free continuum, which
# is what makes the band integral independent of XP's (complicated,
# wavelength-dependent, R ~ 20-100) line-spread function: convolution conserves
# the integral, so a line wholly inside a band contributes the same flux however
# it is smeared. A band edge cutting through a Balmer wing would NOT be safe.
#
# 3550-4000 is omitted entirely -- that is the held-out break. 3360-3540 is the
# blue-of-break band and is the one that matters most for dust; it is narrow
# because XP starts at 3360 A and the break region begins at 3550.
# Edges sit 15+ A inside BOTH the XP range (3360-10200 A) and the model grid
# (3200-9500 A), because tophat() tapers ~10 A beyond each edge and
# common.photometry.project returns NaN for a band that is not fully covered.
# At 3360 and 9500 exactly, the blue and red bands both came back NaN.
XP_BANDS = [('xp_3380_3540', 3380., 3540.),   # Balmer continuum, blueward
            ('xp_4050_4550', 4050., 4550.),   # contains H-delta 4102, H-gamma 4341
            ('xp_4550_5200', 4550., 5200.),   # contains H-beta 4861
            ('xp_5200_5900', 5200., 5900.),
            ('xp_5900_6300', 5900., 6300.),
            ('xp_6300_7100', 6300., 7100.),   # contains H-alpha 6563
            ('xp_7100_7900', 7100., 7900.),
            ('xp_7900_9480', 7900., 9480.)]   # contains the Paschen series

# Gaia XP is externally calibrated to ~1-2%. The quoted per-pixel errors are far
# smaller than that over a whole band, so a calibration floor is applied -- the
# dust constraint is a colour, and an over-tight band error would be read as a
# reddening measurement the calibration cannot support.
XP_CAL_FLOOR = 0.01


def load_xp(star, bands=XP_BANDS, cal_floor=XP_CAL_FLOOR,
            exclude=(BREAK_WINDOW,)):
    """Gaia XP sampled spectrum, integrated into bands -> Observation.

    This is the long-baseline dust lever: 3360-9500 A from space, ~0.036 mag of
    differential extinction per 0.01 mag of E(B-V), against 0.0035 mag for
    NGSL's 3200-3500 A window. Gaia's BROADBAND photometry cannot do this job --
    gaia_g and gaia_rp run past the grid's 9500 A limit, and gaia_bp alone is
    exactly determined by its own free scalar, so it carries no information.

    Bands rather than pixels because XP's LSF is a basis reconstruction, not a
    Gaussian; banding makes the comparison LSF-free (see common.photometry.tophat).
    """
    from common.photometry import tophat
    p = ROOT / 'data' / 'gaia_xp' / f'{star}_xp.csv'
    d = np.loadtxt(p, delimiter=',', skiprows=1)
    w, f, e = d[:, 0], d[:, 1], d[:, 2]

    keep = [b for b in bands
            if not any(b[2] > lo and b[1] < hi for lo, hi in (exclude or []))]
    filters = [tophat(*b) for b in keep]

    flux, unc = [], []
    for nm, lo, hi in keep:
        s = (w >= lo) & (w <= hi) & np.isfinite(f)
        # band flux as a maggie-like quantity is not needed: integrate f_lambda
        # through the same tophat the model will see, via the same code path
        from common.photometry import project as phot_project
        val = phot_project(w, f, [tophat(nm, lo, hi)])[0]
        # error: quadrature sum over the band, then the calibration floor
        rel = (np.sqrt(np.sum(e[s] ** 2)) / np.sum(f[s])) if s.sum() else np.nan
        flux.append(val)
        unc.append(abs(val) * max(rel, cal_floor))

    flux, unc = np.array(flux), np.array(unc)
    return Observation(
        name='xp', star=star, flux=flux, uncertainty=unc,
        mask=np.isfinite(flux), filters=filters,
        resolution=(None,),          # a band integral needs no LSF -- the point
        calibration=('scalar',),     # colours only
        rv_fixed=0.0,
        meta=dict(names=[b[0] for b in keep],
                  bands=keep, cal_floor=cal_floor,
                  excluded=[b[0] for b in bands if b not in keep],
                  file=str(p.relative_to(ROOT))))


def load_all(star, **kw):
    """-> [ngsl, xsl, phot] for one star, skipping any that is unavailable."""
    out = []
    for fn, label in ((load_ngsl, 'ngsl'), (load_xsl, 'xsl'),
                      (load_photometry, 'phot')):
        try:
            out.append(fn(star, **kw.get(label, {})))
        except (KeyError, FileNotFoundError, OSError) as exc:
            print(f'  {star}: no {label} ({type(exc).__name__}: {exc})')
    return out
