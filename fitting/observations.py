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
from common.lines import hydrogen_lines
from common.lsf import NGSL_MOFFAT_BETA

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
        # MEASURED against XSL -- shape as well as width. The single Gaussian
        # at R=600 was splitting the difference between a sharp core and a
        # heavy tail it could not represent; see common/lsf.py.
        resolution=('ngsl_moffat', NGSL_MOFFAT_BETA),
        calibration=('scalar',),        # shape trusted, absolute level is (R/d)^2
        rv_fixed=None,
        meta=dict(offset_px=row['offset_px'], dataqual=row['dataqual'],
                  snr_ceiling=snr_ceiling, file=row['ngsl_file']))


# XSL sits redward of the models by a small, measured amount
# (explore/xsl_line_offsets.py -> data/xsl_line_offsets.csv). Decomposed over
# 12 stars x 6 isolated lines the offset is:
#
#   grand mean      +4.21 km/s   common to every star and every line
#   star-to-star     2.00 km/s
#   line-to-line     1.13 km/s
#   unexplained      2.30 km/s
#
# The CONSTANT dominates, so it is applied as a zero point to every star and a
# per-star departure is allowed on top -- clipped, because the star-to-star
# term is only 2 km/s and the per-star means are measured from as few as one
# usable line. An unclipped per-star value would be fitting noise: HD164967
# (RUWE = 8.32, an astrometric binary) comes out at +10.1 km/s from a single
# line, which is not a radial velocity measurement.
#
# Sign: the measurement is obs - model, so the MODEL is shifted redward by this
# amount to meet the data, which is what a positive `rv` does in predict().
# For scale, 4.2 km/s is 0.4 of an XSL pixel and a seventh of its resolution
# element, so this matters for line-centre residuals and not for line depths.
XSL_RV_ZEROPOINT = 4.21
XSL_RV_STAR_MAX = 2.0


def xsl_rv(star, path=None, zeropoint=XSL_RV_ZEROPOINT,
           star_max=XSL_RV_STAR_MAX):
    """-> velocity in km/s to apply to the MODEL when comparing with XSL."""
    p = Path(path or ROOT / 'data' / 'xsl_line_offsets.csv')
    if not p.exists():
        return zeropoint
    vals, allv = [], []
    for r in csv.DictReader(open(p)):
        if r['blended'] == 'yes':
            continue                   # a blend measures a line ratio, not a shift
        v = float(r['offset_kms'])
        allv.append(v)
        if r['star'] == star:
            vals.append(v)
    if not vals or not allv:
        return zeropoint
    delta = float(np.mean(vals) - np.mean(allv))
    return zeropoint + float(np.clip(delta, -star_max, star_max))


def xsl_resolution_segments():
    """XSL resolving power per arm.

    XSL quotes sigma(v), NOT FWHM: 13 km/s UVB, 11 VIS, 16 NIR. FWHM = 2.3548
    sigma, so the UVB is R ~ 9800 -- reading the quoted number as a FWHM gives
    an answer 2.35x wrong.
    """
    return [(lo, hi, C_KMS / (2.3548 * sig)) for _, lo, hi, sig in ARMS]


# --- XSL fit regions ------------------------------------------------------
#
# XSL is fitted only inside named windows, not across its whole range. Two kinds:
#
#   BALMER   H-alpha, H-beta, H-gamma, H-delta, each +/- XSL_BALMER_HALFWIDTH
#            with the CORE masked. These are the dust-immune Teff / log g
#            diagnostic -- a locally normalised profile cannot be changed by a
#            smooth reddening law -- and they are why XSL is in the analysis.
#
#            H-epsilon and higher orders are excluded. They blend into one
#            another, so a local continuum is not defined for them: measured at
#            XSL resolution the wing of H8 does not return to within 2% of the
#            continuum until 122 A from centre, against 37-41 A for these four.
#            They also sit inside the held-out break window.
#
#   METAL    from data/xsl_metal_windows.csv, chosen by measured [M/H]
#            sensitivity (explore/metal_sensitivity.py) rather than by
#            reputation. That measurement is why Mg I b is NOT used: at
#            ~10,000 K magnesium is largely ionised, so Mg I b 5167 changes by
#            -0.090 in depth against Mg II 4481 at -0.121, and both trail the
#            Fe II blends. Ca II H and K are excluded despite being the most
#            sensitive features in the optical: they sit in the held-out window
#            and carry an interstellar component, so they would bias [M/H] in
#            the same direction as the reddening and look self-consistent.
#
# The windows are a UNION over 9000 / 10000 / 11000 K because the ranking is
# NOT stable in Teff -- at 9000 K only 11-17 of the reference top 30 survive,
# with Spearman -0.04 to 0.22 -- while it is stable in log g. The sample spans
# 8759-10885 K, so one temperature's optimum is wrong at the ends of it. The
# union costs little (a window where the line is weak simply contributes little)
# whereas omitting one loses a star's metallicity constraint outright.

# Core masked because the observed Balmer cores carry a flux excess of ~10% of
# the line EW relative to these LTE models -- almost certainly NLTE in hydrogen,
# which the code does not treat for H. Measured core half-width (50% depth) runs
# 0.8 A at H-alpha to 4.7 A at H-delta, so 6 A covers it with margin. For a fast
# rotator this should grow: v sin i = 200 km/s adds 2.9 A at H-gamma.
XSL_BALMER_HALFWIDTH = 50.0     # wing merges into continuum by 37-41 A
XSL_CORE_MASK = 6.0
XSL_BALMER_ORDER = 1            # local continuum per Balmer window
XSL_METAL_ORDER = 3             # one polynomial per arm across the metal windows
XSL_ARM_SPLIT = 5600.0          # UVB / VIS


# Regions dropped from XSL fitting no matter which window contains them.
#
# A metal-window edge is not enough on its own: the 4127-4137 window lies wholly
# inside H-delta's +/-50 A window, so narrowing it changed nothing -- those
# pixels are masked in through H-delta's red wing regardless. Anything the models
# get wrong has to be excluded explicitly.
MODEL_BAD_REGIONS = [
    (4401.0, 4407.0, 'predicted (K13) O I 3p 5P -> 16s 5S; the level is '
                     'dissolved at photospheric density, see CAVEATS.md'),
    (4824.0, 4831.0, 'predicted (K13) O I 3p 3P -> 15d 3D and -> 16s 3S, same'),
    (4120.0, 4127.0, 'Fe II 4123.8 / 4125.9 poorly predicted: residual reaches '
                     '3.7-5.6% against ~1.5% across the rest of the window'),
    (4137.5, 4143.0, 'red edge of the Si II window, residual 2.8%'),
]


def subtract_intervals(interval, blocked, min_width=1.0):
    """(lo, hi) minus a list of blocked ranges -> the surviving pieces."""
    pieces = [list(interval)]
    for blo, bhi in sorted(blocked):
        out = []
        for lo, hi in pieces:
            if bhi <= lo or blo >= hi:
                out.append([lo, hi])
                continue
            if blo > lo:
                out.append([lo, min(blo, hi)])
            if bhi < hi:
                out.append([max(bhi, lo), hi])
        pieces = out
    return [(lo, hi) for lo, hi in pieces if hi - lo >= min_width]


# Metal windows actually FITTED, by feature centre. The sensitivity ranking says
# which features carry [M/H]; this says which of them the models can be trusted
# to reproduce, which is a different question and has to be answered against the
# data. Everything else stays available for prediction plots -- a feature the
# models get wrong is worth LOOKING at and must not be allowed to drive the fit.
XSL_METAL_KEEP = {
    4550.7: 'Fe II',
    4535.3: 'Ti II + Fe II',
    4134.2: 'Fe II + Si II 4129',
}

# The auto-generated window for the 4134 feature (4126.6-4141.6) is centred
# redward of the lines that actually matter, so it is overridden. In vacuum the
# interesting lines are Si II 4129.22 and 4132.06 -- the Si II 4128/4131 doublet
# in the usual AIR naming, a signature of late-B/early-A stars -- plus Fe II
# 4129.90 and Fe I 4133.22. Shifting blueward centres the window on those four.
#
# Trimmed to 4127-4137 after looking at the residual bin by bin. The model is
# poorly predictive at both edges: 4123-4126 reaches 3.7-5.6% (Fe II 4123.82 and
# 4125.95) and 4138 reaches 2.8%, while 4127-4137 holds all four target lines and
# the inter-line continuum stays within ~1.5%.
XSL_METAL_WINDOW = {4134.2: (4127.0, 4137.0)}
#
# NOTE this window lies wholly inside H-delta's +/-50 A window, so those pixels
# are fitted under H-DELTA's local continuum rather than the metal polynomial.
# That is the right way round: the model supplies the Stark wing and the
# polynomial is only a slow correction to it, whereas giving these lines their
# own continuum would have it fight the wing shape.

# Rejected, with the reason, so the choice is reviewable rather than implicit.
XSL_METAL_REJECT = {
    4410.1: 'model over-absorbs at 4403.4 by 21.6%: three PREDICTED (K13) O I '
            'lines 3p 5P -> 16s 5S, upper level 503 cm^-1 below the O I limit',
    4827.3: 'same pathology: predicted (K13) O I 3p 3P -> 15d 3D and -> 16s 3S, '
            'upper levels 488-496 cm^-1 below the limit; the Cr II line that '
            'leads the blend is real but cannot be separated from them',
    4287.6: 'Ti II, deeper than any grid [M/H] can produce, which conflicts '
            'with the Fe II windows -- unresolved, so it cannot carry weight',
    4390.9: 'not needed once 4535/4550 are in; revisit if [M/H] is unconstrained',
    4183.0: 'as 4390.9',
}


def xsl_metal_windows(path=None, which='selected'):
    """-> [(lo, hi)] metal windows.

    which='selected' keeps only XSL_METAL_KEEP; 'all' returns every window from
    data/xsl_metal_windows.csv (explore/metal_sensitivity.py), which is what the
    prediction plots want.
    """
    p = Path(path or ROOT / 'data' / 'xsl_metal_windows.csv')
    if not p.exists():
        return []
    wins = [(float(r['lo']), float(r['hi'])) for r in csv.DictReader(open(p))]
    if which == 'all':
        return wins
    if which != 'selected':
        raise ValueError(f"which must be 'selected' or 'all', got {which!r}")
    out = []
    for lo, hi in wins:
        for c in XSL_METAL_KEEP:
            if lo <= c <= hi:
                out.append(XSL_METAL_WINDOW.get(c, (lo, hi)))
                break
    return sorted(out)


def xsl_fit_windows(core_mask=XSL_CORE_MASK, half_width=XSL_BALMER_HALFWIDTH,
                    metals='selected'):
    """-> (balmer_windows, metal_windows), each a list of (lo, hi)."""
    from common.lines import BALMER, line_windows
    bal = line_windows(sorted(BALMER.values()), half_width, core_mask)
    met = xsl_metal_windows(which=metals) if metals else []
    return bal, met


def load_xsl(star, exclude=None, core_mask=XSL_CORE_MASK,
             half_width=XSL_BALMER_HALFWIDTH, metals='selected',
             drop_bad=True):
    """XSL DR3 spectrum: vacuum, rest-frame, fitted only in named windows.

    Rest-frame means the RV is already removed, so it is FIXED at 0 -- unlike
    NGSL. The continuum is marginalised per segment, which is what removes
    continuum shape and leaves the line profiles.
    """
    row = _sample_row(star)
    w, f, e, hdr = xsl_load(row['xslid'])          # nm->A, air->vacuum
    ok = (np.isfinite(w) & np.isfinite(f) & np.isfinite(e) & (f > 0) & (e > 0)
          & (w >= MODEL_RANGE[0]) & (w <= MODEL_RANGE[1]))

    bal, met = xsl_fit_windows(core_mask, half_width, metals)
    inwin = np.zeros_like(w, bool)
    for lo, hi in bal + met:
        inwin |= (w >= lo) & (w <= hi)
    ok &= inwin
    for lo, hi in (exclude or []):
        ok &= ~((w >= lo) & (w <= hi))

    # Calibration segments. Each states the pixels it APPLIES to separately
    # from its polynomial DOMAIN: a Balmer line gets a local continuum over its
    # own +/-half_width, while the metal windows are 5-35 A wide, cannot each
    # support a continuum, and so share one polynomial per arm evaluated across
    # the whole arm. Applicability and domain must be separate or the arm
    # polynomial would also claim the Balmer pixels and the design matrix goes
    # rank-deficient.
    from common.lines import BALMER
    from fitting.calibration import check_segments
    segs = []
    for lam in sorted(BALMER.values()):
        segs.append(dict(name=f'balmer_{lam:.0f}', order=XSL_BALMER_ORDER,
                         domain=(lam - half_width, lam + half_width),
                         ranges=[(lam - half_width, lam - core_mask),
                                 (lam + core_mask, lam + half_width)]))
    balmer_span = [(lam - half_width, lam + half_width)
                   for lam in BALMER.values()]
    for arm, lo, hi in (('UVB', 0.0, XSL_ARM_SPLIT),
                        ('VIS', XSL_ARM_SPLIT, np.inf)):
        # Metal windows in this arm with the Balmer windows SUBTRACTED, not
        # dropped. Discarding any window that touched a Balmer span threw away
        # the whole 4382.5-4421.8 window -- which contains the single most
        # [M/H]-sensitive feature in the spectrum -- because it clipped the
        # H-gamma window by 9 A at one end.
        rng = [x for l, h in met if lo <= 0.5 * (l + h) < hi
               for x in subtract_intervals((l, h), balmer_span)]
        # `rng` counts RANGES, not pixels -- comparing it against the polynomial
        # order dropped the whole metal segment whenever few windows survived.
        # calibration.design_matrix already skips a segment with too few pixels.
        if rng:
            segs.append(dict(name=f'metal_{arm}', order=XSL_METAL_ORDER,
                             domain=(min(l for l, _ in rng),
                                     max(h for _, h in rng)),
                             ranges=rng))
    check_segments(segs)

    # the mask must match what the segments actually cover
    covered = np.zeros_like(w, bool)
    for seg in segs:
        for lo_, hi_ in seg['ranges']:
            covered |= (w >= lo_) & (w <= hi_)
    ok &= covered

    # ... minus the regions the models are known to get wrong. Gated, because
    # those regions must not drive a FIT but are exactly what a prediction plot
    # exists to show: drop_bad=False keeps them visible.
    if drop_bad:
        for lo_, hi_, _why in MODEL_BAD_REGIONS:
            ok &= ~((w >= lo_) & (w <= hi_))

    return Observation(
        name='xsl', star=star, wavelength=w, flux=f, uncertainty=e, mask=ok,
        resolution=('R_segments', xsl_resolution_segments()),
        calibration=('segments', segs),
        rv_fixed=xsl_rv(star),
        meta=dict(xslid=row['xslid'], n_balmer=len(bal), n_metal=len(met),
                  core_mask=core_mask, half_width=half_width,
                  bad_regions=[(a_, b_) for a_, b_, _ in MODEL_BAD_REGIONS],
                  rv_kms=xsl_rv(star),
                  segments=[(s['name'], s['order'], len(s['ranges']))
                            for s in segs]))


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

# --- NGSL as synthetic photometry ----------------------------------------
#
# NGSL is used TWICE in different roles, and the roles must not overlap:
#
#   load_ngsl_bands()  the continuum, collapsed into synthetic bands. This is
#                      the dust and continuum-shape constraint.
#   the break region   held out entirely, never conditioned on, predicted.
#
# and the hydrogen lines are used from NEITHER -- XSL resolves them ~16x better
# for the same stars, and the NGSL cores carry the known NLTE excess.
#
# So every NGSL pixel enters the analysis at most once. That is what makes
# banded NGSL photometry legitimate where Gaia XP was not: earlier the plan was
# to take the dust constraint from XP, but XP and NGSL disagree in COLOUR by up
# to 4.8% with only 1.2% star-to-star scatter -- an instrumental difference
# worth 0.040 mag in E(B-V), 8x the error budget (explore/xp_vs_ngsl.py).
# Using NGSL for both roles avoids having to decide which instrument is right.
#
# Band edges need no line-free placement here, unlike the XP bands: the data is
# already at R=600 and the model is broadened to R=600, so both sides carry the
# same LSF and there is no leakage mismatch to dodge. The only requirements are
# to avoid the held-out break and the hydrogen lines.
NGSL_H_MASK_A = 20.0        # half-width dropped around every H line
NGSL_BAND_WIDTH = 400.0     # long stretches are split into bands this wide
NGSL_BAND_MIN = 80.0        # a band narrower than this is not worth carrying
NGSL_BLUE_WIDTH = 165.0     # finer bands blueward of the break (see below)


# Bands run from just inside the grid's blue edge to just BLUEWARD OF THE
# PASCHEN BREAK (8206 A = 911.7635 x 9), which keeps 89.5% of the 3220-9480
# dust lever arm (0.0348 vs 0.0389 mag per 0.01 mag of E(B-V)) because CCM89 is
# nearly flat redward of 8000 A. The 10.5% that is given up buys the entire
# Paschen region back as a SECOND untouched prediction, on the same footing as
# the Balmer break.
#
# Inset from MODEL_RANGE at the blue end because tophat() tapers ~10 A past each
# edge and project() refuses a partially covered band: a band starting at
# exactly 3200 is silently dropped for a spectrum starting at 3201, which cost
# the bluest and most important band the first time.
PASCHEN_LIMIT = 8205.9
NGSL_BAND_RANGE = (3220.0, 8180.0)

# Held out alongside the Balmer break, and predicted rather than fitted.
PASCHEN_WINDOW = (8180.0, 9500.0)


def ngsl_band_edges(h_mask=NGSL_H_MASK_A, break_window=None,
                    width=NGSL_BAND_WIDTH, min_width=NGSL_BAND_MIN,
                    wrange=NGSL_BAND_RANGE, blue_width=NGSL_BLUE_WIDTH):
    """-> [(name, lo, hi)] bands covering `wrange` minus H lines and the break.

    Deterministic: the hydrogen line positions are analytic (Rydberg), so this
    needs no model spectrum and cannot drift with the grid.

    Anything blueward of the break gets `blue_width` instead of `width`. That
    stretch -- 3220-3550 A, the Balmer continuum below the series limit -- is
    where the extinction curve is steepest and where the dust lever therefore
    lives, so it is worth resolving as a SHAPE rather than collapsing to a
    single point. One band there gives one colour against the red; two give an
    internal colour across the steepest part of CCM89 as well.
    """
    if break_window is None:
        break_window = BREAK_WINDOW
    lines = hydrogen_lines(wrange[0] - 200.0, wrange[1] + 200.0, series=(2, 3))
    blocked = [(lam - h_mask, lam + h_mask) for lam in lines] + [tuple(break_window)]
    blocked.sort()

    # complement of the blocked intervals within wrange
    free, cur = [], wrange[0]
    for lo, hi in blocked:
        if hi <= cur:
            continue
        if lo > cur:
            free.append((cur, min(lo, wrange[1])))
        cur = max(cur, hi)
        if cur >= wrange[1]:
            break
    if cur < wrange[1]:
        free.append((cur, wrange[1]))

    out = []
    for lo, hi in free:
        if hi - lo < min_width:
            continue
        w = blue_width if hi <= break_window[0] else width
        n = max(1, int(round((hi - lo) / w)))
        edges = np.linspace(lo, hi, n + 1)
        for a, b in zip(edges[:-1], edges[1:]):
            if b - a >= min_width:
                out.append((f'ngsl_{a:.0f}_{b:.0f}', float(a), float(b)))
    return out


def load_ngsl_bands(star, bands=None, cal_floor=0.01, **kw):
    """NGSL collapsed into synthetic bands -> the continuum/dust constraint.

    The 3200-3540 A band matters most: it is blueward of the Balmer series
    limit, so it holds no hydrogen lines, and it is the longest lever NGSL has
    on reddening. It sits ~2.8% below the SYNTHE continuum from metal
    blanketing, which is part of the continuum shape under test rather than a
    reason to exclude it.

    `cal_floor` is NGSL's spectrophotometric accuracy (~1-3%), not its photon
    noise. Band-integrating thousands of pixels drives the statistical error far
    below the calibration error, and the dust constraint is a colour, so an
    over-tight band error would be read as a reddening measurement the
    calibration cannot support.
    """
    from common.photometry import tophat, project as phot_project
    ng = load_ngsl(star, **kw)
    if bands is None:
        bands = ngsl_band_edges()
    m = ng.mask
    wl, fl, er = ng.wavelength[m], ng.flux[m], ng.uncertainty[m]

    keep, flux, unc = [], [], []
    for nm, lo, hi in bands:
        if wl[0] > lo or wl[-1] < hi:
            continue                      # band not covered by this spectrum
        val = phot_project(wl, fl, [tophat(nm, lo, hi)])[0]
        if not np.isfinite(val):
            continue
        s = (wl >= lo) & (wl <= hi)
        rel = np.sqrt(np.sum(er[s] ** 2)) / np.sum(fl[s]) if s.sum() else np.nan
        keep.append((nm, lo, hi))
        flux.append(val)
        unc.append(abs(val) * max(rel, cal_floor))

    flux, unc = np.array(flux), np.array(unc)
    return Observation(
        name='ngsl_bands', star=star, flux=flux, uncertainty=unc,
        mask=np.isfinite(flux), filters=[tophat(*b) for b in keep],
        resolution=('ngsl_moffat', NGSL_MOFFAT_BETA),
        calibration=('scalar',),
        rv_fixed=None,
        meta=dict(names=[b[0] for b in keep], bands=keep,
                  cal_floor=cal_floor, h_mask=NGSL_H_MASK_A,
                  break_window=tuple(BREAK_WINDOW)))


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


def conditioning_set(star, xsl_kw=None, band_kw=None):
    """The observations a fit may condition on: banded NGSL + XSL lines.

    Deliberately NOT the NGSL spectrum itself. NGSL enters once, as bands; its
    hydrogen lines are left to XSL, which resolves them ~16x better for these
    same stars; and the break region is held out by `heldout`. Returning the
    legal set from one place is what keeps that separation enforceable rather
    than conventional -- double-counting NGSL would quietly shrink the
    uncertainty on exactly the parameter the break prediction is most sensitive
    to.
    """
    out = [load_ngsl_bands(star, **(band_kw or {}))]
    try:
        out.append(load_xsl(star, **(xsl_kw or {})))
    except (FileNotFoundError, KeyError) as exc:
        print(f'  {star}: no XSL ({type(exc).__name__}) -- line constraint absent')
    return out


def heldout(star, window=None, **kw):
    """The NGSL spectrum inside the break window: predicted, never fitted."""
    window = tuple(window or BREAK_WINDOW)
    o = load_ngsl(star, **kw)
    o.mask &= (o.wavelength >= window[0]) & (o.wavelength <= window[1])
    o.meta['role'] = 'held out'
    o.meta['window'] = window
    return o
