"""Per-grating wavelength recalibration of NGSL spectra against the models.

Two corrections, applied in order:

1. AIR -> VACUUM. NGSL v2 wavelengths are in AIR. Established by
   cross-correlating each star against its own ATLAS12 model in seven windows
   from 3300 to 9100 A: as delivered the required shift runs -0.2 to -3.2 A and
   tracks the air-vacuum curve; after converting, the mean shift drops from
   -1.56 to +0.15 A and the maximum from 3.18 to 0.85 A.

   (An earlier test using only 3300-4150 A concluded vacuum. That is the worst
   possible window: the air-vacuum offset there is ~1 A, comparable to the
   residual below, so the two are not separable in it.)

2. PER-GRATING RESIDUAL, LINEAR IN WAVELENGTH where warranted. NGSL took no
   wavecals with the stellar exposures -- the delivery readme states zero
   points were derived per spectrum from the positions of strong stellar
   features -- and the gratings were reduced separately, so each carries its
   own error.

   After conversion G430L shows a clean monotonic ramp (+1.11 A at 3500 A
   falling to +0.08 A at 5350 A) rather than an offset. A constant is therefore
   wrong for it: fitting one at 4200-5600 A, where the residual is ~0, leaves
   ~1 A uncorrected at the Balmer break. Such a linear-in-lambda residual is
   what a DIFFERENT AIR-VACUUM CONVENTION produces -- Edlen (1953/1966) vs
   Ciddor (1996), or different assumed T/P for the air index -- so it is
   modelled as a line, not chased as a zero point.

   G750L scatters window-to-window without a clean trend (and its 7900-8600 A
   window straddles the Paschen limit, where the fit tracks model line
   positions rather than calibration), so it gets a robust constant instead.

APPLYING vs MEASURING. This module applies `data/ngsl_wavecal.csv`;
`explore/ngsl_wavecal_fit.py` measures and writes it. Same split as
`common/lsf.py` applying what `explore/ngsl_lsf.py` measures.
"""
import csv
import sys
import warnings
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common.lines import air_to_vac

ROOT = Path(__file__).resolve().parent.parent
LAM_REF = 4000.0            # pivot for the linear term, near the Balmer break

# name, lo, hi, polynomial degree, fitting sub-windows.
# Shared with explore/ngsl_wavecal_fit.py, which fits into this shape.
SEGMENTS = [
    ('G230LB', 1675.0, 3058.0, 0, [(2600., 3050.)]),
    ('G430L', 3058.0, 5647.0, 1, [(3300., 3700.), (3700., 4100.), (4100., 4600.),
                                  (4600., 5100.), (5100., 5600.)]),
    ('G750L', 5647.0, 10198.0, 0, [(5800., 6500.), (6500., 7200.), (7200., 7900.),
                                   (7900., 8600.), (8600., 9100.)]),
]


class MissingWavecal(UserWarning):
    """A star has no fitted correction, so it gets air->vacuum only."""


# G230LB (1675-3058 A) cannot be calibrated this way and never will be: the
# model grid starts at 3200 A, so there is nothing to cross-correlate against.
# It is left uncorrected deliberately rather than given a fitted number, so a
# missing row for it is the expected state and must not warn.
UNCALIBRATABLE = {'G230LB'}

_WARNED = set()


def apply_wavecal(wave, star, table, quiet=False):
    """Air->vacuum, then remove the fitted per-grating residual.

    The fit returns s such that the model at (wm + s) matches the observation,
    so observed features sit s redward of truth and the data are corrected by
    SUBTRACTING the fitted s(lambda).

    A star with no entry gets air->vacuum and NOTHING ELSE, which is a silent
    ~0.8 A error at the Balmer break. This used to pass unnoticed: the table was
    fitted for four stars of the superseded sample and never regenerated, so
    nine of the thirteen sample stars were uncorrected for weeks while every
    figure and residual looked plausible. It now warns once per star and
    grating. Pass quiet=True only when the caller has already established that
    an uncorrected star is acceptable.
    """
    w = air_to_vac(wave)
    out = w.copy()
    for name, lo, hi, _, _ in SEGMENTS:
        seg = (w >= lo) & (w < hi)
        c = table.get((star, name))
        if c is None:
            if (seg.any() and not quiet and name not in UNCALIBRATABLE
                    and (star, name) not in _WARNED):
                _WARNED.add((star, name))
                warnings.warn(
                    f'no wavelength calibration for {star} {name}: '
                    f'{int(seg.sum())} pixels get air->vacuum only. '
                    f'Refit with explore/ngsl_wavecal_fit.py --star {star}',
                    MissingWavecal, stacklevel=2)
            continue
        out[seg] = w[seg] - (c['a'] + c['b'] * (w[seg] - LAM_REF))
    return out


def load_table(path=None):
    path = Path(path or ROOT / 'data' / 'ngsl_wavecal.csv')
    if not path.exists():
        return {}
    rd = csv.DictReader(open(path))
    if not rd.fieldnames or 'a_A' not in rd.fieldnames:
        return {}                      # stale/foreign format: ignore it
    return {(r['star'], r['grating']):
            dict(a=float(r['a_A']), b=float(r['b_A_per_A']))
            for r in rd if r['a_A'] != ''}
