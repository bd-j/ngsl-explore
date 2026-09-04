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

Writes data/ngsl_wavecal.csv
"""
import csv
import sys
from pathlib import Path

import numpy as np
from astropy.io import fits

sys.path.insert(0, str(Path(__file__).resolve().parent))
from make_atlas_model import hnu_to_flam

ROOT = Path(__file__).resolve().parent.parent
LAM_REF = 4000.0            # pivot for the linear term, near the Balmer break
SHIFTS = np.arange(-6.0, 6.01, 0.02)
MIN_CORR = 0.5

# name, lo, hi, polynomial degree, fitting sub-windows
SEGMENTS = [
    ('G230LB', 1675.0, 3058.0, 0, [(2600., 3050.)]),
    ('G430L', 3058.0, 5647.0, 1, [(3300., 3700.), (3700., 4100.), (4100., 4600.),
                                  (4600., 5100.), (5100., 5600.)]),
    ('G750L', 5647.0, 10198.0, 0, [(5800., 6500.), (6500., 7200.), (7200., 7900.),
                                   (7900., 8600.), (8600., 9100.)]),
]


def air_to_vac(w):
    s = 1e4 / np.asarray(w, dtype=float)
    n = 1 + 0.05792105 / (238.0185 - s * s) + 0.00167917 / (57.362 - s * s)
    return np.asarray(w, dtype=float) * n


def _shift(wo, fo, wm, fm, lo, hi):
    """Cross-correlation shift (A) aligning model to observation, or nan.

    Returns nan when the window is not covered by both spectra (the models
    start at 3200 A, so G230LB cannot be calibrated this way), when the peak
    correlation is poor, or when the peak sits at the search boundary.
    """
    lo, hi = max(lo, float(wm[0])), min(hi, float(wm[-1]))
    if hi - lo < 100.0:
        return np.nan, np.nan
    m = (wo > lo) & (wo < hi) & np.isfinite(fo) & (fo > 0)
    if m.sum() < 80:
        return np.nan, np.nan
    x, y = wo[m], fo[m]
    yn = y / np.polyval(np.polyfit(x, y, 3), x)
    cc = np.empty_like(SHIFTS)
    for i, s in enumerate(SHIFTS):
        mi = np.interp(x, wm + s, fm)
        cc[i] = np.corrcoef(yn, mi / np.polyval(np.polyfit(x, mi, 3), x))[0, 1]
    j = int(np.argmax(cc))
    if cc[j] < MIN_CORR or j in (0, len(SHIFTS) - 1):
        return np.nan, float(cc[j])
    return float(SHIFTS[j]), float(cc[j])


def fit_star(obs_file, model_spec, lsf):
    """-> {grating: dict(a, b, n, rms)} after the air->vacuum conversion."""
    d = fits.getdata('data/spectra/' + obs_file)
    wo, fo = air_to_vac(d['WAVELENGTH'].astype(float)), d['FLUX'].astype(float)
    wm, hnu, _ = np.loadtxt(model_spec, unpack=True)
    fm = lsf(wm, hnu_to_flam(wm, hnu))
    out = {}
    for name, _, _, deg, wins in SEGMENTS:
        xs, ys = [], []
        for lo, hi in wins:
            s, _ = _shift(wo, fo, wm, fm, lo, hi)
            if np.isfinite(s):
                xs.append((lo + hi) / 2 - LAM_REF)
                ys.append(s)
        if not xs:
            out[name] = None
            continue
        xs, ys = np.array(xs), np.array(ys)
        if deg == 1 and len(xs) >= 3:
            b, a = np.polyfit(xs, ys, 1)
            rms = float(np.std(ys - (a + b * xs)))
        else:
            a, b = float(np.median(ys)), 0.0      # median: robust to one bad window
            rms = float(np.std(ys - a))
        out[name] = dict(a=float(a), b=float(b), n=len(xs), rms=rms)
    return out


def apply_wavecal(wave, star, table):
    """Air->vacuum, then remove the fitted per-grating residual.

    _shift returns s such that the model at (wm + s) matches the observation,
    so observed features sit s redward of truth and the data are corrected by
    SUBTRACTING the fitted s(lambda).
    """
    w = air_to_vac(wave)
    out = w.copy()
    for name, lo, hi, _, _ in SEGMENTS:
        c = table.get((star, name))
        if c is None:
            continue
        seg = (w >= lo) & (w < hi)
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


if __name__ == '__main__':
    from plot_model_vs_obs import broaden_ngsl
    cat = [r for r in csv.DictReader(open(ROOT / 'data' / 'balmer_candidates.csv'))
           if r.get('selected') != 'no']
    rows = []
    print(f'Fitted wavelength correction (after air->vacuum), pivot {LAM_REF:.0f} A')
    print(f'{"star":<10}{"grating":<8}{"a (A)":>8}{"b (A/1000A)":>13}{"n":>4}{"rms":>7}')
    for r in cat:
        spec = ROOT / 'models' / 'work' / f'{r["star"]}.spec'
        if not spec.exists():
            continue
        for g, c in fit_star(r['file'], spec, broaden_ngsl).items():
            if c is None:
                rows.append(dict(star=r['star'], grating=g, a_A='', b_A_per_A='',
                                 n_windows=0, fit_rms_A=''))
                continue
            rows.append(dict(star=r['star'], grating=g, a_A=round(c['a'], 4),
                             b_A_per_A=round(c['b'], 7), n_windows=c['n'],
                             fit_rms_A=round(c['rms'], 3)))
            print(f'{r["star"]:<10}{g:<8}{c["a"]:>8.2f}{c["b"]*1000:>13.3f}'
                  f'{c["n"]:>4}{c["rms"]:>7.2f}')
    with open(ROOT / 'data' / 'ngsl_wavecal.csv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=['star', 'grating', 'a_A', 'b_A_per_A',
                                           'n_windows', 'fit_rms_A'])
        w.writeheader()
        w.writerows(rows)
    print(f'\n{len(rows)} rows -> data/ngsl_wavecal.csv')
