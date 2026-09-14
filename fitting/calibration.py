"""The linear nuisance terms, solved in closed form.

Every calibration policy in observations.py is LINEAR in its coefficients:

    model_i = M_i * sum_k c_k B_k(lambda_i)

- 'scalar'    is one column, B_0 = 1.  This is the (R/d)^2 of the star.
- 'poly', n   is n+1 Chebyshev columns across the whole fitted range.
- 'segments'  is an independent Chebyshev per named wavelength segment, each
              zero outside its own segment. This is what XSL needs: a Balmer
              window wants its own local continuum, while the metal windows are
              5-35 A wide and cannot support one each, so they share a low-order
              polynomial per arm.
- 'none'      is no columns at all.

So chi^2 is quadratic in c and the best-fit coefficients come from one weighted
least-squares solve -- no sampling, no optimiser. That is the whole reason the
flux normalisation was never a sampled parameter in this project, and the same
algebra now covers XSL's continuum polynomial and the photometric scalar.

Writing it once is the point: NGSL's normalisation, XSL's continuum and the
photometric scale are then provably the same operation, so a convention error
cannot apply to one dataset and not another.

Marginalising rather than profiling these coefficients adds a -0.5*ln|A| term
(A = B^T C^-1 B) and needs a proper prior on c; that belongs with the likelihood.
This module stops at the point estimate, which is what a residual plot needs.
"""
import numpy as np


def segment_pixels(w, seg):
    """Pixels a calibration segment APPLIES to.

    Separate from the polynomial's domain on purpose. The XSL metal windows are
    scattered across a whole arm and share one polynomial, so the domain is the
    arm while the applicable pixels are only the windows -- and they must
    exclude the Balmer windows, which carry their own local continuum. Defining
    a segment by its domain alone made those overlap, giving pixels columns from
    two segments and a rank-deficient design matrix.
    """
    w = np.asarray(w, float)
    sel = np.zeros_like(w, bool)
    for lo, hi in seg['ranges']:
        sel |= (w >= lo) & (w <= hi)
    return sel


def segment_x(w, seg):
    """Chebyshev variable on [-1, 1] over the segment's DOMAIN."""
    lo, hi = seg['domain']
    return 2.0 * (np.asarray(w, float) - lo) / (hi - lo) - 1.0


def check_segments(segs):
    """Raise if two segments claim the same wavelength -- see segment_pixels."""
    flat = [(lo, hi, i) for i, s in enumerate(segs) for lo, hi in s['ranges']]
    flat.sort()
    for (lo1, hi1, i1), (lo2, hi2, i2) in zip(flat, flat[1:]):
        if lo2 < hi1 and i1 != i2:
            raise ValueError(f'calibration segments {i1} and {i2} overlap at '
                             f'{lo2:.1f}-{min(hi1, hi2):.1f} A')


def usable(obs, model):
    """Masked pixels where the MODEL is also finite and positive.

    The observation mask cannot know about the model, and a non-finite model
    entry is not hypothetical: a filter whose bandpass runs past the grid's
    9500 A limit comes back NaN by design (common.photometry.project refuses to
    integrate over missing spectrum), and one NaN row makes the whole
    least-squares solve fail with 'SVD did not converge'. Drop those rows here so
    the failure is a reported exclusion rather than a crash.
    """
    mod = np.asarray(model, float)
    return np.asarray(obs.mask, bool) & np.isfinite(mod) & (mod > 0)


def design_matrix(obs, model):
    """-> (n_usable, n_coeff) matrix whose columns multiply the model.

    Built only on the usable pixels, so a masked region cannot influence the
    continuum solution.
    """
    m = usable(obs, model)
    mod = np.asarray(model, float)[m]
    kind = obs.calibration[0]

    if kind == 'none':
        return mod[:, None], 0
    if kind == 'scalar':
        return mod[:, None], 1
    if kind == 'segments':
        if obs.wavelength is None:
            raise ValueError("'segments' calibration needs a wavelength axis")
        w = np.asarray(obs.wavelength, float)[m]
        cols = []
        for seg in obs.calibration[1]:
            inseg = segment_pixels(w, seg)
            order = int(seg['order'])
            if inseg.sum() <= order:
                continue          # not enough points to define this continuum
            x = segment_x(w, seg)
            for k in range(order + 1):
                c = np.zeros_like(w)
                c[inseg] = mod[inseg] * np.polynomial.chebyshev.chebval(
                    x[inseg], np.eye(order + 1)[k])
                cols.append(c)
        if not cols:
            raise ValueError(f'{obs.star}/{obs.name}: no usable calibration segment')
        return np.vstack(cols).T, len(cols)

    if kind == 'poly':
        order = int(obs.calibration[1])
        if obs.wavelength is None:
            raise ValueError("'poly' calibration needs a wavelength axis")
        w = np.asarray(obs.wavelength, float)[m]
        # Chebyshev on [-1, 1] across the fitted range. Scaling to the range
        # actually used keeps the columns well conditioned; using raw Angstroms
        # would make the high orders numerically degenerate.
        lo, hi = w.min(), w.max()
        x = 2.0 * (w - lo) / (hi - lo) - 1.0
        cols = [mod * np.polynomial.chebyshev.chebval(
            x, np.eye(order + 1)[k]) for k in range(order + 1)]
        return np.vstack(cols).T, order + 1
    raise ValueError(f'unknown calibration kind: {kind!r}')


def solve(obs, model, rcond=None):
    """Best-fit calibration coefficients -> (calibrated_model_full, coeffs).

    The returned model is full length (unmasked pixels included) so it can be
    plotted against the whole spectrum; only masked pixels entered the solve.
    """
    m = usable(obs, model)
    B, ncoeff = design_matrix(obs, model)
    if ncoeff == 0:
        return np.asarray(model, float), np.array([1.0])
    if m.sum() < ncoeff:
        raise ValueError(f'{obs.star}/{obs.name}: {m.sum()} usable points for '
                         f'{ncoeff} calibration coefficients')

    y = np.asarray(obs.flux, float)[m]
    ivar = 1.0 / np.asarray(obs.uncertainty, float)[m] ** 2
    sw = np.sqrt(ivar)
    c, *_ = np.linalg.lstsq(B * sw[:, None], y * sw, rcond=rcond)

    # Re-evaluate the polynomial on ALL pixels, not just the masked ones.
    full = np.asarray(model, float)
    kind = obs.calibration[0]
    if kind == 'scalar':
        return full * c[0], c
    if kind == 'poly':
        w = np.asarray(obs.wavelength, float)
        wm = w[m]
        lo, hi = wm.min(), wm.max()
        x = 2.0 * (w - lo) / (hi - lo) - 1.0
        poly = np.polynomial.chebyshev.chebval(x, c)
        return full * poly, c
    if kind == 'segments':
        w = np.asarray(obs.wavelength, float)
        wm = w[m]
        out = np.full_like(full, np.nan)
        i = 0
        for seg in obs.calibration[1]:
            order = int(seg['order'])
            if segment_pixels(wm, seg).sum() <= order:
                continue
            n = order + 1
            sel = segment_pixels(w, seg)
            x = segment_x(w, seg)
            out[sel] = full[sel] * np.polynomial.chebyshev.chebval(x[sel], c[i:i + n])
            i += n
        return out, c
    raise ValueError(kind)


def residual(obs, calibrated):
    """Fractional residual (obs - model) / model, NaN outside the mask."""
    out = np.full(np.shape(obs.flux), np.nan)
    m = usable(obs, calibrated)
    cal = np.asarray(calibrated, float)
    out[m] = (np.asarray(obs.flux, float)[m] - cal[m]) / cal[m]
    return out


def chi2(obs, calibrated):
    m = usable(obs, calibrated)
    r = (np.asarray(obs.flux, float)[m] - np.asarray(calibrated, float)[m])
    return float(np.sum((r / np.asarray(obs.uncertainty, float)[m]) ** 2))
