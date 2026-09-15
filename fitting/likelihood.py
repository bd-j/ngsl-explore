"""chi^2 -> log likelihood: marginalising the linear nuisance terms.

fitting/calibration.py PROFILES the calibration coefficients -- it returns the
best-fit c and stops. That understates the uncertainty on everything the
coefficients are covariant with, which matters most for XSL, whose segmented
continuum carries 12 of them.

Marginalising instead: chi^2 is quadratic in c, so

    chi2(c) = chi2_min + (c - c_hat)^T A (c - c_hat),    A = B^T C^-1 B

and with a flat prior on c the integral is analytic:

    ln L = -0.5 * chi2_min - 0.5 * ln|A| + (k/2) ln(2 pi)

The -0.5 ln|A| term is the Occam factor: a calibration with more freedom, or
one whose basis the data constrains loosely, pays for it.

A TRAP IN THE FLAT PRIOR, which is why the convention is a parameter here and
not a hardcoded choice. For a pure 'scalar' calibration c0 is (R/d)^2 times a
grey factor, and the model's absolute normalisation is arbitrary -- these are
surface fluxes and the distance is unknown. Rescale the model by alpha and
chi2_min is unchanged, but ln|A| shifts by 2 ln alpha, so ln L shifts by
-ln alpha. Since surface flux runs roughly as T^4, a FLAT prior on c0 therefore
imposes a spurious Teff-dependent penalty of about -4 ln(T) across the grid --
an artefact of the parameterisation, not evidence about the star.

    coeff_prior='flat'        flat in c. Standard, and the only closed form
                              when k > 1. Carries the alpha dependence above,
                              so it is safe for comparing nodes only if the
                              caller accepts that tilt or supplies a real prior.
    coeff_prior='scale_free'  flat in ln c0, valid ONLY for k = 1 (a pure
                              scalar). Then the marginal is
                              -ln c0_hat - 0.5 ln A, and since A -> alpha^2 A
                              while c0_hat -> c0_hat / alpha the two shifts
                              cancel EXACTLY. It needs the fitted c0_hat, which
                              the scan stores as `scalar`.

There is deliberately no scale-free option for k > 1. Splitting |A| between an
amplitude and a shape subspace is not a matter of scaling ln|A| by (k-1)/k --
an earlier version of this file did exactly that, which is not the marginal of
anything. Doing it properly means reparameterising c = c0 * u and integrating
the shape block, and XSL's 12 segmented coefficients have no single amplitude
to factor out. Until that is worked through, XSL uses the flat prior and the
tilt is a known, documented systematic rather than a silent fudge.

The scan stores chi^2, ln|A| and the pixel count rather than ln L for the same
reason this module is separate: the error model is still an open question, and
a stored likelihood would freeze a convention that a re-run of the scan (1.7 h)
would be needed to change.
"""
import numpy as np

from fitting.calibration import design_matrix, usable

LN2PI = float(np.log(2.0 * np.pi))


def curvature(obs, model):
    """-> (ln|A|, k) for A = B^T C^-1 B, the calibration's Fisher matrix.

    Returns (nan, 0) when there are no calibration coefficients, and raises
    nothing -- a singular A is reported as -inf so a caller can drop the node.
    """
    m = usable(obs, model)
    B, k = design_matrix(obs, model)
    if k == 0:
        return float('nan'), 0
    ivar = 1.0 / np.asarray(obs.uncertainty, float)[m] ** 2
    Bw = B * np.sqrt(ivar)[:, None]
    sign, logdet = np.linalg.slogdet(Bw.T @ Bw)
    return (float(logdet) if sign > 0 else float('-inf')), int(k)


def scale_profile(chi2, n, floor=1.0):
    """Analytic profile over an error-scale s, with s floored.

    Scaling every uncertainty by s gives ln L = -0.5*chi2/s^2 - n*ln(s), whose
    maximum is at s_hat^2 = chi2/n. So node ranking is by n*ln(chi2/n) rather
    than by chi2 -- which is what removes `lnerr` from the parameter vector.

    The floor exists because inflation may WIDEN an error bar and never shrink
    one: the NGSL bands sit at chi^2/n ~ 0.14, and rescaling that to 1 would
    deflate their errors by 2.7x, claiming exactly the precision the 1%
    calibration floor was set to disclaim.

    -> (s_hat, lnL up to a constant independent of the model)
    """
    chi2 = np.asarray(chi2, float)
    n = int(n)
    if n <= 0:
        return np.nan, np.full(chi2.shape, np.nan) if chi2.ndim else np.nan
    s = np.maximum(float(floor), np.sqrt(chi2 / n))
    return s, -0.5 * chi2 / s ** 2 - n * np.log(s)


def lnlike(chi2, n, lndetA=None, k=0, scale='profile', floor=1.0,
           marginalize=True, coeff_prior='flat', c0=None):
    """chi^2 (+ curvature) -> ln L, up to an additive constant.

    chi2 and lndetA may be arrays of any matching shape, so a whole stored
    conditional curve converts in one call.

    scale   'profile' analytic error-scale profile (see scale_profile)
            'fixed'   trust the quoted uncertainties as they stand
    """
    chi2 = np.asarray(chi2, float)
    if scale == 'profile':
        _, ln = scale_profile(chi2, n, floor)
    elif scale == 'fixed':
        ln = -0.5 * chi2
    else:
        raise ValueError(f'unknown scale: {scale!r}')

    if marginalize and k:
        if lndetA is None:
            raise ValueError('marginalize=True needs lndetA')
        lndetA = np.asarray(lndetA, float)
        if coeff_prior == 'flat':
            ln = ln - 0.5 * lndetA + 0.5 * k * LN2PI
        elif coeff_prior == 'scale_free':
            if k != 1:
                raise ValueError(
                    "coeff_prior='scale_free' is only derived for k = 1; "
                    f'this observation has k = {k}. See the module docstring.')
            if c0 is None:
                raise ValueError("coeff_prior='scale_free' needs c0 "
                                 '(the fitted scalar; the scan stores it)')
            c0 = np.asarray(c0, float)
            with np.errstate(divide='ignore', invalid='ignore'):
                ln = ln - np.log(np.abs(c0)) - 0.5 * lndetA + 0.5 * LN2PI
        else:
            raise ValueError(f'unknown coeff_prior: {coeff_prior!r}')
    return ln


def combine(legs, weights=None):
    """Add ln L from independent legs, each already on its own scale.

    Legs are independent datasets (NGSL bands, XSL lines), so their ln L simply
    add. `weights` is provided only to down-weight a leg deliberately and is
    NOT a way to fix a leg that is miscalibrated -- that is what the error-scale
    profile is for.
    """
    out = None
    for i, l in enumerate(legs):
        w = 1.0 if weights is None else float(weights[i])
        out = w * np.asarray(l, float) if out is None else out + w * np.asarray(l, float)
    return out
