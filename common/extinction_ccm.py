"""Cardelli, Clayton & Mathis (1989) extinction law, optical/NIR only.

Implemented here rather than pulled in as a dependency: only the optical-NIR
branch (1.1 <= x <= 3.3 /um, i.e. 3030-9090 A) is needed, plus the IR branch
below it. The UV branch (x > 3.3) is NOT implemented and raises, so it cannot
be used silently outside its validity range.

CCM89 eq. 1: A(lambda)/A(V) = a(x) + b(x)/R_V, with A(V) = R_V * E(B-V).
"""
import numpy as np

R_V_DEFAULT = 3.1


def ccm89_alav(wave_A, r_v=R_V_DEFAULT):
    """A(lambda)/A(V) for wavelengths in Angstroms."""
    x = 1e4 / np.asarray(wave_A, dtype=float)          # inverse microns
    a = np.zeros_like(x)
    b = np.zeros_like(x)

    ir = (x >= 0.3) & (x < 1.1)
    a[ir] = 0.574 * x[ir] ** 1.61
    b[ir] = -0.527 * x[ir] ** 1.61

    opt = (x >= 1.1) & (x <= 3.3)
    y = x[opt] - 1.82
    a[opt] = (1 + 0.17699 * y - 0.50447 * y**2 - 0.02427 * y**3
              + 0.72085 * y**4 + 0.01979 * y**5 - 0.77530 * y**6
              + 0.32999 * y**7)
    b[opt] = (1.41338 * y + 2.28305 * y**2 + 1.07233 * y**3
              - 5.38434 * y**4 - 0.62251 * y**5 + 5.30260 * y**6
              - 2.09002 * y**7)

    bad = (x < 0.3) | (x > 3.3)
    if np.any(bad):
        raise ValueError(
            f'CCM89 here covers 3030-33000 A only; got '
            f'{np.min(wave_A[bad]):.0f}-{np.max(wave_A[bad]):.0f} A. '
            'The UV branch is deliberately not implemented.')
    return a + b / r_v


def deredden(wave_A, flux, ebv, r_v=R_V_DEFAULT):
    """Remove extinction: returns the intrinsic (brighter) flux."""
    a_v = r_v * ebv
    a_lam = a_v * ccm89_alav(wave_A, r_v)
    return flux * 10.0 ** (0.4 * a_lam)


def redden(wave_A, flux, ebv, r_v=R_V_DEFAULT):
    """Apply extinction to a model."""
    a_lam = r_v * ebv * ccm89_alav(wave_A, r_v)
    return flux * 10.0 ** (-0.4 * a_lam)
