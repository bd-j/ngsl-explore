"""Filter projection, so model and observed photometry are computed identically.

Thin wrapper on sedpy. The point of routing both sides through one function is
that a filter curve, a zero point or an integration convention can then only be
wrong for both at once -- which shows up as a constant offset rather than as a
spurious colour, and a colour is exactly what the dust term reads.

Units: f_lambda in erg/s/cm^2/A throughout, wavelengths in vacuum Angstroms.
Maggies rather than magnitudes are the working quantity -- they are linear in
flux, so an uncertainty stays symmetric and a zero flux is representable.
"""
import numpy as np

# AB zero point: m_AB = -2.5 log10(f_nu / 3631 Jy), and a maggie is
# 10**(-0.4 m_AB), i.e. f_nu / 3631 Jy.
AB_ZP_JY = 3631.0


def filter_set(names):
    """-> sedpy FilterSet for a list of filter names (see sedpy's data/filters)."""
    from sedpy.observate import load_filters
    return load_filters(list(names))


def project(wave_A, flam, filters):
    """Synthetic photometry of a spectrum -> maggies, one per filter.

    Returns NaN for any filter whose transmission is not fully covered by the
    spectrum, rather than silently integrating over the missing part. That guard
    matters here: the model grid stops at 9500 A, so gaia_g and gaia_rp are only
    partly covered and would otherwise come back quietly too faint.
    """
    from sedpy.observate import getSED
    wave_A = np.asarray(wave_A, float)
    out = np.full(len(filters), np.nan)
    mags = getSED(wave_A, np.asarray(flam, float), filterlist=filters)
    mags = np.atleast_1d(mags)
    for i, f in enumerate(filters):
        lo, hi = f.wavelength[f.transmission > 1e-3 * f.transmission.max()][[0, -1]]
        if wave_A[0] <= lo and wave_A[-1] >= hi:
            out[i] = 10.0 ** (-0.4 * mags[i])
    return out


def coverage(wave_A, filters):
    """-> list of (name, covered, lo, hi) so an uncovered filter is visible."""
    wave_A = np.asarray(wave_A, float)
    rows = []
    for f in filters:
        t = f.transmission > 1e-3 * f.transmission.max()
        lo, hi = f.wavelength[t][[0, -1]]
        rows.append((f.name, bool(wave_A[0] <= lo and wave_A[-1] >= hi),
                     float(lo), float(hi)))
    return rows


def mag_to_maggies(mag, mag_err=None):
    """AB magnitude -> (maggies, maggies_err). Errors propagate as the
    logarithmic derivative, which is exact for a symmetric magnitude error."""
    m = 10.0 ** (-0.4 * np.asarray(mag, float))
    if mag_err is None:
        return m, None
    return m, m * 0.4 * np.log(10.0) * np.asarray(mag_err, float)


def overlaps(filters, window):
    """Filters whose bandpass touches `window` -- used to keep a held-out
    region out of the conditioning set. The Balmer break cannot be predicted
    from data that was used to constrain the model."""
    lo, hi = window
    out = []
    for f in filters:
        t = f.transmission > 1e-3 * f.transmission.max()
        a, b = f.wavelength[t][[0, -1]]
        if b >= lo and a <= hi:
            out.append(f.name)
    return out
