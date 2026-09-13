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

    Returns NaN for any filter whose bandpass is not FULLY covered by the
    spectrum. This guard is not redundant with sedpy: `obj_counts_hires` has its
    "Source spectrum does not span filter" assertion commented out, so it clips
    the integration to whatever is covered and returns a perfectly plausible
    number. Measured: a source truncated at 9500 A gives gaia_rp 0.040 mag too
    faint, silently. The model grid stops at 9500 A and the entire dust argument
    is a colour, so a 0.04 mag error in one band is exactly the failure that
    must not pass quietly.
    """
    from sedpy.observate import getSED
    wave_A = np.asarray(wave_A, float)
    out = np.full(len(filters), np.nan)
    maggies = np.atleast_1d(getSED(wave_A, np.asarray(flam, float),
                                   filterlist=filters, linear_flux=True))
    for i, f in enumerate(filters):
        lo, hi = f.wavelength[f.transmission > 1e-3 * f.transmission.max()][[0, -1]]
        if wave_A[0] <= lo and wave_A[-1] >= hi:
            out[i] = maggies[i]
    return out


def tophat(name, lo, hi, taper=1.0, step=1.0):
    """A rectangular bandpass as a sedpy Filter, for banding up a spectrum.

    Integrating a spectrum over a band is how the Gaia XP spectra are used
    here without needing XP's line-spread function. Convolution conserves the
    integral, so a band whose EDGES sit in smooth continuum gives the same
    answer whatever the LSF is -- even if the band contains a line, as long as
    the line and its wings are wholly inside. That is why the band edges below
    are placed in line-free stretches rather than at round numbers.

    `taper` softens the edge by a few Angstroms so the transmission is not a
    literal step; a hard edge makes the band integral sensitive to how much
    flux the LSF moves across it.
    """
    from sedpy.observate import Filter
    w = np.arange(lo - 10 * taper, hi + 10 * taper, step)
    t = (0.5 * (1 + np.tanh((w - lo) / taper))
         * 0.5 * (1 + np.tanh((hi - w) / taper)))
    return Filter(kname=name, data=(w, t))


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
