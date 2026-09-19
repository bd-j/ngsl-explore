"""Line positions and window definitions, in one place.

`hydrogen_lines` previously existed TWICE -- in `fitting/fit.py` and again in
`explore/plot_ngsl_vs_model.py` -- which is the same drift that once let the
fitter and the comparison figures disagree about the NGSL LSF. One definition,
imported everywhere.

All wavelengths are VACUUM Angstroms, matching ATLAS12/SYNTHE output and the
converted observations.
"""
import numpy as np

RYDBERG_A = 911.7635        # the Lyman limit; series limit of m is m^2 times it

BALMER_LIMIT = 4 * RYDBERG_A        # 3646.1
PASCHEN_LIMIT = 9 * RYDBERG_A       # 8205.9


def air_to_vac(w):
    """Air -> vacuum Angstroms, Ciddor (1996) via the IAU standard inverse.

    Good to <1e-4 A over this project's range. Every library here is delivered
    in AIR (NGSL, XSL, UVES-POP) while ATLAS12/SYNTHE output is VACUUM, so this
    runs on the way in from all three. It previously existed THREE times, once
    per loader -- identical, but three places to fix and three to get wrong.
    """
    s = 1e4 / np.asarray(w, dtype=float)
    n = 1 + 0.05792105 / (238.0185 - s * s) + 0.00167917 / (57.362 - s * s)
    return np.asarray(w, dtype=float) * n


def hydrogen_lines(wmin, wmax, series=(2, 3), nmax=40):
    """Hydrogen line centres in [wmin, wmax], from the Rydberg formula.

    series=2 is Balmer, 3 is Paschen. Analytic, so this cannot drift with the
    model grid.
    """
    out = []
    for m in series:
        for n in range(m + 1, nmax):
            lam = RYDBERG_A / (1.0 / m ** 2 - 1.0 / n ** 2)
            if wmin <= lam <= wmax:
                out.append(lam)
    return np.array(out)


def balmer_member(n):
    """Vacuum wavelength of the Balmer line n -> 2, from the Rydberg formula.

    BALMER below names only the four members with a usable local continuum.
    The high orders have no name worth writing down and their positions are
    exactly what the formula gives, so they are DERIVED rather than tabulated:
    against XSL the core minimum of H10 lands within 0.4 A of this value in all
    twelve sample stars, which is a tenth of an XSL resolution element.
    """
    return RYDBERG_A / (0.25 - 1.0 / float(n) ** 2)


def balmer_neighbourhood(n):
    """(lo, hi) for Balmer member n: the midpoints of the gaps to n+1 and n-1.

    From H-epsilon up the lines blend into one another -- the wing of H8 does
    not return to within 2% of the continuum until 122 A from centre -- so
    there is no line-free continuum to run a window out to, and the +/-50 A
    used for H-alpha through H-delta would swallow several neighbours whole.
    The midpoint of each gap is the one boundary the line positions themselves
    define, and it lands on the interline maximum: measured on XSL, the
    3785.4-3817.8 window around H10 peaks at both of its own edges.

    The window is ASYMMETRIC (-13.6 / +18.8 A at H10) because the series
    crowds blueward. Forcing it symmetric would either cross into the bluer
    neighbour's core or throw away half of the usable red wing.
    """
    lam = balmer_member(n)
    return 0.5 * (balmer_member(n + 1) + lam), 0.5 * (lam + balmer_member(n - 1))


# Named Balmer members. Only these four are used as XSL fit windows: H-epsilon
# and higher orders blend into each other, so a local continuum is not defined
# for them, and they sit inside the held-out break window besides. Measured at
# XSL resolution, the wing of H8 does not return to within 2% of the continuum
# until 122 A from centre, against 37-41 A for these four.
BALMER = {'Halpha': 6564.6, 'Hbeta': 4862.7, 'Hgamma': 4341.7, 'Hdelta': 4102.9}

# Interstellar features. These are absorption along the sightline, not in the
# star, so they must never enter a stellar fit -- and Ca II K is the single most
# [M/H]-sensitive feature in the optical at these temperatures, which makes it
# the most dangerous rather than the most useful: an ISM line read as stellar
# metallicity would bias [M/H] in the same direction as the reddening, and the
# result would look self-consistent.
ISM_LINES = {'Ca II K': 3934.8, 'Ca II H': 3969.6,
             'Na I D2': 5891.6, 'Na I D1': 5897.6}


def line_windows(centres, half_width, core_mask=0.0):
    """-> [(lo, hi)] windows around each centre, optionally with the core cut.

    Returns up to two intervals per line when `core_mask` is non-zero: the blue
    wing and the red wing, with the core between them dropped. The cores of the
    Balmer lines carry a flux excess of ~10% of the line equivalent width
    relative to these LTE models -- almost certainly NLTE in hydrogen, which the
    code does not treat for H -- so fitting them would drag Teff and log g to
    absorb physics the models are missing.
    """
    out = []
    for lam in np.atleast_1d(centres):
        if core_mask > 0:
            out.append((lam - half_width, lam - core_mask))
            out.append((lam + core_mask, lam + half_width))
        else:
            out.append((lam - half_width, lam + half_width))
    return out
