"""Identify which species dominates a spectral feature, from the Kurucz list.

Needed because line identifications cannot be assigned from memory. An earlier
version of docs/PLAN.md carried labels like "Fe II / Ti II blend" that were
written from expectation rather than derived, which is not a result.

Ranking by log gf alone is badly wrong at these temperatures. The most numerous
species in the Kurucz list over 4100-4900 A are Co I, V I and Nb I -- tens of
thousands of predicted transitions each -- and at 10,000 K essentially all of
that material is ionised, so those lines are absent from a real A-star spectrum.
The ranking therefore needs the full Saha-Boltzmann weight:

    log S = log(gf) + log(eps_element) + log(f_ionstage) - theta * E_low(eV)

with theta = 5040/T, evaluated at the line-forming depth of the actual model
atmosphere (T and N_e read from the .atm file at tau_5000 ~ 2/3).

Partition functions are the crude part: constant values per species, good to a
factor of ~2. That is tolerable here because the Saha exponential spans many
orders of magnitude between stages, so the stage assignment is not in doubt;
the ordering WITHIN a stage depends only on gf and E_low, which are exact.
"""
import numpy as np

LINELIST = '/Users/bjohnson/Projects/atlas12/data/gfallvac08oct17.dat'

EV_PER_CM1 = 1.239841984e-4
K_EV = 8.617333262e-5

ELEMENT = {1: 'H', 2: 'He', 6: 'C', 7: 'N', 8: 'O', 11: 'Na', 12: 'Mg',
           13: 'Al', 14: 'Si', 16: 'S', 20: 'Ca', 21: 'Sc', 22: 'Ti',
           23: 'V', 24: 'Cr', 25: 'Mn', 26: 'Fe', 27: 'Co', 28: 'Ni',
           29: 'Cu', 30: 'Zn', 38: 'Sr', 39: 'Y', 40: 'Zr', 41: 'Nb',
           56: 'Ba', 57: 'La', 58: 'Ce', 60: 'Nd', 62: 'Sm', 63: 'Eu'}

ROMAN = ['I', 'II', 'III', 'IV']

# First and second ionization potentials, eV (NIST).
IP = {1: (13.598, np.inf), 2: (24.587, 54.418), 6: (11.260, 24.383),
      7: (14.534, 29.601), 8: (13.618, 35.121), 11: (5.139, 47.286),
      12: (7.646, 15.035), 13: (5.986, 18.829), 14: (8.152, 16.346),
      16: (10.360, 23.338), 20: (6.113, 11.872), 21: (6.561, 12.800),
      22: (6.828, 13.576), 23: (6.746, 14.618), 24: (6.767, 16.486),
      25: (7.434, 15.640), 26: (7.902, 16.199), 27: (7.881, 17.084),
      28: (7.640, 18.169), 29: (7.726, 20.292), 30: (9.394, 17.964),
      38: (5.695, 11.030), 39: (6.217, 12.224), 40: (6.634, 13.130),
      41: (6.759, 14.320), 56: (5.212, 10.004), 57: (5.577, 11.060),
      58: (5.539, 10.850), 60: (5.525, 10.730), 62: (5.644, 11.070),
      63: (5.670, 11.240)}

# Partition functions near 10,000 K, order of magnitude only (see module note).
U = {(26, 0): 30., (26, 1): 47., (22, 0): 30., (22, 1): 55.,
     (24, 0): 12., (24, 1): 8., (25, 0): 8., (25, 1): 10.,
     (28, 0): 32., (28, 1): 12., (23, 0): 40., (23, 1): 35.,
     (27, 0): 30., (27, 1): 30., (12, 0): 1., (12, 1): 2.,
     (14, 0): 9., (14, 1): 6., (20, 0): 2., (20, 1): 2.,
     (21, 0): 12., (21, 1): 18., (13, 0): 6., (13, 1): 1.,
     (11, 0): 2., (11, 1): 1., (56, 0): 2., (56, 1): 2.,
     (38, 0): 2., (38, 1): 2., (39, 0): 10., (39, 1): 12.,
     (40, 0): 25., (40, 1): 30., (41, 0): 30., (41, 1): 40.}
U_DEFAULT = 20.0

SAHA_CONST = 2.414e15        # (2 pi m_e k / h^2)^{3/2} in cgs, per cm^3 K^-3/2


def atmosphere_point(atm_path, tau=2.0 / 3.0):
    """-> (T, N_e) at a given tau_5000 in an ATLAS12 .atm file."""
    rows = []
    for line in open(atm_path):
        p = line.split()
        if len(p) >= 11:
            try:
                rows.append([float(x) for x in p[:11]])
            except ValueError:
                pass
    a = np.array(rows)
    i = int(np.argmin(np.abs(a[:, 10] - tau)))
    return float(a[i, 1]), float(a[i, 3])


def abundances(atm_path):
    """-> {Z: log eps relative to H} from the .atm ABUNDANCE TABLE."""
    import re
    out, reading = {}, False
    for line in open(atm_path):
        if 'ABUNDANCE TABLE' in line:
            reading = True
            continue
        if reading:
            if not re.match(r'\s*\d', line):
                break
            for m in re.finditer(r'(\d+)\s*([A-Za-z]{1,2})\s*(-?\d+\.\d+)', line):
                z, val = int(m.group(1)), float(m.group(3))
                out[z] = val if val < 0 else np.log10(max(val, 1e-30))
    return out


def ion_fraction(z, stage, T, ne):
    """Fraction of element z in `stage` (0 = neutral), from two-stage Saha."""
    if z not in IP:
        return 1.0 if stage == 0 else 0.0
    n = np.array([1.0, 0.0, 0.0])
    for s in (0, 1):
        chi = IP[z][s]
        if not np.isfinite(chi):
            break
        u_lo = U.get((z, s), U_DEFAULT)
        u_hi = U.get((z, s + 1), U_DEFAULT)
        ratio = (2.0 * SAHA_CONST * T ** 1.5 * (u_hi / u_lo)
                 * np.exp(-chi / (K_EV * T)) / max(ne, 1.0))
        n[s + 1] = n[s] * ratio
    n /= n.sum()
    return float(n[stage]) if stage < len(n) else 0.0


def read_lines(lo_A, hi_A, path=LINELIST):
    """Kurucz lines in [lo, hi] Angstroms -> list of (lam_A, loggf, z, stage, elow_eV).

    The file is wavelength-sorted in NANOMETRES, so the scan stops early.
    """
    lo_nm, hi_nm = lo_A / 10.0, hi_A / 10.0
    out = []
    with open(path, 'r', errors='replace') as fh:
        for line in fh:
            try:
                wl = float(line[0:11])
            except ValueError:
                continue
            if wl < lo_nm:
                continue
            if wl > hi_nm:
                break
            try:
                gf = float(line[11:18])
                code = float(line[18:24])
                elow = float(line[24:36])
            except ValueError:
                continue
            z = int(code)
            stage = int(round((code - z) * 100))
            out.append((wl * 10.0, gf, z, stage, abs(elow) * EV_PER_CM1))
    return out


# Two species closer than this in log strength are reported as a blend: the
# constant partition functions are good to about a factor of 2, so a smaller
# difference is not resolvable by this method and claiming one dominates would
# be over-reading it.
BLEND_DEX = 0.3


def species_label(lo_A, hi_A, T, ne, eps, path=LINELIST):
    """-> a label like 'Fe II' or 'Cr II + O I' for a wavelength window."""
    sp = dominant_species(lo_A, hi_A, T, ne, eps, path=path, top=2)
    if not sp:
        return '?'
    if len(sp) > 1 and abs(sp[1][1]) <= BLEND_DEX:
        return f'{sp[0][0]} + {sp[1][0]}'
    return sp[0][0]


def dominant_species(lo_A, hi_A, T, ne, eps, path=LINELIST, top=2):
    """-> [(label, relative_strength)] for the strongest species in a window."""
    lines = read_lines(lo_A, hi_A, path)
    theta = 5040.0 / T
    best = {}
    for lam, gf, z, stage, elow in lines:
        if z not in ELEMENT or stage > 2:
            continue
        f_ion = ion_fraction(z, stage, T, ne)
        if f_ion <= 0:
            continue
        s = gf + eps.get(z, -12.0) + np.log10(f_ion) - theta * elow
        key = (z, stage)
        if key not in best or s > best[key][0]:
            best[key] = (s, lam)
    if not best:
        return []
    ranked = sorted(best.items(), key=lambda kv: -kv[1][0])
    smax = ranked[0][1][0]
    return [(f'{ELEMENT[z]} {ROMAN[stage]}', float(s - smax), float(lam))
            for (z, stage), (s, lam) in ranked[:top]]


def strong_lines(lo_A, hi_A, T, ne, eps, n=4, within=2.0, path=LINELIST):
    """-> [(lam, 'Fe II', rel_strength)] for the n strongest lines in a window.

    For annotating a panel: which lines are actually there, so the figure names
    them instead of leaving the reader to guess from a window centre.
    `within` keeps only lines this many dex of the strongest.
    """
    theta = 5040.0 / T
    rows = []
    for lam, gf, z, stage, elow in read_lines(lo_A, hi_A, path):
        if z not in ELEMENT or stage > 2:
            continue
        f_ion = ion_fraction(z, stage, T, ne)
        if f_ion <= 0:
            continue
        rows.append((gf + eps.get(z, -12.0) + np.log10(f_ion) - theta * elow,
                     lam, f'{ELEMENT[z]} {ROMAN[stage]}'))
    if not rows:
        return []
    rows.sort(reverse=True)
    smax = rows[0][0]
    out = []
    for s, lam, sp in rows:
        if s < smax - within:
            break
        if any(abs(lam - l) < 0.25 for l, _, _ in out):
            continue          # same feature, already labelled
        out.append((float(lam), sp, float(s - smax)))
        if len(out) >= n:
            break
    return out
