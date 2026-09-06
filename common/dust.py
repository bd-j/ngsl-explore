"""Reddening estimates and distances, shared by selection and fitting.

Three sources, in decreasing order of trustworthiness for a given star:

  fitted      a value measured along the actual sightline (UVES-POP publishes
              one per star; MILES too). Use it.
  tomographic a 3D map evaluated at the star's distance. Not yet wired up --
              Bayestar19 is queryable through the Argonaut API without a map
              download but is PS1-based and covers Dec > -30 only; Lallement /
              Vergely and Edenhofer cover the south more coarsely.
  map column  SFD98 / SF11 from IRSA. This is the TOTAL column through the
              whole Galactic dust layer and therefore an UPPER BOUND for a star
              inside it -- for our stars at 137-317 pc it overshoots by up to
              8x. Never use it as a central value; use it to bound a prior.
"""
import re
import urllib.request

import numpy as np

IRSA = ('https://irsa.ipac.caltech.edu/cgi-bin/DUST/nph-dust'
        '?locstr={ra}+{dec}+equ+j2000')

# Intrinsic (B-V)_0 by spectral type, Pecaut & Mamajek (2013) / Fitzgerald (1970)
BV0 = {'B7': -0.13, 'B8': -0.11, 'B9': -0.07, 'B9.5': -0.05, 'A0': 0.00,
       'A1': 0.03, 'A2': 0.06, 'A3': 0.09, 'A5': 0.15, 'A7': 0.20}


def irsa_ebv(ra, dec, timeout=90):
    """-> (SFD98, SF11) total-column E(B-V) at the reference pixel."""
    req = urllib.request.Request(IRSA.format(ra=ra, dec=dec),
                                 headers={'User-Agent': 'Mozilla/5.0'})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        x = r.read().decode('utf-8', 'replace')

    def grab(tag):
        m = re.search(r'<' + tag + r'>\s*([\d.]+)\s*\(mag\)\s*</' + tag + r'>', x)
        return float(m.group(1)) if m else np.nan
    return grab('refPixelValueSFD'), grab('refPixelValueSandF')


def ebv_photometric(sptype, b_minus_v):
    """E(B-V) = (B-V)_obs - (B-V)_0 from the spectral type, or ''."""
    if b_minus_v in ('', None) or not sptype:
        return ''
    m = re.match(r'([OBAFGKM]\d?(?:\.\d)?)', str(sptype).strip())
    if not m or m.group(1) not in BV0:
        return ''
    return round(float(b_minus_v) - BV0[m.group(1)], 3)


def simbad_info(name):
    """-> dict(sp_type, otype, plx, dist_pc, binary, peculiar) or empty."""
    from astroquery.simbad import Simbad
    sb = Simbad()
    for f in ('sp_type', 'otype', 'plx_value'):
        try:
            sb.add_votable_fields(f)
        except Exception:
            pass
    try:
        q = sb.query_object(name)
        if q is None:
            return {}
        cn = {c.lower(): c for c in q.colnames}
        g = lambda k: q[cn[k]][0] if k in cn else None
        sp = str(g('sp_type') or '').strip()
        ot = str(g('otype') or '').strip()
        plx = g('plx_value')
        try:
            plx = float(plx)
        except (TypeError, ValueError):
            plx = np.nan
        return dict(
            sp_type=sp, otype=ot,
            dist_pc=(round(1000.0 / plx, 1) if np.isfinite(plx) and plx > 0 else ''),
            binary=is_binary(ot, sp),
            peculiar=is_peculiar(sp) or is_peculiar_otype(ot))
    except Exception:
        return {}


# SIMBAD object types implying multiplicity. Parameter cuts do NOT catch these:
# HD162630 looked entirely normal on abundance and rotation and is an SB.
BINARY_OTYPES = {'SB*', '**', 'EB*', 'SB', 'Al*', 'El*', 'RS*', 'bC*', 'Sy*',
                 'CV*', 'XB*', 'LXB', 'HXB'}


def is_binary(otype, sptype=''):
    if str(otype).strip() in BINARY_OTYPES:
        return True
    return '+' in str(sptype)          # e.g. 'A1V+A2Vm'


# SIMBAD object types implying chemical peculiarity or a magnetic field.
# a2* is the Alpha2 CVn class: magnetic Ap stars with abundance patches.
PECULIAR_OTYPES = {'a2*', 'Ap*', 'rC*', 'HB*'}


def is_peculiar_otype(otype):
    return str(otype).strip() in PECULIAR_OTYPES


def is_peculiar(sptype):
    """Ap/Bp/Am and kin: abundance patches or anomalous metals, so a
    scaled-solar model atmosphere does not describe them."""
    s = str(sptype)
    if any(k in s for k in ('Ap', 'Bp', 'pec', 'Cr', 'Eu', 'Sr', 'Si', 'Hg', 'Mn')):
        return True
    # 'p' after a subtype or luminosity class: A0Vp, B9IIIp, A2p, A1VpSrCr.
    # No trailing \b -- the peculiarity suffix is often followed by the
    # element list, as in A1VpSrCr, which slipped through with one.
    if re.search(r'(?:[IV]|\d)p', s):
        return True
    # metallic-lined: A1V+A2Vm, A0mA1Va, kA2hA5mA7
    return bool(re.search(r'[IV]m\b|m[AF]\d|^[AF]\d?m', s))
