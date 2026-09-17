"""Select the fitted sample: stars observed by BOTH NGSL and XSL near 10,000 K.

The two libraries together are what make the Balmer-break test possible, and
each supplies something the other cannot:

  NGSL   space-based spectrophotometry. The continuum slope blueward of the
         Balmer limit (3200-3500 A) is nearly line free in an A star and is
         measured from space, so it constrains dust without the atmospheric
         transmission problem that makes the near-UV hard from the ground.
  XSL    ~16x the resolving power. Locally normalised line profiles CANNOT be
         altered by a smooth reddening law -- a degree-4 polynomial absorbs
         CCM89 over a 1100 A window to 3e-5 -- so XSL carries Teff, log g,
         [M/H] and v sin i information that is immune to the dust degeneracy.

So the sample is the INTERSECTION, not the union: a star in only one library
loses one of the two levers. That is why HD040573, a mainstay of the earlier
NGSL-only comparison, is not here -- XSL never observed it.

Selection, in order:

  Teff window     9000-11000 K, satisfied by EITHER catalog. The two disagree
                  by up to 2900 K for the same star (HD164257: XSL 10885 vs
                  NGSL 7977), which is the entire reason parameters are fitted
                  here rather than adopted. Requiring both to agree would let
                  one catalog's systematic define the sample.
  peculiarity     Ap/Bp/Am and magnetic stars are REJECTED. A scaled-solar
                  model atmosphere does not describe a star with abundance
                  patches, so it cannot test one. Screened on SIMBAD spectral
                  type and object type, plus the abundance/rotation signature
                  (very low v sin i with strongly super-solar [Fe/H]).
  binarity        FLAGGED, NOT REJECTED -- a secondary sample. A composite
                  spectrum has a wavelength-dependent flux ratio, which
                  distorts continuum shape specifically, so binaries test the
                  method differently rather than uselessly.

Reddening from every available source is recorded but NONE is trusted to the
precision this analysis needs (0.0134 mag of predicted D_Balmer per 0.01 mag of
E(B-V)). SFD98/SF11 are total Galactic columns and so upper bounds only -- they
overshoot by up to 8x for stars this close. They are here to bound a prior, not
to centre one.

Writes data/sample.csv -- rejected rows are KEPT, with the reason, so the
selection history survives.
"""
import csv
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common.sample import as_float
from common.dust import (irsa_ebv, simbad_info, ebv_photometric, is_binary,
                         is_peculiar, is_peculiar_otype, is_horizontal_branch)

ROOT = Path(__file__).resolve().parent.parent
TEFF_WINDOW = (9000.0, 11000.0)
MATCH_ARCSEC = 5.0

# The model grid (grid/build_grid.py). A star outside it cannot be fitted by a
# node scan -- the scan piles up on the boundary instead of failing, so the
# coverage is recorded per star and checked rather than discovered later.
GRID = dict(teff=(8500.0, 11500.0), logg=(3.0, 5.0), mh=(-0.5, 0.3))

# Am/Ap signature that a spectral type may miss: sharp lines plus strong metals.
PEC_VSINI_MAX, PEC_FEH_MIN = 30.0, 0.40


def galactic_b(ra, dec):
    """Galactic latitude in degrees (J2000 pole), for judging the map column.

    An SFD/SF11 value is only a useful upper bound if it is small. In the plane
    it is the column through the whole Galaxy and can reach several magnitudes,
    which bounds nothing -- HD174240 sits at b = +0.8 deg and returns SF11 =
    5.00, verified by re-query rather than assumed to be a parse error.
    """
    ra_ngp, dec_ngp = np.radians(192.85948), np.radians(27.12825)
    a, d = np.radians(ra), np.radians(dec)
    return np.degrees(np.arcsin(np.sin(dec_ngp) * np.sin(d)
                                + np.cos(dec_ngp) * np.cos(d) * np.cos(a - ra_ngp)))


# Above this the map column bounds nothing useful and is reported as such.
EBV_MAP_USELESS = 0.5


def sep_arcsec(ra1, dec1, ra2, dec2):
    """Great-circle separation, degrees in -> arcsec out."""
    r1, d1, r2, d2 = map(np.radians, (ra1, dec1, ra2, dec2))
    c = np.sin(d1) * np.sin(d2) + np.cos(d1) * np.cos(d2) * np.cos(r1 - r2)
    return np.degrees(np.arccos(np.clip(c, -1.0, 1.0))) * 3600.0


def crossmatch():
    """-> [(xsl_row, ngsl_row, sep_arcsec)] one entry per XSL star."""
    ngsl = list(csv.DictReader(open(ROOT / 'data' / 'ngsl_catalog.csv')))
    nra = np.array([as_float(r['ra']) for r in ngsl])
    ndec = np.array([as_float(r['dec']) for r in ngsl])

    first = {}
    for r in csv.DictReader(open(ROOT / 'data' / 'xsl_all.csv')):
        first.setdefault(r['star'], r)          # XSL repeats stars per epoch

    out = []
    for x in first.values():
        ra, dec = as_float(x['ra']), as_float(x['dec'])
        if not np.isfinite(ra):
            continue
        s = sep_arcsec(ra, dec, nra, ndec)
        j = int(np.argmin(s))
        if s[j] < MATCH_ARCSEC:
            out.append((x, ngsl[j], float(s[j])))
    return out, len(ngsl), len(first)


# Gaia RUWE above this is treated as a binarity detection. RUWE is a statement
# about the astrometric fit of THIS star, so unlike a SIMBAD object type it
# cannot simply be missing from a catalog -- and a missing object type is
# exactly how a spectroscopic binary got into this project before (HD162630
# looked normal on abundance and rotation and is an SB).
RUWE_BINARY = 1.4


def gaia_flags():
    """-> {star: dict(ruwe, non_single_star, binary)} from data/gaia_sample.csv.

    Optional: the file is produced by explore/fetch_gaia.py, which reads THIS
    script's output for its positions. So the pipeline is two-pass -- run
    build_sample, fetch_gaia, then build_sample again to fold the astrometry in.
    Whether it was available is printed, so a one-pass run cannot be mistaken
    for a vetted one.
    """
    p = ROOT / 'data' / 'gaia_sample.csv'
    if not p.exists():
        return {}
    out = {}
    for r in csv.DictReader(open(p)):
        ruwe = as_float(r.get('ruwe', ''))
        nss = str(r.get('non_single_star', '0')).strip()
        out[r['star']] = dict(
            ruwe=ruwe, non_single_star=nss,
            binary=(np.isfinite(ruwe) and ruwe > RUWE_BINARY)
                   or nss not in ('0', '', '--', 'None'))
    return out


def miles_ebv():
    """Fitted E(B-V) from MILES, keyed by NGSL target name, where it exists."""
    p = ROOT / 'data' / 'ngsl_crossmatch.csv'
    if not p.exists():
        return {}
    out = {}
    for r in csv.DictReader(open(p)):
        key = r.get('target') or r.get('star') or ''
        v = as_float(r.get('miles_ebv', ''))
        if key and np.isfinite(v):
            out[key] = v
    return out


def main():
    matches, n_ngsl, n_xsl = crossmatch()
    print(f'NGSL {n_ngsl} spectra, XSL {n_xsl} unique stars')
    print(f'  positional matches within {MATCH_ARCSEC:.0f}": {len(matches)}')

    cand = [(x, n, s) for x, n, s in matches
            if (TEFF_WINDOW[0] <= as_float(x['teff']) <= TEFF_WINDOW[1])
            or (TEFF_WINDOW[0] <= as_float(n['teff']) <= TEFF_WINDOW[1])]
    print(f'  either Teff in {TEFF_WINDOW[0]:.0f}-{TEFF_WINDOW[1]:.0f} K: '
          f'{len(cand)} stars\n')

    mil = miles_ebv()
    gaia = gaia_flags()
    print('  Gaia astrometry: ' + (f'{len(gaia)} stars from data/gaia_sample.csv'
                                   if gaia else
                                   'NOT AVAILABLE -- binarity is from SIMBAD '
                                   'only; run explore/fetch_gaia.py then rerun'))
    print()
    rows = []
    for x, n, sep in sorted(cand, key=lambda t: -as_float(t[0]['teff'])):
        name = n['target']
        info = simbad_info(x['simbad'] or x['star'])
        sptype = info.get('sp_type', '') or n.get('sptype', '')
        otype = info.get('otype', '')

        try:
            sfd, sf11 = irsa_ebv(as_float(x['ra']), as_float(x['dec']))
        except Exception as exc:
            print(f'  ! IRSA failed for {name}: {type(exc).__name__}')
            sfd = sf11 = np.nan

        bv = as_float(n['bmag']) - as_float(n['vmag'])
        ebv_phot = ebv_photometric(sptype, round(bv, 3) if np.isfinite(bv) else '')

        t_xsl, t_ngsl = as_float(x['teff']), as_float(n['teff'])
        g_xsl, g_ngsl = as_float(x['logg']), as_float(n['logg'])
        z_xsl, z_ngsl = as_float(x['feh']), as_float(n['logz'])

        gf = gaia.get(name, {})
        binary_simbad = is_binary(otype, sptype)
        binary_gaia = bool(gf.get('binary'))
        binary = binary_simbad or binary_gaia
        hb = is_horizontal_branch(otype)
        pec_type = is_peculiar(sptype) or is_peculiar_otype(otype)
        # Abundance signature, which a spectral type can miss entirely. CAVEATS
        # states it as low v sin i AND strongly super-solar [Fe/H]; no v sin i
        # is available for the XSL stars, so on abundance alone this is the
        # weaker half of that test -- and a COMPOSITE spectrum fakes high
        # metallicity too. So it rejects only when nothing says binary; with a
        # binary flag the abundance is more likely dilution than peculiarity,
        # and the star goes to the secondary sample where that can be examined.
        pec_abund = np.isfinite(z_xsl) and z_xsl >= PEC_FEH_MIN
        pec = pec_type or (pec_abund and not binary)

        def inside(v, lo_hi):
            lo, hi = lo_hi
            return '' if not np.isfinite(v) else ('yes' if lo <= v <= hi else 'no')

        # grid coverage judged on BOTH catalogs -- either being outside is worth
        # knowing, since the scan has no way to extrapolate
        out_of_grid = [nm for nm, vals, key in
                       (('teff', (t_xsl, t_ngsl), 'teff'),
                        ('logg', (g_xsl, g_ngsl), 'logg'),
                        ('mh', (z_xsl, z_ngsl), 'mh'))
                       if any(np.isfinite(v) and not (GRID[key][0] <= v <= GRID[key][1])
                              for v in vals)]

        reject = ''
        if pec:
            why = 'sptype/otype' if pec_type else f'[Fe/H]_xsl={z_xsl:+.2f}'
            reject = f'chemically peculiar ({why})'
        tier = 'primary'
        if reject:
            tier = 'rejected'
        elif binary:
            tier = 'secondary'          # kept, flagged: composite spectrum

        rows.append(dict(
            star=name, xsl_name=x['star'], xslid=x['xslid'],
            ngsl_file=n['file'], xsl_file=x['filename'],
            ra=round(as_float(x['ra']), 5), dec=round(as_float(x['dec']), 5),
            match_arcsec=round(sep, 2),
            teff_xsl=('' if not np.isfinite(t_xsl) else round(t_xsl)),
            teff_ngsl=('' if not np.isfinite(t_ngsl) else round(t_ngsl)),
            logg_xsl=('' if not np.isfinite(g_xsl) else round(g_xsl, 2)),
            logg_ngsl=('' if not np.isfinite(g_ngsl) else round(g_ngsl, 2)),
            mh_xsl=('' if not np.isfinite(z_xsl) else round(z_xsl, 2)),
            mh_ngsl=('' if not np.isfinite(z_ngsl) else round(z_ngsl, 2)),
            sptype=sptype, otype=otype, dist_pc=info.get('dist_pc', ''),
            vmag=n['vmag'], bmag=n['bmag'],
            offset_px=n['offset_px'], dataqual=n['dataqual'],
            ebv_sfd98=('' if not np.isfinite(sfd) else round(sfd, 4)),
            ebv_sf11=('' if not np.isfinite(sf11) else round(sf11, 4)),
            ebv_phot=ebv_phot,
            ebv_miles=mil.get(name, ''),
            gal_b=round(galactic_b(as_float(x['ra']), as_float(x['dec'])), 2),
            ebv_map_useful=('no' if (np.isfinite(sf11) and sf11 > EBV_MAP_USELESS)
                            else 'yes'),
            binary='yes' if binary else 'no',
            binary_evidence=','.join(
                [s for s, c in (('simbad', binary_simbad),
                                ('gaia', binary_gaia)) if c]),
            ruwe=(round(gf['ruwe'], 2) if np.isfinite(gf.get('ruwe', np.nan))
                  else ''),
            peculiar='yes' if pec else 'no',
            horizontal_branch='yes' if hb else 'no',
            high_feh='yes' if pec_abund else 'no',
            out_of_grid=','.join(out_of_grid),
            tier=tier, selected='no' if reject else 'yes', reject_reason=reject,
            xsl_comment=x['comment']))

    with open(ROOT / 'data' / 'sample.csv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    print(f'{"star":<11}{"T_xsl":>7}{"T_ngsl":>7}{"g_xsl":>6}{"z_xsl":>6}  '
          f'{"sptype":<12}{"otyp":<5}{"b":>6}{"SF11":>7}{"Eph":>6}  '
          f'{"tier":<10}{"flags / out-of-grid"}')
    for r in rows:
        flags = [k for k, c in (('HB', 'horizontal_branch'),
                                ('hiZ', 'high_feh')) if r[c] == 'yes']
        if r['binary'] == 'yes':
            # quote RUWE only when RUWE is the evidence -- printing a clean
            # RUWE next to BIN(simbad) reads as if 0.92 were the reason
            ev = r['binary_evidence']
            if 'gaia' in ev and r['ruwe']:
                ev = ev.replace('gaia', f'RUWE={r["ruwe"]}')
            flags.insert(0, f'BIN({ev})')
        if r['ebv_map_useful'] == 'no':
            flags.append('map-useless')
        if r['out_of_grid']:
            flags.append('grid:' + r['out_of_grid'])
        print(f'{r["star"]:<11}{str(r["teff_xsl"]):>7}{str(r["teff_ngsl"]):>7}'
              f'{str(r["logg_xsl"]):>6}{str(r["mh_xsl"]):>6}  '
              f'{r["sptype"][:12]:<12}{r["otype"][:5]:<5}'
              f'{r["gal_b"]:>6.1f}{str(r["ebv_sf11"]):>7}{str(r["ebv_phot"]):>6}  '
              f'{r["tier"]:<10}{" ".join(flags)}'
              + (f'   <- {r["reject_reason"]}' if r['reject_reason'] else ''))

    tally = {}
    for r in rows:
        tally[r['tier']] = tally.get(r['tier'], 0) + 1
    print(f'\n{len(rows)} stars -> data/sample.csv')
    for k in ('primary', 'secondary', 'rejected'):
        print(f'  {k:<10}{tally.get(k, 0)}')
    ng = [r['star'] for r in rows if r['out_of_grid'] and r['selected'] == 'yes']
    if ng:
        print(f'  outside the model grid (scan will hit a boundary): {", ".join(ng)}')


if __name__ == '__main__':
    main()
