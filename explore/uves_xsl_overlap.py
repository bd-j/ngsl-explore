"""Stars observed by both UVES-POP and XSL DR3, with the parameters from each.

Both are VLT libraries, so unlike the NGSL x UVES-POP overlap (13 stars, none
in the Balmer window -- explore/uves_ngsl_overlap.py) there is no hemisphere
penalty here. The overlap is nonetheless SMALLER: 9 stars. The two programs
chose nearly disjoint targets -- UVES-POP took bright nearby stars for a
high-resolution atlas, XSL took a stellar-population grid weighted to cool
giants, the bulge and the Magellanic Clouds.

Neither catalog carries the other's identifiers and UVES-POP names its
brightest targets by proper name (Achernar, Betelgeuse, ...), so the match is
positional: SIMBAD resolves the UVES-POP names, and those coordinates are
matched against the XSL DR3 positions in data/xsl_all.csv.

Parameter provenance -- the two sets are independent, which is the point of
putting them side by side:

  XSL       Arentsen et al. (2019), who derived Teff, log g and [Fe/H] for
            the XSL stars -- the DR3 table itself carries only names and
            filenames (see explore/xsl_astars.py, which fetches them). Errors
            are that catalog's, carried through unchanged.
  UVES-POP  VOXAstro-SL re-reduction (2023ApJS..266...11B), full-spectrum fit
            against PHOENIX, with formal errors and v sin i. The quoted errors
            are formal fit errors only -- e_teff runs to ~1 K -- and are not an
            uncertainty budget. Do not read them as accuracy.

Sign convention: every d_* column is UVES-POP minus XSL, so a positive d_teff
means UVES-POP runs hotter. This is asserted against an injected pair below and
the check runs every time.

The UVES-POP `v` column is the radial velocity in km/s, NOT the V magnitude;
that was established in explore/uves_ngsl_overlap.py against SIMBAD (tracks RV
to 3.7 km/s, misses V by up to 114 mag) and is carried through as rv_uves_kms.

Writes data/uves_xsl_overlap.csv

Caveat: the hottest star in the overlap is 7595 K, so like the NGSL overlap
this one contains no A star in the Balmer window and does nothing for the
break. Its use is as an independent parameter-scale check.
"""
import csv
import sys

import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord
from astroquery.simbad import Simbad

# The 9 matches all fall below 0.1"; the next-nearest pair is at 430" (beta Crv
# vs HD 109443, two different stars). Any cut in that four-decade gap gives the
# same 9, so the answer does not depend on where the line goes. 5" is the
# radius build_sample.py uses for XSL and is kept here.
#
# The separations are two orders of magnitude tighter than the NGSL match,
# where four stars landed between 1.7" and 7.3". That is expected and is a
# check in itself: NGSL positions come from the FITS headers at the epoch of
# observation, so proper motion shows up, whereas both catalogs here carry
# J2000 catalog positions.
RADIUS = 5 * u.arcsec
# Near-misses out to here are printed so the rejected ones stay visible.
NEAR = 600 * u.arcsec

# UVES-POP's open-cluster targets carry the library's own running numbers.
# SIMBAD knows them under the Cl* designations below.
CLUSTERS = {
    'IC2391': ('IC 2391', (130.05, -53.03), ['SHJM', 'PP', 'VXR']),
    'NGC6475': ('NGC 6475', (268.4471, -34.8411), ['JJ']),
}
# An unresolved cluster entry is only safe to drop if XSL has nothing in that
# field. Enforced below, not assumed: this margin is the distance from the
# field centre inside which an XSL star would force the entry to be resolved.
CLUSTER_CLEAR_DEG = 3.0


def fnum(x, nd=None):
    """Catalog string -> rounded float, or '' when absent."""
    try:
        v = float(x)
    except (TypeError, ValueError):
        return ''
    if not np.isfinite(v):
        return ''
    return round(v, nd) if nd is not None else v


def delta(a_uves, b_xsl, nd):
    """UVES-POP minus XSL, or '' when either side is missing."""
    if '' in (a_uves, b_xsl):
        return ''
    return round(a_uves - b_xsl, nd)


# Sign check with an injected pair: XSL 6000 K, UVES-POP 6100 K must give +100.
assert delta(6100.0, 6000.0, 1) == 100.0, 'd_* must be UVES-POP minus XSL'
assert delta(3.0, 4.0, 2) == -1.0, 'd_* must be UVES-POP minus XSL'


def fmt(v, nd):
    """Fixed-decimal for the printed table, blank when the value is absent."""
    return '' if v == '' else f'{v:.{nd}f}'


def note(x, uv):
    """Why a row's parameters may not be comparable. Recorded, not eyeballed."""
    n = []
    if not fnum(x['teff']):
        n.append('xsl_unfitted')
    if not fnum(uv['teff']):
        n.append('uves_unfitted')
    g = fnum(uv['logg'])
    if g != '' and g < 0.5:
        n.append(f'uves_logg={g:.2f}_at_grid_edge')
    if x['comment'].strip():
        n.append(f'xsl_comment={x["comment"].strip()}')
    return n


def resolve_uves(names):
    """-> (ra, dec, main_id) arrays. Cluster running numbers are retried under
    their Cl* designations rather than being dropped."""
    t = Simbad().query_objects(names)
    ra = np.asarray(t['ra'], float)
    dec = np.asarray(t['dec'], float)
    mid = [str(v) for v in t['main_id']]
    print(f'  plain names: {np.isfinite(ra).sum()}/{len(names)}')

    todo = [i for i, r in enumerate(ra) if not np.isfinite(r)]
    for pre, (spaced, _, cats) in CLUSTERS.items():
        idx = [i for i in todo if names[i].startswith(pre)]
        for cat in cats:
            idx = [i for i in idx if not np.isfinite(ra[i])]
            if not idx:
                break
            # 'IC2391-0001' -> 'Cl* IC 2391 SHJM 1'; the NGC 6475 numbers
            # already carry their JJ prefix, so it is not doubled.
            num = [names[i].split('-', 1)[1] for i in idx]
            variants = [f'Cl* {spaced} {n if n.startswith(cat) else cat + " " + n.lstrip("0")}'
                        for n in num]
            tc = Simbad().query_objects(variants)
            # Unresolved variants come back masked; fill to nan up front rather
            # than letting float() do it one warning at a time.
            c_ra = np.ma.filled(np.ma.masked_invalid(tc['ra']), np.nan)
            c_dec = np.ma.filled(np.ma.masked_invalid(tc['dec']), np.nan)
            for k, i in enumerate(idx):
                if np.isfinite(c_ra[k]):
                    ra[i], dec[i] = float(c_ra[k]), float(c_dec[k])
                    mid[i] = str(tc['main_id'][k])
    return ra, dec, mid


def check_clusters(names, ra, xsl_coord):
    """Unresolved cluster entries are droppable only if XSL has nothing in the
    field. Asserted, so a future XSL release in one of these fields fails loudly
    instead of silently losing a match."""
    for pre, (spaced, (cra, cdec), _) in CLUSTERS.items():
        left = [n for n, r in zip(names, ra)
                if n.startswith(pre) and not np.isfinite(r)]
        d = SkyCoord(cra * u.deg, cdec * u.deg).separation(xsl_coord).deg
        print(f'  {spaced}: {len(left)} entries unresolved; nearest XSL star '
              f'{d.min():.2f} deg from the field centre, '
              f'{(d < CLUSTER_CLEAR_DEG).sum()} within {CLUSTER_CLEAR_DEG:.0f} deg')
        assert not left or d.min() > CLUSTER_CLEAR_DEG, (
            f'XSL has a star within {CLUSTER_CLEAR_DEG} deg of {spaced} while '
            f'{len(left)} UVES-POP entries there are unresolved -- resolve them '
            f'before trusting the match')


def main():
    uves = list(csv.DictReader(open('data/uves_pop_all.csv')))

    # XSL repeats a star once per epoch; collapse to unique stars but keep the
    # epoch count, since a repeat is a free repeatability check on the spectra.
    first, epochs = {}, {}
    for r in csv.DictReader(open('data/xsl_all.csv')):
        first.setdefault(r['star'], r)
        epochs.setdefault(r['star'], []).append(r['xslid'])
    xsl = list(first.values())
    n_spec = sum(len(v) for v in epochs.values())
    print(f'UVES-POP {len(uves)} stars, XSL {n_spec} spectra of {len(xsl)} stars')

    x_c = SkyCoord([float(r['ra']) for r in xsl] * u.deg,
                   [float(r['dec']) for r in xsl] * u.deg)

    print('Resolving UVES-POP names in SIMBAD...')
    names = [r['name'] for r in uves]
    u_ra, u_dec, u_id = resolve_uves(names)
    ok = np.isfinite(u_ra) & np.isfinite(u_dec)
    print(f'  resolved {ok.sum()}/{len(uves)}')
    check_clusters(names, u_ra, x_c)

    u_c = SkyCoord(u_ra[ok] * u.deg, u_dec[ok] * u.deg)
    u_rows = [r for r, k in zip(uves, ok) if k]
    u_ids = [i for i, k in zip(u_id, ok) if k]

    idx, d2d, _ = u_c.match_to_catalog_sky(x_c)
    hit = d2d < RADIUS
    print(f'\nMatched within {RADIUS}: {hit.sum()} stars')

    rows = []
    for i in np.where(hit)[0]:
        uv, x = u_rows[i], xsl[idx[i]]
        t_x, t_u = fnum(x['teff'], 1), fnum(uv['teff'], 1)
        g_x, g_u = fnum(x['logg'], 3), fnum(uv['logg'], 3)
        z_x, z_u = fnum(x['feh'], 3), fnum(uv['fe_h'], 3)
        rows.append(dict(
            xsl_name=x['star'], uves_name=uv['name'], simbad_id=u_ids[i],
            match_arcsec=round(float(d2d[i].arcsec), 3),
            ra=round(float(x['ra']), 6), dec=round(float(x['dec']), 6),
            n_xsl_epochs=len(epochs[x['star']]),
            xslid=';'.join(epochs[x['star']]),
            teff_xsl=t_x, e_teff_xsl=fnum(x['e_teff'], 1),
            teff_uves=t_u, e_teff_uves=fnum(uv['e_teff'], 1),
            d_teff=delta(t_u, t_x, 1),
            logg_xsl=g_x, e_logg_xsl=fnum(x['e_logg'], 3),
            logg_uves=g_u, e_logg_uves=fnum(uv['e_logg'], 3),
            d_logg=delta(g_u, g_x, 2),
            feh_xsl=z_x, e_feh_xsl=fnum(x['e_feh'], 3),
            feh_uves=z_u, e_feh_uves=fnum(uv['e_fe_h'], 3),
            d_feh=delta(z_u, z_x, 2),
            afe_uves=fnum(uv['a_fe'], 2),
            vsini_uves=fnum(uv['vsini'], 1), e_vsini_uves=fnum(uv['e_vsini'], 1),
            rv_uves_kms=fnum(uv['v'], 2),
            xsl_file=x['filename'], uves_spec_url=uv['spec_url'],
            notes=';'.join(note(x, uv))))

    rows.sort(key=lambda r: r['xsl_name'])
    with open('data/uves_xsl_overlap.csv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print('-> data/uves_xsl_overlap.csv')

    print(f'\n{"XSL":<12}{"UVES":<12}{"sep":>6}{"ep":>3}  {"SIMBAD":<12}'
          f'{"T_xsl":>7}{"T_uves":>8}{"dT":>7}'
          f'{"g_xsl":>7}{"g_uv":>7}{"dg":>7}'
          f'{"z_xsl":>7}{"z_uv":>7}{"dz":>7}{"vsini":>7}')
    for r in rows:
        c = {k: str(v) for k, v in r.items()}
        print(f'{c["xsl_name"]:<12}{c["uves_name"]:<12}{r["match_arcsec"]:>6.2f}'
              f'{r["n_xsl_epochs"]:>3}  {c["simbad_id"]:<12}'
              f'{c["teff_xsl"]:>7}{c["teff_uves"]:>8}{c["d_teff"]:>7}'
              f'{fmt(r["logg_xsl"], 2):>7}{fmt(r["logg_uves"], 2):>7}{c["d_logg"]:>7}'
              f'{fmt(r["feh_xsl"], 2):>7}{fmt(r["feh_uves"], 2):>7}{c["d_feh"]:>7}'
              f'{c["vsini_uves"]:>7}')

    # The positional match, confirmed a second and independent way: resolve the
    # XSL identifier too and require the two main_ids to be the same object.
    # This is what rules out a chance alignment, which 0.1" already makes
    # unlikely but does not exclude.
    tx = Simbad().query_objects([r['xsl_name'] for r in rows])
    agree = [str(a['main_id']) == r['simbad_id'] for a, r in zip(tx, rows)]
    print(f'\nIdentifier check: {sum(agree)}/{len(rows)} pairs resolve to the '
          f'same SIMBAD main_id')
    for a, r in zip(agree, rows):
        if not a:
            print(f'  DISAGREES: {r["xsl_name"]} / {r["uves_name"]}')

    d = np.array([r['d_teff'] for r in rows if r['d_teff'] != ''], float)
    print(f'\nTeff, UVES-POP minus XSL, {len(d)} stars: median {np.median(d):+.0f} K, '
          f'scatter (MAD) {np.median(np.abs(d - np.median(d))):.0f} K, '
          f'range {d.min():+.0f} to {d.max():+.0f} K')

    near = (d2d >= RADIUS) & (d2d < NEAR)
    print(f'\nNear-misses {RADIUS} to {NEAR} -- separate stars, not matched:')
    for i in np.where(near)[0]:
        print(f'  {u_rows[i]["name"]} / {xsl[idx[i]]["star"]}  {d2d[i].arcsec:.0f}"')


if __name__ == '__main__':
    sys.exit(main())
