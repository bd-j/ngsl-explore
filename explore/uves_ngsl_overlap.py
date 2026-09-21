"""Stars observed by both NGSL v2 and UVES-POP, with the parameters from each.

Neither catalog carries the other's identifiers, and UVES-POP names its brightest
targets by proper name (Achernar, Arcturus, ...), so the match is positional:
SIMBAD resolves the 406 UVES-POP names, and those coordinates are matched
against the NGSL v2 FITS-header positions.

Parameter provenance -- the two sets are independent and were fitted very
differently, which is the point of putting them side by side:

  NGSL     Castelli (2004) models on Victoria-Regina isochrones, fitted to the
           NGSL v1 spectra (aaareadme.txt section 4). No uncertainties quoted;
           Teff is printed to 1 K but is not good to that.
  UVES-POP VOXAstro-SL re-reduction (2023ApJS..266...11B), full-spectrum fit
           against PHOENIX, with formal errors and v sin i. The quoted errors
           are formal fit errors only -- e_teff runs to ~2 K -- and are not an
           uncertainty budget.

The library's `v` column is the radial velocity in km/s, NOT the V magnitude.
The check is in the code below and runs every time: across the 12 stars here
with a value it tracks the SIMBAD RV to 3.7 km/s and misses V by up to 114 mag.
It is carried through as rv_uves_kms.

Writes data/uves_ngsl_overlap.csv

Caveat: UVES-POP is a southern library (VLT) and NGSL is all-sky, so the
overlap is small -- 13 stars, none of them A stars in the Balmer window.
"""
import csv
import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord
from astroquery.simbad import Simbad

# 13 pairs fall below 8", the next is the 36" visual pair HD36959/HD36960 --
# two different stars -- and the one after that is at 1952". Any cut in that
# gap gives the same 13, so the answer does not depend on where the line goes.
#
# 10" rather than the 5" build_sample.py uses for XSL, because the residual
# separations here are proper motion between the NGSL header epoch and SIMBAD's
# J2000 positions, and UVES-POP's targets are nearby bright stars. The four
# widest matches are exactly the four fastest movers -- 171 Pup 1.72"/yr at
# 7.33", 61 Vir 1.51"/yr at 3.20", eps Eri 0.98"/yr at 1.75", 10 Tau 0.54"/yr
# at 2.11" -- which is what confirms them as epoch offsets and not mismatches.
# At 5" 171 Pup would be dropped for moving.
RADIUS = 10 * u.arcsec
# Near-misses out to here are printed so the rejected ones stay visible.
NEAR = 120 * u.arcsec


def fnum(x, nd=None):
    """Catalog string -> rounded float, or '' when absent."""
    try:
        v = float(x)
    except (TypeError, ValueError):
        return ''
    if not np.isfinite(v):
        return ''
    return round(v, nd) if nd is not None else v


def note(ng, uv):
    """Why a row's parameters may not be comparable. Recorded, not eyeballed."""
    n = []
    if not fnum(ng['teff']):
        n.append(f'ngsl_unfitted(fit_quality={ng["fit_quality"]})')
    if not fnum(uv['teff']):
        n.append('uves_unfitted')
    g = fnum(uv['logg'])
    if g != '' and g < 0.5:
        n.append(f'uves_logg={g:.2f}_at_grid_edge')
    if ng['dataqual'] != 'good':
        n.append(f'ngsl_dataqual={ng["dataqual"]}')
    return n


ngsl = list(csv.DictReader(open('data/ngsl_catalog.csv')))
uves = list(csv.DictReader(open('data/uves_pop_all.csv')))
print(f'NGSL {len(ngsl)} stars, UVES-POP {len(uves)} stars')

ng_c = SkyCoord([float(r['ra']) for r in ngsl] * u.deg,
                [float(r['dec']) for r in ngsl] * u.deg)

print('Resolving UVES-POP names in SIMBAD...')
sim = Simbad().query_objects([r['name'] for r in uves])
u_ra = np.asarray(sim['ra'], float)
u_dec = np.asarray(sim['dec'], float)
u_id = [str(x) for x in sim['main_id']]
ok = np.isfinite(u_ra) & np.isfinite(u_dec)
print(f'  resolved {ok.sum()}/{len(uves)}')

# The 49 that do not resolve are open-cluster members under the library's own
# running numbers (46 in IC 2391, 3 in NGC 6475). Neither field is reachable:
# the nearest NGSL star to IC 2391 is 6.4 deg away, to NGC 6475 12.2 deg, so
# none of them can be an NGSL target and dropping them costs nothing.
unresolved = [r['name'] for r, k in zip(uves, ok) if not k]
assert all(n.startswith(('IC2391', 'NGC6475')) for n in unresolved), unresolved
print(f'  unresolved (cluster running numbers, both fields >6 deg from any '
      f'NGSL star): {len(unresolved)}')

u_c = SkyCoord(u_ra[ok] * u.deg, u_dec[ok] * u.deg)
u_rows = [r for r, k in zip(uves, ok) if k]
u_ids = [i for i, k in zip(u_id, ok) if k]

idx, d2d, _ = u_c.match_to_catalog_sky(ng_c)
hit = d2d < RADIUS
print(f'\nMatched within {RADIUS}: {hit.sum()} stars')

rows = []
for i in np.where(hit)[0]:
    uv, ng = u_rows[i], ngsl[idx[i]]
    t_ng, t_uv = fnum(ng['teff']), fnum(uv['teff'])
    g_ng, g_uv = fnum(ng['logg']), fnum(uv['logg'])
    z_ng, z_uv = fnum(ng['logz']), fnum(uv['fe_h'])
    rows.append(dict(
        ngsl_target=ng['target'], uves_name=uv['name'], simbad_id=u_ids[i],
        match_arcsec=round(float(d2d[i].arcsec), 2),
        ra=round(float(ng['ra']), 6), dec=round(float(ng['dec']), 6),
        sptype=ng['sptype'], vmag=fnum(ng['vmag']), bmag=fnum(ng['bmag']),
        teff_ngsl=t_ng, teff_uves=fnum(uv['teff'], 1),
        e_teff_uves=fnum(uv['e_teff'], 1),
        logg_ngsl=g_ng, logg_uves=fnum(uv['logg'], 3),
        e_logg_uves=fnum(uv['e_logg'], 3),
        mh_ngsl=z_ng, feh_uves=fnum(uv['fe_h'], 3),
        e_feh_uves=fnum(uv['e_fe_h'], 3),
        afe_uves=fnum(uv['a_fe'], 2),
        vsini_uves=fnum(uv['vsini'], 1), e_vsini_uves=fnum(uv['e_vsini'], 1),
        rv_uves_kms=fnum(uv['v'], 2),
        d_teff=round(t_uv - t_ng, 1) if '' not in (t_ng, t_uv) else '',
        d_logg=round(g_uv - g_ng, 2) if '' not in (g_ng, g_uv) else '',
        d_mh=round(z_uv - z_ng, 2) if '' not in (z_ng, z_uv) else '',
        ngsl_file=ng['file'], ngsl_dataqual=ng['dataqual'],
        ngsl_fit_quality=ng['fit_quality'], uves_spec_url=uv['spec_url'],
        notes=';'.join(note(ng, uv))))

rows.sort(key=lambda r: r['ngsl_target'])
with open('data/uves_ngsl_overlap.csv', 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0]))
    w.writeheader()
    w.writerows(rows)
print('-> data/uves_ngsl_overlap.csv')

print(f'\n{"NGSL":<10}{"UVES":<10}{"sep":>6}  {"SpT":<8}'
      f'{"Teff N":>8}{"Teff U":>9}{"dT":>8}'
      f'{"logg N":>8}{"logg U":>8}{"[M/H] N":>9}{"[Fe/H] U":>9}{"vsini":>7}')
for r in rows:
    c = {k: str(v) for k, v in r.items()}
    print(f'{c["ngsl_target"]:<10}{c["uves_name"]:<10}{r["match_arcsec"]:>6.2f}  '
          f'{c["sptype"]:<8}{c["teff_ngsl"]:>8}{c["teff_uves"]:>9}'
          f'{c["d_teff"]:>8}{c["logg_ngsl"]:>8}{c["logg_uves"]:>8}'
          f'{c["mh_ngsl"]:>9}{c["feh_uves"]:>9}{c["vsini_uves"]:>7}')

# The `v` column: confirm against SIMBAD that it is RV, not the V magnitude.
sv = Simbad()
sv.add_votable_fields('V', 'rvz_radvel')
chk = sv.query_objects([r['ngsl_target'] for r in rows])
drv, dv = [], []
for r, c in zip(rows, chk):
    if r['rv_uves_kms'] == '':
        continue
    drv.append(abs(r['rv_uves_kms'] - float(c['rvz_radvel'])))
    dv.append(abs(r['rv_uves_kms'] - float(c['V'])))
print(f"\nUVES-POP 'v' column vs SIMBAD, {len(drv)} stars: "
      f'max|v - RV| = {max(drv):.1f} km/s, max|v - V| = {max(dv):.0f} mag '
      f'-> it is the radial velocity in km/s, not the V magnitude')

near = (d2d >= RADIUS) & (d2d < NEAR)
print(f'\nNear-misses {RADIUS} to {NEAR} -- separate stars, not matched:')
for i in np.where(near)[0]:
    print(f'  {u_rows[i]["name"]} / {ngsl[idx[i]]["target"]}  '
          f'{d2d[i].arcsec:.1f}" -- distinct components of a visual pair')
