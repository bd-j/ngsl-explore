"""Gaia DR3 photometry and XP spectra for the sample, as an independent dust lever.

Why Gaia: the break test needs E(B-V) measured, not adopted, and measured from
data that is statistically independent of the break itself. Gaia supplies that
from space, over almost exactly the NGSL range, with no slit and no atmosphere.

Two products:

  BP/RP photometry   integrated G, BP, RP, plus RUWE and the parallax.
  XP sampled spectra 336-1020 nm at 2 nm sampling, externally calibrated to
                     ~1-2%. This is the useful one -- an independent,
                     absolutely-calibrated SED across the whole optical, whose
                     long lever arm (0.036 mag of differential extinction per
                     0.01 mag of E(B-V), against 0.0035 for the 3200-3500 A
                     NGSL window) is what makes the dust term small.

TWO SOURCES, deliberately:

  VizieR (I/355/gaiadr3) for the cone search. The ESA TAP endpoint took 115 s
  per cone search and then degraded to failing three retries in a row; VizieR
  answers the same query in 1.1 s. VizieR also returns `Source` as int64, so the
  19-digit source_id survives exactly -- through ESA TAP it came back float64
  rounded (HD147550 as ...768000 against a true ...767872), which would have
  made the XP retrieval point at a different star.

  ESA DataLink for the XP spectra, because only the Gaia archive serves them.
  That part was never slow.

TWO CAVEATS, recorded per star rather than assumed away:

  * BRIGHT-STAR SYSTEMATICS. The sample runs V = 5.5-9.1. Gaia photometry and XP
    calibration degrade for the brightest stars (different gate/window schemes
    below G ~ 6, and XP is not published for the very brightest), so XPcont,
    XPsamp, G and the BP/RP excess factor are all stored and must be checked
    before an XP spectrum is believed.
  * EPOCH. Library coordinates are J2000/epoch 2000; DR3 is epoch 2016.0, and
    two of these stars are SIMBAD `PM*`. The cone is kept TIGHT (5") because a
    wide one invites picking up a neighbour; within it the counterpart is chosen
    by brightness, since for an A star G ~ V to ~0.1 mag. Observed offsets came
    out 0.1-2.7". Anything unmatched is reported, never resolved by widening.

RUWE is the better binarity test for stars this bright: it is a statement about
the astrometric fit of this star, so it cannot simply be absent from a catalog
the way an object type can. A spectroscopic binary got into this project once
already through a clean SIMBAD type (HD162630).

Writes data/gaia_sample.csv and data/gaia_xp/<star>_xp.csv
"""
import csv
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

ROOT = Path(__file__).resolve().parent.parent
XP_DIR = ROOT / 'data' / 'gaia_xp'
CATALOG = 'I/355/gaiadr3'
CONE_ARCSEC = 5.0
G_MINUS_V_MAX = 0.6         # reject a cone hit this far from the catalog V
RUWE_BINARY = 1.4

# XP sampled spectra come on a FIXED grid -- 336-1020 nm in 2 nm steps, 343
# points -- and the table carries no wavelength column, so it is reconstructed
# here and the length checked against it.
XP_NM_MIN, XP_NM_STEP, XP_N = 336.0, 2.0, 343
# XP flux is W m^-2 nm^-1; NGSL and the models are erg s^-1 cm^-2 A^-1.
# 1 W m^-2 nm^-1 = 1e7 erg/s * 1e-4 cm^-2 * 1e-1 A^-1 = 1e2 of them.
XP_FLUX_TO_CGS = 100.0

VIZ_COLS = ['Source', 'RA_ICRS', 'DE_ICRS', 'Gmag', 'BPmag', 'RPmag',
            'e_Gmag', 'e_BPmag', 'e_RPmag', 'BP-RP', 'RUWE', 'Plx', 'e_Plx',
            'pmRA', 'pmDE', 'NSS', 'XPcont', 'XPsamp', 'E(BP/RP)']


def fnum(x):
    try:
        v = float(x)
        return v if np.isfinite(v) else np.nan
    except (TypeError, ValueError):
        return np.nan


def col(row, name):
    """Value of a column that may be absent from the VizieR response."""
    try:
        return fnum(row[name])
    except (KeyError, TypeError, ValueError):
        return np.nan


def yes(row, name):
    v = col(row, name)
    return bool(np.isfinite(v) and v > 0)


def as_array(cell):
    """Array from an XP table cell.

    With retrieval_type='XP_SAMPLED' and format='csv', astroquery hands back ONE
    row whose flux cell is the literal string '(1.2E-14, 1.1E-14, ...)' rather
    than an array column. Accept either, so this keeps working if that changes.
    """
    if isinstance(cell, (str, np.str_)):
        return np.array([float(v) for v in
                         str(cell).strip().strip('()').split(',') if v.strip()])
    return np.atleast_1d(np.asarray(cell, dtype=float))


def sep_arcsec_arr(ra, dec, ras, decs):
    """Great-circle separation of one position from an array, in arcsec."""
    r1, d1 = np.radians(ra), np.radians(dec)
    r2, d2 = np.radians(ras), np.radians(decs)
    c = np.sin(d1) * np.sin(d2) + np.cos(d1) * np.cos(d2) * np.cos(r1 - r2)
    return np.degrees(np.arccos(np.clip(c, -1.0, 1.0))) * 3600.0


def photometry(sample):
    """Cone-search VizieR's Gaia DR3 for each star -> list of records."""
    import astropy.units as u
    from astropy.coordinates import SkyCoord
    from astroquery.vizier import Vizier

    viz = Vizier(row_limit=-1, columns=VIZ_COLS)
    print(f'\n{"star":<11}{"V":>6}{"G":>7}{"BP":>7}{"RP":>7}{"RUWE":>6}'
          f'{"plx":>7}{"d(pc)":>7}{"sep":>6}  xp  flags')

    rows = []
    for s in sample:
        ra, dec, vmag = fnum(s['ra']), fnum(s['dec']), fnum(s['vmag'])
        try:
            res = viz.query_region(SkyCoord(ra, dec, unit='deg'),
                                   radius=CONE_ARCSEC * u.arcsec,
                                   catalog=CATALOG)
        except Exception as exc:
            print(f'{s["star"]:<11}  QUERY FAILED: {type(exc).__name__} '
                  f'{str(exc)[:90]}')
            continue
        if not len(res):
            print(f'{s["star"]:<11}  no Gaia source within {CONE_ARCSEC:.0f}"')
            continue
        t = res[0]
        seps = sep_arcsec_arr(ra, dec,
                             np.array([fnum(v) for v in t['RA_ICRS']]),
                             np.array([fnum(v) for v in t['DE_ICRS']]))

        # brightness, not position, picks the counterpart: for an A star G ~ V
        g = np.array([fnum(v) for v in t['Gmag']])
        ok = np.isfinite(g) & (np.abs(g - vmag) < G_MINUS_V_MAX)
        if not ok.any():
            print(f'{s["star"]:<11}  no cone hit within {G_MINUS_V_MAX} mag of '
                  f'V={vmag:.2f}; nearest G={np.nanmin(g):.2f}')
            continue
        i = int(np.where(ok)[0][np.argmin(np.abs(g[ok] - vmag))])
        r = t[i]

        plx, plxe = col(r, 'Plx'), col(r, 'e_Plx')
        ruwe = col(r, 'RUWE')
        nss = col(r, 'NSS')
        flags = []
        if np.isfinite(ruwe) and ruwe > RUWE_BINARY:
            flags.append(f'RUWE={ruwe:.2f}>{RUWE_BINARY}')
        if np.isfinite(nss) and nss > 0:
            flags.append(f'NSS={nss:.0f}')
        if vmag < 6.0:
            flags.append('bright: check XP calibration')

        rec = dict(star=s['star'], source_id=str(int(r['Source'])),
                   sep_arcsec=round(float(seps[i]), 2), vmag=vmag,
                   gmag=col(r, 'Gmag'), bpmag=col(r, 'BPmag'),
                   rpmag=col(r, 'RPmag'), e_gmag=col(r, 'e_Gmag'),
                   e_bpmag=col(r, 'e_BPmag'), e_rpmag=col(r, 'e_RPmag'),
                   bp_rp=col(r, 'BP-RP'), bp_rp_excess=col(r, 'E(BP/RP)'),
                   parallax=plx, parallax_error=plxe,
                   dist_pc=(round(1000.0 / plx, 1)
                            if np.isfinite(plx) and plx > 0 else ''),
                   pmra=col(r, 'pmRA'), pmdec=col(r, 'pmDE'),
                   ruwe=ruwe,
                   non_single_star=('' if not np.isfinite(nss)
                                    else f'{nss:.0f}'),
                   has_xp_continuous='yes' if yes(r, 'XPcont') else 'no',
                   has_xp_sampled='yes' if yes(r, 'XPsamp') else 'no',
                   ruwe_binary='yes' if (np.isfinite(ruwe)
                                         and ruwe > RUWE_BINARY) else 'no',
                   gaia_flags='; '.join(flags))
        rows.append(rec)
        print(f'{rec["star"]:<11}{vmag:>6.2f}{rec["gmag"]:>7.3f}'
              f'{rec["bpmag"]:>7.3f}{rec["rpmag"]:>7.3f}{ruwe:>6.2f}'
              f'{plx:>7.2f}{str(rec["dist_pc"]):>7}{rec["sep_arcsec"]:>6.1f}  '
              f'{"Y" if rec["has_xp_continuous"] == "yes" else "-"}'
              f'{"Y" if rec["has_xp_sampled"] == "yes" else "-"}'
              f'  {rec["gaia_flags"]}')
    return rows


def xp_spectra(rows):
    """Retrieve XP sampled spectra from the ESA archive -> count written."""
    from astroquery.gaia import Gaia

    want = [r for r in rows if r['has_xp_sampled'] == 'yes']
    print(f'\n{len(want)} stars have XP sampled; retrieving from ESA DataLink')
    XP_DIR.mkdir(parents=True, exist_ok=True)
    got = 0
    for r in want:
        out = XP_DIR / f'{r["star"]}_xp.csv'
        try:
            d = Gaia.load_data(ids=[int(r['source_id'])],
                               data_release='Gaia DR3',
                               retrieval_type='XP_SAMPLED',
                               data_structure='raw', format='csv')
        except Exception as exc:
            print(f'  {r["star"]:<11} FAILED {type(exc).__name__}: '
                  f'{str(exc)[:110]}')
            continue
        if not d:
            print(f'  {r["star"]:<11} no payload returned')
            continue
        key = next(iter(d))
        tab = d[key][0] if isinstance(d[key], list) else d[key]

        # The payload is ONE row whose flux/flux_error cells hold the whole
        # array. Written straight out that gives a 1-row file with the arrays
        # stringified and no wavelength axis at all -- so unpack explicitly.
        flux, ferr = as_array(tab['flux'][0]), as_array(tab['flux_error'][0])
        if flux.size != XP_N:
            print(f'  {r["star"]:<11} UNEXPECTED length {flux.size} '
                  f'(expected {XP_N}); skipped')
            continue
        got_id = str(tab['source_id'][0]).strip()
        if got_id != r['source_id']:
            print(f'  {r["star"]:<11} SOURCE MISMATCH: asked '
                  f'{r["source_id"]}, got {got_id}; skipped')
            continue
        # Gaia is a space instrument with no air path, so these are VACUUM
        # wavelengths -- unlike NGSL, XSL and UVES-POP, no conversion needed.
        wave_A = (XP_NM_MIN + XP_NM_STEP * np.arange(XP_N)) * 10.0
        with open(out, 'w', newline='') as fh:
            wr = csv.writer(fh)
            wr.writerow(['wavelength_A_vacuum', 'flam', 'flam_err'])
            for wv, fl, fe in zip(wave_A, flux * XP_FLUX_TO_CGS,
                                  ferr * XP_FLUX_TO_CGS):
                wr.writerow([f'{wv:.1f}', f'{fl:.6e}', f'{fe:.6e}'])
        got += 1
        print(f'  {r["star"]:<11} {flux.size:>4} pts  '
              f'{wave_A[0]:.0f}-{wave_A[-1]:.0f} A  '
              f'median S/N {np.median(flux / ferr):>5.0f}  '
              f'-> {out.relative_to(ROOT)}')
    print(f'\n{got}/{len(want)} XP spectra -> data/gaia_xp/')
    return got


def main():
    sample = list(csv.DictReader(open(ROOT / 'data' / 'sample.csv')))
    print(f'{len(sample)} stars from data/sample.csv; '
          f'VizieR {CATALOG}, {CONE_ARCSEC:.0f}" cone')

    rows = photometry(sample)
    if not rows:
        raise SystemExit('no Gaia matches at all -- nothing written')
    with open(ROOT / 'data' / 'gaia_sample.csv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(f'\n{len(rows)}/{len(sample)} matched -> data/gaia_sample.csv')
    missing = sorted({s['star'] for s in sample} - {r['star'] for r in rows})
    if missing:
        print(f'  UNMATCHED: {", ".join(missing)}')

    xp_spectra(rows)


if __name__ == '__main__':
    main()
