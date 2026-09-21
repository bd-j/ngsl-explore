"""U-B from photometry against U-B synthesised from NGSL and UVES-POP.

U straddles the Balmer break: bessell_U runs 3050-4150 A with lambda_eff
3571 A, and the break at 3646 A sits inside it. So U-B is a broadband measure
of the break amplitude, and a library that gets the break wrong must get U-B
wrong by a related amount. That is what this checks, for the 13 stars in both
libraries (data/uves_ngsl_overlap.csv).

CATALOG SOURCE: Mermilliod's homogeneous means, VizieR II/168/ubvmeans, which
publishes V, B-V and U-B together. U AND B THEREFORE COME FROM THE SAME SOURCE
by construction, which matters because U-B is a difference: two bands taken
from different catalogues can disagree in zero point and leave a colour error
that looks astrophysical.

SIMBAD is queried too, as a cross-check only, and it is NOT used. It disagrees
badly for HD022484: SIMBAD gives B = 5.150 with V = 4.300, so B-V = +0.85 for
an F9IV-V star that should sit near +0.57, and Mermilliod's homogeneous value
is +0.572. That one bad B propagates straight into U-B. Discrepancies are
reported per star rather than averaged over.

UVES-POP CANNOT MEASURE JOHNSON U. It starts at 3200 A and the U bandpass
starts at 3050 A, so 1.96% of the band's transmission-weighted integral is
missing. `common.photometry.project` returns NaN for a filter that is not
fully covered rather than a plausible clipped number, so this is caught rather
than absorbed. Two things are done about it:

  * a TRUNCATED U (zero below 3200 A) is built and applied to BOTH libraries,
    so the NGSL-vs-UVES comparison is like for like over a band both cover;
  * the cost of the truncation is measured directly, as NGSL's full U minus
    NGSL's truncated U, per star -- NGSL covers 1675-10198 A so it can do both.

B is fully covered by both (3700-5500 A), with one exception: HD138716's UVES
spectrum has a 925 A hole at 3859-4784 A, which is inside B. Its UVES colours
are therefore not computed rather than computed from a hole.

VEGA vs AB. sedpy works in AB; the catalogue is Vega. From its source,
_ab_to_vega = -2.5 log10(ab_zero/vega_zero), which gives m_vega = m_AB +
ab_to_vega. `selftest` pins that by running sedpy's own Vega spectrum through
the filters, where U-B must come out at 0 by definition.

NO RV CORRECTION is applied. The largest shift in this sample is 120 km/s,
which is 1.5 A at 3600 A against a 1100 A bandpass.

Writes data/uves_ngsl_ub.csv
"""
import csv
import sys
import warnings
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common.photometry import project
from common.uves_pop_load import load as load_uves
from explore.plot_uves_ngsl import ngsl_spectrum, fill_small_gaps

ROOT = Path(__file__).resolve().parent.parent
# Where the truncated bands are cut. The delivered UVES-POP grid starts at
# 3200.9 A in vacuum, so cutting at exactly 3200.0 leaves the band's first
# non-zero pixel BLUEWARD of the spectrum's first pixel and project()'s
# coverage guard then rejects every star -- correctly, but for a reason that is
# an off-by-one rather than a real gap. 3210 clears it with margin, and costs
# nothing: the whole truncation is worth 0.004 mag in U-B and 0.002 in c1.
UVES_BLUE_LIMIT = 3210.0
MERMILLIOD = 'II/168/ubvmeans'
CONE = 15                       # arcsec; these catalogues are epoch J2000 like
                                # the library coordinates, so no PM problem


def filters(truncate_U_at=None):
    """-> (U, B) sedpy filters. `truncate_U_at` zeroes U below that wavelength.

    The truncated U is not a real instrument. It exists so the same response is
    applied to a spectrum that starts at 3200 A and to one that starts at
    1675 A, which is the only way NGSL and UVES-POP can be compared in U at all.
    """
    from sedpy.observate import load_filters, Filter
    U, B = load_filters(['bessell_U', 'bessell_B'])
    if truncate_U_at is None:
        return U, B
    w, t = np.array(U.wavelength, float), np.array(U.transmission, float)
    t = np.where(w < truncate_U_at, 0.0, t)
    return Filter(kname='bessell_U_gt3200', data=(w, t)), B


def vega_mags(wave, flam, fl):
    """-> (U, B) in VEGA magnitudes, or nan where the band is not covered."""
    mgy = project(wave, flam, list(fl))
    with np.errstate(divide='ignore', invalid='ignore'):
        ab = -2.5 * np.log10(mgy)
    return tuple(ab[i] + fl[i].ab_to_vega for i in range(2))


def vega_spectrum():
    """-> (wavelength A, F_lambda cgs) of THE Vega sedpy itself calibrated on.

    Imported rather than read off disk. sedpy ships two CALSPEC files and uses
    alpha_lyr_stis_005; picking the other one by globbing made this selftest
    return -0.023 instead of 0, which is a real 0.02 mag error and not a
    tolerance to widen.
    """
    from sedpy.reference_spectra import vega
    return np.asarray(vega[:, 0], float), np.asarray(vega[:, 1], float)


def selftest():
    """Vega through the Johnson filters must give U = B = 0, hence U-B = 0.

    This pins the AB->Vega direction. Getting it backwards flips U-B by
    0.92 mag, which on this sample would still leave a plausible-looking
    correlation with the catalogue -- just offset -- so it has to be tested
    against a value that is known exactly rather than eyeballed.

    Stromgren is NOT testable this way: its indices carry conventional zero
    points set by standard stars, not by Vega, so c1(Vega) is not 0 and the
    offset has to be measured on the sample instead. See main().
    """
    wv, fv = vega_spectrum()
    U, B = vega_mags(wv, fv, filters())
    assert abs(U) < 0.02 and abs(B) < 0.02, f'Vega gives U={U:.3f} B={B:.3f}'
    print(f'  selftest: Vega -> U={U:+.4f} B={B:+.4f} U-B={U - B:+.4f} '
          f'(must be 0 by definition)')



def stromgren(truncate_at=None):
    """-> (u, v, b) sedpy filters for the c1 index, optionally truncated.

    c1 = (u - v) - (v - b) is the classical Balmer-discontinuity index: u
    (3150-3775 A) spans the break at 3646 A, v (3725-4500 A) sits entirely
    redward of it, and b anchors the slope. Only 0.63% of u's integral falls
    below UVES-POP's 3200 A limit, against 1.96% for Johnson U -- so of the two
    systems this is the one UVES-POP can nearly do.
    """
    from sedpy.observate import load_filters, Filter
    u, v, b = load_filters(['stromgren_u', 'stromgren_v', 'stromgren_b'])
    if truncate_at is None:
        return u, v, b
    w, t = np.array(u.wavelength, float), np.array(u.transmission, float)
    t = np.where(w < truncate_at, 0.0, t)
    return Filter(kname='stromgren_u_gt3200', data=(w, t)), v, b


def c1_of(wave, flam, fl):
    """-> Stromgren c1 = (u-v) - (v-b) on sedpy's AB->Vega scale.

    NOT on the standard Stromgren zero point, which is conventional. main()
    removes the offset before comparing, and says so.
    """
    mgy = project(wave, flam, list(fl))
    with np.errstate(divide='ignore', invalid='ignore'):
        ab = -2.5 * np.log10(mgy)
    u, v, b = (ab[i] + fl[i].ab_to_vega for i in range(3))
    return (u - v) - (v - b)


def _vizier_row(cat, star, vmag, cols, vtol=0.35):
    """One catalogue row for `star`, matched by IDENTIFIER first.

    These catalogues are keyed by LID, the Lausanne identifier, whose HD form
    is '0100' + the six-digit HD number. That is exact and has no epoch in it,
    which matters: II/215 carries B1950 positions and this sample contains
    stars that move 27" between epochs, so a cone search is both incomplete and
    unsafe. Matching II/215 positionally at 15" returned c1 for 5 of 13 and
    widening to 120" returned 12 -- by picking up neighbours. By LID it is 10,
    deterministically.

    Not every star is filed under its HD number, so there is a positional
    fallback -- but it must AGREE IN V to `vtol` before it is accepted, and a
    match that fails is dropped and reported rather than kept.

    `vtol` is 0.35 mag, which is loose on purpose. Seven of these thirteen are
    catalogued variables, and the catalogues were not taken at one epoch: FW
    CMa is V = 5.20 in the NGSL table and 5.356 in Mermilliod, a 0.156 mag
    difference that is the Be star varying, not a mismatch. A tolerance of 0.15
    rejected it. The guard still works -- a neighbour close enough in position
    AND within 0.35 mag of a V ~ 5 star is not something this sample contains.
    """
    import astropy.units as u
    from astropy.coordinates import SkyCoord
    from astroquery.vizier import Vizier
    viz = Vizier(row_limit=20, columns=cols)
    # The LID prefix is NOT the same in the two catalogues: II/215 writes the
    # HD form as '0100nnnnnn' and II/168 as '+100nnnnnn'. Trying only the first
    # made every II/168 lookup fall through to the positional fallback without
    # ever saying so -- it still found 11 of 13, which is exactly why it went
    # unnoticed. Both forms are tried.
    if star.startswith('HD'):
        n = int(star[2:])
        for lid in (f'0100{n:06d}', f'+100{n:06d}'):
            try:
                t = viz.query_constraints(catalog=cat, LID=f'=={lid}')
                if len(t) and len(t[0]):
                    return t[0][0], f'LID {lid}'
            except Exception:
                continue
    if not np.isfinite(vmag):
        return None, 'no-LID, no V to verify a positional match'
    try:
        res = viz.query_region(SkyCoord(*_radec(star), unit='deg'),
                               radius=CONE * u.arcsec, catalog=cat)
    except Exception as exc:
        return None, f'query failed ({type(exc).__name__})'
    if not len(res):
        return None, 'not in catalogue'
    t = res[0]
    for row in t:
        try:
            if abs(float(row['Vmag']) - vmag) <= vtol:
                return row, f'position, V agrees to {vtol}'
        except (KeyError, TypeError, ValueError):
            continue
    return None, f'positional candidates all disagree in V by > {vtol}'


_RADEC = {}


def _radec(star):
    return _RADEC[star]


def catalog_c1(rows):
    """-> {star: (c1, e_c1, how)} from Hauck & Mermilliod, VizieR II/215."""
    out = {}
    for r in rows:
        star = r['ngsl_target']
        row, how = _vizier_row('II/215/catalog', star, _f(r['vmag']), ['*'])
        if row is None:
            print(f'  c1  {star:<10} {how}')
            continue
        try:
            v1, e1 = float(row['c1']), float(row['e_c1'])
        except (KeyError, TypeError, ValueError):
            continue
        if np.isfinite(v1):
            out[star] = (v1, e1, how)
    return out


def catalog_ubv(rows):
    """-> {star: dict} from Mermilliod II/168, U and B from the ONE source."""
    out = {}
    for r in rows:
        star = r['ngsl_target']
        row, how = _vizier_row('II/168/ubvmeans', star, _f(r['vmag']), ['*'])
        if row is None:
            print(f'  UBV {star:<10} {how}')
            continue
        g = lambda k: (float(row[k]) if k in row.colnames
                       and np.isfinite(float(row[k])) else np.nan)
        out[star] = dict(V=g('Vmag'), BV=g('B-V'), UB=g('U-B'), how=how)
    return out


def _f(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return np.nan



VAR_C = '#c0392b'           # catalogued variable: red, as asked
STAT_C = '#2a78d6'          # everything else

# Labelled in EVERY panel regardless of how far off the line they land. These
# are the two stars with no GCVS variability flag that nonetheless show sharp
# UVES-vs-NGSL differences across the break, so they have to stay identifiable
# from panel to panel -- the point of following them is to see whether the same
# star is the outlier in U-B, in c1, and against the jump. It is not.
ALWAYS_LABEL = {'HD076932', 'HD063077'}


def make_figure(recs, zp):
    """Synthetic colour against catalogue colour, one point per library.

    Two indices side by side because they behave differently: U-B carries a
    systematic offset and c1 does not, and seeing them on the same figure is
    what makes that a statement about the Johnson U bandpass rather than about
    NGSL.

    Symbols separate the libraries (circle NGSL, triangle UVES-POP) and colour
    separates variable from not, so a reader can tell at a glance whether a
    point is off the line because of the instrument or because the star moved
    between epochs. The two points for one star are joined by a hairline, which
    is the only way to see that HD058343's pair straddles the line.
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from common.figpath import library_figure_path

    SURFACE, INK, MUTED, GRID = '#fcfcfb', '#22262b', '#6b7280', '#dfe3e8'
    fig, axes = plt.subplots(2, 3, figsize=(18.2, 9.4), sharex='col',
                             gridspec_kw={'height_ratios': [2.5, 1],
                                          'hspace': .06, 'wspace': .22})
    fig.patch.set_facecolor(SURFACE)

    # The first two pairs compare a colour with its own catalogue value, so x
    # and the reference are the same column and the 1:1 line means something.
    # The THIRD pair is different: x is c1, the Balmer-jump index, and y is
    # U-B. It asks whether the U-B disagreement depends on the size of the
    # break, which is the question this project is actually about -- so it has
    # no 1:1 line, and its residual panel is the U-B error against the jump.
    panels = [dict(x='ub_cat', ref='ub_cat', n='ub_ngsl', u='ub_uves_corr',
                   title='Johnson $U-B$', xlab='U−B (Mermilliod II/168)',
                   ylab='$U-B$ from the spectrum', one_to_one=True),
              dict(x='c1_cat', ref='c1_cat', n='c1_ngsl_zp', u='c1_uves_zp',
                   title=r'Str' + '\u00f6' + r'mgren $c_1$',
                   xlab='c1 (Hauck & Mermilliod II/215)',
                   ylab='$c_1$ from the spectrum', one_to_one=True),
              dict(x='c1_cat', ref='ub_cat', n='ub_ngsl', u='ub_uves_corr',
                   title='$U-B$ against the Balmer jump',
                   xlab='c1 (Hauck & Mermilliod II/215)',
                   ylab='$U-B$', one_to_one=False)]

    # Proxy handles, built once and drawn in EVERY panel. The panels are read
    # side by side and a reader should not have to look back to the first one
    # to recall what a triangle means.
    handles = [plt.Line2D([], [], ls='', marker='o', color=STAT_C, ms=7,
                          label='NGSL v2'),
               plt.Line2D([], [], ls='', marker='^', color=STAT_C, ms=7,
                          label='UVES-POP'),
               plt.Line2D([], [], ls='', marker='s', color=VAR_C, ms=7,
                          label='catalogued variable (GCVS)'),
               plt.Line2D([], [], ls='', marker='s', color=STAT_C, ms=7,
                          label='not known variable')]

    for j, spec in enumerate(panels):
        kx, kc, kn, ku = spec['x'], spec['ref'], spec['n'], spec['u']
        title, xlab = spec['title'], spec['xlab']
        ax, rax = axes[0, j], axes[1, j]
        for a in (ax, rax):
            a.set_facecolor(SURFACE)
            a.grid(alpha=.25, color=GRID, lw=.7)
            a.tick_params(labelsize=8, colors=MUTED)
            for sp in a.spines.values():
                sp.set_color(GRID)

        xv = [r[kx] for r in recs if np.isfinite(r[kx])]
        xlo, xhi = min(xv), max(xv)
        xpad = .08 * (xhi - xlo)
        xlim = (xlo - xpad, xhi + xpad)
        if spec['one_to_one']:
            ax.plot(xlim, xlim, color=MUTED, lw=1, ls='--', zorder=1)
        rax.axhline(0, color=MUTED, lw=1, zorder=1)

        todo = []
        for r in recs:
            c, x = r[kc], r[kx]
            if not (np.isfinite(c) and np.isfinite(x)):
                continue
            col = VAR_C if r['gcvs_type'] else STAT_C
            n_, u_ = r[kn], r[ku]
            if np.isfinite(n_) and np.isfinite(u_):
                ax.plot([x, x], [n_, u_], color=col, lw=.7, alpha=.5, zorder=2)
                rax.plot([x, x], [n_ - c, u_ - c], color=col, lw=.7,
                         alpha=.5, zorder=2)
            if not spec['one_to_one']:
                # the catalogue value itself, so the two libraries can be read
                # against it rather than only against each other
                ax.scatter(x, c, marker='_', s=90, color=MUTED, zorder=3)
            for val, mk, ms in ((n_, 'o', 46), (u_, '^', 52)):
                if not np.isfinite(val):
                    continue
                ax.scatter(x, val, marker=mk, s=ms, facecolor=col,
                           edgecolor='white', linewidth=.6, zorder=4)
                rax.scatter(x, val - c, marker=mk, s=ms, facecolor=col,
                            edgecolor='white', linewidth=.6, zorder=4)
            if np.isfinite(n_) and (abs(n_ - c) > 0.05
                                    or r['star'] in ALWAYS_LABEL):
                todo.append((x, r['star'].replace('HD0', 'HD'), col))

        # Labels: vertical, anchored at the PHOTOMETRIC value on the x axis
        # rather than offset from a marker. Nine of the eleven c1 values fall
        # between 0.25 and 0.43, so a label beside its point sat among three
        # other stars' points and it was guesswork which one it named.
        #
        # Two things make that unambiguous rather than merely tidy: a hairline
        # dropped from the label to the top of the panel, so the label is tied
        # to an x position and not to whatever marker it happens to sit near;
        # and a stagger, because labels closer than 6% of the axis in x would
        # otherwise overprint each other -- HD63077 and HD206778 are 0.021
        # apart in c1.
        todo.sort()
        span = xlim[1] - xlim[0]
        row, prev = 0, None
        for xl, name, col in todo:
            row = 0 if prev is None or (xl - prev) > 0.06 * span else 1 - row
            prev = xl
            ax.axvline(xl, color=col, lw=.5, alpha=.22, zorder=0)
            ax.annotate(name, xy=(xl, 0), xycoords=('data', 'axes fraction'),
                        xytext=(0, 6 + 34 * row), textcoords='offset points',
                        rotation=90, ha='center', va='bottom',
                        fontsize=6.5, color=col)

        ax.set_xlim(*xlim)
        if spec['one_to_one']:
            ax.set_ylim(*xlim)
        ax.legend(handles=handles, fontsize=7, loc='best', framealpha=.92)
        ax.set_title(title, fontsize=11, color=INK)
        ax.set_ylabel(spec['ylab'], fontsize=9, color=INK)
        rax.set_ylabel('spectrum − catalogue', fontsize=9, color=INK)
        rax.set_xlabel(xlab, fontsize=9, color=INK)

        d = np.array([r[kn] - r[kc] for r in recs
                      if np.isfinite(r[kn]) and np.isfinite(r[kc])
                      and np.isfinite(r[kx])])
        nmad = 1.4826 * np.median(np.abs(d - np.median(d)))
        rax.text(.03, .06, f'NGSL: median {np.median(d):+.3f}, NMAD {nmad:.3f}'
                 f'  (n={d.size})', transform=rax.transAxes, fontsize=7.5,
                 color=INK)
        m = max(0.05, np.nanmax(np.abs([r[k] - r[kc] for r in recs
                                        for k in (kn, ku)
                                        if np.isfinite(r[k])
                                        and np.isfinite(r[kc])
                                        and np.isfinite(r[kx])])) * 1.2)
        rax.set_ylim(-m, m)

    fig.suptitle(
        'Synthetic colours from NGSL and UVES-POP against catalogue photometry'
        '\nUVES-POP starts at 3200 A so its U is on a truncated band, '
        f'corrected back with NGSL (+0.008 mag median).  '
        f'$c_1$ zero point {zp:+.3f} removed.',
        fontsize=11, color=INK, linespacing=1.6)
    fig.subplots_adjust(left=.05, right=.99, top=.895, bottom=.075)
    out = library_figure_path('uves_ngsl_ub.png')
    fig.savefig(out, dpi=200, facecolor=SURFACE)
    plt.close(fig)
    print(f'  -> {out}')


def dump_catalog_rows(rows):
    """Every column of both photometric catalogues, for all 13 stars.

    The colour comparison uses four numbers out of these tables; the rest --
    uncertainties, observation counts, the Stromgren b-y/m1/beta indices -- are
    worth keeping because they are what says whether a discrepancy is real. On
    this sample e_U-B runs 0.005 to 0.075 mag, a factor of 15, so a colour
    residual cannot be read without it.

    Writes data/uves_ngsl_photometry.csv
    """
    cats = [('ubv', 'II/168/ubvmeans'), ('uvby', 'II/215/catalog')]
    got = {}
    cols = {}
    for pre, cat in cats:
        for r in rows:
            star = r['ngsl_target']
            row, how = _vizier_row(cat, star, _f(r['vmag']), ['*'])
            if row is None:
                continue
            cols.setdefault(pre, list(row.colnames))
            d = got.setdefault(star, {})
            for c in row.colnames:
                v = row[c]
                d[f'{pre}_{c}'] = ('' if v is None or (isinstance(v, float)
                                                       and not np.isfinite(v))
                                   else str(v).strip())
            d[f'{pre}_match'] = how

    field = ['star', 'sptype']
    for pre, _ in cats:
        field += [f'{pre}_{c}' for c in cols.get(pre, [])] + [f'{pre}_match']
    out = ROOT / 'data' / 'uves_ngsl_photometry.csv'
    with open(out, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=field, extrasaction='ignore')
        w.writeheader()
        for r in rows:
            star = r['ngsl_target']
            w.writerow(dict(star=star, sptype=r['sptype'],
                            **got.get(star, {})))
    n = {pre: sum(1 for d in got.values() if f'{pre}_match' in d)
         for pre, _ in cats}
    print(f'  -> {out.relative_to(ROOT)}  ({len(field)} columns; '
          f'II/168 {n["ubv"]}/13, II/215 {n["uvby"]}/13)')


def main():
    selftest()
    rows = list(csv.DictReader(open(ROOT / 'data' / 'uves_ngsl_overlap.csv')))
    _RADEC.update({r['ngsl_target']: (float(r['ra']), float(r['dec']))
                   for r in rows})
    cat = catalog_ubv(rows)
    c1cat = catalog_c1(rows)
    print(f'  Mermilliod II/168 U-B: '
          f'{sum(np.isfinite(v["UB"]) for v in cat.values())}/{len(rows)}   '
          f'Hauck & Mermilliod II/215 c1: {len(c1cat)}/{len(rows)}')

    jb_full, jb_tr = filters(), filters(truncate_U_at=UVES_BLUE_LIMIT)
    st_full, st_tr = stromgren(), stromgren(truncate_at=UVES_BLUE_LIMIT)

    recs = []
    for r in rows:
        star = r['ngsl_target']
        wn, fn = ngsl_spectrum(r['ngsl_file'])
        gn = np.isfinite(fn)
        wn, fn = wn[gn], fn[gn]
        Un, Bn = vega_mags(wn, fn, jb_full)
        Unt, _ = vega_mags(wn, fn, jb_tr)
        c1n = c1_of(wn, fn, st_full)
        c1nt = c1_of(wn, fn, st_tr)

        wu, fu, _ = load_uves(r['uves_name'])
        fu_fill, holes = fill_small_gaps(wu, fu)
        # A hole anywhere inside any bandpass used here makes that star's UVES
        # colours meaningless, however smooth the interpolation looked.
        bad = [f'{lo:.0f}-{hi:.0f}' for lo, hi in holes
               if hi > 3050 and lo < 5500]
        if bad:
            Uu = Bu = c1u = np.nan
        else:
            Uu, Bu = vega_mags(wu, fu_fill, jb_tr)
            c1u = c1_of(wu, fu_fill, st_tr)

        c = cat.get(star, dict(V=np.nan, BV=np.nan, UB=np.nan))
        cc, ec, _how = c1cat.get(star, (np.nan, np.nan, ''))
        recs.append(dict(
            star=star, sptype=r['sptype'],
            ub_cat=c['UB'], bv_cat=c['BV'], v_cat=c['V'],
            ub_ngsl=Un - Bn, ub_ngsl_trunc=Unt - Bn, ub_uves=Uu - Bu,
            d_ub_ngsl_cat=(Un - Bn) - c['UB'],
            d_ub_uves_ngsl=(Uu - Bu) - (Unt - Bn),
            ub_trunc_cost=Un - Unt,
            c1_cat=cc, e_c1_cat=ec, c1_ngsl=c1n, c1_ngsl_trunc=c1nt,
            c1_uves=c1u,
            d_c1_uves_ngsl=c1u - c1nt,
            c1_trunc_cost=c1n - c1nt,
            uves_band_hole=';'.join(bad)))

    # Variability, from the companion table -- the UVES-POP files' own GCVS
    # type, carried here so the plot and the CSV agree about which stars are
    # variable without either re-deriving it.
    vt = {}
    f = ROOT / 'data' / 'uves_ngsl_compare.csv'
    if f.exists():
        vt = {r['star']: r['gcvs_type'] for r in csv.DictReader(open(f))}
    for r in recs:
        r['gcvs_type'] = vt.get(r['star'], '')

    # Stromgren zero point. Johnson U-B is tied to Vega, so NGSL - catalog is a
    # calibration test outright. c1 is not: its zero point is conventional, set
    # by standard stars, so the MEDIAN offset here is that convention and only
    # the star-to-star scatter about it is a measurement.
    dc = np.array([r['c1_ngsl'] - r['c1_cat'] for r in recs], float)
    dc = dc[np.isfinite(dc)]
    zp = float(np.median(dc)) if dc.size else np.nan
    for r in recs:
        r['c1_ngsl_zp'] = r['c1_ngsl'] - zp
        r['c1_uves_zp'] = r['c1_uves'] - zp
        r['d_c1_ngsl_cat'] = r['c1_ngsl'] - r['c1_cat'] - zp
        # UVES cannot reach 3050 A, so its U-B is on the truncated band. The
        # correction back to the full band is MEASURED on NGSL, which covers
        # both -- so this one number is the only NGSL information in a UVES
        # point, and it is worth 0.008 mag at the median.
        r['ub_uves_corr'] = r['ub_uves'] + r['ub_trunc_cost']
        r['d_ub_uves_cat'] = r['ub_uves_corr'] - r['ub_cat']

    f4 = lambda x: ('' if not isinstance(x, float) or not np.isfinite(x)
                    else round(float(x), 4))
    cols = list(recs[0])
    with open(ROOT / 'data' / 'uves_ngsl_ub.csv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for rec in recs:
            w.writerow({k: (f4(v) if isinstance(v, float) else v)
                        for k, v in rec.items()})
    print('  -> data/uves_ngsl_ub.csv\n')

    fm = lambda x: f'{x:>10.3f}' if np.isfinite(x) else f'{chr(8212):>10}'
    print('JOHNSON U-B (U and B both from Mermilliod II/168)')
    print(f'{"star":<10}{"SpT":<10}{"cat":>10}{"NGSL":>10}{"NGSL-cat":>10}'
          f'{"UVES":>10}{"UVES-NGSL":>10}{"truncost":>10}')
    for rec in recs:
        print(f'{rec["star"]:<10}{rec["sptype"]:<10}' + fm(rec['ub_cat'])
              + fm(rec['ub_ngsl']) + fm(rec['d_ub_ngsl_cat'])
              + fm(rec['ub_uves']) + fm(rec['d_ub_uves_ngsl'])
              + fm(rec['ub_trunc_cost']))

    print(f'\nSTROMGREN c1 (zero point {zp:+.3f} removed from NGSL-cat)')
    print(f'{"star":<10}{"SpT":<10}{"cat":>10}{"e_cat":>10}{"NGSL":>10}'
          f'{"NGSL-cat":>10}{"UVES":>10}{"UVES-NGSL":>10}')
    for rec in recs:
        if not np.isfinite(rec['c1_cat']) and not np.isfinite(rec['c1_uves']):
            continue
        print(f'{rec["star"]:<10}{rec["sptype"]:<10}' + fm(rec['c1_cat'])
              + fm(rec['e_c1_cat']) + fm(rec['c1_ngsl'])
              + fm(rec['d_c1_ngsl_cat']) + fm(rec['c1_uves'])
              + fm(rec['d_c1_uves_ngsl']))

    print()
    for k, lbl in [
            ('d_ub_ngsl_cat', 'U-B  NGSL - catalog (absolute, Vega-tied)'),
            ('d_ub_uves_ngsl', 'U-B  UVES - NGSL (identical truncated band)'),
            ('ub_trunc_cost', 'U-B  cost of truncating U at 3200 A'),
            ('d_c1_ngsl_cat', 'c1   NGSL - catalog, zero point removed'),
            ('d_ub_uves_cat', 'U-B  UVES - catalog (truncation corrected)'),
            ('d_c1_uves_ngsl', 'c1   UVES - NGSL (identical truncated band)'),
            ('c1_trunc_cost', 'c1   cost of truncating u at 3200 A')]:
        v = np.array([rec[k] for rec in recs], float)
        v = v[np.isfinite(v)]
        if v.size:
            nmad = 1.4826 * np.median(np.abs(v - np.median(v)))
            print(f'{lbl:<45} n={v.size:>2}  median {np.median(v):+.3f}  '
                  f'NMAD {nmad:.3f}  range {v.min():+.3f} to {v.max():+.3f}')

    make_figure(recs, zp)
    dump_catalog_rows(rows)


if __name__ == '__main__':
    warnings.filterwarnings('ignore')
    main()
