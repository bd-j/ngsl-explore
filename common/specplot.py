"""Shared spectrum-panel drawing, so the figures cannot drift apart.

`explore/check_predict.py` and `explore/plot_metal_lines.py` both draw an XSL
window with the observation, one or more models and the masked regions shaded.
They had grown separate copies of that code, which is the same drift that once
let the fitter and the comparison figures disagree about the NGSL LSF -- and it
matters here for a specific reason: the convention that a MASKED pixel is drawn
faded and a FITTED one solid is what stops an exclusion being misread as missing
data. If one figure keeps that convention and the other quietly loses it, the
two say different things about the same spectrum.
"""
import numpy as np

OBS_C, MOD_C, MOD2_C = '#2a78d6', '#eb6834', '#7a3fa8'
BAND_C, HELD_C = '#2a78d6', '#c0392b'
SURFACE, INK, MUTED, GRIDC = '#fcfcfb', '#22262b', '#6b7280', '#dfe3e8'
MODEL_COLORS = (MOD_C, MOD2_C)


def gapped(y, keep):
    """`y` with everything outside `keep` set to NaN, so a line BREAKS there.

    Selecting the kept pixels instead -- ax.plot(w[keep], y[keep]) -- joins the
    last point before a masked stretch to the first point after it, and
    matplotlib draws that join as a straight segment right across the gap. On a
    Balmer panel that segment runs flat across the masked core at continuum
    level, which is the most misleading thing the figure could draw there: it
    looks like an observation that has no line in it. NaN breaks the path and
    leaves the gap empty, where the faded full-range pass shows the real
    profile.
    """
    out = np.asarray(y, float).copy()
    out[~np.asarray(keep, bool)] = np.nan
    return out


def style(ax):
    """House axis styling."""
    ax.set_facecolor(SURFACE)
    ax.grid(alpha=.25, color=GRIDC, lw=.7)
    ax.tick_params(labelsize=8, colors=MUTED)
    for s in ax.spines.values():
        s.set_color(GRIDC)


def spectrum_panel(ax, obs, models, lo, hi, title=None, band=None,
                   obs_label='XSL (solid = fitted)', shade_masked=True,
                   legend=False, fontsize=9):
    """Draw one wavelength window: observation, models, and what is excluded.

    obs      an Observation with .wavelength, .flux, .mask
    models   [(label, calibrated_flux_array)], same length as obs.flux
    band     optional (lo, hi) to tint, e.g. the feature the panel is about

    The observation is drawn twice: faded across the whole range and solid only
    where `obs.mask` is True. Masked stretches are tinted, so a deliberate
    exclusion looks deliberate rather than like a gap in the data.
    """
    style(ax)
    w = np.asarray(obs.wavelength, float)
    inr = (w > lo) & (w < hi)
    if inr.sum() < 3:
        ax.set_visible(False)
        return False
    fit = np.asarray(obs.mask, bool) & inr

    if band is not None:
        ax.axvspan(band[0], band[1], color=MUTED, alpha=.06, lw=0)
    if shade_masked:
        gap = inr & ~np.asarray(obs.mask, bool)
        if gap.any():
            # tint each contiguous masked run, not one span across all of them
            idx = np.where(gap)[0]
            for run in np.split(idx, np.where(np.diff(idx) > 1)[0] + 1):
                if run.size:
                    ax.axvspan(w[run[0]], w[run[-1]], color=HELD_C,
                               alpha=.09, lw=0)

    flux = np.asarray(obs.flux, float)
    ax.plot(w[inr], flux[inr], color=OBS_C, lw=.8, alpha=.30)
    ax.plot(w[inr], gapped(flux, fit)[inr], color=OBS_C, lw=1.15,
            label=obs_label)
    for (lab, m), c in zip(models, MODEL_COLORS):
        # same convention as the observation: faded across the whole window,
        # solid only where it was actually fitted. Without the faded pass the
        # model simply vanished inside a masked Balmer core, so the panel showed
        # a hole in the middle of the line it exists to display.
        mv = np.asarray(m, float)
        ax.plot(w[inr], mv[inr], color=c, lw=.8, alpha=.30)
        ax.plot(w[inr], gapped(mv, fit)[inr], color=c, lw=1.0, label=lab)

    ax.set_xlim(lo, hi)
    if title:
        ax.set_title(title, fontsize=fontsize, color=INK)
    ax.set_xlabel(r'$\lambda$ [$\AA$, vacuum]', fontsize=8, color=INK)
    ax.set_ylabel(r'F$_\lambda$', fontsize=8, color=INK)
    if legend:
        ax.legend(fontsize=7, loc='lower left', framealpha=.92)
    return True


def offscale_note(ax, values_pct, ylim, color=HELD_C):
    """Say so when a fixed residual range hides points off the top.

    The break and Balmer residual panels all use FIXED ranges, so that a 0.5%
    residual and a 15% one cannot look identical between stars. The cost is
    that a panel can imply the residual stayed inside a range it left, so
    every one of them carries this note instead.
    """
    v = np.asarray(values_pct, float)
    v = v[np.isfinite(v)]
    n_out = int(np.sum((v < ylim[0]) | (v > ylim[1])))
    if not n_out:
        return 0
    ax.text(0.99, 0.04,
            f'{n_out} px outside ±{ylim[1]:.0f}% (to {np.max(np.abs(v)):.0f}%)',
            transform=ax.transAxes, ha='right', va='bottom', fontsize=6.5,
            color=color, bbox=dict(fc='white', ec='none', alpha=.8, pad=1.5))
    return n_out


def residual_panel(ax, w, models, lo, hi, bands=(), ylim=None, marks=(),
                   ylabel='(data − model)/model [%]', xlabel=True,
                   note=True):
    """Fractional residuals over one window, fitted and predicted drawn apart.

    models  [(resid, used, pred[, median_masks])] per model curve, in
            MODEL_COLORS order. `resid` is a FRACTION and is plotted as a
            percentage; `used` are the pixels that set the model's scale and
            `pred` the ones that did not. They are drawn differently for the
            same reason the spectrum panels fade a masked pixel: a residual on
            a pixel that helped fit the scale is not a prediction, and a panel
            that draws the two alike invites reading one as the other.

            `median_masks` is a LIST of masks, each getting its own dotted
            median drawn across its own wavelength span. It defaults to
            [pred], which is right for a panel holding one feature -- but a
            Balmer panel quotes the median over the masked CORE while `pred`
            may also hold a bad region 40 A away, and a panel whose dotted
            line is not the number in its title is worse than no line.
    bands   [(lo, hi)] to tint -- the masked cores, for a Balmer panel. A list
            rather than one range because a panel may hold two members.
    marks   [(lambda, label)] vertical markers, e.g. an interstellar line

    """
    style(ax)
    w = np.asarray(w, float)
    for blo, bhi in bands:
        ax.axvspan(blo, bhi, color=HELD_C, alpha=.09, lw=0, zorder=0)
    for lam, label in marks:
        if lo < lam < hi:
            ax.axvline(lam, color=MUTED, ls=':', lw=.9, alpha=.8)
            ax.text(lam, 0.96, label, transform=ax.get_xaxis_transform(),
                    rotation=90, va='top', ha='center', fontsize=6, color=MUTED,
                    bbox=dict(fc=SURFACE, ec='none', alpha=.75, pad=.6))
    ax.axhline(0, color=MUTED, lw=1)

    shown = []
    for entry, c in zip(models, MODEL_COLORS):
        resid, used, pred = entry[:3]
        r = 100 * np.asarray(resid, float)
        inr = (w > lo) & (w < hi)
        # gapped, not selected: the scale-setting wings sit on BOTH sides of
        # the core, so selecting them would draw a segment straight across it
        ax.plot(w[inr], gapped(r, used)[inr], color=c, lw=.7, alpha=.35)
        ax.plot(w[inr], gapped(r, pred)[inr], color=c, lw=1.2)
        for mm in (entry[3] if len(entry) > 3 else [pred]):
            sel = inr & np.asarray(mm, bool) & np.isfinite(r)
            if sel.any():
                ax.plot([w[sel].min(), w[sel].max()],
                        [float(np.median(r[sel]))] * 2, color=c, ls=':',
                        lw=1.4, zorder=4)
        shown.append(r[inr])

    ax.set_xlim(lo, hi)
    if ylim is not None:
        ax.set_ylim(*ylim)
        if note and shown:
            offscale_note(ax, np.concatenate(shown), ylim)
    ax.set_ylabel(ylabel, fontsize=8, color=INK)
    # `xlabel` drops the axis TEXT only, never the tick labels: every member in
    # explore/plot_balmer_lines.py sits in a different window, so a panel whose
    # wavelengths cannot be read is not worth drawing. Only the repeated
    # 'lambda [A, vacuum]' under ten stacked panels is worth dropping.
    if xlabel:
        ax.set_xlabel(r'$\lambda$ [$\AA$, vacuum]', fontsize=8, color=INK)
