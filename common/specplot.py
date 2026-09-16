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

    ax.plot(w[inr], np.asarray(obs.flux, float)[inr], color=OBS_C, lw=.8,
            alpha=.30)
    ax.plot(w[fit], np.asarray(obs.flux, float)[fit], color=OBS_C, lw=1.15,
            label=obs_label)
    for (lab, m), c in zip(models, MODEL_COLORS):
        # same convention as the observation: faded across the whole window,
        # solid only where it was actually fitted. Without the faded pass the
        # model simply vanished inside a masked Balmer core, so the panel showed
        # a hole in the middle of the line it exists to display.
        mv = np.asarray(m, float)
        ax.plot(w[inr], mv[inr], color=c, lw=.8, alpha=.30)
        ax.plot(w[fit], mv[fit], color=c, lw=1.0, label=lab)

    ax.set_xlim(lo, hi)
    if title:
        ax.set_title(title, fontsize=fontsize, color=INK)
    ax.set_xlabel(r'$\lambda$ [$\AA$, vacuum]', fontsize=8, color=INK)
    ax.set_ylabel(r'F$_\lambda$', fontsize=8, color=INK)
    if legend:
        ax.legend(fontsize=7, loc='lower left', framealpha=.92)
    return True
