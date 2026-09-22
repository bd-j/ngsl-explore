# Stale results, dead ends and abandoned work

Everything here is **wrong, superseded or abandoned**. Nothing in this file
describes how the project works now. It exists for one reason: these numbers and
designs are still findable — in old figures, in commit messages, in notebooks,
in `explore/superseded/` — and a reader who meets one needs to know it is dead
and what replaced it.

The current state lives in [DATA.md](DATA.md), [GRID.md](GRID.md),
[FITTING.md](FITTING.md) and [LSF.md](LSF.md). Those four describe things as
they are and do not re-litigate how they got there.

For traps that are still live — ways the comparison goes wrong *today* — see
[CAVEATS.md](CAVEATS.md). A trap is not a dead end.

---

## Retracted numbers

| number | where it came from | why it is wrong | what replaced it |
|---|---|---|---|
| **R = 600 ± 40**, constant in velocity | a single Gaussian fitted to NGSL against XSL (`superseded/ngsl_lsf_from_xsl.py`) | a one-parameter profile forced to stand in for a core plus a halo. It inflates to split the difference, and drifts with wavelength as the halo's weight changes — which is what made the profile look constant in *velocity* when its core is constant in *Ångströms* | Moffat, core 3.54 Å (G430L) / 8.38 Å (G750L), β = 1.52 — [LSF.md](LSF.md) |
| **R = 939** | the tabulated STIS LSF, read as a resolving power | the tabulated profile is not Gaussian; a FWHM read off it is not an R | as above |
| **R = 665** | pixel sampling | same error, different route | as above |
| **R = 950** | printed by `common/uves_pop_load.py` as an aside | inherited from R = 939 | the line now names the Moffat core instead of quoting an R |
| **"a winged profile removes 86% of the core excess"** | `superseded/ngsl_core_excess.py` | the core excess is degenerate with the effective **width**, not diagnostic of the **shape**. Changing only the Moffat core at HD194453's ML node runs it from +2.68% at 3.54 Å to −4.09% at 7.00 Å — any profile can be tuned to zero it | measured against XSL rather than tuned: median **+0.42%** over the eight in-grid stars, both signs, consistent with zero — and only quotable alongside the profile in use ([LSF.md](LSF.md)) |
| `NGSL_R_MEASURED = 600` | `common/lsf.py` | as above. **Removed from the code**, not kept for reference — it is neither a width nor a resolution | `common.lsf.to_ngsl_pixels` |
| core 4.02 Å (G430L), 8.34 Å (G750L), β = 1.6 | `common/lsf.py`, adopted constants | **not reproducible from the committed code.** The committed `ngsl_lsf_shape.csv` gives a median core of 4.096 Å and β = 1.556, and no script in the repository fitted G750L at all | `ngsl_lsf.py` now writes every fitted parameter for every star, grating, profile and sampling convention to a CSV, and the adopted numbers are a stated reduction of it |

---

## Retracted and revised claims

### The Balmer cores are NLTE in hydrogen
Revised twice, then settled. The entry first read that the filled Balmer cores
were "most likely NLTE in hydrogen". That did not survive: at XSL's R ~ 9800 the
same models fit the **full Hγ profile, core included**, for the metal-rich
stars, and a physical NLTE core deficit cannot be present at NGSL's resolution
and absent at R = 9800.

The replacement claim — that a winged profile "removes 86% of the core excess" —
was also wrong, for the reason in the table above.

A third revision concerned how far the conclusion reached. "The core excess is
the instrument profile, not NLTE" was reached with a width that *nulled* the
excess, so the number attached to it was never meaningful.

Current statement in [CAVEATS.md](CAVEATS.md) and [LSF.md](LSF.md).

### The experiment is dust-limited
An early draft of FITTING.md called the problem dust-limited. It is
**Teff-limited**: dD/dTeff = −0.017 mag per 100 K against dD/dE(B−V) = +0.0027
per 0.01 mag at fixed Teff. The correction mattered because it changes which leg
of the fit deserves the work — the XSL line-profile leg is the critical path,
not the dust lever.

### Gaia XP validates the NGSL calibration
An earlier version of the XP section read the median flux ratio of **1.011**
(XP against NGSL for HD194453, 4000–9000 Å) as validating both calibrations.
That was too strong a conclusion from a median: it validates the mean *level*
and the unit conversion, and says nothing about the *colour* — which is the
whole dust lever. XP and NGSL disagree in colour by up to 4.8%.

### The NLTE question is re-opened
An earlier version of the core-excess section generalised HD194453's value to
the sample and concluded the NLTE question was re-opened. Eight stars say it is
not. This is the specific failure of publishing a one-star result as a sample
conclusion.

### UVES-POP's coverage gaps were attributed to XSL

DATA.md's XSL entry listed "a dichroic gap at 5750-5844 A, inter-order gaps
every ~150 A redward of 8515 A, and one star (HD162678) missing 3859-4779 A
outright" as **XSL's** coverage gaps, and said CAVEATS.md explained "why the
blue has no gaps and the red does" for XSL.

All three belong to UVES-POP. XSL has 486 points in 5750-5844 A and 610 in
8515-8690 A, and HD162678 is not in XSL at all — it is a UVES-POP star. The
echelle order-width argument in CAVEATS.md is about UVES-POP's spectrograph
and says nothing about X-shooter.

Replaced by a checked statement of XSL's actual coverage in
[DATA.md](DATA.md) and the UVES-POP gap survey in [CAVEATS.md](CAVEATS.md).

### The blue of a UVES-POP spectrum has 100% coverage in every star

CAVEATS.md's echelle order-width entry concluded that because orders overlap
in the blue, "the Balmer region has 100% coverage in every star while the red
end is riddled with holes."

True of ORDER gaps, and false as written. A missing spectral SETTING is a
different failure that owes nothing to order spacing, and three of the 22
UVES-POP spectra on disk have one in the blue: HD162678 (3859-4779 A),
HD138716 (3859-4784 A) and **Betelgeuse (3201-3753 A, which removes the
Balmer break outright)**. The same entry's own table already listed
HD162678's hole, so the document contradicted itself.

Replaced by the survey in [CAVEATS.md](CAVEATS.md): every spectrum has holes,
13-15 of them, 5.9-40% of pixels NaN.

### `HB*` is a peculiarity type
Horizontal-branch stars were originally rejected as chemically peculiar,
conflating two different problems. A field HB star at 9000–11000 K sits *below*
the ~11500 K Grundahl jump where radiative levitation starts, so a scaled-solar
atmosphere still describes it. What actually disqualifies these stars is low
log g and low [M/H] — a grid-coverage question, now recorded as such.
Reclassifying it is why HD143459, HD074721 and HD128801 are in the sample.

---

## Superseded designs

### The emcee whole-spectrum fitter (`fitting/fit.py`, deleted)
The earlier design fitted the whole NGSL spectrum with emcee over free Teff,
log g, [M/H], E(B−V), v sin i, `inst`, RV and an error scale, with the break
**included** in the fit and a `DustPrior` doing the work of separating Teff from
reddening.

Two things killed it:

* **A Gaussian dust prior does not bite.** With 1466 pixels and a free error
  scale the likelihood formally measures E(B−V) to ~0.001, so N(0.00, 0.02) came
  back 2σ out at 0.038.
* **Fitting the break and then reporting its residual answers a different
  question** from the one the project asks. Hence the conditioning / held-out
  design in [FITTING.md](FITTING.md).

The file was deleted rather than kept, because it held a second copy of the
forward model. That is the drift that once turned a request for R = 600 into
R = 83 — the fitter reimplemented its own kernels, `kind='R'` ignored the grid
spacing, and because `inst` was free nothing crashed. The kernels now delegate
to `common.lsf`.

### UVES-POP as a fitted library
UVES-POP is no longer fitted. It shares **no star** with NGSL, so it can only
ever be a separate sample rather than a cross-check on the same object, and its
continuum normalisation is too uncertain for a break measurement.
`explore/plot_uves_vs_model.py` still runs and its figures are in
`figures/explore_libraries/`, but nothing downstream depends on them.

### The earlier NGSL + UVES-POP sample
The sample below was used for the *catalog-parameter comparison* figures, before
the fit was set up. The current sample — NGSL ∩ XSL, 13 stars — is in
[DATA.md](DATA.md).

**NGSL**

| star | Teff | log g | [M/H] | E(B-V) | slit offset | notes |
|---|---|---|---|---|---|---|
| HD194453 | 10241 | 3.9 | +0.0 | -0.01 | **0.00 px** | primary target |
| HD040573 | 10200 | 4.2 | -0.4 | 0.06 | 0.36 px | classification and gravity agree |
| HD128801 | 10123 | 3.7 | -1.9 | 0.03 | 0.58 px | metal-poor comparison; in MILES |
| HD143459 | 9878 | 3.6 | -0.6 | 0.04 | 0.09 px | horizontal-branch star |
| ~~HD147550~~ | 10074 | 3.9 | -0.0 | **0.125** | 0.24 px | dropped: reddened |

HD194453 was the primary target because its **slit offset is 0.00 px** — the
wavelength-dependent slit-throughput correction, the dominant systematic on
break *shape*, is essentially null for it.

**UVES-POP**

| star | Teff | log g | [Fe/H] | v sin i | E(B-V) | S/N | notes |
|---|---|---|---|---|---|---|---|
| HD162678 | 9908 | 3.53 | +0.03 | 37 | 0.077 | 169 | slowest rotator; best for profiles |
| HD188294 | 11016 | 4.04 | +0.05 | 182 | 0.040 | 299 | best S/N and lowest reddening |
| HD162393 | 9955 | 4.05 | -0.55 | 142 | 0.065 | 132 | metal-poor comparison |
| ~~HD162817~~ | 10153 | 3.66 | +0.11 | 65 | **0.102** | 176 | dropped: reddened |
| ~~HD162630~~ | 10494 | 3.92 | -0.18 | 46 | 0.01 | 226 | dropped: spectroscopic binary |

UVES-POP spectra were **dereddened** with CCM89 (R_V = 3.1) using the library's
own fitted E(B-V), and models rotationally broadened with each star's published
`v sin i`. The two samples are **disjoint**, so they were independent tests
rather than a repeat measurement; between them they spanned 9878–11016 K,
log g 3.5–4.2, [M/H] −1.9 to +0.1 and rotation 37–182 km/s.

---

## Abandoned conventions

### Truncating the LSF kernel at a fixed number of FWHM
The kernel was cut at **40 × FWHM**. That gave the G750L kernel a 335 Å reach
against G430L's 142 Å — the same profile behaving differently in the two
gratings for no instrumental reason. It is now cut in **detector pixels**
(`NGSL_TRUNC_PX = 15`), the unit the instrument works in, which also keeps the
reach comparable to the tabulated STIS profiles that stop at ±20 px.

### A line-free sub-window for the LSF fit
A 5100–5647 Å window was tried. Four of nine stars walked their fit to the
guard, the survivors scattering by ±0.95 Å: the red end of G430L has almost no
features in an A star. Every sub-window is now anchored on a Balmer line.

### The blue continuum as a median at ~3570 Å
An early version of the break metric took the blue continuum as a median centred
at ~3570 Å rather than extrapolating it to 3646 Å, so a sloping SED did not
cancel. Values shifted by 0.02–0.09 mag when fixed. See
`common/balmer_metric.py`.

### Normalising in a window that contains Hδ
An early version scaled fluxes in a 4000–4200 Å window, which *contains* Hδ at
4102.9 Å. Each break panel is now normalised in windows bracketing that break,
with every Balmer and Paschen line masked ±20 Å.

### Writing the smoothed R = 10,000 CSV for every grid node
`make_model.py` writes `models/<name>_R10000.csv` for single-star use, and the
grid suppresses it (`build_grid.py` passes `--no-csv`). Over the full grid it
would have been **~19 GB and ~2.4 hours** of pure waste; 385 such files
(4.2 GB) were written before this was noticed and have been deleted. The
reasoning is in [GRID.md](GRID.md).

---

## Superseded scripts

Three are kept in [`explore/superseded/`](../explore/superseded/) with a README
recording what each got wrong:

* **`ngsl_lsf_from_xsl.py`** fitted Gaussians and reported R = 600. Correct as a
  Gaussian fit; not a line spread function.
* **`ngsl_lsf_shape.py`** compared profile families on a single 400 Å window and
  selected the winner on leftover Balmer core excess — a criterion this project
  had already retracted as degenerate with width.
* **`ngsl_core_excess.py`** produced the retracted "a winged profile removes 86%
  of the core excess" claim.

All three are superseded by `explore/ngsl_lsf.py`, which does the whole
measurement in one script and writes every fitted parameter to CSV.
