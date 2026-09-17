# Superseded LSF scripts

Everything in here was replaced by **`explore/ngsl_lsf.py`**, which does the
whole NGSL line-spread-function measurement in one place. See
[docs/LSF.md](../../docs/LSF.md) for the result.

These are kept because the project's record of *how a conclusion was reached*
is part of the work, and because two of them were wrong in ways worth being
able to point at. They are not maintained, they are not run by the pipeline,
and their outputs now go to `data/superseded/` so a stale CSV cannot be picked
up by mistake. Nothing imports them.

| file | what it did | why it was replaced |
|---|---|---|
| `ngsl_lsf_from_xsl.py` | Gaussian width vs XSL, 3 stars, 5 windows → `ngsl_lsf_measured.csv` | Gaussian only, so it measured a compromise between a core and a halo and reported it as a resolution (R = 600). Correct as a Gaussian fit; not a line spread function. |
| `ngsl_lsf_shape.py` | profile family comparison on one 400 A window → `ngsl_lsf_shape.csv` | Fitted a single control window (4200–4600 A) and could not produce the per-grating widths that were nonetheless quoted from it. Its `best` profile was selected on leftover Balmer core excess — a criterion this project has since retracted. |
| `ngsl_core_excess.py` | Balmer core excess, NGSL vs model and NGSL vs degraded XSL → `ngsl_core_excess.csv` | Its model-free half is now a standard output of `ngsl_lsf.py` (`core_mean`, `core_peak`, per-line in `ngsl_lsf_lines.csv`) across 9 stars rather than 8, and for every profile rather than two. Its model-dependent half produced the retracted "86%" claim. |

## What was wrong, specifically

Worth keeping in view, because both errors were invisible in the outputs.

**A selection criterion the project had already retracted was still in the
code.** `ngsl_lsf_shape.py` ends by choosing the profile with the smallest
leftover Balmer core excess. That excess is degenerate with the kernel's
effective *width* and says nothing about its *shape* — a plain Gaussian spans
+16.8% to −1.2% across plausible widths — which `docs/CAVEATS.md` stated while
the script went on ranking profiles by it.

**The adopted numbers were not reproducible from the committed code.**
`common/lsf.py` carried `NGSL_MOFFAT_LSF = [..., 4.02, 8.34]` and
`NGSL_MOFFAT_BETA = 1.6`. Recomputing from the committed `ngsl_lsf_shape.csv`
gives a median Moffat core of **4.096 A** and **beta = 1.556**, and no script in
the repository fitted G750L at all, so **8.34 A had no producing code**. The
three-window average that 4.02 came from was assembled by hand.

That is the reason `ngsl_lsf.py` writes every fitted parameter, for every star,
grating, profile and sampling convention, to `data/ngsl_lsf.csv` — the adopted
constants are a documented reduction of that table, and re-running the script
regenerates them.
