# Nickel 1-m Automation — Project Overview

## Purpose

Lick Observatory's Nickel 1-m Telescope requires observers to manually
focus the telescope and verify pointing at the start of each night — a
repetitive, ~10-minute-per-iteration procedure performed several times per
session. This repository automates that startup sequence using the Keck
Task Library (KTL) Python API, which is the control interface for the
Nickel's mechanisms (secondary focus, pointing) and its science camera.

The project was started by undergraduate student Scott Hakoda (Utah Tech
University) in summer 2025, mentored by Kyle Westfall and Will Deich
(UCO Santa Cruz). Most of the repository is unchanged since then. The one
file under active revision is `scripts/focus.py`, which Kyle Westfall began
refactoring from the original `scripts/automate.py` approach.

> **Note on naming:** the production code directory was originally named
> `ktl/`, which turned out to collide with Python's import resolution for
> the real `ktl` package (a directory literally named `ktl` sitting on
> `sys.path` can shadow the real module). It was renamed to `scripts/` to
> avoid that, and may be renamed again — check the current top-level
> directory listing rather than assuming `scripts/` is permanent.

## Repository state at a glance

The repo mixes three tiers of code that are useful to distinguish when
reading it:

| Tier | Location | Status |
|---|---|---|
| Production / control code | `scripts/` | Talks to the real telescope via `ktl.cache(...)`. `automate.py` is last year's working version; `focus.py` is the in-progress refactor superseding it. |
| Prototype | `practice/` | Earlier developmental versions of the same ideas (focus sequencing, photometry, quadratic fitting), built and tested against **simulated** FITS images rather than real hardware. |
| Learning scratch work | `practice-ktl/` | Small standalone scripts and notes used to learn the `ktl` API's `cache`/`callback`/`monitor`/`waitFor` mechanics. Not meant to run as part of the pipeline. |

`obs_images/` is not code — it's real logs and notes from actual focus runs
at the telescope (`focusing.notes`, plus captured KTL/terminal output from
live sessions), useful as ground truth for what real keyword behavior looks
like.

## `scripts/` — production code

This is the code intended to actually run at the telescope (each script
starts with the `#! @KPYTHON@` shebang used by the Keck-style Python
environment).

- **`automate.py`** (631 lines) — The original end-to-end focus-finding
  driver. Defines a `Keyword` class that opens and monitors all the KTL
  keywords needed (`POCSECPA/D/LK`, `nickucam` exposure keywords, `EVENT`),
  an `Event` class that performs one "change focus → expose → measure"
  step, and both an automatic curve-following search (`curve_finder`,
  `curve_helper`) and a manual fixed-range search
  (`manual_focus_finder`). Fits a quadratic to FWHM-vs-focus via
  `quadratic.py` and plots results live via a `Grid` helper. Also supports
  `--refit` (recompute the curve from a saved `focus_data.ecsv`, omitting
  bad observations) and `--reevaluate` (re-run photometry on an old focus
  sequence's images).
- **`focus.py`** (564 lines) — Kyle's refactor of the above, restructured
  around small single-responsibility KTL wrapper classes (`Focus`,
  `ExposurePath`, `ExposureConfig`, `Exposure`) and a `FocusSequence` base
  class with three concrete strategies: `GridFocusSequence` (fixed step
  grid), `AutomatedFocusSequence` (adaptive curve-following, ported from
  `automate.py`'s logic but as stateful step methods), and
  `ArchiveFocusSequence` (replay previously-taken exposures for
  reprocessing instead of re-observing, driven by `--obsnum`). This is the
  file actively being worked on; several pieces are explicitly unfinished
  (see below). It can now run without a live `ktl` connection when
  reprocessing an archived sequence: `FocusSequence` only constructs its
  hardware-control objects (`Focus`, `Exposure`) if `ktl` imported
  successfully, and the `--obsnum` archive path builds exposure filenames
  from explicit `--datadir`/`--prefix`/`--suffix` arguments instead of
  reading the recording directory/prefix/suffix from KTL keywords. This
  path has been confirmed working end-to-end against a real archived
  exposure sequence, without `ktl` installed.
- **`photometry.py`** (513 lines) — Image-quality measurement used by
  `focus.py`. Iteratively detects sources and estimates background
  (`find_sources`, using `photutils.segmentation` + sigma-clipping),
  characterizes source shape via image moments to get FWHM
  (`evaluate_shape`, `evaluate_sources`, `moment2d`), and picks a target
  source either by brightness, flux-weighted mean, or proximity to
  user-given coordinates (`image_quality`). This is a cleaner rewrite of
  the photometry approach prototyped in `practice/`.
- **`move_to_target.py`** — Given the telescope's current RA/Dec (read
  from KTL) and a catalog of known pointing/focus stars
  (`point_focus.txt`), finds the nearest star via `SkyCoord.separation`
  and commands the telescope to slew to it (`POCRAD`/`POCDECD`).
- **`quadratic.py`** — Tiny shared helper: fit a quadratic to
  (focus, FWHM) points and return its vertex (the best-focus estimate).
- **`test_photometry.py`** — Ad hoc manual test script (not pytest-style —
  it calls its test function at import time) that runs `photometry.py`
  against a real FITS frame from a hardcoded local path.
- **`point_focus.txt`**, **`focus_data.ecsv`** — Data files: a static
  catalog of standard pointing/focus stars, and an example saved focus
  sequence (ECSV table of ObsNum/Focus/FWHM/Centroid) of the kind
  `focus.py`'s `ArchiveFocusSequence` / `automate.py`'s `--refit` consume.

## `practice/` — prototype work

Earlier iteration of the same pipeline, built and exercised against
synthetic data so it didn't require telescope time:

- `simulation.py` / `ktl_simulation.py` generate synthetic Moffat-PSF FITS
  images with a specified FWHM.
- `practice/focus.py`'s `Image` class calls the simulator directly in
  place of a real camera exposure.
- `practice/automate.py` (242 lines) is a thinner predecessor of
  `scripts/automate.py`, using the simulated `Image` instead of real
  photometry/KTL calls.
- `practice/photometry.py`, `DAO_find.py`, `fwhm_practice.py` are rougher,
  more exploratory drafts of source-detection/FWHM logic (including
  `photutils.DAOStarFinder`-based attempts) that predate the moment-based
  approach finalized in `scripts/photometry.py`.
- `practice/quadratic.py` is nearly identical to `scripts/quadratic.py`.
- `plot_practice.py` is scratch matplotlib/gridspec layout experimentation.
- `data/nickel/`, `images/`, `images/test_images/` hold sample FITS
  files (mostly synthetic, plus a few real test frames) used only as
  inputs while developing the code above.

## `practice-ktl/` — KTL API learning exercises

- `ktl_behavior.py`, `pie_behavior.py` — small scripts probing how
  `ktl.cache().callback()/.monitor()/.waitFor()` behave, one against a
  toy `pie` service and one against real Nickel keyword names, used to
  work out polling vs. event-driven patterns.
- `ktl.txt` — scratch pseudocode outlining the focus loop; effectively an
  early design sketch of what became `scripts/focus.py` / `automate.py`.
- `ktl.focus.keywords` — a dump of real KTL keyword definitions
  (`POCSECPD`, `EVENT`, `RECORD`, etc.) kept as an API reference.

## Dependencies

`requirements.txt` currently lists `astropy`, `IPython`, `matplotlib`,
`numpy`, `photutils`. The `ktl` package itself (the observatory's telescope
control middleware) is not pip-installable and is assumed to be provided
by the Lick/Keck software environment. Note `photometry.py` also imports
`scipy.ndimage`, which is not currently listed in `requirements.txt`.

`photometry.py` uses the `photutils.segmentation` API names introduced in a
recent `photutils` release (`n_sigma`, `n_pixels`, `n_labels`,
`data_masked` — replacing the older `nsigma`/`npixels`/`nlabels`/`data_ma`).
This effectively requires `photutils` 2.x; `requirements.txt` does not
currently pin a minimum version, so installing an older `photutils` will
fail with an unexpected-keyword or missing-attribute error rather than an
obviously version-related one.

## Known rough edges / open items

These are things worth knowing before extending the code — not all are
bugs, but they mark where the codebase is unfinished or where the two
"production" scripts have diverged:

- **`automate.py` and `photometry.py` are out of sync.** `automate.py`
  does `from photometry import photometry, Grid`, but the current
  `scripts/photometry.py` defines `image_quality` (and no `Grid` class) —
  `photometry.py` was rewritten for `focus.py`'s API and `automate.py`
  was never updated to match. `automate.py` should be treated as
  historical/reference rather than runnable as-is.
- **`focus.py --refit` is a stub** — it immediately raises
  `NotImplementedError`, and the code after it (reading a saved ECSV,
  honoring `--omit`) is unreachable dead code left in place for when it's
  implemented.
- **`focus.py`'s TODOs**: writing the output data file and plotting
  results at the end of `main()` are marked but not implemented; a
  `# TODO` also flags uncertainty about how to correctly check the
  `EXPREC` keyword state, and another about the exact keyword name for
  starting an exposure (`'StartX'`, pending confirmation with Will
  Deich).
- **Hardcoded assumptions**: `automate.py`'s `detect_outliers` hardcodes a
  1024-pixel image size for cutout bounds (flagged in-line as
  "USING HARDCODED 1024 FOR NOW").
- **`test_photometry.py`** references an absolute path
  (`/Users/westfall/Work/Lick/Nickel/science_camera/obs/2025-08-14/...`)
  specific to one machine, and runs as a script rather than via
  `pytest`.
- **The production directory has been renamed once already** (`ktl/` →
  `scripts/`, to avoid the import collision noted above) and may move
  again; anything that hardcodes the path (docs, shell aliases, cron jobs)
  should be easy to update.

## Suggested next steps

- Decide whether `scripts/automate.py` should be retired/archived now that
  `focus.py` covers its functionality with a cleaner design, to avoid the
  stale-import trap above.
- Finish the TODOs in `focus.py` (`--refit`, output-file writing, results
  plotting) so it reaches parity with what `automate.py` could already do.
- Add `scipy` to `requirements.txt`, and pin `photutils>=2.0` given the
  API-name dependency noted above.
- Consider whether `practice/` and `practice-ktl/` should be archived out
  of the main tree (e.g. moved under a clearly-labeled `archive/` or
  removed from version control) now that `scripts/` has superseded them,
  to reduce confusion for new contributors.
