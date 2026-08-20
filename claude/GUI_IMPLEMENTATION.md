# Nickel Focus/Pointing GUI — Implementation Log

This document tracks progress implementing the plan in `GUI_DESIGN.md`:
what's been done in each phase, any issues hit along the way, and any
places the implementation ended up diverging from the original design.
Once GUI development is complete, `GUI_DESIGN.md` should be revised to
match what was actually built, using this log as the record of what
changed and why.

Each phase gets its own section below, added as it's completed.

## Phase 1: Model refactor

**Design doc reference:** §9 Phased plan, item 1 ("Model refactor"); test
approach per §8.

### Summary

Refactored `scripts/focus.py` so the sequence-stepping engine is
decoupled from CLI-specific plotting, as required before any Controller
code (CLI or future GUI) can drive it:

- Added a `StepResult` dataclass — `index`, `focus_value`, `exposure`,
  `fwhm`, `frame`, `stamp`, `centroid`, `is_outlier` — as the shared unit
  of output for one sequence step.
- Split `FocusSequence.execute()`'s `while` loop into a `step()`
  generator (hardware + photometry only, no plotting) and a thin
  `execute()` that iterates `step()` and feeds each `StepResult` to
  `FocusPlot`. `step()` is the same seam a future GUI worker thread will
  drive (§4.3 of the design doc).
- While touching this code, also tightened `Path`/`str` consistency:
  `ExposurePath.previous` and `ExposurePath.for_obsnum` now return
  `pathlib.Path` instead of `str`, `StepResult.exposure` is typed `Path`,
  and `main()`'s archive-mode file list (`--obsnum`/`--datadir`/
  `--prefix`/`--suffix`) is built as `Path` objects throughout. This also
  fixed a latent bug: the "missing files" `FileNotFoundError` message
  used to call `str.join()` on what would have been a list of `Path`
  objects.
- New docstrings use Sphinx cross-references (`:class:`numpy.ndarray``,
  `:func:`FocusPlot.is_outlier``, etc.) per the project's documentation
  convention.
- Added a top-level `tests/` directory (§8): `conftest.py` (adds
  `scripts/` to `sys.path`, forces the `Agg` matplotlib backend, and a
  `focus_sweep` fixture generating small synthetic Gaussian-source FITS
  frames on a known FWHM-vs-focus quadratic), `test_quadratic.py`,
  `test_focus_sequence.py`, and `test_focus_cli.py`.
- `requirements.txt` updated to add `scipy` and `pytest` (both were
  already-used/needed but missing from the file); `pytest` installed into
  the `nickel` virtualenv.

### Deviations from `GUI_DESIGN.md`

- §4.1's sketch of `StepResult` included an `obsnum` field. The actual
  codebase doesn't track a literal integer observation number separate
  from the exposure's filename — `ExposurePath` reads `OBSNUM` from KTL
  only while live-exposing, and `ArchiveFocusSequence` never retains it
  as its own field. `StepResult.exposure` (a `Path`) is used instead;
  consumers that want a display label derive one from
  `exposure.stem`, same as `FocusPlot` already did. Revisit this naming
  when `GUI_DESIGN.md` gets its post-implementation revision.

### Issues encountered

- The synthetic-data test for `execute()`'s fitted best-FWHM initially
  failed: `photometry.py`'s moment-based FWHM is computed over a
  thresholded segmentation mask, which systematically underestimates a
  Gaussian source's true (analytic) FWHM because the mask truncates the
  profile's wings before the second moment is computed. This is a
  property of the existing photometry algorithm, not a refactor bug. Fixed
  the test to check self-consistency instead — comparing `execute()`'s
  fitted vertex FWHM against a direct `image_quality()` call on the frame
  nearest the true best focus — rather than against the synthetic image's
  analytic input FWHM.
- Re-ran the exact real-data command used for manual testing during
  Phase 0 (`focus.py 345 5 -n 4 --obsnum 2165` against the real
  `n2165`-`n2168` sequence) after the refactor: output unchanged
  (`Best focus: 356.1`, `Expected sigma: 3.3 pixels`), confirming no
  behavior regression on real data.

### Testing

12 tests added, all passing in the `nickel` virtualenv (which has no
`ktl` installed):

- `test_quadratic.py`: quadratic fit/vertex correctness and input
  validation.
- `test_focus_sequence.py`: `ArchiveFocusSequence` never constructs
  hardware objects; `step()` yields one correctly-indexed `StepResult`
  per exposure and matches the sequence's own bookkeeping; `step()` can
  be driven one exposure at a time (as a worker thread would); `execute()`
  recovers the known best focus from the synthetic sweep; `execute()`
  with `plot=True` works under the headless `Agg` backend;
  `execute(goto=True)` raises without a `ktl` connection;
  `fit_best_focus()` rejects fewer than 3 points.
- `test_focus_cli.py`: `focus.main()`'s archive-mode CLI path runs
  end-to-end against the synthetic fixture.

### Next

Phase 2 (§9): GUI skeleton + archive mode, plus interactive source
selection and `reanalyze()` (§5.6), all `ktl`-free.
