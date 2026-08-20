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

## Phase 2: GUI skeleton + archive mode

**Design doc reference:** §9 Phased plan, item 2.

**Status:** In progress — sub-phase 1 complete. Broken into sub-phases
below before implementation starts, both because it's a large chunk of
work and because of a
constraint not fully spelled out in `GUI_DESIGN.md`: Qt must stay an
*optional* dependency, exactly like `ktl`. `scripts/focus.py` (the CLI)
must keep working with no Qt installed at all — nothing in `scripts/`
should ever import PySide6, directly or indirectly. That's naturally true
as long as `gui/` only ever imports *from* `scripts/` and never the other
way around, but it also means: `gui/`'s own modules should fail with one
clear, friendly message if PySide6 isn't installed (mirroring the
`ktl is None` pattern) rather than a confusing traceback; the GUI's Qt
dependency belongs in its own requirements file, not the base
`requirements.txt`; and the test suite must keep collecting and passing
on a Qt-free environment, with GUI-specific tests skipped (not erroring)
when PySide6 is absent.

### Sub-phases

1. **`FocusSequence.reanalyze()`** (Qt-free, lives in `scripts/focus.py`).
   Re-runs `photometry.image_quality()` over a sequence's already-collected
   `exposures`/`observed_focus`, in place, with a new `method` — the
   Model-layer piece §5.6's interactive source selection needs, and the
   one piece of Phase 2 with no Qt involvement at all. Testable
   immediately with the existing `tests/` setup.
2. **GUI scaffolding and the optional-Qt boundary.** `gui/__init__.py`;
   a single shim module (e.g. `gui/qt.py`) that imports PySide6 once,
   raising a clear `ImportError` if it's missing, so every other `gui/`
   module gets its Qt classes from that shim rather than importing
   PySide6 directly; the same `scripts/`-on-`sys.path` setup `tests/conftest.py`
   already does; a new `requirements-gui.txt` holding `PySide6` (kept out
   of `requirements.txt`); a `tests/gui/` subtree with its own
   `conftest.py` that calls `pytest.importorskip('PySide6')` so the rest
   of the suite is unaffected if it's not installed. Deliverable: a
   `gui/main.py` that opens an empty window, proving the scaffolding and
   the optional-import boundary both work, before any real panel exists.
3. **`SequenceWorker`** (Qt-dependent). The `QThread` engine from §4.3:
   drives `FocusSequence.step()`/`reanalyze()` off the GUI thread, emits
   `stepComplete`/`sequenceFinished`/`sequenceFailed`, honors "Stop only
   between steps," and handles the early-stop-with-<3-points case. Tested
   against the synthetic `ArchiveFocusSequence` fixture with no real
   widgets involved.
4. **`ImagePanel`.** Embedded-matplotlib view (§3 phase 1 choice): pan via
   scroll bars, wheel zoom, recenter, a small stretch selector, the
   outlier-highlight box, the focus-value-sorted exposure drop-down
   (§5.2), and the click-to-select-source signal (§5.6) — gating on
   whether a sequence is running is the Controller's job, but the panel
   needs an enable/disable hook for it.
5. **`FocusCurvePanel`.** Embedded-matplotlib FWHM-vs-focus scatter with
   the best-fit quadratic/vertex and distinctly-styled outlier points
   (§5.3).
6. **`FocusControlPanel`.** Archive-mode sequence configuration
   (datadir/prefix/suffix/obsnum/start/step/nstep), the current photometry
   method (brightest/weighted/selected-source, read-only in this phase
   except via an image click), Start/Stop, live status/log, and the
   best-focus result display (§5.4). Grid/Automated sequence-type options
   are present in the layout but disabled with a "requires live ktl —
   Phase 3" note, so the widget shape doesn't need to change later.
7. **`MainWindow` + `Controller` wiring.** Lay out the three panels per
   §5.1; the Controller owns the current sequence, wires
   `SequenceWorker`'s signals to the panels, implements the
   hardware-exclusivity state machine from §4.3 (even though there's no
   real hardware yet, so Phase 3 doesn't need to re-architect it), and
   wires image clicks through to `reanalyze()`.
8. **End-to-end smoke test + log update.** A `tests/gui/` test that opens
   `MainWindow` offscreen, runs an archive-mode sequence against the
   synthetic fixture, and confirms the curve/best-focus populate; a
   source-selection → reanalysis smoke test. Update this document with
   what shipped, deviations, and issues, same as Phase 1.

### Sub-phase 1 results: `FocusSequence.reanalyze()`

Added `reanalyze(method='brightest')` to `FocusSequence` (right after
`step()`/before `execute()`), sharing `step()`'s general shape (a
generator yielding `StepResult`s) but iterating over the sequence's
already-collected `observed_focus`/`exposures` directly instead of
calling `continue_sequence()`/`step_focus()`/`take_exposure()` — so it
works unmodified on any subclass (`Archive`, `Grid`, `Automated`) since it
never touches hardware, and degrades cleanly to an empty generator when
nothing's been collected yet. Also had `step()` and `reanalyze()` both set
`self.method`, giving the sequence a live record of "the method currently
in effect" — previously declared in `__init__` but never actually set or
read anywhere; needed for the GUI to know what to show for the current
photometry method and for §5.6's "the newly-selected source persists for
future exposures" behavior.

No deviations from the design doc for this sub-phase.

**Testing:** 4 tests added to `test_focus_sequence.py`: `reanalyze()`
updates `img_quality`/`source_stamps`/`centroids` in place without
touching `observed_focus`/`exposures`/`step_iter`; it's a no-op on a
sequence with nothing collected yet; it updates `self.method`; and
reanalyzing with a coordinate-based `method` pointing at the fixture's
single source recovers the same measurement as `'brightest'` (the
synthetic fixture only has one source per frame, so this doesn't exercise
actually picking a *different* source among several — worth keeping in
mind if a future test wants to check that specifically, e.g. with a
multi-source synthetic frame). All 16 tests pass (12 from Phase 1 + 4
new).

### Sub-phase 2 results: GUI scaffolding and the optional-Qt boundary

Added the `gui/` package and its supporting test infrastructure:

- `gui/qt.py` — the single place PySide6 gets imported. Every other
  `gui` module is expected to get `QtCore`/`QtGui`/`QtWidgets` from here
  (`from gui.qt import QtWidgets`) rather than importing PySide6 directly.
  Catches the broad `ImportError` (not the narrower `ModuleNotFoundError`)
  so the clear, custom message is guaranteed regardless of exactly how the
  import fails.
- `gui/__init__.py` — only adds `scripts/` to `sys.path` (mirroring
  `tests/conftest.py`). Deliberately does *not* import `gui.qt`, so
  `import gui` alone never requires Qt to be installed; only importing a
  submodule that actually needs it does.
- `gui/main.py` — the sub-phase's deliverable: `build_app()`/
  `build_window()`/`main()`, opening an empty, titled `QMainWindow`. Split
  into separate functions (rather than one `main()`) specifically so
  tests can construct and inspect the window without calling the
  blocking `app.exec()`.
- `requirements-gui.txt` — holds `PySide6`, kept separate from
  `requirements.txt` so installing the base file never pulls in Qt.
- `tests/gui/conftest.py` — calls `pytest.importorskip('PySide6')` so the
  whole `tests/gui/` subtree skips as a unit if Qt isn't installed; also
  sets `QT_QPA_PLATFORM=offscreen` (so GUI tests never need a real
  display) and provides a session-scoped `qapp` fixture.
- `tests/test_optional_qt_boundary.py` — a cheap, Qt-free regression
  guard that scans `scripts/*.py` for any reference to a Qt binding name
  and fails if it finds one. Runs unconditionally, as part of the base
  suite.

### Issues encountered

- The `nickel` environment already has `pytest-qt` installed (unrelated
  to this project's requirements files). `pytest-qt` registers a plugin
  that unconditionally imports Qt in `pytest_configure`, *before* any
  test or conftest-level `importorskip` runs — so simulating a Qt-free
  environment by hiding `PySide6` (via `sys.modules['PySide6'] = None`)
  while `pytest-qt` is active makes pytest itself crash with an
  `INTERNALERROR`, not a clean skip. This isn't a bug in this repo's
  code: a genuinely Qt-free environment (only `requirements.txt`
  installed, no `requirements-gui.txt`) would never have `pytest-qt`
  installed either, so its plugin would never load and this can't happen.
  Confirmed by re-running the same simulation with `-p no:pytest-qt`
  (i.e., modeling the actual target environment): 17 passed, and the
  entire `tests/gui/` subtree (3 tests) skipped as one unit, exactly as
  intended. **Caveat for later:** if `pytest-qt` ends up genuinely useful
  for a future sub-phase (e.g., `qtbot` for simulating clicks or waiting
  on signals), it must be added to `requirements-gui.txt`, specifically
  so it's never present without Qt itself, and any dedicated Qt-free CI
  leg must be provisioned from `requirements.txt` alone. Not adopted yet
  — this sub-phase's own tests use a small hand-written `qapp` fixture
  instead, since nothing here needed more than that.

No other deviations from the design doc for this sub-phase.

**Testing:** 4 tests added: `tests/gui/test_main.py` (the empty window
opens and its title is set) and `tests/gui/test_qt_shim.py` (`gui.qt`
exposes the expected names; simulating `PySide6` being absent via
`sys.modules['PySide6'] = None` makes `import gui.qt` raise a clear
`ImportError` mentioning PySide6), plus `tests/test_optional_qt_boundary.py`
(1 test, Qt-free). All 20 tests pass with PySide6 installed; 17 pass (3
cleanly skipped) with it simulated as absent.

### Sub-phase 3 results: `SequenceWorker`

Added `gui/model/__init__.py` and `gui/model/sequence_worker.py`, a
`QThread` subclass driving a `focus.FocusSequence`:

- `SequenceWorker(sequence, method='brightest', mode='step')` —
  `mode='step'` drives `sequence.step(method=...)` (a live/archive run);
  `mode='reanalyze'` drives `sequence.reanalyze(method=...)` (re-measure
  what's already been collected, no new exposures). Both modes share the
  same signal contract, since a real workflow should be able to refresh
  its best-focus result after a reanalysis the same way it does after a
  normal run.
- Signals: `stepComplete(object)` (one `StepResult` per yield),
  `sequenceFinished(float, float)` (`best_focus, best_fwhm`, once),
  `sequenceFailed(str)` (once, human-readable) — the latter two are
  mutually exclusive per run.
- `request_stop()` sets a plain flag checked only *after* a `StepResult`
  is yielded, matching §4.3's "Stop only between steps, never
  mid-exposure" decision.
- Beyond what the design doc's pseudocode showed, wrapped the whole
  step/reanalyze loop in a broad `except Exception` (not just the
  narrower `fit_best_focus()` call) so a genuine mid-sequence failure
  (e.g., photometry raising because a frame has no detectable source)
  reports `sequenceFailed` instead of the thread dying silently, which
  Qt would otherwise not surface as a crash -- just a `QThread` that
  quietly stopped emitting anything, which would leave a Controller
  waiting forever for a signal that's never coming. Kept the "too few
  points for a fit" case as its own narrower `except ValueError` around
  just `fit_best_focus()`, so that specific message stays as
  clear/specific as the design doc intended.

No deviations from the design doc's signal contract; the broader
exception handling above is treated as a robustness addition, not a
deviation, since it doesn't change any documented behavior in the
success/early-stop paths.

**Testing:** 4 tests in `tests/gui/test_sequence_worker.py`, using a
small helper that starts the worker and pumps a local `QEventLoop` until
`QThread.finished` fires (the built-in completion signal, independent of
which of `sequenceFinished`/`sequenceFailed` was emitted), with all test
signal connections made `DirectConnection` so callbacks run synchronously
on the worker thread -- the simplest way to collect results
deterministically without any real widgets involved. Covers: an invalid
`mode` is rejected at construction; `mode='step'` runs the full synthetic
sequence and recovers the known best focus; calling `request_stop()`
from inside a `stepComplete` handler after 2 steps halts the worker at
exactly 2 (not the full 5) and emits `sequenceFailed` (not
`sequenceFinished`) with a message explaining why; `mode='reanalyze'`
re-measures every exposure and updates `sequence.method`. All 24 tests
pass (20 from Phases 1-2 sub-phases 1-2, + 4 new); reran the new tests
several times back-to-back to check for threading-related flakiness --
none observed.

### Sub-phase 4 results: `ImagePanel`

Added `gui/views/__init__.py` and `gui/views/image_panel.py`
(`ImagePanel(QWidget)`), plus pinned matplotlib's Qt backend to PySide6
(`os.environ.setdefault('QT_API', 'pyside6')` in `gui/qt.py`, since
`ImagePanel` is the first module to embed a `FigureCanvasQTAgg`, and
matplotlib's own auto-detection shouldn't be left to guess which binding
to use).

Covers §5.2/§5.6's baseline: an embedded matplotlib canvas showing the
current `StepResult`'s frame with the measured source boxed (red, or
yellow + "Outlier centroid" when `is_outlier`); a stretch selector
(`ZScale`/`Min-Max`, both via `astropy.visualization`); a "Recenter"
button; a single exposure drop-down, always kept sorted by focus value
(via `bisect.bisect_right`, independent of insertion/acquisition order)
and showing both the exposure's filename stem and its focus value (e.g.
`"n2165 — Focus 345.0"` — see Phase 1's `obsnum` deviation note for why
it's the filename stem rather than a literal observation number); a
`sourceSelected(float, float)` signal emitted on click, gated by
`set_selection_enabled()`, which the Controller (sub-phase 7) will drive.

**Deviations from `GUI_DESIGN.md`:**

- **Pan via scroll bars, implemented differently than the literal phrase
  suggests.** Rather than wrapping the canvas in a `QScrollArea` and
  resizing the canvas widget itself to make native scrollbars appear
  (fragile with a matplotlib-rendered canvas, since it means resizing the
  `Figure` in physical inches/DPI terms to simulate a zoom), `ImagePanel`
  uses two real `QScrollBar` widgets whose range/page-step are computed
  from the current zoom level and whose `valueChanged` sets the
  matplotlib axes' `xlim`/`ylim`. This satisfies the requirement (real
  scroll bars drive panning) without the more brittle approach.
- **The "center on measured source" recenter variant was not
  implemented** — the design doc listed it as a "possibly," and
  "Recenter" here only ever fits the whole frame. Worth adding later if
  it turns out to matter in practice.
- **Switching exposures (via the drop-down or a new result arriving)
  always resets the zoom/pan to fit-to-view**, rather than preserving
  whatever zoom/pan the user had. Not explicitly specified either way in
  the design doc; chosen because leftover pan/zoom state from a
  previously-viewed frame seemed more likely to confuse than help,
  especially across frames of different sizes. Worth revisiting if it
  turns out to be annoying in practice (e.g., wanting to stay zoomed in
  on the same region while stepping through several exposures).

**Testing:** 6 tests in `tests/gui/test_image_panel.py`: the drop-down
stays sorted by focus value regardless of the order results are added in,
and "latest" tracks *acquisition* order rather than sorted position (used
a deliberately scrambled, non-monotonic add order to exercise this, the
same scenario `AutomatedFocusSequence` would produce live); manually
selecting an older entry survives a subsequent `add_result()`; clicking
is a no-op until `set_selection_enabled(True)`, then emits the clicked
data coordinates (simulated by constructing a `matplotlib.backend_bases.MouseEvent`
directly and calling `panel._on_click(event)`, rather than driving a real
Qt mouse event through pixel-to-data coordinate translation, which would
test matplotlib's own event system more than this panel's logic); the
outlier box color follows `StepResult.is_outlier` (checked via
`dataclasses.replace()` on a real fixture result, comparing
`patch.get_edgecolor()` against `matplotlib.colors.to_rgba(...)`); zooming
opens up the scroll ranges and `recenter()` collapses them back to zero;
`reset()` clears the dropdown and current result. All 30 tests pass (24
from earlier sub-phases + 6 new).

### Sub-phase 5 results: `FocusCurvePanel`

Added `gui/views/focus_curve_panel.py` (`FocusCurvePanel(QWidget)`): an
embedded-matplotlib scatter of FWHM vs. focus, redrawn on every
`add_result()`. Normal points are blue circles labeled "Measured";
outlier points (`StepResult.is_outlier`) are gold "x" markers labeled
"Outlier" — distinguished by both color and marker shape, per §5.3's
decision. Once 3+ points have been added, overlays the best-fit
quadratic (red curve) and its vertex (a green point labeled
`"Best focus: {value}"`), using the same `quadratic.fit_quadratic`/
`vertex` Model functions the CLI's `FocusPlot` and `execute()` already
use — this panel only owns the plotting, not the math. Wrapped the fit in
a narrow `except (ValueError, ZeroDivisionError)` so a degenerate fit
(e.g., collinear points) skips drawing the curve instead of crashing the
panel; not explicitly called out in the design doc, treated as the same
kind of defensive addition as `SequenceWorker`'s broad exception guard in
sub-phase 3.

Mirrors `ImagePanel`'s `add_result()`/`reset()` shape (no
`show_result()` here, since every point stays on screen at once — there's
no single "current" point the way there's a single current exposure).

No deviations from the design doc for this sub-phase.

**Testing:** 5 tests in `tests/gui/test_focus_curve_panel.py`: points
accumulate and a scatter appears; the fitted curve only appears once 3+
points exist; the plotted vertex's x-coordinate (found by its
`"Best focus"`-prefixed label) matches the synthetic fixture's known best
focus; outlier and normal points are drawn as separate, distinctly
labeled series; `reset()` clears everything. All 35 tests pass (30 from
earlier sub-phases + 5 new).

### Sub-phase 6 results: `FocusControlPanel`

Added `gui/views/focus_control_panel.py` (`FocusControlPanel(QWidget)`).
Purely a View: it exposes configuration and emits signals for requested
actions (`startRequested`, `stopRequested`, `methodChanged`,
`moveToBestFocusRequested`); it never constructs a `FocusSequence` or a
`SequenceWorker` itself. Grouped into: Sequence Type (Archive checked and
enabled, Grid/Automated present but disabled with a "requires live ktl —
Phase 3" tooltip); Archive Configuration (datadir + Browse…, prefix,
suffix, obsnum, start focus, step size, number of steps —
`get_archive_config()` returns these as a `dict` with `datadir` as a
`Path`, matching Phase 1's Path-everywhere decision); Exposure Settings
(exposure time/speed/binning, all disabled — not applicable until Phase 3
takes real exposures); Photometry Method (a Brightest/Weighted combo the
user can pick from directly, `methodChanged`, plus a separate read-only
`method_label` driven only by `set_method()`, which is how "selected
source" ever gets displayed — only `ImagePanel`'s click, routed through
the Controller, can put it in that state, matching §5.4's "read-only
except via an image click"); Start/Stop; a Status group (a status line,
a live step/focus/FWHM label, and a scrolling log); and a Result group
(best-focus/FWHM display and the "Move to Best Focus" button).

**Deviations from `GUI_DESIGN.md`:**

- **"Move to Best Focus" is unconditionally disabled in this phase**,
  not just "disabled until a result exists" as §5.4 describes for the
  finished feature. There's no hardware to move yet, so enabling it once
  a result exists would let the user click a button that does nothing
  meaningful. The confirmation-dialog logic (`_on_move_to_best_focus_clicked`,
  `moveToBestFocusRequested`) is fully implemented and tested now,
  matching the same "build the shape now, flip it on in Phase 3" approach
  already used for Grid/Automated and exposure settings — Phase 3 should
  only need to remove this button's blanket disable, not add new logic.

No other deviations for this sub-phase.

**Testing:** 14 tests in `tests/gui/test_focus_control_panel.py`,
covering: only Archive is selectable; exposure settings are disabled;
`get_archive_config()` round-trips exactly what the form fields hold;
the method combo/`methodChanged`/`get_selected_method()` agree with each
other; `set_method()` renders both a plain method name and an `(x, y)`
tuple correctly; Start/Stop signals fire (Stop needed `set_running(True)`
first, since a disabled `QPushButton` doesn't emit `clicked()` at all —
caught this from a real test failure, not anticipated in advance);
`set_running()`/`set_stopping()` toggle the right widgets; `update_step()`,
`show_best_focus()`, `show_failure()`, and `reset()` all update the
right labels/log; the move-to-best-focus confirmation flow (using
`monkeypatch` on `QMessageBox.question`, since it's normally a blocking
modal dialog) correctly no-ops with no result yet, no-ops on "No", and
emits on "Yes"; the Browse… button (`monkeypatch` on
`QFileDialog.getExistingDirectory`) updates the directory field, and
leaves it alone if the dialog is canceled (empty string). All 49 tests
pass (35 from earlier sub-phases + 14 new).

### Next

Phase 2, sub-phase 7: `MainWindow` + `Controller` wiring.
