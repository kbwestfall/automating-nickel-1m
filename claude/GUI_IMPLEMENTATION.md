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

**Status:** Complete (all 8 sub-phases). Broken into sub-phases below
before implementation started, both because it's a large chunk of work
and because of a
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

### Sub-phase 7 results: `MainWindow` + `Controller` wiring

Added `gui/views/main_window.py` (lays out the three panels per §5.1 in
nested `QSplitter`s — purely structural, no wiring) and
`gui/controller.py` (`Controller(QObject)`), and updated `gui/main.py` to
build the real `MainWindow` + `Controller` instead of sub-phase 2's empty
placeholder window.

`Controller` owns `self.sequence` (the current `focus.FocusSequence`, or
`None`) and `self.worker` (the current `SequenceWorker`, or `None`), and
implements §4.3's hardware-exclusivity state machine via `_set_running()`,
which both `FocusControlPanel.set_running()` (Start/Stop/config-field
states) and `ImagePanel.set_selection_enabled()` key off of — the same
single flag drives both, so they can't drift out of sync.
`start_archive_sequence()` builds the target-focus values the same way
`focus.py`'s own CLI archive-mode path does: constructing a
`GridFocusSequence` purely for its focus-value arithmetic (never touches
hardware, since `ktl` isn't connected) to get the file list to hand to
`ArchiveFocusSequence`. `reanalyze()` and `_on_source_selected()` wire
`ImagePanel.sourceSelected` through to a `mode='reanalyze'`
`SequenceWorker` run.

**A real integration gap surfaced while wiring this up, not anticipated
in the design doc or earlier sub-phases:** `ImagePanel.add_result()` and
`FocusCurvePanel.add_result()` only ever *appended* — reanalyzing an
already-collected exposure would have inserted a second, duplicate entry
for the same file rather than updating the existing one. Added
`update_result()` to both panels (matched by `StepResult.exposure`,
falling back to `add_result()` if the exposure isn't already shown), and
had `Controller._on_step_complete()` branch on `self.worker.mode` to call
the right one. This is exactly the "replacing its stored
FWHM/stamp/centroid" behavior §5.6 describes for reanalysis — it just
hadn't been built into the panels yet, since sub-phases 4-5 only ever
exercised the append path.

No other deviations from the design doc for this sub-phase.

**Testing:** 9 tests in `tests/gui/test_controller.py`, driving a real
`MainWindow` + `Controller` pair (using the same `QEventLoop`-pumping
pattern as sub-phase 3's `SequenceWorker` tests, waiting on
`controller.worker.finished`): a full archive run populates both panels
and the result display with the known best focus; `_set_running(True)`
takes effect synchronously before the worker thread does any real work
(no race to test around); hardware exclusivity blocks a second start
(tested with a fake truthy `worker` sentinel, not real timing);
misconfigured paths report a failure without ever starting a worker;
`stop()`/`reanalyze()`'s no-op guards (nothing running, no sequence
loaded); `set_method()` updates both the controller and the display; and
-- the key regression test for the gap above -- clicking a source after a
full run re-measures every exposure and leaves the panels at the *same*
result count, not double. Also added 3 tests each to
`test_image_panel.py`/`test_focus_curve_panel.py` directly for
`update_result()`'s in-place-replace/redisplay-if-current/fall-back-to-add
behavior. All 63 tests pass; reran the full GUI subtree twice more
back-to-back with no flakiness observed.

### Sub-phase 8 results: end-to-end smoke test + Phase 2 close-out

Added `tests/gui/test_app_smoke.py` with two tests:

- A genuine entry-point smoke test going through `gui.main.build_app()`/
  `build_window()` + `Controller` + `window.show()` — everything prior
  sub-phases tested constructed `MainWindow`/`Controller` directly, never
  the actual `gui/main.py` path a user running `python -m gui.main` would
  take.
- A GUI/Model equivalence check: rather than re-testing `focus.main()`
  itself (already covered by `test_focus_cli.py`), this builds the same
  sequence directly — the same `GridFocusSequence`-for-focus-values +
  `ArchiveFocusSequence` + `execute()` calls both `focus.py`'s CLI and
  `Controller.start_archive_sequence()` make under the hood — and
  confirms the GUI's fitted best focus matches exactly. This is the
  concrete version of §8's testing philosophy: the CLI and GUI share the
  same underlying `FocusSequence` code because of the Phase 1 `step()`
  refactor, so this test is what actually proves that sharing holds at
  the GUI's integration point, rather than assuming it.

**No deviations for this sub-phase.** Final Phase 2 check: reran the full
suite with `PySide6` simulated absent (and `pytest-qt` explicitly
disabled, matching a real Qt-free environment — see sub-phase 2's
caveat) — 17 tests pass, the entire `tests/gui/` subtree (now 6 files)
skips as one unit, confirming the optional-Qt boundary still holds with
the complete `gui/` package in place, not just the scaffolding from
sub-phase 2.

**Testing:** 2 new tests. All 65 tests pass with PySide6 installed; 17
pass (48 cleanly skipped) with it simulated as absent.

### Phase 2 summary

All 8 sub-phases complete. What exists now:

- `gui/qt.py`, `gui/__init__.py` — the optional-Qt import boundary.
- `gui/model/sequence_worker.py` — `SequenceWorker`, the `QThread` engine.
- `gui/views/image_panel.py`, `focus_curve_panel.py`,
  `focus_control_panel.py`, `main_window.py` — the three panels and their
  layout.
- `gui/controller.py`, `gui/main.py` — wiring and the entry point.
- `requirements-gui.txt` — `PySide6`, kept separate from `requirements.txt`.
- `tests/gui/` (`conftest.py` + 8 test files) +
  `tests/test_optional_qt_boundary.py` — 48 GUI tests, all `ktl`-free,
  all passing offscreen.

Functionally, the GUI can today: configure and run an archive/replay
sequence against files on disk; watch the image, stamp-implicit outlier
box, and FWHM curve update live as it runs; browse prior exposures by a
focus-value-sorted drop-down without disrupting a running sequence;
click a source to reanalyze the whole sequence with a new photometry
target, in place; and Stop cleanly between steps. It cannot yet: talk to
real hardware (Grid/Automated sequence types and exposure settings are
present but disabled), take a single exposure, or move the telescope to
the best focus (the button and its confirmation dialog exist but are
unconditionally disabled) — all Phase 3.

Two real gaps were found and fixed only because tests were written
alongside each piece rather than after the fact: the `pytest-qt`/Qt-free
interaction (sub-phase 2) and the append-vs-replace bug in
`ImagePanel`/`FocusCurvePanel` that reanalysis would have silently
duplicated every point through (sub-phase 7). Both are documented above
in their sub-phases as issues/gaps, not just as fixes.

## Phase 3: Live mode

**Design doc reference:** §9 Phased plan, item 3.

**Status:** Planning — broken into sub-phases below before implementation
starts, for the same reason as Phase 2 (large chunk of work with natural
component boundaries), plus one constraint that changes the shape of the
plan itself: **this development machine has no real telescope/camera
hardware access.**

That's more than a testing inconvenience. `GridFocusSequence`/
`AutomatedFocusSequence` already construct fine without `ktl` (Phase 1's
design), but the moment `step_focus()` calls `self._focus.set_to(...)`,
that's `None.set_to(...)` — an `AttributeError` — since `self._focus` is
`None` without `ktl`. So without hardware, none of Phase 3's live-mode
code can be run at all here, not even manually, let alone tested.

The plan's answer to that: a lightweight `FakeFocus`/`FakeExposure` test
double matching the real classes' public interface, injected via
`monkeypatch.setattr(focus, 'Focus', FakeFocus)` — not via new
constructor parameters on `FocusSequence`, honoring the earlier decision
that it shouldn't need injection hooks. That unlocks real, meaningful
automated coverage of `GridFocusSequence`/`AutomatedFocusSequence`'s
stepping logic for the first time (they've had zero coverage since Phase
1 — only `ArchiveFocusSequence` was testable before, since it never
touches hardware). It's a genuine benefit beyond working around the
hardware gap: it's the first real exercise of `AutomatedFocusSequence`'s
adaptive curve-following branches.

What the fake-hardware tests *can't* do: validate against the real KTL
keyword protocol (naming, timing, actual hardware quirks). That already
has open uncertainty flagged in the existing code (the `EXPREC`
state-check TODO, the `'StartX'` exposure-start name pending confirmation
with Will Deich) and can only be resolved at the telescope. Phase 3's
fake-hardware testing validates the *software logic* — stepping,
stopping, exposure configuration, error handling, UI wiring — and that
distinction should stay explicit in each sub-phase's results, not get
blurred into "tested" without qualification.

**A design simplification found while scoping this phase:** `execute()`'s
`goto=True` path is broken -- it calls `self.measure_fwhm(...)`, a method
that doesn't exist anywhere in the codebase, and separately never even
appends its verification exposure to `self.exposures` before trying to
read it back. Rather than just patching that bug in place, the fix is to
build a proper `take_single_exposure(focus_value, method, **exp_kwargs)`
primitive on `FocusSequence` -- move to a focus value, take one exposure,
measure it via `image_quality()`, return a `StepResult` -- since
"move to best focus" is exactly one such call at the fitted best-focus
value, and the same primitive is what §5.5's single-exposure workflow
needs too. Building it once, early, turns both "Move to Best Focus" and
"Take Single Exposure" into thin GUI wiring on top of already-tested
Model logic, rather than two separate pieces of Model work.

### Sub-phases

1. **Fake hardware test double** (`tests/` — Qt-free). `FakeFocus`/
   `FakeExposure` matching `focus.Focus`/`focus.Exposure`'s public
   interface, monkeypatched in for tests. Unlocks real coverage of
   `GridFocusSequence`/`AutomatedFocusSequence`'s stepping logic.
2. **`take_single_exposure()` Model primitive**, tested against sub-phase
   1's fake hardware. `execute()`'s `goto=True` path is rewritten to call
   it, fixing the `measure_fwhm()`/exposure-list bugs described above.
3. **Exposure-settings plumbing.** `SequenceWorker` currently never calls
   `sequence._exposure.cfg.configure(...)` the way the old `execute()`
   did; live mode needs exposure time/speed/binning to actually take
   effect before stepping, regardless of which feature triggers a real
   exposure.
4. **Enable Grid/Automated sequence types in the GUI.** Lift
   `FocusControlPanel`'s blanket-disable on those radios and the exposure
   settings group; add `AutomatedFocusSequence`'s max-steps field;
   `Controller` gets a unified `start_sequence()` branching on the
   selected type instead of `start_archive_sequence()` alone.
5. **Enable "Move to Best Focus."** Lift its blanket-disable; wire
   `moveToBestFocusRequested` to `sequence.take_single_exposure(best_focus,
   method=...)` from sub-phase 2 — genuinely simple now that the hard
   part is already built and tested.
6. **Single-exposure workflow (§5.5).** The GUI action (focus-value
   entry, hardware-exclusivity gating for a standalone action, built on
   the same `take_single_exposure()` primitive) plus the "add to existing
   sequence / start new sequence" choice specific to this workflow.
7. **End-to-end fake-hardware smoke test + log update**, closing out
   Phase 3.

### Sub-phase 1 results: fake hardware test double

Added `tests/fake_hardware.py` (`FakeFocus`, `FakeExposurePath`,
`FakeExposureConfig`, `FakeExposure` — matching `focus.Focus`/
`focus.Exposure`'s public interface, with `expose()` synthesizing a small
Gaussian-source FITS frame instead of talking to a camera) and a
`fake_hardware` fixture in `tests/conftest.py` that monkeypatches
`focus.ktl` (to anything not `None`) and `focus.Focus`/`focus.Exposure`
(to return the shared fakes), following the same
`fwhm_min + curvature*(focus - best_focus)**2` relationship as the
existing `focus_sweep` fixture. This is injected by monkeypatching the
module-level names, not new constructor parameters on `FocusSequence`,
per the earlier decision that it shouldn't need injection hooks.

**A real, previously-undetectable bug was found immediately:**
`AutomatedFocusSequence.__init__` set `self.step = step` (the step-size
argument) — but `step` is also the name of the generator method
`FocusSequence.step()` introduced in Phase 1. The instance attribute
shadowed the inherited method entirely, so `seq.step(method=...)` tried
to call a `float`. This has been broken since Phase 1 and could not have
been caught before now: `AutomatedFocusSequence` was never testable
without live hardware, and the CLI has never exercised it either (`main()`
only ever calls `.execute()`, and even that path was never manually
run against a real automated sequence during Phase 0/1 testing — only
Grid/Archive paths were). Fixed by renaming the attribute to
`step_size` throughout the class; confirmed nothing else in the
codebase referenced the old name.

This is exactly the kind of gap the design doc's testing strategy (§8)
was meant to catch, and exactly why building this fixture came first in
Phase 3 rather than being deferred: had the GUI wiring for Grid/Automated
sequences (sub-phase 4) been built and "tested" only by inspection or
against `ArchiveFocusSequence`, this bug would have shipped invisibly and
only surfaced at the telescope.

**Testing:** 5 tests in `tests/test_live_focus_sequence.py` (Qt-free,
alongside the existing `test_focus_sequence.py`, not under `tests/gui/`,
since none of this touches Qt): the injected fake hardware is really
what `FocusSequence` uses; `GridFocusSequence` commands its exact target
focus values in order and recovers the known best focus; a focus value
below the valid range correctly raises; `AutomatedFocusSequence`'s
adaptive search converges near the known best focus from both above and
below it (exercising both the `direction = 1` and `direction = -1`
branches of `step_focus()` for the first time). All 70 tests pass;
reconfirmed the Qt-free boundary holds with the new fixture in place (22
pass, `tests/gui/` still skips as one unit).

### Sub-phase 2 results: `take_single_exposure()` Model primitive

Added `FocusSequence.take_single_exposure(focus_value, method='brightest',
**exp_kwargs)`: moves to `focus_value`, exposes, measures via
`image_quality()`, and returns a `StepResult` with `index=0` — one
iteration of `step()`'s body, but without touching `step_iter`,
`observed_focus`, `exposures`, or any other sequence bookkeeping, since
it isn't part of a larger sequence. Raises the same
"ktl not connected"-style `ValueError` as `execute(goto=True)` used to if
`self._focus`/`self._exposure` are `None`.

`execute()`'s `goto=True` path is now just:

```python
if goto:
    self.take_single_exposure(best_focus, method=method)
```

replacing the old three lines that called the nonexistent
`measure_fwhm()` and separately never appended their own verification
exposure to `self.exposures` before trying to read it back (the
`self._focus is None` check moved into `take_single_exposure()` itself,
so it isn't duplicated).

No deviations from the plan for this sub-phase — this was scoped and
agreed in the Phase 3 planning discussion itself (see the Phase 3
introduction above), not discovered mid-implementation.

**Testing:** 4 tests in `tests/test_take_single_exposure.py` (Qt-free):
calling it on a sequence with no hardware (`ArchiveFocusSequence`) raises
clearly; against the fake hardware fixture, it moves the focus, returns
a correctly-populated `StepResult`, and leaves sequence bookkeeping
untouched; an out-of-range focus value raises; and -- the regression test
for the fix -- `GridFocusSequence(...).execute(goto=True, ...)` against
fake hardware now runs to completion and actually moves the telescope to
the fitted best focus, instead of crashing. All 74 tests pass; Qt-free
boundary reconfirmed (26 pass, `tests/gui/` still skips as one unit).

### Sub-phase 3 results: exposure-settings plumbing for `SequenceWorker`

`SequenceWorker.__init__` gained an `exp_kwargs` parameter (a `dict` of
`record`/`speed`/`binning`/`exptime`, defaulting to `None` and stored as
`{}` -- an explicit `is None` check, not `exp_kwargs or {}`, per your
note). `run()` now applies it via
`sequence._exposure.cfg.configure(**self.exp_kwargs)` before the
step/reanalyze loop, but only when `mode == 'step'` and
`sequence._exposure is not None` — mirroring exactly what the old
CLI-only `execute()` did (`if self._exposure is not None:
self._exposure.cfg.configure(**exp_kwargs)`), so archive sequences (no
exposure hardware at all) and `mode='reanalyze'` (no new exposures)
silently ignore it rather than erroring or doing something meaningless.
Wrapped in the same `try/except` as the stepping loop, so a bad exposure
setting (e.g. an invalid speed/binning value, which `ExposureConfig.configure()`
already validates and raises on) reports `sequenceFailed` instead of
crashing the thread silently.

No deviations from the plan for this sub-phase.

**Testing:** 4 new tests in `tests/gui/test_sequence_worker.py` (8
total in that file now): `exp_kwargs` defaults to `{}`, not `None`;
passing it against an `ArchiveFocusSequence` (no exposure hardware) is a
harmless no-op rather than an `AttributeError` on `None.cfg`; against
the `fake_hardware` fixture, a `GridFocusSequence` run actually applies
the given exposure time/speed before stepping (checked via
`seq._exposure.cfg.exptime`/`.speed` after the run); and `exp_kwargs` is
correctly ignored in `mode='reanalyze'` (confirmed `cfg.exptime` stays
`None` even when `exp_kwargs` is passed). All 78 tests pass; Qt-free
boundary reconfirmed (26 pass, `tests/gui/` still skips as one unit).

### Sub-phase 4 results: enable Grid/Automated sequence types in the GUI

`FocusControlPanel`: lifted the blanket-disable on `grid_radio`/
`automated_radio` and the exposure-settings group; added a
`maxsteps_spin` field; added `_on_sequence_type_changed()`, wired to all
three radios' `toggled` signals (and called once at the end of `__init__`
for the correct initial state), which enables/disables fields per
selected type -- Archive gets the datadir/prefix/suffix/obsnum fields and
`nstep`, no exposure settings; Grid gets `nstep` and exposure settings,
no archive fields; Automated gets `maxsteps` (not `nstep`) and exposure
settings, no archive fields. `set_running(False)` now calls
`_on_sequence_type_changed()` to restore the *correct* per-type state
rather than just re-enabling everything uniformly. Renamed
`get_archive_config()` to `get_sequence_config()` (now returns every
field regardless of type -- the Controller reads what's relevant) and
added `get_sequence_type()` and `get_exposure_config()`.

`Controller`: renamed `start_archive_sequence()` to `start_sequence()`,
which now branches on `get_sequence_type()` to build an
`ArchiveFocusSequence` (unchanged logic, extracted into
`_build_archive_sequence()`), a `GridFocusSequence`, or an
`AutomatedFocusSequence`, and passes `get_exposure_config()` through to
`SequenceWorker` as `exp_kwargs` for the two live types.
`SequenceWorker`/`_start_worker` gained the `exp_kwargs` pass-through.
Added one explicit check after constructing a live sequence: if
`sequence._focus`/`._exposure` are `None` (no `ktl` connection), report a
clear "no ktl connection is available for a live sequence" failure
*before* ever starting a worker thread, rather than letting the worker
fail later with a cryptic `AttributeError` on `None.set_to(...)` (still
caught by `SequenceWorker`'s existing broad exception guard from
sub-phase 3, but with a much less useful message).

**A real test gap was caught while reviewing this sub-phase, not by a
failure:** the first pass of new tests only exercised the *widget
enabling/disabling* logic in `FocusControlPanel` and `Controller`'s
config-branching in isolation -- nothing actually drove a Grid or
Automated sequence through the full `Controller` against the
`fake_hardware` fixture end-to-end. Added
`test_start_sequence_runs_grid_against_fake_hardware` and
`test_start_sequence_runs_automated_against_fake_hardware` to close
that gap; both passed immediately, but the absence would have meant
sub-phase 4's actual deliverable -- live sequences work through the
GUI -- was never really proven, only its surrounding scaffolding.

No deviations from the plan otherwise.

**Testing:** updated `test_focus_control_panel.py` (renamed/replaced the
old "only Archive is enabled" test with `test_all_sequence_types_are_selectable`,
`test_get_sequence_type`, `test_archive_fields_enabled_only_for_archive`
(checks all three types' field states), `test_get_sequence_config_reflects_form_fields`,
`test_get_exposure_config_reflects_form_fields`, and
`test_set_running_restores_the_correct_per_type_state`); renamed
`Controller.start_archive_sequence()` call sites to `start_sequence()`
throughout `test_controller.py`/`test_app_smoke.py`; added the two
fake-hardware integration tests above plus
`test_start_sequence_without_ktl_reports_clear_failure` (confirms the
new upfront check fires with no worker ever starting, using the real
"no ktl" state of this dev machine -- no fixture needed for that one).
All 84 tests pass; Qt-free boundary reconfirmed (26 pass, `tests/gui/`
still skips as one unit).

### Sub-phase 5 results: enable "Move to Best Focus"

`SequenceWorker` gained a third mode, `'single'`, alongside `'step'`/
`'reanalyze'`: it calls `sequence.take_single_exposure(focus_value,
method=method, **exp_kwargs)` once and emits a new
`singleExposureFinished(object)` signal with the resulting `StepResult`,
instead of driving a generator. It never emits `stepComplete`/
`sequenceFinished` for this mode -- there's no sequence of steps to
report and no curve to fit, matching the earlier design insight that a
confirmation exposure is a single measurement, not sequence data.
Failures (e.g. no hardware, out-of-range focus) go through the existing
`sequenceFailed` signal rather than a new one, since "report this
human-readable failure message" behaves identically regardless of mode.
`mode='single'` requires a `focus_value` constructor argument; omitting
it raises `ValueError` immediately, the same fail-fast style as the
existing mode-string validation.

`FocusControlPanel.show_best_focus()` gained a `can_move` argument
(default `False`) that enables/disables `move_to_best_focus_button`
accordingly and updates its tooltip; `reset()` now explicitly disables
the button again (it should never carry over between runs). Added
`show_confirmation(result)`, which reports a completed "Move to Best
Focus" in the status label and log, mirroring `update_step`/
`show_failure`'s style. The button is also added to `set_running(True)`'s
force-disable list -- hardware exclusivity applies to it exactly like
every other action.

`Controller` gained `move_to_best_focus(focus_value)`, connected to
`moveToBestFocusRequested`: a no-op if something is already running or
no sequence has ever run, otherwise starts a `mode='single'` worker with
the current exposure settings. `_on_sequence_finished` now computes
`can_move` from whether the just-finished sequence actually has hardware
(`sequence._focus`/`._exposure` not `None`) -- the same check already
used to reject starting a live sequence type without `ktl` -- and passes
it to `show_best_focus()`. New `_on_single_exposure_finished` adds the
confirmation result to `ImagePanel` (so the confirmation frame is what's
displayed) and calls `show_confirmation()`; it is *not* added to
`FocusCurvePanel`, since re-fitting the quadratic around a confirmation
point taken at (approximately) the already-fitted vertex would only
muddy which points were actual sequence data.

**One real, pre-existing bug surfaced by the new integration test, not
by inspection:** `set_running(False)` unconditionally cleared
`status_label` (to erase a stale "Stopping — waiting..." message once a
stop actually completes). That happens to run in response to the
worker's `finished` signal, which is queued *after* `sequenceFailed`/
`singleExposureFinished` for the same run -- so any status message set
by those (a failure, or now a move-to-best-focus confirmation) was
being wiped out a moment after being shown. No earlier test caught this
because none exercised a real worker failure or the new confirmation
message through to completion; `test_missing_files_reports_failure`
only exercises the synchronous config-validation failure path, which
never starts a worker at all. Fixed by only clearing `status_label` when
it currently holds the stale "Stopping..." text, leaving a genuine
failure or confirmation message alone.

No deviations from the plan otherwise.

**Testing:** `tests/gui/test_sequence_worker.py` gained
`test_single_mode_requires_a_focus_value`,
`test_single_mode_moves_and_exposes` (against `fake_hardware`, checks
the fake focus was moved, `exp_kwargs` applied, and
`singleExposureFinished` fired with no `stepComplete`/`sequenceFinished`),
and `test_single_mode_reports_failure_without_hardware` (an
`ArchiveFocusSequence` has no hardware, so `take_single_exposure` raises
and `sequenceFailed` fires instead of crashing the worker thread).
`tests/gui/test_focus_control_panel.py`: replaced the old
"stays disabled regardless" test with
`test_show_best_focus_updates_result_and_gates_move_button` (checks both
`can_move` states) and added `test_show_confirmation_updates_status_and_log`;
updated `test_set_running_toggles_widget_states`/
`test_reset_clears_status_log_and_result` to also check the move button.
`tests/gui/test_controller.py` gained
`test_move_to_best_focus_is_a_noop_without_a_sequence`,
`test_move_to_best_focus_blocks_while_something_is_running`, and
`test_move_to_best_focus_runs_against_fake_hardware` (runs a real Grid
sequence to completion against `fake_hardware`, then emits
`moveToBestFocusRequested` with the fitted best focus and confirms the
fake hardware was actually commanded to it, the image panel gained the
confirmation exposure, and the status message survived). This last test
is what caught the `status_label`-clearing bug above. All 91 tests pass;
Qt-free boundary reconfirmed (26 pass, `tests/gui/` still skips as one
unit).

### Sub-phase 6 results: the standalone single-exposure workflow (§5.5)

New "Single Exposure" group in `FocusControlPanel`: a focus-value spin
box, a "Take Single Exposure" button, a label reporting the pending
result (or "No pending exposure"), and an "Add to Existing Sequence"
button (disabled until there's both a pending result and a loaded
sequence to add it to -- mirroring `move_to_best_focus_button`'s
`can_move`-gating pattern from sub-phase 5). Deliberately did **not**
add a separate "Start New Sequence" button/signal: per §5.5, starting a
new sequence while a single exposure is pending should simply discard
it, and that's exactly what the *existing* "Start" button already does
once `Controller.start_sequence()` is taught to clear pending state --
introducing a second action that does the same thing as "Start" would
just be two names for one behavior.

`SequenceWorker`'s `mode='single'` (built in sub-phase 5) is reused
as-is -- no changes needed there, since "take a single exposure" is
identical hardware-wise whether it's a best-focus confirmation or a
standalone confirmation shot; only what the Controller does with the
result differs.

`Controller` gained:
- `pending_result` / `_standalone_sequence` state: a not-yet-committed
  `StepResult` and the throwaway `focus.FocusSequence()` used to take
  it. A *throwaway* sequence, not `self.sequence`, deliberately -- a
  standalone exposure has no relationship to whatever sequence (if any)
  is currently loaded, and using a fresh `focus.FocusSequence()` just
  for its auto-detected `_focus`/`_exposure` handles (never subclassed,
  since `take_single_exposure()`/`reanalyze()` don't touch
  `continue_sequence()`/`step_focus()`) means it's never confused with a
  real Grid/Automated/Archive sequence.
- `take_single_exposure(focus_value)`: hardware-exclusivity- and
  `ktl`-connection-gated (identical style to `start_sequence()`'s
  upfront check), starts a `mode='single'` worker against the throwaway
  sequence.
- `add_pending_to_sequence()`: manually appends `pending_result`'s data
  into `self.sequence`'s bookkeeping arrays (the same fields
  `step()`/`reanalyze()` populate), recomputes `is_outlier` against the
  now-larger `centroids` list via `focus.FocusPlot.is_outlier` (it's
  always `False` fresh out of `take_single_exposure()`, since that has
  no sequence context to compare against), and adds the result to
  `FocusCurvePanel` (refitting the curve) -- `ImagePanel` already has it
  from when the exposure was taken, so only `update_result()` to reflect
  the corrected outlier styling.
- `_on_single_exposure_finished` now distinguishes "Move to Best Focus"
  from a standalone exposure by checking `worker.sequence is
  self.sequence` (true only for the sub-phase 5 case, which reuses the
  loaded sequence's hardware handles) rather than adding a redundant
  purpose flag -- a standalone exposure's worker always runs against
  the throwaway sequence, which is never `self.sequence`.
- `reanalyze()` now targets `self.sequence` if loaded, else
  `_standalone_sequence` -- implementing §5.6's "reanalyze() runs
  against that one exposure" requirement for a standalone exposure with
  no sequence loaded, via the *same* `mode='reanalyze'` worker path
  already used for real sequences, with no new Model code.
- `_on_step_complete` distinguishes reanalyzing the standalone target
  (via the same `worker.sequence is` identity check) to update
  `pending_result` and `ImagePanel` in place, but -- consistent with a
  standalone exposure never appearing on the curve until committed --
  does *not* touch `FocusCurvePanel` for it.
- `start_sequence()` now calls a new `_clear_pending()` after a sequence
  successfully starts (not on a failed attempt -- a bad config shouldn't
  silently discard a pending exposure the user hasn't decided about
  yet), matching §5.5's "not automatically counted... discarded when a
  new sequence starts."

No deviations from the plan otherwise, though the "Start New Sequence"
simplification above is a deliberate divergence from the design doc's
literal wording of two separate buttons -- worth folding back into
`GUI_DESIGN.md` during the end-of-development reconciliation pass.

**Testing:** `tests/gui/test_focus_control_panel.py` gained tests for
the new widgets: signal emission for both new buttons,
`show_pending_exposure`/`clear_pending_exposure`, `reset()` clearing
pending state, and `set_running()` locking out/restoring the new
widgets (mirroring the `move_to_best_focus_button` pattern -- re-enables
on stop except `add_to_sequence_button`, which stays gated).
`tests/gui/test_controller.py` gained hardware-exclusivity and
no-`ktl` tests for `take_single_exposure()`, an integration test against
`fake_hardware` confirming a standalone exposure is taken, displayed,
held pending, and does *not* touch `sequence`/the curve panel; tests for
`add_pending_to_sequence()` being a no-op without a pending result or a
loaded sequence and, positively, actually committing into a loaded
Grid sequence's data (checked against real fake-hardware-driven
exposure counts and curve-panel point counts); a test confirming
interactive source selection reanalyzes a standalone pending exposure
in place (fresh measurement, no duplicate in `ImagePanel`, still absent
from the curve); and two tests distinguishing that starting a new
sequence discards a pending exposure only on success, not on a failed
config. All 106 tests pass; Qt-free boundary reconfirmed (26 pass,
`tests/gui/` still skips as one unit).

### Sub-phase 7 results: end-to-end fake-hardware smoke test + Phase 3 close-out

Added `tests/gui/test_phase3_smoke.py` with one test,
`test_phase3_full_live_workflow`, that drives a single `Controller`
through every live-mode action built across sub-phases 1-6 back to back,
in the order GUI_DESIGN.md §5.5's own "typical use" narrative describes
it: confirm the field with a standalone exposure, mark a star on it via
a click (reanalyzing that one exposure since nothing else is loaded
yet), start a Grid sequence that measures the marked star, move to its
best focus, take another standalone exposure, mark a *different* source
(now reanalyzing the *loaded sequence* instead, since it takes priority
once one exists), commit the pending exposure into the sequence, and
finally discard one last pending exposure by starting a fresh Automated
sequence. Every prior sub-phase's tests exercise these actions
individually or in pairs; nothing before this checked that all of them
compose correctly in one continuous session.

**This actually caught a real bug in the test itself, not the
implementation:** the first draft of this test assumed clicking a
source while *both* a loaded sequence and a pending standalone exposure
exist would reanalyze the standalone one (matching §5.6's "no sequence
loaded" wording taken too literally). Running it immediately failed an
identity check -- `pending_result` hadn't changed -- which is actually
`Controller.reanalyze()` working exactly as implemented and tested in
sub-phase 6: it always prefers `self.sequence` over `_standalone_sequence`
when both exist. The design doc's §5.6 sentence describes only the
literally-nothing-loaded case; it doesn't say what happens when both
exist simultaneously, which sub-phase 6 resolved by giving the loaded
sequence priority. Fixed the test to assert the actual (already
deliberate, already covered by unit tests) priority rule instead of
rewriting the Controller to match a misreading -- a good example of a
smoke test needing correcting against already-settled behavior, not the
other way around.

No deviations from the plan otherwise. Reran the full suite with
`PySide6` simulated absent (`-p no:pytest-qt`, matching a real Qt-free
environment) -- 26 pass, the entire `tests/gui/` subtree (now 9 files)
skips as one unit (`conftest.py`'s `pytest.importorskip('PySide6')`
fails at collection, so pytest reports it as a single skip rather than
one per test), confirming the optional-Qt boundary holds with the
complete Phase 3 feature set in place.

**Testing:** 1 new end-to-end test. All 107 tests pass with PySide6
installed; 26 pass (`tests/gui/` skipped as one unit) with it simulated
as absent.

### Phase 3 summary

All 7 sub-phases complete. What exists now, on top of Phase 2:

- `tests/fake_hardware.py` + the `fake_hardware` fixture in
  `tests/conftest.py` -- the test double standing in for real
  telescope/camera hardware on a development machine with no access to
  either, used throughout Phase 3.
- `focus.FocusSequence.take_single_exposure()` -- the shared primitive
  behind both "Move to Best Focus" and the standalone single-exposure
  workflow, per the reassessment that a confirmation move is just one
  step of a sequence taken in isolation.
- `SequenceWorker`'s `mode='single'` and `singleExposureFinished` signal.
- `FocusControlPanel`: Grid/Automated are now fully selectable and
  configured (exposure settings, `maxsteps`); "Move to Best Focus" is
  live, gated on the finished sequence actually having hardware to move;
  a new "Single Exposure" group drives the standalone workflow.
- `Controller`: `start_sequence()` branches over all three sequence
  types against real (or fake) hardware; `move_to_best_focus()`,
  `take_single_exposure()`, and `add_pending_to_sequence()` round out
  the live-mode action set; `reanalyze()` transparently targets whatever
  is actually loaded -- a real sequence, or a standalone pending
  exposure if nothing else is.

Functionally, the GUI can now do everything GUI_DESIGN.md scoped for
v1 except real telescope/camera I/O, which is untestable on this
machine by construction -- every action has been proven against
`fake_hardware` instead, which validates the software logic (stepping,
exclusivity, bookkeeping, worker threading) but not the real KTL
keyword protocol itself; that can only be checked at the telescope.

Two real gaps were caught only because tests were written alongside
each piece rather than assumed correct: `set_running(False)` silently
erasing a just-shown failure/confirmation message (sub-phase 5), and
the missing Grid/Automated-through-Controller integration test
(sub-phase 4). This sub-phase's own smoke test additionally caught a
mistaken assumption in *its own first draft*, corrected against
already-settled, already-tested Controller behavior rather than the
other way around -- see above.

Remaining before GUI development as a whole is considered finished:
reconciling `GUI_DESIGN.md` with the as-built implementation (explicitly
deferred until now, per earlier discussion) -- most notably §5.5's
"Start New Sequence" button, which sub-phase 6 deliberately never built
as a separate action since the existing "Start" button already does the
same thing.

## Post-Phase-3 revisions

Fixes to real usage problems found by actually running the GUI on a
physical display, ahead of Phase 4 -- something offscreen/fake-hardware
testing can prove correct in isolation but can't surface (an
`availableGeometry()`/layout-minimum-size interaction only shows up
against a screen's *actual* resolution).

### Window opens larger than the screen

**Problem:** `MainWindow` unconditionally resized to a fixed 1200x800.
On a smaller/lower-resolution screen, that's larger than the available
desktop area, and since the initial position left the bottom-right
resize handle off screen, the window couldn't be resized smaller either
-- a genuine dead end, not just an inconvenience.

**Fix:** `MainWindow._size_to_screen()` now sizes to the smaller of
1200x800 and 90% of `screen().availableGeometry()` (excluding
docks/taskbars), then centers the window -- so the resize handles are
always reachable regardless of dock position. Also wrapped
`FocusControlPanel` in a `QScrollArea` in `MainWindow.__init__`: without
it, the fix above didn't actually work, because the panel's stacked
group boxes (Sequence Type, Sequence Configuration, Exposure Settings,
Photometry Method, Single Exposure, Status, Result) impose a combined
minimum-height layout constraint taller than 720px (90% of the offscreen
test screen's 800px height), and Qt's layout engine silently overrides
`resize()` to respect a widget's minimum size hint -- so the window
still opened oversized until the panel itself no longer bubbled its
full minimum height up to the window. A `QScrollArea`'s own minimum size
hint doesn't inherit its contents', which is what actually lets the
outer window shrink; the control panel remains fully usable, just
scrollable if the window is small.

**Testing:** `tests/gui/test_main.py` gained
`test_build_window_fits_within_the_screen`, asserting the built,
shown window's size never exceeds `screen().availableGeometry()`. This
test is what caught the layout-minimum-size issue above -- the naive
resize-based fix alone still failed it. All 108 tests pass; Qt-free
boundary reconfirmed.

### Default window size

Changed the preferred size passed to `_size_to_screen()` from 1200x800
to 1200x600 (2:1), per explicit request -- no functional change beyond
the aspect ratio.

### `ImagePanel` usability fixes

Four fixes from actually using the GUI on a real display:

- **Excess whitespace around the image:** matplotlib reserves default
  margins around an `Axes` for ticks/labels/title -- all of which
  `ImagePanel` already suppresses -- so those margins were pure wasted
  space. Added `figure.subplots_adjust(left=0, right=1, bottom=0,
  top=1)` once in `__init__` (a `Figure`-level setting, unaffected by
  `ax.clear()` in `_render`, so it doesn't need repeating there) to let
  the image fill the whole canvas.
- **Zoom didn't track the cursor:** `_on_wheel_zoom` called the public
  `set_zoom()`, which recomputes scrollbar *ranges* for the new zoom
  level but leaves their *value* (the view's top-left corner) alone --
  so the view always zoomed toward that fixed corner, not wherever the
  cursor was. Noticeable with any scroll input, but especially with a
  trackpad's many small high-frequency scroll events, where the drift
  compounds visibly. Added `_zoom_at(xdata, ydata, new_zoom)`, used only
  by the wheel handler: it captures the cursor's data point as a
  fraction of the *current* view before changing `_zoom`, then
  recomputes the scrollbar values so that same fraction -- and thus the
  same data point -- lands back under the cursor at the new zoom level.
  `set_zoom()` itself is unchanged (still used by `recenter()`, where
  there's no cursor position to anchor to anyway).
- **Redundant title:** removed `ax.set_title(...)` from `_render` -- the
  exposure/focus value it showed is already in the drop-down above the
  image, per the request.
- **Default colormap:** changed `imshow`'s hardcoded `cmap` from
  `'viridis'` to `'gray'`, matching the conventional grayscale display
  for astronomical images.

**Testing:** added `test_wheel_zoom_keeps_the_cursor_point_fixed` to
`tests/gui/test_image_panel.py`, which scrolls in on a point away from
the image center and confirms its fractional position within the view
is essentially unchanged before/after -- this is what would have caught
the original corner-anchored behavior. No test was added for the
colormap/title/whitespace changes -- they're cosmetic rendering details
with no behavior to regress. All 109 tests pass; Qt-free boundary
reconfirmed.

### Source selection: click → 'm' key

**Problem:** on macOS, clicking an unfocused window both refocuses it
*and* delivers that same click to whatever's underneath the cursor. A
GUI window that had lost focus (e.g. after switching to another app)
needed exactly the kind of click that, over `ImagePanel`, used to also
mean "select this source" -- so simply refocusing the window could
silently kick off a reanalysis nobody asked for.

**Fix:** replaced the `button_press_event` connection with
`key_press_event`, and `_on_click` with `_on_key_press`, which only acts
on `event.key == 'm'` -- otherwise unchanged (same
`_selection_enabled`/`inaxes`/`xdata`/`ydata` checks, same
`sourceSelected` emission). The intended gesture is now "hover over the
star, press 'm'" rather than "click the star." Also explicitly set
`canvas.setFocusPolicy(QtCore.Qt.FocusPolicy.StrongFocus)`, needed for
the canvas to actually receive key events at all -- a plain click still
naturally gives it keyboard focus as a side effect (standard Qt
behavior for a `StrongFocus` widget), so no extra step is needed before
pressing 'm' works.

**Testing:** replaced `test_click_emits_signal_only_when_selection_enabled`
with `test_m_key_emits_signal_only_when_selection_enabled`,
`test_other_keys_do_not_select_a_source`, and
`test_plain_click_does_not_select_a_source` -- the last one dispatches a
real `button_press_event` through the canvas's own callback registry
(not just checking a method no longer exists) to prove nothing is
listening for clicks at all anymore. All 111 tests pass; Qt-free
boundary reconfirmed.

### Control panel width balloons and Archive/Replay stops responding

**Reported symptom:** on this machine's real display, hitting "Start"
for an Archive/Replay sequence made the bottom-right control panel's
*contents* dramatically expand in width (a horizontal scrollbar
appeared; buttons got much wider), while the window and splitter
positions stayed the same size, and the sequence did not appear to run.
No error appeared anywhere -- not the terminal, not the Status label.

**Investigation:** confirmed via a headless (`QT_QPA_PLATFORM=offscreen`)
script that `Controller.start_sequence()` and `ArchiveFocusSequence`
itself work fine in isolation (5/5 steps complete, no exception) even
with a deliberately small window, which rules out a functional
regression in the archive-loading/stepping code -- nothing about the
recent 'm'-key or window-sizing changes touches that path. Separately,
confirmed a genuine `QScrollArea.setWidgetResizable(True)` bug in
isolation: if its viewport is squeezed narrower than the scrolled
widget's own minimum width, the widget doesn't get a horizontal
scrollbar as you'd expect -- it snaps to some unrelated, much larger
width (a ~340px-minimum panel jumped to ~640px in testing) and gets
stuck there, matching the reported symptom closely. This is a plausible
consequence of the sub-phase-adjacent "window too large" fix, which
introduced the `QScrollArea` wrapping `FocusControlPanel` in the first
place -- a bare `QSplitter` child (the previous arrangement) can be
compressed below its minimum size hint via handle-dragging without this
failure mode; a resizable `QScrollArea` apparently cannot.

**Caveat -- this is a mitigation, not a confirmed root-cause fix:**
despite reproducing the underlying `QScrollArea` bug in isolation, it
could not be reproduced through the actual embedded `MainWindow`
splitter hierarchy in the headless test environment -- `QMainWindow`
enforces the aggregate minimum size of its central widget on ordinary
resize, so `window.resize()` alone never actually squeezed the scroll
area's viewport below the panel's minimum in testing, even deliberately
far below (splitter `setSizes()` calls that did force it narrow also
didn't trigger the bug there). It's possible the real trigger on a
physical display involves something offscreen mode doesn't reproduce
faithfully (e.g. manually dragging a splitter handle, real font
metrics/DPI, or a timing/event-order difference). Applied the fix
regardless, since it's a correct, low-risk hardening either way: giving
the `QScrollArea` an explicit `setMinimumWidth()` -- a hard floor,
unlike a sizeHint -- equal to the panel's own `minimumSizeHint().width()`
plus a small buffer for the vertical scrollbar, so the splitter can
never ask the scroll area to go narrower than the panel can actually
provide, which is exactly the precondition the bug needs to trigger.

**Testing:** added `test_control_scroll_area_has_a_floor_at_the_panels_minimum_width`,
which checks the configured floor directly rather than trying to
reproduce the elusive rendering bug end-to-end. All 112 tests pass;
Qt-free boundary reconfirmed. **Flagged for the user to confirm** this
actually resolves what they saw on their real display, since it
couldn't be independently verified end-to-end here.

Confirmed by the user: reproducing the exact same steps that triggered
this no longer shows the problem.

### Drop-down selection preserves zoom/center

Picking a different exposure from `ImagePanel`'s drop-down used to
always reset to a fit-to-view zoom via `show_result()`. Per request,
`_on_combo_changed` now calls a new `_show_result_preserving_view()`
instead: it captures the current view's center in data coordinates
(from the scrollbar values and `_view_extent()`), swaps in the new
result's frame and shape, recomputes scrollbar ranges for that frame,
then re-centers the view on the same data point at the same zoom
(clamped to the new frame's valid scroll range). Falls back to the
existing `show_result()` (fit-to-view) if nothing was displayed yet.
`add_result()`'s own call to `show_result()` for a newly-collected
exposure is unchanged -- resetting to fit-to-view still makes sense
there; this only applies to manually browsing already-collected frames.

No deviations. **Testing:** added
`test_dropdown_selection_preserves_zoom_and_center`, which zooms in and
pans, switches to a different drop-down entry, and confirms both the
zoom factor and the rendered view limits (center) are unchanged. All
113 tests pass; Qt-free boundary reconfirmed.

### Focus curve panel: clipped x-axis label

`FocusCurvePanel`'s `Figure` used matplotlib's default fixed subplot
margins, sized for a "normal" standalone plot window -- too little
bottom margin once the panel is squeezed short by the splitter, which
clipped the "Focus Value" x-axis label. Fixed by constructing the
`Figure` with `layout='constrained'`, which recomputes margins on every
draw based on the actual rendered size and label text, rather than a
one-time fixed fraction -- so the label always has room regardless of
how small the panel gets. `ImagePanel` didn't need the equivalent
treatment: it draws no ticks/labels/title at all (removed in an earlier
revision), so it has nothing that constrained layout would need to make
room for.

**Testing:** added `test_xlabel_is_not_clipped_when_the_panel_is_short`,
which resizes the panel to a deliberately short 300x150 and checks the
rendered x-axis label's bounding box stays within the figure's bounds.
Confirmed this fails without the fix (label bbox extends ~29px below
the figure) and passes with it. All 114 tests pass; Qt-free boundary
reconfirmed.

### Focus control panel: no spinner buttons, integer focus values

Two changes to `FocusControlPanel`, per request:

- Every spin box (`obsnum_spin`, `start_spin`, `step_spin`, `nstep_spin`,
  `maxsteps_spin`, `exptime_spin`, `single_focus_spin`) now has
  `setButtonSymbols(QAbstractSpinBox.ButtonSymbols.NoButtons)`, removing
  the increment/decrement arrows -- typing/pasting a value still works
  identically.
- `start_spin` ("Start focus") and `single_focus_spin` ("Focus value",
  single-exposure workflow) changed from `QDoubleSpinBox` to `QSpinBox`
  (range unchanged, 165-500): a focus value is a whole-unit telescope
  position, even though a best-fit focus (from `fit_best_focus`'s
  quadratic vertex) can land on a fractional value. `step_spin` (a step
  *size*, not a target position) and `exptime_spin` (exposure time, not
  a focus value at all) are unaffected -- the request was specifically
  about focus *value* entries, not every numeric field.
- `_on_move_to_best_focus_clicked` now rounds `self._best_focus` via
  Python's `round()` (not truncation) *before* both building the
  confirmation dialog's text and emitting `moveToBestFocusRequested` --
  so the dialog shows the observer exactly the integer value that will
  actually be commanded, and a fractional fit result (e.g. 356.6) never
  reaches `Focus.set_to()`.

No deviations. **Testing:** added
`test_number_entry_boxes_have_no_increment_buttons`,
`test_focus_value_entries_are_integers`,
`test_move_to_best_focus_rounds_up_not_just_truncates` (356.6 -> 357,
chosen specifically to distinguish rounding from truncation, which
would give 356), and `test_move_to_best_focus_dialog_shows_the_rounded_value`.
Updated `test_move_to_best_focus_confirmation_flow` and
`test_take_single_exposure_signal`, whose old expectations assumed
float values. All 118 tests pass; Qt-free boundary reconfirmed.

Follow-up: `step_spin` ("Step size") changed from `QDoubleSpinBox` to
`QSpinBox` (range 1-100) too, per request -- it's expressed in the same
focus-position units as `start_spin`/`single_focus_spin`, so the same
"whole units" reasoning applies. Updated
`test_get_sequence_config_reflects_form_fields` (now asserts an int
`step`) and `test_focus_value_entries_are_integers` (now also covers
`step_spin`). All 118 tests still pass.

While this was in progress, the user made their own matching edit to
`ImagePanel.add_result()`/`update_result()`, changing the drop-down
label's focus-value format from `.1f` to `.0f`. Asked whether other
spots needed the same treatment: found every other place
`FocusControlPanel` displays a `StepResult.focus_value` readout
(`update_step`'s step label/log, `show_confirmation`, and
`show_pending_exposure`'s label/log) still used `.1f`, and updated all
of them to `.0f` for consistency -- these are all the same kind of
per-exposure focus readout as the drop-down label. Deliberately left
`show_best_focus`'s `best_focus`/`best_fwhm` display as `.1f`: that's
the fitted quadratic's continuous optimum, not an entry or readout of
an actual commanded/measured position, and showing its sub-integer
precision is informative (the rounded integer *target* is what's
already shown separately, in the move-to-best-focus confirmation
dialog). Updated three tests
(`test_update_step_updates_label_and_log`,
`test_show_confirmation_updates_status_and_log`,
`test_show_pending_exposure_updates_label_and_gates_add_button`) plus
the two integration tests checking the confirmation status text
(`test_move_to_best_focus_runs_against_fake_hardware`,
`test_phase3_full_live_workflow`) that asserted the old `.1f` text. All
118 tests still pass; Qt-free boundary reconfirmed.

## Phase 4: Focus control panel redesign (tabs)

Per request: replace `FocusControlPanel`'s flat, always-fully-laid-out
structure (every sequence type's fields present at once, shown/hidden
by enabling/disabling) with a `QTabWidget`: **Single, Grid, Auto,
Replay, Log, Help**. Each tab shows only the fields relevant to that
action -- the actual fix for "minimize white space," since disabled
fields still consume layout space but a hidden tab page doesn't.

Design decisions, confirmed with the user before starting:

- **Single tab absorbs "Move to Best Focus" entirely.** Its focus field
  defaults to the most recent fit's best focus (rounded to an integer);
  acquiring a single exposure there with the default value *is* moving
  to best focus, and changing the value first tests any other focus.
  There is no separate confirmation dialog -- reviewing/editing the
  value before pressing Start *is* the confirmation. This removes
  `moveToBestFocusRequested`, `move_to_best_focus_button`, and its
  `QMessageBox` confirmation entirely.
- **"Add to Existing Sequence" is dropped.** A standalone single
  exposure is never folded into another sequence's curve -- simplest,
  and matches the Single tab's narrower new role. This removes
  `addToSequenceRequested`, `pending_result`, `add_to_sequence_button`,
  `show_pending_exposure()`/`clear_pending_exposure()`. Interactive
  source selection ('m' key) reanalyzing a standalone single exposure in
  place when no sequence is loaded (§5.6) is *not* dropped -- that's a
  different feature (marking a source), orthogonal to committing data
  into a sequence -- so `Controller.reanalyze()`'s existing
  self.sequence-else-`_standalone_sequence` targeting is kept as-is.
- **Photometry method and Start/Stop are shared, below the tabs**, not
  duplicated on each one -- they mean the same thing regardless of which
  tab is active. "Start" reads whichever tab is currently selected
  (`get_sequence_type()`/`get_sequence_config()`), except on Single,
  where it emits `takeSingleExposureRequested` directly instead.
- **Live progress lives on the Log tab only** -- no persistent status
  line outside it. The Log tab replaces the old "Status" group
  one-for-one (`status_label`, `step_label`, `log_widget`, unchanged
  internally); there is no more separate "Result" group, since the
  fitted best focus is now conveyed by the Single tab's default value
  (plus the existing "Sequence finished..." log line) rather than a
  dedicated label.
- The tab *bar* is never disabled while something is running (only the
  input widgets within each tab are) -- otherwise the Log tab wouldn't
  be reachable while a sequence runs, defeating the point above.

### Sub-phases

1. **Rewrite `FocusControlPanel`** as the tabbed view described above,
   with its own updated unit tests (`tests/gui/test_focus_control_panel.py`
   rewritten to match -- most existing tests reference removed
   radios/groups and won't survive as-is). No `Controller`/`MainWindow`
   changes yet, so this is verifiable in isolation. Includes small
   shared helpers (`_int_spin`/`_float_spin`/exposure-fields-row
   builder) reused across Single/Grid/Auto tabs, per the "minimize
   replicated code" request -- Grid and Auto in particular are identical
   apart from `nstep` vs. `maxsteps`.
2. **Update `Controller`** to match: drop `move_to_best_focus()`/
   `add_pending_to_sequence()`/their signal connections and the
   `_single_purpose`-style branching in `_on_single_exposure_finished`
   they required (only one single-exposure flow remains); wire
   `show_best_focus()`'s simplified signature (no more `can_move`).
   Update `tests/gui/test_controller.py` accordingly.
3. **Update the end-to-end tests** (`tests/gui/test_app_smoke.py`,
   `tests/gui/test_phase3_smoke.py`) to drive the new tabbed workflow,
   run the full suite plus the Qt-free boundary check, and log results/
   deviations here.

`scripts/focus.py` and `gui/model/sequence_worker.py` need no changes at
all -- this is purely a View/Controller reorganization; `take_single_exposure()`
and `mode='single'` already do exactly what the unified Single tab needs.

### Sub-phase 1 results: rewrite `FocusControlPanel`

`gui/views/focus_control_panel.py` rewritten per the design above.
Structure: `QTabWidget` (Single/Grid/Auto/Replay/Log/Help) +
`_build_method_group()` + Start/Stop row, all inside one `QVBoxLayout`.
Each tab is its own plain `QWidget` with a single `QFormLayout` (or
`QVBoxLayout` for Log/Help) -- no nested `QGroupBox` wrappers, since the
tab title already says what the fields are for; that redundant nesting
was part of the old whitespace problem.

Shared helpers, per the "minimize replicated code" request:
`_int_spin()`/`_float_spin()` (module-level, no `self` needed) wrap
`QSpinBox`/`QDoubleSpinBox` construction with the no-buttons behavior
from the last round of revisions, used by every numeric field across
all four tabs; `_add_exposure_rows(form)` adds the
exptime/speed/binning rows identically for Single/Grid/Auto (the three
tabs that take real exposures), collapsing what would otherwise be
three near-identical copies of that group into one function called
three times.

Per-tab widgets are fully independent instances (`grid_start_spin` vs.
`auto_start_spin` vs. `replay_start_spin`, etc., all separately
named/prefixed) rather than one shared set of fields toggled by type,
as before. This means Grid/Auto/Replay each remember their own
last-used configuration independently, which is arguably a small UX
improvement, not just an implementation detail -- and it's what makes
"each tab shows only its own relevant fields" trivially true rather
than something that needs enable/disable logic to maintain.

New public API on the panel, replacing the old
type-radio-driven one: `get_sequence_type()`/`get_sequence_config()`
read whichever tab is currently active (`self.tabs.currentWidget()`),
returning `None`/`{}` for Single/Log/Help; `get_exposure_config()`
similarly reads Single/Grid/Auto's own fields, `{}` for Replay/Log/Help.
`show_best_focus(best_focus, best_fwhm)` dropped the old `can_move`
parameter (no more "Move to Best Focus" button to gate) and now also
sets `single_focus_spin`'s value to `round(best_focus)` -- the
mechanism behind "Single tab defaults to the best fit." New
`show_single_exposure_result(result)` replaces the old
`show_confirmation()`/`show_pending_exposure()` pair (only one flavor
of single-exposure result exists now). `set_running()` disables every
input widget across all four actionable tabs plus the method combo,
but deliberately never touches `self.tabs.setEnabled()` -- confirmed
directly with `test_tab_bar_stays_enabled_while_running` -- so the Log
tab stays reachable mid-run, per the earlier discussion.

Removed entirely: the sequence-type radios and
`_on_sequence_type_changed()` (superseded by tabs);
`moveToBestFocusRequested`/`move_to_best_focus_button`/
`_on_move_to_best_focus_clicked()` and its confirmation `QMessageBox`;
`addToSequenceRequested`/`add_to_sequence_button`/
`add_pending_to_sequence()`-adjacent view methods
(`show_pending_exposure()`/`clear_pending_exposure()`); the separate
"Result" group (`result_label`) -- its information now lives in the Log
tab's existing log line plus the Single tab's default value.

No `Controller`/`MainWindow` changes yet (as planned) -- `Controller`
still calls the old API, so `tests/gui/test_controller.py`,
`test_app_smoke.py`, and `test_phase3_smoke.py` currently fail (mostly
`AttributeError` on removed panel attributes/methods) when run on their
own; that's expected and is sub-phases 2-3's job. One thing worth
recording: running the *entire* `tests/gui/` directory together in this
intermediate state triggered a PySide6 segfault, not just clean test
failures -- but running `test_controller.py` alone (or any other single
file alone) reproduces only ordinary `AttributeError` failures, and a
standalone script constructing 30 `MainWindow`s back-to-back in a loop
raised nothing at all. This points to instability from combining many
already-broken (mid-refactor) test files in one pytest session rather
than a real bug in the new panel code, but it's not fully explained;
worth keeping an eye on once sub-phases 2-3 restore a fully consistent
test suite -- if the segfault recurs on a *consistent* suite, that
would be a real bug to chase.

**Testing:** `tests/gui/test_focus_control_panel.py` rewritten from
scratch (29 tests) covering: tab order; no-spinner-buttons and
integer-focus-value coverage across all four tabs;
`get_sequence_type()`/`get_sequence_config()`/`get_exposure_config()`
per tab, including the empty-dict cases; Start emitting the right
signal per tab (`startRequested` vs. `takeSingleExposureRequested`) and
being disabled on Log/Help; `set_running()` locking fields/method while
leaving the tab bar itself enabled, and correctly restoring (not
blindly re-enabling) the Start button for whichever tab is active once
stopped; `show_best_focus()` rounding into `single_focus_spin`, and
`reset()` explicitly *not* clearing that default; the Browse dialog;
and a light check that the Help tab actually mentions every other tab
by name. All 29 pass in isolation.

### Sub-phase 2 results: update `Controller`

Removed `move_to_best_focus()`, `add_pending_to_sequence()`, their
signal connections, `pending_result`, and the
`worker.sequence is self.sequence`/`_single_purpose`-style branching in
`_on_single_exposure_finished` that distinguishing "move to best focus"
from "standalone exposure" required -- only one single-exposure flow
exists now, so `_on_single_exposure_finished` unconditionally does what
the old "standalone" branch did: display the result, seed
`_standalone_sequence`'s bookkeeping (kept, for §5.6 reanalyze
targeting -- see the Phase 4 plan above), and report it via the new
`show_single_exposure_result()`. `_on_sequence_finished` now calls
`show_best_focus(best_focus, best_fwhm)` without the removed `can_move`
argument. `start_sequence()` resets `_standalone_sequence = None`
directly (no more `_clear_pending()` helper needed, since there's no
paired View state to clear alongside it anymore -- `reset()` on the
View side already doesn't touch anything Single-tab-related, per
sub-phase 1).

No deviations otherwise. **Testing:** `tests/gui/test_controller.py`
rewritten: `_configure_control_panel` became `_configure_replay_tab`
(switches to `replay_tab`, sets its own fields); every
`grid_radio.setChecked(True)`-style call became
`panel.tabs.setCurrentWidget(panel.grid_tab)`, etc. Dropped all
move-to-best-focus and add-pending-to-sequence tests; renamed/adapted
the standalone-exposure tests
(`test_source_selection_reanalyzes_a_standalone_exposure`,
`test_start_sequence_discards_standalone_sequence_state`) to check
`controller._standalone_sequence` instead of the removed
`pending_result`. Added an assertion to each of the three
sequence-completion tests confirming `single_focus_spin` picks up the
rounded fitted best focus -- the concrete behavior behind "Single tab
defaults to best focus," which no prior test exercised end-to-end
through the Controller (sub-phase 1's tests only checked the View
method in isolation). All 17 tests pass.

Ran `tests/gui/test_focus_control_panel.py` + `test_controller.py`
together (46 tests, all pass) to specifically check for the segfault
noted at the end of sub-phase 1 -- it did not recur. Ran the entire
`tests/gui/` directory: only `test_app_smoke.py` and
`test_phase3_smoke.py` fail now (plain `AttributeError`s on the removed
API, e.g. `datadir_edit`, `pending_result`), exactly the two files
sub-phase 3 will update -- no segfault there either. This is reasonably
strong evidence the earlier crash really was instability from a
mostly-broken suite (many files failing against a half-migrated API at
once) rather than a bug in the new code, though it was never fully
root-caused.

### Sub-phase 3 results: end-to-end tests + Phase 4 close-out

`tests/gui/test_app_smoke.py`: `_configure_control_panel` became
`_configure_replay_tab` (switches to `replay_tab`, sets its own
fields). `test_gui_controller_matches_direct_execute` lost its
assertion target -- `window.control_panel._best_focus` no longer exists,
since there's no separate "Result" state on the panel anymore, only the
rounded value baked into `single_focus_spin`. Comparing against the
rounded integer would have weakened the test's actual purpose (bit-exact
GUI/CLI fit equivalence), so instead it now connects directly to
`controller.worker.sequenceFinished` right after `start_sequence()`
returns (before pumping the event loop) to capture the *raw* fitted
float the same way `test_sequence_worker.py`'s `_run_and_collect` does,
and compares that.

`tests/gui/test_phase3_smoke.py` rewritten around the Single tab's new
role: standalone exposure with no sequence loaded -> mark a source (§5.6
reanalyze-in-place) -> start a Grid sequence with that method -> switch
to the Single tab and confirm its focus value already defaults to the
fitted best focus -> click Start there to actually move to it (no
confirmation dialog needed -- the previous sub-phase already removed
that entirely) -> mark a different source on the now-loaded sequence
(confirming reanalyze-target priority still favors a loaded sequence
over a standalone exposure) -> take one more standalone exposure ->
discard it by starting an Automated sequence. Dropped the old "Add to
Existing Sequence" step entirely, per the Phase 4 design decision.

No deviations. Full-suite and Qt-free-boundary checks: 112 tests pass;
26 pass (`tests/gui/` skipped as one unit) with PySide6 simulated
absent. Re-ran `tests/gui/` three times in a row (86 tests each) to
specifically watch for the sub-phase-1 segfault recurring now that the
whole suite is consistent again -- it did not, in any of the three
runs, reinforcing that it really was an artifact of the mid-refactor
state rather than a latent bug.

Also spot-checked the actual visual outcome the whole redesign was for:
building a real `MainWindow` offscreen, `FocusControlPanel`'s size
dropped from roughly 640x1050 (pre-redesign, with every group always
laid out) to 312x384 -- a concrete, measured confirmation that showing
only the active tab's fields (rather than disabling irrelevant ones in
place) actually solves the "minimize white space" request, not just in
principle.

### Phase 4 summary

All 3 sub-phases complete. `FocusControlPanel` is now a `QTabWidget`
(Single/Grid/Auto/Replay/Log/Help) with a shared photometry-method
selector and Start/Stop row below it. The Single tab absorbed "Move to
Best Focus" entirely (defaults to the last fit's rounded result; no
confirmation dialog); "Add to Existing Sequence" was dropped outright.
`Controller` and `scripts/focus.py`/`SequenceWorker` needed
correspondingly little and no change, respectively -- confirming the
sub-phase 1 prediction that this was purely a View/Controller
reorganization. `claude/GUI_DESIGN.md` still describes the pre-Phase-4
layout (radios, a separate Result group, the Move-to-Best-Focus
confirmation dialog, Add-to-Existing-Sequence) and needs reconciling
against this document during the end-of-GUI-development pass already
on record as deferred.

## Post-Phase-4 revisions

### Remove photometry method selector; report coordinates in the Log instead

Per request: dropped the "Photometry Method" group (`method_combo`,
`method_label`) and `methodChanged` entirely. Measurement always uses
the brightest detected source by default (`self.method = 'brightest'`
at `Controller` construction, unchanged); clicking a source and
pressing 'm' still overrides it to that source's coordinates via
`Controller.set_method()`, which is now just `self.method = method`
with no paired View call. Rather than a persistent "Current method: ..."
label, the coordinates *actually measured* -- `StepResult.centroid`,
which reflects whichever source was used either way -- are now folded
into the existing per-exposure log lines: `update_step()` (every
step/reanalysis) and the new-in-Phase-4 `show_single_exposure_result()`
both append `Source (x, y)` to their text. This is more accurate than
the old display ever was, too -- the old "Selected source (x, y)" text
showed the raw *clicked* point, not the actual detected centroid
`image_quality()` measured nearest to it.

`FocusControlPanel.get_selected_method()`/`set_method()` are gone; the
Help tab's "Method" paragraph became "Source selection," describing the
automatic-brightest default and the 'm'-key override without mentioning
a combo box that no longer exists.

**Testing:** removed `test_method_combo_and_signal` and
`test_set_method_updates_display_for_string_and_coordinates`
(`tests/gui/test_focus_control_panel.py`); dropped the `method_combo`
assertions from the renamed `test_set_running_locks_config_widgets`;
added coordinate assertions to `test_update_step_updates_label_and_log`
and `test_show_single_exposure_result_updates_status_and_log`. In
`tests/gui/test_controller.py`, `test_set_method_updates_state_and_display`
(renamed `test_set_method_updates_state`) and
`test_source_selection_updates_measurements_in_place` lost their
`method_label` assertions -- the latter already covers the coordinate
now living in the log indirectly via the shared `update_step()` path
tested elsewhere. All 110 tests pass; Qt-free boundary reconfirmed.

### Two-column tab layout to compress height further

New `_two_column_row(left_rows, right_rows)` helper lays a tab's fields
out as two side-by-side `QFormLayout`s instead of one stacked column:
exposure parameters (exptime/speed/binning) on the left, the focus
value(s) that define the sequence on the right, per request. Single's
right column is just its one focus field, as expected. `_add_exposure_rows`
(which mutated a `QFormLayout` directly) became `_exposure_field_rows()`,
returning `(label, widget)` pairs instead so it composes with the new
helper -- still one function shared by Single/Grid/Auto, unchanged in
spirit from Phase 4 sub-phase 1. Replay keeps its data-directory row
on its own line on top (unchanged position), with prefix/suffix/obsnum
(file-name construction) on the left and start/step/nstep (focus
values) on the right below it. Log and Help are untouched, per request.

Measured effect: a real `MainWindow`'s `FocusControlPanel` height
dropped from 384px (Phase 4's single-column tabs) to 264px with the
same window/tab content otherwise unchanged -- confirming the
two-column layout achieves real additional vertical compression, not
just a rearrangement.

**Testing:** no new tests -- purely a layout change with no new
behavior (every widget referenced by existing tests keeps the same
name and `panel.tabs`/tab-widget structure; `_config_widgets`,
`get_sequence_config()`, etc. are unaffected since they read widgets by
attribute, not by layout position). Reran the full suite to confirm
nothing broke: all 110 tests pass; Qt-free boundary reconfirmed.

### Move acquisition buttons into their own tabs

Per request: the shared Start/Stop row below the tab widget is gone.
Each actionable tab now owns its own button(s) at the bottom of its
layout (after the `addStretch(1)` already there, so they sit at the
tab's bottom edge regardless of window height -- confirmed directly by
checking button `y()` against tab `height()` once actually displayed):
Single has one "Acquire" (`single_acquire_button`); Grid and Auto each
have "Acquire"/"Interrupt" (`grid_acquire_button`/`grid_interrupt_button`,
`auto_acquire_button`/`auto_interrupt_button`); Replay has one "Load"
(`replay_load_button`). Log and Help have none.

This actually *simplified* the controller-facing logic, not just moved
buttons around: because a `QTabWidget` only lets the user interact with
the currently-visible page, each button unambiguously means "run this
tab's action" the moment it's clickable at all -- there's no longer a
need to track which tab is active to decide what a click *means* (only
`get_sequence_type()`/`get_sequence_config()` still consult
`self.tabs.currentWidget()`, to decide what a click that already
happened should build). This removed `_on_tab_changed()`,
`_update_start_button()`, `_on_start_clicked()`, and the
`currentChanged` connection entirely -- none of them are needed once
enablement no longer depends on which tab is showing, only on whether
something is running. `set_running()` now disables all four acquire
buttons and enables both interrupt buttons uniformly (either Interrupt
button works regardless of which tab is visible, since `stopRequested`
always targets whatever `Controller.worker` actually is, not "whichever
tab's sequence" -- there's only ever one). `Controller` needed **zero**
changes: it only ever depended on the `startRequested`/`stopRequested`/
`takeSingleExposureRequested` signals, never on the button widgets
themselves.

**Testing:** replaced the old shared-button tests in
`tests/gui/test_focus_control_panel.py` with
`test_each_tab_has_the_requested_buttons` (checks the exact button set
per tab, including that Replay's pre-existing "Browse…" button coexists
with "Load," and that Log/Help have none at all),
per-tab signal tests (`test_single_tab_acquire_emits_take_single_exposure_requested`,
`test_grid_tab_acquire_and_interrupt_signals`,
`test_auto_tab_acquire_and_interrupt_signals`,
`test_replay_tab_load_emits_start_requested`), and
`test_set_running_locks_config_widgets_and_acquire_buttons` (all four
acquire buttons disable, both interrupt buttons enable, and vice versa).
Updated the button references in `tests/gui/test_controller.py`
(`replay_load_button`) and `tests/gui/test_phase3_smoke.py` (each
`panel.start_button.click()` replaced with the specific tab's own
button). All 109 tests pass; Qt-free boundary reconfirmed.

### Reanalysis was iteratively replacing curve points instead of clearing first

**Problem:** pressing 'm' to mark a source and reanalyze a loaded
sequence fed each re-measured `StepResult` to
`FocusCurvePanel.update_result()`, which replaces one point in place
and redraws -- so as results streamed back from the worker one at a
time, the plot showed a mix of old and new measurements, with old
points disappearing one-by-one rather than the whole curve clearing
upfront.

**Fix:** `Controller.reanalyze()` now calls `self.window.curve_panel.reset()`
synchronously, before starting the worker, whenever the reanalyze
target is the actively-loaded `sequence` (not a standalone Single-tab
exposure, which never touches the curve panel at all). `_on_step_complete`'s
reanalyze branch now calls `curve_panel.add_result()` instead of
`update_result()`, since the panel starts empty and simply
accumulates fresh points as they arrive -- there's nothing left to
"replace." This made `FocusCurvePanel.update_result()` completely
unused (it was only ever called from this one call site), so it was
removed along with its two dedicated tests, rather than leaving dead
code behind.

**Testing:** added `test_reanalyze_clears_the_curve_panel_immediately_not_iteratively`,
which populates the curve, calls `controller.reanalyze()`, and asserts
`curve_panel._results == []` *before* pumping the event loop at all --
proving the clear happens synchronously rather than as a side effect of
processing results. Verified this fails without the fix (the curve
still holds all 5 old points at that point) and passes with it. All
108 tests pass; Qt-free boundary reconfirmed.

### A long Log message balloons the panel width -- again

While debugging an unrelated user-error (an incorrect obsnum on the
Replay tab), a long failure message surfaced the same symptom as the
"control panel width balloons" bug from Post-Phase-3 -- but a different
root cause this time. Confirmed directly: `status_label.minimumSizeHint()`
grew to ~1989px wide after `show_failure()` was given a realistic long
message (a "missing files" list). A plain `QLabel` without word wrap
reports its *entire single-line text's rendered width* as its minimum
size -- there's nothing a scroll-area width floor can do about the
widget's own natural minimum genuinely growing at runtime; that
mitigation was for a different failure mode (the viewport being
squeezed narrower than a *static* minimum).

**Fix:** `status_label.setWordWrap(True)` and `step_label.setWordWrap(True)`
(the latter proactively -- `update_step()` can also produce a fairly
long line with focus/FWHM/source coordinates all on one line). Wrapping
lets a long message grow the Log tab's *height* instead of the whole
panel's width; height is exactly what the surrounding `QScrollArea`
(from the original window-sizing fix) is already built to absorb.

**Testing:** added `test_status_and_step_labels_wrap_instead_of_widening`,
which feeds `show_failure()` a realistic long, space-separated message
(a single unbroken token wouldn't exercise the fix meaningfully, since
word wrap can only break between words) and asserts
`status_label.minimumSizeHint().width()` stays small. Verified this
fails without `setWordWrap(True)` and passes with it. All 109 tests
pass; Qt-free boundary reconfirmed.

## Phase 5: Polish

Per GUI_DESIGN.md §9, item 5. Per explicit decision, item 4 (evaluating
`pyqtgraph.ImageView` for `ImagePanel` responsiveness) is skipped
entirely -- matplotlib has been responsive enough in practice to not
warrant the swap. Phase 5 is three largely independent pieces of work:

### Sub-phase 1: CLI output-file writing

Implements `scripts/focus.py`'s standing `# TODO: Write the output file
if provided` in `main()`. When `--ofile` is given, write an ECSV table
of the completed sequence's per-exposure results. Column names must
match what the (still-unimplemented) `--refit` path already reads back
(`tbl['FOCUS']`, `tbl['SIGMA']`) for internal consistency, even though
"SIGMA" is really the FWHM-based image-quality metric
(`StepResult.fwhm`) rather than a literal Gaussian sigma -- that naming
predates this work and isn't being changed here. Additional columns
(obsnum/exposure filename, outlier flag) are added alongside for
traceability, since `--refit`'s stub only reads the two it needs and
extra columns are harmless.

Explicitly out of scope: implementing `--refit` itself (stays
`raise NotImplementedError('Not ready to refit.')`) -- the design doc's
Phase 5 bullet only names "output-file writing," not refitting, and
that's a separable, larger feature.

This is a Qt-free, CLI-only change with no GUI involvement -- testable
entirely with `tests/test_focus_cli.py`-style tests already in place.

### Sub-phase 2: GUI settings persistence

"Last-used exposure time, directories, etc." -- persist each tab's own
configuration (Single/Grid/Auto's exposure settings; Grid/Auto/Replay's
focus-sequence fields; Replay's datadir/prefix/suffix/obsnum) plus
`ImagePanel`'s stretch choice, using `QtCore.QSettings` (Qt's built-in
cross-platform persistent-settings mechanism -- an macOS
`~/Library/Preferences` plist, a Windows registry key, or a Linux ini
file, depending on platform, with no new dependency). Saved on window
close (`MainWindow.closeEvent`), restored on construction. Deliberately
*not* saved on every keystroke/change -- these are convenience defaults
for next time, not data that needs crash-safety, and per-field change
wiring across ~20 fields is a lot of code for that marginal benefit.

Uses flat, individually-typed key/value pairs (e.g. `'grid/start'`)
rather than one serialized blob, since `QSettings`'s native backends
don't reliably round-trip nested Python structures across platforms --
each panel gets `get_settings_state()`/`set_settings_state(state)`
methods operating on a flat `dict` of primitives, which `MainWindow`
reads/writes via `QSettings` key-by-key. `set_settings_state()` must
tolerate a state dict missing keys (a fresh install, or an older saved
state from before a field existed).

The Single tab's focus value itself is deliberately *not* persisted --
its default is meant to reflect the most recent fit
(`show_best_focus()`), and restoring a stale persisted value on startup
would fight with that.

**Testing:** pure `get_settings_state()`/`set_settings_state()`
round-trip tests on `FocusControlPanel`/`ImagePanel` (no real
`QSettings` needed); one `MainWindow`-level test using a real but
isolated `QSettings` instance (a distinct organization/application name
or an explicit `IniFormat` pointed at a temp file, cleaned up after the
test) confirming save-then-reconstruct actually restores values end to
end.

### Sub-phase 3: Richer stretch options

`ImagePanel.STRETCHES` currently only has ZScale and Min/Max (both
`LinearStretch`). Add more of `astropy.visualization`'s stretch classes
(already a dependency -- no new package) paired with `ZScaleInterval`,
e.g. `SqrtStretch`, `LogStretch`, `AsinhStretch`, giving more ways to
bring out faint structure or compress a large dynamic range.

**Testing:** confirm each new stretch entry actually produces a valid
`ImageNormalize` and renders without error against a synthetic frame
(mirroring how the existing ZScale/Min-Max entries are implicitly
covered today).

### Sub-phase 1 results: CLI output-file writing

Implemented in `scripts/focus.py`'s `main()`: when `--ofile` is given,
after `seq.execute()` completes, writes an ECSV table with columns
`EXPOSURE` (the exposure file path, as a string), `FOCUS`
(`seq.observed_focus`), `SIGMA` (`seq.img_quality` -- the FWHM-based
image-quality metric; kept as `SIGMA` to match what `--refit`'s
already-written reader expects, per the plan above), and `OUTLIER`.
`OUTLIER` is recomputed via `FocusPlot.is_outlier(seq.centroids[:i+1])`
for each row -- the same *cumulative-up-to-that-point* comparison
`step()` itself uses live, not a single pass over the complete, final
centroid list (which could flag different points as outliers than what
was actually seen/plotted during the run). Prints `Wrote focus data to
{path}` for user feedback, matching the existing `Best focus:`/
`Expected sigma:` print style.

No deviations from the plan. `--refit` remains
`raise NotImplementedError('Not ready to refit.')`, unchanged --
correctly out of scope per the sub-phase boundary.

**Testing:** added `test_cli_writes_output_file_when_requested` to
`tests/test_focus_cli.py`, which runs the CLI archive-mode path with
`--ofile`, reads the written ECSV back with `Table.read()`, and checks
row count, the recorded focus values, and that all four expected
columns exist. Added a one-line check to the existing
`test_cli_archive_mode_runs_end_to_end` confirming *no* "Wrote focus
data" message appears when `--ofile` isn't given. Also manually
verified end-to-end against real synthetic FITS files outside pytest,
confirming the file is valid, readable ECSV. All 110 tests pass;
Qt-free boundary reconfirmed (this sub-phase touches no GUI code at
all).

### Sub-phase 2 results: GUI settings persistence

`FocusControlPanel` gained `get_settings_state()`/`set_settings_state()`
plus a private `_settings_fields()` -- a single dict mapping each
persistable key (e.g. `'grid/start'`) to its widget, used by *both*
methods via `isinstance` dispatch (`QComboBox` -> `currentText()`/
`setCurrentText()`, `QLineEdit` -> `text()`/`setText()`, spin boxes ->
`value()`/`setValue()`). Keeping one field table rather than writing
the getter and setter as two separately-maintained lists is what
guarantees they can't drift apart into covering different fields.
`ImagePanel` gained the same two methods for its one persistable
preference, the stretch choice. The Single tab's focus value is
excluded from `FocusControlPanel`'s field table entirely, per the plan
-- it must keep tracking `show_best_focus()`'s most recent fit, not a
stale session-old default.

`MainWindow` calls `_load_settings()` at the end of `__init__` and
`_save_settings()` from an overridden `closeEvent()`, via
`QtCore.QSettings(_SETTINGS_ORG, _SETTINGS_APP)`. The load side reuses
each panel's *own current* `get_settings_state()` output as
`QSettings.value(key, default)`'s default argument -- Qt uses that
default's Python type to coerce whatever was actually stored, so no
separate key-to-type map is needed; a key that was never saved (fresh
install, or a field added after the user's last save) simply falls
back to the value the widget was already constructed with.

**A real bug caught before it reached anyone, not by a test:** the
first working version wrote to the actual per-machine `QSettings` store
-- an real macOS preferences plist at
`~/Library/Preferences/com.lickobservatory.NickelFocusGUI.plist` --
every time *any* test constructed (and every time one that called
`.close()` on) a `MainWindow`, since `_load_settings()`/`_save_settings()`
run unconditionally with no test-mode guard. Running the test suite
during this sub-phase actually created that file on this machine with
placeholder test values; it was deleted immediately upon discovery.
Fixed by adding an `autouse=True` fixture, `isolate_qsettings`, to
`tests/gui/conftest.py`: it monkeypatches `MainWindow._settings()` to
return a `QSettings` pointed at a per-test-function temp `.ini` file
instead, for every test under `tests/gui/` regardless of whether it
requests the fixture by name. Re-ran the full suite afterward and
confirmed via `ls`/`defaults read` that no real preferences file
reappears. This is the same category of lesson as the segfault noted in
Phase 4 sub-phase 1 -- a new mechanism (there: many widgets under a
half-migrated API; here: persistent state shared across the whole test
session) needs its own isolation, which doesn't exist until something
new actually needs it.

No other deviations. **Testing:** `get_settings_state()`/
`set_settings_state()` round-trip tests on both panels, plus a test
confirming `'single/focus'` is never a key in `FocusControlPanel`'s
state, and a test confirming an unknown/partial state dict only
touches the keys it names, leaving everything else alone (also
covered per-panel for `ImagePanel`'s single key with an unrecognized
stretch name). `tests/gui/test_main.py` gained an end-to-end test
building a window, changing settings across both panels, closing it,
building a *second* window against the same isolated settings file,
and confirming the values actually round-tripped -- plus a test
confirming the Single tab's focus value specifically does *not*
survive that round trip. Manually verified the same save/restore cycle
outside pytest against a real (but manually isolated) ini file, to be
sure this wasn't a self-fulfilling green test suite. All 117 tests
pass; Qt-free boundary reconfirmed; no real settings file created by a
full test run.

### Settings persistence made opt-in, via a new Options tab

Follow-up to sub-phase 2, before starting sub-phase 3: the user asked
whether persistence could be opt-in, since silently writing *any* file
-- even an inert preferences file -- for a user who never asked for it
is unwelcome. New **Options** tab (Single, Grid, Auto, Replay, Log,
**Options**, Help), holding one checkbox for now -- "Remember settings
between sessions" -- and explicitly scoped to grow as more app-wide (as
opposed to per-sequence) preferences become useful later.

The first working version still had a privacy gap: it saved the
checkbox's *own* checked state as an always-written flag, meaning every
user got a settings file the moment they closed the window, whether
they'd ever opted in or not -- exactly the silent-write problem the
request was about, just narrowed to one boolean instead of the whole
configuration. Caught by the user during review, before it shipped.

**Fix:** there is no separately-persisted opt-in flag at all. Whether
anything was ever saved *is* the opt-in signal:
`MainWindow._load_settings()` checks the checkbox if (and only if)
`'control_panel' in settings.childGroups()` -- i.e. a previous session
actually wrote something -- and `_load_settings()` never itself calls
anything that writes (`QSettings.value()`/`.childGroups()` are
read-only; merely constructing/reading a `QSettings` object doesn't
create its backing file). `_save_settings()` only calls `setValue()` at
all if the checkbox is checked; if it's unchecked, it removes the
`control_panel`/`image_panel` groups *only if they already exist*
(guarded the same way), so a session that never opts in never makes a
single mutating `QSettings` call and therefore never creates a file.
Verified manually end-to-end with a real isolated ini file: closing
without opting in leaves no file on disk at all; opting in and setting
a value creates one and it's readable next launch (with the checkbox
auto-checked); unchecking and closing again empties it back out, and
the next launch correctly detects "not opted in" and reverts to code
defaults.

**Testing:** rewrote the sub-phase 2 `MainWindow` tests around explicit
opt-in: `test_settings_default_unchecked_and_nothing_is_ever_written`
(closing without ever checking the box leaves no file at all, checked
via `os.path.exists` on the isolated settings file's actual path, not
just "values didn't come back"),
`test_settings_are_saved_on_close_and_restored_on_open_once_opted_in`
(now explicitly checks the box first), the existing
`test_settings_do_not_restore_a_stale_single_tab_focus_value` (now also
opts in first, since otherwise it wouldn't be testing anything), and a
new `test_unchecking_settings_erases_previously_saved_configuration`
(opt in, save, opt back out, confirm the saved group is actually gone
and a third launch reverts to defaults). Updated tab-order/button-set
tests in `tests/gui/test_focus_control_panel.py` for the new seventh
tab. All 119 tests pass; Qt-free boundary reconfirmed; confirmed via
`ls`/`defaults read` that no real preferences file is created by a full
test run.

### Sub-phase 3 results: richer stretch options

Added three entries to `ImagePanel.STRETCHES`: `Sqrt`, `Log`, `Asinh`,
each pairing `ZScaleInterval` (the same already-good limits `ZScale`
itself uses) with a different `astropy.visualization` stretch curve
instead of `LinearStretch` -- Sqrt/Asinh to bring out faint structure,
Log to compress a large dynamic range. All from `astropy.visualization`,
already a dependency, so no new package.

**Deviation caught by the full-suite run, not by a targeted test:**
the first attempt named these `'ZScale + Sqrt'`/`'ZScale + Log'`/
`'ZScale + Asinh'`, to make the paired interval explicit in the
dropdown. That made `stretch_combo`'s minimum width wide enough to push
`ImagePanel`'s (and thus the whole window's) minimum width to 812px --
just over the 800px offscreen test screen's available width, failing
`test_build_window_fits_within_the_screen` (a Post-Phase-3-era test,
run as part of the routine full-suite pass, not something written for
this sub-phase). This is the same class of problem as the "control
panel width balloons"/"long Log message" bugs from earlier -- a widget
whose minimum size depends on its content quietly grows past a
constraint elsewhere in the layout -- just triggered by wider combo-box
text instead of an unwrapped label this time. Fixed by shortening the
names to single words (`Sqrt`/`Log`/`Asinh`, dropping the redundant
"ZScale + " prefix) rather than doing anything more elaborate to the
layout -- simplest fix, and arguably better anyway: none of the other
combo entries spell out their interval either (`Min/Max` doesn't say
"Min/Max + Linear").

**Testing:** added `test_every_stretch_option_renders_without_error`,
which iterates every key in `STRETCHES`, selects it, calls
`show_result()`, and confirms an image actually got drawn (`ax.get_images()`
non-empty) -- covering the three new entries generically, and any
future ones added the same way, without hardcoding their names. All
120 tests pass; Qt-free boundary reconfirmed.

## Phase 5 summary

All 3 sub-phases complete, plus the opt-in settings revision that
followed sub-phase 2. `scripts/focus.py` gained `--ofile` output
writing; the GUI gained an Options tab whose one setting (for now)
governs opt-in persistence of every other tab's configuration and the
image stretch choice; `ImagePanel` gained three more stretch curves.
`claude/GUI_DESIGN.md` remains unreconciled against the as-built
GUI (tabs, buttons-in-tabs, no method selector, the Options tab, opt-in
persistence) -- still deferred to the end-of-GUI-development pass
already on record.

### Next

Awaiting further direction -- Phase 5 (and the currently-known scope of
GUI_DESIGN.md's phased plan, since Phase 4/pyqtgraph was explicitly
skipped) is complete.

## Code-quality pass (pre-GUI_SUMMARY.md)

Per request, ahead of writing a `GUI_SUMMARY.md` (replacing the earlier
plan to reconcile `GUI_DESIGN.md` itself, which is being left as-is): a
pass over every file in `gui/` for conciseness/modularity, complete
docstring coverage (including trivial/private helpers), and inline
comments explaining Qt mechanics at points a first-time Qt user would
need them.

**Docstrings:** added one to every function/class that lacked one,
confirmed after the fact with a small AST-based script checking every
`gui/` file for `FunctionDef`/`ClassDef` nodes with no docstring (empty
output = none left) -- more reliable than eyeballing across nine files.
`__init__` methods are excluded throughout, per direction: a class's own
docstring already covers its constructor. Multi-line docstrings follow
the existing convention of opening/closing `"""` on their own line;
single-line ones stay inline as long as they fit the file's existing
~98-character line length (checked the same way, across all of `gui/`).

**Conciseness/modularity** -- three concrete, low-risk extractions,
chosen for being exact, unambiguous duplication (not speculative
"might be reused later" abstraction):
- `ImagePanel`: `_zoom_at`/`_show_result_preserving_view` had an
  identical 3-line "clamp and set both scrollbars without triggering
  `_on_scroll`" block. Extracted to `_set_scroll_position(x0, y0)`.
- `FocusControlPanel`: every `_build_*_tab` ended with the same
  `QWidget()` + `setLayout()` + `return` triplet -- extracted to a
  module-level `_tab_widget(layout)`. The Grid/Auto tabs' identical
  Acquire/Interrupt button pair (construct, wire to
  `startRequested`/`stopRequested`, disable Interrupt) became
  `_build_acquire_interrupt_row()`. A smaller `_button_row(*buttons)`
  replaced the repeated "new QHBoxLayout, addWidget each button" used
  by all four actionable tabs. Deliberately did *not* merge the Grid
  and Auto tab-builders further (into one parameterized method) -- the
  residual duplication after these extractions is small, and a shared
  builder would need either `setattr`-based dynamic attribute names or
  a multi-item tuple/namespace return, either of which would make the
  two tabs *harder* to read top-to-bottom than the current
  near-identical-but-plainly-written pair, for a marginal line-count
  win. "Three similar lines beat a premature abstraction" applies just
  as much at nine similar lines.
- `Controller._on_step_complete`: `self.worker is not None and
  self.worker.mode == 'reanalyze'` was computed twice (once inline,
  once via the `reanalyzing_standalone` tuple); factored into one
  `is_reanalyze` local.

**Inline comments for Qt mechanics**, placed once at each concept's
first real occurrence rather than repeated everywhere: signal/slot
`connect()` (`Controller.__init__`), a signal's `.emit` used directly
as another signal's slot (`FocusControlPanel._build_acquire_interrupt_row`),
`QThread.run()` as the actual background-thread entry point plus why
emitting a signal from it is still safe (`SequenceWorker.run`),
overriding a Qt *event handler* (`closeEvent`) as distinct from a
signal connection, `QSplitter` vs. plain layouts and nesting them
(`MainWindow.__init__`), `QScrollArea`'s one-child-widget model,
`QMainWindow.setCentralWidget`, `QTabWidget` page-switching semantics
and stretch factors (`QHBoxLayout.addWidget(..., 1)`), `QFormLayout`'s
label/field pairing, matplotlib's own `mpl_connect` callback registry
as something distinct from Qt signals, and `QLabel.setTextFormat`
turning on HTML interpretation.

No behavior changes anywhere in this pass -- confirmed by running the
full suite (unchanged at 120 tests) and the Qt-free boundary check
(unchanged at 26 pass) after every file, not just once at the end.
