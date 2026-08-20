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

### Next

Phase 3, sub-phase 6: the standalone single-exposure workflow (§5.5) --
a GUI action to take one exposure at a user-chosen focus value
(reusing `take_single_exposure()` and the `mode='single'` worker path
built in this sub-phase), plus a choice to add it to the currently
displayed sequence or start a new one.
