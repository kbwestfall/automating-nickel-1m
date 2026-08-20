# Nickel Focus/Pointing GUI — Design Document

**Status:** Draft — under discussion on the `gui` branch.

## 1. Purpose

Replace the current command-line workflow in `scripts/focus.py` (and
eventually the pointing workflow in `scripts/move_to_target.py`) with an
interactive Qt GUI that lets an observer watch a focus sequence run,
inspect each exposure as it arrives, and see the FWHM-vs-focus curve build
up in real time — without needing to read terminal output or tab over to a
separate matplotlib window.

This document is meant to be edited as decisions get made. Sections marked
**OPEN QUESTION** need an answer before the corresponding part of the
design is final.

## 2. Relationship to existing code

The non-GUI work in `scripts/` already contains most of the domain logic
this GUI needs to drive:

- `scripts/focus.py`: `Focus`, `Exposure`, `ExposureConfig`, `ExposurePath`
  (KTL hardware wrappers), and `FocusSequence` and its subclasses
  (`GridFocusSequence`, `AutomatedFocusSequence`, `ArchiveFocusSequence`)
  — the step-by-step logic of "change focus, expose, measure."
- `scripts/photometry.py`: `image_quality()` — turns a FITS file into a
  FWHM measurement, a cutout stamp, and a centroid.
- `scripts/quadratic.py`: quadratic fit + vertex (best focus).
- The recently-added fallback behavior means all of this can already run
  without a live `ktl` connection (see `ArchiveFocusSequence` and
  `scripts/focus.py`'s `--obsnum`/`--datadir`/`--prefix`/`--suffix`
  handling) — useful both for GUI development away from the telescope and
  for reprocessing archived sequences from within the GUI itself.

**The one piece of `scripts/focus.py` that does not carry over as-is** is
`FocusPlot` and the way `FocusSequence.execute()` drives it: `execute()`
currently owns a blocking `while` loop that calls straight into a
matplotlib-specific plotting object. A Qt GUI needs to:

1. Run that loop somewhere other than the GUI's main thread (see
   §4.3), and
2. Report each step's results (focus value, FWHM, stamp, frame, centroid,
   outlier flag) to *whatever* is currently observing it, rather than to a
   hardcoded matplotlib object.

So the `FocusSequence` classes are the right foundation for the Model
(§4.1), but `execute()` should be refactored to decouple "step the
sequence and measure the result" from "draw the result somewhere."

**Decision:** turn the current `while` loop into a generator — `step()`
yields one `StepResult` at a time — rather than a callback/observer list.
`execute()` becomes a thin wrapper that consumes its own generator (so the
CLI's behavior is unchanged); the GUI's worker thread (§4.3) drives the
same generator directly, one `next()` per loop iteration, checking for a
stop request in between. `FocusPlot` (CLI) and the new Qt views (GUI)
become two interchangeable consumers of the same sequence-stepping code
instead of the GUI forking its own copy of the hardware/photometry logic.

## 3. Technology choices

**Decision — Qt binding: PySide6.** LGPL-licensed, free to use and
redistribute, and it's the binding Qt itself maintains. If the telescope
control-room software stack turns out to be pinned to Qt 5 for some
reason, this can be revisited (PySide2 is the closest fallback).

**Decision — FITS image display widget: staged.**

- **Phase 1:** embed matplotlib via `FigureCanvasQTAgg`, reusing the
  `ImageNormalize`/`ZScaleInterval`/`LinearStretch` display code already
  written for `FocusPlot`. Minimizes new code and gets the rest of the GUI
  (control panel, curve panel, worker-thread plumbing) up and running
  first.
- **Phase 2:** evaluate replacing it with `pyqtgraph.ImageView`, which has
  native pan/zoom/histogram-stretch built in. The primary criterion for
  switching is UI **responsiveness** of the image view — matplotlib
  redraws are comparatively slow, and `pyqtgraph` is built for
  fast/frequent image updates. `ginga` and a hand-rolled `QGraphicsView`
  remain fallback options if `pyqtgraph` doesn't meet that bar, but
  aren't the current plan.

Because `ImagePanel` is a self-contained View component (§4.2) that only
consumes `StepResult.frame`/`stamp` data from the Controller, this swap
should not require touching the Model or Controller — worth confirming
that stays true as `ImagePanel` gets built out, so the phase 2 migration
stays cheap.

## 4. Architecture (Model-View-Controller)

### 4.1 Model

Owns telescope/camera state and sequence logic; **no Qt imports**, so it
stays usable from the existing CLI script and stays unit-testable without
a display.

- Wraps (or *is*) the existing `Focus`, `Exposure`, and `FocusSequence`
  subclasses from `scripts/focus.py`.
- Emits step results as plain data (e.g., a small dataclass:
  `StepResult(obsnum, focus_value, fwhm, frame, stamp, centroid,
  is_outlier)`) via the `step()` generator (§2).
- Should work in two modes, mirroring what `scripts/focus.py` already
  supports:
  - **Live mode** — real `ktl` connection, real exposures.
  - **Archive/replay mode** — reprocess a directory of existing FITS files
    (this is `ArchiveFocusSequence` today), useful both for observers
    reviewing a past sequence and for GUI development without telescope
    access.
- **Single exposure** (§5.5): a Model operation that takes one exposure at
  a caller-specified focus value and runs photometry on it, returning a
  `StepResult`, without going through a `FocusSequence`'s `step()`
  generator or its stepping/stopping bookkeeping. Built from the same
  primitives a sequence step already uses (`Focus.set_to()`,
  `Exposure.expose()`, `photometry.image_quality()`), just invoked once.
- **Reanalyze** (§5.6): given one or more exposures already on disk
  (live-acquired or archived) — a full sequence, or just a single
  standalone exposure (§5.5) with no sequence loaded — rerun
  `image_quality()` on each file with a new `method` (e.g., a
  newly-selected source's coordinates) and replace the stored
  FWHM/stamp/centroid data — no new exposures are taken. For a sequence,
  this is the same "replay existing files through photometry" pattern
  `ArchiveFocusSequence` already implements for stepping, so
  `reanalyze()` should be built by reusing that path (e.g., driving
  `step()` again with `self.exposures`/`self.observed_focus` as the
  archive source and the new `method`) rather than duplicating it; for a
  single standalone exposure it's just one more `image_quality()` call.
  Either way, `reanalyze()` updates data in place (same obsnum(s),
  refreshed FWHM/stamp/centroid) rather than producing a new object, so
  Controller/View references to "the current sequence" or "the current
  exposure" stay valid across a reanalysis.

### 4.2 View

Passive Qt widgets; contain no telescope/photometry logic, only display
state and emit Qt signals for user actions.

- `ImagePanel` — the FITS display widget (§3, §5.2), showing the most
  recent frame by default, or a user-selected one via the exposure
  drop-down (§5.2) — live or archived alike.
- `FocusControlPanel` — sequence configuration and controls (§5.4).
- `FocusCurvePanel` — the FWHM-vs-focus plot (§5.3), likely a
  `FigureCanvasQTAgg` embedding matplotlib so the quadratic-fit drawing
  code in `FocusPlot.update_curve()` can be reused directly.
- `MainWindow` — lays the three panels out per §5.1 and hosts the menu/
  toolbar (open archive directory, choose live vs. archive mode, etc.).

### 4.3 Controller

Mediates between Model and View; this is also where **threading** lives.

Every hardware call in `scripts/focus.py` is blocking (`waitFor(...,
timeout=...)`, `time.sleep`-equivalent polling), and a real exposure takes
several seconds to tens of seconds. Running `FocusSequence.execute()`
directly on the GUI's main/event-loop thread would freeze the entire UI —
no cancel button, no repaint, nothing — for the full duration of the
sequence. The Controller therefore needs to run sequence execution on a
worker thread (`QThread`, or a `QRunnable` on a `QThreadPool`) and marshal
results back to the GUI thread via queued Qt signals:

```
SequenceWorker(QThread)
    stepComplete = Signal(object)      # a StepResult
    sequenceFinished = Signal(float, float)   # best_focus, best_fwhm
    sequenceFailed = Signal(str)
    stopRequested: bool                # checked between steps for a clean Stop
```

The `step()` generator from §2 is what makes this straightforward: the
worker thread's `run()` loop just iterates it directly —

```python
for result in sequence.step():
    self.stepComplete.emit(result)
    if self.stopRequested:
        break
try:
    self.sequenceFinished.emit(*sequence.fit_best_focus(...))
except ValueError as e:
    # e.g., fewer than 3 focus values collected (early Stop, or a very
    # short sequence) — fit_best_focus() raises rather than returning
    # a fit, per the existing CLI code.
    self.sequenceFailed.emit(f'Stopped early — not enough points for a focus fit: {e}')
```

— so the same generator backs both the CLI's `main()` loop and the GUI's
worker thread. **Decision:** if `fit_best_focus()` can't produce a fit
(fewer than 3 points — most likely from an early Stop), the worker emits
`sequenceFailed` with a clear message rather than letting the exception
propagate; `FocusCurvePanel` still shows whatever raw (un-fit) scatter
points were collected, just without a curve or vertex.

**Decision:** "Stop" only aborts *between* steps — it never tries to
interrupt an in-progress exposure. Concretely, `stopRequested` is checked
after a `StepResult` is yielded (as in the loop above) and never during
`Exposure.expose()` itself, so a Stop click lets the current exposure and
measurement finish normally before the sequence halts.

Because that means there can be a real delay (up to one full exposure)
between clicking Stop and the sequence actually stopping, `FocusControlPanel`
(§5.4) needs to make that wait visible rather than leaving the user
wondering whether the click registered: on Stop, immediately disable the
button and switch the status area to something like "Stopping — waiting
for current exposure to finish...", then clear it once
`sequenceFinished`/`sequenceFailed` actually arrives.

**Hardware exclusivity.** `Focus` and `Exposure` represent single, shared
pieces of hardware, so the Controller only ever allows one operation
against them at a time. Concretely, "Start Sequence" and "Take Single
Exposure" (§5.5) are mutually exclusive and disabled while the other is
running, and interactive source selection (§5.6) — which doesn't touch the
hardware at all, only re-runs photometry on files already on disk — is
correspondingly only *enabled* while neither is running (i.e., before a
sequence/single exposure starts, using an already-obtained image, or after
one finishes).

## 5. Layout and functional requirements

### 5.1 Overall layout

```
┌─────────────────────────────┬───────────────────────────────┐
│                             │   FocusCurvePanel               │
│                             │   (FWHM vs. focus, fitted        │
│                             │    quadratic + vertex)           │
│        ImagePanel           ├───────────────────────────────┤
│   (pan / zoom / recenter /  │                                 │
│      stretch controls)      │   FocusControlPanel              │
│                             │   (sequence setup, start/stop,   │
│                             │    live status/log)              │
└─────────────────────────────┴───────────────────────────────┘
```

### 5.2 `ImagePanel` (left)

Baseline functionality, per your description:

- **Pan** via scroll bars.
- **Zoom** (mouse wheel and/or +/- controls; exact interaction TBD with
  the chosen widget, §3).
- **Recenter** (e.g., a button to re-fit the image to the view, and
  possibly a "center on measured source" action tied to the current
  `StepResult.centroid`).
- **Stretch** — at minimum, the ZScale/linear stretch already used by
  `FocusPlot`; likely worth exposing a small set of options (linear, log,
  asinh — all available via `astropy.visualization`) plus manual
  vmin/vmax.

**Decision:** yes to both.

- The outlier-highlight box from `FocusPlot.update_frame()` (red/yellow
  box around the measured source, "Outlier centroid" label when flagged)
  carries over to `ImagePanel`.
- The user can step back through prior exposures via a single **drop-down**
  (not a filmstrip/slider, and not two separate menus) listing every
  exposure taken so far, with each entry showing *both* its observation
  number and its focus value (e.g. "d2165 — Focus 345.0"), so one
  selection identifies the exposure by either value at once.
- The drop-down is **always kept sorted by focus value**, not by
  acquisition order. Since `AutomatedFocusSequence`'s adaptive search
  doesn't acquire focus values monotonically, a new `StepResult` is
  inserted at the position its focus value belongs at, which may be in the
  middle of the list rather than at the end — the visible order can shift
  as a sequence progresses.

This applies during a **live** sequence too, not just in archive/review
mode: the drop-down is browsable while a sequence is running. The
"auto-follow latest" behavior only kicks in if the currently-displayed
exposure *is* the most recently-acquired one (regardless of where it sorts
in the list) — in that case, each new `StepResult` is inserted into the
drop-down and immediately displayed. If the user has manually selected a
different entry, new results are still inserted into the drop-down (so
they don't disappear) but do **not** yank the view away from what the user
is currently looking at.

### 5.3 `FocusCurvePanel` (top right)

- Scatter of FWHM vs. focus value for the sequence so far.
- Best-fit quadratic curve and vertex marker, once ≥3 points exist (same
  logic as `FocusPlot.update_curve()`).
- **Decision:** outlier points (per `FocusPlot.is_outlier()`) are visually
  distinguished from normal points using a different color and/or marker
  symbol (exact styling — e.g. an "x" marker vs. a filled circle, and/or
  the same yellow used in `ImagePanel`'s outlier box — to be picked during
  implementation), consistent with how they're already flagged in the
  image panel and stamp grid.

### 5.4 `FocusControlPanel` (bottom right)

Needs to cover what `scripts/focus.py`'s CLI arguments currently do:

- **Sequence type**: grid (fixed start/step/end or nstep), automated
  (adaptive curve-following with a max-steps cap), or archive/replay
  (existing directory + obsnum + prefix/suffix).
- **Exposure settings**: exposure time, speed, binning (mirrors
  `ExposureConfig`).
- **Photometry method**: brightest / weighted / coordinate-based (mirrors
  `image_quality()`'s `method` argument).
- **Start / Stop** controls, disabled/enabled appropriately while a
  sequence is running.
- **Live feedback**: current step number, current focus value, current
  FWHM, and a running log/status area for warnings or errors (e.g., "no
  source detected," "controller not ready").
- **Result display**: once the sequence completes, the best-fit focus
  value and expected FWHM (today just printed to the terminal).

**Decision:** yes — "Move to best focus" is a separate, explicit action
(not an automatic side effect of a sequence finishing), matching
`FocusSequence.execute(goto=...)`, and it requires a confirmation dialog
before actually commanding `Focus.set_to()`. This is deliberately the
*only* action in the GUI that gets a confirmation dialog — ordinary
sequence steps and "Take Single Exposure" (§5.5) also call `Focus.set_to()`
but don't get one, since those focus values were directly typed/configured
by the user, which already counts as confirmation. "Move to best focus"
is different on two counts: it's a software-computed value the user
hasn't explicitly typed (so `Focus.set_to()`'s existing safety comment,
"Do NOT enable movement via this script!", applies most directly here),
and clicking it is also the observer's explicit signal that they're
satisfied with the sequence's result and ready to act on it — the GUI
should never *assume* success and quietly move the telescope on the
sequence's behalf, since that risks leaving the observer unsure what the
telescope's actual current focus is. The action should be disabled until a
sequence has produced a best-focus result, and the confirmation should
state the specific focus value about to be commanded.

### 5.5 Single-exposure workflow

- A new action in `FocusControlPanel`, e.g. "Take Single Exposure": the
  user specifies a focus value (and uses whatever exposure settings/method
  are currently configured), and the Model's single-exposure operation
  (§4.1) takes one exposure and measures it, without any sequence
  bookkeeping.
- Subject to the hardware-exclusivity rule above: unavailable while a
  sequence is running, and starting it disables "Start Sequence" until it
  completes.
- Typical use: confirm the telescope landed on the intended field/star,
  and — via interactive source selection (§5.6) — mark exactly which star
  should be used for FWHM measurement before a real sequence ever starts.
  For example: an observer slews to a focus field and takes a single
  exposure to confirm the star is visible, loads that frame into the GUI,
  and marks the star they want to use (instead of relying on, say,
  "brightest") — then starts a focus sequence that measures the FWHM of
  whichever source is nearest those marked coordinates at every focus
  value.
- The result is displayed in `ImagePanel` like any other `StepResult`, and
  is eligible for interactive source selection (§5.6) since it counts as
  "an already-obtained image" before a sequence begins.
- Once satisfied, the user explicitly chooses one of:
  - **Add to existing sequence** — appends this exposure's `StepResult`,
    as measured, to the currently-loaded sequence's data; it's inserted
    into the `ImagePanel` drop-down (§5.2) like any other exposure, and
    `FocusCurvePanel` refits/redraws.
  - **Start new sequence** — begins a new
    `GridFocusSequence`/`AutomatedFocusSequence` using the currently
    configured settings (start/step/end, exposure settings, and whatever
    `method` — including a just-marked source — is active). The single
    exposure is **not** automatically counted as one of the new
    sequence's own data points; its role was to confirm the field and/or
    mark the source. If the observer wants it included in the new
    sequence's curve too, they can add it afterward via "Add to existing
    sequence" once that sequence exists.

### 5.6 Interactive source selection

- `ImagePanel` supports **clicking near a source** in the displayed image
  to select it for FWHM measurement. This reuses `photometry.image_quality()`'s
  existing coordinate-based `method` argument — a `(x, y)` tuple already
  selects whichever detected source is nearest those coordinates — so the
  GUI's job is just translating a display click into image-data
  coordinates and setting `method` to that tuple; no new photometry logic
  is needed.
- **Gating:** per the hardware-exclusivity rule above, this is only
  enabled before a sequence/single exposure starts (using an
  already-obtained image) or after one finishes — never while one is
  actively running.
- **No gesture conflict:** since panning is via scroll bars and zoom is via
  the mouse wheel (§5.2's baseline interactions), a plain click on the
  image isn't used for anything else, so it's free to mean "select this
  source" without needing a separate "selection mode" toggle.
- **Reanalysis:** selecting a source automatically reruns the Model's
  `reanalyze()` operation (§4.1) over every exposure currently on display —
  re-measuring each already-acquired FITS file with the newly selected
  coordinates, replacing its stored FWHM/stamp/centroid, and redrawing
  `ImagePanel`, `FocusCurvePanel`, and the stamp grid as applicable. No new
  exposures are taken. This applies whether there's a loaded sequence
  (reanalyze every exposure in it) *or* just a standalone single exposure
  (§5.5) with no sequence loaded yet — in the latter case, `reanalyze()`
  runs against that one exposure, so its own displayed box/stamp/FWHM
  update immediately and the observer gets instant visual confirmation
  that they marked the right star, rather than only affecting future
  exposures.
- Because `reanalyze()` only touches files already on disk (no `ktl`
  calls), it doesn't need the same worker-thread treatment as a live
  sequence for correctness — but running it on a background thread is
  still worth doing if there are enough frames that re-running photometry
  on all of them would visibly stall the UI.

## 6. Pointing workflow — scope for v1

The branch is named for "the focus and pointing process," but the layout
sketch above only covers focus. `scripts/move_to_target.py` (find nearest
pointing star, slew there) is a separate, much simpler workflow today.

**Decision:** v1 is focus-only. Pointing (`scripts/move_to_target.py`)
will be designed and added in a later pass, once the focus GUI's
architecture (Model/View/Controller split, worker-thread pattern, image
panel) is proven out — at that point it should mostly be a matter of
adding a new Model wrapper around `move_to_target.py` and a small control
panel, likely reusing `ImagePanel` to confirm the telescope landed on the
right star.

## 7. Proposed module layout

Sketch, to be refined once §3's technology choices are settled:

```
gui/
    __init__.py
    main.py                 # entry point, builds MainWindow and starts the Qt event loop
    model/
        __init__.py
        sequence_worker.py   # QThread wrapper driving FocusSequence
        step_result.py        # StepResult dataclass
    views/
        __init__.py
        main_window.py
        image_panel.py
        focus_curve_panel.py
        focus_control_panel.py
    controller.py            # wires View signals -> Model calls, Model results -> View slots
```

(This would sit alongside `scripts/`, importing from it rather than
duplicating any hardware or photometry logic.)

## 8. Phased plan

1. **Model refactor**: turn `FocusSequence.execute()`'s loop into a
   `step()` generator (§2), with `execute()` becoming a thin wrapper so
   the CLI script's behavior is unchanged. Do this before or alongside the
   GUI skeleton, since the Controller design depends on it.
2. **Skeleton + archive mode**: `MainWindow` with all three panels wired
   up against `ArchiveFocusSequence` only (no `ktl` dependency), using
   embedded-matplotlib for `ImagePanel` (§3), so the GUI is fully testable
   away from the telescope using data already on hand. Interactive source
   selection and `reanalyze()` (§5.6) fit naturally here too, since they
   only need files already on disk, no live hardware.
3. **Live mode**: wire the Controller's worker thread to real
   `GridFocusSequence`/`AutomatedFocusSequence`, including Start/Stop, the
   safety-relevant confirmations from §5.4, and the single-exposure
   workflow (§5.5), which needs real hardware access.
4. **Image panel responsiveness pass**: evaluate swapping `ImagePanel`'s
   matplotlib backend for `pyqtgraph.ImageView` (§3), driven by measured
   responsiveness rather than switching preemptively.
5. **Polish**: settings persistence (last-used exposure time, directories,
   etc.), output-file writing (this is also an open TODO in the CLI
   script), richer stretch options.
6. **Pointing workflow** (§6), once the above is stable.

## 9. Summary of decisions

No open questions remain from this design pass. Decisions made so far:

- **Qt binding:** PySide6 (§3).
- **Image widget:** embedded matplotlib first, `pyqtgraph.ImageView`
  evaluated later against measured responsiveness (§3).
- **Model refactor:** `FocusSequence.execute()`'s loop becomes a `step()`
  generator, driving both the CLI and the GUI worker thread (§2, §4.3).
- **Stop behavior:** only between steps, never mid-exposure; the control
  panel must show a "stopping — waiting for current exposure" state while
  it waits (§4.3).
- **Image panel:** shows the outlier-highlight box, and lets the user
  browse prior exposures — live or archived — via a single drop-down
  keyed on both obsnum and focus value, always kept sorted by focus value
  (entries insert wherever they belong, not just at the end); browsing an
  older entry doesn't get interrupted by new results arriving, which only
  auto-display if the view is already tracking the most recently-acquired
  exposure (§5.2).
- **FWHM curve:** outlier points get a distinct color/marker (§5.3).
- **Single exposure:** a hardware-exclusive action independent of any
  sequence; its result can be added to an existing sequence's data or
  left standalone while a new sequence is started, but is never
  automatically folded into a new sequence's own data points (§5.5).
- **Interactive source selection:** clicking a source in `ImagePanel` sets
  the coordinate-based photometry `method` and automatically triggers
  `reanalyze()` in place — over a loaded sequence's existing files, or
  just a standalone single exposure if no sequence is loaded — so the
  observer gets immediate visual confirmation; only enabled outside of an
  active sequence/single exposure, since panning/zoom use scroll
  bars/wheel and never conflict with a plain selection click (§5.6).
- **Early stop with too few points:** if Stop leaves fewer than 3 points
  for `fit_best_focus()`, the worker emits `sequenceFailed` with a clear
  message instead of raising; the curve panel still shows the raw points
  collected (§4.3).
- **Pointing scope:** out of v1, revisited once the focus GUI is stable (§6).
- **Move to best focus:** the one action in the GUI that gets a
  confirmation dialog, stating the target focus value, disabled until a
  result exists — both because it's a computed (not user-typed) value, and
  because clicking it is the observer's explicit sign-off that the
  sequence succeeded (§5.4).

Anything not covered above should still be treated as undecided and raised
before being implemented.
