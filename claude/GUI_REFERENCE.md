# GUI Reference

Developer-facing reference for the `gui` package: what it is, how it's put
together, and how its pieces talk to each other. This document tracks the
implementation as it actually exists and should be kept up to date as the
GUI changes (e.g. the pointing workflow, the eventual package restructuring).

It is a companion to, not a replacement for, `claude/GUI_DESIGN.md` (the
original design spec, left as-is) and `claude/GUI_IMPLEMENTATION.md` (the
phase-by-phase build log).

## 1. Summary

`gui` is a PySide6 (Qt) desktop application wrapping `scripts/focus.py`'s
telescope focus-sequence automation for the Nickel 1-m telescope at Lick
Observatory. It lets an observer configure and run a focus sequence (a
fixed grid of focus values, or an adaptive search), inspect each exposure
and its measured FWHM as it comes in, watch the live focus curve, and
replay/reanalyze previously collected exposures from disk -- all without
using the `focus.py` CLI directly.

The GUI is strictly a View+Controller layer on top of the existing Model
code in `scripts/` (`focus.py`, `photometry.py`, `quadratic.py`). It adds
no new domain logic of its own; every hardware call, exposure, fit, and
sequence rule still lives in `scripts/focus.py`, unchanged and still
independently usable from the CLI with no Qt installed at all.

## 2. Directory structure and design

```
gui/
├── __init__.py            adds scripts/ to sys.path; no Qt import
├── qt.py                  the only module that imports PySide6
├── main.py                entry point (`python -m gui.main`)
├── controller.py          wires views to the model; owns app state
├── model/
│   ├── __init__.py
│   └── sequence_worker.py QThread that drives a FocusSequence
└── views/
    ├── __init__.py
    ├── main_window.py       top-level window, layout, settings persistence
    ├── image_panel.py       exposure display/pan/zoom/source selection
    ├── focus_curve_panel.py live FWHM-vs-focus plot
    └── focus_control_panel.py sequence configuration tabs, log, options
```

This is a fairly conventional Model-View-Controller split, with one Qt-
specific wrinkle:

- **Model** is `scripts/focus.py` and friends -- reused as-is, not
  reimplemented. `gui/model/sequence_worker.py` isn't domain logic; it's
  the adapter that runs the Model's blocking, synchronous methods
  (`FocusSequence.step`, `.reanalyze`, `.take_single_exposure`) on a
  background `QThread` so a real exposure (seconds to tens of seconds)
  doesn't freeze the GUI's event loop.
- **View** is `gui/views/` -- three panels (`ImagePanel`,
  `FocusCurvePanel`, `FocusControlPanel`) composed inside `MainWindow`.
  Views only display state and emit Qt signals for user actions
  (button clicks, a key press, a picked drop-down entry); they never
  import `focus` or talk to `SequenceWorker` directly.
- **Controller** (`gui/controller.py`) is the only module that imports
  both the views and the model. It owns the currently active
  `focus.FocusSequence` and `SequenceWorker` (if any), connects each
  view's signals to its own handler methods, and implements the
  "hardware exclusivity" rule: only one operation (a running sequence, a
  single exposure) may be active at a time.

`gui/qt.py` exists so exactly one module fails, with one clear message, if
PySide6 isn't installed -- every other `gui` module imports Qt classes from
here (`from gui.qt import QtWidgets`) rather than importing PySide6
directly. `scripts/` never imports `gui.qt` and must never need to: the
base `requirements.txt` is enough to run the CLI with no Qt involved,
`requirements-gui.txt` adds the one extra dependency the GUI needs.

Views hold no reference to each other or to the Controller. All
cross-panel communication (e.g. a new exposure updating both the image
display and the focus curve) is mediated by the Controller reacting to one
signal and calling methods on multiple views -- no view ever calls another
view's methods directly.

## 3. Key modules

| Module | Purpose |
| --- | --- |
| `gui/qt.py` | Single point of contact with PySide6; the only module allowed to `import PySide6` directly. |
| `gui/main.py` | Builds the `QApplication`, `MainWindow`, and `Controller`; starts the Qt event loop. |
| `gui/controller.py` | Mediates between views and model; owns the active sequence/worker and the hardware-exclusivity state machine. |
| `gui/model/sequence_worker.py` | `QThread` subclass that runs one `FocusSequence` operation (step/reanalyze/single exposure) off the GUI thread, reporting progress via signals. |
| `gui/views/main_window.py` | Top-level `QMainWindow`; lays out the three panels in splitters; persists/restores settings via `QSettings`. |
| `gui/views/image_panel.py` | Displays collected exposures one at a time: pan, zoom, stretch selection, and click-to-select-source ('m' key) for reanalysis. |
| `gui/views/focus_curve_panel.py` | Live scatter plot of FWHM vs. focus value, with the fitted quadratic and its vertex once enough points exist. |
| `gui/views/focus_control_panel.py` | Tabbed sequence configuration (Single/Grid/Auto/Replay), plus Log and Options tabs; emits the signals that request an action. |

### Lines of communication

```
                        ┌─────────────────────────┐
                        │        MainWindow        │
                        │  (layout + settings)     │
                        └───────────┬───────────────┘
                                    │ owns
              ┌─────────────────────┼─────────────────────┐
              ▼                     ▼                      ▼
      ImagePanel            FocusCurvePanel         FocusControlPanel
   (exposure display)       (FWHM vs focus)       (config tabs, log, options)
              │                     ▲                      │
    sourceSelected            add_result /            startRequested /
    (Signal)                  update_result           stopRequested /
              │                     ▲                takeSingleExposureRequested
              ▼                     │                      │
              └─────────────► Controller ◄──────────────────┘
                                    │
                          builds / starts / stops
                                    ▼
                            SequenceWorker (QThread)
                     stepComplete / sequenceFinished /
                     sequenceFailed / singleExposureFinished
                                    │
                                    ▼
                     scripts/focus.py (FocusSequence and
                     subclasses, unmodified Model code)
```

- **Views → Controller**: exclusively via Qt signals
  (`FocusControlPanel.startRequested`/`stopRequested`/
  `takeSingleExposureRequested`, `ImagePanel.sourceSelected`), connected
  once in `Controller.__init__`. Views never call Controller methods
  directly.
- **Controller → Views**: exclusively via direct method calls (e.g.
  `self.window.image_panel.add_result(result)`,
  `self.window.control_panel.show_failure(message)`) -- views expose a
  small "public API used by the Controller" for this, documented at the
  top of each view's method list.
- **Controller → SequenceWorker**: the Controller constructs a new
  `SequenceWorker` per run (`_start_worker`), connects its signals to its
  own `_on_*` handlers, and calls `.start()`. `SequenceWorker.request_stop()`
  is the only method called on a running worker.
  SequenceWorker → Controller: Qt signals only
  (`stepComplete`, `sequenceFinished`, `sequenceFailed`,
  `singleExposureFinished`, plus `QThread`'s built-in `finished`), safely
  delivered across the thread boundary by Qt's queued-connection
  machinery.
- **SequenceWorker → Model**: ordinary synchronous Python calls into
  `focus.FocusSequence` (`.step()`, `.reanalyze()`, `.take_single_exposure()`,
  `.fit_best_focus()`) -- the Model has no idea it's being driven from a
  background thread.

## 4. Usage

```
python -m gui.main
```

Requires PySide6 (`pip install -r requirements-gui.txt`, in addition to
the base `requirements.txt`). Configure a sequence on the Single, Grid,
Auto, or Replay tab and press Acquire/Load; watch progress on the image
panel, focus curve, and Log tab; press Interrupt to stop a running
Grid/Auto sequence between steps. See the in-app Help tab for a short
per-tab usage summary.

Installation and testing instructions will be added here once the repo is
restructured into an installable Python package.

## 5. Version history

This column will eventually hold git tags; until the project starts
tagging releases, it holds the short commit hash current as of each entry.

| Version | Date | Summary |
| --- | --- | --- |
| `c6058e8` | 2026-08-24 | Initial reference document, covering the GUI as of the code-quality/docstring pass following Phase 5. |
