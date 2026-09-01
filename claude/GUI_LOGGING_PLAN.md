# Wire `logging` into the GUI's Log tab (phased)

## Context

The GUI's Log tab (`FocusControlPanel`, built in `_build_log_tab()`) currently
shows a hand-built narration: `Controller` calls methods like `show_failure()`
and `update_step()`, which each format a string and push it into
`self.log_widget` (a `QPlainTextEdit`) via direct `appendPlainText()` calls.
None of this touches Python's `logging` module. Meanwhile, `nickel_focus`
already has a real logging system (`pkg/logger.py`'s `NickelFocusLogger`/
`get_logger()`), instantiated once as the package-level singleton
`nickel_focus.log`, but it's used in exactly one place (`slew.py`) — every
other module (`focus.py`, `photometry.py`, the CLI scripts) still uses raw
`print()`.

The goal (user's own two milestones for this branch):
1. Wire the existing `logging` singleton into the Log tab so log records
   actually appear there, safely across the GUI's worker threads.
2. Consolidate around it — retire the ad hoc string-passing in favor of real
   `log.*()` calls, and finish the `print()` → `logging` migration in the
   domain/CLI layers that was previously deferred (see `project_logging_migration`
   memory: this was blocked on deciding how to handle `capsys`-based CLI tests,
   which this plan resolves by switching them to `caplog`).

Key technical fact driving the design: `FocusWorker`/`SlewWorker` are `QThread`
subclasses that do the real work off the GUI thread and report back only via
Qt signals (`stepComplete`, `focusSequenceFailed`, etc.), which Qt auto-queues
onto the GUI thread before invoking connected slots. Once domain code
(`focus.py`, `photometry.py`) starts calling `log.*()` directly, those calls
can happen **on a worker thread**. A `logging.Handler` must not touch
`log_widget` directly from `emit()` — it must forward through a Qt signal,
mirroring the existing worker-signal pattern, so Qt does the thread marshaling
for us.

Also relevant: `photometry.py`'s `verbose=` diagnostic prints (lines 82-96,
122-137, 236, 247) are **never actually reachable** today — `focus.py`'s
`step()`/`reanalyze()`/`take_single_exposure()` call `image_quality()` without
ever passing `verbose=True` (grep confirms no caller sets it). Migrating these
to `log.debug()` isn't just a style change, it fixes dead diagnostics and
makes them available for free wherever the console/file logger's level is
already turned up (`scriptbase.py`'s existing `-v/--verbosity` flag already
wires `NickelFocusLogger.convert_verbosity_to_logging_level()` for the
CLI — nothing new needed there).

---

## Phase 1 — Thread-safe Qt log handler, wired in additively

Get real log records flowing into the existing Log tab without touching any
other behavior yet.

- **New file** `nickel_focus/gui/log_handler.py`: a `QtLogHandler` that is
  both a `QtCore.QObject` and a `logging.Handler`, with a
  `record_logged = QtCore.Signal(str)`. Its `emit(record)` formats the record
  (see formatter note below) and does `self.record_logged.emit(text)` —
  nothing else. This lives in `gui/`, not `pkg/logger.py`, because `pkg/logger.py`
  must stay import-light and PySide6-free (per `gui/qt.py`'s own docstring:
  Qt is an optional extra, only ever imported through that shim).
- **New formatter** in `pkg/logger.py`: add a `GuiFormatter(logging.Formatter)`
  next to the existing `StreamFormatter`/`FileFormatter`, e.g.
  `"%(levelname)8s | %(message)s"` — no ANSI color codes (they'd render as
  literal escape junk in a `QPlainTextEdit`), and terser than `FileFormatter`'s
  full `filename:funcName:lineno` since this panel is narrow and scrolls live.
- **New widget method** on `FocusControlPanel`: `append_log_line(text)`,
  a one-line wrapper around `self.log_widget.appendPlainText(text)`. Keeps
  `log_widget` encapsulated behind the panel's existing public-method
  convention (same pattern as `update_step`/`show_failure`/etc.).
- **Wire-up in `Controller.__init__`** (`gui/controller.py`), since it's
  already the layer that imports `from nickel_focus import focus, slew`
  (which triggers package init and creates `nickel_focus.log`) and already
  wires View ↔ Model:
  ```python
  from nickel_focus import log
  from nickel_focus.gui.log_handler import QtLogHandler
  ...
  # in __init__, after window/state setup:
  for h in [h for h in log.handlers if isinstance(h, QtLogHandler)]:
      log.removeHandler(h)  # see gotcha below
  self._qt_log_handler = QtLogHandler()
  self._qt_log_handler.setFormatter(GuiFormatter())
  self._qt_log_handler.setLevel(logging.INFO)
  self._qt_log_handler.record_logged.connect(window.control_panel.append_log_line)
  log.addHandler(self._qt_log_handler)
  ```
- **Gotcha to handle explicitly**: `nickel_focus.log` is a **module-level
  singleton** (created once at package import in `__init__.py:15`), not
  per-`Controller`. If `Controller.__init__` unconditionally
  `addHandler()`s, every new `Controller` (e.g. one per GUI test) leaves a
  stale handler pointing at a since-destroyed `log_widget`, accumulating
  duplicate/dangling handlers across the test session or across any future
  window-rebuild. Strip prior `QtLogHandler` instances first (shown above),
  mirroring how `NickelFocusLogger.init()` itself clears old handlers on
  reinit (`pkg/logger.py:178-183`).

This phase is purely additive: `show_failure`/`update_step`/etc. keep working
exactly as before, and log output now also shows up whenever anything calls
`log.*()` (today, just `slew.py:95`).

**Verification**: `nickel_focus_gui --test` and trigger a slew (exercises
`slew.py`'s `log.info` call) — confirm the line appears in the Log tab.

**New unit tests**:
- `tests/test_logger.py` (new, or add to an existing `pkg/logger.py` test
  module if one exists): `test_gui_formatter_plain_output` — feed
  `GuiFormatter().format(record)` a hand-built `LogRecord` per level and
  assert the output matches `"<LEVEL> | <message>"` with no ANSI escape
  sequences (contrast with `StreamFormatter`, which does add them).
- `tests/gui/test_log_handler.py` (new): `test_emit_forwards_formatted_text`
  — attach a bare `QtLogHandler()` (no `Controller`/window involved) to a
  throwaway `logging.Logger`, connect `record_logged` to a plain Python list
  via a lambda, call `logger.info('hello')`, run
  `QtCore.QCoreApplication.processEvents()`, and assert the collected text
  matches what `GuiFormatter` would produce. Keeps this test independent of
  the full GUI wiring so it fails only if the handler itself is wrong.
- `tests/gui/test_controller.py` (extend, using the existing `qapp` fixture
  from `tests/gui/conftest.py`):
  - `test_log_records_appear_in_log_tab` — build a `Controller`, call
    `nickel_focus.log.info('hello')`, `processEvents()`, assert `'hello'` is
    in `window.control_panel.log_widget.toPlainText()`.
  - `test_repeated_controller_construction_does_not_leak_handlers` —
    construct two `Controller`s in sequence against fresh windows (mirroring
    how each test in this file gets its own `Controller`) and assert
    `nickel_focus.log.handlers` contains exactly one `QtLogHandler`
    afterward. This is a regression test for the stale-handler gotcha above;
    without the fix it would either fail this assertion or raise/segfault
    when the first (dangling) handler tries to append into a destroyed
    `log_widget`.

---

## Phase 2 — Consolidate the Log tab around real logging calls

Retire the hand-built narration in favor of the handler from Phase 1 being
the *only* thing that populates `log_widget`.

- In `focus_control_panel.py`, the `show_*`/`update_step` methods
  (lines 946-1024) currently do two things each: update a one-line status
  widget (`status_label`/`step_label` — genuine live UI *state*, not log
  history) **and** call `self.log_widget.appendPlainText(...)` (log
  *history*, now redundant with Phase 1's handler). Strip the second part
  from each of: `update_step`, `show_nearest_target`, `show_slew_result`,
  `show_best_focus`, `show_single_exposure_result`, `show_failure`.
- In `gui/controller.py`, at each call site of those methods, add the
  matching `log.*()` call so the Log tab keeps getting the same information,
  now via the real logger:
  - `show_failure(message)` call sites → also `log.error(message)`
  - `_on_slew_finished`/`show_slew_result` → `log.info(...)`
  - `_on_focus_sequence_finished`/`show_best_focus` → `log.info(...)`
  - `_on_single_exposure_finished`/`show_single_exposure_result` → `log.info(...)`
  - `_find_nearest`/`show_nearest_target` → `log.info(...)`
  - `_on_step_complete`/`update_step` → `log.info(...)` (matches today's
    always-visible behavior; `QtLogHandler`'s default level of INFO from
    Phase 1 already covers this, no change needed there)
- Result: `focus_control_panel.py`'s view methods go back to being pure
  "update this label" code; `log_widget`'s entire content is sourced from one
  place (whatever calls `log.*()`), matching how the file/console handlers
  already work.

**Verification**: re-run the Phase 1 GUI test plus `test_controller.py`'s
existing signal-handling tests (they assert on `show_failure` etc. being
called/status text — update any assertions that checked `log_widget` content
directly via the old `appendPlainText` path, since that content now arrives
asynchronously through the logging handler instead).

**New unit tests**:
- `tests/gui/test_focus_control_panel.py` (extend): for each of `update_step`,
  `show_nearest_target`, `show_slew_result`, `show_best_focus`,
  `show_single_exposure_result`, `show_failure` — assert calling the method
  directly updates `status_label`/`step_label` as before but leaves
  `log_widget` untouched (e.g. `log_widget.toPlainText() == ''` beforehand and
  after). This pins down that the view is now pure label-state and locks in
  the intended behavior change from Phase 1's version of these tests.
- `tests/gui/test_controller.py` (extend), one test per call site added in
  this phase, using `caplog` at the appropriate level:
  - `test_focus_sequence_failure_logs_error` — drive `focus_worker.
    focusSequenceFailed.emit('boom')` (or trigger the equivalent failure
    path directly) and assert a `log.error` record with that message via
    `caplog.records`.
  - `test_slew_failure_logs_error`, `test_slew_finished_logs_info`,
    `test_focus_sequence_finished_logs_info`,
    `test_single_exposure_finished_logs_info`,
    `test_find_nearest_target_logs_info`, `test_step_complete_logs_info` —
    same pattern for each remaining call site, confirming the right level
    and message content.
- `tests/gui/test_controller.py`: `test_log_tab_reflects_worker_signal_end_to_end`
  — an integration-style test in this same file: emit a worker signal (e.g.
  `focusSequenceFailed`), `processEvents()`, and assert the resulting text
  appears in `window.control_panel.log_widget.toPlainText()`, confirming the
  full path (worker signal → Controller → `log.error` → `QtLogHandler` →
  Log tab) works end-to-end now that the direct `appendPlainText` path is
  gone.

---

## Phase 3 — Migrate domain-layer `print()` to `logging`

Now that there's a real consumer (the GUI Log tab) for `log.debug()` output,
finish wiring the domain layer:

- `focus.py` (lines 85, 89, 93, 95, 103, 215): unconditional narration of the
  secondary-mirror focus move → `log.info(...)`.
- `photometry.py`: the `verbose=`-gated diagnostics (lines 82, 89, 96, 122,
  123, 236, 247) → `log.debug(...)`; the two `"Warning: ..."` lines (127, 137)
  → `log.warning(...)`. Since `logging` already filters by level, drop the
  `verbose` parameter from `find_sources`, `evaluate_shape`, `evaluate_sources`,
  `image_quality` (photometry.py:27, 102, 211, 280) — check `focus.py:559`'s
  separate `execute(verbose=True, ...)` parameter too, since it may be a
  distinct, older CLI-only verbosity switch that needs its own resolution
  (either wire it to `log.setLevel` or drop it if `execute()` is otherwise
  superseded by `FocusSequence.step()`).
- Update `nickel_focus/__init__.py`'s module-level `log = get_logger(...)`
  usage pattern: no change needed there, but every module now calling
  `log.*()` should do `from nickel_focus import log` at the top (matching
  `slew.py`'s existing pattern) — per the imports-at-top-of-file convention,
  not an inline import inside functions.

**Verification**: run `nickel_focus_gui --test`, trigger a focus move/step
with `-v 2`-equivalent (DEBUG) logging enabled, confirm the previously-dead
diagnostic messages now appear in both the terminal and the Log tab.

**New unit tests**:
- `tests/test_focus.py` (extend, or new if this module doesn't already have a
  logging-focused test): `test_secondary_focus_move_logs_info` — using
  `caplog` at `INFO`, exercise the secondary-mirror focus-change path and
  assert the migrated messages (e.g. `'Unlocking secondary'`,
  `'Successfully changed focus to ...'`) appear as `log.info` records; also
  assert nothing is written to stdout any more
  (`capsys.readouterr().out == ''`), since this replaces a `print()`-based
  test if one currently exists for this code path.
- `tests/test_photometry.py` (extend, or new): 
  - `test_find_sources_debug_logging` — with `caplog.at_level(logging.DEBUG)`,
    call `find_sources(...)` and assert the previously `verbose`-gated
    messages (`'Updated background: ...'`, `'Updated threshold: ...'`,
    `'Source detection converged after ...'`) now appear unconditionally as
    `log.debug` records — this is the direct regression test proving the
    "dead diagnostics" gap described in Context is fixed.
  - `test_find_sources_silent_at_info_level` — same call with the logger at
    its default `INFO` level and assert none of those messages are captured,
    confirming the messages are gated by level rather than always emitted.
  - `test_evaluate_shape_warns_on_degenerate_sigma` and
    `test_evaluate_shape_warns_on_dissimilar_sigma` — construct inputs that
    hit each of the two degenerate-sigma branches (photometry.py:125-137) and
    assert a `log.warning` record is produced (via `caplog`), replacing the
    removed `verbose`-gated print branches.
  - `test_verbose_parameter_removed` — use `inspect.signature` on
    `find_sources`, `evaluate_shape`, `evaluate_sources`, `image_quality` to
    assert none of them still accept a `verbose` parameter, guarding against
    the flag being silently reintroduced later.
- Whatever resolution is chosen for `focus.py:559`'s `execute(verbose=True, ...)`
  parameter gets a matching test: either
  `test_execute_verbose_controls_log_level` (if wired to `log.setLevel`) or a
  simple confirmation that `execute()` no longer accepts `verbose` (if
  dropped), consistent with whichever approach is decided during
  implementation of this phase.

---

## Phase 4 — Migrate CLI scripts' `print()` and fix their tests

- `scripts/focus.py` (lines 180, 181, 197) and `scripts/slew_to_nearest.py`
  (line 80) → `log.info(...)`. These scripts already call
  `scriptbase.ScriptBase.init_log(args)` (`scripts/scriptbase.py:129-136`),
  which already derives the logging level from `-v/--verbosity` — so once
  these become `log.*()` calls, the existing CLI flag controls them for free.
- This breaks the `capsys`-based assertions in `test_focus_cli.py`,
  `test_slew_to_nearest_cli.py`, and `test_focus_gui_cli.py` (they currently
  do `capsys.readouterr()` and check stdout text). Switch those specific
  assertions to pytest's `caplog` fixture (`caplog.text` / `caplog.records`),
  which is the standard way to assert on `logging` output instead of stdout.
  This is the concrete resolution to the `capsys` blocker noted as deferred
  in project memory.
- `scripts/focus.py`'s own separate `--verbose` flag (`scripts/focus.py:96-97`,
  distinct from `scriptbase`'s `-v/--verbosity`) most likely becomes
  redundant once `photometry.py`'s `verbose` param is gone (Phase 3) — confirm
  its only remaining use and either remove it or fold it into `-v`.

**Verification**: `pytest nickel_focus/tests/test_focus_cli.py
nickel_focus/tests/test_slew_to_nearest_cli.py nickel_focus/tests/test_focus_gui_cli.py`
after switching to `caplog`, plus a manual run of both CLI scripts at a couple
of `-v` levels to confirm console output still matches expectations.

**New unit tests**:
- `test_focus_cli.py::test_cli_archive_mode_runs_end_to_end` and
  `test_cli_writes_output_file_when_requested`, `test_slew_to_nearest_cli.py::
  test_cli_dry_run_does_not_move_the_telescope` and
  `test_cli_search_string_selects_matching_target`, and
  `test_focus_gui_cli.py::test_test_flag_is_suppressed_from_help` — rewritten
  in place to use `caplog` instead of `capsys` wherever they were asserting on
  the migrated `print()` messages specifically (the `test_test_flag_is_suppressed_from_help`
  case only needs to change if it happens to assert on affected output;
  otherwise leave it on `capsys` since it's testing `argparse`'s own help
  text, not application logging).
- New: `test_cli_verbosity_flag_controls_log_level` (in `test_focus_cli.py`
  and/or `test_slew_to_nearest_cli.py`) — run the script's `main()` with
  `-v 0`, `-v 1`, and `-v 2` and use `caplog` to assert that DEBUG-level
  messages (now real, thanks to Phase 3) appear only at `-v 2`, INFO messages
  appear at `-v 1` and `-v 2` but not `-v 0`, and WARNING/ERROR appear at
  every level — this is the first test to actually exercise
  `scriptbase.py`'s existing `-v/--verbosity` → logging-level wiring against
  real emitted messages, since before this migration there was nothing at
  DEBUG/INFO for it to filter.
- A test for whichever resolution is chosen for `scripts/focus.py`'s separate
  `--verbose` flag (`scripts/focus.py:96-97`): either
  `test_verbose_flag_removed_from_parser` (confirms the redundant flag is
  gone from `get_parser()`) or a test confirming it's folded into `-v`,
  matching the Phase 3 resolution of `photometry.py`'s `verbose` parameter.

---

## Notes on scope not included above

Polish ideas (per-level color in `log_widget`, a live level-filter control,
hooking up `NickelFocusLogger`'s existing `FileHandler` support to a
GUI "save log to file" action) are natural follow-ons but out of scope for
these four phases — flagged here so they aren't forgotten, not because they're
required for #1/#2.

---

## Change log

### Phase 1 — done (2026-09-01)

Implemented as planned, plus one unplanned fix discovered while testing:

- Added `GuiFormatter` to `pkg/logger.py`, `QtLogHandler` in new file
  `gui/log_handler.py`, `FocusControlPanel.append_log_line()`, and wired it
  all up in `Controller.__init__` (including the stale-handler-removal
  guard for repeated `Controller` construction).
- Added tests: `tests/test_logger.py`, `tests/gui/test_log_handler.py`, and
  two new tests in `tests/gui/test_controller.py`
  (`test_log_records_appear_in_log_tab`,
  `test_repeated_controller_construction_does_not_leak_handlers`). Also
  added an autouse `_remove_qt_log_handlers` fixture to `tests/gui/conftest.py`
  so a GUI test's handler never lingers on the `nickel_focus.log` singleton
  into a later, non-GUI test module.
- **Unplanned fix**: `NickelFocusLogger.init()` (`pkg/logger.py`) set levels
  on its handlers but never called `self.setLevel(...)` on the logger
  itself. Since a logger's own effective level gates every call *before*
  any handler is even consulted, and nothing in this codebase set a level
  anywhere else in the parent chain either, every logger fell back to the
  root logger's default (`WARNING`) — meaning `log.info()`/`log.debug()`
  calls were silently dropped everywhere in the app, not just the GUI (this
  predates this branch; confirmed via `log.getEffectiveLevel()` returning
  `30` even after `get_logger(level=logging.INFO)`). Fixed by adding
  `self.setLevel(logging.DEBUG)` — the *most permissive* value, not the
  passed-in `level` — so the logger's own gate never filters anything, and
  all real filtering happens per-handler (console at `level`, file at
  `log_file_level`, the GUI handler at its own level). Pinning the logger to
  the passed-in `level` instead would have silently reintroduced the same
  bug for any handler configured more verbose than that value (e.g.
  `log_file_level` more verbose than the console).
- Verified manually end-to-end (`nickel_focus_gui`, offscreen): both
  `log.info` and `log.warning` calls now appear in the Log tab via
  `GuiFormatter`'s plain (non-ANSI) text, while the console still shows
  `StreamFormatter`'s colorized version of the same records.
- Full suite: `209 passed` (`pytest nickel_focus/tests`).

### Phase 2 — done (2026-09-01)

Implemented as planned, no changes to the plan itself:

- Stripped the `self.log_widget.appendPlainText(...)` call out of each of
  `FocusControlPanel.update_step`, `show_nearest_target`, `show_slew_result`,
  `show_best_focus`, `show_single_exposure_result`, `show_failure` — these
  now only update `status_label`/`step_label`, matching the panel's other
  pure-display methods (`set_current_position`, etc.).
- Added a matching `log.*()` call in `Controller` at every corresponding
  call site: `log.error(...)` alongside all 9 `show_failure(...)` call
  sites (in `start_focus_sequence`, `take_single_exposure`, `move_to_target`,
  `_find_nearest`, `_on_focus_sequence_failed`, `_on_slew_failed`);
  `log.info(...)` in `_on_slew_finished`, `_on_focus_sequence_finished`,
  `_on_single_exposure_finished`, `_find_nearest` (success path), and
  `_on_step_complete`. The step/best-focus/single-exposure messages
  intentionally duplicate the same short f-string the View used to build
  internally (rather than introducing a shared formatter/changing method
  signatures to accept pre-formatted text) — each is 1-6 lines, and this
  keeps Phase 2 scoped to what was planned; worth revisiting if this
  duplication grows.
- Result: `log_widget`'s entire content is now sourced from one place
  (`Controller`'s `log.*()` calls via Phase 1's `QtLogHandler`) — verified
  manually (`_on_focus_sequence_failed('smoke test failure')` followed by
  `_on_slew_finished()` produced both an ERROR and an INFO line in the Log
  tab, while `status_label` correctly showed only the latest one).
- Updated 7 existing tests in `tests/gui/test_focus_control_panel.py` whose
  assertions checked `log_widget` content directly after calling a view
  method (now asserting `log_widget` is left untouched instead); renamed a
  few to match (e.g. `test_show_failure_updates_status_and_log` →
  `test_show_failure_updates_status_only`). `test_reset_clears_status_and_log`
  now populates `log_widget` via `append_log_line()` before calling
  `reset()`, so it still meaningfully exercises clearing a non-empty log.
- Added 8 new tests in `tests/gui/test_controller.py`, one `caplog`-based
  test per call site added above, plus
  `test_log_tab_reflects_worker_signal_end_to_end`.
- Full suite: `217 passed` (`pytest nickel_focus/tests`).

### Phase 3 — done (2026-09-01)

Implemented largely as planned, with a few resolutions/scope adjustments
made along the way (flagged individually below) and two planned tests
dropped at the user's request:

- `focus.py` (`Focus.set_to`, `Exposure.expose`): all 6 `print()` calls
  → `log.info(...)`. Added `from nickel_focus import log`.
- `photometry.py`: the `verbose=`-gated diagnostics in `find_sources`
  → `log.debug(...)` (unconditional); the two `"Warning: ..."` branches in
  `evaluate_shape` → `log.warning(...)`. Dropped the now-redundant
  `verbose` parameter from `find_sources`, `evaluate_shape`,
  `evaluate_sources`, `image_quality`. Added `from nickel_focus import log`.
- **Scope addition, done at the user's request**: while removing
  `evaluate_sources`'s `verbose` parameter, its two `if verbose:
  warnings.warn(...)` guards needed *some* resolution (the parameter they
  gated no longer exists). Initially planned to just drop the `if verbose:`
  guard and leave `warnings.warn`; the user asked for `log.warning(...)`
  instead. Applied that consistently to all `warnings.warn` calls actually
  touched by this migration: all 3 in `photometry.py` (`evaluate_sources`)
  and, for the same consistency, the 2 plain (non-verbose-gated)
  `warnings.warn` calls already in `focus.py` (`FocusPlot.add_stamp`,
  `AutomatedFocusSequence.step_focus`) — nothing else in the codebase uses
  `warnings.warn` (checked via grep). Removed the now-unused `import
  warnings` from both files.
- **Resolution**: `focus.py:559`'s separate `FocusSequence.execute(verbose=True,
  ...)` parameter turned out to be dead code -- never referenced in the
  method body, never passed by any caller (checked via grep). Dropped it
  outright rather than wiring it to anything.
- Manually verified the level-gating end-to-end (not just that the calls
  exist): `log.init(level=logging.DEBUG)` then `log.init(level=logging.INFO)`
  around the same `find_sources(...)` call shows the debug diagnostics
  appear (with full `file:func:line` context, since `DebugStreamFormatter`
  kicks in at that level) and then vanish completely -- confirming the
  Phase 1 logger-level fix and this phase's migration compose correctly,
  i.e. `-v 2`-equivalent actually surfaces detail that was previously
  unreachable at any verbosity.
- Added `tests/test_photometry.py` (new) with `caplog`-based tests for the
  `find_sources` debug messages (present at DEBUG, absent at INFO) and the
  `evaluate_shape` warning branches (degenerate/dissimilar sigma). Added
  `tests/test_focus.py` (new) with a minimal fake `ktl` module (just enough
  to drive `Focus.set_to()`'s golden path and its "already at position"
  early-return) to verify its migrated messages log at INFO -- this was
  more work than anticipated, since every existing test fixture replaces
  `focus.Focus` wholesale (`FakeFocus`) rather than exercising its real
  implementation against a faked `ktl`; this is the first test to do the
  latter. Deliberately did not build an equivalent fake for
  `Exposure.expose()`'s one migrated line -- verified manually instead
  (see above) -- to keep this new test-infrastructure investment scoped to
  what the plan actually asked for.
- Two planned tests dropped at the user's explicit request (not needed):
  `test_verbose_parameter_removed` (photometry) and a test confirming
  `execute()` dropped its `verbose` parameter.
- Full suite: `223 passed` (`pytest nickel_focus/tests`).
