"""
End-to-end fake-hardware smoke test closing out Phase 3.

Unlike `test_controller.py`, which tests each live-mode action in
isolation, this drives one `Controller` through the full Phase 3 feature
set against `fake_hardware`, in the order an observer would actually use
it (GUI_DESIGN.md §5.5's own "typical use" narrative): take a standalone
exposure to confirm the field and mark a star via a click, start a real
sequence that measures that marked star, move to its best focus, take
another standalone exposure, mark a *different* source (now that a
sequence is loaded, this should reanalyze the loaded sequence instead of
the by-then-irrelevant pending exposure), commit the pending exposure
into the sequence, and finally discard one last pending exposure by
starting a fresh Automated sequence. If any of these actions don't
compose correctly with each other -- as opposed to each working fine
alone -- this is what would catch it.
"""
import pytest

import focus
from gui.qt import QtCore
from gui.views.main_window import MainWindow
from gui.controller import Controller


def _wait_for_worker(controller, timeout_ms=5000):
    worker = controller.worker
    if worker is None:
        return
    loop = QtCore.QEventLoop()
    worker.finished.connect(loop.quit)
    QtCore.QTimer.singleShot(timeout_ms, loop.quit)
    loop.exec()


def test_phase3_full_live_workflow(qapp, fake_hardware):
    window = MainWindow()
    controller = Controller(window)
    panel = window.control_panel

    # 1. Confirm the field with a standalone exposure -- no sequence is
    #    loaded yet.
    panel.takeSingleExposureRequested.emit(340.)
    _wait_for_worker(controller)

    assert controller.sequence is None, 'no sequence should be loaded yet'
    assert controller.pending_result is not None, 'the standalone exposure should be held pending'
    assert len(window.image_panel._results) == 1, 'the standalone exposure should be displayed'
    assert len(window.curve_panel._results) == 0, 'not sequence data until committed'
    assert not panel.add_to_sequence_button.isEnabled(), \
        'nothing is loaded to add the pending exposure to'

    # 2. Mark a star on it via a click -- with no sequence loaded, this
    #    should reanalyze that one standalone exposure in place (§5.6).
    original_pending = controller.pending_result
    controller._on_source_selected(50., 50.)
    _wait_for_worker(controller)

    assert controller.method == (50., 50.), 'clicking a source should update the active method'
    assert controller.pending_result is not original_pending, \
        'reanalysis should produce a fresh measurement of the standalone exposure'
    assert len(window.image_panel._results) == 1, 'reanalysis should update in place, not duplicate'
    assert len(window.curve_panel._results) == 0, 'still not sequence data'

    # 3. Start a Grid sequence -- using the just-marked source as the
    #    active method -- which discards the pending exposure (it's not
    #    automatically counted as sequence data, per §5.5).
    panel.grid_radio.setChecked(True)
    panel.start_spin.setValue(340.)
    panel.step_spin.setValue(5.)
    panel.nstep_spin.setValue(5)
    panel.exptime_spin.setValue(5.0)
    controller.start_sequence()
    _wait_for_worker(controller)

    assert controller.pending_result is None, 'starting a new sequence should discard the pending exposure'
    assert isinstance(controller.sequence, focus.GridFocusSequence), 'setup: a Grid sequence should be loaded'
    assert controller.sequence.method == (50., 50.), \
        'the sequence should have measured the marked source, not "brightest"'
    assert len(window.image_panel._results) == 5, 'one image per exposure'
    assert len(window.curve_panel._results) == 5, 'one curve point per exposure'
    best_focus = panel._best_focus
    assert best_focus is not None, 'the sequence should have fit a best focus'
    assert panel.move_to_best_focus_button.isEnabled(), \
        'a live sequence should allow moving to best focus'

    # 4. Move to best focus: a confirmation exposure is taken and
    #    displayed, but the sequence/curve data are untouched.
    panel.moveToBestFocusRequested.emit(best_focus)
    _wait_for_worker(controller)

    assert fake_hardware['focus'].current == best_focus, \
        'the fake Focus should have been commanded to the fitted best focus'
    assert len(window.image_panel._results) == 6, 'the confirmation exposure should be displayed'
    assert len(window.curve_panel._results) == 5, 'a confirmation exposure is not sequence data'
    assert f'{best_focus:.0f}' in panel.status_label.text(), \
        'the confirmation should be reported in the status'

    # 5. Take another standalone single exposure -- e.g. a post-hoc
    #    check -- independent of the loaded sequence.
    panel.takeSingleExposureRequested.emit(362.)
    _wait_for_worker(controller)

    assert controller.pending_result is not None, 'the standalone exposure should be held pending'
    assert controller.pending_result.focus_value == 362., \
        'the fake Focus should have been moved to the requested value'
    assert len(controller.sequence.exposures) == 5, \
        'the standalone exposure must not touch the loaded sequence'
    assert len(window.image_panel._results) == 7, 'the standalone exposure should be displayed'
    assert len(window.curve_panel._results) == 5, 'still not sequence data until committed'
    assert panel.add_to_sequence_button.isEnabled(), \
        'a sequence is loaded, so committing the pending exposure should be available'

    # 6. Mark a different source via a click. A sequence is loaded now,
    #    so this reanalyzes *it*, not the by-then-secondary pending
    #    exposure (whichever is "loaded" takes priority -- see
    #    `Controller.reanalyze`).
    pending_before_click = controller.pending_result
    controller._on_source_selected(55., 55.)
    _wait_for_worker(controller)

    assert controller.method == (55., 55.)
    assert len(controller.sequence.exposures) == 5, 'reanalysis updates in place, doesn\'t change count'
    assert len(window.curve_panel._results) == 5, 'the loaded sequence was reanalyzed, not added to'
    assert len(window.image_panel._results) == 7, 'no duplicates from reanalysis'
    assert controller.pending_result is pending_before_click, \
        'the pending standalone exposure should be untouched by reanalyzing the loaded sequence'

    # 7. Commit the pending exposure into the loaded sequence.
    panel.addToSequenceRequested.emit()

    assert controller.pending_result is None, 'the pending result should be cleared once committed'
    assert len(controller.sequence.exposures) == 6, 'the sequence should gain the committed exposure'
    assert len(window.curve_panel._results) == 6, 'the curve panel should gain the committed point'
    assert panel.pending_label.text() == 'No pending exposure'

    # 8. Take one more standalone exposure, then discard it by starting
    #    a fresh Automated sequence -- confirming the earlier sequence's
    #    (now 6-point) data is fully replaced and pending state cleared.
    panel.takeSingleExposureRequested.emit(370.)
    _wait_for_worker(controller)
    assert controller.pending_result is not None, 'setup: a second pending exposure should exist'

    panel.automated_radio.setChecked(True)
    panel.start_spin.setValue(340.)
    panel.step_spin.setValue(5.)
    panel.maxsteps_spin.setValue(12)
    controller.start_sequence()
    _wait_for_worker(controller)

    assert controller.pending_result is None, \
        'starting a new sequence should discard the pending exposure'
    assert isinstance(controller.sequence, focus.AutomatedFocusSequence), \
        'the new sequence should have replaced the old Grid sequence'
    assert panel._best_focus == pytest.approx(fake_hardware['best_focus'], abs=5.), \
        'the adaptive search should converge near the known best focus'
