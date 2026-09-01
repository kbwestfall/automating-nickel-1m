"""
End-to-end fake-hardware smoke test for the full live-mode feature set
(Phase 3, plus the Phase 4 tabbed redesign).

Unlike `test_controller.py`, which tests each live-mode action in
isolation, this drives one `Controller` through the full feature set
against `fake_hardware`, in the order an observer would actually use it
(GUI_DESIGN.md §5.5's own "typical use" narrative, adapted for the
Single tab): take a standalone exposure to confirm the field and mark a
star via a click, start a real sequence that measures that marked star,
move to its best focus via the Single tab's auto-populated default, mark
a *different* source (now that a sequence is loaded, this should
reanalyze the loaded sequence instead of a standalone exposure), take
one more standalone exposure, and finally discard it by starting a fresh
Automated sequence. If any of these actions don't compose correctly with
each other -- as opposed to each working fine alone -- this is what
would catch it.
"""
import pytest

from nickel_focus import focus
from nickel_focus.gui.qt import QtCore
from nickel_focus.gui.views.main_window import MainWindow
from nickel_focus.gui.controller import Controller


def _wait_for_focus_worker(controller, timeout_ms=5000):
    worker = controller.focus_worker
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

    # 1. Confirm the field with a standalone exposure from the Single
    #    tab -- no sequence is loaded yet.
    panel.tabs.setCurrentWidget(panel.single_tab)
    panel.single_focus_spin.setValue(340)
    panel.single_acquire_button.click()
    _wait_for_focus_worker(controller)

    assert controller.focus_sequence is None, 'no sequence should be loaded yet'
    assert controller._standalone_focus_sequence is not None, \
        'the standalone exposure should be held for possible reanalysis'
    assert len(window.image_panel._results) == 1, 'the standalone exposure should be displayed'
    assert len(window.curve_panel._results) == 0, 'not sequence data'

    # 2. Mark a star on it via a click -- with no sequence loaded, this
    #    should reanalyze that one standalone exposure in place (§5.6).
    controller._on_source_selected(50., 50.)
    _wait_for_focus_worker(controller)

    assert controller.method == (50., 50.), 'clicking a source should update the active method'
    assert len(window.image_panel._results) == 1, 'reanalysis should update in place, not duplicate'
    assert len(window.curve_panel._results) == 0, 'still not sequence data'

    # 3. Start a Grid sequence -- using the just-marked source as the
    #    active method -- which discards the standalone exposure state.
    panel.tabs.setCurrentWidget(panel.grid_tab)
    panel.grid_start_spin.setValue(340)
    panel.grid_step_spin.setValue(5)
    panel.grid_nstep_spin.setValue(5)
    panel.grid_exptime_spin.setValue(5.0)
    panel.grid_acquire_button.click()
    _wait_for_focus_worker(controller)

    assert controller._standalone_focus_sequence is None, \
        'starting a new sequence should discard standalone single-exposure state'
    assert isinstance(controller.focus_sequence, focus.GridFocusSequence), \
        'setup: a Grid sequence should be loaded'
    assert controller.focus_sequence.method == (50., 50.), \
        'the sequence should have measured the marked source, not "brightest"'
    assert len(window.image_panel._results) == 5, 'one image per exposure'
    assert len(window.curve_panel._results) == 5, 'one curve point per exposure'

    # 4. Move to best focus: the Single tab's focus value should now
    #    default to the fitted best focus. Acquiring it as-is is "move to
    #    best focus" -- a confirmation exposure is taken and displayed,
    #    but the sequence/curve data are untouched.
    panel.tabs.setCurrentWidget(panel.single_tab)
    best_focus = panel.single_focus_spin.value()
    assert best_focus == pytest.approx(fake_hardware['best_focus'], abs=5.), \
        'the Single tab should default to the fitted best focus'
    panel.single_acquire_button.click()
    _wait_for_focus_worker(controller)

    assert fake_hardware['focus'].current == best_focus, \
        'the fake Focus should have been commanded to the (rounded) fitted best focus'
    assert len(window.image_panel._results) == 6, 'the confirmation exposure should be displayed'
    assert len(window.curve_panel._results) == 5, 'a confirmation exposure is not sequence data'
    assert f'{best_focus:.0f}' in panel.status_label.text(), \
        'the result should be reported in the status'

    # 5. Mark a different source via a click. A sequence is loaded now,
    #    so this reanalyzes *it*, not the by-then-secondary standalone
    #    exposure from step 4 (whichever is "loaded" takes priority --
    #    see `Controller.reanalyze`).
    controller._on_source_selected(55., 55.)
    _wait_for_focus_worker(controller)

    assert controller.method == (55., 55.)
    assert len(controller.focus_sequence.exposures) == 5, \
        "reanalysis updates in place, doesn't add"
    assert len(window.curve_panel._results) == 5, 'the loaded sequence was reanalyzed, not added to'
    assert len(window.image_panel._results) == 6, 'no duplicates from reanalysis'

    # 6. Take one more standalone exposure, then discard it by starting a
    #    fresh Automated sequence -- confirming the earlier sequence's
    #    data is fully replaced and standalone state is cleared.
    panel.single_focus_spin.setValue(370)
    panel.single_acquire_button.click()
    _wait_for_focus_worker(controller)
    assert controller._standalone_focus_sequence is not None, \
        'setup: a standalone exposure should exist'

    panel.tabs.setCurrentWidget(panel.auto_tab)
    panel.auto_start_spin.setValue(340)
    panel.auto_step_spin.setValue(5)
    panel.auto_maxsteps_spin.setValue(12)
    panel.auto_acquire_button.click()
    _wait_for_focus_worker(controller)

    assert controller._standalone_focus_sequence is None, \
        'starting a new sequence should discard standalone single-exposure state'
    assert isinstance(controller.focus_sequence, focus.AutomatedFocusSequence), \
        'the new sequence should have replaced the old Grid sequence'
    assert panel.single_focus_spin.value() == pytest.approx(fake_hardware['best_focus'], abs=5), \
        'the adaptive search should converge near the known best focus'
