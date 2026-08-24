"""
End-to-end smoke test for the GUI entry point (`gui.main`), plus a
GUI/Model equivalence check: the whole point of the `step()` generator
refactor (GUI_DESIGN.md §2) is that the CLI and the GUI drive the same
underlying `focus.FocusSequence` code, so this confirms the GUI's
`Controller` produces the identical fit as calling
`FocusSequence.execute()` directly -- the same call the CLI's `main()`
ultimately makes -- rather than re-testing `focus.main()` itself (already
covered by `test_focus_cli.py`).
"""
import pytest

import focus
import gui.main
from gui.controller import Controller
from gui.qt import QtCore


def _wait_for_worker(controller, timeout_ms=5000):
    worker = controller.worker
    if worker is None:
        return
    loop = QtCore.QEventLoop()
    worker.finished.connect(loop.quit)
    QtCore.QTimer.singleShot(timeout_ms, loop.quit)
    loop.exec()


def _configure_replay_tab(panel, focus_sweep):
    panel.tabs.setCurrentWidget(panel.replay_tab)
    panel.replay_datadir_edit.setText(str(focus_sweep['datadir']))
    panel.replay_prefix_edit.setText(focus_sweep['prefix'])
    panel.replay_suffix_edit.setText(focus_sweep['suffix'])
    panel.replay_obsnum_spin.setValue(focus_sweep['obsnum'])
    panel.replay_start_spin.setValue(int(focus_sweep['focus_values'][0]))
    panel.replay_step_spin.setValue(
        int(focus_sweep['focus_values'][1] - focus_sweep['focus_values'][0]))
    panel.replay_nstep_spin.setValue(len(focus_sweep['focus_values']))


def test_app_opens_and_shows_without_error(qapp):
    app = gui.main.build_app()
    assert app is qapp, 'build_app() should reuse the existing QApplication, not create a second one'

    window = gui.main.build_window()
    controller = Controller(window)
    window.show()
    qapp.processEvents()

    assert window.windowTitle() == 'Nickel Focus GUI'
    assert window.image_panel is not None, 'MainWindow should build a real ImagePanel'
    assert window.curve_panel is not None, 'MainWindow should build a real FocusCurvePanel'
    assert window.control_panel is not None, 'MainWindow should build a real FocusControlPanel'

    window.close()


def test_gui_controller_matches_direct_execute(qapp, focus_sweep):
    # Build the same sequence directly -- the same GridFocusSequence (for
    # its focus-value arithmetic) + ArchiveFocusSequence + execute() calls
    # that both focus.py's CLI and Controller.start_sequence()
    # make -- and confirm the GUI reproduces the identical fit.
    grid = focus.GridFocusSequence(focus_sweep['focus_values'][0],
                                    focus_sweep['focus_values'][1] - focus_sweep['focus_values'][0],
                                    nstep=len(focus_sweep['focus_values']))
    direct_seq = focus.ArchiveFocusSequence(list(grid.target_focus), focus_sweep['files'])
    expected_best_focus, expected_best_fwhm = direct_seq.execute(
        goto=False, plot=False, method='brightest')

    window = gui.main.build_window()
    controller = Controller(window)
    _configure_replay_tab(window.control_panel, focus_sweep)

    controller.start_sequence()
    fitted = []
    controller.worker.sequenceFinished.connect(lambda bf, bfw: fitted.append((bf, bfw)))
    _wait_for_worker(controller)

    assert fitted, 'sequenceFinished should have fired exactly once'
    assert fitted[0][0] == pytest.approx(expected_best_focus), \
        'the GUI Controller should produce the same fit as calling FocusSequence.execute() directly'
