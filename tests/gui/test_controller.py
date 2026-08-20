"""Tests for :mod:`gui.controller`."""
import types

import pytest

import focus
from gui.qt import QtCore
from gui.views.main_window import MainWindow
from gui.controller import Controller


def _configure_control_panel(panel, focus_sweep):
    panel.datadir_edit.setText(str(focus_sweep['datadir']))
    panel.prefix_edit.setText(focus_sweep['prefix'])
    panel.suffix_edit.setText(focus_sweep['suffix'])
    panel.obsnum_spin.setValue(focus_sweep['obsnum'])
    panel.start_spin.setValue(focus_sweep['focus_values'][0])
    panel.step_spin.setValue(focus_sweep['focus_values'][1] - focus_sweep['focus_values'][0])
    panel.nstep_spin.setValue(len(focus_sweep['focus_values']))


def _wait_for_worker(controller, timeout_ms=5000):
    worker = controller.worker
    if worker is None:
        return
    loop = QtCore.QEventLoop()
    worker.finished.connect(loop.quit)
    QtCore.QTimer.singleShot(timeout_ms, loop.quit)
    loop.exec()


def _make_controller(focus_sweep):
    window = MainWindow()
    controller = Controller(window)
    _configure_control_panel(window.control_panel, focus_sweep)
    return window, controller


def test_start_sequence_runs_to_completion_for_archive(qapp, focus_sweep):
    window, controller = _make_controller(focus_sweep)

    controller.start_sequence()
    _wait_for_worker(controller)

    assert controller.worker is None, 'worker should be cleared once finished'
    assert len(window.image_panel._results) == len(focus_sweep['files']), \
        'image panel should have one entry per exposure'
    assert len(window.curve_panel._results) == len(focus_sweep['files']), \
        'curve panel should have one entry per exposure'
    assert window.control_panel._best_focus == pytest.approx(focus_sweep['best_focus'], abs=5.), \
        'control panel should display the fitted best focus'


def test_running_state_is_set_synchronously_before_worker_runs(qapp, focus_sweep):
    window, controller = _make_controller(focus_sweep)

    controller.start_sequence()

    # _set_running(True) happens before worker.start(), so this should
    # already be true without needing to pump the event loop at all.
    assert not window.control_panel.start_button.isEnabled(), \
        'Start should be disabled as soon as a sequence begins'
    assert not window.image_panel._selection_enabled, \
        'image selection should be disabled while a sequence is running'

    _wait_for_worker(controller)
    assert window.control_panel.start_button.isEnabled(), 'Start should re-enable once finished'
    assert window.image_panel._selection_enabled, 'selection should re-enable once finished'


def test_hardware_exclusivity_blocks_second_start(qapp, focus_sweep):
    window, controller = _make_controller(focus_sweep)
    controller.worker = types.SimpleNamespace()  # pretend something is already running

    controller.start_sequence()

    assert controller.sequence is None, \
        'a second start should be a no-op while something is already running'


def test_missing_files_reports_failure(qapp, focus_sweep):
    window, controller = _make_controller(focus_sweep)
    window.control_panel.prefix_edit.setText('does-not-exist-')

    controller.start_sequence()

    assert controller.worker is None, 'a config error should never start a worker'
    assert 'Could not start sequence' in window.control_panel.status_label.text(), \
        'the failure should be reported to the user'


def test_stop_delegates_to_the_worker(qapp):
    window = MainWindow()
    controller = Controller(window)
    stop_calls = []
    controller.worker = types.SimpleNamespace(request_stop=lambda: stop_calls.append(True))

    controller.stop()

    assert stop_calls == [True], 'stop() should call request_stop() on the active worker'
    assert 'Stopping' in window.control_panel.status_label.text(), \
        'stop() should show the "stopping" status'


def test_stop_is_a_noop_with_nothing_running(qapp):
    window = MainWindow()
    controller = Controller(window)
    controller.stop()  # should not raise
    assert window.control_panel.status_label.text() == ''


def test_set_method_updates_state_and_display(qapp):
    window = MainWindow()
    controller = Controller(window)

    controller.set_method('weighted')

    assert controller.method == 'weighted'
    assert window.control_panel.method_label.text() == 'Current method: Weighted'


def test_reanalyze_is_a_noop_without_a_sequence(qapp):
    window = MainWindow()
    controller = Controller(window)
    controller.reanalyze()
    assert controller.worker is None


def test_source_selection_updates_measurements_in_place(qapp, focus_sweep):
    window, controller = _make_controller(focus_sweep)
    controller.start_sequence()
    _wait_for_worker(controller)

    n_before = len(window.image_panel._results)

    controller._on_source_selected(50., 50.)
    _wait_for_worker(controller)

    assert controller.method == (50., 50.), 'clicking a source should update the active method'
    assert window.control_panel.method_label.text() == 'Current method: Selected source (50.0, 50.0)'
    assert len(window.image_panel._results) == n_before, \
        'reanalysis should update existing entries, not add duplicates'
    assert len(window.curve_panel._results) == n_before, \
        'reanalysis should update existing entries, not add duplicates'


def test_start_sequence_runs_grid_against_fake_hardware(qapp, fake_hardware):
    window = MainWindow()
    controller = Controller(window)
    window.control_panel.grid_radio.setChecked(True)
    window.control_panel.start_spin.setValue(340.)
    window.control_panel.step_spin.setValue(5.)
    window.control_panel.nstep_spin.setValue(5)
    window.control_panel.exptime_spin.setValue(7.5)

    controller.start_sequence()
    _wait_for_worker(controller)

    assert len(window.image_panel._results) == 5, 'image panel should have one entry per exposure'
    assert window.control_panel._best_focus == pytest.approx(fake_hardware['best_focus'], abs=5.), \
        'control panel should display the fitted best focus'
    assert fake_hardware['focus'].current == 360., 'the fake Focus should have been commanded'
    assert controller.sequence._exposure.cfg.exptime == 7.5, \
        'exposure settings from the control panel should have been applied'


def test_start_sequence_runs_automated_against_fake_hardware(qapp, fake_hardware):
    window = MainWindow()
    controller = Controller(window)
    window.control_panel.automated_radio.setChecked(True)
    window.control_panel.start_spin.setValue(340.)
    window.control_panel.step_spin.setValue(5.)
    window.control_panel.maxsteps_spin.setValue(12)

    controller.start_sequence()
    _wait_for_worker(controller)

    assert isinstance(controller.sequence, focus.AutomatedFocusSequence)
    assert window.control_panel._best_focus == pytest.approx(fake_hardware['best_focus'], abs=5.), \
        'the adaptive search should converge near the known best focus'


def test_move_to_best_focus_is_a_noop_without_a_sequence(qapp):
    window = MainWindow()
    controller = Controller(window)
    controller.move_to_best_focus(356.)
    assert controller.worker is None, 'no sequence has run yet, so there is nothing to move'


def test_move_to_best_focus_blocks_while_something_is_running(qapp, focus_sweep):
    window, controller = _make_controller(focus_sweep)
    controller.sequence = object()
    controller.worker = types.SimpleNamespace()  # pretend something is already running

    controller.move_to_best_focus(356.)

    assert isinstance(controller.worker, types.SimpleNamespace), \
        'a move-to-best-focus request should be ignored while a sequence is running'


def test_move_to_best_focus_runs_against_fake_hardware(qapp, fake_hardware):
    window = MainWindow()
    controller = Controller(window)
    window.control_panel.grid_radio.setChecked(True)
    window.control_panel.start_spin.setValue(340.)
    window.control_panel.step_spin.setValue(5.)
    window.control_panel.nstep_spin.setValue(5)
    window.control_panel.exptime_spin.setValue(7.5)
    controller.start_sequence()
    _wait_for_worker(controller)
    n_before = len(window.image_panel._results)
    best_focus = window.control_panel._best_focus

    window.control_panel.moveToBestFocusRequested.emit(best_focus)
    _wait_for_worker(controller)

    assert controller.worker is None, 'worker should be cleared once the confirmation finishes'
    assert fake_hardware['focus'].current == best_focus, \
        'the fake Focus should have been commanded to the best-focus value'
    assert len(window.image_panel._results) == n_before + 1, \
        'the confirmation exposure should be added to the image panel'
    assert f'{best_focus:.1f}' in window.control_panel.status_label.text(), \
        'the confirmation should be reported in the status'


def test_start_sequence_without_ktl_reports_clear_failure(qapp):
    # No fake_hardware fixture here: this dev machine has no ktl, so a
    # live sequence type should fail immediately and clearly, without
    # ever starting a worker thread.
    window = MainWindow()
    controller = Controller(window)
    window.control_panel.grid_radio.setChecked(True)
    window.control_panel.nstep_spin.setValue(5)

    controller.start_sequence()

    assert controller.worker is None, 'no worker should start without a ktl connection'
    assert 'no ktl connection' in window.control_panel.status_label.text().lower()
