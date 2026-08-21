"""Tests for :mod:`gui.controller`."""
import types

import pytest

import focus
from gui.qt import QtCore
from gui.views.main_window import MainWindow
from gui.controller import Controller


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
    _configure_replay_tab(window.control_panel, focus_sweep)
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
    best_focus = round(focus_sweep['best_focus'])
    assert window.control_panel.single_focus_spin.value() == pytest.approx(best_focus, abs=5), \
        'the Single tab should default to the fitted best focus'


def test_running_state_is_set_synchronously_before_worker_runs(qapp, focus_sweep):
    window, controller = _make_controller(focus_sweep)

    controller.start_sequence()

    # _set_running(True) happens before worker.start(), so this should
    # already be true without needing to pump the event loop at all.
    assert not window.control_panel.replay_load_button.isEnabled(), \
        'Load should be disabled as soon as a sequence begins'
    assert not window.image_panel._selection_enabled, \
        'image selection should be disabled while a sequence is running'
    assert window.control_panel.tabs.isEnabled(), 'the tab bar should stay reachable while running'

    _wait_for_worker(controller)
    assert window.control_panel.replay_load_button.isEnabled(), 'Load should re-enable once finished'
    assert window.image_panel._selection_enabled, 'selection should re-enable once finished'


def test_hardware_exclusivity_blocks_second_start(qapp, focus_sweep):
    window, controller = _make_controller(focus_sweep)
    controller.worker = types.SimpleNamespace()  # pretend something is already running

    controller.start_sequence()

    assert controller.sequence is None, \
        'a second start should be a no-op while something is already running'


def test_missing_files_reports_failure(qapp, focus_sweep):
    window, controller = _make_controller(focus_sweep)
    window.control_panel.replay_prefix_edit.setText('does-not-exist-')

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


def test_set_method_updates_state(qapp):
    window = MainWindow()
    controller = Controller(window)

    controller.set_method((123.4, 567.8))

    assert controller.method == (123.4, 567.8)


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
    assert len(window.image_panel._results) == n_before, \
        'reanalysis should update existing entries, not add duplicates'
    assert len(window.curve_panel._results) == n_before, \
        'reanalysis should update existing entries, not add duplicates'


def test_start_sequence_runs_grid_against_fake_hardware(qapp, fake_hardware):
    window = MainWindow()
    controller = Controller(window)
    panel = window.control_panel
    panel.tabs.setCurrentWidget(panel.grid_tab)
    panel.grid_start_spin.setValue(340)
    panel.grid_step_spin.setValue(5)
    panel.grid_nstep_spin.setValue(5)
    panel.grid_exptime_spin.setValue(7.5)

    controller.start_sequence()
    _wait_for_worker(controller)

    assert len(window.image_panel._results) == 5, 'image panel should have one entry per exposure'
    assert panel.single_focus_spin.value() == pytest.approx(fake_hardware['best_focus'], abs=5), \
        'the Single tab should default to the fitted best focus'
    assert fake_hardware['focus'].current == 360., 'the fake Focus should have been commanded'
    assert controller.sequence._exposure.cfg.exptime == 7.5, \
        'exposure settings from the control panel should have been applied'


def test_start_sequence_runs_automated_against_fake_hardware(qapp, fake_hardware):
    window = MainWindow()
    controller = Controller(window)
    panel = window.control_panel
    panel.tabs.setCurrentWidget(panel.auto_tab)
    panel.auto_start_spin.setValue(340)
    panel.auto_step_spin.setValue(5)
    panel.auto_maxsteps_spin.setValue(12)

    controller.start_sequence()
    _wait_for_worker(controller)

    assert isinstance(controller.sequence, focus.AutomatedFocusSequence)
    assert panel.single_focus_spin.value() == pytest.approx(fake_hardware['best_focus'], abs=5), \
        'the adaptive search should converge near the known best focus'


def test_start_sequence_without_ktl_reports_clear_failure(qapp):
    # No fake_hardware fixture here: this dev machine has no ktl, so a
    # live sequence type should fail immediately and clearly, without
    # ever starting a worker thread.
    window = MainWindow()
    controller = Controller(window)
    panel = window.control_panel
    panel.tabs.setCurrentWidget(panel.grid_tab)
    panel.grid_nstep_spin.setValue(5)

    controller.start_sequence()

    assert controller.worker is None, 'no worker should start without a ktl connection'
    assert 'no ktl connection' in window.control_panel.status_label.text().lower()


def test_take_single_exposure_is_a_noop_while_something_is_running(qapp):
    window = MainWindow()
    controller = Controller(window)
    controller.worker = types.SimpleNamespace()  # pretend something is already running

    controller.take_single_exposure(340.)

    assert controller._standalone_sequence is None, \
        'a take-single-exposure request should be ignored while something is running'


def test_take_single_exposure_without_ktl_reports_clear_failure(qapp):
    window = MainWindow()
    controller = Controller(window)

    controller.take_single_exposure(340.)

    assert controller.worker is None, 'no worker should start without a ktl connection'
    assert 'no ktl connection' in window.control_panel.status_label.text().lower()


def test_take_single_exposure_runs_against_fake_hardware(qapp, fake_hardware):
    window = MainWindow()
    controller = Controller(window)
    panel = window.control_panel
    panel.tabs.setCurrentWidget(panel.single_tab)
    panel.single_exptime_spin.setValue(6.0)

    panel.takeSingleExposureRequested.emit(345.)
    _wait_for_worker(controller)

    assert controller.worker is None, 'worker should be cleared once the exposure finishes'
    assert fake_hardware['focus'].current == 345., 'the fake Focus should have been commanded'
    assert controller.sequence is None, \
        'a standalone single exposure should not create/replace a loaded sequence'
    assert len(window.image_panel._results) == 1, \
        'the exposure should be displayed in the image panel'
    assert window.curve_panel._results == [], \
        'a standalone exposure should not appear on the focus curve'
    assert f'{345:.0f}' in panel.status_label.text(), \
        'the result should be reported in the status'


def test_source_selection_reanalyzes_a_standalone_exposure(qapp, fake_hardware):
    window = MainWindow()
    controller = Controller(window)
    controller.take_single_exposure(345.)
    _wait_for_worker(controller)
    original_exposure = window.image_panel._current.exposure

    controller._on_source_selected(50., 50.)
    _wait_for_worker(controller)

    assert controller.method == (50., 50.), 'clicking a source should update the active method'
    assert len(window.image_panel._results) == 1, \
        'reanalysis should update the exposure in place, not add a duplicate'
    assert window.image_panel._current.exposure == original_exposure, \
        'still the same exposure, just reanalyzed'
    assert window.curve_panel._results == [], \
        'a standalone reanalyzed exposure still should not appear on the focus curve'


def test_start_sequence_discards_standalone_sequence_state(qapp, fake_hardware):
    window = MainWindow()
    controller = Controller(window)
    controller.take_single_exposure(345.)
    _wait_for_worker(controller)
    assert controller._standalone_sequence is not None, 'setup: a standalone exposure should exist'

    panel = window.control_panel
    panel.tabs.setCurrentWidget(panel.grid_tab)
    panel.grid_start_spin.setValue(340)
    panel.grid_step_spin.setValue(5)
    panel.grid_nstep_spin.setValue(5)
    controller.start_sequence()

    assert controller._standalone_sequence is None, \
        'starting a new sequence should discard standalone single-exposure state'
    _wait_for_worker(controller)
