"""Tests for :mod:`gui.controller`."""
import types

import pytest

from nickel_focus import focus
from nickel_focus.gui.qt import QtCore
from nickel_focus.gui.views.main_window import MainWindow
from nickel_focus.gui.controller import Controller


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


def _wait_for_slew_worker(controller, timeout_ms=5000):
    worker = controller.slew_worker
    if worker is None:
        return
    loop = QtCore.QEventLoop()
    worker.finished.connect(loop.quit)
    QtCore.QTimer.singleShot(timeout_ms, loop.quit)
    loop.exec()


def _write_starlist(tmp_path, lines):
    """Write a small starlist file for a test to point the Slew tab's file field at."""
    path = tmp_path / 'stars.txt'
    path.write_text('\n'.join(lines) + '\n')
    return path


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


def test_reanalyze_clears_the_curve_panel_immediately_not_iteratively(qapp, focus_sweep):
    window, controller = _make_controller(focus_sweep)
    controller.start_sequence()
    _wait_for_worker(controller)
    assert len(window.curve_panel._results) == len(focus_sweep['files']), \
        'setup: the curve should be populated from the initial sequence'

    controller.reanalyze()

    # The clear must happen synchronously, before the worker thread has
    # produced any results -- not get whittled down to zero one point at
    # a time as each re-measured result streams back in.
    assert window.curve_panel._results == [], \
        'the curve should be cleared immediately when reanalysis starts'

    _wait_for_worker(controller)
    assert len(window.curve_panel._results) == len(focus_sweep['files']), \
        'the curve should be fully repopulated once reanalysis completes'


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


def test_controller_without_ktl_shows_placeholder_current_position(qapp):
    # No fake_telescope fixture here: this dev machine has no ktl, so
    # `telescope` should end up None and the live display should show a
    # clear placeholder rather than stale/bogus coordinates.
    window = MainWindow()
    controller = Controller(window)

    assert controller.telescope is None, 'setup: no ktl connection is available'
    assert window.control_panel.slew_current_ra_label.text() == '—'
    assert window.control_panel.slew_current_dec_label.text() == '—'


def test_controller_polls_current_position_from_fake_telescope(qapp, fake_telescope):
    fake_telescope.ra, fake_telescope.dec = 5.5, 30.25

    window = MainWindow()
    Controller(window)

    assert window.control_panel.slew_current_ra_label.text() == '05:30:00.00', \
        'the current-position display should be seeded from the telescope on construction'
    assert window.control_panel.slew_current_dec_label.text() == '+30:15:00.00'


def test_move_to_target_runs_against_fake_telescope(qapp, fake_telescope):
    window = MainWindow()
    controller = Controller(window)

    controller.move_to_target('05:30:00', '+20:15:00')
    _wait_for_slew_worker(controller)

    assert controller.slew_worker is None, 'the slew worker should be cleared once finished'
    assert len(fake_telescope.slew_calls) == 1, \
        'the telescope should have been commanded to slew once'
    commanded_ra, commanded_dec = fake_telescope.slew_calls[0]
    assert commanded_ra.to_string(unit='hourangle', sep=':') == '5:30:00', \
        'the telescope should be commanded to the entered target RA'
    assert 'complete' in window.control_panel.status_label.text().lower(), \
        'a successful move should be reported'


def test_move_to_target_reports_telescope_failure(qapp, fake_telescope):
    fake_telescope.tracking_on = False
    window = MainWindow()
    controller = Controller(window)

    controller.move_to_target('05:30:00', '+20:15:00')
    _wait_for_slew_worker(controller)

    assert controller.slew_worker is None, 'the slew worker should be cleared once finished'
    assert 'Tracking is disabled' in window.control_panel.status_label.text(), \
        'the telescope failure message should be reported'


def test_move_to_target_reports_a_parse_failure_for_bad_coordinates(qapp, fake_telescope):
    window = MainWindow()
    controller = Controller(window)

    controller.move_to_target('not a coordinate', '+20:15:00')

    assert controller.slew_worker is None, 'a bad target should never start a slew worker'
    assert 'Could not parse target coordinates' in window.control_panel.status_label.text()


def test_move_to_target_without_ktl_reports_clear_failure(qapp):
    window = MainWindow()
    controller = Controller(window)

    controller.move_to_target('05:30:00', '+20:15:00')

    assert controller.slew_worker is None, 'no worker should start without a ktl connection'
    assert 'no ktl connection' in window.control_panel.status_label.text().lower()


def test_move_to_target_is_a_noop_while_something_is_running(qapp, fake_telescope):
    window = MainWindow()
    controller = Controller(window)
    controller.worker = types.SimpleNamespace()  # pretend a focus sequence is running

    controller.move_to_target('05:30:00', '+20:15:00')

    assert controller.slew_worker is None, \
        'a move request should be ignored while a focus sequence is running'
    assert fake_telescope.slew_calls == [], 'the telescope should never have been commanded'


def test_start_sequence_is_a_noop_while_slewing(qapp, fake_telescope):
    window = MainWindow()
    controller = Controller(window)
    controller.slew_worker = types.SimpleNamespace()  # pretend a slew is running
    panel = window.control_panel
    panel.tabs.setCurrentWidget(panel.grid_tab)
    panel.grid_nstep_spin.setValue(5)

    controller.start_sequence()

    assert controller.sequence is None, \
        'a focus sequence should be ignored while the telescope is slewing'


def test_find_nearest_pointing_populates_target_fields(qapp, fake_telescope):
    fake_telescope.ra, fake_telescope.dec = 0.5, 29.75  # near the packaged Pointing00 entry
    window = MainWindow()
    controller = Controller(window)

    controller.find_nearest_pointing()

    assert window.control_panel.slew_target_ra_edit.text() != '', \
        'the target RA field should be populated with the found target'
    assert 'Pointing' in window.control_panel.status_label.text(), \
        'the found target name should be reported'


def test_find_nearest_focus_populates_target_fields(qapp, fake_telescope):
    fake_telescope.ra, fake_telescope.dec = 0.5, 28.1  # near the packaged Focusing00 entry
    window = MainWindow()
    controller = Controller(window)

    controller.find_nearest_focus()

    assert window.control_panel.slew_target_ra_edit.text() != '', \
        'the target RA field should be populated with the found target'
    assert 'Focusing' in window.control_panel.status_label.text(), \
        'the found target name should be reported'


def test_find_nearest_object_requires_a_file(qapp, fake_telescope):
    window = MainWindow()
    controller = Controller(window)

    with pytest.raises(ValueError, match='starlist file must be provided'):
        controller.find_nearest_object('', '')


def test_find_nearest_object_searches_the_given_file(qapp, fake_telescope, tmp_path):
    path = _write_starlist(tmp_path, [
        'StarA 01:00:00 +10:00:00 2000.0',
        'StarB 13:00:00 -20:00:00 2000.0',
    ])
    fake_telescope.ra, fake_telescope.dec = 1.0, 10.0
    window = MainWindow()
    controller = Controller(window)

    controller.find_nearest_object(str(path), '')

    assert 'StarA' in window.control_panel.status_label.text(), \
        'the nearest target in the given file should be found and reported'


def test_find_nearest_object_reports_a_missing_file(qapp, fake_telescope):
    window = MainWindow()
    controller = Controller(window)

    controller.find_nearest_object('/does/not/exist.txt', '')

    assert 'Could not find nearest target' in window.control_panel.status_label.text()


def test_find_nearest_without_ktl_reports_clear_failure(qapp):
    window = MainWindow()
    controller = Controller(window)

    controller.find_nearest_pointing()

    assert 'no ktl connection' in window.control_panel.status_label.text().lower()
