"""Tests for :mod:`gui.views.focus_control_panel`."""
from pathlib import Path

import focus
from gui.qt import QtWidgets
from gui.views.focus_control_panel import FocusControlPanel


def _step_result(focus_sweep, index=0):
    seq = focus.ArchiveFocusSequence(focus_sweep['focus_values'], focus_sweep['files'])
    return list(seq.step(method='brightest'))[index]


def test_all_sequence_types_are_selectable(qapp):
    panel = FocusControlPanel()
    assert panel.archive_radio.isChecked(), 'Archive should be the default sequence type'
    for radio in (panel.archive_radio, panel.grid_radio, panel.automated_radio):
        assert radio.isEnabled(), 'all three sequence types should be selectable'


def test_get_sequence_type(qapp):
    panel = FocusControlPanel()
    assert panel.get_sequence_type() == 'archive'

    panel.grid_radio.setChecked(True)
    assert panel.get_sequence_type() == 'grid'

    panel.automated_radio.setChecked(True)
    assert panel.get_sequence_type() == 'automated'


def test_archive_fields_enabled_only_for_archive(qapp):
    panel = FocusControlPanel()
    archive_fields = (panel.datadir_edit, panel.browse_button, panel.prefix_edit,
                      panel.suffix_edit, panel.obsnum_spin)

    assert all(w.isEnabled() for w in archive_fields), 'Archive fields start enabled (default type)'
    assert panel.nstep_spin.isEnabled(), 'Archive uses number-of-steps'
    assert not panel.maxsteps_spin.isEnabled(), 'Archive does not use max-steps'
    assert not panel.exptime_spin.isEnabled(), 'exposure settings are not applicable in archive mode'

    panel.grid_radio.setChecked(True)
    assert not any(w.isEnabled() for w in archive_fields), 'Grid has no archive fields to fill in'
    assert panel.nstep_spin.isEnabled(), 'Grid uses number-of-steps'
    assert not panel.maxsteps_spin.isEnabled(), 'Grid does not use max-steps'
    assert panel.exptime_spin.isEnabled(), 'Grid takes real exposures, so settings apply'

    panel.automated_radio.setChecked(True)
    assert not any(w.isEnabled() for w in archive_fields), 'Automated has no archive fields either'
    assert not panel.nstep_spin.isEnabled(), 'Automated does not use number-of-steps'
    assert panel.maxsteps_spin.isEnabled(), 'Automated uses max-steps instead'
    assert panel.exptime_spin.isEnabled(), 'Automated takes real exposures, so settings apply'


def test_get_sequence_config_reflects_form_fields(qapp):
    panel = FocusControlPanel()
    panel.datadir_edit.setText('/tmp/some/dir')
    panel.prefix_edit.setText('x')
    panel.suffix_edit.setText('.fit')
    panel.obsnum_spin.setValue(1234)
    panel.start_spin.setValue(300.)
    panel.step_spin.setValue(2.5)
    panel.nstep_spin.setValue(7)
    panel.maxsteps_spin.setValue(9)

    config = panel.get_sequence_config()

    assert config == {
        'datadir': Path('/tmp/some/dir'),
        'prefix': 'x',
        'suffix': '.fit',
        'obsnum': 1234,
        'start': 300.,
        'step': 2.5,
        'nstep': 7,
        'maxsteps': 9,
    }, 'get_sequence_config() should reflect exactly what the form fields hold'


def test_get_exposure_config_reflects_form_fields(qapp):
    panel = FocusControlPanel()
    panel.exptime_spin.setValue(12.5)
    panel.speed_combo.setCurrentText('Fast')
    panel.binning_combo.setCurrentText('2,2')

    assert panel.get_exposure_config() == {'exptime': 12.5, 'speed': 'Fast', 'binning': '2,2'}


def test_method_combo_and_signal(qapp):
    panel = FocusControlPanel()
    received = []
    panel.methodChanged.connect(received.append)

    assert panel.get_selected_method() == 'brightest', 'Brightest should be the default'

    panel.method_combo.setCurrentText('Weighted')
    assert panel.get_selected_method() == 'weighted'
    assert received == ['weighted'], 'methodChanged should emit the lowercase method name'


def test_set_method_updates_display_for_string_and_coordinates(qapp):
    panel = FocusControlPanel()

    panel.set_method('brightest')
    assert panel.method_label.text() == 'Current method: Brightest'

    panel.set_method((123.4, 567.8))
    assert panel.method_label.text() == 'Current method: Selected source (123.4, 567.8)'


def test_start_and_stop_signals(qapp):
    panel = FocusControlPanel()
    starts = []
    stops = []
    panel.startRequested.connect(lambda: starts.append(True))
    panel.stopRequested.connect(lambda: stops.append(True))

    panel.start_button.click()
    assert starts == [True], 'clicking Start should emit startRequested'

    # Stop starts disabled (nothing is running yet); a disabled button
    # doesn't emit clicked() at all, so put the panel in a running state
    # first, matching how the Controller will actually use it.
    panel.set_running(True)
    panel.stop_button.click()
    assert stops == [True], 'clicking Stop should emit stopRequested'


def test_set_running_toggles_widget_states(qapp):
    panel = FocusControlPanel()
    panel.show_best_focus(350., 3.0, can_move=True)

    panel.set_running(True)
    assert not panel.start_button.isEnabled(), 'Start should be disabled while running'
    assert panel.stop_button.isEnabled(), 'Stop should be enabled while running'
    assert not panel.datadir_edit.isEnabled(), 'config fields should be locked while running'
    assert not panel.move_to_best_focus_button.isEnabled(), \
        'move-to-best-focus should be locked out while anything else is running'

    panel.set_running(False)
    assert panel.start_button.isEnabled(), 'Start should be re-enabled once stopped'
    assert not panel.stop_button.isEnabled(), 'Stop should be disabled once stopped'
    assert panel.datadir_edit.isEnabled(), 'config fields should unlock once stopped'
    assert not panel.move_to_best_focus_button.isEnabled(), \
        'move-to-best-focus should stay disabled until the next result explicitly re-enables it'


def test_set_running_restores_the_correct_per_type_state(qapp):
    # set_running(False) must not just re-enable everything -- it should
    # restore whatever's correct for the currently-selected sequence
    # type, e.g. Automated has no archive fields and no nstep, even after
    # a run finishes.
    panel = FocusControlPanel()
    panel.automated_radio.setChecked(True)

    panel.set_running(True)
    panel.set_running(False)

    assert not panel.datadir_edit.isEnabled(), 'Automated has no archive fields to restore'
    assert not panel.nstep_spin.isEnabled(), 'Automated does not use number-of-steps'
    assert panel.maxsteps_spin.isEnabled(), 'Automated should have max-steps re-enabled'
    assert panel.exptime_spin.isEnabled(), 'Automated should have exposure settings re-enabled'


def test_set_stopping_disables_stop_and_shows_status(qapp):
    panel = FocusControlPanel()
    panel.set_running(True)

    panel.set_stopping()

    assert not panel.stop_button.isEnabled(), 'Stop should disable itself once clicked'
    assert 'Stopping' in panel.status_label.text(), 'status should explain the wait'


def test_update_step_updates_label_and_log(qapp, focus_sweep):
    result = _step_result(focus_sweep, index=1)
    panel = FocusControlPanel()

    panel.update_step(result, total_expected=5)

    assert 'Step 2/5' in panel.step_label.text(), 'step label should show 1-based step/total'
    assert f'{result.focus_value:.1f}' in panel.step_label.text(), 'step label should show focus'
    assert panel.log_widget.toPlainText().strip(), 'the step should also be appended to the log'


def test_show_best_focus_updates_result_and_gates_move_button(qapp):
    panel = FocusControlPanel()

    panel.show_best_focus(356.1, 3.3, can_move=False)
    assert '356.1' in panel.result_label.text(), 'result label should show the best focus'
    assert '3.3' in panel.result_label.text(), 'result label should show the expected FWHM'
    assert not panel.move_to_best_focus_button.isEnabled(), \
        'move-to-best-focus should stay disabled when the sequence has no hardware to move'

    panel.show_best_focus(356.1, 3.3, can_move=True)
    assert panel.move_to_best_focus_button.isEnabled(), \
        'move-to-best-focus should enable once a live sequence finishes'


def test_show_confirmation_updates_status_and_log(qapp, focus_sweep):
    result = _step_result(focus_sweep)
    panel = FocusControlPanel()

    panel.show_confirmation(result)

    assert f'{result.focus_value:.1f}' in panel.status_label.text(), \
        'status should report the confirmed focus value'
    assert f'{result.fwhm:.2f}' in panel.status_label.text(), \
        'status should report the measured FWHM'
    assert panel.log_widget.toPlainText().strip(), 'the confirmation should also be logged'


def test_show_failure_updates_status_and_log(qapp):
    panel = FocusControlPanel()
    panel.show_failure('Stopped early -- not enough points for a focus fit')

    assert panel.status_label.text() == 'Stopped early -- not enough points for a focus fit'
    assert 'ERROR' in panel.log_widget.toPlainText(), 'failures should be logged too'


def test_reset_clears_status_log_and_result(qapp, focus_sweep):
    result = _step_result(focus_sweep)
    panel = FocusControlPanel()
    panel.update_step(result)
    panel.show_best_focus(350., 3.0, can_move=True)
    panel.show_failure('oops')

    panel.reset()

    assert panel.step_label.text() == 'Step: —'
    assert panel.result_label.text() == 'Best focus: —'
    assert panel.status_label.text() == ''
    assert panel.log_widget.toPlainText() == ''
    assert not panel.move_to_best_focus_button.isEnabled(), \
        'reset should disable move-to-best-focus until the next result arrives'


def test_move_to_best_focus_confirmation_flow(qapp, monkeypatch):
    panel = FocusControlPanel()
    received = []
    panel.moveToBestFocusRequested.connect(received.append)

    # No result yet: nothing should happen, regardless of the dialog.
    panel._on_move_to_best_focus_clicked()
    assert received == [], 'should be a no-op with no best-focus result yet'

    panel.show_best_focus(356.1, 3.3, can_move=True)

    monkeypatch.setattr(
        QtWidgets.QMessageBox, 'question',
        lambda *a, **k: QtWidgets.QMessageBox.StandardButton.No)
    panel._on_move_to_best_focus_clicked()
    assert received == [], 'declining the confirmation should not emit the signal'

    monkeypatch.setattr(
        QtWidgets.QMessageBox, 'question',
        lambda *a, **k: QtWidgets.QMessageBox.StandardButton.Yes)
    panel._on_move_to_best_focus_clicked()
    assert received == [356.1], 'confirming should emit the target focus value'


def test_browse_button_updates_datadir(qapp, monkeypatch):
    panel = FocusControlPanel()

    monkeypatch.setattr(QtWidgets.QFileDialog, 'getExistingDirectory',
                         lambda *a, **k: '/some/chosen/dir')
    panel._on_browse_clicked()
    assert panel.datadir_edit.text() == '/some/chosen/dir'

    # Canceling the dialog (empty string) should leave the field alone.
    monkeypatch.setattr(QtWidgets.QFileDialog, 'getExistingDirectory', lambda *a, **k: '')
    panel._on_browse_clicked()
    assert panel.datadir_edit.text() == '/some/chosen/dir', 'canceling should not clear the field'
