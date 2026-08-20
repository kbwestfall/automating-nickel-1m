"""Tests for :mod:`gui.views.focus_control_panel`."""
from pathlib import Path

import focus
from gui.qt import QtWidgets
from gui.views.focus_control_panel import FocusControlPanel


def _step_result(focus_sweep, index=0):
    seq = focus.ArchiveFocusSequence(focus_sweep['focus_values'], focus_sweep['files'])
    return list(seq.step(method='brightest'))[index]


def test_only_archive_sequence_type_is_enabled(qapp):
    panel = FocusControlPanel()
    assert panel.archive_radio.isChecked(), 'Archive should be the default sequence type'
    assert panel.archive_radio.isEnabled(), 'Archive should be selectable in this phase'
    assert not panel.grid_radio.isEnabled(), 'Grid requires live ktl, not available yet'
    assert not panel.automated_radio.isEnabled(), 'Automated requires live ktl, not available yet'


def test_exposure_settings_are_disabled(qapp):
    panel = FocusControlPanel()
    assert not panel.exptime_spin.isEnabled(), 'exposure time is not applicable in archive mode'
    assert not panel.speed_combo.isEnabled(), 'speed is not applicable in archive mode'
    assert not panel.binning_combo.isEnabled(), 'binning is not applicable in archive mode'


def test_get_archive_config_reflects_form_fields(qapp):
    panel = FocusControlPanel()
    panel.datadir_edit.setText('/tmp/some/dir')
    panel.prefix_edit.setText('x')
    panel.suffix_edit.setText('.fit')
    panel.obsnum_spin.setValue(1234)
    panel.start_spin.setValue(300.)
    panel.step_spin.setValue(2.5)
    panel.nstep_spin.setValue(7)

    config = panel.get_archive_config()

    assert config == {
        'datadir': Path('/tmp/some/dir'),
        'prefix': 'x',
        'suffix': '.fit',
        'obsnum': 1234,
        'start': 300.,
        'step': 2.5,
        'nstep': 7,
    }, 'get_archive_config() should reflect exactly what the form fields hold'


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

    panel.set_running(True)
    assert not panel.start_button.isEnabled(), 'Start should be disabled while running'
    assert panel.stop_button.isEnabled(), 'Stop should be enabled while running'
    assert not panel.datadir_edit.isEnabled(), 'config fields should be locked while running'

    panel.set_running(False)
    assert panel.start_button.isEnabled(), 'Start should be re-enabled once stopped'
    assert not panel.stop_button.isEnabled(), 'Stop should be disabled once stopped'
    assert panel.datadir_edit.isEnabled(), 'config fields should unlock once stopped'


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


def test_show_best_focus_updates_result_but_not_move_button(qapp):
    # Per this phase's scope, "Move to Best Focus" stays disabled
    # regardless of whether a result exists -- there's no hardware to
    # move yet. Only Phase 3 should change that.
    panel = FocusControlPanel()
    panel.show_best_focus(356.1, 3.3)

    assert '356.1' in panel.result_label.text(), 'result label should show the best focus'
    assert '3.3' in panel.result_label.text(), 'result label should show the expected FWHM'
    assert not panel.move_to_best_focus_button.isEnabled(), \
        'move-to-best-focus stays disabled in this phase regardless of a result existing'


def test_show_failure_updates_status_and_log(qapp):
    panel = FocusControlPanel()
    panel.show_failure('Stopped early -- not enough points for a focus fit')

    assert panel.status_label.text() == 'Stopped early -- not enough points for a focus fit'
    assert 'ERROR' in panel.log_widget.toPlainText(), 'failures should be logged too'


def test_reset_clears_status_log_and_result(qapp, focus_sweep):
    result = _step_result(focus_sweep)
    panel = FocusControlPanel()
    panel.update_step(result)
    panel.show_best_focus(350., 3.0)
    panel.show_failure('oops')

    panel.reset()

    assert panel.step_label.text() == 'Step: —'
    assert panel.result_label.text() == 'Best focus: —'
    assert panel.status_label.text() == ''
    assert panel.log_widget.toPlainText() == ''


def test_move_to_best_focus_confirmation_flow(qapp, monkeypatch):
    panel = FocusControlPanel()
    received = []
    panel.moveToBestFocusRequested.connect(received.append)

    # No result yet: nothing should happen, regardless of the dialog.
    panel._on_move_to_best_focus_clicked()
    assert received == [], 'should be a no-op with no best-focus result yet'

    panel.show_best_focus(356.1, 3.3)

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
