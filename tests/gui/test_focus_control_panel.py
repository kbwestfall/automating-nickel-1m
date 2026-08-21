"""Tests for :mod:`gui.views.focus_control_panel`."""
from pathlib import Path

import focus
from gui.qt import QtWidgets
from gui.views.focus_control_panel import FocusControlPanel


def _step_result(focus_sweep, index=0):
    seq = focus.ArchiveFocusSequence(focus_sweep['focus_values'], focus_sweep['files'])
    return list(seq.step(method='brightest'))[index]


def test_tabs_exist_in_the_requested_order(qapp):
    panel = FocusControlPanel()
    titles = [panel.tabs.tabText(i) for i in range(panel.tabs.count())]
    assert titles == ['Single', 'Grid', 'Auto', 'Replay', 'Log', 'Help']


def test_number_entry_boxes_have_no_increment_buttons(qapp):
    panel = FocusControlPanel()
    spin_boxes = (
        panel.single_focus_spin, panel.single_exptime_spin,
        panel.grid_start_spin, panel.grid_step_spin, panel.grid_nstep_spin,
        panel.grid_exptime_spin,
        panel.auto_start_spin, panel.auto_step_spin, panel.auto_maxsteps_spin,
        panel.auto_exptime_spin,
        panel.replay_obsnum_spin, panel.replay_start_spin, panel.replay_step_spin,
        panel.replay_nstep_spin,
    )
    for spin_box in spin_boxes:
        assert spin_box.buttonSymbols() == QtWidgets.QAbstractSpinBox.ButtonSymbols.NoButtons, \
            f'{spin_box} should not show increment/decrement buttons'


def test_focus_value_entries_are_integers(qapp):
    # A focus value/step size is a whole-unit telescope position, even
    # though a best-fit focus can land on a fractional value (that gets
    # rounded before it's ever offered as the Single tab's default --
    # see the show_best_focus tests below).
    panel = FocusControlPanel()
    integer_entries = (
        panel.single_focus_spin,
        panel.grid_start_spin, panel.grid_step_spin,
        panel.auto_start_spin, panel.auto_step_spin,
        panel.replay_start_spin, panel.replay_step_spin,
    )
    for spin_box in integer_entries:
        assert isinstance(spin_box, QtWidgets.QSpinBox), f'{spin_box} should be an integer entry'


def test_get_sequence_type(qapp):
    panel = FocusControlPanel()

    panel.tabs.setCurrentWidget(panel.single_tab)
    assert panel.get_sequence_type() is None
    panel.tabs.setCurrentWidget(panel.grid_tab)
    assert panel.get_sequence_type() == 'grid'
    panel.tabs.setCurrentWidget(panel.auto_tab)
    assert panel.get_sequence_type() == 'automated'
    panel.tabs.setCurrentWidget(panel.replay_tab)
    assert panel.get_sequence_type() == 'archive'
    panel.tabs.setCurrentWidget(panel.log_tab)
    assert panel.get_sequence_type() is None
    panel.tabs.setCurrentWidget(panel.help_tab)
    assert panel.get_sequence_type() is None


def test_get_sequence_config_for_grid(qapp):
    panel = FocusControlPanel()
    panel.tabs.setCurrentWidget(panel.grid_tab)
    panel.grid_start_spin.setValue(300)
    panel.grid_step_spin.setValue(2)
    panel.grid_nstep_spin.setValue(7)

    assert panel.get_sequence_config() == {'start': 300, 'step': 2, 'nstep': 7}


def test_get_sequence_config_for_auto(qapp):
    panel = FocusControlPanel()
    panel.tabs.setCurrentWidget(panel.auto_tab)
    panel.auto_start_spin.setValue(310)
    panel.auto_step_spin.setValue(3)
    panel.auto_maxsteps_spin.setValue(9)

    assert panel.get_sequence_config() == {'start': 310, 'step': 3, 'maxsteps': 9}


def test_get_sequence_config_for_replay(qapp):
    panel = FocusControlPanel()
    panel.tabs.setCurrentWidget(panel.replay_tab)
    panel.replay_datadir_edit.setText('/tmp/some/dir')
    panel.replay_prefix_edit.setText('x')
    panel.replay_suffix_edit.setText('.fit')
    panel.replay_obsnum_spin.setValue(1234)
    panel.replay_start_spin.setValue(300)
    panel.replay_step_spin.setValue(2)
    panel.replay_nstep_spin.setValue(7)

    assert panel.get_sequence_config() == {
        'datadir': Path('/tmp/some/dir'),
        'prefix': 'x',
        'suffix': '.fit',
        'obsnum': 1234,
        'start': 300,
        'step': 2,
        'nstep': 7,
    }


def test_get_sequence_config_is_empty_outside_grid_auto_replay(qapp):
    panel = FocusControlPanel()
    for tab in (panel.single_tab, panel.log_tab, panel.help_tab):
        panel.tabs.setCurrentWidget(tab)
        assert panel.get_sequence_config() == {}


def test_get_exposure_config_for_single_grid_auto(qapp):
    panel = FocusControlPanel()

    panel.tabs.setCurrentWidget(panel.single_tab)
    panel.single_exptime_spin.setValue(1.5)
    panel.single_speed_combo.setCurrentText('Fast')
    panel.single_binning_combo.setCurrentText('2,2')
    assert panel.get_exposure_config() == {'exptime': 1.5, 'speed': 'Fast', 'binning': '2,2'}

    panel.tabs.setCurrentWidget(panel.grid_tab)
    panel.grid_exptime_spin.setValue(2.5)
    panel.grid_speed_combo.setCurrentText('Slow')
    panel.grid_binning_combo.setCurrentText('1,1')
    assert panel.get_exposure_config() == {'exptime': 2.5, 'speed': 'Slow', 'binning': '1,1'}

    panel.tabs.setCurrentWidget(panel.auto_tab)
    panel.auto_exptime_spin.setValue(3.5)
    panel.auto_speed_combo.setCurrentText('Fast')
    panel.auto_binning_combo.setCurrentText('4,4')
    assert panel.get_exposure_config() == {'exptime': 3.5, 'speed': 'Fast', 'binning': '4,4'}


def test_get_exposure_config_is_empty_for_replay_log_help(qapp):
    panel = FocusControlPanel()
    for tab in (panel.replay_tab, panel.log_tab, panel.help_tab):
        panel.tabs.setCurrentWidget(tab)
        assert panel.get_exposure_config() == {}


def test_start_button_enabled_only_on_actionable_tabs(qapp):
    panel = FocusControlPanel()
    for tab in (panel.single_tab, panel.grid_tab, panel.auto_tab, panel.replay_tab):
        panel.tabs.setCurrentWidget(tab)
        assert panel.start_button.isEnabled(), f'{tab} should allow Start'
    for tab in (panel.log_tab, panel.help_tab):
        panel.tabs.setCurrentWidget(tab)
        assert not panel.start_button.isEnabled(), f'{tab} has nothing to start'


def test_start_click_on_single_tab_emits_take_single_exposure_requested(qapp):
    panel = FocusControlPanel()
    panel.tabs.setCurrentWidget(panel.single_tab)
    panel.single_focus_spin.setValue(356)
    starts, singles = [], []
    panel.startRequested.connect(lambda: starts.append(True))
    panel.takeSingleExposureRequested.connect(singles.append)

    panel.start_button.click()

    assert starts == [], 'Single tab should not emit startRequested'
    assert singles == [356], 'Single tab should emit the configured focus value'


def test_start_click_on_grid_tab_emits_start_requested(qapp):
    panel = FocusControlPanel()
    panel.tabs.setCurrentWidget(panel.grid_tab)
    starts, singles = [], []
    panel.startRequested.connect(lambda: starts.append(True))
    panel.takeSingleExposureRequested.connect(singles.append)

    panel.start_button.click()

    assert starts == [True], 'Grid tab should emit startRequested'
    assert singles == [], 'Grid tab should not emit takeSingleExposureRequested'


def test_stop_signal(qapp):
    panel = FocusControlPanel()
    stops = []
    panel.stopRequested.connect(lambda: stops.append(True))

    # Stop starts disabled (nothing is running yet); a disabled button
    # doesn't emit clicked() at all, so put the panel in a running state
    # first, matching how the Controller will actually use it.
    panel.set_running(True)
    panel.stop_button.click()
    assert stops == [True], 'clicking Stop should emit stopRequested'


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


def test_set_running_locks_config_widgets_and_method(qapp):
    panel = FocusControlPanel()
    panel.tabs.setCurrentWidget(panel.grid_tab)

    panel.set_running(True)
    assert not panel.start_button.isEnabled(), 'Start should be disabled while running'
    assert panel.stop_button.isEnabled(), 'Stop should be enabled while running'
    assert not panel.grid_start_spin.isEnabled(), 'config fields should be locked while running'
    assert not panel.method_combo.isEnabled(), 'method should be locked while running'

    panel.set_running(False)
    assert panel.start_button.isEnabled(), 'Start should be re-enabled once stopped (Grid is active)'
    assert not panel.stop_button.isEnabled(), 'Stop should be disabled once stopped'
    assert panel.grid_start_spin.isEnabled(), 'config fields should unlock once stopped'
    assert panel.method_combo.isEnabled(), 'method should unlock once stopped'


def test_set_running_restores_start_button_per_current_tab(qapp):
    panel = FocusControlPanel()
    panel.tabs.setCurrentWidget(panel.log_tab)

    panel.set_running(True)
    panel.set_running(False)

    assert not panel.start_button.isEnabled(), \
        'Start should reflect the Log tab (nothing to start), not just re-enable blindly'


def test_tab_bar_stays_enabled_while_running(qapp):
    # The Log tab must stay reachable while a sequence runs.
    panel = FocusControlPanel()
    panel.set_running(True)
    assert panel.tabs.isEnabled(), 'the tab bar itself should never be disabled'


def test_switching_tabs_while_running_does_not_enable_start(qapp):
    panel = FocusControlPanel()
    panel.tabs.setCurrentWidget(panel.grid_tab)
    panel.set_running(True)

    panel.tabs.setCurrentWidget(panel.single_tab)

    assert not panel.start_button.isEnabled(), \
        'Start should stay disabled while running, regardless of tab switches'


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
    assert f'{result.focus_value:.0f}' in panel.step_label.text(), 'step label should show focus'
    assert panel.log_widget.toPlainText().strip(), 'the step should also be appended to the log'


def test_show_best_focus_updates_log_and_single_tab_default(qapp):
    panel = FocusControlPanel()

    panel.show_best_focus(356.6, 3.3)

    assert '356.6' in panel.status_label.text(), 'status should show the best focus'
    assert '3.3' in panel.status_label.text(), 'status should show the expected FWHM'
    assert 'ERROR' not in panel.log_widget.toPlainText()
    assert panel.single_focus_spin.value() == 357, \
        'the Single tab default should be the best focus, rounded to the nearest integer'


def test_show_single_exposure_result_updates_status_and_log(qapp, focus_sweep):
    result = _step_result(focus_sweep)
    panel = FocusControlPanel()

    panel.show_single_exposure_result(result)

    assert f'{result.focus_value:.0f}' in panel.status_label.text(), \
        'status should report the exposure focus value'
    assert f'{result.fwhm:.2f}' in panel.status_label.text(), \
        'status should report the measured FWHM'
    assert panel.log_widget.toPlainText().strip(), 'the result should also be logged'


def test_show_failure_updates_status_and_log(qapp):
    panel = FocusControlPanel()
    panel.show_failure('Stopped early -- not enough points for a focus fit')

    assert panel.status_label.text() == 'Stopped early -- not enough points for a focus fit'
    assert 'ERROR' in panel.log_widget.toPlainText(), 'failures should be logged too'


def test_reset_clears_status_and_log(qapp, focus_sweep):
    result = _step_result(focus_sweep)
    panel = FocusControlPanel()
    panel.update_step(result)
    panel.show_failure('oops')

    panel.reset()

    assert panel.step_label.text() == 'Step: —'
    assert panel.status_label.text() == ''
    assert panel.log_widget.toPlainText() == ''


def test_reset_does_not_touch_the_single_tab_default(qapp):
    panel = FocusControlPanel()
    panel.show_best_focus(356.6, 3.3)

    panel.reset()

    assert panel.single_focus_spin.value() == 357, \
        "reset() clears status/log, but shouldn't discard the Single tab's best-focus default"


def test_browse_button_updates_datadir(qapp, monkeypatch):
    panel = FocusControlPanel()

    monkeypatch.setattr(QtWidgets.QFileDialog, 'getExistingDirectory',
                         lambda *a, **k: '/some/chosen/dir')
    panel._on_browse_clicked()
    assert panel.replay_datadir_edit.text() == '/some/chosen/dir'

    # Canceling the dialog (empty string) should leave the field alone.
    monkeypatch.setattr(QtWidgets.QFileDialog, 'getExistingDirectory', lambda *a, **k: '')
    panel._on_browse_clicked()
    assert panel.replay_datadir_edit.text() == '/some/chosen/dir', 'canceling should not clear the field'


def test_help_tab_has_reminder_text(qapp):
    panel = FocusControlPanel()
    labels = panel.help_tab.findChildren(QtWidgets.QLabel)
    assert labels, 'the Help tab should contain some descriptive text'
    text = ' '.join(label.text() for label in labels)
    for keyword in ('Single', 'Grid', 'Auto', 'Replay', 'Log'):
        assert keyword in text, f'Help text should mention {keyword}'
