"""Smoke tests for :mod:`nickel_focus.gui.launcher`'s scaffolding."""
import os

from nickel_focus.gui import launcher
from nickel_focus.gui.qt import QtCore, QtWidgets


def test_build_window_opens_without_error(qapp):
    window = launcher.build_window()
    window.show()
    qapp.processEvents()

    assert window.windowTitle() == 'Nickel Focus GUI', 'window title not set as expected'

    window.close()


def test_build_window_fits_within_the_screen(qapp):
    # A fixed 1200x800 window can be larger than the available screen
    # (excluding docks/taskbars), pushing its resize handles off screen
    # and leaving the user stuck -- the window should never open larger
    # than what's actually available.
    window = launcher.build_window()
    window.show()
    qapp.processEvents()

    available = window.screen().availableGeometry()
    assert window.width() <= available.width(), 'window should not be wider than the screen'
    assert window.height() <= available.height(), 'window should not be taller than the screen'

    window.close()


def test_control_scroll_area_has_a_floor_at_the_panels_minimum_width(qapp):
    # QScrollArea.setWidgetResizable(True) has a real bug: if its
    # viewport ever gets squeezed narrower than the scrolled widget's own
    # minimum width, the widget doesn't get a horizontal scrollbar as
    # you'd expect -- it snaps to some unrelated, much larger width (seen
    # in practice: a ~340px-minimum panel jumping to ~640px) and gets
    # stuck there, buttons and all. The mitigation is giving the scroll
    # area itself an explicit minimum width -- a hard floor, unlike a
    # sizeHint -- that the splitter can't ask it to go below, so that
    # code path never triggers. This confirms that floor is actually in
    # place and covers the panel's current minimum.
    window = launcher.build_window()
    window.show()
    qapp.processEvents()

    control_scroll = window.centralWidget().widget(1).widget(1)
    assert isinstance(control_scroll, QtWidgets.QScrollArea), \
        'setup: expected the control panel to be wrapped in a QScrollArea'
    assert control_scroll.minimumWidth() >= window.control_panel.minimumSizeHint().width(), \
        "the scroll area's minimum width must never be less than the panel's own minimum"

    window.close()


def test_settings_default_unchecked_and_nothing_is_ever_written(qapp):
    # The `isolate_qsettings` autouse fixture (see conftest.py) points
    # MainWindow._settings() at a per-test ini file, but doesn't create
    # it -- so this also verifies that merely opening/closing a window
    # without opting in creates no settings file at all, not just that
    # values don't come back.
    window = launcher.build_window()
    settings_path = window._settings().fileName()
    assert not window.control_panel.remember_settings_checkbox.isChecked(), \
        'persistence should default to opted out'

    window.control_panel.grid_start_spin.setValue(311)
    window.close()

    assert not os.path.exists(settings_path), \
        'closing without opting in should never create a settings file'

    window2 = launcher.build_window()
    assert not window2.control_panel.remember_settings_checkbox.isChecked(), \
        'still opted out -- there is nothing to have detected'
    assert window2.control_panel.grid_start_spin.value() != 311, \
        'nothing should have been restored'
    window2.close()


def test_settings_are_saved_on_close_and_restored_on_open_once_opted_in(qapp):
    # Both windows here share the same isolated settings file, since
    # `isolate_qsettings` fixes it for the life of one test.
    window = launcher.build_window()
    window.control_panel.remember_settings_checkbox.setChecked(True)
    window.control_panel.grid_start_spin.setValue(311)
    window.control_panel.tabs.setCurrentWidget(window.control_panel.replay_tab)
    window.control_panel.replay_datadir_edit.setText('/persisted/dir')
    window.image_panel.stretch_combo.setCurrentText('Min/Max')

    window.close()
    window2 = launcher.build_window()

    assert window2.control_panel.remember_settings_checkbox.isChecked(), \
        'opting in should itself be detected and restored'
    assert window2.control_panel.grid_start_spin.value() == 311, \
        'the Grid start focus should have been restored'
    assert window2.control_panel.replay_datadir_edit.text() == '/persisted/dir', \
        'the Replay data directory should have been restored'
    assert window2.image_panel._stretch_name == 'Min/Max', \
        'the stretch preference should have been restored'

    window2.close()


def test_settings_do_not_restore_a_stale_single_tab_focus_value(qapp):
    window = launcher.build_window()
    window.control_panel.remember_settings_checkbox.setChecked(True)
    window.control_panel.single_focus_spin.setValue(499)
    window.close()

    window2 = launcher.build_window()

    assert window2.control_panel.single_focus_spin.value() != 499, \
        "the Single tab's focus default should not be restored from a previous session"
    window2.close()


def test_unchecking_settings_erases_previously_saved_configuration(qapp):
    window = launcher.build_window()
    window.control_panel.remember_settings_checkbox.setChecked(True)
    window.control_panel.grid_start_spin.setValue(311)
    window.close()

    window2 = launcher.build_window()
    assert window2.control_panel.grid_start_spin.value() == 311, 'setup: opting in should persist'
    window2.control_panel.remember_settings_checkbox.setChecked(False)
    window2.close()

    settings_path = window2._settings().fileName()
    assert os.path.exists(settings_path), \
        'the file itself may still exist (Qt does not delete an emptied ini file), but...'
    settings = QtCore.QSettings(settings_path, QtCore.QSettings.Format.IniFormat)
    assert 'control_panel' not in settings.childGroups(), \
        'its contents should be erased once the user opts back out'

    window3 = launcher.build_window()
    assert not window3.control_panel.remember_settings_checkbox.isChecked(), \
        'opting back out should itself be detected on the next launch'
    assert window3.control_panel.grid_start_spin.value() != 311, \
        'the erased configuration should not come back'
    window3.close()
