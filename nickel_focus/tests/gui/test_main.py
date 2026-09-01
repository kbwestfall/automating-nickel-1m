"""Smoke tests for :mod:`nickel_focus.gui.launcher`'s scaffolding."""
import os

import pytest

from nickel_focus.gui import launcher
from nickel_focus.gui.qt import QtCore


def test_build_window_opens_without_error(qapp):
    window = launcher.build_window()
    window.show()
    qapp.processEvents()

    assert window.windowTitle() == 'Nickel Focus GUI', 'window title not set as expected'

    window.close()


def test_build_window_fits_within_the_screen(qapp):
    # A fixed 1200x600 window can be larger than the available screen
    # (excluding docks/taskbars), pushing its resize handles off screen
    # and leaving the user stuck -- `_size_to_screen` should never
    # request a size larger than what's actually available. But no
    # request can shrink the window below what its content genuinely
    # requires (`minimumSizeHint`) -- e.g. the offscreen test platform's
    # synthetic screen is a mere 800px wide, narrower than this app's
    # real, content-driven minimum width -- so the bound checked here is
    # whichever of the two is larger, rather than a hard "always fits
    # this one screen" expectation with no real content behind it. A
    # regression in the clamping logic itself (e.g. `_size_to_screen`
    # never being called) would still show up here, since the window
    # would then exceed *both* bounds, not just the screen-specific one.
    window = launcher.build_window()
    window.show()
    qapp.processEvents()

    available = window.screen().availableGeometry()
    assert window.width() <= max(available.width(), window.minimumSizeHint().width()), \
        'window should not be wider than the screen, unless its content cannot shrink further'
    assert window.height() <= max(available.height(), window.minimumSizeHint().height()), \
        'window should not be taller than the screen, unless its content cannot shrink further'

    window.close()


def test_control_panel_sits_directly_in_the_splitter(qapp):
    # Each tab page scrolls internally instead (see
    # `FocusControlPanel._tab_widget`) -- wrapping the whole panel in
    # another, outer scroll area on top of that would be redundant and
    # would let the tab bar itself scroll out of view along with a tall
    # page's content, which is exactly what an outer scroll area used to
    # do here.
    window = launcher.build_window()
    window.show()
    qapp.processEvents()

    right = window.centralWidget().widget(1)
    assert right.widget(1) is window.control_panel, \
        'the control panel should sit directly in the splitter, not wrapped in a scroll area'

    window.close()


def test_control_panel_does_not_snap_to_an_unrelated_width_when_squeezed(qapp):
    # QScrollArea.setWidgetResizable(True) has a real bug: if its
    # viewport is ever squeezed narrower than the scrolled widget's own
    # minimum width, the widget doesn't get a horizontal scrollbar as
    # you'd expect -- it snaps to some unrelated, much larger width (seen
    # in practice: a ~340px-minimum panel jumping to ~640px) and gets
    # stuck there. That risk no longer applies at this level now that the
    # control panel isn't wrapped in an outer scroll area (see
    # test_control_panel_sits_directly_in_the_splitter) -- this confirms
    # squeezing it instead just shrinks it toward its own minimum.
    window = launcher.build_window()
    window.show()
    qapp.processEvents()

    central = window.centralWidget()
    central.setSizes([central.width() - 50, 50])
    qapp.processEvents()

    panel = window.control_panel
    assert panel.width() < panel.sizeHint().width() / 2, \
        'squeezing below its sizeHint should shrink the panel toward its floor, not snap wider'

    window.close()


def test_control_panel_gets_a_typical_tab_sized_initial_share(qapp):
    # `right`'s very first layout hands the control panel just enough
    # height for a typical (non-Help) tab (see
    # `FocusControlPanel.preferred_height_excluding_help`) and gives the
    # curve panel everything else, rather than the panel claiming most of
    # the window the way an unconfigured splitter's default share would.
    window = launcher.build_window()
    window.show()
    qapp.processEvents()

    right = window.centralWidget().widget(1)
    panel = window.control_panel
    assert right.sizes()[1] == pytest.approx(panel.preferred_height_excluding_help(), abs=5), \
        "the control panel's initial height share should match its preferred height"
    for tab_name in ('slew_tab', 'single_tab', 'grid_tab', 'auto_tab', 'replay_tab',
                      'log_tab', 'options_tab'):
        tab = getattr(panel, tab_name)
        assert right.sizes()[1] >= tab.widget().minimumSizeHint().height(), \
            f'{tab_name} should fit within the initial share without scrolling'

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
