"""Tests for :mod:`gui.views.image_panel`."""
import dataclasses

import pytest
from matplotlib.backend_bases import KeyEvent, MouseEvent
from matplotlib.colors import to_rgba

import focus
from gui.views.image_panel import ImagePanel


def _step_results(focus_sweep):
    seq = focus.ArchiveFocusSequence(focus_sweep['focus_values'], focus_sweep['files'])
    return list(seq.step(method='brightest'))


def test_add_result_keeps_dropdown_sorted_by_focus_value(qapp, focus_sweep):
    results = _step_results(focus_sweep)
    # results are naturally ascending in focus (340..360); add them out of
    # order to confirm sorting doesn't depend on acquisition order.
    scrambled = [results[2], results[0], results[4], results[1], results[3]]

    panel = ImagePanel()
    for r in scrambled:
        panel.add_result(r)

    combo_focus_values = [r.focus_value for r in panel._results]
    assert combo_focus_values == sorted(r.focus_value for r in results), \
        'dropdown should stay sorted by focus value regardless of insertion order'
    assert panel.exposure_combo.count() == len(results), 'dropdown item count should match'

    # "Latest" tracks acquisition order: the last-added result (focus 355)
    # should be displayed, even though it isn't the highest focus value.
    assert panel._current is scrambled[-1], \
        'panel should auto-follow the most recently *added* result, not the sorted-last one'


def test_manual_selection_is_not_disrupted_by_new_results(qapp, focus_sweep):
    results = _step_results(focus_sweep)
    panel = ImagePanel()
    for r in results[:3]:
        panel.add_result(r)

    panel.exposure_combo.setCurrentIndex(0)
    assert panel._current is results[0], 'setup: manual selection should change the displayed result'

    panel.add_result(results[3])
    assert panel._current is results[0], \
        'adding a new result should not disrupt a manually-selected older entry'
    assert results[3] in panel._results, 'the new result should still be listed in the dropdown'


def test_m_key_emits_signal_only_when_selection_enabled(qapp, focus_sweep):
    results = _step_results(focus_sweep)
    panel = ImagePanel()
    panel.add_result(results[0])

    received = []
    panel.sourceSelected.connect(lambda x, y: received.append((x, y)))

    event = KeyEvent('key_press_event', panel.canvas, 'm')
    event.inaxes = panel.ax
    event.xdata, event.ydata = 42.0, 17.0

    panel._on_key_press(event)
    assert received == [], "'m' should be ignored while selection is disabled"

    panel.set_selection_enabled(True)
    panel._on_key_press(event)
    assert received == [(42.0, 17.0)], \
        "'m' should emit the data coordinates under the cursor once enabled"


def test_other_keys_do_not_select_a_source(qapp, focus_sweep):
    results = _step_results(focus_sweep)
    panel = ImagePanel()
    panel.add_result(results[0])
    panel.set_selection_enabled(True)

    received = []
    panel.sourceSelected.connect(lambda x, y: received.append((x, y)))

    other_key = KeyEvent('key_press_event', panel.canvas, 'a')
    other_key.inaxes = panel.ax
    other_key.xdata, other_key.ydata = 42.0, 17.0
    panel._on_key_press(other_key)
    assert received == [], "only the 'm' key should select a source"


def test_plain_click_does_not_select_a_source(qapp, focus_sweep):
    # A plain click must never select a source on its own -- on macOS,
    # the click that refocuses an inactive window is delivered to
    # whatever's underneath the cursor, so a click-to-select would
    # trigger reanalysis just from bringing the window forward.
    results = _step_results(focus_sweep)
    panel = ImagePanel()
    panel.add_result(results[0])
    panel.set_selection_enabled(True)

    received = []
    panel.sourceSelected.connect(lambda x, y: received.append((x, y)))

    # Dispatch through the canvas's own callback registry (matching how
    # an actual click arrives), rather than calling a handler directly --
    # this proves nothing is listening for clicks anymore, not just that
    # some particular method is gone.
    event = MouseEvent('button_press_event', panel.canvas, x=10, y=10, button=1)
    event.inaxes = panel.ax
    event.xdata, event.ydata = 42.0, 17.0
    panel.canvas.callbacks.process('button_press_event', event)

    assert received == [], 'a plain click should never select a source'


def test_outlier_box_color_reflects_is_outlier(qapp, focus_sweep):
    results = _step_results(focus_sweep)
    panel = ImagePanel()

    panel.show_result(dataclasses.replace(results[0], is_outlier=False))
    colors = [p.get_edgecolor() for p in panel.ax.patches]
    assert to_rgba('red') in colors, 'a non-outlier result should be boxed in red'

    panel.show_result(dataclasses.replace(results[0], is_outlier=True))
    colors = [p.get_edgecolor() for p in panel.ax.patches]
    assert to_rgba('yellow') in colors, 'an outlier result should be boxed in yellow'


def test_dropdown_selection_preserves_zoom_and_center(qapp, focus_sweep):
    results = _step_results(focus_sweep)
    panel = ImagePanel()
    for r in results:
        panel.add_result(r)

    panel.set_zoom(2.0)
    panel.h_scroll.setValue(10)
    panel.v_scroll.setValue(15)
    zoom_before = panel._zoom
    xlim_before, ylim_before = panel.ax.get_xlim(), panel.ax.get_ylim()

    panel.exposure_combo.setCurrentIndex(0)

    assert panel._current is results[0], 'setup: the drop-down selection should have changed'
    assert panel._zoom == zoom_before, \
        'switching frames via the drop-down should not change the zoom factor'
    assert panel.ax.get_xlim() == pytest.approx(xlim_before), \
        'switching frames via the drop-down should keep the same view center'
    assert panel.ax.get_ylim() == pytest.approx(ylim_before), \
        'switching frames via the drop-down should keep the same view center'


def test_zoom_and_recenter_affect_scroll_range(qapp, focus_sweep):
    results = _step_results(focus_sweep)
    panel = ImagePanel()
    panel.show_result(results[0])

    assert panel.h_scroll.maximum() == 0, 'a fit-to-view zoom should start with nothing to scroll'
    assert panel.v_scroll.maximum() == 0, 'a fit-to-view zoom should start with nothing to scroll'

    panel.set_zoom(2.0)
    assert panel.h_scroll.maximum() > 0, 'zooming in should open up horizontal scroll range'
    assert panel.v_scroll.maximum() > 0, 'zooming in should open up vertical scroll range'

    panel.recenter()
    assert panel.h_scroll.maximum() == 0, 'recenter should return to a fit-to-view zoom'


def test_update_result_replaces_in_place(qapp, focus_sweep):
    results = _step_results(focus_sweep)
    panel = ImagePanel()
    for r in results:
        panel.add_result(r)

    updated = dataclasses.replace(results[2], fwhm=results[2].fwhm + 10., is_outlier=True)
    panel.update_result(updated)

    assert len(panel._results) == len(results), \
        'update_result() should not add a new entry for an already-known exposure'
    assert panel.exposure_combo.count() == len(results), \
        'update_result() should not add a new dropdown entry either'
    assert panel._results[2].fwhm == updated.fwhm, 'the stored measurement should be replaced'


def test_update_result_redisplays_if_it_is_the_current_one(qapp, focus_sweep):
    results = _step_results(focus_sweep)
    panel = ImagePanel()
    for r in results:
        panel.add_result(r)  # panel is now showing results[-1] (the latest added)

    updated = dataclasses.replace(results[-1], is_outlier=True)
    panel.update_result(updated)

    assert panel._current is updated, 'updating the currently-displayed result should redisplay it'


def test_update_result_falls_back_to_add_for_an_unknown_exposure(qapp, focus_sweep):
    results = _step_results(focus_sweep)
    panel = ImagePanel()
    panel.add_result(results[0])

    panel.update_result(results[1])

    assert len(panel._results) == 2, \
        'update_result() for an exposure not yet shown should behave like add_result()'


def test_wheel_zoom_keeps_the_cursor_point_fixed(qapp, focus_sweep):
    results = _step_results(focus_sweep)
    panel = ImagePanel()
    panel.show_result(results[0])  # 101x101 frame, fit-to-view (zoom=1)

    # Scroll in centered on a point away from the image center; if the
    # zoom weren't anchored to the cursor, the view would instead stay
    # anchored at the (0, 0) top-left corner and this point would drift
    # away from the cursor's screen position.
    cursor_x, cursor_y = 20.0, 80.0
    event = MouseEvent('scroll_event', panel.canvas, x=10, y=10, button='up')
    event.inaxes = panel.ax
    event.xdata, event.ydata = cursor_x, cursor_y

    panel._on_wheel_zoom(event)

    assert panel._zoom > 1.0, 'scrolling up should zoom in'
    xlim, ylim = panel.ax.get_xlim(), panel.ax.get_ylim()
    assert xlim[0] <= cursor_x <= xlim[1], 'the cursor location should remain within view'
    assert ylim[0] <= cursor_y <= ylim[1], 'the cursor location should remain within view'
    # The cursor's fractional position within the view should be
    # essentially unchanged (it started at fraction (0.2, 0.8) of the
    # full 0-101 frame).
    frac_x = (cursor_x - xlim[0]) / (xlim[1] - xlim[0])
    frac_y = (cursor_y - ylim[0]) / (ylim[1] - ylim[0])
    assert frac_x == pytest.approx(0.2, abs=0.02), 'the cursor should stay over the same data point'
    assert frac_y == pytest.approx(0.8, abs=0.02), 'the cursor should stay over the same data point'


def test_reset_clears_everything(qapp, focus_sweep):
    results = _step_results(focus_sweep)
    panel = ImagePanel()
    for r in results:
        panel.add_result(r)

    panel.reset()

    assert panel.exposure_combo.count() == 0, 'reset should clear the dropdown'
    assert panel._current is None, 'reset should clear the currently displayed result'
    assert panel._results == [], 'reset should clear the stored results'


def test_every_stretch_option_renders_without_error(qapp, focus_sweep):
    result = _step_results(focus_sweep)[0]
    panel = ImagePanel()

    for name in panel.STRETCHES:
        panel.stretch_combo.setCurrentText(name)
        panel.show_result(result)
        assert panel._stretch_name == name
        assert panel.ax.get_images(), f'{name} should have actually drawn an image'


def test_get_and_set_settings_state_round_trip(qapp):
    panel = ImagePanel()
    panel.stretch_combo.setCurrentText('Min/Max')

    state = panel.get_settings_state()
    assert state == {'stretch': 'Min/Max'}

    fresh = ImagePanel()
    fresh.set_settings_state(state)
    assert fresh._stretch_name == 'Min/Max'
    assert fresh.stretch_combo.currentText() == 'Min/Max'


def test_set_settings_state_ignores_an_unknown_stretch_name(qapp):
    panel = ImagePanel()
    original = panel._stretch_name

    panel.set_settings_state({'stretch': 'NotARealStretch'})

    assert panel._stretch_name == original, 'an unrecognized stretch name should be ignored'
