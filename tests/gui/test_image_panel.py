"""Tests for :mod:`gui.views.image_panel`."""
import dataclasses

import pytest
from matplotlib.backend_bases import MouseEvent
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


def test_click_emits_signal_only_when_selection_enabled(qapp, focus_sweep):
    results = _step_results(focus_sweep)
    panel = ImagePanel()
    panel.add_result(results[0])

    received = []
    panel.sourceSelected.connect(lambda x, y: received.append((x, y)))

    event = MouseEvent('button_press_event', panel.canvas, x=10, y=10, button=1)
    event.inaxes = panel.ax
    event.xdata, event.ydata = 42.0, 17.0

    panel._on_click(event)
    assert received == [], 'click should be ignored while selection is disabled'

    panel.set_selection_enabled(True)
    panel._on_click(event)
    assert received == [(42.0, 17.0)], 'click should emit the clicked data coordinates once enabled'


def test_outlier_box_color_reflects_is_outlier(qapp, focus_sweep):
    results = _step_results(focus_sweep)
    panel = ImagePanel()

    panel.show_result(dataclasses.replace(results[0], is_outlier=False))
    colors = [p.get_edgecolor() for p in panel.ax.patches]
    assert to_rgba('red') in colors, 'a non-outlier result should be boxed in red'

    panel.show_result(dataclasses.replace(results[0], is_outlier=True))
    colors = [p.get_edgecolor() for p in panel.ax.patches]
    assert to_rgba('yellow') in colors, 'an outlier result should be boxed in yellow'


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


def test_reset_clears_everything(qapp, focus_sweep):
    results = _step_results(focus_sweep)
    panel = ImagePanel()
    for r in results:
        panel.add_result(r)

    panel.reset()

    assert panel.exposure_combo.count() == 0, 'reset should clear the dropdown'
    assert panel._current is None, 'reset should clear the currently displayed result'
    assert panel._results == [], 'reset should clear the stored results'
