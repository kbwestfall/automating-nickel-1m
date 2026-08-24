"""Tests for :mod:`gui.views.focus_curve_panel`."""
import dataclasses

import pytest

import focus
from gui.views.focus_curve_panel import FocusCurvePanel


def _step_results(focus_sweep):
    seq = focus.ArchiveFocusSequence(focus_sweep['focus_values'], focus_sweep['files'])
    return list(seq.step(method='brightest'))


def test_add_result_accumulates_points(qapp, focus_sweep):
    results = _step_results(focus_sweep)
    panel = FocusCurvePanel()

    for r in results[:2]:
        panel.add_result(r)

    assert len(panel._results) == 2, 'panel should track every added result'
    assert panel.ax.collections, 'a scatter plot should exist once results are added'


def test_curve_and_vertex_appear_once_three_points_exist(qapp, focus_sweep):
    results = _step_results(focus_sweep)
    panel = FocusCurvePanel()

    for r in results[:2]:
        panel.add_result(r)
    assert not panel.ax.lines, 'fewer than 3 points: no fitted curve should be drawn yet'

    panel.add_result(results[2])
    assert panel.ax.lines, 'with 3+ points, the fitted quadratic curve should be drawn'


def test_curve_fit_recovers_known_best_focus(qapp, focus_sweep):
    results = _step_results(focus_sweep)
    panel = FocusCurvePanel()
    for r in results:
        panel.add_result(r)

    vertex_x = None
    for coll in panel.ax.collections:
        if coll.get_label().startswith('Best focus'):
            vertex_x = coll.get_offsets()[0][0]

    assert vertex_x is not None, 'a labeled best-focus point should be plotted'
    assert vertex_x == pytest.approx(focus_sweep['best_focus'], abs=5.), \
        'plotted vertex should match the known best focus'


def test_outlier_points_are_drawn_as_a_distinct_series(qapp, focus_sweep):
    results = _step_results(focus_sweep)
    normal = [dataclasses.replace(r, is_outlier=False) for r in results[:2]]
    outlier = dataclasses.replace(results[2], is_outlier=True)

    panel = FocusCurvePanel()
    for r in normal:
        panel.add_result(r)
    panel.add_result(outlier)

    labels = [coll.get_label() for coll in panel.ax.collections]
    assert 'Measured' in labels, 'normal points should be drawn as their own series'
    assert 'Outlier' in labels, 'outlier points should be drawn as a visually distinct series'


def test_xlabel_is_not_clipped_when_the_panel_is_short(qapp, focus_sweep):
    # A fixed-margin figure can have its bottom label clipped once the
    # panel is squeezed shorter than that margin needs -- constrained
    # layout is what's supposed to prevent that.
    results = _step_results(focus_sweep)
    panel = FocusCurvePanel()
    for r in results:
        panel.add_result(r)
    panel.resize(300, 150)  # deliberately shorter than a fixed margin would assume
    panel.show()
    qapp.processEvents()

    panel.figure.canvas.draw()
    renderer = panel.figure.canvas.get_renderer()
    xlabel_bbox = panel.ax.xaxis.label.get_window_extent(renderer=renderer)

    assert xlabel_bbox.y0 >= panel.figure.bbox.y0, \
        "the x-axis label should stay within the figure's bounds, not get clipped off the bottom"


def test_reset_clears_the_plot(qapp, focus_sweep):
    results = _step_results(focus_sweep)
    panel = FocusCurvePanel()
    for r in results:
        panel.add_result(r)

    panel.reset()

    assert panel._results == [], 'reset should clear the stored results'
    assert not panel.ax.collections, 'reset should clear the scatter points'
    assert not panel.ax.lines, 'reset should clear the fitted curve'
