"""
Tests for the ktl-free parts of :mod:`focus`: the
:func:`focus.FocusSequence.step`/:func:`focus.FocusSequence.execute`
sequence-stepping engine, exercised via :class:`focus.ArchiveFocusSequence`
so no :class:`focus.Focus`/:class:`focus.Exposure` hardware objects are
ever constructed.
"""
from pathlib import Path

import numpy as np
import pytest

import focus
from photometry import image_quality


def _make_sequence(focus_sweep):
    return focus.ArchiveFocusSequence(focus_sweep['focus_values'], focus_sweep['files'])


def test_archive_sequence_has_no_hardware(focus_sweep):
    seq = _make_sequence(focus_sweep)
    assert seq._focus is None
    assert seq._exposure is None


def test_step_yields_one_result_per_exposure(focus_sweep):
    seq = _make_sequence(focus_sweep)

    results = list(seq.step(method='brightest'))

    assert len(results) == len(focus_sweep['files'])
    for i, result in enumerate(results):
        assert result.index == i
        assert result.focus_value == focus_sweep['focus_values'][i]
        assert isinstance(result.exposure, Path)
        assert result.exposure == focus_sweep['files'][i]
        assert result.fwhm > 0
        assert result.frame.shape == (101, 101)
        assert result.stamp.size > 0
        assert isinstance(result.is_outlier, (bool, np.bool_))

    # The sequence's own bookkeeping should match what step() yielded.
    assert seq.observed_focus == focus_sweep['focus_values']
    assert seq.img_quality == [r.fwhm for r in results]
    assert seq.step_iter == len(focus_sweep['files'])


def test_step_can_be_driven_one_at_a_time(focus_sweep):
    # Mirrors how a future GUI worker thread would drive the sequence:
    # advancing it one step at a time (so it can check a stop flag between
    # steps) rather than consuming it in one `for` loop.
    seq = _make_sequence(focus_sweep)
    stepper = seq.step(method='brightest')

    first = next(stepper)
    assert first.index == 0
    assert seq.step_iter == 1

    second = next(stepper)
    assert second.index == 1
    assert seq.step_iter == 2


def test_execute_recovers_known_best_focus(focus_sweep):
    seq = _make_sequence(focus_sweep)

    best_focus, best_fwhm = seq.execute(goto=False, plot=False, method='brightest')

    assert seq.plot is None
    assert best_focus == pytest.approx(focus_sweep['best_focus'], abs=5.)

    # The moment-based FWHM in photometry.py is measured over a
    # thresholded source mask, which is expected to be systematically
    # biased relative to the Gaussian profile's analytic FWHM (truncating
    # the wings before computing the second moment underestimates it) --
    # so check self-consistency against a direct measurement of the frame
    # taken nearest the true best focus, rather than against the analytic
    # input value.
    middle = focus_sweep['focus_values'].index(focus_sweep['best_focus'])
    _, _, _, direct_fwhm, _, _ = image_quality(focus_sweep['files'][middle],
                                                method='brightest')
    assert best_fwhm == pytest.approx(direct_fwhm, rel=0.3)


def test_execute_with_plot_does_not_require_a_display(focus_sweep):
    # `plot=True` must still work under the non-interactive Agg backend
    # conftest.py selects, so this exercises FocusPlot without needing a
    # real display.
    seq = _make_sequence(focus_sweep)

    seq.execute(goto=False, plot=True, method='brightest')

    assert seq.plot is not None


def test_execute_with_goto_requires_hardware(focus_sweep):
    seq = _make_sequence(focus_sweep)
    with pytest.raises(ValueError, match='ktl not connected'):
        seq.execute(goto=True, plot=False)


def test_fit_best_focus_requires_at_least_three_points():
    with pytest.raises(ValueError):
        focus.FocusSequence.fit_best_focus([340., 345.], [4., 3.])
