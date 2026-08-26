"""
Tests for :func:`focus.FocusSequence.take_single_exposure`, and the fix
it enables to :func:`focus.FocusSequence.execute`'s previously-broken
``goto=True`` path.

See claude/GUI_IMPLEMENTATION.md, Phase 3 sub-phase 2.
"""
import pytest

from nickel_focus import focus


def test_take_single_exposure_requires_hardware(focus_sweep):
    seq = focus.ArchiveFocusSequence(focus_sweep['focus_values'], focus_sweep['files'])
    with pytest.raises(ValueError, match='ktl not connected'):
        seq.take_single_exposure(350., method='brightest')


def test_take_single_exposure_moves_and_measures(fake_hardware):
    seq = focus.FocusSequence()

    result = seq.take_single_exposure(345., method='brightest')

    assert fake_hardware['focus'].current == 345., 'take_single_exposure should move the focus'
    assert result.index == 0, "a standalone exposure isn't part of a larger sequence"
    assert result.focus_value == 345., 'result should report the focus value it was taken at'
    assert result.fwhm > 0, 'measured FWHM should be positive'

    # No sequence bookkeeping should be touched -- this is a standalone
    # operation, not a step.
    assert seq.observed_focus == [], 'take_single_exposure should not append to observed_focus'
    assert seq.exposures == [], 'take_single_exposure should not append to exposures'
    assert seq.step_iter == 0, 'take_single_exposure should not advance step_iter'


def test_take_single_exposure_out_of_range_raises(fake_hardware):
    seq = focus.FocusSequence()
    with pytest.raises(ValueError, match='out of range'):
        seq.take_single_exposure(50., method='brightest')  # below the 165 minimum


def test_execute_with_goto_now_works_against_fake_hardware(fake_hardware):
    # Before this sub-phase, execute(goto=True) called the nonexistent
    # measure_fwhm() and never appended its verification exposure to
    # self.exposures before trying to read it back -- this is the
    # regression test proving that path is fixed.
    seq = focus.GridFocusSequence(340., 5., nstep=5)

    best_focus, best_fwhm = seq.execute(goto=True, plot=False, method='brightest')

    assert best_focus == pytest.approx(fake_hardware['best_focus'], abs=5.), \
        'the fitted best focus should still be correct'
    assert fake_hardware['focus'].current == best_focus, \
        'goto=True should move the telescope to the fitted best focus'
