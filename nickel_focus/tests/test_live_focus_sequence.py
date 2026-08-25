"""
Tests for :class:`focus.GridFocusSequence`/
:class:`focus.AutomatedFocusSequence`, driven against the fake hardware
in :mod:`fake_hardware` (see the ``fake_hardware`` fixture in
``conftest.py``) rather than real ``ktl``/telescope/camera hardware.

This is the first automated coverage either of these two classes has
had: :mod:`test_focus_sequence` only exercises
:class:`focus.ArchiveFocusSequence`, which never touches
:class:`focus.Focus`/:class:`focus.Exposure` at all. These tests
validate the *software logic* -- stepping, focus-range checks, the
adaptive curve-following algorithm -- not the real KTL keyword protocol,
which can only be checked at the telescope.
"""
import pytest

from nickel_focus import focus


def test_hardware_is_really_injected(fake_hardware):
    seq = focus.GridFocusSequence(340., 5., nstep=3)
    assert seq._focus is fake_hardware['focus'], 'FocusSequence should use the injected fake Focus'
    assert seq._exposure is not None, 'FocusSequence should construct a fake Exposure'


def test_grid_sequence_commands_expected_focus_values(fake_hardware):
    seq = focus.GridFocusSequence(340., 5., nstep=5)

    results = list(seq.step(method='brightest'))

    assert [r.focus_value for r in results] == [340., 345., 350., 355., 360.], \
        'GridFocusSequence should command exactly its target_focus values, in order'
    assert fake_hardware['focus'].current == 360., \
        'the fake Focus should end at the last commanded value'

    best_focus, best_fwhm = seq.fit_best_focus(seq.observed_focus, seq.img_quality)
    assert best_focus == pytest.approx(fake_hardware['best_focus'], abs=5.), \
        'fitting a grid run against the fake hardware should recover the known best focus'


def test_grid_sequence_rejects_out_of_range_focus(fake_hardware):
    seq = focus.GridFocusSequence(160., 5., nstep=3)  # 160 is below the 165 minimum
    with pytest.raises(ValueError, match='out of range'):
        list(seq.step(method='brightest'))


def test_automated_sequence_converges_near_best_focus(fake_hardware):
    seq = focus.AutomatedFocusSequence(340., 5., maxsteps=12)

    results = list(seq.step(method='brightest'))

    assert 3 <= len(results) <= 12, 'the adaptive search should finish within its step budget'
    best_focus, best_fwhm = seq.fit_best_focus(seq.observed_focus, seq.img_quality)
    assert best_focus == pytest.approx(fake_hardware['best_focus'], abs=5.), \
        'the adaptive curve-following search should converge near the known best focus'


def test_automated_sequence_starting_past_best_focus_still_converges(fake_hardware):
    # Starting on the *other* side of the true best focus (350) exercises
    # the "direction = -1" branch of step_focus(), not just the default
    # "search upward first" case above.
    seq = focus.AutomatedFocusSequence(365., 5., maxsteps=12)

    list(seq.step(method='brightest'))

    best_focus, best_fwhm = seq.fit_best_focus(seq.observed_focus, seq.img_quality)
    assert best_focus == pytest.approx(fake_hardware['best_focus'], abs=5.), \
        'the adaptive search should converge regardless of which side it starts on'
