"""Tests for :mod:`gui.model.focus_worker`."""
import pytest

from nickel_focus import focus
from nickel_focus.gui.qt import QtCore
from nickel_focus.gui.model.focus_worker import FocusWorker


def _make_sequence(focus_sweep):
    return focus.ArchiveFocusSequence(focus_sweep['focus_values'], focus_sweep['files'])


def _run_and_collect(worker, timeout_ms=5000):
    """
    Start `worker`, pump the event loop until it finishes, and return the
    lists of ``(stepComplete, focusSequenceFinished, focusSequenceFailed,
    singleExposureFinished)`` payloads it emitted.

    Connections use ``DirectConnection`` so the callbacks run
    synchronously on the worker thread as each signal is emitted --
    the simplest way to deterministically collect results in a test with
    no real widgets involved. Real GUI code relies on Qt's automatic
    queued connections instead, to safely marshal updates back to the GUI
    thread.
    """
    steps, finished, failed, single = [], [], [], []
    direct = QtCore.Qt.ConnectionType.DirectConnection
    worker.stepComplete.connect(steps.append, direct)
    worker.focusSequenceFinished.connect(lambda *args: finished.append(args), direct)
    worker.focusSequenceFailed.connect(failed.append, direct)
    worker.singleExposureFinished.connect(single.append, direct)

    loop = QtCore.QEventLoop()
    worker.finished.connect(loop.quit)
    QtCore.QTimer.singleShot(timeout_ms, loop.quit)
    worker.start()
    loop.exec()

    return steps, finished, failed, single


def test_invalid_mode_is_rejected(focus_sweep):
    seq = _make_sequence(focus_sweep)
    with pytest.raises(ValueError, match='mode'):
        FocusWorker(seq, mode='bogus')


def test_step_mode_runs_full_sequence(focus_sweep, qapp):
    seq = _make_sequence(focus_sweep)
    worker = FocusWorker(seq, method='brightest', mode='step')

    steps, finished, failed, single = _run_and_collect(worker)

    assert len(steps) == len(focus_sweep['files']), 'should emit one stepComplete per exposure'
    assert not failed, f'sequence should not fail: {failed}'
    assert len(finished) == 1, 'focusSequenceFinished should fire exactly once'
    best_focus, best_fwhm = finished[0]
    assert best_focus == pytest.approx(focus_sweep['best_focus'], abs=5.), \
        'worker should recover the known best focus'


def test_request_stop_halts_between_steps(focus_sweep, qapp):
    seq = _make_sequence(focus_sweep)
    worker = FocusWorker(seq, method='brightest', mode='step')

    steps = []

    def on_step(result):
        steps.append(result)
        if len(steps) == 2:
            worker.request_stop()

    finished, failed = [], []
    direct = QtCore.Qt.ConnectionType.DirectConnection
    worker.stepComplete.connect(on_step, direct)
    worker.focusSequenceFinished.connect(lambda *args: finished.append(args), direct)
    worker.focusSequenceFailed.connect(failed.append, direct)

    loop = QtCore.QEventLoop()
    worker.finished.connect(loop.quit)
    QtCore.QTimer.singleShot(5000, loop.quit)
    worker.start()
    loop.exec()

    assert len(steps) == 2, 'worker should stop after the requested step, not run to completion'
    assert not finished, \
        'too few points remain for a fit, so focusSequenceFinished should not fire'
    assert len(failed) == 1, 'focusSequenceFailed should fire once for the early stop'
    assert 'not enough points' in failed[0].lower(), 'failure message should explain why'


def test_reanalyze_mode_updates_measurements(focus_sweep, qapp):
    seq = _make_sequence(focus_sweep)
    seq.execute(goto=False, plot=False, method='brightest')

    worker = FocusWorker(seq, method='weighted', mode='reanalyze')
    steps, finished, failed, single = _run_and_collect(worker)

    assert len(steps) == len(focus_sweep['files']), 'reanalyze should emit one result per exposure'
    assert not failed, f'reanalysis should not fail: {failed}'
    assert len(finished) == 1, 'focusSequenceFinished should fire once after reanalysis'
    assert seq.method == 'weighted', 'reanalyze should update the sequence method'


def test_exp_kwargs_defaults_to_an_empty_dict(focus_sweep):
    seq = _make_sequence(focus_sweep)
    worker = FocusWorker(seq)
    assert worker.exp_kwargs == {}, 'omitting exp_kwargs should default to an empty dict, not None'


def test_exp_kwargs_is_a_noop_for_archive_mode(focus_sweep, qapp):
    # ArchiveFocusSequence has no exposure hardware at all (_exposure is
    # None); passing exp_kwargs must not try to call .cfg.configure() on
    # it, which would be an AttributeError on None.
    seq = _make_sequence(focus_sweep)
    worker = FocusWorker(seq, mode='step', exp_kwargs={'exptime': 12.5})

    steps, finished, failed, single = _run_and_collect(worker)

    assert not failed, f'exp_kwargs should be harmless when there is no exposure hardware: {failed}'
    assert len(finished) == 1


def test_exp_kwargs_configures_the_exposure_before_stepping(fake_hardware, qapp):
    seq = focus.GridFocusSequence(340., 5., nstep=5)
    worker = FocusWorker(seq, mode='step', exp_kwargs={'exptime': 12.5, 'speed': 'Fast'})

    steps, finished, failed, single = _run_and_collect(worker)

    assert not failed, f'sequence should not fail: {failed}'
    assert seq._exposure.cfg.exptime == 12.5, 'exposure time should be configured before stepping'
    assert seq._exposure.cfg.speed == 'Fast', 'speed should be configured before stepping'


def test_exp_kwargs_are_ignored_in_reanalyze_mode(fake_hardware, qapp):
    seq = focus.GridFocusSequence(340., 5., nstep=5)
    seq.execute(goto=False, plot=False, method='brightest')
    assert seq._exposure.cfg.exptime is None, 'setup: nothing configured yet'

    worker = FocusWorker(seq, mode='reanalyze', exp_kwargs={'exptime': 12.5})
    steps, finished, failed, single = _run_and_collect(worker)

    assert not failed, f'reanalysis should not fail: {failed}'
    assert seq._exposure.cfg.exptime is None, \
        'reanalyze never takes new exposures, so exp_kwargs should be ignored'


def test_single_mode_requires_a_focus_value(focus_sweep):
    seq = _make_sequence(focus_sweep)
    with pytest.raises(ValueError, match='focus_value'):
        FocusWorker(seq, mode='single')


def test_single_mode_moves_and_exposes(fake_hardware, qapp):
    seq = focus.GridFocusSequence(340., 5., nstep=5)
    worker = FocusWorker(seq, method='brightest', mode='single', focus_value=356.,
                          exp_kwargs={'exptime': 8.0})

    steps, finished, failed, single = _run_and_collect(worker)

    assert not failed, f'single exposure should not fail: {failed}'
    assert not steps, 'mode=single has no generator, so stepComplete should never fire'
    assert not finished, \
        'mode=single never fits a curve, so focusSequenceFinished should never fire'
    assert len(single) == 1, 'singleExposureFinished should fire exactly once'
    assert single[0].focus_value == 356., 'the fake hardware should have been moved to focus_value'
    assert seq._exposure.cfg.exptime == 8.0, 'exp_kwargs should be applied before exposing'


def test_single_mode_reports_failure_without_hardware(focus_sweep, qapp):
    # ArchiveFocusSequence has no _focus/_exposure at all, so
    # take_single_exposure() should raise and the worker should report it
    # via focusSequenceFailed rather than crashing the thread.
    seq = _make_sequence(focus_sweep)
    worker = FocusWorker(seq, mode='single', focus_value=356.)

    steps, finished, failed, single = _run_and_collect(worker)

    assert not single, 'singleExposureFinished should not fire on failure'
    assert len(failed) == 1, 'focusSequenceFailed should fire once'
    assert 'move to best focus' in failed[0].lower(), 'failure message should explain the action'
