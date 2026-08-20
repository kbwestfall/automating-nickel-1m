"""Tests for :mod:`gui.model.sequence_worker`."""
import pytest

import focus
from gui.qt import QtCore
from gui.model.sequence_worker import SequenceWorker


def _make_sequence(focus_sweep):
    return focus.ArchiveFocusSequence(focus_sweep['focus_values'], focus_sweep['files'])


def _run_and_collect(worker, timeout_ms=5000):
    """
    Start `worker`, pump the event loop until it finishes, and return the
    lists of ``(stepComplete, sequenceFinished, sequenceFailed)`` payloads
    it emitted.

    Connections use ``DirectConnection`` so the callbacks run
    synchronously on the worker thread as each signal is emitted --
    the simplest way to deterministically collect results in a test with
    no real widgets involved. Real GUI code relies on Qt's automatic
    queued connections instead, to safely marshal updates back to the GUI
    thread.
    """
    steps, finished, failed = [], [], []
    direct = QtCore.Qt.ConnectionType.DirectConnection
    worker.stepComplete.connect(steps.append, direct)
    worker.sequenceFinished.connect(lambda *args: finished.append(args), direct)
    worker.sequenceFailed.connect(failed.append, direct)

    loop = QtCore.QEventLoop()
    worker.finished.connect(loop.quit)
    QtCore.QTimer.singleShot(timeout_ms, loop.quit)
    worker.start()
    loop.exec()

    return steps, finished, failed


def test_invalid_mode_is_rejected(focus_sweep):
    seq = _make_sequence(focus_sweep)
    with pytest.raises(ValueError, match='mode'):
        SequenceWorker(seq, mode='bogus')


def test_step_mode_runs_full_sequence(focus_sweep, qapp):
    seq = _make_sequence(focus_sweep)
    worker = SequenceWorker(seq, method='brightest', mode='step')

    steps, finished, failed = _run_and_collect(worker)

    assert len(steps) == len(focus_sweep['files']), 'should emit one stepComplete per exposure'
    assert not failed, f'sequence should not fail: {failed}'
    assert len(finished) == 1, 'sequenceFinished should fire exactly once'
    best_focus, best_fwhm = finished[0]
    assert best_focus == pytest.approx(focus_sweep['best_focus'], abs=5.), \
        'worker should recover the known best focus'


def test_request_stop_halts_between_steps(focus_sweep, qapp):
    seq = _make_sequence(focus_sweep)
    worker = SequenceWorker(seq, method='brightest', mode='step')

    steps = []

    def on_step(result):
        steps.append(result)
        if len(steps) == 2:
            worker.request_stop()

    finished, failed = [], []
    direct = QtCore.Qt.ConnectionType.DirectConnection
    worker.stepComplete.connect(on_step, direct)
    worker.sequenceFinished.connect(lambda *args: finished.append(args), direct)
    worker.sequenceFailed.connect(failed.append, direct)

    loop = QtCore.QEventLoop()
    worker.finished.connect(loop.quit)
    QtCore.QTimer.singleShot(5000, loop.quit)
    worker.start()
    loop.exec()

    assert len(steps) == 2, 'worker should stop after the requested step, not run to completion'
    assert not finished, 'too few points remain for a fit, so sequenceFinished should not fire'
    assert len(failed) == 1, 'sequenceFailed should fire once for the early stop'
    assert 'not enough points' in failed[0].lower(), 'failure message should explain why'


def test_reanalyze_mode_updates_measurements(focus_sweep, qapp):
    seq = _make_sequence(focus_sweep)
    seq.execute(goto=False, plot=False, method='brightest')

    worker = SequenceWorker(seq, method='weighted', mode='reanalyze')
    steps, finished, failed = _run_and_collect(worker)

    assert len(steps) == len(focus_sweep['files']), 'reanalyze should emit one result per exposure'
    assert not failed, f'reanalysis should not fail: {failed}'
    assert len(finished) == 1, 'sequenceFinished should fire once after reanalysis'
    assert seq.method == 'weighted', 'reanalyze should update the sequence method'


def test_exp_kwargs_defaults_to_an_empty_dict(focus_sweep):
    seq = _make_sequence(focus_sweep)
    worker = SequenceWorker(seq)
    assert worker.exp_kwargs == {}, 'omitting exp_kwargs should default to an empty dict, not None'


def test_exp_kwargs_is_a_noop_for_archive_mode(focus_sweep, qapp):
    # ArchiveFocusSequence has no exposure hardware at all (_exposure is
    # None); passing exp_kwargs must not try to call .cfg.configure() on
    # it, which would be an AttributeError on None.
    seq = _make_sequence(focus_sweep)
    worker = SequenceWorker(seq, mode='step', exp_kwargs={'exptime': 12.5})

    steps, finished, failed = _run_and_collect(worker)

    assert not failed, f'exp_kwargs should be harmless when there is no exposure hardware: {failed}'
    assert len(finished) == 1


def test_exp_kwargs_configures_the_exposure_before_stepping(fake_hardware, qapp):
    seq = focus.GridFocusSequence(340., 5., nstep=5)
    worker = SequenceWorker(seq, mode='step', exp_kwargs={'exptime': 12.5, 'speed': 'Fast'})

    steps, finished, failed = _run_and_collect(worker)

    assert not failed, f'sequence should not fail: {failed}'
    assert seq._exposure.cfg.exptime == 12.5, 'exposure time should be configured before stepping'
    assert seq._exposure.cfg.speed == 'Fast', 'speed should be configured before stepping'


def test_exp_kwargs_are_ignored_in_reanalyze_mode(fake_hardware, qapp):
    seq = focus.GridFocusSequence(340., 5., nstep=5)
    seq.execute(goto=False, plot=False, method='brightest')
    assert seq._exposure.cfg.exptime is None, 'setup: nothing configured yet'

    worker = SequenceWorker(seq, mode='reanalyze', exp_kwargs={'exptime': 12.5})
    steps, finished, failed = _run_and_collect(worker)

    assert not failed, f'reanalysis should not fail: {failed}'
    assert seq._exposure.cfg.exptime is None, \
        'reanalyze never takes new exposures, so exp_kwargs should be ignored'
