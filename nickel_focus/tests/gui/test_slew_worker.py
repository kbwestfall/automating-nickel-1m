"""Tests for :mod:`gui.model.slew_worker`."""
from astropy.coordinates import Angle

from nickel_focus.gui.model.slew_worker import SlewWorker
from nickel_focus.gui.qt import QtCore
from nickel_focus.tests.fake_hardware import FakeTelescopePointing


def _run_and_collect(worker, timeout_ms=5000):
    """
    Start ``worker``, pump the event loop until it finishes, and return
    the lists of ``(slewFinished, slewFailed)`` payloads it emitted.

    See `test_focus_worker._run_and_collect` for why a
    ``DirectConnection`` is used here.
    """
    finished, failed = [], []
    direct = QtCore.Qt.ConnectionType.DirectConnection
    worker.slewFinished.connect(lambda: finished.append(True), direct)
    worker.slewFailed.connect(failed.append, direct)

    loop = QtCore.QEventLoop()
    worker.finished.connect(loop.quit)
    QtCore.QTimer.singleShot(timeout_ms, loop.quit)
    worker.start()
    loop.exec()

    return finished, failed


def test_slew_finished_is_emitted_on_success(qapp):
    telescope = FakeTelescopePointing()
    worker = SlewWorker(telescope, Angle(30., unit='deg'), Angle(10., unit='deg'))

    finished, failed = _run_and_collect(worker)

    assert finished == [True], 'a successful slew should emit slewFinished exactly once'
    assert failed == [], 'a successful slew should never emit slewFailed'
    assert telescope.slew_calls == [(Angle(30., unit='deg'), Angle(10., unit='deg'))], \
        'the worker should have called slew_to with the given ra/dec'


def test_slew_failed_is_emitted_when_movement_is_disabled(qapp):
    telescope = FakeTelescopePointing(movement_allowed=False)
    worker = SlewWorker(telescope, Angle(30., unit='deg'), Angle(10., unit='deg'))

    finished, failed = _run_and_collect(worker)

    assert finished == [], 'a failed slew should never emit slewFinished'
    assert len(failed) == 1, 'a failed slew should emit slewFailed exactly once'
    assert 'disabled' in failed[0], 'the failure message should explain why the slew failed'
