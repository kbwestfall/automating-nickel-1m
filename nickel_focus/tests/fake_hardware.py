"""
Lightweight stand-ins for :class:`focus.Focus`/:class:`focus.Exposure`,
used to test :class:`focus.GridFocusSequence`/
:class:`focus.AutomatedFocusSequence` without any real ``ktl`` connection
or telescope/camera hardware.

See ``claude/GUI_IMPLEMENTATION.md``, Phase 3 sub-phase 1, for why this
exists: without it, live sequences can't be run or tested at all on a
machine with no hardware access, since ``self._focus``/``self._exposure``
would just be ``None``. These fakes are injected by monkeypatching
``focus.Focus``/``focus.Exposure``/``focus.ktl`` directly (see the
``fake_hardware`` fixture in ``conftest.py``), not through new
constructor parameters on :class:`focus.FocusSequence` -- keeping faith
with the earlier decision that it shouldn't need injection hooks for
this.

They validate the *software logic* (stepping, focus-range checks,
exposure configuration) -- they cannot validate the real KTL keyword
protocol, which can only be checked at the telescope.
"""
from astropy.io import fits


class FakeFocus:
    """Stand-in for :class:`focus.Focus`, with no ktl/hardware dependency."""

    def __init__(self, start=340.):
        self.current = start

    def set_to(self, focus_value):
        if focus_value < 165 or focus_value > 500:
            raise ValueError(f'Requested focus ({focus_value}) is out of range (165-500).')
        self.current = focus_value


class FakeExposurePath:
    """Stand-in for :class:`focus.ExposurePath`."""

    def __init__(self):
        self.previous = None


class FakeExposureConfig:
    """Stand-in for :class:`focus.ExposureConfig`; just records what it's told."""

    def __init__(self):
        self.record = None
        self.speed = None
        self.binning = None
        self.exptime = None

    def configure(self, record=None, speed=None, binning=None, exptime=None):
        if record is not None:
            self.record = record
        if speed is not None:
            self.speed = speed
        if binning is not None:
            self.binning = binning
        if exptime is not None:
            self.exptime = exptime


class FakeExposure:
    """
    Stand-in for :class:`focus.Exposure`. Instead of talking to real
    camera hardware, :func:`expose` synthesizes a FITS frame for whatever
    focus value ``focus_obj.current`` holds at the time it's called, via
    ``make_frame``, and writes it to ``directory``.

    Parameters
    ----------
    focus_obj : FakeFocus
        The (shared) fake focus object whose current value determines
        what frame gets synthesized -- this must be the same object a
        `focus.FocusSequence` is using as its ``_focus``, so that a
        `step_focus()` call immediately before `expose()` is reflected.
    make_frame : callable
        ``make_frame(focus_value) -> numpy.ndarray``, e.g. a closure over
        :func:`conftest.gaussian_frame` with a chosen focus-vs-FWHM
        relationship.
    directory : pathlib.Path
        Where to write each synthesized exposure.
    """

    def __init__(self, focus_obj, make_frame, directory):
        self._focus = focus_obj
        self._make_frame = make_frame
        self._directory = directory
        self._count = 0
        self.path = FakeExposurePath()
        self.cfg = FakeExposureConfig()

    def expose(self, record=None, speed=None, binning=None, exptime=None):
        self.cfg.configure(record=record, speed=speed, binning=binning, exptime=exptime)
        data = self._make_frame(self._focus.current)
        self._count += 1
        path = self._directory / f'fake{self._count:04d}.fits'
        fits.writeto(path, data, overwrite=True)
        self.path.previous = path
