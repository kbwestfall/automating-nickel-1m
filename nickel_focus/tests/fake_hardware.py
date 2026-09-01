"""
Lightweight stand-ins for :class:`focus.Focus`/:class:`focus.Exposure` and
:class:`slew.NickelTelescopePointing`, used to test the sequences/scripts
built on top of them without any real ``ktl`` connection or
telescope/camera hardware.

See ``claude/GUI_IMPLEMENTATION.md``, Phase 3 sub-phase 1, for why this
exists: without it, live sequences can't be run or tested at all on a
machine with no hardware access, since ``self._focus``/``self._exposure``
would just be ``None``. These fakes are injected by monkeypatching
``focus.Focus``/``focus.Exposure``/``focus.ktl``/``slew.NickelTelescopePointing``
directly (see the ``fake_hardware``/``fake_telescope`` fixtures in
``conftest.py``), not through new constructor parameters on
:class:`focus.FocusSequence` -- keeping faith with the earlier decision
that it shouldn't need injection hooks for this.

They validate the *software logic* (stepping, focus-range checks,
exposure configuration, target selection) -- they cannot validate the
real KTL keyword protocol, which can only be checked at the telescope.
"""
from astropy.coordinates import SkyCoord
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


class FakeTelescopePointing:
    """
    Stand-in for :class:`slew.NickelTelescopePointing`, with no
    ktl/hardware dependency.

    `slew_to` reproduces the real class's four checks/failure modes (in
    the same order, raising the same messages), each individually
    toggleable via a constructor flag, so callers such as
    :class:`nickel_focus.scripts.slew_to_nearest.NickelSlewToNearest` can
    be exercised against both the success path and each failure path.

    Parameters
    ----------
    ra : float
        Initial right ascension, hours. Public, and meant to be read or
        set directly by a test both before and after a `slew_to` call.
    dec : float
        Initial declination, degrees. Public, for the same reason as
        ``ra``.
    movement_allowed : bool
        Whether `slew_to` passes the "movement disabled" check.
    tracking_on : bool
        Whether `slew_to` passes the "tracking disabled" check.
    ready : bool
        Whether `slew_to` passes the "not ready for a new target" check.
    reaches_target : bool
        Whether `slew_to` passes the final "failed to reach target"
        check.
    """

    def __init__(self, ra=0.0, dec=0.0, movement_allowed=True, tracking_on=True,
                 ready=True, reaches_target=True):
        self.ra = ra
        self.dec = dec
        self.movement_allowed = movement_allowed
        self.tracking_on = tracking_on
        self.ready = ready
        self.reaches_target = reaches_target
        self.slew_calls = []

    @property
    def current(self):
        """
        Stand-in for :attr:`slew.NickelTelescopePointing.current`.

        Returns
        -------
        astropy.coordinates.SkyCoord
            The current RA/Dec, built from ``self.ra``/``self.dec``.
        """
        return SkyCoord(ra=self.ra, dec=self.dec, unit=('hourangle', 'deg'))

    def slew_to(self, ra, dec):
        """
        Stand-in for :meth:`slew.NickelTelescopePointing.slew_to`.

        Every call is recorded in ``self.slew_calls`` (as the
        ``(ra, dec)`` arguments given), and, if it doesn't raise,
        ``self.ra``/``self.dec`` are updated to reflect the new pointing.

        Parameters
        ----------
        ra
            Right ascension of the target; anything accepted by
            `~astropy.coordinates.SkyCoord`'s ``ra`` argument.
        dec
            Declination of the target; anything accepted by
            `~astropy.coordinates.SkyCoord`'s ``dec`` argument.

        Raises
        ------
        ValueError
            Raised if ``self.movement_allowed``, ``self.tracking_on``,
            ``self.ready``, or ``self.reaches_target`` is ``False``,
            with the same messages `slew.NickelTelescopePointing.slew_to`
            would raise for the corresponding real failure.
        """
        if not self.movement_allowed:
            raise ValueError('Telescope movement is disabled!')
        if not self.tracking_on:
            raise ValueError('Tracking is disabled!')
        if not self.ready:
            raise ValueError('Telescope is not ready to move to a new target!')

        target_coo = SkyCoord(ra=ra, dec=dec, unit=('hourangle', 'deg'))
        self.slew_calls.append((ra, dec))

        if not self.reaches_target:
            raise ValueError('Telescope failed to make it to target within 5 min.')

        self.ra = target_coo.ra.to('hourangle').value
        self.dec = target_coo.dec.to('deg').value
