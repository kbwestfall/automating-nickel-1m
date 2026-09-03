"""
Single point of contact with the :mod:`ktl` package.

Everything else in this package should get :mod:`ktl` from here (``from
nickel_focus import ktl``) rather than importing it directly, so there's exactly
one place that decides whether it's available, with one warning emitted at most
once (when :mod:`~nickel_focus` is first imported) rather than once per module
that would otherwise each run their own ``try``/``except``.

Unlike :mod:`~nickel_focus.gui.qt` (which hard-fails if PySide6 isn't installed,
since the ``gui`` extra requires it), :mod:`ktl` being unavailable is not an
error here: ``scripts/``'s archive/replay CLI paths and
:class:`~nickel_focus.focus.ArchiveFocusSequence` must keep working with no live
:mod:`ktl` connection at all.  Code that actually talks to hardware
(:class:`~nickel_focus.focus.Focus`, :class:`~nickel_focus.focus.Exposure`,
:class:`~nickel_focus.slew.NickelTelescopePointing`) checks ``ktl is None``
itself and raises :class:`RuntimeError` there instead.
"""
import warnings

try:
    import ktl
except ModuleNotFoundError:
    warnings.warn('Unable to import ktl package.  Functionality will be limited.')
    ktl = None

__all__ = ['ktl']
