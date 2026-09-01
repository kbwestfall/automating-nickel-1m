"""
Single point of contact with the `ktl` package.

Everything else in this package should get `ktl` from here (`from
nickel_focus import ktl`) rather than importing it directly, so there's
exactly one place that decides whether it's available, with one warning
emitted at most once (when `nickel_focus` is first imported) rather than
once per module that would otherwise each run their own `try`/`except`.

Unlike `gui.qt` (which hard-fails if PySide6 isn't installed, since the
`gui` extra requires it), `ktl` being unavailable is not an error here:
`scripts/`'s archive/replay CLI paths and `focus.ArchiveFocusSequence`
must keep working with no live `ktl` connection at all. Code that
actually talks to hardware (`focus.Focus`, `focus.Exposure`,
`slew.NickelTelescopePointing`) checks `ktl is None` itself and raises
`RuntimeError` there instead.
"""
import warnings

try:
    import ktl
except ModuleNotFoundError:
    warnings.warn('Unable to import ktl package.  Functionality will be limited.')
    ktl = None

__all__ = ['ktl']
