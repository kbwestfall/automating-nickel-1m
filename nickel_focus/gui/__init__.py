"""
Qt-based GUI for the Nickel focus-finding workflow.

Unlike :mod:`~nickel_focus.scripts` (the CLI), which must keep working with no
Qt installed at all, this package requires a Qt binding -- see
:mod:`~nickel_focus.gui.qt` for the one place that dependency is actually
imported. Importing :mod:`~nickel_focus.gui` alone doesn't force a Qt
dependency; only importing a submodule that actually needs it (via
:mod:`~nickel_focus.gui.qt`) does.
"""
