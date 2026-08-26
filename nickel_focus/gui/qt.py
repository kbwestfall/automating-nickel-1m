"""
Single point of contact with the Qt binding used by this GUI.

Every other :mod:`gui` module should get its Qt classes from here (e.g.
``from gui.qt import QtWidgets``) rather than importing PySide6 directly,
so there's exactly one place that fails, with one clear message, if Qt
isn't installed. `scripts/` (the CLI) never imports this module and must
never need to -- installing the base package alone is enough to run the
CLI with no Qt involved at all; the ``gui`` extra in ``pyproject.toml``
holds the additional dependency this module needs.
"""
import os

try:
    from PySide6 import QtCore, QtGui, QtWidgets
except ImportError as e:
    raise ImportError(
        'The GUI requires PySide6, which is not installed in this environment.  Install it '
        'with, e.g., "pip install .[gui]" (see also the base "dependencies" in pyproject.toml '
        'for the dependencies shared with the CLI, which do not include Qt).'
    ) from e

# Pin matplotlib's Qt backend (matplotlib.backends.backend_qtagg) to PySide6,
# so views that embed a matplotlib canvas (e.g. ImagePanel) don't depend on
# whatever binding matplotlib's own auto-detection happens to prefer.
os.environ.setdefault('QT_API', 'pyside6')

__all__ = ['QtCore', 'QtGui', 'QtWidgets']
