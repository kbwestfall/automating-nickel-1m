"""
Qt-based GUI for the Nickel focus-finding workflow.

Unlike `scripts/` (the CLI), which must keep working with no Qt
installed at all, this package requires a Qt binding -- see :mod:`gui.qt`
for the one place that dependency is actually imported.

Importing :mod:`gui` adds `scripts/` to :obj:`sys.path`, so the shared
Model code (`focus`, `photometry`, `quadratic`) can be imported the same
flat way `scripts/focus.py` imports its own siblings, without requiring
Qt itself -- that way, `import gui` alone doesn't force a Qt dependency;
only importing a submodule that actually needs it (via :mod:`gui.qt`)
does.
"""
import sys
from pathlib import Path

_SCRIPTS_DIR = Path(__file__).resolve().parent.parent / 'scripts'
if str(_SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS_DIR))
