"""
Regression guard for the optional-Qt boundary: `scripts/` (the CLI) must
never *import* a Qt binding, directly or indirectly, so it keeps working
in an environment with no Qt installed at all. Deliberately doesn't need
Qt (or even `gui/`) to run, so it always executes as part of the base
suite.

Checked via the AST rather than a raw text search, so that merely
*mentioning* a Qt binding's name -- e.g. in a comment or docstring
explaining why a signal connection must be kept alive, as
`scripts/focus_gui.py` does -- doesn't trip the guard. What actually
matters is whether the module would require Qt to be installed just to
import it; `focus_gui.py` legitimately needs to talk about Qt while
deferring its real `nickel_focus.gui` imports to inside a function body,
so this check only walks `import`/`from ... import` statements, not
arbitrary text.
"""
import ast
from pathlib import Path

_SCRIPTS_DIR = Path(__file__).resolve().parent.parent / 'scripts'
_QT_MARKERS = ('PySide6', 'PySide2', 'PyQt6', 'PyQt5')


def _imported_qt_bindings(path):
    """Return which of `_QT_MARKERS` `path` imports, directly or via a submodule."""
    tree = ast.parse(path.read_text(), filename=str(path))
    found = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            names = [alias.name for alias in node.names]
        elif isinstance(node, ast.ImportFrom):
            names = [node.module] if node.module else []
        else:
            continue
        for name in names:
            root = name.split('.')[0]
            if root in _QT_MARKERS:
                found.add(root)
    return found


def test_scripts_do_not_reference_qt():
    offenders = {
        path.name: bindings for path in _SCRIPTS_DIR.glob('*.py')
        if (bindings := _imported_qt_bindings(path))
    }
    assert not offenders, f'scripts/ must not import a Qt binding, found in: {offenders}'
