"""
Regression guard for the optional-Qt boundary: `scripts/` (the CLI) must
never reference a Qt binding, directly or indirectly, so it keeps working
in an environment with no Qt installed at all. Deliberately doesn't need
Qt (or even `gui/`) to run, so it always executes as part of the base
suite.
"""
from pathlib import Path

_SCRIPTS_DIR = Path(__file__).resolve().parent.parent / 'scripts'
_QT_MARKERS = ('PySide6', 'PySide2', 'PyQt6', 'PyQt5')


def test_scripts_do_not_reference_qt():
    offenders = [
        path.name for path in _SCRIPTS_DIR.glob('*.py')
        if any(marker in path.read_text() for marker in _QT_MARKERS)
    ]
    assert not offenders, f'scripts/ must not reference a Qt binding, found in: {offenders}'
