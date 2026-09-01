"""
Tests for :class:`nickel_focus.scripts.focus_gui.NickelFocusGUI`'s argument
parsing. Deliberately doesn't call :func:`NickelFocusGUI.main` -- that
builds a real Qt window and blocks on the event loop -- so, unlike
`focus_gui.py` itself, this needs no Qt import at all (`get_parser`/
`parse_args` never reach the function-local imports in `main`).
"""
from nickel_focus.scripts.focus_gui import NickelFocusGUI


def test_test_flag_defaults_to_false():
    args = NickelFocusGUI.parse_args([])
    assert args.test is False, '--test should default to off'


def test_test_flag_can_be_set():
    args = NickelFocusGUI.parse_args(['--test'])
    assert args.test is True, '--test should be settable on the command line'


def test_test_flag_is_suppressed_from_help(capsys):
    NickelFocusGUI.get_parser().print_help()
    captured = capsys.readouterr()
    assert '--test' not in captured.out, \
        '--test is a developer-only escape hatch and should not appear in --help'
