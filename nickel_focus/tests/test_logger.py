"""Tests for :mod:`pkg.logger`."""
import logging

from nickel_focus.pkg.logger import GuiFormatter


def _make_record(level, message):
    return logging.LogRecord(
        name='nickel_focus', level=level, pathname=__file__, lineno=1, msg=message, args=(),
        exc_info=None)


def test_gui_formatter_has_no_color_escape_codes():
    record = _make_record(logging.INFO, 'hello')
    text = GuiFormatter().format(record)
    assert '\x1b' not in text, 'GUI formatter output must not contain ANSI escape codes'


def test_gui_formatter_includes_level_and_message():
    record = _make_record(logging.WARNING, 'careful now')
    text = GuiFormatter().format(record)
    assert 'WARNING' in text, 'formatted line should include the level name'
    assert 'careful now' in text, 'formatted line should include the message'
