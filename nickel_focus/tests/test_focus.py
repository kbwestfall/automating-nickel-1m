"""Tests for the ktl-backed parts of :mod:`focus`."""
import logging

from nickel_focus import focus


class _FakeKeyword:
    """A minimal stand-in for a ktl keyword: read/write a value, always ready."""

    def __init__(self, value=None):
        self.value = value

    def read(self):
        return self.value

    def write(self, value):
        self.value = value

    def waitFor(self, condition, timeout=None):
        return True


class _FakeKtl:
    """
    A minimal stand-in for the ``ktl`` module, just enough to drive
    :meth:`focus.Focus.set_to` through its golden path: every keyword is
    always ready/allowed, except ``POCSECPD``, which starts away from
    ``secpd_start`` so a `set_to` call doesn't take the early "already at
    that position" return.
    """

    def __init__(self, secpd_start=300.):
        self._secpd_start = secpd_start
        self._keywords = {}

    def cache(self, system, name):
        if name not in self._keywords:
            start = self._secpd_start if name == 'POCSECPD' else 0.
            self._keywords[name] = _FakeKeyword(start)
        return self._keywords[name]


def test_set_to_logs_info(monkeypatch, caplog):
    monkeypatch.setattr(focus, 'ktl', _FakeKtl(secpd_start=300.))
    f = focus.Focus()

    with caplog.at_level(logging.INFO, logger='nickel_focus'):
        f.set_to(350.)

    messages = [r.getMessage() for r in caplog.records if r.name == 'nickel_focus']
    assert any('Unlocking secondary' in m for m in messages), \
        'set_to() should log that it is unlocking the secondary'
    assert any('Successfully changed focus to 350' in m for m in messages), \
        'set_to() should log the successful focus change'


def test_set_to_logs_info_when_already_at_position(monkeypatch, caplog):
    monkeypatch.setattr(focus, 'ktl', _FakeKtl(secpd_start=350.))
    f = focus.Focus()

    with caplog.at_level(logging.INFO, logger='nickel_focus'):
        f.set_to(350.)

    messages = [r.getMessage() for r in caplog.records if r.name == 'nickel_focus']
    assert any('No change needed' in m for m in messages), \
        'set_to() should log when the focus is already at the requested position'
