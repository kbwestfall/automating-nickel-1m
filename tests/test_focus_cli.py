"""
End-to-end test of :func:`focus.main`'s archive/replay path
(``--obsnum``/``--datadir``/``--prefix``/``--suffix``), exercised without
any ``ktl`` connection.
"""
import sys

import focus


def test_cli_archive_mode_runs_end_to_end(focus_sweep, capsys, monkeypatch):
    argv = [
        'focus.py',
        str(focus_sweep['focus_values'][0]), '5',
        '-n', str(len(focus_sweep['focus_values'])),
        '--obsnum', str(focus_sweep['obsnum']),
        '--datadir', str(focus_sweep['datadir']),
        '--prefix', focus_sweep['prefix'],
        '--suffix', focus_sweep['suffix'],
        '--no-plot',
    ]
    monkeypatch.setattr(sys, 'argv', argv)

    focus.main()

    captured = capsys.readouterr()
    assert 'Best focus:' in captured.out, 'CLI should print the fitted best focus'
    assert 'Expected sigma:' in captured.out, 'CLI should print the expected sigma'
