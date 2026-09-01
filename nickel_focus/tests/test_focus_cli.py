"""
End-to-end test of :class:`nickel_focus.scripts.focus.NickelFocus`'s
archive/replay path (``--obsnum``/``--datadir``/``--prefix``/``--suffix``),
exercised without any ``ktl`` connection.
"""
import logging

import pytest
from astropy.table import Table

from nickel_focus import log
from nickel_focus.scripts.focus import NickelFocus


def test_cli_archive_mode_runs_end_to_end(focus_sweep, caplog):
    argv = [
        str(focus_sweep['focus_values'][0]), '5',
        '-n', str(len(focus_sweep['focus_values'])),
        '--obsnum', str(focus_sweep['obsnum']),
        '--datadir', str(focus_sweep['datadir']),
        '--prefix', focus_sweep['prefix'],
        '--suffix', focus_sweep['suffix'],
        '--no-plot',
    ]

    with caplog.at_level(logging.INFO, logger='nickel_focus'):
        NickelFocus.main(NickelFocus.parse_args(argv))

    messages = [r.getMessage() for r in caplog.records if r.name == 'nickel_focus']
    assert any('Best focus:' in m for m in messages), 'CLI should log the fitted best focus'
    assert any('Expected sigma:' in m for m in messages), 'CLI should log the expected sigma'
    assert not any('Wrote focus data' in m for m in messages), \
        'no output file was requested, so nothing should be written'


def test_cli_writes_output_file_when_requested(focus_sweep, caplog, tmp_path):
    ofile = tmp_path / 'focus_data.ecsv'
    argv = [
        str(focus_sweep['focus_values'][0]), '5',
        '-n', str(len(focus_sweep['focus_values'])),
        '--obsnum', str(focus_sweep['obsnum']),
        '--datadir', str(focus_sweep['datadir']),
        '--prefix', focus_sweep['prefix'],
        '--suffix', focus_sweep['suffix'],
        '--no-plot',
        '--ofile', str(ofile),
    ]

    with caplog.at_level(logging.INFO, logger='nickel_focus'):
        NickelFocus.main(NickelFocus.parse_args(argv))

    messages = [r.getMessage() for r in caplog.records if r.name == 'nickel_focus']
    assert any(f'Wrote focus data to {ofile}' in m for m in messages), \
        'the CLI should confirm where the output file was written'
    assert ofile.is_file(), 'the output file should actually have been written'

    tbl = Table.read(ofile, format='ecsv')
    assert len(tbl) == len(focus_sweep['focus_values']), 'one row per exposure'
    assert list(tbl['FOCUS']) == pytest.approx(focus_sweep['focus_values']), \
        'the output table should record the observed focus values'
    assert 'SIGMA' in tbl.colnames, 'the output table should record the image-quality metric'
    assert 'EXPOSURE' in tbl.colnames, 'the output table should record which file each row is from'
    assert 'OUTLIER' in tbl.colnames, 'the output table should record the outlier flag'


def test_cli_verbosity_flag_controls_console_level():
    # caplog can't observe this: it captures independently of whatever
    # handlers the app configures, which is exactly the point of caplog,
    # but it means it can't tell us whether -v actually changed what a
    # user watching the console would see. Check the console handler's
    # level directly instead.
    NickelFocus.init_log(NickelFocus.parse_args(['300', '5', '-v', '0']))
    assert log.sh.level == logging.WARNING, '-v 0 should set the console handler to warning'

    NickelFocus.init_log(NickelFocus.parse_args(['300', '5', '-v', '1']))
    assert log.sh.level == logging.INFO, \
        '-v 1 (the default) should set the console handler to info'

    NickelFocus.init_log(NickelFocus.parse_args(['300', '5', '-v', '2']))
    assert log.sh.level == logging.DEBUG, '-v 2 should set the console handler to debug'
