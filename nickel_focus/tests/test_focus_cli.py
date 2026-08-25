"""
End-to-end test of :class:`nickel_focus.scripts.focus.NickelFocus`'s
archive/replay path (``--obsnum``/``--datadir``/``--prefix``/``--suffix``),
exercised without any ``ktl`` connection.
"""
import pytest
from astropy.table import Table

from nickel_focus.scripts.focus import NickelFocus


def test_cli_archive_mode_runs_end_to_end(focus_sweep, capsys):
    argv = [
        str(focus_sweep['focus_values'][0]), '5',
        '-n', str(len(focus_sweep['focus_values'])),
        '--obsnum', str(focus_sweep['obsnum']),
        '--datadir', str(focus_sweep['datadir']),
        '--prefix', focus_sweep['prefix'],
        '--suffix', focus_sweep['suffix'],
        '--no-plot',
    ]

    NickelFocus.main(NickelFocus.parse_args(argv))

    captured = capsys.readouterr()
    assert 'Best focus:' in captured.out, 'CLI should print the fitted best focus'
    assert 'Expected sigma:' in captured.out, 'CLI should print the expected sigma'
    assert 'Wrote focus data' not in captured.out, \
        'no output file was requested, so nothing should be written'


def test_cli_writes_output_file_when_requested(focus_sweep, capsys, tmp_path):
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

    NickelFocus.main(NickelFocus.parse_args(argv))

    captured = capsys.readouterr()
    assert f'Wrote focus data to {ofile}' in captured.out, \
        'the CLI should confirm where the output file was written'
    assert ofile.is_file(), 'the output file should actually have been written'

    tbl = Table.read(ofile, format='ecsv')
    assert len(tbl) == len(focus_sweep['focus_values']), 'one row per exposure'
    assert list(tbl['FOCUS']) == pytest.approx(focus_sweep['focus_values']), \
        'the output table should record the observed focus values'
    assert 'SIGMA' in tbl.colnames, 'the output table should record the image-quality metric'
    assert 'EXPOSURE' in tbl.colnames, 'the output table should record which file each row is from'
    assert 'OUTLIER' in tbl.colnames, 'the output table should record the outlier flag'
