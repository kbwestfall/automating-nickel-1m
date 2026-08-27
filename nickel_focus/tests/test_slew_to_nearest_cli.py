"""
End-to-end test of
:class:`nickel_focus.scripts.slew_to_nearest.NickelSlewToNearest`,
exercised without any ``ktl`` connection or telescope hardware (see the
``fake_telescope`` fixture in ``conftest.py``).
"""
from astropy.coordinates import SkyCoord
import pytest

from nickel_focus import slew
from nickel_focus.scripts.slew_to_nearest import NickelSlewToNearest


def _write_starlist(tmp_path, lines):
    """
    Write a small starlist file for a test to point ``-f``/``--starlist``
    at.

    Parameters
    ----------
    tmp_path
        Directory to write the file into (typically pytest's ``tmp_path``
        fixture).
    lines
        The starlist's data lines (standard format; see
        `nickel_focus.starlist`), without trailing newlines.

    Returns
    -------
    pathlib.Path
        Path to the written file.
    """
    path = tmp_path / 'stars.txt'
    path.write_text('\n'.join(lines) + '\n')
    return path


def test_cli_dry_run_does_not_move_the_telescope(fake_telescope, tmp_path, capsys):
    path = _write_starlist(tmp_path, ['StarA 01:00:00 +10:00:00 2000.0'])
    fake_telescope.ra, fake_telescope.dec = 1.0, 10.0

    NickelSlewToNearest.main(
        NickelSlewToNearest.parse_args(['-f', str(path), '--dry_run'])
    )

    captured = capsys.readouterr()
    assert 'Nearest object is: StarA' in captured.out, 'CLI should report the nearest target'
    assert fake_telescope.slew_calls == [], 'dry-run mode should never command a slew'


def test_cli_slews_to_the_nearest_target(fake_telescope, tmp_path):
    path = _write_starlist(tmp_path, [
        'StarA 01:00:00 +10:00:00 2000.0',
        'StarB 13:00:00 -20:00:00 2000.0',
    ])
    fake_telescope.ra, fake_telescope.dec = 1.0, 10.0
    telescope_coo = SkyCoord(ra=fake_telescope.ra, dec=fake_telescope.dec,
                              unit=('hourangle', 'deg'))
    _, expected_ra, expected_dec = slew.find_nearest_target(telescope_coo, file=path)

    NickelSlewToNearest.main(NickelSlewToNearest.parse_args(['-f', str(path)]))

    assert len(fake_telescope.slew_calls) == 1, 'the telescope should be commanded to slew once'
    commanded_ra, commanded_dec = fake_telescope.slew_calls[0]
    assert commanded_ra.deg == pytest.approx(expected_ra.deg, abs=1e-10), (
        'the telescope should be commanded to the RA that find_nearest_target selected -- '
        'regression check for the RA-degrees-vs-hours unit bug'
    )
    assert commanded_dec.deg == pytest.approx(expected_dec.deg, abs=1e-10), (
        'the telescope should be commanded to the Dec that find_nearest_target selected'
    )


def test_cli_search_string_selects_matching_target(fake_telescope, tmp_path, capsys):
    path = _write_starlist(tmp_path, [
        'Pointing01 05:00:00 +30:00:00 2000.0',
        'Focusing01 05:00:00 +30:00:00 2000.0',
    ])
    fake_telescope.ra, fake_telescope.dec = 5.0, 30.0

    NickelSlewToNearest.main(
        NickelSlewToNearest.parse_args(['-f', str(path), '-s', 'Focusing', '--dry_run'])
    )

    captured = capsys.readouterr()
    assert 'Nearest object is: Focusing01' in captured.out, (
        '--search should restrict the candidate list to matching names, even though '
        'Pointing01 is an equally close non-match'
    )


def test_cli_raises_when_search_string_matches_nothing(fake_telescope, tmp_path):
    path = _write_starlist(tmp_path, ['StarA 01:00:00 +10:00:00 2000.0'])
    fake_telescope.ra, fake_telescope.dec = 1.0, 10.0

    argv = ['-f', str(path), '-s', 'Nonexistent']
    with pytest.raises(ValueError, match='no object names containing'):
        NickelSlewToNearest.main(NickelSlewToNearest.parse_args(argv))


@pytest.mark.parametrize('flag,message', [
    ('movement_allowed', 'Telescope movement is disabled!'),
    ('tracking_on', 'Tracking is disabled!'),
    ('ready', 'Telescope is not ready to move to a new target!'),
    ('reaches_target', 'Telescope failed to make it to target within 5 min.'),
])
def test_cli_propagates_telescope_failures(fake_telescope, tmp_path, flag, message):
    path = _write_starlist(tmp_path, ['StarA 01:00:00 +10:00:00 2000.0'])
    fake_telescope.ra, fake_telescope.dec = 1.0, 10.0
    setattr(fake_telescope, flag, False)

    with pytest.raises(ValueError, match=message):
        NickelSlewToNearest.main(NickelSlewToNearest.parse_args(['-f', str(path)]))
