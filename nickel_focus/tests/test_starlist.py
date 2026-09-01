"""Tests for :mod:`starlist`."""
from importlib import resources

from astropy import units
from astropy.coordinates import FK5
import pytest

from nickel_focus import starlist


def _fields(directive=None):
    return starlist.compile_data_directive(directive)


def _entries(lines):
    """Parse raw lines into `StarlistEntry` objects, bypassing `parse_starlist`
    (which now combines entries into a `Table` before returning)."""
    fields = _fields()
    return [
        starlist.parse_data_line(line, fields, lineno=i)
        for i, line in enumerate(lines, start=1)
    ]


class TestStandardFormat:
    """Section I: the default (no ``!Data`` directive) line format."""

    @pytest.mark.parametrize('line', [
        'obj1a 12 34 56 1 2 3 2000.0',
        'obj1b 12.58222222 1 2 3 2000.0',
        'obj1c 12 34.9333333 1 2 3 2000.0',
        'obj1d 12 34 56 1.034166667 2000.0',
        'obj1e 12 34 56 1 2.05 2000.0',
    ])
    def test_decimal_and_sexagesimal_forms_agree(self, line):
        entry = starlist.parse_data_line(line, _fields())
        assert entry.ra.to('deg').value == pytest.approx(188.733333333, abs=1e-6), (
            'decimal and sexagesimal RA forms should parse to the same angle'
        )
        assert entry.dec.deg == pytest.approx(1.034166667, abs=1e-6), (
            'decimal and sexagesimal Dec forms should parse to the same angle'
        )
        assert entry.equinox == 2000.0, 'equinox should be parsed as a float'
        assert entry.frame == 'fk5', 'equinox 2000.0 with no B/J prefix should default to FK5'

    def test_colon_separated_position(self):
        entry = starlist.parse_data_line('star 12:34:56 -01:02:03 2000.0', _fields())
        assert entry.ra.to('hourangle').value == pytest.approx(12.582222, abs=1e-5), (
            'colon-separated RA should parse to the expected hour angle'
        )
        assert entry.dec.deg == pytest.approx(-1.034166667, abs=1e-6), (
            'colon-separated negative Dec should parse to the expected angle'
        )

    def test_space_separated_negative_sign(self):
        entry = starlist.parse_data_line('star 12 34 56 - 1 2 3 2000.0', _fields())
        assert entry.dec.deg == pytest.approx(-1.034166667, abs=1e-6), (
            'a space-separated leading minus sign should negate the Dec value'
        )

    def test_keyword_value_fields(self):
        entry = starlist.parse_data_line(
            'star 12 34 56 1 2 3 2000.0 pmra=12.3 pmdec=-4.5 pmepoch=1995 pri=2 Vmag=15.1 '
            'a trailing comment',
            _fields(),
        )
        assert entry.pm_ra == 12.3, 'pmra= keyword should set pm_ra'
        assert entry.pm_dec == -4.5, 'pmdec= keyword should set pm_dec'
        assert entry.pm_epoch == 1995.0, 'pmepoch= keyword should set pm_epoch'
        assert entry.priority == 2.0, 'pri= keyword should set priority'
        assert entry.mag == {'V': 15.1}, 'Vmag= keyword should set a V-band magnitude'
        assert entry.comment == 'a trailing comment', (
            'text after the last key=value pair should become the comment'
        )

    def test_backward_compatible_bare_magnitude(self):
        entry = starlist.parse_data_line('star 12 34 56 1 2 3 2000.0 15.036 a comment', _fields())
        assert entry.mag == {'mag': 15.036}, (
            'a bare number immediately after the equinox should be treated as a magnitude'
        )
        assert entry.comment == 'a comment', (
            'text after the bare magnitude should become the comment'
        )

    def test_non_numeric_token_after_equinox_is_not_a_magnitude(self):
        entry = starlist.parse_data_line('star 12 34 56 1 2 3 2000.0 not a magnitude', _fields())
        assert entry.mag == {}, 'a non-numeric token after the equinox should not be a magnitude'
        assert entry.comment == 'not a magnitude', (
            'a non-numeric token after the equinox should become part of the comment'
        )

    @pytest.mark.parametrize('equinox_token,frame', [
        ('2000.0', 'fk5'),
        ('1950.0', 'fk4'),
        ('1975.0', 'fk4'),
        ('J2000.0', 'fk5'),
        ('B1950.0', 'fk4'),
    ])
    def test_equinox_prefix_and_default_frame(self, equinox_token, frame):
        entry = starlist.parse_data_line(f'star 12 34 56 1 2 3 {equinox_token}', _fields())
        assert entry.frame == frame, f'equinox {equinox_token!r} should imply frame {frame!r}'

    def test_missing_required_field_raises(self):
        with pytest.raises(ValueError):
            starlist.parse_data_line('star 12 34 56 1 2 3', _fields())


class TestCommentDirective:
    """Section II: ``!Comment`` directives and blank lines."""

    def test_default_comment_pattern(self):
        patterns = starlist.compile_comment_patterns()
        assert starlist.is_comment_line('  # a comment', patterns), (
            'a line starting with optional whitespace then # should be a comment'
        )
        assert starlist.is_comment_line('', patterns), 'a blank line should always be a comment'
        assert not starlist.is_comment_line('star 12 34 56 1 2 3 2000.0', patterns), (
            'a normal data line should not be treated as a comment'
        )

    def test_custom_comment_patterns(self):
        patterns = starlist.compile_comment_patterns('{^#} {Object.*RA}')
        assert starlist.is_comment_line('#junk', patterns), (
            'a custom pattern matching lines starting with # should mark them as comments'
        )
        assert starlist.is_comment_line('Object   RA   Dec', patterns), (
            'a custom pattern matching "Object...RA" should mark that line as a comment'
        )
        assert not starlist.is_comment_line('star 12 34 56 1 2 3 2000.0', patterns), (
            'a normal data line should not match either custom comment pattern'
        )


class TestDataDirective:
    """Section III: ``!Data`` directives."""

    def test_fixed_width_name_and_hms_shorthand(self):
        fields = _fields('{name %20} ra_hms dec_dms equinox mag keyval {comment *}')
        entry = starlist.parse_data_line(
            'Some Object Name    00:55:16 +01:01:58 2000.0 15.036 pmra=12.3 nice star',
            fields,
        )
        assert entry.name == 'Some Object Name', (
            'a %20 fixed-width field should capture the object name, trailing spaces stripped'
        )
        assert entry.ra.to('hourangle').value == pytest.approx(0.921111, abs=1e-5), (
            'the ra_hms shorthand should parse a colon-separated RA correctly'
        )
        assert entry.dec.deg == pytest.approx(1.032778, abs=1e-5), (
            'the dec_dms shorthand should parse a colon-separated Dec correctly'
        )
        assert entry.mag == {'mag': 15.036}, 'mag should be parsed after the ra/dec/equinox fields'
        assert entry.pm_ra == 12.3, 'pmra= should still be parsed via the keyval field'
        assert entry.comment == 'nice star', 'trailing text should become the comment'

    def test_mag_before_literal_equinox_default(self):
        fields = _fields('name ra_h ra_m ra_s dec_d dec_m dec_s mag {equinox 2000.0} {comment *}')
        entry = starlist.parse_data_line(
            'XXX92.412 00 55 16 +01 01 58  15.036  a comment', fields,
        )
        assert entry.mag == {'mag': 15.036}, 'mag should be read from the line before equinox'
        assert entry.equinox == 2000.0, 'the literal {equinox 2000.0} default should be used'
        assert entry.comment == 'a comment', 'text after mag should become the comment'

    def test_skip_field(self):
        fields = _fields('name ra_hms dec_dms {equinox 2000.0} {skip %11} mag {comment *}')
        # After the dec field, the 11-character skip field consumes the leading
        # space plus a 10-character token, landing exactly on the next field.
        entry = starlist.parse_data_line(
            'XX92.412 00:55:16 +01:01:58 yadda12345 15.036 a comment', fields,
        )
        assert entry.mag == {'mag': 15.036}, 'mag should be read correctly after the skip field'
        assert entry.comment == 'a comment', 'comment should be read correctly after mag'

    def test_empty_directive_restores_default(self):
        assert _fields('') == _fields(None) == _fields(), (
            'an empty or None !Data directive should restore the default field list'
        )

    def test_fieldformat_forbidden_on_hms_shorthand(self):
        with pytest.raises(ValueError):
            _fields('name {ra_hms %s} dec_dms equinox')

    def test_orphaned_subfield_raises(self):
        with pytest.raises(ValueError):
            _fields('name ra_m ra_s dec_d dec_m dec_s equinox')

    def test_unrecognized_fieldname_is_captured_as_extra(self):
        # A field name outside the documented set (here mimicking the doc's own
        # 'epoch' vs. 'equinox' example 4) is captured rather than rejected,
        # as long as the required 'equinox' field is supplied elsewhere.
        fields = _fields('name ra_hms dec_dms {foo bar} equinox mag {comment *}')
        entry = starlist.parse_data_line('star 00:55:16 +01:01:58 2000.0 15.036 a comment', fields)
        assert entry.extra == {'foo': 'bar'}, (
            'an unrecognized field name should be captured in entry.extra, not rejected'
        )
        assert entry.equinox == 2000.0, 'the equinox field should still be parsed normally'

    def test_missing_equinox_field_raises(self):
        fields = _fields('name ra_hms dec_dms {epoch 2000.0} mag {comment *}')
        with pytest.raises(ValueError):
            starlist.parse_data_line('star 00:55:16 +01:01:58 2000.0 15.036 a comment', fields)


class TestParseStarlist:
    """
    End-to-end parsing of a full file into a `Table`, including mid-file
    directives.
    """

    def test_directives_comments_and_blank_lines(self):
        text = [
            '!Comment {^#} {^junk}',
            '# header comment',
            'junk line to skip',
            '',
            'star1 12:34:56 -01:02:03 2000.0 mag=15.2 pri=1 a nice star',
            '!Data name ra_h ra_m ra_s dec_d dec_m dec_s equinox',
            'star2 1 2 3 -4 5 6 1950',
        ]
        table = starlist.parse_starlist(text)
        assert list(table['name']) == ['star1', 'star2'], (
            'comment lines, junk lines, and blank lines should all be skipped'
        )

        star1 = table[0]
        assert star1['mag'] == {'mag': 15.2}, 'mag=15.2 keyval should be parsed for star1'
        assert star1['priority'] == 1.0, 'pri=1 keyval should be parsed for star1'
        assert star1['comment'] == 'a nice star', 'trailing text should be the comment for star1'

        star2 = table[1]
        assert star2['src_equinox'] == 1950.0, (
            'star2 should be read using the mid-file !Data directive, giving equinox=1950'
        )
        assert star2['src_frame'] == 'fk4', 'equinox 1950 should imply the FK4 frame for star2'

    def test_empty_file_returns_empty_table(self):
        table = starlist.parse_starlist(['# nothing but a comment', ''])
        assert len(table) == 0, (
            'a file with only comment/blank lines should produce an empty table'
        )
        expected_columns = [
            'name', 'ra', 'dec', 'src_equinox', 'src_frame', 'pm_ra', 'pm_dec', 'pm_epoch',
            'priority', 'comment', 'mag', 'extra'
        ]
        assert list(table.colnames) == expected_columns, (
            'an empty table should still have all the expected columns'
        )

    def test_reads_from_path(self):
        path = resources.files('nickel_focus') / 'data' / 'point_focus.txt'
        table = starlist.parse_starlist(path)
        assert len(table) == 24, 'point_focus.txt should yield one row per data line'
        assert table['name'][0] == 'Pointing00', 'the first row should be the first data line'
        assert 'ICRS' in table.meta['frame'], 'the default output frame should be ICRS'

    def test_frame_and_equinox_are_forwarded(self):
        table = starlist.parse_starlist(
            ['star 12 34 56 1 2 3 2000.0'], frame='fk5', equinox='J2050.0',
        )
        expected = _entries(['star 12 34 56 1 2 3 2000.0'])[0].coord.transform_to(
            FK5(equinox='J2050.0')
        )
        assert table['ra'][0] == pytest.approx(expected.ra.deg, abs=1e-10), (
            'parse_starlist should forward frame/equinox to entries_to_table for RA'
        )
        assert table['dec'][0] == pytest.approx(expected.dec.deg, abs=1e-10), (
            'parse_starlist should forward frame/equinox to entries_to_table for Dec'
        )


class TestEntriesToTable:
    """`entries_to_table`: combining entries into a single-frame Table."""

    def test_default_icrs_frame(self):
        entries = _entries(['star 12 34 56 1 2 3 2000.0'])
        table = starlist.entries_to_table(entries)
        expected = entries[0].coord.icrs
        assert table['ra'].unit == units.deg, 'the ra column should be in degrees'
        assert table['dec'].unit == units.deg, 'the dec column should be in degrees'
        assert table['ra'][0] == pytest.approx(expected.ra.deg, abs=1e-10), (
            "the default frame should match astropy's own ICRS transform for RA"
        )
        assert table['dec'][0] == pytest.approx(expected.dec.deg, abs=1e-10), (
            "the default frame should match astropy's own ICRS transform for Dec"
        )
        assert 'ICRS' in table.meta['frame'], 'table.meta should record the ICRS frame'

    def test_mixed_source_frames_share_one_target_frame(self):
        entries = _entries(['star1950 12 34 56 1 2 3 1950.0', 'star2000 12 34 56 1 2 3 2000.0'])
        table = starlist.entries_to_table(entries)
        assert list(table['src_frame']) == ['fk4', 'fk5'], (
            "src_frame should record each entry's original frame, not the output frame"
        )
        # Same sexagesimal RA/Dec but different source equinoxes should
        # precess to different ICRS positions.
        assert table['ra'][0] != pytest.approx(table['ra'][1], abs=1e-6), (
            'FK4 1950 and FK5 2000 sources with the same sexagesimal RA should precess '
            'to different ICRS positions'
        )

    def test_explicit_frame_and_equinox(self):
        entries = _entries(['star 12 34 56 1 2 3 2000.0'])
        table = starlist.entries_to_table(entries, frame='fk5', equinox='J2050.0')
        expected = entries[0].coord.transform_to(FK5(equinox='J2050.0'))
        assert table['ra'][0] == pytest.approx(expected.ra.deg, abs=1e-10), (
            "an explicit frame/equinox should match astropy's own transform for RA"
        )
        assert table['dec'][0] == pytest.approx(expected.dec.deg, abs=1e-10), (
            "an explicit frame/equinox should match astropy's own transform for Dec"
        )

    def test_equinox_ignored_for_frames_without_it(self):
        entries = _entries(['star 12 34 56 1 2 3 2000.0'])
        # icrs has no equinox attribute; passing one should be silently
        # ignored rather than raising.
        with_equinox = starlist.entries_to_table(entries, equinox='J2050.0')
        without_equinox = starlist.entries_to_table(entries)
        assert with_equinox['ra'][0] == without_equinox['ra'][0], (
            'an unsupported equinox kwarg should be silently ignored for the icrs frame'
        )

    def test_other_fields_are_preserved(self):
        entries = _entries(['star 12 34 56 1 2 3 2000.0 pmra=1.5 pri=2 mag=15.0 a comment'])
        table = starlist.entries_to_table(entries)
        assert table['name'][0] == 'star', 'the name column should be preserved'
        assert table['pm_ra'][0] == 1.5, 'the pm_ra column should be preserved'
        assert table['priority'][0] == 2.0, 'the priority column should be preserved'
        assert table['mag'][0] == {'mag': 15.0}, 'the mag column should be preserved'
        assert table['comment'][0] == 'a comment', 'the comment column should be preserved'
        assert table['src_equinox'][0] == 2000.0, 'the src_equinox column should be preserved'
        assert table['src_frame'][0] == 'fk5', 'the src_frame column should be preserved'

    def test_empty_entries_returns_empty_table(self):
        table = starlist.entries_to_table([])
        assert len(table) == 0, 'an empty entries list should produce a zero-row table'
        assert table['ra'].unit == units.deg, (
            'the empty table should still have a degree-unit ra column'
        )
        assert table['dec'].unit == units.deg, (
            'the empty table should still have a degree-unit dec column'
        )
        assert 'ICRS' in table.meta['frame'], (
            'the empty table should still record the target frame'
        )

    def test_unrecognized_frame_name_raises(self):
        entries = _entries(['star 12 34 56 1 2 3 2000.0'])
        with pytest.raises(ValueError):
            starlist.entries_to_table(entries, frame='not_a_frame')
