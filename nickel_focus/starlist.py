"""
Parser for Lick Observatory "starlist" files, as read by the ``coords`` program
used at the 1-m (Nickel) and 3-m (Shane) telescopes.

The format is documented at
https://mthamilton.ucolick.org/techdocs/telescopes/starlistFormat.html; this
module implements the standard line format (section I), the ``!Comment``
directive (section II), and the ``!Data`` directive (section III) described
there.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

from astropy import units
from astropy.coordinates import Angle, SkyCoord, frame_transform_graph
from astropy.table import Table

_TOKEN_RE = re.compile(r'\{[^}]*\}|\S+')
_KEYVAL_RE = re.compile(r'^(?P<key>[A-Za-z0-9_]+)=(?P<val>.+)$')


@dataclass
class StarlistEntry:
    """
    A single target parsed from a starlist file.

    Parameters
    ----------
    name
        The object's name.
    ra
        Right ascension.
    dec
        Declination.
    equinox
        The equinox value (e.g.  ``2000.0``).
    frame
        Reference frame implied by the equinox: ``'fk4'`` or ``'fk5'``.
    pm_ra
        Proper motion in RA, mas/yr.
    pm_dec
        Proper motion in Dec, mas/yr.
    pm_epoch
        Epoch of the proper motion; defaults to the equinox if unset.
    mag
        Magnitudes, keyed by band letter (or ``'mag'`` if unlabeled).
    priority
        Priority of the object, from a ``pri=`` keyword field.
    comment
        Any trailing comment text on the line.
    extra
        Any other ``key=value`` pairs, or otherwise-named ``!Data`` fields, not
        otherwise recognized.
    line
        The original raw line text this entry was parsed from.
    lineno
        The line number this entry was parsed from, if known.
    """
    name: str
    ra: Angle
    dec: Angle
    equinox: float
    frame: str
    pm_ra: float = 0.0
    pm_dec: float = 0.0
    pm_epoch: float | None = None
    mag: dict = field(default_factory=dict)
    priority: float | None = None
    comment: str = ''
    extra: dict = field(default_factory=dict)
    line: str = ''
    lineno: int | None = None

    @property
    def coord(self) -> SkyCoord:
        """
        The target position as an :class:`astropy.coordinates.SkyCoord`, in this
        entry's own reference frame and equinox (rather than any common frame).

        Returns
        -------
        astropy.coordinates.SkyCoord
            The position.
        """
        prefix = 'B' if self.frame == 'fk4' else 'J'
        return SkyCoord(
            ra=self.ra, dec=self.dec, frame=self.frame, equinox=f'{prefix}{self.equinox}'
        )


class _LineCursor:
    """
    Tracks a read position while sequentially consuming fields from a starlist
    data line.

    Parameters
    ----------
    line
        The raw data line to be parsed, in its entirety.
    """

    def __init__(self, line: str):
        self.line = line
        self.pos = 0

    def peek_token(self) -> str | None:
        """
        Look ahead at the next whitespace-delimited token without consuming it.

        Returns
        -------
        str or None
            The next token, or ``None`` if only whitespace (or nothing) remains
            in the line.
        """
        m = re.match(r'\s*(\S+)', self.line[self.pos:])
        return None if m is None else m.group(1)

    def take_token(self) -> str:
        """
        Consume and return the next whitespace-delimited token.

        Returns
        -------
        str
            The consumed token.

        Raises
        ------
        ValueError
            Raised if only whitespace (or nothing) remains in the line.
        """
        m = re.match(r'\s*(\S+)', self.line[self.pos:])
        if m is None:
            raise ValueError(f'Expected another field in line: {self.line!r}')
        self.pos += m.end()
        return m.group(1)

    def take_fixed(self, width: int) -> str:
        """
        Consume and return a fixed-width slice of the line, counted from the
        current position without first skipping leading whitespace.

        Parameters
        ----------
        width
            The number of characters to consume.

        Returns
        -------
        str
            The consumed slice, exactly ``width`` characters long (including any
            leading/trailing whitespace it contains).

        Raises
        ------
        ValueError
            Raised if fewer than ``width`` characters remain in the line.
        """
        if self.pos + width > len(self.line):
            raise ValueError(
                f'Line is too short to extract a {width}-character field: {self.line!r}'
            )
        value = self.line[self.pos:self.pos + width]
        self.pos += width
        return value

    def take_rest(self) -> str:
        """
        Consume and return everything remaining in the line.

        Returns
        -------
        str
            The remaining line text, stripped of leading/trailing whitespace.
        """
        value = self.line[self.pos:].strip()
        self.pos = len(self.line)
        return value


def _consume_angle(
    cursor: _LineCursor, unit: units.Unit, colon_required: bool, has_minsec: bool
) -> Angle:
    """
    Consume 1-3 whitespace-delimited tokens describing an hours- or
    degrees-based angle, honoring the decimal-fraction and colon shorthand rules
    in section I/III of the format document.

    Parameters
    ----------
    cursor
        The line cursor to consume tokens from.
    unit
        The angular unit (``units.hourangle`` or ``units.deg``) that a bare
        (non-colon, non-decimal) leading value is expressed in.
    colon_required
        If ``True`` (the ``ra_hms``/``ra_dms``/``dec_dms`` shorthands), the
        single token consumed must be colon-separated into exactly three parts
        (``h:m:s``/``d:m:s``); a plain or decimal value is rejected.
    has_minsec
        If ``True``, a bare leading value (no colon, no decimal point) is
        followed by a separate minutes token and, unless that minutes token is
        itself decimal, a separate seconds token.  If ``False`` (e.g.  a
        standalone ``ra_d``/``dec_d`` field with no declared ``_m``/``_s``
        partner fields), a bare leading value is the whole angle and no further
        tokens are consumed.

    Returns
    -------
    astropy.coordinates.Angle
        The parsed angle, in ``unit``.

    Raises
    ------
    ValueError
        Raised if ``colon_required`` is set and the token either has no colons
        or does not split into exactly three colon-separated parts.
    """
    token = cursor.take_token()

    sign = 1.0
    match token:
        case '+' | '-':
            sign = -1.0 if token == '-' else 1.0
            token = cursor.take_token()
        case _ if token[:1] in ('+', '-'):
            sign = -1.0 if token[0] == '-' else 1.0
            token = token[1:]

    if ':' in token:
        parts = token.split(':')
        if colon_required and len(parts) != 3:
            raise ValueError(f'Expected both minutes and seconds in {token!r}.')
        primary = float(parts[0])
        minutes = float(parts[1]) if len(parts) > 1 else 0.0
        seconds = float(parts[2]) if len(parts) > 2 else 0.0
        return Angle(sign * (primary + minutes / 60. + seconds / 3600.), unit)

    if colon_required:
        raise ValueError(f'Expected a colon-separated {unit} field, got {token!r}.')

    if '.' in token:
        # The fractional part supplies any remaining minutes/seconds.
        return Angle(sign * float(token), unit)

    primary = float(token)
    if not has_minsec:
        return Angle(sign * primary, unit)

    m_token = cursor.take_token()
    if '.' in m_token:
        return Angle(sign * (primary + float(m_token) / 60.), unit)

    s_token = cursor.take_token()
    return Angle(sign * (primary + float(m_token) / 60. + float(s_token) / 3600.), unit)


def _parse_equinox(token: str) -> tuple[float, str]:
    """
    Parse an ``equinox`` field.

    Parameters
    ----------
    token
        The equinox value, optionally prefixed with ``B`` or ``J`` to select the
        reference frame explicitly (e.g.  ``'J2000.0'``, ``'B1950.0'``), or a
        bare numeric value (e.g.  ``'2000.0'``).

    Returns
    -------
    value : float
        The numeric equinox value.
    frame : str
        The implied reference frame, ``'fk4'`` or ``'fk5'``: taken from an
        explicit ``B``/``J`` prefix if present, otherwise ``'fk4'`` if ``value
        <= 1975`` and ``'fk5'`` otherwise.

    Raises
    ------
    ValueError
        Raised if ``token`` does not match the expected ``[BJ]?<number>``
        format.
    """
    m = re.match(r'^(?P<prefix>[BJ])?(?P<value>[+-]?\d+\.?\d*)$', token, re.IGNORECASE)
    if m is None:
        raise ValueError(f'Cannot parse equinox from {token!r}.')
    value = float(m.group('value'))
    match m.group('prefix'):
        case None:
            frame = 'fk4' if value <= 1975 else 'fk5'
        case p if p.upper() == 'B':
            frame = 'fk4'
        case _:
            frame = 'fk5'
    return value, frame


def _parse_keyval_token(token: str, entry: StarlistEntry) -> None:
    """
    Apply one ``key=value`` token to ``entry``, updating it in place.

    Parameters
    ----------
    token
        A single ``key=value`` token (e.g.  ``'pmra=12.3'``, ``'Vmag=15.1'``).
        Recognized keys are ``pmra``, ``pmdec``, ``pmepoch``, and ``pri``; a
        band-letter-prefixed magnitude such as ``Vmag`` or ``Jmag``; a bare
        ``mag``; or a single band letter such as ``V`` or ``J``. Anything else
        is stored verbatim in ``entry.extra``.
    entry
        The :class:`~nickel_focus.starlist.StarlistEntry` to update.
    """
    m = _KEYVAL_RE.match(token)
    key, val = m.group('key'), m.group('val')
    lk = key.lower()
    match lk:
        case 'pmra':
            entry.pm_ra = float(val)
        case 'pmdec':
            entry.pm_dec = float(val)
        case 'pmepoch':
            entry.pm_epoch = float(val)
        case 'pri':
            entry.priority = float(val)
        case 'mag':
            entry.mag['mag'] = float(val)
        case _ if lk.endswith('mag') and len(lk) > 3:
            entry.mag[key[:-3]] = float(val)
        case _ if len(key) == 1 and key.isalpha():
            entry.mag[key] = float(val)
        case _:
            entry.extra[key] = val


def _resolve_generic(cursor: _LineCursor, fmt: str | None) -> str:
    """
    Resolve a non-angle, non-keyword field's value given its ``fieldformat``.

    Parameters
    ----------
    cursor
        The line cursor to consume from.
    fmt
        The field's ``fieldformat``, as produced by
        :func:`~nickel_focus.starlist._parse_field_list`: ``None`` or ``'%s'``
        for a whitespace-delimited token, ``'%N'`` (``N`` an integer) for a
        fixed-width slice, ``'*'`` for the rest of the line, or any other
        string, which is a literal default value and consumes nothing from the
        line.

    Returns
    -------
    str
        The resolved field value.
    """
    match fmt:
        case None | '%s':
            return cursor.take_token()
        case '*':
            return cursor.take_rest()
        case _ if (m := re.match(r'^%(\d+)$', fmt)) is not None:
            return cursor.take_fixed(int(m.group(1))).strip()
        case _:
            return fmt


def _parse_field_list(directive: str) -> list[tuple[str, str | None]]:
    """
    Tokenize a ``!Data`` directive into ``(fieldname, fieldformat)`` pairs.

    Parameters
    ----------
    directive
        The directive text (or default line format string) to tokenize, e.g.
        ``'name ra_h ra_m ra_s dec_d dec_m dec_s equinox mag keyval {comment
        *}'``. Each whitespace-separated token is either a bare field name
        (format defaults to ``None``) or a brace-delimited ``{fieldname
        fieldformat}`` group.

    Returns
    -------
    list
        The ``(fieldname, fieldformat)`` pairs, in the order they appear in
        ``directive``.
    """
    fields = []
    for tok in _TOKEN_RE.findall(directive):
        if tok.startswith('{') and tok.endswith('}'):
            inner = tok[1:-1].strip()
            name, _, fmt = inner.partition(' ')
            stripped_fmt = fmt.strip()
            fields.append((name, stripped_fmt if stripped_fmt != '' else None))
        else:
            fields.append((tok, None))
    return fields


def _collapse_angle_fields(fields: list[tuple[str, str | None]]) -> list:
    """
    Collapse the consecutive ``ra_h ra_m ra_s`` / ``dec_d dec_m dec_s`` triplets
    (or the ``ra_hms``/``ra_dms``/``dec_dms`` shorthands) that a ``!Data``
    directive must always specify together into single pseudo fields that
    :func:`~nickel_focus.starlist.parse_data_line` knows how to consume as a
    unit.

    Parameters
    ----------
    fields
        The ``(fieldname, fieldformat)`` pairs returned by
        :func:`~nickel_focus.starlist._parse_field_list`.

    Returns
    -------
    list
        ``fields``, with any RA/Dec sexagesimal sub-fields collapsed into single
        ``('__ra__', info)``/``('__dec__', info)`` pseudo-fields, where ``info``
        is a dict with keys ``'unit'``, ``'colon_required'``, and
        ``'has_minsec'`` as consumed by
        :func:`~nickel_focus.starlist._consume_angle`. All other fields are
        passed through unchanged.

    Raises
    ------
    ValueError
        Raised if a ``fieldformat`` is given for one of the RA/Dec primary
        fieldnames (not allowed, per the format document), or if a
        minutes/seconds sub-field (``ra_m``, ``ra_s``, ``dec_m``, or ``dec_s``)
        appears without its primary field immediately preceding it.
    """
    ra_primary_units = {'ra_h': units.hourangle, 'ra_d': units.deg}
    dec_primary_units = {'dec_d': units.deg}

    collapsed = []
    i, n = 0, len(fields)
    while i < n:
        name, fmt = fields[i]
        match name:
            case _ if name in ra_primary_units or name in dec_primary_units:
                if fmt is not None:
                    raise ValueError(f'A fieldformat cannot be used with {name!r}.')
                is_ra = name in ra_primary_units
                m_name, s_name = ('ra_m', 'ra_s') if is_ra else ('dec_m', 'dec_s')
                has_minsec = (
                    i + 2 < n and fields[i + 1][0] == m_name and fields[i + 2][0] == s_name
                )
                unit = ra_primary_units[name] if is_ra else dec_primary_units[name]
                pseudo = '__ra__' if is_ra else '__dec__'
                collapsed.append(
                    (pseudo, {'unit': unit, 'colon_required': False, 'has_minsec': has_minsec})
                )
                i += 3 if has_minsec else 1
            case 'ra_hms' | 'ra_dms' | 'dec_dms':
                if fmt is not None:
                    raise ValueError(f'A fieldformat cannot be used with {name!r}.')
                pseudo = '__dec__' if name == 'dec_dms' else '__ra__'
                unit = units.hourangle if name == 'ra_hms' else units.deg
                collapsed.append(
                    (pseudo, {'unit': unit, 'colon_required': True, 'has_minsec': True})
                )
                i += 1
            case 'ra_m' | 'ra_s' | 'dec_m' | 'dec_s':
                raise ValueError(
                    f'Field {name!r} must immediately follow its primary RA/Dec field.'
                )
            case _:
                collapsed.append((name, fmt))
                i += 1
    return collapsed


def compile_data_directive(directive: str | None = None) -> list:
    """
    Parse a ``!Data`` directive into the field-specification list consumed by
    :func:`~nickel_focus.starlist.parse_data_line`.

    Parameters
    ----------
    directive
        The text following ``!Data`` on a directive line.  If ``None`` or blank,
        the default line format is used: ``name ra_h ra_m ra_s dec_d dec_m dec_s
        equinox mag keyval {comment *}`` ("Example 1" in the format document).

    Returns
    -------
    list
        A list of ``(fieldname, fieldformat)`` tuples, with any RA/Dec
        sexagesimal sub-fields collapsed into single pseudo-fields.
    """
    if directive is None or directive.strip() == '':
        text = 'name ra_h ra_m ra_s dec_d dec_m dec_s equinox mag keyval {comment *}'
    else:
        text = directive
    return _collapse_angle_fields(_parse_field_list(text))


def compile_comment_patterns(directive: str | None = None) -> list[str]:
    """
    Parse a ``!Comment`` directive's list of regex patterns.

    Parameters
    ----------
    directive
        The text following ``!Comment`` on a directive line.  If ``None`` or
        blank, the built-in ``coords`` default (``^[ \\t]*#``, i.e.  any line
        beginning with optional whitespace then ``#``) is used.

    Returns
    -------
    list
        Regex patterns; a line is a comment if it matches any of them.
    """
    if directive is None or directive.strip() == '':
        return [r'^[ \t]*#']
    return [
        tok[1:-1] if tok.startswith('{') and tok.endswith('}') else tok
        for tok in _TOKEN_RE.findall(directive)
    ]


def is_comment_line(line: str, patterns: list[str]) -> bool:
    """
    Determine whether a line should be treated as a comment.

    Parameters
    ----------
    line
        The raw line text to test.
    patterns
        Regular expressions, as returned by
        :func:`~nickel_focus.starlist.compile_comment_patterns`. ``line`` is a
        comment if it matches (via :func:`re.search`) any of them.

    Returns
    -------
    bool
        ``True`` if ``line`` is blank, or matches one of ``patterns``.
    """
    return line.strip() == '' or any(re.search(p, line) for p in patterns)


def parse_data_line(line: str, fields: list, lineno: int | None = None) -> StarlistEntry:
    """
    Parse one starlist data line into a
    :class:`~nickel_focus.starlist.StarlistEntry`.

    Parameters
    ----------
    line
        The raw line text (not a comment or directive line).
    fields
        The field specification, as returned by
        :func:`~nickel_focus.starlist.compile_data_directive`.
    lineno
        Optional line number, recorded on the entry for error reporting.

    Returns
    -------
    StarlistEntry
        The parsed target.

    Raises
    ------
    ValueError
        Raised if the line cannot be parsed, or does not supply all of object
        name, RA, Dec, and equinox.
    """
    cursor = _LineCursor(line)
    entry = StarlistEntry(
        name='', ra=None, dec=None, equinox=None, frame=None, line=line, lineno=lineno
    )

    for name, fmt in fields:
        match name:
            case '__ra__':
                entry.ra = _consume_angle(
                    cursor, fmt['unit'], fmt['colon_required'], fmt['has_minsec']
                )
            case '__dec__':
                entry.dec = _consume_angle(
                    cursor, fmt['unit'], fmt['colon_required'], fmt['has_minsec']
                )
            case 'name':
                entry.name = _resolve_generic(cursor, fmt)
            case 'equinox':
                entry.equinox, entry.frame = _parse_equinox(_resolve_generic(cursor, fmt))
            case 'mag':
                token = cursor.peek_token()
                if token is not None and re.match(r'^[+-]?(\d+\.?\d*|\.\d+)$', token) is not None:
                    entry.mag['mag'] = float(cursor.take_token())
            case 'keyval':
                token = cursor.peek_token()
                while token is not None and _KEYVAL_RE.match(token) is not None:
                    _parse_keyval_token(cursor.take_token(), entry)
                    token = cursor.peek_token()
            case 'comment':
                entry.comment = _resolve_generic(cursor, fmt)
            case 'skip':
                _resolve_generic(cursor, fmt)
            case _:
                entry.extra[name] = _resolve_generic(cursor, fmt)

    if entry.name == '' or entry.ra is None or entry.dec is None or entry.equinox is None:
        raise ValueError(f'Line {lineno}: missing required name/RA/Dec/equinox field(s): {line!r}')

    return entry


def entries_to_table(
    entries: list[StarlistEntry], frame: str | object = 'icrs', **frame_attrs
) -> Table:
    """
    Combine a list of :class:`~nickel_focus.starlist.StarlistEntry` objects into
    a single :class:`astropy.table.Table`, with every position transformed to
    one common reference frame.

    Each entry can carry its own equinox/reference frame (``'fk4'`` or
    ``'fk5'``, set line-by-line from the starlist's ``equinox`` field), so
    positions cannot simply be stacked into a table as-is.  This function
    transforms each entry's :class:`astropy.coordinates.SkyCoord`
    (:attr:`~nickel_focus.starlist.StarlistEntry.coord`) into ``frame``
    individually before combining them, which sidesteps that mismatch.

    Parameters
    ----------
    entries
        The entries to combine, e.g.  as returned by
        :func:`~nickel_focus.starlist.parse_data_line`.
    frame
        The frame to transform every position into: a frame name recognized by
        :attr:`astropy.coordinates.frame_transform_graph` (e.g.  ``'icrs'``,
        ``'fk5'``, ``'fk4'``), or an already-constructed frame instance/class
        accepted by :meth:`astropy.coordinates.SkyCoord.transform_to`.
    **frame_attrs
        Frame attributes used to construct the target frame when ``frame`` is
        given as a name, e.g.  ``equinox='J2000.0'``. Ignored if ``frame`` is
        passed as an instance, or if the named frame does not accept a given
        attribute (as is the case for ``equinox`` with the default ``'icrs'``,
        which has no free equinox parameter).

    Returns
    -------
    astropy.table.Table
        One row per entry, with ``name``, ``ra``/``dec`` (deg, in ``frame``),
        ``pm_ra``, ``pm_dec``, ``pm_epoch``, ``priority``, ``comment``, ``mag``,
        and ``extra`` columns, plus ``src_equinox``/ ``src_frame`` recording
        each entry's original equinox/frame before transformation.  The frame
        actually used is recorded in ``Table.meta['frame']``.

    Raises
    ------
    ValueError
        Raised if ``frame`` is a string not recognized as a coordinate frame
        name.
    """
    if isinstance(frame, str):
        frame_cls = frame_transform_graph.lookup_name(frame)
        if frame_cls is None:
            raise ValueError(f'Unrecognized coordinate frame name: {frame!r}.')
        supported = {k: v for k, v in frame_attrs.items() if k in frame_cls.frame_attributes}
        target_frame = frame_cls(**supported)
    else:
        target_frame = frame

    if len(entries) == 0:
        names = [
            'name', 'ra', 'dec', 'src_equinox', 'src_frame', 'pm_ra', 'pm_dec', 'pm_epoch',
            'priority', 'comment', 'mag', 'extra'
        ]
        dtype = [str, float, float, float, str, float, float, float, float, str, object, object]
        table = Table(names=names, dtype=dtype)
        table['ra'].unit = units.deg
        table['dec'].unit = units.deg
        table.meta['frame'] = str(target_frame)
        return table

    ra = units.Quantity([e.coord.transform_to(target_frame).ra for e in entries])
    dec = units.Quantity([e.coord.transform_to(target_frame).dec for e in entries])

    table = Table({
        'name': [e.name for e in entries],
        'ra': ra.to(units.deg),
        'dec': dec.to(units.deg),
        'src_equinox': [e.equinox for e in entries],
        'src_frame': [e.frame for e in entries],
        'pm_ra': [e.pm_ra for e in entries],
        'pm_dec': [e.pm_dec for e in entries],
        'pm_epoch': [e.pm_epoch for e in entries],
        'priority': [e.priority for e in entries],
        'comment': [e.comment for e in entries],
        'mag': [e.mag for e in entries],
        'extra': [e.extra for e in entries],
    })
    table.meta['frame'] = str(target_frame)
    return table


def parse_starlist(
    source: str | Path | list, frame: str | object = 'icrs', **frame_attrs
) -> Table:
    """
    Parse a Lick Observatory starlist file into a single-frame table of targets.

    Parameters
    ----------
    source
        Path to the starlist file, or an iterable of lines (e.g., an open file
        handle or list of strings).
    frame
        The frame to transform every position into, forwarded to
        :func:`~nickel_focus.starlist.entries_to_table`: a frame name recognized
        by :attr:`astropy.coordinates.frame_transform_graph` (e.g.  ``'icrs'``,
        ``'fk5'``, ``'fk4'``), or an already-constructed frame instance/class
        accepted by :meth:`astropy.coordinates.SkyCoord.transform_to`.
    **frame_attrs
        Frame attributes forwarded to
        :func:`~nickel_focus.starlist.entries_to_table`, used to construct the
        target frame when ``frame`` is given as a name, e.g.
        ``equinox='J2000.0'``.

    Returns
    -------
    astropy.table.Table
        One row per data line in the file, in file order, as built by
        :func:`~nickel_focus.starlist.entries_to_table`. Blank lines,
        ``!Comment``-matched lines, and directive lines themselves are skipped.
    """
    if isinstance(source, (str, Path)):
        lines = Path(source).read_text().splitlines()
    else:
        lines = [str(line).rstrip('\n') for line in source]

    directive_re = re.compile(r'^!(?P<directive>Comment|Data)\b\s*(?P<args>.*)$')
    comment_patterns = compile_comment_patterns()
    fields = compile_data_directive()
    entries = []

    for lineno, line in enumerate(lines, start=1):
        directive = directive_re.match(line)
        if directive is not None:
            if directive.group('directive') == 'Comment':
                comment_patterns = compile_comment_patterns(directive.group('args'))
            else:
                fields = compile_data_directive(directive.group('args'))
            continue

        if is_comment_line(line, comment_patterns):
            continue

        entries.append(parse_data_line(line, fields, lineno=lineno))

    return entries_to_table(entries, frame=frame, **frame_attrs)
