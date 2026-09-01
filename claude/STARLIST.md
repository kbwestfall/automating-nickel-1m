# Starlist Parser — Design Notes

**Status:** Implemented, tested, not yet wired into `slew.py`.

## Purpose

Lick's `coords` program (used at the 1-m Nickel and 3-m Shane telescopes)
reads target catalogs in a "starlist" format, documented at
https://mthamilton.ucolick.org/techdocs/telescopes/starlistFormat.html.
That format is more than the simple whitespace-delimited layout used by
the existing `nickel_focus/data/point_focus.txt` — it also supports a
standard-format line with optional keyword fields (proper motion,
magnitude, priority), a `!Comment` directive for filtering junk lines, and
a `!Data` directive that lets a file redefine the column order/format
entirely (fixed-width fields, literal defaults, key=value matching, etc.).

`nickel_focus/starlist.py` implements this format, with tests in
`nickel_focus/tests/test_starlist.py`.

## Design decisions

- **One engine, not two.** Rather than writing a separate code path for
  "standard format" versus "custom `!Data` format," the standard format is
  itself just the default field list (`name ra_h ra_m ra_s dec_d dec_m
  dec_s equinox mag keyval {comment *}` — literally "Example 1" in the
  format doc). A single generic line-consumption engine
  (`parse_data_line`) drives both cases by walking a field spec
  left-to-right over the line. This also unifies the format doc's
  backward-compatible "bare magnitude after equinox" rule with the general
  `mag` field type, since both reduce to "peek the next token, consume it
  only if numeric."
- **RA/Dec triplets are pseudo-fields.** `ra_h ra_m ra_s`, `dec_d dec_m
  dec_s`, and the `ra_hms`/`ra_dms`/`dec_dms` shorthands are collapsed at
  directive-compile time (`compile_data_directive` →
  `_collapse_angle_fields`) into single `__ra__`/`__dec__` pseudo-fields.
  Their consumption — 1 vs. 3 whitespace tokens, depending on whether a
  decimal point or colon appears in the data — is inherently coupled and
  can't be decided field-by-field.
- **Output type.** Each data line becomes a `StarlistEntry` dataclass
  (name, `Angle` ra/dec, equinox, `'fk4'`/`'fk5'` frame, proper motion, a
  `mag` dict keyed by band letter, priority, comment, and an `extra` dict
  for anything unrecognized), with a `.coord` property that builds an
  `astropy.coordinates.SkyCoord` in the correct frame/equinox. A whole
  file's worth of entries is then combined into a single
  `astropy.table.Table` by `entries_to_table` (see below), which
  `parse_starlist` now calls internally — see "Unified-frame output."
- **Unrecognized `!Data` field names are captured, not rejected.** They
  land in `entry.extra` rather than raising. This matters in practice: the
  format doc's own "Example 4" uses a field name `epoch` where the
  surrounding text elsewhere insists on `equinox` — almost certainly a
  documentation typo. Being lenient here means a file with an unexpected
  field name doesn't hard-crash the parser; it still correctly fails
  validation afterward if `equinox` itself was never actually supplied.

## Verification

Checked against:

- All five decimal/colon-equivalent RA/Dec forms from section I of the
  format doc (`obj1a`–`obj1e`), confirming they agree to ~1e-7 deg.
- The worked `!Data` examples in section III (fixed-width name field +
  `ra_hms`/`dec_dms` shorthand; magnitude before a literal equinox default;
  a fixed-width `skip` field) — including the `epoch`/`equinox`
  inconsistency in Example 4, which now surfaces as a clear `ValueError`
  rather than a silent misparse.
- The real `nickel_focus/data/point_focus.txt` file, end-to-end through
  `parse_starlist`.
- The full test suite (157 tests as of the latest change, 37 of them in
  `test_starlist.py`) passes.

## Unified-frame output: `entries_to_table` / `parse_starlist`

Added to resolve the second follow-up noted below. `entries_to_table(entries,
frame='icrs', **frame_attrs)` combines a list of `StarlistEntry` objects
into a single `astropy.table.Table`, sidestepping the FK4/FK5-mixing
problem by transforming each entry's `.coord` into the target frame
*individually* (via `SkyCoord.transform_to`) before stacking the results,
rather than trying to concatenate `SkyCoord`s in mismatched frames.

- `frame` defaults to `'icrs'` and can be any name recognized by
  `astropy.coordinates.frame_transform_graph` (e.g. `'fk5'`, `'fk4'`), or
  an already-constructed frame instance.
- Extra keyword arguments (e.g. `equinox='J2000.0'`) become frame
  attributes used to build the target frame when `frame` is given as a
  name; they're silently ignored if the named frame doesn't accept them —
  which is why the default `'icrs'` case satisfies "ICRS, J2000.0" without
  needing an equinox at all (ICRS has no free equinox parameter; it's
  already effectively J2000).
- The returned table carries `ra`/`dec` (deg, in `frame`) plus each
  entry's other fields (`pm_ra`, `pm_dec`, `pm_epoch`, `priority`,
  `comment`, `mag`, `extra`) and `src_equinox`/`src_frame` recording the
  *original* per-entry equinox/frame before transformation, so the
  original starlist metadata isn't lost.
- Passing an empty list returns an empty `Table` with the right columns
  and units rather than raising — needed because `parse_starlist` (below)
  now always builds a table, even for an all-comment/blank file.

Verified against `astropy`'s own `SkyCoord.icrs`/`.transform_to()` as an
independent oracle, on both a single entry and a pair of entries sharing
identical sexagesimal RA/Dec but different (1950/FK4 vs. 2000/FK5)
equinoxes — confirming they land at distinct ICRS positions once
precessed, rather than being (incorrectly) treated as identical.

**`parse_starlist` now returns a `Table` directly**, rather than a list of
`StarlistEntry`. It gained the same `frame`/`**frame_attrs` parameters as
`entries_to_table` and simply forwards them:
`entries_to_table(entries, frame=frame, **frame_attrs)` on the entries it
collects internally, right before returning. Callers no longer need to
chain the two functions themselves; `parse_data_line`/`compile_data_directive`
/etc. are still exposed for anyone who wants a raw `StarlistEntry` list
(e.g. the test suite builds one directly with a small `_entries()` helper
to test `entries_to_table` in isolation).

## Style refactor (module cleanup pass)

A follow-up pass tightened up `starlist.py` itself, with no behavior
changes beyond the `parse_starlist` return type above:

- **`from astropy import units`** replaces `import astropy.units as u`
  throughout (and in the tests), spelling out `units.deg` /
  `units.hourangle` instead of `u.deg` / `u.hourangle`.
- **Fewer module-level constants.** `DEFAULT_COMMENT_PATTERNS`,
  `DEFAULT_DATA_FIELDS`, and several single-use compiled regexes
  (`_FLOAT_RE`, `_EQUINOX_RE`, `_FIXED_WIDTH_RE`, `_DIRECTIVE_RE`) and
  lookup dicts (`_RA_PRIMARY_UNITS`, `_DEC_PRIMARY_UNITS`) were each
  referenced from exactly one function, so they're now defined locally
  right where they're used instead of cluttering module scope. `_TOKEN_RE`
  and `_KEYVAL_RE` stayed at module level since they're genuinely shared
  across more than one function.
- **`match`/`case` instead of `if`/`elif`/`else` chains** in
  `_parse_keyval_token` (dispatching on the keyword name), `_resolve_generic`
  (dispatching on `fieldformat` kind, using a guard clause with a walrus
  assignment for the `%N` fixed-width case), `_collapse_angle_fields`
  (dispatching on field name), `parse_data_line`'s per-field loop (the
  main `__ra__`/`__dec__`/`name`/`equinox`/... dispatch), and the sign- and
  B/J-prefix-detection branches in `_consume_angle` / `_parse_equinox`.
- **No truthy/falsy checks on non-boolean values.** Every conditional now
  spells out exactly what it's testing instead of relying on Python's
  truthiness rules: `directive is not None and directive.strip() != ''`
  rather than `directive and directive.strip()`; `if m is None:` rather
  than `if not m:` for a regex `Match`/`None` result (including inside a
  `match`/`case` guard, via `(m := re.match(...)) is not None`);
  `entry.name == ''` rather than `not entry.name`; `len(entries) == 0`
  rather than `not entries`. `and`/`or` between genuinely boolean
  expressions (comparisons, `isinstance`, `.startswith()`, a parameter
  documented/typed as `bool`) were left alone, since those already are
  explicit boolean tests.
- **No imports inside function bodies.** `test_starlist.py` had a
  `from importlib import resources` local to one test function; moved to
  the top of the file with the rest of the imports. (One exception exists
  in this codebase: a script's `main()`/CLI entry point, as in
  `scripts/focus.py`, may still defer heavy imports for `--help` speed —
  not applicable to this module.)
- **Line wrapping**: a def signature, call, or conditional that exceeds
  99 characters now wraps with the opening parenthesis alone on its own
  line and the contents at a flat 4-space indent, closing parenthesis on
  its own line back at the statement's indent — not the previous "hanging
  indent aligned under the opening delimiter" style. Anything that
  actually fits on one line at 99 characters or fewer was joined back
  onto one line rather than left wrapped (checked in both `starlist.py`
  and `test_starlist.py`). `compile_data_directive`'s wrapped ternary
  (`x if cond else y \` with a backslash continuation) became a plain
  `if`/`else` statement instead, avoiding the backslash entirely.
- **Full docstrings on every function, including private ones.** Every
  function in `starlist.py` — `_consume_angle`, `_parse_equinox`,
  `_parse_keyval_token`, `_resolve_generic`, `_parse_field_list`,
  `_collapse_angle_fields`, `is_comment_line`, and the `_LineCursor`
  methods — now has a numpy-style docstring explaining each argument,
  matching the level of detail already used for the public functions.
  One exception: `_LineCursor.__init__` itself stays undocumented; its
  `line` parameter is documented in the *class* docstring instead. The
  `StarlistEntry` dataclass followed the same rule — its per-field
  attribute docstrings were consolidated into one `Parameters` section
  on the class docstring, since a dataclass's fields are its (implicit)
  constructor's parameters.
- **Assert messages in the test suite.** Every `assert` in
  `test_starlist.py` now carries a message stating the behavior being
  checked (e.g. `'equinox 2000.0 with no B/J prefix should default to
  FK5'`), so a failure is legible from the pytest output alone. A
  `pytest.raises(...)` block with no bare `assert` inside doesn't need
  one.

## Open items / possible follow-ups

- **Not yet wired into `slew.py`.** `slew.py`'s `find_nearest_target` still
  reads `point_focus.txt` via `np.loadtxt` and builds a `SkyCoord` array
  manually. It could be switched to `starlist.parse_starlist`, which now
  returns a ready-to-use `Table` directly and would add support for
  comments, custom column order, and proper motion/priority metadata in
  that catalog — but this hasn't been done, pending a decision on whether
  `find_nearest_target`'s string-array-based API should change.
