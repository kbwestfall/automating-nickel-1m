#!/usr/bin/env python3
"""
Reflow every docstring in one or more Python files to a target width
(80 columns by default), mimicking vim's ``gqip``.

Numpydoc structure is respected: section headers and their underlines
(e.g. ``Parameters`` / ``----------``), field-list entries, and bulleted
list items (``- ``/``* ``, whether separated by blank lines or adjacent)
are each reflowed independently at their own indentation, never merged
with a neighboring one. Sphinx inline-markup spans (e.g.
``:meth:`~a.b.C.d```) and other unbroken tokens (URLs, etc.) are never
split, since they contain no whitespace and ``textwrap`` only breaks on
whitespace.

By default this only reports what would change; pass ``--apply`` to
write the changes.
"""
import argparse
import ast
import textwrap
from pathlib import Path


def is_underline(line, prev_text_len):
    """
    Whether `line` looks like a numpydoc section-header underline (e.g.
    ``----------`` under ``Parameters``).

    Parameters
    ----------
    line : str
        The candidate underline text.
    prev_text_len : int
        Length of the preceding (stripped) header line, which a real
        underline is expected to roughly match.

    Returns
    -------
    bool
        True if `line` is a run of a single repeated underline-style
        character, long enough to plausibly underline `prev_text_len`
        characters of header text.
    """
    stripped = line.strip()
    if not stripped:
        return False
    chars = set(stripped)
    return (
        len(chars) == 1 and chars <= set('-=~^"\'')
        and len(stripped) >= max(3, prev_text_len - 2)
    )


def reflow_single_line(text, indent, width):
    """
    Reflow a docstring body that has no internal newline at all -- the
    compact one-line style, e.g. a docstring reading just ``Build a Foo
    from bar.`` with no line break between its opening and closing
    quotes.

    Parameters
    ----------
    text : str
        The docstring's text (no surrounding triple quotes).
    indent : int
        The column the docstring's own ``\"\"\"`` starts at in the
        source.
    width : int
        Target total column width (indentation included).

    Returns
    -------
    str
        The new body text (no surrounding triple quotes). Stays a
        single physical line if it still fits with the quotes and
        indent included; otherwise switches to this codebase's
        multi-line convention (opening quotes alone, indented body,
        closing quotes alone at the same indent).
    """
    indent_str = ' ' * indent
    if len(f'{indent_str}"""{text}"""') <= width:
        return text
    wrapped = textwrap.fill(
        text, width=width, initial_indent=indent_str, subsequent_indent=indent_str,
        break_long_words=False, break_on_hyphens=False, fix_sentence_endings=True,
    )
    return '\n' + wrapped + '\n' + indent_str


def reflow_body(body, indent, width=80):
    """
    Reflow a docstring body (the text between the triple quotes, NOT
    including them), preserving numpydoc structure.

    Parameters
    ----------
    body : str
        The raw docstring content, exactly as `ast` gives it via the
        string constant's value (i.e. no surrounding triple quotes).
    indent : int
        The column the docstring's own ``\"\"\"`` starts at in the
        source. The docstring text itself carries no indentation
        information for a single-line docstring, and using the true
        source column is more robust than guessing from the body's own
        lines even when it is multi-line.
    width : int
        Target total column width (indentation included).

    Returns
    -------
    str
        The reflowed body, in the same style (single-line or
        multi-line) as best fits the target width.
    """
    if '\n' not in body:
        return reflow_single_line(body, indent, width)

    lines = body.split('\n')

    # The body typically ends with a whitespace-only line that exists only
    # to indent the closing triple-quote -- preserve it verbatim rather
    # than let it fall through as a blank-paragraph separator.
    trailing_indent_line = None
    if lines and lines[-1] != '' and lines[-1].strip() == '':
        trailing_indent_line = lines[-1]
        lines = lines[:-1]

    def flush_paragraph(buf, indent, out):
        """Reflow one paragraph's worth of already-stripped lines in `buf`
        (all at the same `indent`) and append the result to `out`."""
        if not buf:
            return
        text = ' '.join(s.strip() for s in buf)
        pad = ' ' * indent
        wrapped = textwrap.fill(
            text, width=width, initial_indent=pad, subsequent_indent=pad,
            break_long_words=False, break_on_hyphens=False,
            fix_sentence_endings=True,
        )
        out.extend(wrapped.split('\n'))

    out = []
    i = 0
    n = len(lines)
    while i < n:
        line = lines[i]
        if not line.strip():
            out.append('')
            i += 1
            continue

        indent = len(line) - len(line.lstrip(' '))
        stripped = line.strip()

        # Section header + underline (e.g. "Parameters" / "----------").
        if i + 1 < n and is_underline(lines[i + 1], len(stripped)):
            out.append(line)
            out.append(lines[i + 1])
            i += 2
            continue

        # Bulleted list item: "- " or "* " marker: always its own block.
        if stripped[:2] in ('- ', '* '):
            marker_indent = indent
            cont_indent = marker_indent + 2
            buf = [stripped]
            i += 1
            while i < n and lines[i].strip():
                li = len(lines[i]) - len(lines[i].lstrip(' '))
                nstripped = lines[i].strip()
                if li <= marker_indent and nstripped[:2] in ('- ', '* '):
                    break
                if li < cont_indent:
                    break
                buf.append(nstripped)
                i += 1
            text = ' '.join(buf)
            wrapped = textwrap.fill(
                text, width=width, initial_indent=' ' * marker_indent,
                subsequent_indent=' ' * cont_indent,
                break_long_words=False, break_on_hyphens=False,
                fix_sentence_endings=True,
            )
            out.extend(wrapped.split('\n'))
            continue

        # Otherwise: accumulate a run of contiguous non-blank lines at
        # exactly this indentation into one reflow-able paragraph.
        buf = [stripped]
        cur_indent = indent
        i += 1
        while i < n and lines[i].strip():
            li = len(lines[i]) - len(lines[i].lstrip(' '))
            nstripped = lines[i].strip()
            if li != cur_indent:
                break
            if nstripped[:2] in ('- ', '* '):
                break
            if i + 1 < n and is_underline(lines[i + 1], len(nstripped)):
                break
            buf.append(nstripped)
            i += 1
        flush_paragraph(buf, cur_indent, out)

    if trailing_indent_line is not None:
        out.append(trailing_indent_line)

    return '\n'.join(out)


def find_docstring_nodes(tree):
    """
    Find every real docstring `ast.Constant` node in a parsed module:
    the module's own docstring, plus one for every class/function/method
    whose first statement is a bare string literal.

    Parameters
    ----------
    tree : ast.Module
        The parsed source, from `ast.parse`.

    Returns
    -------
    list of ast.Constant
        Each node's `.value` is the docstring text; its `.lineno`,
        `.col_offset`, `.end_lineno`, and `.end_col_offset` locate it in
        the source.
    """
    nodes = []

    def maybe_add(container):
        body = getattr(container, 'body', None)
        if not body:
            return
        first = body[0]
        if (isinstance(first, ast.Expr) and isinstance(first.value, ast.Constant)
                and isinstance(first.value.value, str)):
            nodes.append(first.value)

    maybe_add(tree)
    for node in ast.walk(tree):
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            maybe_add(node)
    return nodes


def process_file(path, width=80, dry_run=True):
    """
    Reflow every docstring in the Python file at `path`.

    Parameters
    ----------
    path : str or pathlib.Path
        The file to process.
    width : int
        Target total column width (indentation included).
    dry_run : bool
        If True (the default), only print what would change; if False,
        write the changes to `path`.

    Returns
    -------
    int
        The number of docstrings changed (or that would change, in a
        dry run).
    """
    with open(path, 'r') as f:
        src = f.read()
    tree = ast.parse(src)

    replacements = []
    for node in find_docstring_nodes(tree):
        old_text = ast.get_source_segment(src, node)
        if old_text is None or not (old_text.startswith('"""') or old_text.startswith('r"""')):
            continue
        prefix = 'r' if old_text.startswith('r"""') else ''
        inner = old_text[len(prefix) + 3:-3]
        new_inner = reflow_body(inner, node.col_offset, width=width)
        if new_inner == inner:
            continue
        new_text = f'{prefix}"""{new_inner}"""'
        replacements.append((
            node.lineno, node.col_offset, node.end_lineno, node.end_col_offset,
            old_text, new_text,
        ))

    if dry_run:
        for (sl, sc, el, ec, old, new) in replacements:
            print(f'--- {path}, lines {sl}-{el} ---')
            print('BEFORE:')
            print(old)
            print('AFTER:')
            print(new)
            print()
        print(f'{len(replacements)} docstring(s) would change in {path}')
        return len(replacements)

    # Apply replacements from bottom to top so earlier offsets stay valid.
    src_lines = src.split('\n')
    by_position = sorted(replacements, key=lambda r: (r[0], r[1]), reverse=True)
    for (sl, sc, el, ec, old, new) in by_position:
        if sl == el:
            line = src_lines[sl - 1]
            src_lines[sl - 1] = line[:sc] + new + line[ec:]
        else:
            first_line = src_lines[sl - 1][:sc] + new
            src_lines[sl - 1:el] = first_line.split('\n')
    with open(path, 'w') as f:
        f.write('\n'.join(src_lines))
    print(f'{len(replacements)} docstring(s) changed in {path}')
    return len(replacements)


def iter_python_files(paths):
    """
    Expand a list of file/directory arguments into the ``.py`` files
    they name, recursing into directories.

    Parameters
    ----------
    paths : list of str
        Files and/or directories, as given on the command line.

    Returns
    -------
    list of pathlib.Path
        Every ``.py`` file found, in sorted order, with `__pycache__`
        contents excluded.
    """
    files = []
    for raw in paths:
        p = Path(raw)
        if p.is_dir():
            files.extend(
                f for f in p.rglob('*.py') if '__pycache__' not in f.parts
            )
        else:
            files.append(p)
    return sorted(set(files))


def main():
    """Parse command-line arguments and reflow the requested files."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        'paths', nargs='+',
        help='Python files and/or directories (searched recursively) to reflow.'
    )
    parser.add_argument(
        '--width', type=int, default=80,
        help='Target column width for reflowed docstring text.'
    )
    parser.add_argument(
        '--apply', action='store_true',
        help='Write changes to disk. Without this, only reports what would change.'
    )
    args = parser.parse_args()

    total = 0
    for path in iter_python_files(args.paths):
        total += process_file(path, width=args.width, dry_run=not args.apply)
    if not args.apply:
        print(f'\n{total} docstring(s) total would change (dry run; pass --apply to write).')
    else:
        print(f'\n{total} docstring(s) total changed.')


if __name__ == '__main__':
    main()
