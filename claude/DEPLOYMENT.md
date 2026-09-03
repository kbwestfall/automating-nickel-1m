# Deployment: kroot Makefiles and CVS sync

This is the working document for actually building the kroot/CVS deployment
path, superseding the more speculative `claude/DEPLOYMENT_PLAN.md` (kept for
history/context) now that several open questions there have been answered.
This doc tracks concrete decisions as they're made, starting with the
Makefile tree.

## Decisions so far

- **The repo root now has its own top Makefile** (`DIRS = nickel_focus`),
  matching Brad's original top-Makefile pattern directly — Brad's APF
  reference example has been renamed to `Makefile.example` for reference
  rather than being adapted in place. Only this `/Makefile` and the
  `nickel_focus/` tree beneath it are meant for CVS; `dev/`, `claude/`,
  `.github/`, `pyproject.toml`, `.gitignore`, `Makefile.example`, etc. are
  not part of the deployed tree (see `MANIFEST.in`, below, which also keeps
  most of these out of the pip sdist/wheel for the same reason).
- **The Makefile tree mirrors the directory tree**, and this is done
  entirely by hand via explicit `DIRS =` directives at each level — there's
  no automatic recursive-copy behavior. `nickel_focus/Makefile` sets
  `DIRS = config data gui pkg scripts tests`; `gui/Makefile` further sets
  `DIRS = model views`; `tests/Makefile` sets `DIRS = gui`.
- **Data-like files (`.toml` configs, `point_focus.txt`) are declared via
  `LIBFILES`/`LIBSUB`, not `DATAFILES`/`DATASUB`.** Brad's `DATA*`
  directives install to a separate directory tree from the code; the
  preference for now is to keep these files alongside the code they belong
  to (so relative-path lookups from within the package keep working the
  same way pre- and post-install). `DATAFILES`/`DATASUB` are unused in this
  package for now.
- **`focus_data.ecsv` and `data/starlistFormat.html` are excluded** from
  the deployed tree — the former looks like a runtime example/output
  rather than static reference data, the latter a dev-facing reference doc.
  Only `data/point_focus.txt` is deployed.
- **Directories whose entire `LIBFILES` list is "every `.py` file here" use
  a wildcard** (`LIBFILES = $(wildcard *.py)`) instead of an explicit
  enumeration, mirroring the existing `SINFILES = $(wildcard *.sin)`
  pattern from Brad's example. This applies to `nickel_focus/Makefile`
  itself and every subdirectory Makefile except `config/` (explicit list,
  since only two of its files are meant to ship) and `data/` (explicit,
  since two of its three files are deliberately excluded).
- **`tests/` and `tests/gui/` are included** in the mirrored tree (their
  own `LIBFILES = $(wildcard *.py)`), even though the test suite isn't
  needed at runtime — kept for consistency/completeness of the mirrored
  package for now.
- **The GUI script is part of the deployment.** `nickel_focus_gui.sin` is
  installed the same way as the other two.
- **`.sin` scripts live in their own `nickel_focus/sin/` directory**, not
  in `scripts/` — kept separate from the importable modules they wrap, so
  `scripts/` stays a plain, `.sin`-free library directory. Only the
  `@KPYTHON@` token is needed for now (no other `@VAR@` tokens, no
  `@@VAR@@` data-path tokens).
- **`sin/Makefile`'s `BINFILES = $(FILES)`** (derived from
  `SINFILES = $(wildcard *.sin)` / `FILES = $(SINFILES:%.sin=%)`) is kept
  as a deliberate choice, not an oversight — programmatically stripping
  `.sin` from every file in the directory is preferred over hand-listing
  the names again. If this turns out to not play well with the real
  `Make.rules`, revisit then.
- **Pip-packaging cleanliness is handled via a root `MANIFEST.in`**
  (`global-exclude Makefile`, `global-exclude *.sin`, `prune .github/`,
  `prune dev/`, `prune claude/`, `exclude .gitignore`,
  `exclude Makefile.example`) rather than a `pyproject.toml`
  `exclude-package-data` table. This keeps the kroot-only build artifacts
  and non-package directories out of the sdist/wheel that `pip install .`
  produces, without touching `pyproject.toml` itself.

## Makefile tree (as written)

| Directory | Role | Installs |
|---|---|---|
| `/` (repo root) | router | `DIRS = nickel_focus` |
| `nickel_focus/` | router + leaf | `LIBFILES = $(wildcard *.py)` → `LIBSUB = nickel_focus/`; `DIRS = config data gui pkg scripts sin tests` |
| `nickel_focus/config/` | leaf | `LIBFILES = nickucam.toml nscicam.toml` → `nickel_focus/config/` |
| `nickel_focus/data/` | leaf | `LIBFILES = point_focus.txt` → `nickel_focus/data/` |
| `nickel_focus/gui/` | router + leaf | `LIBFILES = $(wildcard *.py)` → `nickel_focus/gui/`; `DIRS = model views` |
| `nickel_focus/gui/model/` | leaf | `LIBFILES = $(wildcard *.py)` → `nickel_focus/gui/model/` |
| `nickel_focus/gui/views/` | leaf | `LIBFILES = $(wildcard *.py)` → `nickel_focus/gui/views/` |
| `nickel_focus/pkg/` | leaf | `LIBFILES = $(wildcard *.py)` → `nickel_focus/pkg/` (includes `version.py`, see below) |
| `nickel_focus/scripts/` | leaf | `LIBFILES = $(wildcard *.py)` → `nickel_focus/scripts/` |
| `nickel_focus/sin/` | leaf | `SINFILES = $(wildcard *.sin)`, `FILES`/`BINFILES` derived from it → `BINSUB = nickel` |
| `nickel_focus/tests/` | router + leaf | `LIBFILES = $(wildcard *.py)` → `nickel_focus/tests/`; `DIRS = gui` |
| `nickel_focus/tests/gui/` | leaf | `LIBFILES = $(wildcard *.py)` → `nickel_focus/tests/gui/` |

Three `.sin` files live in `nickel_focus/sin/`, each a thin wrapper reusing
the existing `ScriptBase.entry_point()` machinery (no logic duplicated):
`nickel_focus.sin` → `NickelFocus`, `nickel_focus_gui.sin` →
`NickelFocusGUI`, `nickel_slew_to_nearest.sin` → `NickelSlewToNearest`.

## Still open / unverified

- **Combining `DIRS =` with local `LIBFILES =` in the same Makefile**
  (used at `nickel_focus/`, `gui/`, and `tests/`, which all have both
  local files and subdirectories) is an assumption — Brad's one example
  Makefile only showed a pure-leaf directory, and the email's "top"
  Makefile only showed pure `DIRS`-routing with no local files. The
  repo-root `/Makefile` is now a clean pure-`DIRS` router (`DIRS =
  nickel_focus`, no local files), so it doesn't exercise this combination
  either — only `nickel_focus/`, `gui/`, and `tests/` do. Needs confirming
  against the real `Make.rules`, or just trying a real build against
  `$(LROOT)` and seeing whether both directives are honored together.
- **`BINSUB = nickel`** is a placeholder pending naming-convention
  confirmation (Brad's APF example used `master` for its `Main` script).
  Easy to rename later.
- **`nickel_focus/pkg/version.py` regeneration without `pip`.** Answering
  the "does such a tool exist" question: yes — `setuptools_scm` (already a
  `[build-system] requires` dependency) exposes this independent of any
  `pip install`/build step. Since `pyproject.toml` already configures
  `[tool.setuptools_scm] version_file = "nickel_focus/pkg/version.py"`,
  running `python -m setuptools_scm` from the repo root (with
  `setuptools_scm` installed in whatever *dev* Python runs the sync step —
  not `kpython`) recomputes the version from the current git tag and
  rewrites that file, with no other build/install action attached. Given
  every CVS-synced version will be a tagged commit, the plan is to run
  this as a step in the (not-yet-written) CVS sync script, right before
  mirroring `nickel_focus/` into CVS, so the exported `version.py` always
  reflects the tag rather than a stale/dev version string.
- **`point_focus.txt` lookup still uses `importlib.resources`**
  (`nickel_focus/slew.py:141`, `resources.files('nickel_focus') / 'data' /
  'point_focus.txt'`). You noted this will likely need to change to a
  `__file__`-relative lookup instead, since the kroot/CVS deployment is a
  plain files-on-`PYTHONPATH` layout rather than a proper installed
  distribution with package metadata that `importlib.resources` depends on.
  Not changed yet — flagged as a follow-up once you've finished assessing
  which files move where.
- **CVS sync script** itself: still blocked on the CVS-side specifics
  ("details to come later" per your note) — mechanism will be written from
  scratch once those are available. It will need to run the `version.py`
  regeneration step above before mirroring `nickel_focus/`.
- You mentioned reassessing which files currently in `nickel_focus/` should
  move to better isolate what is/isn't synced — the table above reflects
  the directory layout *as it exists today*; expect it to need updating
  once that assessment is done.

## Verification (once buildable against a real `$(LROOT)`)

- Confirm each directory's files actually land under the expected
  `LIBSUB`/`BINSUB` path and the installed `nickel_focus` tree is still a
  valid, importable package (same nested structure as the source tree).
- Confirm the `.sin` → installed-binary step preserves the literal
  `#! @KPYTHON@` line pre-substitution, and that post-substitution each
  script actually runs and reaches its `entry_point()`.
- Confirm none of this affects `pip install .` / the existing pytest/tox
  suite, which remain the dev/CI path.
