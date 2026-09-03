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
- **`nickel_focus/pkg/version.py` regeneration without `pip` — done.**
  `tools/write_version.py` answers the "does such a tool exist" question:
  yes — `setuptools_scm` (already a `[build-system] requires` dependency,
  now also added to the `dev` extra in `pyproject.toml`) exposes this
  independent of any `pip install`/build step, via its `get_version()` API
  rather than its bare CLI (verified: `python -m setuptools_scm` alone only
  *prints* the version, it does **not** write `version_file` — that only
  happens through the API, called with the file's configured
  `version_file` explicitly). The script reads
  `[tool.setuptools_scm] version_file` straight out of `pyproject.toml`
  (via stdlib `tomllib`) so the path can't drift out of sync, and warns
  (without failing) if the computed version isn't a clean tag, since every
  CVS-synced version is expected to be one. Verified end-to-end against
  this repo. Note `nickel_focus/pkg/version.py` is `.gitignore`d (never
  tracked by git in the first place — it's a pure build artifact), so
  regenerating it is always safe to do freely. This becomes a step the
  (not-yet-written) CVS sync script runs right before mirroring
  `nickel_focus/` into CVS.
- **`point_focus.txt` lookup — done.** `nickel_focus/slew.py`'s
  `find_nearest_target` now locates it via
  `Path(__file__).resolve().parent / 'data' / 'point_focus.txt'` instead of
  `importlib.resources.files('nickel_focus') / ...`, since the kroot/CVS
  deployment is a plain files-on-`PYTHONPATH` layout rather than a proper
  installed distribution with package metadata that `importlib.resources`
  depends on. Verified against the full non-GUI test suite (89 passed).
  `test_starlist.py`'s own unrelated `resources.files(...)` call (used just
  to locate a file for a `parse_starlist` test) was left alone — it only
  ever runs in a normal pip/pytest context.
- **`tools/deploy.sh` drafted** (see its own section below) — the CVS-sync
  step inside it (step 3) is still nominal placeholder commands, pending
  the CVS-side specifics ("details to come later" per your note).
- You mentioned reassessing which files currently in `nickel_focus/` should
  move to better isolate what is/isn't synced — the table above reflects
  the directory layout *as it exists today*; expect it to need updating
  once that assessment is done.

## `tools/deploy.sh`

A draft end-to-end deployment script, per your requested workflow: (1)
update the deployment-side git checkout to the latest `main`, (2)
regenerate `version.py` (via `tools/write_version.py`), (3) sync
`nickel_focus/` into the CVS working copy, (4) `make install` from there.
**Not tested end-to-end** — can't be, on this machine — treat it as a
starting point to finish with Brad/Will, not a ready-to-run tool.

Design choices made, flagged here in case they need revisiting:

- **`GIT_CHECKOUT`/`CVS_CHECKOUT`** are placeholder paths (env-var
  overridable, default to obvious `/TBD/...` strings) — real values depend
  on how the deployment machine is actually laid out.
- **Step 1 uses `git fetch` + `git reset --hard origin/main`** (forcibly
  discarding any local changes on that checkout) rather than a plain
  `git pull`, per your choice — safer/more reproducible for a checkout
  nobody should be editing directly, at the cost of being destructive if
  that assumption is ever wrong. It prints a warning (via `git status
  --short`) before doing so, so the discarded changes are at least visible
  in the log.
- **Tag enforcement is a warning, not a hard failure.** The script checks
  `git describe --tags --exact-match HEAD` and warns (but continues) if
  HEAD isn't exactly a tagged commit, matching the same
  warn-don't-block choice already made in `tools/write_version.py`. If you
  and Brad/Will decide CVS syncs must never happen from an untagged
  commit, this should become a hard `exit 1` instead.
- **Step 3's commands are nominal, not final** — mirror `nickel_focus/` +
  the top `Makefile` into `$CVS_CHECKOUT` via `rsync --delete`, restore the
  executable bit on `*.sin` files (CVS doesn't track POSIX permissions the
  way git does), then use the standard `cvs -qn update` dry-run idiom
  (`?` = untracked-by-CVS file on disk → `cvs add`; `!` = CVS-tracked file
  missing from disk → `cvs remove`) before a single `cvs commit`. This is a
  reasonable, standard-idiom starting point, not a verified-correct
  mechanism — expected to be directly edited once the real CVS conventions
  are confirmed.
- **`set -euxo pipefail`** (note the `-x`) — every command is echoed as
  it runs, at your request, since this can't be tested here and will need
  to be debugged live by whoever runs it first.

## Verification (once buildable against a real `$(LROOT)`)

- Confirm each directory's files actually land under the expected
  `LIBSUB`/`BINSUB` path and the installed `nickel_focus` tree is still a
  valid, importable package (same nested structure as the source tree).
- Confirm the `.sin` → installed-binary step preserves the literal
  `#! @KPYTHON@` line pre-substitution, and that post-substitution each
  script actually runs and reaches its `entry_point()`.
- Confirm none of this affects `pip install .` / the existing pytest/tox
  suite, which remain the dev/CI path.
