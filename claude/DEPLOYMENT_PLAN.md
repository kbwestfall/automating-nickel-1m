# Integrate `nickel_focus` into the Lick/kroot CVS deployment system

*This is a tentative planning document, meant to be reviewed with Brad
Holden before any implementation begins. It captures what can be inferred
from the material provided so far, a concrete-but-adjustable proposal, and
— most importantly — the open questions that need answers before this can
be finalized. Most of these should be answerable by Brad; a few (flagged
below) may ultimately need input from Will Deich instead, given his role in
the KTL/`@KPYTHON@` side of things noted in `claude/PROJECT.md`.*

## Context

Development happens in this git repository; production deployment at the
telescope happens through a separate, pre-existing Lick/Keck-style (`kroot`)
build system whose source of truth is CVS, not git. This is per an email
from Brad Holden describing the Makefile conventions (an example Makefile
he sent is now at the repo root, `/Makefile`). Two things are needed:

1. **A one-way git → CVS sync**, run manually by a person, that pushes the
   current state of this git repo into the CVS tree that the kroot build
   system actually deploys from.
2. **Kroot-convention Makefiles and `.sin` scripts** for `nickel_focus`, so
   the existing top-level `kroot` build machinery (`Make.defs`/`Make.rules`,
   `@VAR@` preprocessing) can build and install it the same way it does
   every other instrument's software (e.g. APF, per Brad's example
   Makefile).

This supersedes the earlier `bin/`-directory plan (scrapped) — the
`#! @KPYTHON@` / thin-wrapper script *content* designed there is still
correct and reusable, it just needs to become `.sin` files living where the
kroot Makefile convention expects them, rather than a bespoke `bin/`
directory.

## What we know about the kroot build system

From Brad's email and the example Makefile now committed at the repo root
(`/Makefile`):

- **Top Makefile** (per-package, one per kroot module):
  ```make
  ifndef LROOT
      LROOT = /usr/local/lick
  endif
  include $(LROOT)/etc/Make.defs
  DIRS = main utils
  include $(LROOT)/etc/Make.rules
  ```
  Defines `LROOT` (default `/usr/local/lick`) if not already set, pulls in
  shared boilerplate, then just recurses into each directory named in
  `DIRS`.

- **Per-directory Makefile** (the one currently at repo root is Brad's APF
  example, not yet adapted for us) declares, per directory:
  - `LIBFILES` / `LIBSUB` — Python modules to install as importable library
    code, placed under `LIBSUB` (a subdirectory of the shared Python lib
    path).
  - `BINFILES` / `BINSUB` — installed executables; `BINSUB` puts them in a
    *subdirectory* of `bin/` specifically so they are **not** on everyone's
    `PATH` by default (APF's example: `BINFILES = Main`, `BINSUB = master`).
  - `DATASUB` / `DATAFILES` — data files and their install subdirectory
    (APF's example: `DATASUB = apf/main`).
  - `SINFILES = $(wildcard *.sin)` — every `.sin` file in the directory is
    picked up automatically; `FILES = $(SINFILES:%.sin=%)` is the resulting
    installed name (`Main.sin` → installed as `Main`).

- **`.sin` files**: by UCO/Keck convention, an executable's source is named
  `whatever.sin` with a first line like `#! @KPYTHON@`. A preprocessor
  replaces `@VAR@`-delimited tokens with real values at build time (e.g.
  `@KPYTHON@` → the local kpython install, something like
  `/opt/kroot/rel/default/bin/kpython`). This is the same shebang-token
  convention already seen (unused) in `dev/automate.py` and a couple of
  `nickel_focus/scripts/*.py` files.

- A separate **`@@VAR@@`** (double-@) convention exists for referencing data
  file install locations from within code — mechanism and exact usage not
  yet confirmed (see open questions).

## Proposed adaptation (tentative)

1. **`.sin` scripts replace the earlier `bin/` scripts.** Same thin-wrapper
   content as designed previously — reusing `nickel_focus.scripts.*`'s
   existing `ScriptBase.entry_point()` machinery, no logic duplicated:

   ```python
   #! @KPYTHON@

   from nickel_focus.scripts.focus import NickelFocus


   def main():
       NickelFocus.entry_point()


   if __name__ == "__main__":
       main()
   ```

   as `nickel_focus.sin`, `nickel_focus_gui.sin`, `nickel_slew_to_nearest.sin`
   (paralleling `NickelFocusGUI`/`nickel_focus.scripts.focus_gui` and
   `NickelSlewToNearest`/`nickel_focus.scripts.slew_to_nearest`).

2. **Adapt the Makefiles** to this repo's actual layout instead of APF's:
   - Top Makefile: replace `DIRS = main utils` with whatever subdirectory
     breakdown we settle on for this repo (see open questions — APF's
     `main`/`utils` split doesn't map cleanly onto our single nested
     `nickel_focus/` package).
   - A per-directory Makefile (tentatively living in/under `nickel_focus/`)
     declaring:
     - `LIBFILES` = every library module (`__init__.py`, `focus.py`,
       `photometry.py`, `quadratic.py`, `slew.py`, `starlist.py`,
       `pkg/*.py`, `gui/*.py`, `gui/model/*.py`, `gui/views/*.py`,
       `scripts/*.py` — 24 files today, excluding `tests/`), **assuming**
       `LIBFILES`/`LIBSUB` preserves relative subdirectory structure so the
       installed tree is still a valid, importable nested package (open
       question below — this is the single biggest structural unknown).
     - `LIBSUB` = wherever the shared kpython `PYTHONPATH` root is.
     - `SINFILES` picked up automatically from the three `.sin` files above.
     - `BINSUB` = a Nickel-specific bin subdirectory (mirroring APF's
       `BINSUB = master`), exact name TBD.
     - `DATAFILES` = `config/nscicam.toml`, `config/nickucam.toml`,
       `data/point_focus.txt`, `data/starlistFormat.html` (excluding
       `data/focus_data.ecsv`, which looks like a runtime example/output
       rather than static deployed data — to confirm).
     - `DATASUB` = a Nickel-specific data subdirectory (mirroring APF's
       `DATASUB = apf/main`), exact name TBD.

3. **`pip`/`pyproject.toml` packaging is unaffected.** `[project.scripts]`
   and a normal `pip install .` remain exactly as they are today, for
   dev-laptop and CI use (including the archive/replay `ArchiveFocusSequence`
   path that's designed to run without a live `ktl`/kpython environment).
   The kroot/CVS path becomes the **sole production deployment mechanism**,
   replacing the scrapped `bin/`-directory idea — there is no collision to
   manage since these are now two entirely separate toolchains (`pip` vs.
   kroot `make`) rather than two mechanisms fighting over the same install
   step.

4. **One-way git → CVS sync script.** A script (kept in this git repo,
   e.g. `dev/sync_to_cvs.sh` or similar — naming/location TBD), run manually
   by a human, that:
   - Exports the current git `HEAD` tree (`git archive` or equivalent) for
     whatever subset of the repo is meant to live in CVS (see open question
     on scope).
   - Reconciles that export against a local CVS working copy: `cvs add` for
     new files, `cvs remove` for files no longer present, then a single
     `cvs commit`.
   - Needs explicit handling for the executable bit on `.sin` files (CVS
     doesn't track POSIX permissions the way git does) and probably a
     dry-run/diff mode before committing, given it's a manual, human-run
     step rather than an automated pipeline.

## Open questions for Brad (critical — need answers before finalizing)

**CVS mechanics**
- Where does the CVS repository/module for this software live (`CVSROOT`,
  module name), and is there already a local CVS working-copy checkout
  convention to reuse, or does the sync script need to create one?
- Does the one-way sync mirror the *entire* git repo into CVS, or only the
  kroot-deployable subset (i.e., just what ends up under `LIBFILES`/
  `SINFILES`/`DATAFILES`, excluding `dev/`, `claude/`, `.github/`, tests,
  `pyproject.toml`/pip-only tooling)?
- Any existing precedent/script from another instrument's repo we should
  follow instead of writing this from scratch?

**Makefile / directory structure**
- Does `LIBFILES`/`Make.rules` preserve relative subdirectory paths (so a
  nested package like `nickel_focus/gui/views/main_window.py` installs to
  `<LIBSUB>/nickel_focus/gui/views/main_window.py` and stays a valid,
  importable package), or does it flatten everything into one directory?
  This determines whether `nickel_focus`'s current nested-package layout
  can be deployed as-is, or needs restructuring for the CVS side.
- What's the naming convention other instruments use for `DIRS`/`BINSUB`/
  `DATASUB` (APF uses `master`/`apf/main`) — is there a Nickel-specific
  precedent we should match (e.g. `nickel`)?
- Should the top Makefile's `DIRS` be a single directory (e.g. just
  `nickel_focus`) given this is one nested package rather than APF's
  flat `main`/`utils` split, or does Brad want a different top-level
  reorganization to match the `main`/`utils` convention more literally?

**Runtime environment**
- Are `nickel_focus`'s third-party dependencies (astropy, numpy,
  matplotlib, photutils, scipy, IPython, and PySide6 for the GUI) already
  available in the shared kpython environment, or does deployment need to
  provision them separately? This is untouched by the Makefile/LIBFILES
  mechanism, which only places our own source files.
- `nickel_focus/pkg/version.py` is currently generated by `setuptools_scm`
  during a `pip`/build step. With no `pip` involved in the kroot deployment
  path, how should this file be produced/maintained instead (static,
  manually bumped; generated by a small step in the sync script; something
  else)?
- Should the GUI script (`nickel_focus_gui.sin`, depending on PySide6) be
  part of this deployment at all, or is it dev/laptop-only for now?

**`.sin` preprocessing details** *(may need Will rather than Brad — this
overlaps the KTL/`kpython` side of things)*
- Besides `@KPYTHON@`, are there other `@VAR@` tokens we should expect to
  need (e.g. anything for locating installed data files from within a
  script)?
- What does the double-`@@VAR@@` convention specifically provide, and would
  any of our code (e.g. locating `point_focus.txt` at runtime, which today
  is just a relative/hardcoded path) need to use it?
- Confirming what's available in the shared kpython environment (its `ktl`
  package, and the third-party dependencies noted above) is really a
  question for Will, even though Brad may also know the answer.

## Verification (once the above is resolved)

- Dry-run the sync script against a scratch/throwaway CVS module first, not
  the real deployment tree.
- Build via the adapted Makefiles in a sandbox with a fake `LROOT` (or
  against the real one read-only) and confirm `SINFILES`→installed files
  keep the literal `#! @KPYTHON@` line intact, and that installed
  `LIBFILES` form a still-importable `nickel_focus` package.
- Confirm `pip install .` / existing pytest/tox suite are completely
  unaffected by any of these additions.
