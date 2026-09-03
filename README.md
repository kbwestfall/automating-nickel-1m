# nickel_focus

Acquisition and focusing tools for the Nickel 1-m Telescope at Lick
Observatory.

Scott Hakoda
Utah Tech University

Site: University of California Observatories, Santa Cruz, California
Mentors: Kyle Westfall, Will Deich

## Overview

Lick Observatory's Nickel 1-m Telescope, located atop Mount Hamilton, is
used by a wide variety of observers, including faculty, staff, and
students both within and outside the University of California
Observatories. At the start of each observing session, observers must
perform a set of startup and calibration procedures, including focusing
and adjusting the telescope's pointing. Focusing in particular has
historically been a manual, ~10-minute-per-iteration process, repeated
several times per night — tedious, error-prone, and a meaningful source of
lost telescope time.

This package automates that process using the Keck Task Library (KTL) API,
the control interface for the Nickel's mechanisms (secondary focus,
pointing) and its science camera. Given the telescope's current position,
it selects and slews to the nearest known focus or pointing star, takes a
series of exposures across a range of focus values, measures each
exposure's image quality (FWHM), and fits a curve to the results to
determine the optimal focus. A Qt-based GUI is also available for
interactive use.

## Installation

`nickel_focus` requires Python 3.12–3.14. The `ktl` package (Keck's
telescope-control middleware) is not pip-installable — it's provided
separately by the observatory's `kpython` environment. `nickel_focus` is
designed to degrade gracefully without it: importing the package still
works, and code paths that don't need a live telescope connection (e.g.
re-analyzing an already-recorded focus sequence) still run.

For most uses — development, testing on a laptop, or reprocessing archived
data — install directly from a checkout of this repository:

```console
pip install .
```

Optional extras add functionality on top of the base install:

| Extra | Adds |
|---|---|
| `gui` | The Qt-based focusing GUI (`nickel_focus_gui`), via PySide6 |
| `test` | The pytest-based test suite |
| `gui_test` | GUI-specific test dependencies (`pytest-qt`) |
| `dev` | All of the above, plus `setuptools_scm` for the maintainer tools under `tools/` |

For example, to install with GUI support:

```console
pip install ".[gui]"
```

or, for a full development setup:

```console
pip install ".[dev]"
```

This installs the three command-line scripts described below directly
onto your `PATH` (via the environment's usual `bin`/`Scripts` directory).
Note this is separate from how the package reaches the telescope itself —
see **Maintenance** below.

## Scripts

Three command-line scripts are installed. All three share a common set of
logging options: `-v`/`--verbosity` (0 = warnings/errors only, 1 = default,
adds informational messages, 2 = adds debug messages), `--log_file` (write
a log file; use `default` for an automatically-named one), and
`--log_level` (verbosity for the log file specifically, if different from
the console). Pass `--help` to any of them for the full, current option
list.

### `nickel_focus`

The main focus-finding driver. Given a starting focus value and a step
size, it either follows an automated curve-fitting search or steps
through a fixed grid, taking an exposure and measuring its image quality
at each step, then fits a quadratic to focus-vs-FWHM to find the optimum.

```console
# Automated search starting at focus 350, step size 5
nickel_focus 350 5

# Fixed grid from focus 350 to 400 in steps of 5
nickel_focus 350 5 400

# Fixed grid with an explicit number of steps instead of an end value
nickel_focus 350 5 --nstep 10

# Save the measured focus curve to a file for later reference
nickel_focus 350 5 --ofile focus_data.ecsv
```

Useful options include `--exptime`/`-t` (exposure time in seconds),
`--binning`/`-b` (`1,1`, `2,2`, or `4,4`), `--speed`/`-s` (`Slow` or
`Fast`), and `--no-plot` (disable the live focus-curve plot). A sequence
can also be re-analyzed from previously-recorded exposures, without a live
`ktl` connection, using `--obsnum` together with `--datadir`, `--prefix`,
and `--suffix` to locate the archived files. (`--refit`/`--omit`, for
refitting a saved curve after excluding bad points, are accepted but not
yet implemented.)

### `nickel_focus_gui`

Launches a Qt-based GUI for the same focusing functionality — slewing,
single exposures, grid sequences, and automated sequences — with live
plots of the current frame, per-exposure source stamps, and the evolving
focus curve. Requires the `gui` extra (PySide6). Takes no arguments beyond
the shared logging options above.

```console
nickel_focus_gui
```

### `nickel_slew_to_nearest`

Slews the telescope to the nearest known pointing or focus star, or the
nearest match in a user-supplied starlist.

```console
# Slew to the nearest star in the default pointing/focus catalog
nickel_slew_to_nearest

# Restrict to the pointing stars specifically, without actually moving
nickel_slew_to_nearest --search Pointing --dry_run

# Use a custom starlist file instead of the default catalog
nickel_slew_to_nearest --starlist my_targets.txt
```

## Maintenance

> **Notional.** The description below reflects the deployment approach as
> currently designed and is expected to change as it's finalized with
> Brad Holden and Will Deich; see `claude/DEPLOYMENT.md` in this
> repository for the current, evolving detail.

Development happens in this git repository. At the telescope, however,
`nickel_focus` runs under the observatory's own `kpython` environment and
is maintained at the system level by observatory staff — observers use
that system-level install directly, rather than a `pip`-managed
environment.

That system install goes through Lick/Keck's existing "kroot"-style build
system, whose source of truth is CVS rather than git. The plan is:

- A one-way, manually-triggered sync (`tools/deploy.sh`) mirrors this
  repository's `nickel_focus/` tree into the CVS-managed working copy that
  the kroot build actually installs from.
- A tree of kroot-convention Makefiles under `nickel_focus/` (mirroring
  its directory structure) declare which files are importable library code
  versus installed executables, following the same conventions used by
  other instruments' software at Lick.
- The three command-line scripts are installed from `.sin` template files
  under `nickel_focus/sin/`, each starting with a `#! @KPYTHON@` shebang
  that gets substituted with the real `kpython` interpreter path at
  install time — the standard UCO/Keck convention for this kind of script.
- Because no `pip` build ever happens along this path,
  `nickel_focus/pkg/version.py` (normally generated by `setuptools_scm` as
  a side effect of a `pip` build) is instead regenerated directly via
  `tools/write_version.py`, run as a step in the deploy script.

`pip install .` and the scripts it installs (see above) remain the
supported path for development, CI, and any other pip-managed use — the
kroot/CVS path above is specific to the telescope's own system
installation.
