# nickel_focus

> **NOTE**: This repo has been moved to
> https://github.com/UCObservatories/nickel_focus .  This fork of the
> original version of the code is now archived.

Acquisition and focusing tools for the Nickel 1-m Telescope at Lick
Observatory.

**Original Author**: Scott Hakoda (Utah Tech University; Akamai Workforce Initiative 2025)

**Site**: University of California Observatories, Santa Cruz, California

**Mentors**: Kyle Westfall, Will Deich

## Overview

Lick Observatory's Nickel 1-m Telescope, located atop Mount Hamilton, is used by
a wide variety of observers, including faculty, staff, and students both within
and outside the University of California Observatories. At the start of each
observing session, observers must confirm the telescope pointing and focus the
optics.  Focusing in particular has historically been a tedious,
~10-minute-per-iteration process that can be error-prone and a meaningful source
of lost telescope time.

This package automates that process using the Keck Task Library (KTL) API, the
control interface for the Nickel's mechanisms (secondary focus, pointing) and
its science camera. Given the telescope's current position, it selects and slews
to the nearest known focus or pointing star, takes a series of exposures across
a range of focus values, measures each exposure's image quality (FWHM), and fits
a curve to the results to determine the optimal focus. A Qt-based GUI is also
available for interactive use.

## Installation

Although provided as a python package, the code is installed and maintained by
observatory staff.  Local installation is possible but of limited use for
observers.  For completeness, we describe installation here.

`nickel_focus` requires Python 3.12–3.14. The `ktl` package (Keck's
telescope-control middleware) is not pip-installable — it's provided separately
by the observatory's `kpython` environment. `nickel_focus` is designed to
degrade gracefully without it: importing the package still works, and code paths
that don't need a live telescope connection (e.g.  re-analyzing an
already-recorded focus sequence) still run.

For most uses (development, testing on a laptop, or reprocessing archived
data), the installation process is as follows:

```console
# Create a python environment (if it doesn't already exist)
cd ~
python3 -m venv nickel
# Activate the environment
source ~/nickel/bin/activate
# Clone the repository
git clone https://github.com/kbwestfall/automating-nickel-1m.git
# Install using pip
cd automating-nickel-1m
pip install -e .
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
pip install -e ".[gui]"
```

or, for a full development setup:

```console
pip install -e ".[dev]"
```

This installs the three command-line scripts described below.  Note this
installation is separate from how the package reaches the telescope itself — see
**Maintenance** below.

## Scripts

Three command-line scripts are installed. All three share a common set of
logging options:

 - `-v`/`--verbosity`: Set the verbosity level: 0 = warnings/errors only, 1 =
   default, adds informational messages, 2 = adds debug messages.

 - `--log_file`: All logging messages are also passed to a log file; set to
   `default` for an automatically generated file name.

 - `--log_level`: Set the verbosity for the log file specifically, if different
   from the console.

 - `--help`: Provide the full list of command-line options.

### `nickel_slew_to_nearest`

Slews the telescope to the nearest known pointing or focus star, or the nearest
match in a user-supplied starlist.

Examples:

```console
# Slew to the nearest star in the default pointing/focus catalog
nickel_slew_to_nearest

# Restrict to the pointing stars specifically, without actually moving
nickel_slew_to_nearest --search Pointing --dry_run

# Use a custom starlist file instead of the default catalog
nickel_slew_to_nearest --starlist my_targets.txt
```

### `nickel_focus`

The main focus-finding driver. Given a starting focus value and a step size, it
either follows an automated curve-fitting search or steps through a fixed grid,
taking an exposure and measuring its image quality at each step, then fits a
quadratic to focus-vs-FWHM to find the optimum.

Examples:

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

Useful options include:

 -  `--exptime`/`-t`: Set the exposure time in seconds

 - `--binning`/`-b`: Set the binning (`1,1`, `2,2`, or `4,4`)

 - `--speed`/`-s`: Set the read speed (`Slow` or `Fast`)
 
 - `--no-plot`: Disable the live focus-curve plot; just report the best focus value.

A sequence can also be re-analyzed from previously-recorded exposures, without a
live `ktl` connection, using `--obsnum` together with `--datadir`, `--prefix`,
and `--suffix` to locate the archived files.

### `nickel_focus_gui`

Launch a Qt-based GUI to perform the same focusing functionality — slewing,
single exposures, grid sequences, and automated sequences — with live plots of
the current frame, per-exposure source stamps, and the evolving focus curve.
Requires the `gui` extra dependencies (PySide6). Takes no arguments beyond the
shared logging options above.

```console
nickel_focus_gui
```

## Maintenance

> **Notional.** The description below reflects the deployment approach as
> currently designed and is expected to change as it's finalized with
> Brad Holden and Will Deich; see `claude/DEPLOYMENT.md` in this
> repository for the current, evolving detail.

Development happens in this git repository.  Installation is done via

```console
pip install -e ".[dev]"
```

The package includes continuous integration tests performed via `tox` and every time changes are pushed to the GitHub repository.

Tests can be executed locally using `pytest`:

```console
cd automating-nickel-1m
pytest -W ignore
```

or via `tox` (see `tox.ini` for the versions that are tested)

```console
tox -e 3.12-test-qt
```

> **IMPORTANT**: None of the tests actually exercise the `ktl` connections
> established to manipulate the telescope hardware.  All tests wire in fake
> hardware for this purpose, but it represents a gap in actual, real-world
> testing.

At the telescope, `nickel_focus` runs under the observatory's own `kpython`
environment and is maintained at the system level by observatory staff —
observers use that system-level install directly, rather than a `pip`-managed
environment.

The maintenance workflow is as follows:

 - Maintainers edit their local install of the git repository.
 - Updates are merged to the remote GitHub repository.
 - New versions deployed via the observatory should ideally be tagged.
 - To deploy a new version, maintainers execute the `tools/deploy.sh` bash
   script on the telescope host computer.  This script:
    - Update the deployment-side git checkout to the latest `main`.  This is a
      hard reset, meaning any changes local to the host computer will be lost
      (i.e., the GitHub version is the authoritative version of the code).
    - Regenerate `nickel_focus/pkg/version.py` from the current git state using
     `tools/write_version.py`.
    - Sync the updated `nickel_focus/` tree into the CVS working copy.
    - Run `make install` from the CVS working copy to install the package into
      `$(LROOT)`.

