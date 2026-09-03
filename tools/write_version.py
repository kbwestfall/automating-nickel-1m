#!/usr/bin/env python3
"""
Regenerate ``nickel_focus/pkg/version.py`` from the current git tag/commit,
without requiring a ``pip install`` of this package.

``setuptools_scm`` (already a ``[build-system] requires`` dependency in
``pyproject.toml``) normally writes this file as a side effect of a pip
build. This script calls its API directly instead, so the file can be
regenerated on demand -- e.g. immediately before mirroring a tagged commit
to CVS, where no pip build ever happens.

Requires ``setuptools_scm`` to be installed in whatever Python runs this
script (``pip install setuptools_scm``); it is not a runtime dependency of
``nickel_focus`` itself.
"""
import sys
import tomllib
import warnings
from pathlib import Path

import setuptools_scm

REPO_ROOT = Path(__file__).resolve().parent.parent


def _setuptools_scm_config():
    """
    Read the ``[tool.setuptools_scm]`` table from the repo's
    ``pyproject.toml``.

    Returns
    -------
    dict
        The ``[tool.setuptools_scm]`` table, e.g. containing the configured
        ``version_file`` path (relative to the repo root).
    """
    pyproject = REPO_ROOT / "pyproject.toml"
    with open(pyproject, "rb") as f:
        config = tomllib.load(f)
    return config.get("tool", {}).get("setuptools_scm", {})


def main():
    """
    Regenerate the version file configured in ``pyproject.toml``'s
    ``[tool.setuptools_scm]`` table, and print the computed version. Warns
    (but does not fail) if the computed version doesn't look like a clean
    tagged release, since CVS syncs are expected to only ever happen from
    tagged commits.
    """
    scm_config = _setuptools_scm_config()
    version_file = scm_config.get("version_file")
    if version_file is None:
        raise RuntimeError(
            "No [tool.setuptools_scm] version_file configured in pyproject.toml; "
            "nothing to write."
        )

    # `root` is already absolute, so `version_file` (relative to `root`, per
    # pyproject.toml's own convention) resolves correctly on its own;
    # `relative_to` is only passed because setuptools_scm asserts it's not
    # None internally, not because its value is actually used here -- which
    # is exactly what triggers the (harmless, in this case) warning below.
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore", message="absolute root path .* overrides relative_to .*"
        )
        version = setuptools_scm.get_version(
            root=REPO_ROOT,
            relative_to=REPO_ROOT / "pyproject.toml",
            version_file=version_file,
        )

    if "dev" in version or "dirty" in version:
        print(
            f"Warning: computed version '{version}' does not look like a clean "
            "tagged release.",
            file=sys.stderr,
        )

    print(f"Wrote {version_file}: {version}")


if __name__ == "__main__":
    main()
