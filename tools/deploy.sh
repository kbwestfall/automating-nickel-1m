#!/usr/bin/env bash
#
# Deploy the latest nickel_focus release to the telescope's kroot/CVS-
# managed system installation.
#
# Workflow:
#   1. Update the deployment-side git checkout to the latest `main`.
#   2. Regenerate nickel_focus/pkg/version.py from the current git state.
#   3. Sync the updated nickel_focus/ tree into the CVS working copy.
#   4. Run `make install` from the CVS working copy to install the
#      package into $(LROOT).
#
# NOTE: step 3's commands are nominal placeholders -- the exact CVS
# conventions (working-copy location, add/remove/commit mechanics) are
# still TBD pending confirmation from Brad Holden. See
# claude/DEPLOYMENT.md and claude/DEPLOYMENT_PLAN.md for context. This
# script has not been tested end-to-end; treat it as a draft to finish
# with Brad/Will, not as a ready-to-run tool.

set -euxo pipefail

# --- Configuration (TBD) --------------------------------------------------

# Path to this repo's deployment-side git clone. Expected to always sit on
# the `main` branch; this script forcibly resets it to match origin/main,
# discarding any local changes there.
GIT_CHECKOUT="${GIT_CHECKOUT:-/TBD/path/to/nickel_focus/git/checkout}"

# Path to the CVS working copy that `make install` runs from -- i.e. the
# CVS-side mirror of this repo's nickel_focus/ tree (plus the top-level
# Makefile). Assumed to already be a valid CVS checkout (its CVS/Root
# metadata already points at the right CVSROOT). TBD once Brad's CVS
# conventions are confirmed.
CVS_CHECKOUT="${CVS_CHECKOUT:-/TBD/path/to/nickel_focus/cvs/checkout}"

# ---------------------------------------------------------------------------

echo "==> Step 1: updating $GIT_CHECKOUT to the latest main"
cd "$GIT_CHECKOUT"

if [ -n "$(git status --porcelain)" ]; then
    echo "WARNING: $GIT_CHECKOUT has local changes that are about to be discarded:" >&2
    git status --short >&2
fi

git fetch origin main
git checkout main
git reset --hard origin/main

if ! git describe --tags --exact-match HEAD >/dev/null 2>&1; then
    echo "WARNING: HEAD is not exactly a tagged commit; CVS syncs are expected" >&2
    echo "         to only happen from tagged releases." >&2
fi

echo "==> Step 2: regenerating nickel_focus/pkg/version.py"
# Requires setuptools_scm to be installed in whatever python3 is on PATH
# here (e.g. via `pip install .[dev]`) -- it is not a runtime dependency
# of nickel_focus itself. See tools/write_version.py.
python3 "$GIT_CHECKOUT/tools/write_version.py"

echo "==> Step 3: syncing nickel_focus/ into the CVS working copy"

rsync -a --delete "$GIT_CHECKOUT/nickel_focus/" "$CVS_CHECKOUT/nickel_focus/"
rsync -a "$GIT_CHECKOUT/Makefile" "$CVS_CHECKOUT/Makefile"

cd "$CVS_CHECKOUT"

# Restore the executable bit on .sin files -- CVS does not track POSIX
# permissions the way git does, so rsync alone won't preserve this.
find nickel_focus -name '*.sin' -exec chmod +x {} \;

# `cvs -qn update` (dry run) reports '?' for files on disk that CVS
# doesn't know about yet, and '!' for files CVS tracks that are now
# missing from disk.
cvs -qn update | awk '/^\?/ {print $2}' | xargs -r cvs add
cvs -qn update | awk '/^!/ {print $2}' | xargs -r cvs remove

cvs commit -m "Deploy nickel_focus $(git -C "$GIT_CHECKOUT" describe --tags)"

echo "==> Step 4: running make install"
cd "$CVS_CHECKOUT"
make install

echo "==> Deployment complete."
