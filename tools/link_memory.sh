#!/bin/sh
# Point Claude Code's per-project memory store at the copy tracked in this repo.
#
# Claude Code keeps project memory in ~/.claude/projects/<slug>/memory, where
# <slug> is the project's absolute path with '/' replaced by '-'. That is
# outside the repo and under no version control, so the research narrative --
# campaign history, retractions, measured constants -- had no history and no
# backup. This replaces that directory with a symlink to notes/memory/ so the
# repo is the single source of truth; nothing syncs, there is only one copy.
#
# Safe to re-run. Refuses to clobber an existing non-empty real directory.
set -e

ROOT=$(cd "$(dirname "$0")/.." && pwd)
SRC="$ROOT/notes/memory"
SLUG=$(printf '%s' "$ROOT" | tr '/' '-')
DEST="$HOME/.claude/projects/$SLUG/memory"

[ -d "$SRC" ] || { echo "no $SRC in this checkout"; exit 1; }

if [ -L "$DEST" ]; then
    cur=$(readlink "$DEST")
    if [ "$cur" = "$SRC" ]; then
        echo "already linked: $DEST -> $SRC"
        exit 0
    fi
    echo "re-pointing $DEST (was -> $cur)"
    rm "$DEST"
elif [ -e "$DEST" ]; then
    if [ -n "$(ls -A "$DEST" 2>/dev/null)" ]; then
        echo "REFUSING: $DEST exists and is not empty."
        echo "Merge it into $SRC by hand, then remove it and re-run:"
        ls -A "$DEST" | sed 's/^/  /'
        exit 1
    fi
    rmdir "$DEST"
fi

mkdir -p "$(dirname "$DEST")"
ln -s "$SRC" "$DEST"
echo "linked: $DEST -> $SRC"
