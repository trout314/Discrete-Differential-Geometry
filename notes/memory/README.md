# Project memory

Claude Code's per-project memory for this repository, tracked here rather than
left in `~/.claude/`.

`MEMORY.md` is the index — one line per memory, loaded at the start of every
session. Each other file is one fact or one campaign's state, with frontmatter
giving its name, a one-line description used for recall, and a type
(`user` / `feedback` / `project` / `reference`).

## Why these live in the repo

Claude Code writes memory to `~/.claude/projects/<slug>/memory`, outside the
working tree and under no version control. That directory held the entire
research narrative — campaign history, gotchas, retractions, measured constants
whose runs took days — with no history, no diffs, and no backup. Moving it here
puts it under the same version control as the code it describes, so a claim and
the commit that established (or retracted) it sit side by side.

## Hooking it up on a new machine or a fresh clone

The `~/.claude` path is **not** created by cloning. Run:

    tools/link_memory.sh

which replaces `~/.claude/projects/<slug>/memory` with a symlink to this
directory. There is exactly one copy of every file and nothing to sync. The
script is safe to re-run and refuses to clobber an existing non-empty store
(if you already have memories there, merge them in by hand first).

Note `<slug>` is the project's absolute path with `/` replaced by `-`, so it
differs per machine and per checkout location; the script computes it.

## Reading these critically

They are a working record, not a set of established results. Specifically:

- **They reflect what was true when written.** Several have been retracted or
  amended in place — see `flicker-catalysis.md` (retracted outright) and the
  partial retraction inside `worm-sampler-program.md`. If a memory names a
  file, function or flag, verify it still exists before acting on it.
- **Claims of exhaustiveness deserve suspicion.** Three separate enumerations
  in this project turned out to have coverage defects
  (`smatrix-sweep-running.md` records all three for its sweep alone). Prefer
  the code's own docstring, which is versioned alongside the implementation.
- **Load-bearing corrections are duplicated into code docstrings** precisely so
  they survive independently of this directory.
