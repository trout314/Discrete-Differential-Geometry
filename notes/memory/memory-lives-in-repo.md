---
name: memory-lives-in-repo
description: Project memory is TRACKED at notes/memory/ in the repo (symlinked from ~/.claude); run tools/link_memory.sh on a new machine; the repo is PUBLIC
metadata: 
  node_type: memory
  type: project
  originSessionId: c395a1df-09fb-4d44-bb80-ddc781dd7a5b
  modified: 2026-08-03T00:50:43.709Z
---

Since 2026-08-02, this project's Claude Code memory is version-controlled with
the code, at **`notes/memory/`** in the repo.
`~/.claude/projects/<slug>/memory` is a **symlink** to it, so there is exactly
one copy and nothing to sync — writing a memory the normal way lands it in the
working tree as an untracked/modified file, ready to commit.

**On a new machine or a fresh clone the link does not exist** (cloning does not
create `~/.claude`). Run:

    tools/link_memory.sh

Idempotent; refuses to clobber an existing non-empty store rather than
overwrite memories written elsewhere.

**Commit memory changes.** They are ordinary tracked files now, so an
uncommitted memory is just an uncommitted file — it shows in `git status` and
is lost with the working tree. Offer to commit them alongside the code change
they describe, so a claim and the commit that established or retracted it sit
together.

**The repo is PUBLIC** (`trout314/Discrete-Differential-Geometry`). The corpus
was scanned before publishing — no credentials, home paths, emails or external
URLs — but everything written here from now on is public the moment it is
pushed. Keep it to research content.

**GOTCHA that nearly ate this.** `.gitignore` carried a bare `memory/` rule
commented "Claude Code local memory store", matching any directory of that
name. Committing without checking would have added the README and tooling but
silently NOT the 59 memory files, leaving `git status` clean afterwards. Caught
with `git check-ignore`, not by trusting `git add`. If memory files ever stop
appearing in `git status`, suspect an ignore rule first.

See `notes/memory/README.md` for how to read the corpus critically (it is a
working record, not established results — several entries are retracted or
amended in place, and claims of exhaustiveness in this project have a bad
track record: see [[smatrix-sweep-running]] and the partial retraction in
[[worm-sampler-program]]).
