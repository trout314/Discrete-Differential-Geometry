---
name: run-scripts-with-caffeinate
description: Always launch scripts/compute on this machine wrapped in caffeinate (prevents idle-sleep killing runs)
metadata: 
  node_type: memory
  type: feedback
  originSessionId: 1ed7f5d0-df5c-4552-b65e-5060180eafe2
  modified: 2026-07-21T14:09:24.740Z
---

Always run scripts on this machine wrapped in `caffeinate` (Aaron's standing
instruction, 2026-07-21).

**Why:** this machine idle-sleeps and the sleep KILLS running/background
processes. A 5-hour production MCMC run (the [[five-illegal-knot]] study) was
killed ~12 min in by overnight idle-sleep, wasting the launch. caffeinate holds
a power assertion that prevents idle system sleep for the wrapped command's
lifetime.

**How to apply:** prefix launches with `caffeinate -i`, e.g.
`caffeinate -i zsh launcher.sh` or `caffeinate -i python3 script.py`. Do this for
ANY backgrounded or multi-minute compute (and it's harmless on quick ones, so
default to it). caffeinate runs the command as a child and releases the assertion
when it exits. Verify it's active with
`pmset -g assertions | grep PreventUserIdleSystemSleep` (value 1 = held).
CAVEAT: `-i` stops IDLE sleep only — closing the laptop lid can still force sleep
(tell Aaron to keep the lid open / stay on AC for long runs). Add `-s` to also
prevent system sleep on AC.
