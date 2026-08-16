#!/usr/bin/env python3
"""MOVED: tcp_melt now lives in the research package as ``ddg_lab.tcp_melt``
(2026-08 cleanup, Phase 1b -- see notes/CLEANUP_PLAN.md).

Import shim + CLI dispatch: legacy ``from tcp_melt import ...`` (with this
directory on sys.path) and ``python scripts/tcp_melt.py`` both keep working; new code
should import ``ddg_lab.tcp_melt``.
"""
import importlib
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "python"))
if __name__ == "__main__":
    import runpy
    runpy.run_module("ddg_lab.tcp_melt", run_name="__main__", alter_sys=True)
else:
    sys.modules[__name__] = importlib.import_module("ddg_lab.tcp_melt")
