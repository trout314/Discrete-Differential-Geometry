#!/usr/bin/env python3
"""MOVED: sixhundred_cell now lives in the research package as ``ddg_lab.sixhundred_cell``
(2026-08 cleanup, Phase 1b -- see notes/CLEANUP_PLAN.md).

Import shim + CLI dispatch: legacy ``from sixhundred_cell import ...`` (with this
directory on sys.path) and ``python scripts/sixhundred_cell.py`` both keep working; new code
should import ``ddg_lab.sixhundred_cell``.
"""
import importlib
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "python"))
if __name__ == "__main__":
    import runpy
    runpy.run_module("ddg_lab.sixhundred_cell", run_name="__main__", alter_sys=True)
else:
    sys.modules[__name__] = importlib.import_module("ddg_lab.sixhundred_cell")
