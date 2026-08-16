#!/usr/bin/env python3
"""MOVED: tip_retract_search now lives in the research package as ``ddg_lab.tip_retract_search``
(2026-08 cleanup, Phase 1b -- see notes/CLEANUP_PLAN.md).

Import shim + CLI dispatch: legacy ``from tip_retract_search import ...`` (with this
directory on sys.path) and ``python scripts/defect_dynamics/tip_retract_search.py`` both keep working; new code
should import ``ddg_lab.tip_retract_search``.
"""
import importlib
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "python"))
if __name__ == "__main__":
    import runpy
    runpy.run_module("ddg_lab.tip_retract_search", run_name="__main__", alter_sys=True)
else:
    sys.modules[__name__] = importlib.import_module("ddg_lab.tip_retract_search")
