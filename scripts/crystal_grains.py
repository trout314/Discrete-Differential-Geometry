#!/usr/bin/env python3
"""MOVED: crystal_grains now lives in the package as
``discrete_differential_geometry.crystal_grains`` (2026-08 cleanup, Phase 1a).

The covering-map crystalline-grain detector -- the single authoritative
crystalline/defect identifier. This shim keeps legacy
``from crystal_grains import ...`` working for scripts that put scripts/ on
sys.path, and keeps the CLI runnable:

    python scripts/crystal_grains.py STATE.mfd
    python scripts/crystal_grains.py STATE.mfd --ref c15 r --min-size 20
    python scripts/crystal_grains.py 'data/melt_test/*.mfd' --json out/grains.json
"""
import importlib
import os
import sys

sys.path.insert(0, os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "python"))
_mod = importlib.import_module("discrete_differential_geometry.crystal_grains")
sys.modules[__name__] = _mod

if __name__ == "__main__":
    _mod.main()
