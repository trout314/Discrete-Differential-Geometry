#!/usr/bin/env python3
"""MOVED: link_classes now lives in the package as
``discrete_differential_geometry.link_classes`` (2026-08 cleanup, Phase 1a).

This shim keeps legacy ``from link_classes import ...`` working for scripts
that put tools/ on sys.path, and keeps the CLI runnable:

    python tools/link_classes.py --enumerate 6
    python tools/link_classes.py data/tcp_reference/T3_R_m2_N1360.mfd
    python tools/link_classes.py --verify <seed.mfd>
"""
import importlib
import os
import sys

sys.path.insert(0, os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "python"))
_mod = importlib.import_module("discrete_differential_geometry.link_classes")
sys.modules[__name__] = _mod

if __name__ == "__main__":
    _mod.main()
