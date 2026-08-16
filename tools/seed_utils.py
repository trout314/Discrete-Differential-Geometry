"""MOVED: seed_utils now lives in the package as
``discrete_differential_geometry.seed_utils`` (2026-08 cleanup, Phase 1a).

This shim keeps legacy ``from seed_utils import ...`` working for scripts
that put tools/ on sys.path; new code should import the package module.
The sys.modules swap makes ``import seed_utils`` yield the real package
module, so every attribute (including private ones) stays reachable.
"""
import importlib
import os
import sys

sys.path.insert(0, os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "python"))
sys.modules[__name__] = importlib.import_module(
    "discrete_differential_geometry.seed_utils")
