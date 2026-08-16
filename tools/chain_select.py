"""MOVED: chain_select now lives in the package as
``discrete_differential_geometry.chain_select`` (2026-08 cleanup, Phase 1a).

This shim keeps legacy ``from chain_select import ...`` working for scripts
that put tools/ on sys.path; new code should import the package module.
"""
import importlib
import os
import sys

sys.path.insert(0, os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "python"))
sys.modules[__name__] = importlib.import_module(
    "discrete_differential_geometry.chain_select")
