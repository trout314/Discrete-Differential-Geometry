"""Smoke tests for the ddg_lab research package (Phase 1b of the cleanup).

Every module must import cleanly (no sys.path hacks, no argv/env reads at
module level), and the shims at the old script locations must hand back the
very same module objects. Real invariant tests for this machinery arrive
with the Phase-3 validator conversion; this is the bitrot guard.
"""
import importlib
import os
import sys

import pytest

MODULES = ["defect_state", "event_replay", "crystal_flicker", "worm_moves",
           "worm_helix", "worm_slide", "dressed_generators", "worm_deg4_slide",
           "f0_worm", "link_planner", "fk_moves", "tip_retract_search",
           "dopant_pairs", "tcp_melt", "graft_signature", "steiner_geodesic",
           "heat_geodesic", "sixhundred_cell", "grid_sweep"]

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


@pytest.mark.parametrize("name", MODULES)
def test_imports(name):
    importlib.import_module(f"ddg_lab.{name}")


def test_shims_hand_back_the_package_module():
    for p in ("scripts", "tools", "scripts/defect_dynamics"):
        sys.path.insert(0, os.path.join(_ROOT, p))
    try:
        import worm_slide, grid_sweep, tcp_melt          # via the shims
        assert worm_slide is importlib.import_module("ddg_lab.worm_slide")
        assert grid_sweep is importlib.import_module("ddg_lab.grid_sweep")
        assert tcp_melt is importlib.import_module("ddg_lab.tcp_melt")
    finally:
        for _ in range(3):
            sys.path.pop(0)


def test_worm_moves_on_a_crystal():
    """Cheap invariant spot check: move-site census on pristine A15."""
    import numpy as np
    from discrete_differential_geometry.tcp_reference import build_t3_triangulation
    from ddg_lab import worm_moves as wm
    F = np.asarray(build_t3_triangulation("a15", 2)[0])
    faces, edeg, vedges = wm.build_tables(F)
    # pristine TCP: every edge degree in {5, 6}
    assert set(edeg.values()) <= {5, 6}
    # every 2->3 site is a face whose opposite apexes are non-adjacent
    sites = list(wm.two_three_sites(F, faces, edeg))
    assert len(sites) > 0
