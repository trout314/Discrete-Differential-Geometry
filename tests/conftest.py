"""Shared fixtures: reference-crystal files for tests that exercise the
validator scripts at their real data paths.

Reference crystals under data/tcp_reference/ are gitignored, deterministic
outputs of scripts/tcp_reference.py -- if one is missing (fresh clone) it is
rebuilt here, exactly as the regenerate-references workflow does. The mgas
snapshots are NOT reproducible; tests needing them skip when absent.
"""
import os
import subprocess
import sys

import pytest

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _ensure_ref(struct, m, fname):
    path = os.path.join(REPO, "data", "tcp_reference", fname)
    if not os.path.exists(path):
        r = subprocess.run(
            [sys.executable, os.path.join(REPO, "scripts", "tcp_reference.py"),
             struct, "-m", str(m)],
            cwd=REPO, capture_output=True, text=True, timeout=600)
        assert r.returncode == 0 and os.path.exists(path), \
            f"could not rebuild {fname}:\n{r.stdout}\n{r.stderr}"
    return path


@pytest.fixture(scope="session")
def r_m2_ref():
    """The R m=2 reference crystal several validators hard-code."""
    return _ensure_ref("r", 2, "T3_R_m2_N7248.mfd")


@pytest.fixture(scope="session")
def a15_m3_ref():
    return _ensure_ref("a15", 3, "T3_A15_m3_N1242.mfd")


@pytest.fixture(scope="session")
def lam35r_snapshot():
    """The lam=0.35 R-melt snapshot some regressions are pinned to.
    Not regenerable -- skip when absent."""
    path = os.path.join(REPO, "data", "mgas", "lam35r_snap15000.mfd")
    if not os.path.exists(path):
        pytest.skip("data/mgas/lam35r_snap15000.mfd not present "
                    "(unreproducible snapshot)")
    return path
