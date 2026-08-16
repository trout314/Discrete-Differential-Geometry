"""The latent validator scripts, run as pytest (Phase 3 of the cleanup).

Each test runs a validator script in a subprocess from the repo root with
reduced-but-meaningful parameters. Scripts that assert / exit nonzero on
failure are checked by return code; scripts that only PRINT their verdict
are checked against their exact verdict strings.

Kept as scripts + wrapped (rather than transcribed) so the campaign-scale
invocations documented in each script keep working unchanged.
"""
import os
import subprocess
import sys

import pytest

from conftest import REPO


def run_script(relpath, *args, timeout=600):
    r = subprocess.run([sys.executable, os.path.join(REPO, relpath), *args],
                       cwd=REPO, capture_output=True, text=True,
                       timeout=timeout)
    return r


def check(r, *expect_stdout):
    assert r.returncode == 0, f"rc={r.returncode}\n{r.stdout[-2000:]}\n{r.stderr[-2000:]}"
    for s in expect_stdout:
        assert s in r.stdout, f"missing {s!r} in:\n{r.stdout[-2000:]}"


# --------------------------------------------------------------------------
# fast: exact, self-contained batteries
# --------------------------------------------------------------------------

def test_fpkmc_m0_formulas():
    """SymPy verification of the FPKMC chain formulas (exact)."""
    check(run_script("scripts/defect_dynamics/fpkmc_m0_derivations.py"),
          "ALL PASS")


def test_transport_battery(r_m2_ref):
    """Quaternion/holonomy toolkit: exact equality battery T1-T5."""
    check(run_script("scripts/defect_dynamics/transport_tests.py"),
          "ALL TRANSPORT TESTS PASS (exact)")


def test_transport_injectivity(r_m2_ref):
    """Wilson rotation determines the path (mod backtracking) at depth 5."""
    check(run_script("scripts/defect_dynamics/transport_injectivity_check.py",
                     "--maxlen", "5"),
          "INJECTIVE at this depth")


def test_samerung_slide_exactly_free(r_m2_ref):
    """Non-local slide: dS = c*dQ exact; same-rung hops cost exactly zero."""
    check(run_script("scripts/defect_dynamics/samerung_validate.py"),
          "max |dS - c*dQ| over all targets: 0.000e+00",
          "max |dS| on same-rung hops: 0.000e+00")


# --------------------------------------------------------------------------
# slow: sampler-scale runs / pinned-snapshot regressions
# --------------------------------------------------------------------------

@pytest.mark.slow
def test_worm_cycles_path_independence(a15_m3_ref):
    """Apply/undo DFS over worm cycles on A15 m=3 (walks=1)."""
    r = run_script("scripts/defect_dynamics/worm_cycles.py",
                   a15_m3_ref, "--walks", "1")
    check(r, "creation-orbit representatives")


@pytest.mark.slow
def test_hb_kernel_detailed_balance(r_m2_ref):
    """FPKMC V3b: numerical DB certification of the HB kernel (10 pairs)."""
    check(run_script("scripts/defect_dynamics/fpkmc_v3b_db.py",
                     "--pairs", "10"),
          "PASS")


@pytest.mark.slow
def test_contract_split_pair_test():
    """Labeled pair test for the D-side contract/split channel (1 pair)."""
    check(run_script("scripts/validate_contract_split.py",
                     "--pair-test", "--trials", "250000", "--pairs", "1"),
          "OVERALL: PASS")


@pytest.mark.slow
def test_worm_enum_vs_oracle(lam35r_snapshot):
    """D worm_enum vs the Python oracle: superset + dS multiset inclusion.

    Also the regression guard for the ddg_sampler_worm_at binding signature
    (a stale 10th argtype broke this call until 2026-08-16)."""
    check(run_script("scripts/defect_dynamics/worm_crossval2.py"),
          "PASS: D superset + dS multiset inclusion at every anchor")


@pytest.mark.slow
def test_fusion_two_move(lam35r_snapshot):
    """The accidental 2-move fusion on the D core still validates."""
    check(run_script("scripts/defect_dynamics/fusion_verify.py"))
