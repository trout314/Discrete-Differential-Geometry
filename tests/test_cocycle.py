"""The T^3 integer-cocycle validation ladder as pytest (ddg.cocycle).

Pytest port of the scripts/cocycle_check.py ladder, on a small A15 m=2
crystal with light churn:

  0. a tampered cocycle must be REJECTED at enable
  1. the pristine crystal cocycle is accepted
  2. zero-churn harmonic gauge recovers the crystal's own coordinates
  3. closedness is exact after churn
  4. the fundamental-cycle winding lattice is still (scale*m) Z^3 after churn
  5. harmonic gauge + consumers run on the churned state
"""
import numpy as np
import pytest

import discrete_differential_geometry as ddg
from discrete_differential_geometry import cocycle as coc
from discrete_differential_geometry.tcp_reference import (
    build_t3_triangulation, reference_frac_positions)

SCALE = 10 ** 6
M_SUP = 2
SWEEPS = 3


def wrapped_rms(a, b, M):
    """RMS deviation between torus position sets after removing the global
    translation (wrap-aware)."""
    d = a - b
    d = d - d[0]
    d -= M * np.round(d / M)
    d -= d.mean(axis=0)
    d -= M * np.round(d / M)
    return float(np.sqrt((d ** 2).sum(axis=1).mean()))


@pytest.fixture(scope="module")
def ladder():
    fac, n_verts = build_t3_triangulation("a15", M_SUP)
    frac = reference_frac_positions("a15", M_SUP)
    mfd = ddg.Manifold(3, np.asarray(fac).tolist())
    qbar = mfd.mean_degree(1)
    params = ddg.SamplerParams(
        num_facets_target=mfd.num_facets, num_facets_coef=0.1,
        hinge_degree_target=qbar, num_hinges_coef=2.0,
        hinge_degree_target_coef=0.55 * qbar / 6.0)
    ddg.set_random_seed(3)
    s = ddg.ManifoldSampler(mfd, params)

    edges = np.asarray(s.manifold.simplices(1))
    omega = coc.build_from_positions(edges, frac, M_SUP, scale=SCALE)
    out = dict(n_verts=n_verts, frac=frac, edges=edges, omega=omega,
               edge_len=float(np.linalg.norm(omega, axis=1).mean()))

    # rung 0: tampered cocycle rejected at enable
    bad = omega.copy()
    bad[0, 0] += 1
    try:
        s.enable_cocycle(edges, bad)
        out["tamper_rejected"] = False
    except RuntimeError:
        out["tamper_rejected"] = True

    # rung 1: pristine cocycle accepted
    s.enable_cocycle(edges, omega)

    # rung 3 setup: churn, then audit
    s.run(sweeps=SWEEPS)
    s.check_cocycle()                       # raises if closedness broken
    out["accepted"] = s.get_stats().total_accepted
    out["e2"], out["w2"] = s.read_cocycle()
    return out


def test_rung0_tampered_rejected(ladder):
    assert ladder["tamper_rejected"]


def test_rung2_harmonic_round_trip(ladder):
    n, frac = ladder["n_verts"], ladder["frac"]
    edges, omega = ladder["edges"], ladder["omega"]
    _, phi, _ = coc.harmonic_gauge(edges, omega, n)
    pos, _, _ = coc.tree_positions(edges, omega, n)
    rms = wrapped_rms(pos - phi, SCALE * frac, SCALE * M_SUP)
    assert rms < 0.3 * ladder["edge_len"]


def test_rung3_churn_happened(ladder):
    # check_cocycle already passed inside the fixture; churn must be real
    assert ladder["accepted"] > 0


def test_rung4_winding_lattice(ladder):
    e2, w2 = ladder["e2"], ladder["w2"]
    n2 = int(e2.max()) + 1
    _, cyc2, _ = coc.tree_positions(e2, w2, n2)
    basis = coc.lattice_basis(cyc2)
    assert len(basis) == 3
    det = abs(int(basis[0, 0]) * int(basis[1, 1]) * int(basis[2, 2]))
    assert det == (SCALE * M_SUP) ** 3


def test_rung5_consumers_run(ladder):
    e2, w2 = ladder["e2"], ladder["w2"]
    n2 = int(e2.max()) + 1
    oh2, phi2, gram2 = coc.harmonic_gauge(e2, w2, n2)
    ev = np.linalg.eigvalsh(gram2 / np.trace(gram2) * 3)
    assert np.all(ev > 0)                   # positive-definite cell shape
    _, _, S_nem = coc.nematic_q(oh2)
    assert -0.5 - 1e-9 <= S_nem <= 1.0 + 1e-9
