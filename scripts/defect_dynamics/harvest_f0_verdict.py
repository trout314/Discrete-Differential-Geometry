"""Thermodynamic verdict on an f0 harvest: continued sampling of the
before/after states under the SAME action; compare late-window mean
objective / pin gap / n_ill. (Provisional single-chain late-window
means -- not a certified stationarity claim.)

Usage:
  python scripts/defect_dynamics/harvest_f0_verdict.py A.mfd B.mfd \
      [ETARGET [CIMP [F3T]]]
Env: SWEEPS (default 400), STRIDE (default 10).

RESULT 2026-07-30 (quench_down5q_wOFF vs its 14-removal harvest, 3000
sweeps): quenched f0=1536 late<S>=262+-12 gap +10.2 n_ill 63; harvested
f0=1522 late<S>=357+-20 gap -0.01 n_ill 189 -- the gap-closed sector is
defect-financed (~9 illegal edges per vacancy), stationary, no anneal.
"""
import os, sys
import numpy as np

_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, os.path.join(_ROOT, "python"))
import discrete_differential_geometry as ddg

ETARGET = float(sys.argv[3]) if len(sys.argv) > 3 else 5.1067907
CIMP = float(sys.argv[4]) if len(sys.argv) > 4 else 0.7
F3T = int(sys.argv[5]) if len(sys.argv) > 5 else 8704
SWEEPS = int(os.environ.get("SWEEPS", "400"))
STRIDE = int(os.environ.get("STRIDE", "10"))

for fn in (sys.argv[1], sys.argv[2]):
    m = ddg.Manifold.load(fn, 3)
    p = ddg.SamplerParams(
        num_facets_target=F3T, num_facets_coef=0.1,
        hinge_degree_target=ETARGET, num_hinges_coef=1.0,
        hinge_degree_variance_coef=0.0, codim3_degree_variance_coef=0.0,
        hinge_degree_target_coef=0.0, codim3_degree_target_coef=0.0)
    s = ddg.ManifoldSampler(m, p)
    s.set_n6_potential(0.0, CIMP, tilt=[0.0] * 5)
    ddg.set_random_seed(777)
    S, G, NI = [], [], []
    for i in range(SWEEPS // STRIDE):
        s.run(sweeps=STRIDE)
        fv = [int(x) for x in s.manifold.f_vector]
        pairs, degs = s.manifold.illegal_edges()
        S.append(s.current_objective)
        G.append(fv[1] - 6.0 * fv[3] / ETARGET)
        NI.append(len(degs))
    S, G, NI = map(np.array, (S, G, NI))
    h = len(S) // 2
    print(f"{os.path.basename(fn):40s}: f0={fv[0]}  "
          f"late<S>={S[h:].mean():8.2f} +-{S[h:].std():5.2f}  "
          f"late<gap>={G[h:].mean():+6.2f}  "
          f"late<n_ill>={NI[h:].mean():6.1f}  "
          f"(first/last S {S[0]:.1f}/{S[-1]:.1f})", flush=True)
print("verdict run done", flush=True)
