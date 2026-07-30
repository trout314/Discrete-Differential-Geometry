"""Post-refactor regression: D worm_enum vs Python oracle on lam35r_snap15000.

Validated invariant from the port: at every anchor, D's candidate set is a
SUPERSET of the oracle's (D has no radius-2 patch cutoff), and every oracle
dS value appears among D's trial dS values.  Checks counts + dS multiset
inclusion at all 121 anchors.
"""
import os, sys
from collections import Counter

_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, os.path.join(_ROOT, "python"))
sys.path.insert(0, os.path.join(_ROOT, "scripts/defect_dynamics"))
import numpy as np
import discrete_differential_geometry as ddg
import worm_deg4_slide as W

SNAP = os.path.join(_ROOT, "data/mgas/lam35r_snap15000.mfd")
LAM, ESTAR, ZLEG, CIMP = 0.35, 5.105025, 0.6, 1.0

m = ddg.Manifold.load(SNAP, 3)
p = ddg.SamplerParams(
    num_facets_target=m.num_facets, hinge_degree_target=ESTAR,
    num_facets_coef=0.1, num_hinges_coef=0.0,
    hinge_degree_variance_coef=0.0, codim3_degree_variance_coef=0.0,
    hinge_degree_target_coef=LAM * ESTAR / 6.0)
s = ddg.ManifoldSampler(m, p)
s.set_n6_potential(ZLEG * LAM, CIMP * LAM, tilt=[0.0] * 5)
L = W.Live(s.manifold)
anchors = sorted(L.deg4())
print(f"{len(anchors)} anchors")

bad = 0
no, nd = 0, 0
for a in anchors:
    oc = W.candidates(L, a)
    land, dsd = s.worm_enum(a[0], a[1])
    no += len(oc)
    nd += len(land)
    if len(land) < len(oc):
        print(f"  FAIL {a}: D {len(land)} < oracle {len(oc)}")
        bad += 1
        continue
    co = Counter(round(c[2], 8) for c in oc)
    cd = Counter(round(float(x), 8) for x in dsd)
    miss = co - cd
    if miss:
        print(f"  FAIL {a}: oracle dS not in D: {dict(miss)} "
              f"(oracle {len(oc)}, D {len(land)})")
        bad += 1
print(f"oracle candidates: {no}, D candidates: {nd}")
print("PASS: D superset + dS multiset inclusion at every anchor"
      if bad == 0 else f"FAIL at {bad} anchors")
