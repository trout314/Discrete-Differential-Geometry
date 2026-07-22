#!/usr/bin/env python3
"""Hunt for the constrained-liquid window: a coupling scale where FK defects are
PRESENT and MOBILE (knot gas equilibrates by its own dynamics) while the state
stays strongly constrained (legal-edge fraction high).

Sweep the overall coupling LAM multiplying the run5h objective (edge pin at
native+bump with hinge_degree_target_coef, plus the n6 FK push zleg/cimp);
volume pin fixed. At each LAM: burn, then sample the illegal-vertex set every
DT sweeps; report density (nill, ncomp), constraint (legaledge/legalvert), and
MOBILITY = Jaccard overlap J(lag) of the illegal set (frozen -> plateau at
~0.3-0.6; mobile -> decays to the random baseline ~ nill/V).

Usage: mobility_sweep.py CELL out.jsonl [lam1 lam2 ...]
"""
import json
import os
import sys
import time
from collections import Counter
from itertools import combinations

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
for p in ("../../python", "../../scripts"):
    sys.path.insert(0, os.path.join(_HERE, p))
import discrete_differential_geometry as ddg
from fk_skeleton import edges_from_facets
from dopant_pairs import vertex_classes

CELL = sys.argv[1]
OUT = sys.argv[2]
LAMS = [float(x) for x in sys.argv[3:]] or [1.0, 0.7, 0.5, 0.35, 0.25, 0.18, 0.12, 0.08]
BUMP, ZLEG, CIMP = 8e-4, 0.6, 1.0
BURN, SPAN, DT = 1500, 3000, 50
LAGS = [50, 150, 300, 600, 1200]

ref = ddg.Manifold.load(CELL, 3)
native = float(edges_from_facets(ref.facets())[1].mean())
out = open(OUT, "w")

for lam in LAMS:
    t0 = time.time()
    ddg.set_random_seed(int(1000 * lam) + 7)
    m = ddg.Manifold.load(CELL, 3)
    et = native + BUMP
    params = ddg.SamplerParams(
        num_facets_target=ref.num_facets, num_facets_coef=0.1,
        hinge_degree_target=et, num_hinges_coef=0.0,
        hinge_degree_variance_coef=0.0, codim3_degree_variance_coef=0.0,
        hinge_degree_target_coef=lam * et / 6.0)
    s = ddg.ManifoldSampler(m, params)
    s.set_n6_potential(ZLEG * lam, CIMP * lam, tilt=[0.0] * 5)
    v = s.manifold
    s.run(sweeps=BURN)

    sets, nill, ncomp, le, lv = [], [], [], [], []
    for k in range(SPAN // DT):
        s.run(sweeps=DT)
        fac = np.asarray(v.facets())
        eu, edeg, V = edges_from_facets(fac)
        n6, imp, adj = vertex_classes(fac)
        ill = frozenset(int(x) for x in np.nonzero(imp > 0)[0])
        sets.append(ill)
        nill.append(len(ill))
        le.append(float(np.mean((edeg == 5) | (edeg == 6))))
        lv.append(float(np.mean(imp == 0)))
        seen = set()
        nc = 0
        for x in ill:
            if x in seen:
                continue
            nc += 1
            st = [x]; seen.add(x)
            while st:
                u = st.pop()
                for w in adj[u]:
                    if imp[w] > 0 and w not in seen:
                        seen.add(w); st.append(w)
        ncomp.append(nc)

    Vn = V
    jac = {}
    for lag in LAGS:
        step = lag // DT
        vals = []
        for i in range(len(sets) - step):
            a, b = sets[i], sets[i + step]
            if a or b:
                vals.append(len(a & b) / len(a | b))
        jac[lag] = float(np.mean(vals)) if vals else None
    base = float(np.mean(nill)) / (2 * Vn) if np.mean(nill) else 0.0
    row = dict(lam=lam, V=Vn, nill=float(np.mean(nill)), ncomp=float(np.mean(ncomp)),
               legaledge=float(np.mean(le)), legalvert=float(np.mean(lv)),
               jaccard={str(k): v for k, v in jac.items()}, jac_random=base,
               secs=round(time.time() - t0, 1))
    out.write(json.dumps(row) + "\n"); out.flush()
    js = " ".join(f"J{k}={v:.2f}" if v is not None else f"J{k}=-"
                  for k, v in jac.items())
    print(f"lam={lam:5.2f}: nill={row['nill']:7.1f} ncomp={row['ncomp']:5.1f} "
          f"legaledge={row['legaledge']:.3f} legalvert={row['legalvert']:.3f}  "
          f"{js}  (rand {base:.3f})  [{row['secs']}s]", flush=True)
out.close()
print("done")
