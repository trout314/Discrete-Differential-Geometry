#!/usr/bin/env python3
"""Pin-release test: continue a run5h final snapshot with BUMP=0 (edge-degree
target back to native) and everything else identical (zleg=0.6, cimp=1.0,
nhinge=6.0, VDV/HDV off). Watch whether the knot multimers anneal away
(thermodynamically stabilized by the bump) or persist (kinetically trapped).

args: start.mfd out_prefix sweeps seed
Logs n_illegal + complex members every 150 sweeps (jsonl, flushed); saves a
resume .mfd every 2500 sweeps.
"""
import json
import os
import sys
import time

import numpy as np

_R = "/Users/atrout/Desktop/Discrete-Differential-Geometry"
for p in ("python", "scripts", "tools"):
    sys.path.insert(0, os.path.join(_R, p))
import discrete_differential_geometry as ddg
from dopant_pairs import vertex_classes
from fk_skeleton import edges_from_facets

START, OUT = sys.argv[1], sys.argv[2]
SWEEPS, SEED = int(sys.argv[3]), int(sys.argv[4])
CELL = os.path.join(os.path.dirname(os.path.dirname(OUT)), "r_m4.mfd")
SAMPLE, SAVE_EVERY = 150, 2500
ZLEG, CIMP, NHINGE, LAM = 0.6, 1.0, 6.0, 1.0


def comps_of(fac):
    n6, imp, adj = vertex_classes(fac)
    V = len(n6)
    lab = None                    # dense == report space is fine here
    illv = [v for v in range(V) if imp[v] > 0]
    seen, comps = set(), []
    for s0 in illv:
        if s0 in seen:
            continue
        st, comp = [s0], []
        seen.add(s0)
        while st:
            u = st.pop(); comp.append(u)
            for w in adj[u]:
                if imp[w] > 0 and w not in seen:
                    seen.add(w); st.append(w)
        comps.append(sorted(int(x) for x in comp))   # py ints (json-safe)
    return comps


ddg.set_random_seed(SEED)
ref = ddg.Manifold.load(CELL, 3)
native = float(edges_from_facets(ref.facets())[1].mean())   # 5.104225
m = ddg.Manifold.load(START, 3)
params = ddg.SamplerParams(
    num_facets_target=ref.num_facets, num_facets_coef=0.1,
    hinge_degree_target=native, num_hinges_coef=NHINGE,     # BUMP = 0
    hinge_degree_variance_coef=0.0, codim3_degree_variance_coef=0.0,
    hinge_degree_target_coef=LAM * native / 6.0)
s = ddg.ManifoldSampler(m, params)
s.set_n6_potential(ZLEG * LAM, CIMP * LAM, tilt=[0.0] * 5)
v = s.manifold

out = open(OUT + ".release.jsonl", "w")
fac0 = np.asarray(v.facets())
c0 = comps_of(fac0)
sizes0 = sorted((len(c) for c in c0), reverse=True)
print(f"{os.path.basename(OUT)}: start n_illegal={sum(sizes0)} sizes={sizes0} "
      f"(target=native={native:.6f})", flush=True)
out.write(json.dumps({"sweep": 0, "n_illegal": sum(sizes0),
                      "sizes": sizes0,
                      "members": [c for c in c0]}) + "\n")
out.flush()

t0 = time.time()
done = 0
while done < SWEEPS:
    s.run(sweeps=SAMPLE)
    done += SAMPLE
    fac = np.asarray(v.facets())
    comps = comps_of(fac)
    sizes = sorted((len(c) for c in comps), reverse=True)
    out.write(json.dumps({"sweep": done, "t": round(time.time() - t0, 1),
                          "n_illegal": sum(sizes), "sizes": sizes,
                          "members": [c for c in comps]}) + "\n")
    out.flush()
    if done % SAVE_EVERY == 0:
        v.save(OUT + f"_rel{done}.mfd")
        print(f"{os.path.basename(OUT)} sweep{done} "
              f"({time.time()-t0:.0f}s): n_illegal={sum(sizes)} "
              f"sizes={sizes[:8]}", flush=True)
out.close()
v.save(OUT + "_release_final.mfd")
print(f"{os.path.basename(OUT)} DONE {done} sweeps in {time.time()-t0:.0f}s",
      flush=True)
