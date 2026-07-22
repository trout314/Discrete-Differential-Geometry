#!/usr/bin/env python3
"""Mobile-gas production chain (m4): equilibrate the constrained knot liquid at
coupling scale LAM (the mobility window found by mobility_sweep, lam ~ 0.35) and
record everything the OCP/Stillinger-Lovett test needs:

  * ts.jsonl every TS sweeps: illegal complexes (members; sizes) -> mobility /
    turnover / carrier census at m4;
  * snapshot .mfd + .cocycle.npz every SNAP sweeps -> S_knot(k), charge S(k),
    charge-budget audit;
  * achieved mean edge degree each TS (pin satisfaction -> neutrality budget).

args: cell lam bump zleg cimp burn span ts snap mcell seed out [start] [startcoc]
  optional start/startcoc: resume/init from a snapshot .mfd + its .cocycle.npz
  (cell still supplies num_facets_target and the native degree).
"""
import json
import os
import sys
import time
from collections import Counter
from itertools import combinations

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
for p in ("../../python", "../../scripts", "../../tools"):
    sys.path.insert(0, os.path.join(_HERE, p))
import discrete_differential_geometry as ddg
from discrete_differential_geometry import cocycle as coc
from cocycle_check import reference_frac_positions
from fk_skeleton import edges_from_facets
from dopant_pairs import vertex_classes

a = sys.argv
CELL = a[1]; LAM = float(a[2]); BUMP = float(a[3])
ZLEG = float(a[4]); CIMP = float(a[5])
BURN = int(a[6]); SPAN = int(a[7]); TS = int(a[8]); SNAP = int(a[9])
MCELL = int(a[10]); SEED = int(a[11]); OUT = a[12]
START = a[13] if len(a) > 13 else None
STARTCOC = a[14] if len(a) > 14 else None

ddg.set_random_seed(SEED)
ref = ddg.Manifold.load(CELL, 3)
native = float(edges_from_facets(ref.facets())[1].mean())
et = native + BUMP
m = ddg.Manifold.load(START if START else CELL, 3)
params = ddg.SamplerParams(
    num_facets_target=ref.num_facets, num_facets_coef=0.1,
    hinge_degree_target=et, num_hinges_coef=0.0,
    hinge_degree_variance_coef=0.0, codim3_degree_variance_coef=0.0,
    hinge_degree_target_coef=LAM * et / 6.0)
s = ddg.ManifoldSampler(m, params)
s.set_n6_potential(ZLEG * LAM, CIMP * LAM, tilt=[0.0] * 5)
v = s.manifold
if STARTCOC:
    e0, w0, _ = coc.load_cocycle(STARTCOC)
    s.enable_cocycle(np.asarray(e0), np.asarray(w0))
else:
    edges = np.asarray(v.simplices(1))
    s.enable_cocycle(edges, coc.build_from_positions(
        edges, reference_frac_positions("r", MCELL), MCELL))

log = open(OUT + ".ts.jsonl", "w")
t0 = time.time()
s.run(sweeps=BURN)
print(f"[{os.path.basename(OUT)}] burned {BURN} ({time.time()-t0:.0f}s) "
      f"lam={LAM} target={et:.6f}", flush=True)

done = BURN
nsnap = 0
while done - BURN < SPAN:
    s.run(sweeps=TS)
    done += TS
    fac = np.asarray(v.facets())
    eu, edeg, V = edges_from_facets(fac)
    n6, imp, adj = vertex_classes(fac)
    lab = np.unique(fac)
    illv = [i for i in range(V) if imp[i] > 0]
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
        comps.append(sorted(int(lab[x]) for x in comp))
    sizes = sorted((len(c) for c in comps), reverse=True)
    log.write(json.dumps({
        "sweep": done, "t": round(time.time() - t0, 1),
        "n_illegal": int(sum(sizes)), "sizes": sizes, "members": comps,
        "mean_edeg": float(edeg.mean()),
        "legaledge": float(np.mean((edeg == 5) | (edeg == 6))),
        "legalvert": float(np.mean(imp == 0))}) + "\n")
    log.flush()
    if (done - BURN) % SNAP == 0:
        nsnap += 1
        stem = f"{OUT}_snap{done}"
        v.save(stem + ".mfd")
        e1, w1 = s.read_cocycle()
        coc.save_cocycle(stem + ".cocycle.npz", e1, w1, sweeps=done)
        try:
            s.check_cocycle()
            drift = ""
        except Exception as e:
            drift = f" COCYCLE-DRIFT {e}"
        print(f"[{os.path.basename(OUT)}] sw{done} ({time.time()-t0:.0f}s): "
              f"nill={sum(sizes)} ncomp={len(sizes)} "
              f"<edeg>={edeg.mean():.6f} snap{nsnap}{drift}", flush=True)
log.close()
print(f"[{os.path.basename(OUT)}] DONE {done} sweeps, {nsnap} snapshots, "
      f"{time.time()-t0:.0f}s", flush=True)
