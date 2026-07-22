#!/usr/bin/env python3
"""Produce a local-FK glass: melt an R crystal (free sampling, n6 off), then
quench under the FK-legality n6 potential (VDV/HDV off) at the flat edge target.

The FK n6 potential is zero iff every vertex is a Frank-Kasper coordination
(Z12/14/15/16) -- a purely LOCAL constraint with no crystalline registry. From a
fully disordered melt (no seed), a fast quench should heal toward FK-legal while
staying amorphous. We log FK-legality (edge/vertex) AND crystallinity
(crystal_grains registry fraction) separately: glass = high FK-legal + low
crystalline; failure = crystallizes or stays a disclination fluid.

args: cell out zleg cimp melt_sweeps quench_sweeps mcell seed
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
from discrete_differential_geometry import cocycle as coc
from cocycle_check import reference_frac_positions
from fk_skeleton import edges_from_facets, vertex_class_census
from dopant_pairs import vertex_classes
from defect_census import defect_cores

a = sys.argv
CELL, OUT = a[1], a[2]
ZLEG, CIMP = float(a[3]), float(a[4])
MELT, QUENCH, MCELL, SEED = int(a[5]), int(a[6]), int(a[7]), int(a[8])
LAM = 1.0
REF_FACETS = np.asarray(ddg.Manifold.load(CELL, 3).facets())

ddg.set_random_seed(SEED)
ref = ddg.Manifold.load(CELL, 3)
native = float(edges_from_facets(ref.facets())[1].mean())
m = ddg.Manifold.load(CELL, 3)
# Pin BOTH volume (f3) and mean edge degree (num_hinges) -> fixes f0 via Euler
# (f0 = f1 - f3 on T^3), so V is held while the degree VARIANCE is free to melt.
# hinge_degree_target_coef (per-edge pin) MUST be 0: it blocks melting and fights
# FK-legality (penalizes the legal degree-6 edges). The n6 potential is the only
# FK force. Variance penalties (VDV/HDV) off per project convention.
params = ddg.SamplerParams(
    num_facets_target=ref.num_facets, num_facets_coef=0.1,
    hinge_degree_target=native, num_hinges_coef=2.0,
    hinge_degree_variance_coef=0.0, codim3_degree_variance_coef=0.0,
    hinge_degree_target_coef=0.0)
s = ddg.ManifoldSampler(m, params)
v = s.manifold
edges = np.asarray(v.simplices(1))
s.enable_cocycle(edges, coc.build_from_positions(
    edges, reference_frac_positions("r", MCELL), MCELL))

log = open(OUT + ".glass.jsonl", "w")


def census():
    fac = np.asarray(v.facets())
    eu, edeg, V = edges_from_facets(fac)
    fz, _ = vertex_class_census(eu, edeg, V)
    n6, imp, adj = vertex_classes(fac)
    return dict(V=V, nfac=int(v.num_facets), edeg=float(edeg.mean()),
                legaledge=float(np.mean((edeg == 5) | (edeg == 6))),
                legalvert=float(1 - fz["impure"]),
                nill=int((imp > 0).sum()),
                edv=float(edeg.var()))


def logrow(phase, sweep, extra=None):
    r = census()
    r.update(phase=phase, sweep=sweep)
    if extra:
        r.update(extra)
    log.write(json.dumps(r) + "\n")
    log.flush()
    cs = f" cryst={extra['cryst']:.3f}" if extra else ""
    print(f"[{os.path.basename(OUT)}] {phase} sw{sweep}: nill={r['nill']} "
          f"legaledge={r['legaledge']:.3f} legalvert={r['legalvert']:.3f} "
          f"<edeg>={r['edeg']:.4f} edv={r['edv']:.3f}{cs}", flush=True)


t0 = time.time()
s.set_n6_potential(0.0, 0.0)                      # MELT
logrow("start", 0)
for k in range(MELT // 250):
    s.run(sweeps=250)
    logrow("melt", (k + 1) * 250)

s.set_n6_potential(ZLEG * LAM, CIMP * LAM, tilt=[0.0] * 5)   # QUENCH
done = MELT
for k in range(QUENCH // 250):
    s.run(sweeps=250)
    done += 250
    extra = None
    if done % 1000 == 0:
        fac = np.asarray(v.facets())
        dv, _, _, _, gs = defect_cores(fac, REF_FACETS, "r", 30)
        extra = {"cryst": 1 - len(dv) / len(np.unique(fac)),
                 "largest_grain": int(gs[0]) if gs else 0}
    logrow("quench", done, extra)

v.save(OUT + "_glass.mfd")
e1, w1 = s.read_cocycle()
coc.save_cocycle(OUT + "_glass.cocycle.npz", e1, w1, sweeps=done)
try:
    s.check_cocycle()
    drift = "ok"
except Exception as e:
    drift = f"DRIFT {e}"
print(f"[{os.path.basename(OUT)}] DONE {done} sweeps {time.time()-t0:.0f}s  "
      f"cocycle={drift}", flush=True)
