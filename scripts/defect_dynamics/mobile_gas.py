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
  optional slideprob=P (anywhere in argv): enable the D-side knot-slide move
  type with probability P of proposing a slide once the unified proposal
  lands on a degree-3 edge (0 = off, the default). Slides are clean
  (species-preserving) by construction, so acceptance is plain Metropolis on
  the exact action change; they run inside the D sampler and participate in
  objective tracking and cocycle updates like any other move type.
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
import defect_state as ds

SLIDEPROB = next((float(x.split("=", 1)[1]) for x in sys.argv
                  if x.startswith("slideprob=")), 0.0)
a = [x for x in sys.argv if not x.startswith("slideprob=")]
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
if SLIDEPROB:
    s.set_slide_prob(SLIDEPROB)
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
    # Census via the shared incremental implementation (defect_state), which
    # reproduces the historical imp>0 fields EXACTLY -- verified field for
    # field against the block this replaced -- while adding the broadened
    # counts (illegal-edge incidence OR non-FK coordination) alongside rather
    # than redefining any existing column.
    dstate = ds.DefectState(v)
    c = ds.census(dstate, members=True)
    comps = c["members"]
    # per-complex centroids from the raw cocycle tree lift (the validated
    # pass1 position protocol), min-imaged about each complex's first member
    # so a complex straddling a lift seam can't split; aligned with `members`.
    cents = []
    if comps:
        e1, w1 = s.read_cocycle()
        e1 = np.asarray(e1)
        pos = coc.tree_positions(e1, np.asarray(w1), int(e1.max()) + 1)[0]
        pbox = 1.0e6 * MCELL         # raw units per box period (1e6/cell)
        for cc in comps:
            p0 = pos[cc[0]].astype(float)
            rel = (pos[cc].astype(float) - p0 + pbox / 2) % pbox - pbox / 2
            cents.append([round(float(x), 1)
                          for x in (p0 + rel.mean(0)) % pbox])
    rec = {
        "sweep": done, "t": round(time.time() - t0, 1),
        # --- historical columns, unchanged in meaning
        "n_illegal": c["n_illegal"], "sizes": c["sizes"], "members": comps,
        "cents": cents,
        "mean_edeg": c["mean_edeg"],
        "legaledge": c["legaledge"],
        "legalvert": c["legalvert"],
        # --- broadened defect definition (see defect_state)
        "n_defect_broad": c["n_defect_broad"],
        "sizes_broad": c["sizes_broad"],
        "n_nonfk_all_legal": c["n_nonfk_all_legal"],
        "n_nonfk_n6_1": c["n_nonfk_n6_1"],
        "n_nonfk_n6_ge5": c["n_nonfk_n6_ge5"],
        "legalvert_fk": c["legalvert_fk"]}
    if SLIDEPROB:
        st, sa = s.slide_stats()          # cumulative; difference per block
        rec["slide_tries"] = st
        rec["slide_accepts"] = sa
    log.write(json.dumps(rec) + "\n")
    log.flush()
    if (done - BURN) % SNAP == 0:
        nsnap += 1
        stem = f"{OUT}_snap{done}"
        v.save(stem + ".mfd")
        e1, w1 = s.read_cocycle()
        # canonicalize before persisting: the .mfd writer renumbers labels to
        # rank order; raw live labels (holes from 1<->4 churn) would skew the
        # saved pair (the "phantom detachment" postmortem, 2026-07-24).
        e1s, w1s = coc.canonicalize_labels(np.asarray(e1), np.asarray(w1))
        coc.save_cocycle(stem + ".cocycle.npz", e1s, w1s, sweeps=done)
        try:
            s.check_cocycle()
            drift = ""
        except Exception as e:
            drift = f" COCYCLE-DRIFT {e}"
        # strong guard (lam35 postmortem): check_cocycle can pass while the
        # tracked cocycle DETACHES from the manifold (edge set goes stale, no
        # error). Compare edge sets directly; a mismatch poisons cents + resume.
        eset = {tuple(sorted(e)) for e in np.asarray(v.simplices(1)).tolist()}
        cset = {tuple(sorted(e)) for e in np.asarray(e1).tolist()}
        if eset != cset:
            drift += f" COCYCLE-DETACHED({len(eset ^ cset)} edges)"
        print(f"[{os.path.basename(OUT)}] sw{done} ({time.time()-t0:.0f}s): "
              f"nill={c['n_illegal']} ncomp={len(c['sizes'])} "
              f"nonFK={c['n_nonfk_all_legal']} "
              f"<edeg>={c['mean_edeg']:.6f} snap{nsnap}{drift}", flush=True)
log.close()
print(f"[{os.path.basename(OUT)}] DONE {done} sweeps, {nsnap} snapshots, "
      f"{time.time()-t0:.0f}s", flush=True)
