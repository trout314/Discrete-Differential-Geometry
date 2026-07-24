#!/usr/bin/env python3
"""Regenerate a T^3 integer cocycle for a mostly-crystalline state from scratch.

For a state whose cocycle record is lost or corrupt (see memory: cocycle
detachment, lam35 snap14000+): the manifold itself still fixes the physical
chart, because ~95% of its vertices develop single-valuedly onto the reference
crystal (crystal_grains covering map). Method:

  1. find_grains() -> per-tet vertex correspondence into the reference torus;
  2. every edge inside an assigned tet gets its displacement from the
     reference images (min-image in the reference box; exact closure on any
     triangle contained in one tet, single-valued within a grain);
  3. BFS-integrate those displacements -> a lift p(v) of every reachable
     vertex to the cover; unreachable (defect-pocket) vertices are attached to
     a positioned neighbor and relaxed (mean of positioned neighbors, a few
     Jacobi sweeps);
  4. build_from_positions(edges, p, MCELL) -> exactly-closed omega with the
     physical cohomology class (the registry chart); wrap indicators handle
     torus-crossing edges.

The regenerated omega is a FRESH GAUGE: valid for all coordinate observables
and future tracking, but absolute winding bookkeeping is not continuous with
the state's lost history.

Usage: regen_cocycle.py STATE.mfd REF.mfd MCELL OUT.cocycle.npz
"""
import os
import sys
from collections import Counter, defaultdict, deque

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
for p in ("../../python", "../../scripts"):
    sys.path.insert(0, os.path.join(_HERE, p))
import discrete_differential_geometry as ddg
from discrete_differential_geometry import cocycle as coc
from cocycle_check import reference_frac_positions
from crystal_grains import build_struct, find_grains, ref_index
from fk_skeleton import edges_from_facets

STATE, REF, MCELL, OUT = sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4]
BOX = float(MCELL)

fac = np.asarray(ddg.Manifold.load(STATE, 3).facets())
ref_fac = np.asarray(ddg.Manifold.load(REF, 3).facets())
st = build_struct(fac)
refs = {"r": build_struct(ref_fac)}
idx = ref_index(refs)
grain_of_tet, sig_of_tet, phase_of_grain = find_grains(st, refs, idx)
rp = np.asarray(reference_frac_positions("r", MCELL), float)   # (Vref, 3) cells

# ---- 2. edge displacements from assigned-tet reference images ----
def mim(d):
    return d - np.round(d / BOX) * BOX

disp = {}                      # (u, v) u<v -> displacement u->v (cells)
n_conflict = 0
gmain = Counter(grain_of_tet.values()).most_common(1)[0][0]
for t, g in grain_of_tet.items():
    if g != gmain:
        continue                # one grain = one frame; skip re-seed remnants
    sig = sig_of_tet[t]
    vs = list(sig.keys())
    for i in range(4):
        for j in range(i + 1, 4):
            u, w = vs[i], vs[j]
            d = mim(rp[sig[w]] - rp[sig[u]])
            key, dd = ((u, w), d) if u < w else ((w, u), -d)
            if key in disp:
                if not np.allclose(disp[key], dd, atol=1e-9):
                    n_conflict += 1
            else:
                disp[key] = dd

eu, edeg, V = edges_from_facets(fac)
edges = np.asarray(sorted(tuple(sorted(e)) for e in eu))
lab = np.unique(fac)
print(f"{os.path.basename(STATE)}: V={len(lab)} E={len(edges)}; "
      f"grain edges with displacement: {len(disp)} "
      f"({100 * len(disp) / len(edges):.2f}%), conflicts: {n_conflict}")

# ---- 3. BFS lift over grain edges, then relax the pockets ----
adj_d = defaultdict(list)                       # grain edges (with disp)
for (u, w), d in disp.items():
    adj_d[u].append((w, d))
    adj_d[w].append((u, -d))
adj_all = defaultdict(list)                     # all edges
for u, w in edges:
    adj_all[u].append(w)
    adj_all[w].append(u)

pos = {}
for root in lab:
    if int(root) in pos or int(root) not in adj_d:
        continue
    if pos:
        break                                    # only lift the main component
    pos[int(root)] = np.zeros(3)
    q = deque([int(root)])
    while q:
        u = q.popleft()
        for w, d in adj_d[u]:
            if w not in pos:
                pos[w] = pos[u] + d
                q.append(w)
print(f"BFS-lifted {len(pos)}/{len(lab)} vertices via grain edges")

# attach the rest to any positioned neighbor, then Jacobi-relax the pockets
todo = [int(v) for v in lab if int(v) not in pos]
changed = True
while changed:
    changed = False
    for v in todo:
        if v in pos:
            continue
        nb = [w for w in adj_all[v] if w in pos]
        if nb:
            pos[v] = np.mean([mim(pos[w] - pos[nb[0]]) + pos[nb[0]]
                              for w in nb], axis=0)
            changed = True
assert len(pos) == len(lab), f"unpositioned vertices: {len(lab) - len(pos)}"
for _ in range(50):                              # relax pocket interiors
    for v in todo:
        nb = adj_all[v]
        pos[v] = np.mean([mim(pos[w] - pos[v]) + pos[v] for w in nb], axis=0)

# ---- 4. exactly-closed omega via the standard wrap construction ----
frac = np.zeros((int(lab.max()) + 1, 3))
for v, p in pos.items():
    frac[v] = p
omega = coc.build_from_positions(edges, frac, BOX)
coc.save_cocycle(OUT, edges, omega, sweeps=-1)
print(f"wrote {OUT}")

# ---- validation: enable on the state, check, winding lattice ----
m = ddg.Manifold.load(STATE, 3)
params = ddg.SamplerParams(num_facets_target=m.num_facets)
s = ddg.ManifoldSampler(m, params)
s.enable_cocycle(edges, omega)
s.check_cocycle()
p3, cyc, _ = coc.tree_positions(edges, omega, int(lab.max()) + 1)
W = coc.lattice_basis(np.asarray(cyc))
print("check_cocycle: OK;  winding-lattice basis (raw units):")
print(W)
