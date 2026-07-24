#!/usr/bin/env python3
"""Web-constrained curvature transport: jellium-neutralized flows on the six-web.

Physics (see notes/CONVENTIONS.md; discussion 2026-07-24): on the legal manifold
curvature IS web density (ebar = 5 + f6), so the question "can local moves anneal
the defect gas into the medium?" becomes a transport problem ON the six-web:

  * charges: per-COMPLEX cell-mean curvature charge Q_c, deposited on the
    complex's on-web members (off-web core vertices carry none -- the web docks
    at the complex surface);
  * jellium: the compensating uniform background source on every legal web
    vertex ("each legal region can absorb ~Q_tot/V6 of curvature demand by
    local composition shift");
  * L1 (min-cost) flow with Euclidean edge lengths = Wasserstein-1 optimal
    transport of curvature demand along the web.  Path decomposition gives the
    charge-weighted TRANSPORT-DISTANCE spectrum: distances ~ nearest-complex
    baseline => demand locally absorbable (guided local moves well-posed);
    a fat tail => demand misplaced at long wavelength (no local protocol fast).
  * L2 (min-norm) flow = discrete Coulomb field of the same charges (lsmr on
    the incidence matrix; underdetermined consistent => minimum-norm solution).

Chart consistency is verified (manifold edge set == cocycle edge set) before
any coordinate is trusted -- detached charts (see memory: cocycle-detachment)
are refused.

Usage: web_transport.py LABEL:snap.mfd [LABEL:snap.mfd ...] [--quantum QRAD]
       [--json OUT.json]
Needs sibling .cocycle.npz per snapshot.  Charges quantized to --quantum rad
(default 1e-3) for the integer min-cost solver; L2 uses exact real charges.
"""
import argparse
import heapq
import json
import os
import sys
from collections import defaultdict

import networkx as nx
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import lsmr

_HERE = os.path.dirname(os.path.abspath(__file__))
for p in ("../../python", "../../scripts"):
    sys.path.insert(0, os.path.join(_HERE, p))
import discrete_differential_geometry as ddg
from discrete_differential_geometry import cocycle as coc
from discrete_differential_geometry.vertex_fields import FIELDS
from fk_skeleton import edges_from_facets
from dopant_pairs import vertex_classes

CELL = 1.0e6                      # raw chart units per unit cell


def analyze(label, path, quantum):
    fac = np.asarray(ddg.Manifold.load(path, 3).facets())
    eu, edeg, Vn = edges_from_facets(fac)
    n6v, imp, adj = vertex_classes(fac)
    lab = np.unique(fac)

    # --- chart, with the detachment guard ---
    edges_c, omega, _ = coc.load_cocycle(path[:-4] + ".cocycle.npz")
    edges_c, omega = coc.canonicalize_labels(np.asarray(edges_c),
                                             np.asarray(omega))
    eset = {tuple(sorted(e)) for e in eu.tolist()}
    cset = {tuple(sorted(e)) for e in edges_c.tolist()}
    if eset != cset:
        raise RuntimeError(f"{label}: chart DETACHED ({len(eset ^ cset)} edge "
                           "mismatches) -- regenerate with regen_cocycle.py")
    frac, basis = coc.torus_positions(fac, edges_c, omega)
    P = np.abs(np.diag(basis))

    # --- six-web graph with Euclidean lengths (cells) ---
    six = np.nonzero(edeg == 6)[0]
    eu6 = eu[six]
    d = frac[eu6[:, 1]] - frac[eu6[:, 0]]
    d -= np.round(d)                              # min image (frac period 1)
    elen = np.linalg.norm(d @ basis, axis=1) / CELL          # cells
    web_verts = sorted({int(v) for e in eu6 for v in e})
    widx = {v: i for i, v in enumerate(web_verts)}
    V6, E6 = len(web_verts), len(eu6)

    # --- complexes and their charges (cell-mean reference) ---
    qR = FIELDS["curvature_charge"](fac)
    dq = qR - qR.mean()
    lab_to_i = {int(l): i for i, l in enumerate(lab)}
    illv = [i for i in range(len(lab)) if imp[i] > 0]
    seen, comps = set(), []
    for s0 in illv:
        if s0 in seen:
            continue
        st, comp = [s0], []
        seen.add(s0)
        while st:
            u = st.pop()
            comp.append(u)
            for w in adj[u]:
                if imp[w] > 0 and w not in seen:
                    seen.add(w)
                    st.append(w)
        comps.append(comp)
    Qc = [float(dq[c].sum()) for c in comps]
    onweb = [[int(lab[j]) for j in c if int(lab[j]) in widx] for c in comps]

    # --- integer charge deposit (quantum units) + jellium background ---
    b = np.zeros(V6, dtype=np.int64)
    sink_of = {}                                  # web vertex idx -> complex id
    skipped = 0
    for ci, (mem, q) in enumerate(zip(onweb, Qc)):
        qi = int(round(q / quantum))
        if not mem or qi == 0:
            skipped += 1 if qi != 0 else 0
            continue
        share, rem = divmod(abs(qi), len(mem))
        sgn = 1 if qi > 0 else -1
        for k, v in enumerate(mem):
            amt = sgn * (share + (1 if k < rem else 0))
            b[widx[v]] += amt
            if amt:
                sink_of[widx[v]] = ci
    # jellium must balance PER WEB COMPONENT (flux cannot cross components)
    parent = list(range(V6))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for a_, b_ in eu6:
        ra, rb = find(widx[int(a_)]), find(widx[int(b_)])
        if ra != rb:
            parent[ra] = rb
    comp_of = np.array([find(i) for i in range(V6)])
    legal_mask = np.array([imp[lab_to_i[v]] == 0 for v in web_verts])
    ncomp_web = len(set(comp_of.tolist()))
    for c in set(comp_of.tolist()):
        members = np.nonzero(comp_of == c)[0]
        tot = int(b[members].sum())
        if tot == 0:
            continue
        lm = [int(i) for i in members if legal_mask[i]] or [int(i) for i in members]
        base, rem = divmod(-tot, len(lm))
        for k, i in enumerate(lm):
            b[i] += base + (1 if k < rem else 0)
    assert b.sum() == 0

    # --- L1 min-cost flow (network simplex; integer costs in milli-cells) ---
    G = nx.DiGraph()
    for i in range(V6):
        G.add_node(i, demand=int(-b[i]))          # nx: demand = wanted inflow
    cost = np.maximum(1, np.round(elen * 1000).astype(int))
    for k, (a_, b_) in enumerate(eu6):
        i, j = widx[int(a_)], widx[int(b_)]
        G.add_edge(i, j, weight=int(cost[k]))
        G.add_edge(j, i, weight=int(cost[k]))
    flow = nx.min_cost_flow(G)
    w1 = nx.cost_of_flow(G, flow) / 1000.0 * quantum   # rad*cells

    # --- path stripping -> charge-weighted transport-distance spectrum ---
    lenof = {}
    for k, (a_, b_) in enumerate(eu6):
        i, j = widx[int(a_)], widx[int(b_)]
        lenof[(i, j)] = lenof[(j, i)] = elen[k]
    res = {u: dict(vs) for u, vs in flow.items()}
    sup = b.copy()
    paths = []                                     # (amount_int, length_cells, sink)
    srcs = [i for i in range(V6) if sup[i] > 0]
    for s in srcs:
        while sup[s] > 0:
            u, walk, plen, seenv = s, [s], 0.0, {s}
            while sup[u] >= 0:
                nxt = None
                for v, f in res[u].items():
                    if f > 0:
                        nxt = v
                        break
                if nxt is None or nxt in seenv:    # dead end / cycle guard
                    nxt = None
                    break
                plen += lenof[(u, nxt)]
                u = nxt
                walk.append(u)
                seenv.add(u)
            if nxt is None and sup[u] >= 0:
                break                              # cannot extend (numerical)
            amt = min(sup[s], -sup[u],
                      min(res[walk[i]][walk[i + 1]]
                          for i in range(len(walk) - 1)))
            for i in range(len(walk) - 1):
                res[walk[i]][walk[i + 1]] -= amt
            sup[s] -= amt
            sup[u] += amt
            paths.append((int(amt), plen, sink_of.get(u, -1)))

    # --- baseline: charge-weighted distance to NEAREST complex (Dijkstra) ---
    adj6 = defaultdict(list)
    for k, (a_, b_) in enumerate(eu6):
        i, j = widx[int(a_)], widx[int(b_)]
        adj6[i].append((j, elen[k]))
        adj6[j].append((i, elen[k]))
    dist = np.full(V6, np.inf)
    pq = []
    for i in sink_of:
        dist[i] = 0.0
        heapq.heappush(pq, (0.0, i))
    while pq:
        du, u = heapq.heappop(pq)
        if du > dist[u]:
            continue
        for v, L in adj6[u]:
            if du + L < dist[v]:
                dist[v] = du + L
                heapq.heappush(pq, (du + L, v))
    wsup = np.abs(b).astype(float)                 # background charge (either
    wsup[list(sink_of)] = 0.0                      # polarity), complexes excluded
    base_mean = float((dist * wsup).sum() / wsup.sum())

    # --- L2 min-norm (Coulomb) flow with exact real charges ---
    br = np.zeros(V6)
    for ci, (mem, q) in enumerate(zip(onweb, Qc)):
        if mem:
            for v in mem:
                br[widx[v]] += q / len(mem)
    for c in set(comp_of.tolist()):
        members = np.nonzero(comp_of == c)[0]
        lm = members[legal_mask[members]]
        lm = lm if len(lm) else members
        br[lm] -= br[members].sum() / len(lm)
    rows = np.concatenate([[widx[int(a_)] for a_, _ in eu6],
                           [widx[int(b_)] for _, b_ in eu6]])
    cols = np.concatenate([np.arange(E6), np.arange(E6)])
    vals = np.concatenate([-np.ones(E6), np.ones(E6)])
    B = sp.csr_matrix((vals, (rows, cols)), shape=(V6, E6))
    phi = lsmr(B, br, atol=1e-12, btol=1e-12, maxiter=20000)[0]
    resid = np.abs(B @ phi - br).max()

    # --- report ---
    amts = np.array([a for a, _, _ in paths], float)
    lens = np.array([L for _, L, _ in paths])
    order = np.argsort(lens)
    cw = np.cumsum(amts[order]) / amts.sum()

    def qtile(f):
        return float(lens[order][np.searchsorted(cw, f)])

    tmean = float((amts * lens).sum() / amts.sum())
    print(f"\n== {label}  ({os.path.basename(path)})")
    print(f"  web: V6={V6} E6={E6} components={ncomp_web}; complexes {len(comps)} "
          f"(on-web {sum(1 for m in onweb if m)}), Q_tot={sum(Qc):+.2f} rad, "
          f"quantum {quantum * 1000:.1f} mrad, skipped(no-web/zero) {skipped}")
    print(f"  W1 = {w1:.1f} rad*cells over {len(paths)} paths; stripped "
          f"{amts.sum() * quantum:.2f} rad of demand")
    print(f"  transport distance (charge-weighted): mean {tmean:.2f} cells, "
          f"median {qtile(0.5):.2f}, p90 {qtile(0.9):.2f}, max {lens.max():.2f}")
    print(f"  nearest-complex baseline (charge-weighted mean): "
          f"{base_mean:.2f} cells  -> transport/baseline = "
          f"{tmean / base_mean:.2f}")
    print(f"  L2 Coulomb: ||phi||^2={phi @ phi:.4f}, max|phi|="
          f"{np.abs(phi).max():.4f} rad, resid {resid:.2e}")
    return dict(label=label, path=path, V6=V6, E6=E6, n_cx=len(comps),
                Q_tot=sum(Qc), w1=w1, t_mean=tmean, t_med=qtile(0.5),
                t_p90=qtile(0.9), t_max=float(lens.max()),
                base_mean=base_mean, ratio=tmean / base_mean,
                l2_energy=float(phi @ phi), l2_max=float(np.abs(phi).max()))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("states", nargs="+", help="LABEL:snapshot.mfd")
    ap.add_argument("--quantum", type=float, default=1e-3,
                    help="charge quantum in rad for the integer solver")
    ap.add_argument("--json", default=None)
    args = ap.parse_args()
    out = []
    for arg in args.states:
        lab_, path = arg.split(":", 1)
        out.append(analyze(lab_, path, args.quantum))
    if args.json:
        with open(args.json, "w") as f:
            json.dump(out, f, indent=1)
        print(f"\nwrote {os.path.abspath(args.json)}")


if __name__ == "__main__":
    main()
