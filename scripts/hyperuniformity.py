#!/usr/bin/env python3
"""Hyperuniformity test: is the curvature field more uniform than random?

Two estimators of large-scale fluctuation suppression, both against a
SHUFFLED baseline (same graph, same windows, field values permuted over
vertices -- kills spatial correlations, keeps the marginal distribution and
the window geometry, so the ratio isolates spatial organization):

1. Window variance: BFS balls grown to EXACTLY m vertices (fixed-size window,
   sidestepping the fluctuating-ball-volume problem); Var over random centers
   of the window charge Q = sum q(v), vs m. Thermal: Var/Var_shuffle -> const;
   hyperuniform: ratio decays as a power of m.
2. Spectral: graph-Laplacian eigenmodes, S(lambda_n) = |phi_n . dq|^2 with
   k ~ sqrt(lambda); S/S_shuffle -> 0 as k -> 0 iff hyperuniform.

Fields: q_R(v) = sum_{e ni v} delta_e / 2 (curvature charge; global total
pinned by the edge pin -- that kills only the strict k=0 mode, the test is the
APPROACH), and q_V(v) = D_v/4 (volume charge; total pinned by the facet pin).
"""
import os, sys, glob, json, time
import numpy as np
from scipy.sparse import coo_matrix, identity
from scipy.sparse.linalg import eigsh

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_ROOT, "python"))
from discrete_differential_geometry import Manifold

THETA = float(np.arccos(1.0 / 3.0))
RNG = np.random.default_rng(0)
NCENT, NSHUF = 300, 4

ENSEMBLES = [
    ("ext2e-3_N1e4",  "S3_N1e4_1e-1_ED5p0043_2_VDVs_2e-3",   6),
    ("ext2e-3_N31623","S3_N31623_1e-1_ED5p0043_2_VDVs_2e-3", 4),
    ("ext2e-3_N1e5",  "S3_N1e5_1e-1_ED5p0043_2_VDVs_2e-3",   3),
    ("glass8e-3_N31623","S3_N31623_1e-1_ED5p0043_4_VDVs_8e-3", 4),
    ("hub_b0_N1e4",   "S3_N1e4_1e-1_ED5p1043_2",             4),
]

def load_fields(path):
    m = Manifold.load(path, 3)
    F = np.asarray(m.facets(), np.int64)
    lab, inv = np.unique(F, return_inverse=True)
    T = inv.reshape(F.shape); V = len(lab)
    deg = np.bincount(T.ravel(), minlength=V)
    pr = [(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)]
    epairs = np.sort(np.vstack([T[:, [i, j]] for i, j in pr]), axis=1)
    eu, ecnt = np.unique(epairs, axis=0, return_counts=True)
    delta = 2*np.pi - THETA*ecnt
    qR = np.zeros(V)
    np.add.at(qR, eu[:, 0], delta/2); np.add.at(qR, eu[:, 1], delta/2)
    qV = deg/4.0
    # adjacency lists
    adj = [[] for _ in range(V)]
    for a, b in eu:
        adj[a].append(b); adj[b].append(a)
    adj = [np.array(x) for x in adj]
    return qR, qV, adj, eu, V

def bfs_order(adj, V, src, mmax):
    """Vertices in BFS order from src (shells shuffled), length mmax."""
    seen = np.zeros(V, bool); seen[src] = True
    order = [src]; frontier = [src]
    while len(order) < mmax and frontier:
        nxt = []
        for u in frontier:
            for w in adj[u]:
                if not seen[w]:
                    seen[w] = True; nxt.append(w)
        RNG.shuffle(nxt)
        order.extend(nxt); frontier = nxt
    return np.array(order[:mmax])

def window_variances(q, orders, mgrid):
    """Var over centers of prefix sums of q along each BFS order."""
    sums = np.array([np.cumsum(q[o])[mgrid - 1] for o in orders])  # (nc, nm)
    return sums.var(axis=0)

results = {}
for tag, stem, reps in ENSEMBLES:
    t0 = time.time()
    paths = sorted(glob.glob(f"{_ROOT}/seeds/{stem}_s0*.mfd"))[:reps]
    acc = None
    for path in paths:
        qR, qV, adj, eu, V = load_fields(path)
        mmax = V // 6
        mgrid = np.unique(np.geomspace(8, mmax, 12).astype(int))
        orders = [bfs_order(adj, V, int(RNG.integers(V)), mmax)
                  for _ in range(NCENT)]
        row = {}
        for fname, q in (("R", qR), ("V", qV)):
            vr = window_variances(q, orders, mgrid)
            vs = np.mean([window_variances(RNG.permutation(q), orders, mgrid)
                          for _ in range(NSHUF)], axis=0)
            row[fname] = (vr, vs)
        # spectral (first path per ensemble only -- the expensive part)
        if acc is None:
            ii = np.r_[eu[:, 0], eu[:, 1]]; jj = np.r_[eu[:, 1], eu[:, 0]]
            A = coo_matrix((np.ones(len(ii)), (ii, jj)), shape=(V, V)).tocsr()
            L = coo_matrix((A.sum(1).A1, (range(V), range(V))),
                           shape=(V, V)).tocsr() - A
            k = min(41, V - 2)
            lam, phi = eigsh(L, k=k, sigma=0, which="LM")
            o = np.argsort(lam); lam, phi = lam[o], phi[:, o]
            lam, phi = lam[1:], phi[:, 1:]              # drop constant mode
            spec = {}
            for fname, q in (("R", qR), ("V", qV)):
                dq = q - q.mean()
                Sr = (phi.T @ dq)**2
                Ss = np.mean([(phi.T @ (RNG.permutation(dq)))**2
                              for _ in range(8)], axis=0)
                spec[fname] = (lam.tolist(), (Sr/Ss).tolist())
        acc = acc or {"mgrid": mgrid.tolist(), "R": [], "V": []}
        for fname in ("R", "V"):
            vr, vs = row[fname]
            acc[fname].append((vr/vs).tolist())
    out = {"mgrid": acc["mgrid"], "spec": spec}
    for fname in ("R", "V"):
        M = np.array(acc[fname])
        out[fname] = dict(ratio_mean=M.mean(0).tolist(),
                          ratio_sem=(M.std(0, ddof=1)/np.sqrt(len(M))).tolist()
                          if len(M) > 1 else (M.std(0)*0).tolist())
    results[tag] = out
    r = out["R"]["ratio_mean"]
    print(f"[{tag:>18}] V-window ratio R: m={acc['mgrid'][0]}:{r[0]:.3f} ... "
          f"m={acc['mgrid'][-1]}:{r[-1]:.3f}   ({time.time()-t0:.0f}s)", flush=True)

OUT = os.path.join(_ROOT, "out"); os.makedirs(OUT, exist_ok=True)
json.dump(results, open(os.path.join(OUT,"hyperuniformity.json"), "w"), indent=1)
print("wrote out/hyperuniformity.json")
