#!/usr/bin/env python3
"""Minimum m^2 coefficient that keeps a defect gas DISPERSED.

The barrier stopping a flicker from entering a defect complex is 100% the
m^2 anti-clustering term (see notes/memory/flight-contact-barrier.md: pins
and geometry contribute exactly zero). We want that barrier as low as
possible so defects can collide and exchange momentum -- but not so low that
equilibrium collapses to one clump of high-degree edges. This finds the
crossover.

DESIGN

* Fixed defect population. nfc >= 30 pins f_3 so hard that 3->2 never
  accepts: the k fliers are exactly conserved and only their ARRANGEMENT
  can change. So this measures clustering at fixed particle number, which is
  the question, rather than a creation/annihilation equilibrium.

* TWO-SIDED. Every cimp is run from a CLUMPED start (all k fliers mutually
  adjacent, build_blob) and from a DISPERSED start (k fliers scattered with a
  minimum pairwise separation). If the two agree at late times, that is
  evidence of equilibrium; if they do not, the state is kinetically arrested
  and no single-start "stationarity" test can tell the difference -- exactly
  the false-stationarity trap of notes/memory/ecmc-blob-ab.md, where arm A
  passed a late-window slope test on all four observables while sitting at
  spread 3.57 with equilibrium near 4.0.

* FAST TRANSPORT ON. Defects are essentially immobile under local Pachner
  moves (notes/memory/defect-travel.md: median max excursion exactly 0.000
  cells for nv<=5). Without a transport channel the dispersed start simply
  cannot clump and the clumped start cannot spread, and every cimp would look
  "stable" from both sides. The diffusive nonlocal slide -- the fastest
  channel measured (16-43x the plain dynamics) -- is therefore ON throughout.
  It is f-preserving, so it pays no pin cost and does not perturb the
  ensemble it is sampling.

OBSERVABLES (all center-free, so they cannot depend on where the blob began)

    ncomp        number of connected defect complexes; k = fully dispersed,
                 1 = the single clump the user wants to avoid
    top1         size of the largest complex, in vertices
    max_edeg     largest edge degree anywhere -- "high degree edges" is the
                 concrete failure mode
    mean_sep     mean pairwise graph distance between complex centroids
    m2           sum_v m(v)^2, the clustering energy itself

Usage:
    python scripts/defect_dynamics/cimp_scan.py --cimps 0 0.05 0.1 0.2 0.35 0.5 \
        --sweeps 6000 --seeds 3 --out data/ecmc_ab/cimp_scan.json
"""
import argparse
import json
import os
import sys
from collections import Counter, deque

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(os.path.dirname(_HERE))
sys.path.insert(0, os.path.join(_ROOT, "python"))
sys.path.insert(0, _HERE)

import discrete_differential_geometry as ddg
from discrete_differential_geometry import Manifold, ManifoldSampler, SamplerParams
from discrete_differential_geometry.recording import _components
# Blob-construction helpers, absorbed from ecmc_ab.py when the ECMC A/B
# experiment was archived (2026-08 cleanup, Phase 2).
def _face_apex_map(m):
    f2a = {}
    for t in np.asarray(m.facets()):
        t = tuple(sorted(int(x) for x in t))
        for i in range(4):
            f2a.setdefault(t[:i] + t[i + 1:], []).append(t[i])
    return f2a


def _edge_set(m):
    return {tuple(sorted(map(int, e))) for e in m.simplices(1)}


def build_blob(m, k, rng, radius=4):
    """Insert k vertex-adjacent 2->3 fliers inside a graph ball; return
    (list of (face, chord), center vertex)."""
    # ball on the pristine crystal
    edges = _edge_set(m)
    V = int(np.asarray(m.facets()).max()) + 1
    adj = [[] for _ in range(V)]
    for a, b in edges:
        adj[a].append(b)
        adj[b].append(a)
    center = int(rng.integers(V))
    dist = np.full(V, -1, np.int32)
    dist[center] = 0
    dq = deque([center])
    while dq:
        u = dq.popleft()
        if dist[u] >= radius:
            continue
        for w in adj[u]:
            if dist[w] < 0:
                dist[w] = dist[u] + 1
                dq.append(w)
    ball = {v for v in range(V) if 0 <= dist[v] <= radius}

    placed, support = [], set()
    for _ in range(50 * k):                       # retry budget
        if len(placed) >= k:
            break
        f2a = _face_apex_map(m)
        edges = _edge_set(m)
        faces = [(f, ap) for f, ap in f2a.items() if len(ap) == 2]
        rng.shuffle(faces)
        for f, ap in faces:
            sup = set(f) | set(map(int, ap))
            if not sup <= ball:
                continue
            if placed and not (sup & support):
                continue                          # must touch the blob
            chord = tuple(sorted(map(int, ap)))
            if chord in edges:
                continue                          # 2->3 invalid
            m.do_bistellar_move(list(f), list(ap))
            placed.append((f, chord))
            support |= sup
            break
        else:
            break                                 # no candidate face at all
    if len(placed) < k:
        raise RuntimeError(f"only placed {len(placed)}/{k} fliers")
    return placed, center

REF = os.path.join(_ROOT, "data/tcp_reference/T3_C15_m3_N3672.mfd")


def scatter_once(m, k, rng, sep):
    """Place k fliers on `m` with pairwise support separation >= sep.

    Returns the placed list, or None if it could not fit all k. Moves COMMIT
    as they go, so a failed attempt leaves `m` dirty -- the caller must pass a
    freshly loaded manifold per attempt (see build_scattered)."""
    edges0 = _edge_set(m)
    V = int(np.asarray(m.facets()).max()) + 1
    adj = [[] for _ in range(V)]
    for a, b in edges0:
        adj[a].append(b)
        adj[b].append(a)

    def ball(src, r):
        dist = {s: 0 for s in src}
        dq = deque(src)
        while dq:
            u = dq.popleft()
            if dist[u] >= r:
                continue
            for w in adj[u]:
                if w not in dist:
                    dist[w] = dist[u] + 1
                    dq.append(w)
        return set(dist)

    placed, blocked = [], set()
    for _ in range(200 * k):
        if len(placed) >= k:
            return placed
        f2a = _face_apex_map(m)
        edges = _edge_set(m)
        faces = [(f, ap) for f, ap in f2a.items() if len(ap) == 2]
        rng.shuffle(faces)
        for f, ap in faces:
            sup = set(f) | set(map(int, ap))
            if sup & blocked:
                continue
            chord = tuple(sorted(map(int, ap)))
            if chord in edges:
                continue
            if not m.has_bistellar_move(list(f), list(ap)):
                continue
            m.do_bistellar_move(list(f), list(ap))
            placed.append((f, chord))
            blocked |= ball(sup, sep)
            break
        else:
            break
    return placed if len(placed) >= k else None


def build_scattered(load, k, rng, min_sep=2):
    """Dispersed counterpart to build_blob, which REQUIRES each new flier to
    touch the existing clump and so can only make a connected start.

    `load` is a zero-arg callable returning a FRESH manifold: separation is a
    CEILING and we back off until k fit, and each attempt needs a clean
    manifold because the placements commit as they go. C15 m3 has only V=648
    vertices with mean graph distance ~4, so blocking a radius-r ball around
    each 5-vertex support saturates it quickly -- measured: r=4 fits 2 fliers,
    r=3 fits 4, r=2 fits 10, so k=12 needs r<=1.

    Returns (manifold, placed, achieved_sep). Raises rather than quietly
    placing fewer -- a silently empty scattered start reads as "no defects",
    which looks like a result rather than a bug (it did: the first run
    reported ncomp 0.00 across the board)."""
    for sep in range(int(min_sep), -1, -1):
        m = load()
        placed = scatter_once(m, k, rng, sep)
        if placed is not None:
            return m, placed, sep
    raise RuntimeError(f"build_scattered: cannot place {k} fliers even at "
                       f"separation 0")


def observe(view):
    """Center-free clustering observables."""
    pairs, degs = view.illegal_edges()
    pairs = np.asarray(pairs).reshape(-1, 2)
    degs = np.asarray(degs)
    if not len(pairs):
        return {"n_ill": 0, "n3": 0, "ncomp": 0, "top1": 0, "m2": 0,
                "max_edeg": 0, "mean_sep": 0.0}
    mv = Counter()
    for a, b in pairs:
        mv[int(a)] += 1
        mv[int(b)] += 1
    comps = _components(pairs)

    # mean pairwise graph distance between complex centroids, via BFS from one
    # representative vertex of each complex on the CURRENT 1-skeleton
    E = np.asarray(view.simplices(1))
    V = int(E.max()) + 1
    adj = [[] for _ in range(V)]
    for a, b in E:
        adj[int(a)].append(int(b))
        adj[int(b)].append(int(a))
    reps = []
    seen = set()
    for a, b in pairs:
        c = None
        for v in (int(a), int(b)):
            if v not in seen:
                c = v
                break
        if c is not None:
            reps.append(c)
            seen.add(c)
    reps = reps[:len(comps)] if comps else []
    ds = []
    for i, r in enumerate(reps):
        dist = np.full(V, -1, np.int32)
        dist[r] = 0
        dq = deque([r])
        while dq:
            u = dq.popleft()
            for w in adj[u]:
                if dist[w] < 0:
                    dist[w] = dist[u] + 1
                    dq.append(w)
        for q in reps[i + 1:]:
            if dist[q] >= 0:
                ds.append(int(dist[q]))
    return {"n_ill": int(len(degs)), "n3": int((degs == 3).sum()),
            "ncomp": len(comps), "top1": int(comps[0]) if comps else 0,
            "m2": int(sum(c * c for c in mv.values())),
            "max_edeg": int(degs.max()),
            "mean_sep": float(np.mean(ds)) if ds else 0.0}


def run_one(start, cimp, seed, cfg):
    rng = np.random.default_rng(seed)
    ddg.set_random_seed(int(rng.integers(2 ** 31)))
    if start == "clumped":
        m = Manifold.load(REF, 3)
        build_blob(m, cfg["k_fliers"], rng, radius=cfg["radius"])
    else:
        m, _, _ = build_scattered(lambda: Manifold.load(REF, 3),
                                  cfg["k_fliers"], rng,
                                  min_sep=cfg["min_sep"])
    fv = [int(x) for x in m.f_vector]
    estar = 6.0 * fv[3] / fv[1]
    s = ManifoldSampler(m, SamplerParams(
        num_facets_target=fv[3], num_facets_coef=cfg["nfc"],
        hinge_degree_target=estar, num_hinges_coef=cfg["k1"],
        hinge_degree_variance_coef=0.0, codim3_degree_variance_coef=0.0,
        hinge_degree_target_coef=0.0))
    s.set_n6_potential(0.0, cimp, None, cfg.get("imp_offset", 0),
                       cfg.get("imp_lin", 0.0))
    s.set_nonlocal_slide_prob(cfg["nslide_prob"], cfg["nslide_max_step"])

    series = [{"sw": 0, **observe(s.manifold)}]
    nch = cfg["sweeps"] // cfg["chunk"]
    for ci in range(nch):
        s.run(sweeps=cfg["chunk"])
        series.append({"sw": (ci + 1) * cfg["chunk"], **observe(s.manifold)})
    return series


def late(series, frac=0.5):
    cut = series[-1]["sw"] * frac
    tail = [r for r in series if r["sw"] >= cut]
    return {k: float(np.mean([r[k] for r in tail]))
            for k in ("ncomp", "top1", "m2", "max_edeg", "mean_sep", "n_ill")}


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--cimps", type=float, nargs="+",
                    default=[0.0, 0.05, 0.1, 0.2, 0.35, 0.5])
    ap.add_argument("--seeds", type=int, default=3)
    ap.add_argument("--seed0", type=int, default=5000)
    ap.add_argument("--sweeps", type=int, default=6000)
    ap.add_argument("--chunk", type=int, default=100)
    ap.add_argument("--k-fliers", type=int, default=12)
    ap.add_argument("--radius", type=int, default=4)
    ap.add_argument("--min-sep", type=int, default=4)
    ap.add_argument("--imp-lin", type=float, default=0.0,
                   help="pure chemical potential on impure edges: adds "
                        "imp_lin*m, i.e. 2*imp_lin per impure edge, with NO "
                        "arrangement dependence")
    ap.add_argument("--imp-offset", type=int, default=0,
                   help="flat foot of the impurity quadratic: "
                        "V(m) = cimp * max(0, m - offset)^2. 0 = bare quadratic.")
    ap.add_argument("--nfc", type=float, default=30.0)
    ap.add_argument("--k1", type=float, default=1.0)
    ap.add_argument("--nslide-prob", type=float, default=0.15)
    ap.add_argument("--nslide-max-step", type=int, default=30)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    cfg = vars(args)

    out = {}
    print(f"[V(m) = {args.imp_lin}*m + cimp*max(0, m - {args.imp_offset})^2]")
    print(f"{'cimp':>6} {'start':>10} {'ncomp':>7} {'top1':>7} {'max_ed':>7} "
          f"{'mean_sep':>9} {'m2':>7} {'n_ill':>7}")
    for cimp in args.cimps:
        out[str(cimp)] = {}
        for start in ("clumped", "dispersed"):
            runs = []
            for i in range(args.seeds):
                try:
                    runs.append(run_one(start, cimp, args.seed0 + i, cfg))
                except Exception as ex:                 # noqa: BLE001
                    print(f"  FAIL cimp={cimp} {start} seed {i}: {ex}")
            if not runs:
                continue
            L = [late(r) for r in runs]
            agg = {k: (float(np.mean([x[k] for x in L])),
                       float(np.std([x[k] for x in L], ddof=1)
                             / np.sqrt(len(L))) if len(L) > 1 else 0.0)
                   for k in L[0]}
            out[str(cimp)][start] = {"late": agg,
                                     "series": [r for r in runs]}
            print(f"{cimp:>6.2f} {start:>10} {agg['ncomp'][0]:>7.2f} "
                  f"{agg['top1'][0]:>7.2f} {agg['max_edeg'][0]:>7.2f} "
                  f"{agg['mean_sep'][0]:>9.2f} {agg['m2'][0]:>7.1f} "
                  f"{agg['n_ill'][0]:>7.1f}", flush=True)
        # two-sided agreement: do the clumped and dispersed starts meet?
        g = out[str(cimp)]
        if "clumped" in g and "dispersed" in g:
            dn = abs(g["clumped"]["late"]["ncomp"][0]
                     - g["dispersed"]["late"]["ncomp"][0])
            se = (g["clumped"]["late"]["ncomp"][1] ** 2
                  + g["dispersed"]["late"]["ncomp"][1] ** 2) ** 0.5
            g["two_sided_gap_ncomp"] = dn
            g["two_sided_z"] = dn / se if se > 0 else float("inf")
            print(f"       -> two-sided gap in ncomp = {dn:.2f}"
                  f" ({'MET' if dn < 1.0 else 'NOT met -- arrested?'})")

    if args.out:
        os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
        with open(args.out, "w") as f:
            json.dump({"scan": out, "cfg": cfg}, f, indent=1, default=float)
        print(f"\nwrote {os.path.abspath(args.out)}")


if __name__ == "__main__":
    main()
