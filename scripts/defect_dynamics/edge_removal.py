"""Targeted single-edge removal: the elementary f-sector surgery primitive.

Goal: given a target edge e = (u, v), find and (optionally) apply a short
move sequence that deletes e from the triangulation -- drive d(e) down to 3
with base-2->3s (each lowers d(e) by 1), then consume it with the axis 3->2.
A vertex removal is Z-4 of these plus a final 4->1; the general bilocal
proposer's "deep half-move" alphabet is built from exactly this operation.

Method: staged depth-limited DFS with exact rollback (the machinery that
cracked the vertex-collapse problem where pure greed dead-ends). Per stage,
search sequences of <= depth moves confined to the closed stars of u and v,
committing the best sequence that lowers d(e) (ties by dS); repeat. The
stage alphabet includes cage moves that do not touch e -- these unblock the
apex-adjacency dead-ends that stop greedy descent.

Survey mode (CLI): sample random target edges stratified by initial degree,
attempt removal with rollback (state restored after every attempt), and
report success rate, move count, net dS and barrier by starting degree --
the cost-of-surgery measurement for the bilocal half-move census.
"""
import argparse
import os
import random
import sys
import time
from collections import defaultdict
from itertools import combinations

_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, os.path.join(_ROOT, "python"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import discrete_differential_geometry as ddg     # noqa: E402
import worm_deg4_slide as W                      # noqa: E402


def _star_verts(L, vv):
    out = set()
    for t in L.v2t[vv]:
        out.update(t)
    return out


def _moves_for(L, e):
    """Stage alphabet for target edge e: direct reducers + cage moves in
    star(u) u star(v). Yields (center, cocenter) bistellar moves."""
    u, v = e
    region = _star_verts(L, u) | _star_verts(L, v)
    out = []
    # 3->2 on deg-3 edges inside the region (cleanup; includes e at d=3)
    for ed, d in list(L.edeg.items()):
        if d != 3 or not (ed[0] in region and ed[1] in region):
            continue
        tets = [t for t in L.v2t[ed[0]] if ed[1] in t]
        if len(tets) != 3:
            continue
        lk = sorted({x for t in tets for x in t} - set(ed))
        if len(lk) == 3:
            out.append((tuple(ed), tuple(lk)))
    # 2->3 on faces within the region (direct reducers = faces containing e)
    faces = set()
    for vv in (u, v):
        for t in L.v2t[vv]:
            for f in combinations(t, 3):
                if set(f) <= region:
                    faces.add(tuple(sorted(f)))
    for f in faces:
        tets = [t for t in L.v2t[f[0]] if set(f) <= set(t)]
        if len(tets) != 2:
            continue
        apex = sorted((set(tets[0]) | set(tets[1])) - set(f))
        if len(apex) != 2 or tuple(sorted(apex)) in L.edeg:
            continue
        out.append((f, tuple(apex)))
    return out


def remove_edge(L, obj, e, depth=3, nodecap=20000, keep=True):
    """Attempt to remove edge e from the live state.

    L: worm_deg4_slide.Live (driver-attached for exact objective tracking);
    obj: callable returning the current objective; keep: commit on success
    (False = measure only, always roll back).

    Returns dict: ok, moves (list of (center, cocenter)), n_moves, ds_net,
    barrier (max S - S_start along the committed path), nodes, stages.
    """
    e = tuple(sorted(e))
    o0 = obj()
    all_moves = []
    smax = o0
    total_nodes = 0
    stages = 0

    def gone():
        return e not in L.edeg

    while not gone():
        d0 = L.edeg[e]
        best = {"d": d0, "ds": float("inf"), "seq": None}
        nodes = [0]
        o_base = obj()

        def dfs(dep, seq):
            if nodes[0] > nodecap:
                return
            if gone() or L.edeg.get(e, 0) < best["d"]:
                dd = 0 if gone() else L.edeg[e]
                ds = obj() - o_base
                if dd < best["d"] or ds < best["ds"]:
                    best.update(d=dd, ds=ds, seq=list(seq))
            if dep == 0 or gone():
                return
            for cen, coc in _moves_for(L, e):
                nodes[0] += 1
                try:
                    L.do(cen, coc)
                except Exception:
                    continue
                seq.append((cen, coc))
                dfs(dep - 1, seq)
                seq.pop()
                L.do(coc, cen)
        dfs(depth, [])
        total_nodes += nodes[0]
        stages += 1
        if best["seq"] is None:
            break
        for cen, coc in best["seq"]:
            L.do(cen, coc)
            all_moves.append((cen, coc))
            smax = max(smax, obj())

    ok = gone()
    result = {"ok": ok, "moves": list(all_moves), "n_moves": len(all_moves),
              "ds_net": obj() - o0, "barrier": smax - o0,
              "nodes": total_nodes, "stages": stages}
    if not ok or not keep:
        for cen, coc in reversed(all_moves):
            L.do(coc, cen)
    return result


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--state", required=True)
    ap.add_argument("--etarget", type=float, required=True)
    ap.add_argument("--cimp", type=float, default=0.7)
    ap.add_argument("--hinges-coef", type=float, default=1.0)
    ap.add_argument("--facets-target", type=int, required=True)
    ap.add_argument("--n", type=int, default=30, help="target edges to try")
    ap.add_argument("--depth", type=int, default=3)
    ap.add_argument("--nodecap", type=int, default=20000)
    ap.add_argument("--seed", type=int, default=11)
    args = ap.parse_args()

    random.seed(args.seed)
    m = ddg.Manifold.load(args.state, 3)
    p = ddg.SamplerParams(
        num_facets_target=args.facets_target, num_facets_coef=0.1,
        hinge_degree_target=args.etarget, num_hinges_coef=args.hinges_coef,
        hinge_degree_variance_coef=0.0, codim3_degree_variance_coef=0.0,
        hinge_degree_target_coef=0.0, codim3_degree_target_coef=0.0)
    s = ddg.ManifoldSampler(m, p)
    s.set_n6_potential(0.0, args.cimp, tilt=[0.0] * 5)
    v = s.manifold
    L = W.Live(v, driver=s.do_bistellar_move)
    obj = lambda: s.current_objective    # noqa: E731
    o0 = obj()
    print(f"[edge_removal] state={os.path.basename(args.state)} "
          f"S0={o0:.3f} depth={args.depth} cap={args.nodecap}", flush=True)

    # stratified targets by initial degree
    by_d = defaultdict(list)
    for e, d in L.edeg.items():
        by_d[d].append(e)
    degrees = sorted(by_d)
    per = max(1, args.n // len(degrees))
    targets = []
    for d in degrees:
        random.shuffle(by_d[d])
        targets += [(d, e) for e in by_d[d][:per]]
    targets = targets[:args.n]

    stats = defaultdict(lambda: {"try": 0, "ok": 0, "moves": [], "ds": [],
                                 "bar": [], "sec": []})
    for d0, e in targets:
        t0 = time.time()
        r = remove_edge(L, obj, e, depth=args.depth,
                        nodecap=args.nodecap, keep=False)
        dt = time.time() - t0
        st = stats[d0]
        st["try"] += 1
        if r["ok"]:
            st["ok"] += 1
            st["moves"].append(r["n_moves"])
            st["ds"].append(r["ds_net"])
            st["bar"].append(r["barrier"])
        st["sec"].append(dt)
        drift = obj() - o0
        print(f"  d0={d0} e={e}: {'OK' if r['ok'] else 'FAIL'} "
              f"moves={r['n_moves']} dS={r['ds_net']:+.2f} "
              f"barrier={r['barrier']:+.2f} nodes={r['nodes']} "
              f"({dt:.0f}s, drift {drift:+.1e})", flush=True)

    print("\n=== survey by initial degree ===", flush=True)
    import numpy as np
    for d in sorted(stats):
        st = stats[d]
        if st["ok"]:
            print(f"  d0={d}: {st['ok']}/{st['try']} removed; "
                  f"moves med {np.median(st['moves']):.0f}; "
                  f"dS med {np.median(st['ds']):+.2f} "
                  f"[{min(st['ds']):+.2f},{max(st['ds']):+.2f}]; "
                  f"barrier med {np.median(st['bar']):+.2f}; "
                  f"{np.median(st['sec']):.0f}s", flush=True)
        else:
            print(f"  d0={d}: 0/{st['try']} removed "
                  f"({np.median(st['sec']):.0f}s median)", flush=True)


if __name__ == "__main__":
    main()
