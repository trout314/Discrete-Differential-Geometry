#!/usr/bin/env python3
"""Anatomy of the pin-release survivor chains (graph-based; no coords needed).

For each survivor complex in the release finals:
  * order its intra-complex deg-4 edges into paths (the +2 disclination lines);
  * caps (path endpoints): six-edge count n6, full incident edge-degree profile,
    and adjacency to the native six-edge web;
  * interior chain vertices: same;
  * straightness: end-to-end 1-skeleton distance vs contour length;
  * off-chain members: how they attach to the chain;
  * far-crystal baseline: per-class six-edge stats for comparison.
"""
import json
import sys
from collections import Counter, defaultdict, deque

import numpy as np

sys.path.insert(0, "python"); sys.path.insert(0, "scripts")
from discrete_differential_geometry import Manifold
from discrete_differential_geometry.vertex_fields import edges_and_degrees
from dopant_pairs import vertex_classes

SP = sys.argv[1]


def bfs_dist(adj, src, dst, cap=30):
    if src == dst:
        return 0
    seen = {src}
    q = deque([(src, 0)])
    while q:
        u, d = q.popleft()
        for w in adj[u]:
            if w == dst:
                return d + 1
            if w not in seen and d + 1 < cap:
                seen.add(w)
                q.append((w, d + 1))
    return cap


for tag in ("above10", "above13", "below3"):
    rows = [json.loads(l) for l in open(f"{SP}/release/{tag}.release.jsonl")]
    fin = rows[-1]["members"]
    fac = np.asarray(Manifold.load(f"{SP}/release/{tag}_release_final.mfd", 3).facets())
    eu, ecnt, deg, V = edges_and_degrees(fac)
    n6, imp, adj = vertex_classes(fac)
    # per-vertex incident edge-degree profile
    inc_deg = defaultdict(Counter)
    for (a, b), c in zip(eu, ecnt):
        inc_deg[int(a)][int(c)] += 1
        inc_deg[int(b)][int(c)] += 1

    print(f"\n════ {tag} ════")
    # far-crystal baseline (legal vertices): six-edge count distribution
    legal = [v for v in range(V) if imp[v] == 0]
    cls_counts = Counter(int(n6[v]) for v in legal)
    print("crystal baseline n6 census (legal): "
          + "  ".join(f"Z{12 + (0 if k == 0 else k)}[n6={k}]:{v}"
                      for k, v in sorted(cls_counts.items())))

    for comp in fin:
        cset = set(comp)
        e4 = [(int(a), int(b)) for (a, b), c in zip(eu, ecnt)
              if c == 4 and int(a) in cset and int(b) in cset]
        # build deg-4 subgraph, order into paths
        g4 = defaultdict(list)
        for a, b in e4:
            g4[a].append(b)
            g4[b].append(a)
        caps = [v for v in g4 if len(g4[v]) == 1]
        paths = []
        used = set()
        for c0 in caps:
            if c0 in used:
                continue
            path = [c0]
            used.add(c0)
            cur, prev = c0, None
            while True:
                nxt = [w for w in g4[cur] if w != prev]
                if not nxt:
                    break
                prev, cur = cur, nxt[0]
                path.append(cur)
                used.add(cur)
                if len(g4[cur]) == 1:
                    break
            paths.append(path)
        onchain = set(g4)
        off = cset - onchain
        print(f"\n  survivor size {len(comp)}: {len(e4)} deg-4 edges -> "
              f"{len(paths)} path(s), {len(off)} off-chain members")
        for p in paths:
            L = len(p) - 1
            d = bfs_dist(adj, p[0], p[-1])
            print(f"    path len {L} edges, end-to-end skeleton dist {d} "
                  f"(straightness {d / max(L, 1):.2f})")
            for role, v in [("capA", p[0]), ("capB", p[-1])]:
                prof = inc_deg[v]
                print(f"      {role} v{v}: n6={int(n6[v])} illegal={bool(imp[v])} "
                      f" edges by deg: " + " ".join(f"{k}:{prof[k]}"
                                                    for k in sorted(prof)))
            mids = p[1:-1]
            if mids:
                mn6 = [int(n6[v]) for v in mids]
                m4 = [inc_deg[v][4] for v in mids]
                print(f"      interior ({len(mids)}): n6 vals {mn6}, "
                      f"deg4-incidence {m4}")
        if off:
            for v in sorted(off):
                # how does an off-chain illegal member attach?
                touch = [w for w in adj[v] if w in onchain]
                prof = inc_deg[v]
                print(f"      off-chain v{v}: n6={int(n6[v])} "
                      f"touches {len(touch)} chain verts; edges by deg: "
                      + " ".join(f"{k}:{prof[k]}" for k in sorted(prof)))
