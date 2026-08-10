#!/usr/bin/env python3
"""Can an edge contraction / vertex split ever map an FK-legal state to another?

The 2<->3 Pachner move provably cannot (every one is born with a degree-3
edge), which is why FK-legal states are isolated points for the sampler
(notes/memory/fk-move-search.md). The contract/split channel is a genuinely
different move -- it changes f0 and it is a composite, block-scale operation --
so the question has to be asked again from scratch.

This script answers it two ways.

ANALYTIC.  Contract edge uv of degree r, ring x_1..x_r (the link cycle).
The r tetrahedra around uv are deleted, u and v merge to w.  Exactly three
kinds of edge change degree:

    ring edges   (x_i, x_i+1)       deg -> deg - 1
    spoke pairs  (u,x_i),(v,x_i)    merge, deg -> deg(u,x_i) + deg(v,x_i) - 4
    edge uv                          deleted

Everything else is untouched.  Starting FK-legal (all degrees in {5,6}) and
demanding the result be edge-legal forces

    ring edges were 6  (6 -> 5),  spoke pairs were 5 + 5  (-> 6),

since 5+6-4 = 7 and 6+6-4 = 8 are both illegal.  But then ALL r merged spokes
at w have degree 6, so n6(w) >= r, and r = deg(uv) must itself have been legal,
i.e. r in {5,6}.  Hence n6(w) >= 5: w is a Z17-or-higher hub, not a
Frank-Kasper coordination.  So:

    NO edge contraction maps FK-legal -> FK-legal,
    and by reversibility no vertex split does either.

Edge legality and vertex FK-ness are traded against each other, exactly the
failure mode the old VDV anneal hit (all-{5,6} edges, but hubs instead of
Z12/14/15/16).

NUMERIC.  The above rests on the three degree-change rules, so this script
verifies them against a literal contraction (relabel v -> u, drop degenerate
facets, recount) on real reference crystals, then censuses every contractible
edge of each crystal for the resulting damage:

    n_ill_e   illegal edges created (degree outside {5,6})
    n_nonfk   non-FK vertices created (m > 0, or n6 outside {0,2,3,4})

The minimum over all sites is the illegality BUDGET a worm-style sampler must
be allowed to spend for this channel to connect legal states at all.

Usage:
    python scripts/defect_dynamics/fk_channel_census.py                 # table
    python scripts/defect_dynamics/fk_channel_census.py --crystal a15 --verify
"""
import argparse
import math
import os
import sys
from collections import Counter, defaultdict

_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(_ROOT, "python"))
import numpy as np
import discrete_differential_geometry as ddg

REF = os.path.join(_ROOT, "data", "tcp_reference")
E_FLAT = 2.0 * math.pi / math.acos(1.0 / 3.0)

# name -> (file, mean coordination Zbar of the ideal structure)
LIBRARY = {
    "a15":   ("T3_A15_m3_N1242.mfd",    13.5),
    "c15":   ("T3_C15_m3_N3672.mfd",    40.0 / 3.0),
    "c14":   ("T3_C14_m3_N1836.mfd",    40.0 / 3.0),
    "c36":   ("T3_C36_m4_N8704.mfd",    40.0 / 3.0),
    "sigma": ("T3_SIGMA_m3_N4644.mfd",  13.466666666666667),
    "z":     ("T3_Z_m6_N8640.mfd",      13.428571428571429),
    "mu":    ("T3_MU_m4_N14208.mfd",    13.384615384615385),
    "r":     ("T3_R_m2_N7248.mfd",      13.397129186602871),
    "p":     ("T3_P_m3_N8640.mfd",      13.466666666666667),
    "delta": ("T3_DELTA_m3_N8640.mfd",  13.466666666666667),
}

LEGAL = (5, 6)
FK_N6 = (0, 2, 3, 4)


# --- combinatorial helpers (pure python over the facet array) ---------------

def edge_degrees(facets):
    """dict frozenset({u,v}) -> #facets containing the edge."""
    deg = Counter()
    for f in facets:
        a, b, c, d = f
        deg[(a, b)] += 1
        deg[(a, c)] += 1
        deg[(a, d)] += 1
        deg[(b, c)] += 1
        deg[(b, d)] += 1
        deg[(c, d)] += 1
    return deg


def key(u, v):
    return (u, v) if u < v else (v, u)


def star_map(facets):
    """vertex -> list of facet indices."""
    st = defaultdict(list)
    for i, f in enumerate(facets):
        for x in f:
            st[x].append(i)
    return st


def link_pieces(facets, star, v):
    """(vertex set, edge set) of the link of vertex v."""
    vs, es = set(), set()
    for i in star[v]:
        t = [x for x in facets[i] if x != v]
        vs.update(t)
        es.add(key(t[0], t[1]))
        es.add(key(t[0], t[2]))
        es.add(key(t[1], t[2]))
    return vs, es


def ring_of(facets, star, u, v):
    """(ring vertices, ring edges) = the link of edge uv, or None if uv absent."""
    rv, re = set(), set()
    n = 0
    for i in star[u]:
        f = facets[i]
        if v in f:
            n += 1
            x, y = [z for z in f if z != u and z != v]
            rv.add(x)
            rv.add(y)
            re.add(key(x, y))
    return (rv, re, n) if n else None


def link_condition(facets, star, u, v, rv, re):
    """lk(u) cap lk(v) == lk(uv): the contraction yields a manifold."""
    vu, eu = link_pieces(facets, star, u)
    vv, ev = link_pieces(facets, star, v)
    vu.discard(v)
    vv.discard(u)
    eu = {e for e in eu if v not in e}
    ev = {e for e in ev if u not in e}
    return (vu & vv) == rv and (eu & ev) == re


def n6_and_m(deg, incident):
    """(n6, m) from a vertex's incident-edge degree list."""
    n6 = sum(1 for d in incident if d >= 6)
    m = sum(1 for d in incident if d not in LEGAL)
    return n6, m


# --- the census -------------------------------------------------------------

def contraction_damage(facets, star, deg, adj, u, v):
    """Damage of contracting uv, computed locally. None if not contractible.

    Returns (n_ill_e, n_nonfk, r) -- illegal edges and non-FK vertices in the
    resulting complex, assuming the input was fully FK-legal.
    """
    ring = ring_of(facets, star, u, v)
    if ring is None:
        return None
    rv, re, r = ring
    if len(rv) != r or len(re) != r:
        return None                      # link of uv is not a simple r-cycle
    if not link_condition(facets, star, u, v, rv, re):
        return None

    # new degrees on the affected edges
    new = {}
    for e in re:
        new[e] = deg[e] - 1                       # ring edge
    for x in rv:
        new[("w", x)] = deg[key(u, x)] + deg[key(v, x)] - 4   # merged spoke

    # unaffected edges at w (neighbours of exactly one of u, v)
    others = {}
    for x in adj[u]:
        if x != v and x not in rv:
            others[x] = deg[key(u, x)]
    for x in adj[v]:
        if x != u and x not in rv:
            others[x] = deg[key(v, x)]

    # --- illegal-edge delta over the affected edges -------------------------
    after_ill = sum(1 for d in new.values() if d not in LEGAL)
    before_ill = sum(1 for e in re if deg[e] not in LEGAL)
    before_ill += sum(1 for x in rv
                      for e in (key(u, x), key(v, x)) if deg[e] not in LEGAL)
    before_ill += 0 if deg[key(u, v)] in LEGAL else 1

    # --- non-FK vertex delta over the affected vertices ---------------------
    def bad_before(x):
        n6, m = n6_and_m(deg, [deg[key(x, y)] for y in adj[x]])
        return 1 if (m > 0 or n6 not in FK_N6) else 0

    before_bad = bad_before(u) + bad_before(v) + sum(bad_before(x) for x in rv)

    inc_w = [new[("w", x)] for x in rv] + list(others.values())
    n6, m = n6_and_m(deg, inc_w)
    after_bad = 1 if (m > 0 or n6 not in FK_N6) else 0
    for x in rv:
        inc = [new[("w", x)]]
        for y in adj[x]:
            if y in (u, v):
                continue
            e = key(x, y)
            inc.append(new.get(e, deg[e]))
        n6, m = n6_and_m(deg, inc)
        if m > 0 or n6 not in FK_N6:
            after_bad += 1
    return after_ill - before_ill, after_bad - before_bad, r


def hub_splits(facets, w):
    """All vertex splits of w that are edge-legal AND restore full FK-legality.

    Cutting link(w) along a simple cycle gamma sends
        ambient (x_i, x_i+1) along gamma:   deg -> deg + 1
        spoke   (w, x):                     deg -> (a, b) with a + b = deg + 4
    so edge-legality forces every gamma edge to have had degree 5 and every
    spoke on gamma to have had degree 6 (6 + 4 = 10 = 5 + 5; a degree-5 spoke
    would need 9 = 5 + 4, illegal).  Hence gamma runs entirely along w's
    six-web spokes -- which is why only a hub (n6 >= 5) can be split at all.

    Returns a list of (gamma, new_facets).
    """
    star = star_map(facets)
    deg = edge_degrees(facets)
    tris = [tuple(x for x in facets[i] if x != w) for i in star[w]]
    lk_adj = defaultdict(set)
    for t in tris:
        for a, b in ((t[0], t[1]), (t[0], t[2]), (t[1], t[2])):
            lk_adj[a].add(b)
            lk_adj[b].add(a)

    six = sorted(x for x in lk_adj if deg[key(w, x)] >= 6)
    if len(six) < 5:
        return []
    six_set = set(six)

    # gamma = a simple cycle of length 5 or 6 among w's six-web spokes, riding
    # only ambient degree-5 link edges.  The new edge has degree |gamma|, so
    # only lengths 5 and 6 are legal.
    cycles = []

    def walk(path, seen, start):
        if len(path) in (5, 6):
            if start in lk_adj[path[-1]] and deg[key(path[-1], start)] == 5:
                cycles.append(tuple(path))
        if len(path) >= 6:
            return
        for y in sorted(lk_adj[path[-1]] & six_set):
            if y in seen or y < start or deg[key(path[-1], y)] != 5:
                continue
            walk(path + [y], seen | {y}, start)

    for s in six:
        walk([s], {s}, s)
    # each cycle is found twice (both traversal directions); canonicalise by
    # rotating to the minimum vertex and taking the smaller direction
    uniq = {}
    for c in cycles:
        i = c.index(min(c))
        fwd = c[i:] + c[:i]
        rev = (fwd[0],) + tuple(reversed(fwd[1:]))
        uniq[min(fwd, rev)] = c

    out = []
    for gamma in uniq.values():
        nf = apply_split(facets, w, gamma, tris)
        if nf is None:
            continue
        nd = edge_degrees(nf)
        if any(d not in LEGAL for d in nd.values()):
            continue
        nadj = defaultdict(set)
        for (a, b) in nd:
            nadj[a].add(b)
            nadj[b].add(a)
        ok = True
        for x in nadj:
            n6, mm = n6_and_m(nd, [nd[key(x, y)] for y in nadj[x]])
            if mm > 0 or n6 not in FK_N6:
                ok = False
                break
        if ok:
            out.append((gamma, nf))
    return out


def apply_split(facets, w, gamma, tris):
    """Split w along the link cycle gamma; the fresh vertex is max+1.

    gamma separates link(w) (an S^2) into two open disks; triangles on side A
    keep w, triangles on side B take the new vertex w2, and each gamma edge
    (x, y) gains the tetrahedron (w, w2, x, y).
    """
    gset = set(gamma)
    gedges = {key(gamma[i], gamma[(i + 1) % len(gamma)]) for i in range(len(gamma))}
    # 2-colour the link triangles across gamma: adjacency across link EDGES
    # that are not gamma edges.
    by_edge = defaultdict(list)
    for i, t in enumerate(tris):
        for a, b in ((t[0], t[1]), (t[0], t[2]), (t[1], t[2])):
            by_edge[key(a, b)].append(i)
    side = {}
    stack = [0]
    side[0] = 0
    while stack:
        i = stack.pop()
        t = tris[i]
        for a, b in ((t[0], t[1]), (t[0], t[2]), (t[1], t[2])):
            e = key(a, b)
            if e in gedges:
                continue
            for j in by_edge[e]:
                if j not in side:
                    side[j] = side[i]
                    stack.append(j)
    if len(side) != len(tris):
        # gamma did not separate into exactly two pieces in one flood: colour
        # the remainder and require exactly two components
        comp = 1
        for i in range(len(tris)):
            if i in side:
                continue
            comp += 1
            if comp > 2:
                return None
            stack = [i]
            side[i] = 1
            while stack:
                j = stack.pop()
                t = tris[j]
                for a, b in ((t[0], t[1]), (t[0], t[2]), (t[1], t[2])):
                    e = key(a, b)
                    if e in gedges:
                        continue
                    for k in by_edge[e]:
                        if k not in side:
                            side[k] = 1
                            stack.append(k)
        if comp != 2:
            return None
    if len(set(side.values())) != 2:
        return None

    w2 = max(max(f) for f in facets) + 1
    tri_side = {tris[i]: side[i] for i in range(len(tris))}
    out = []
    for f in facets:
        if w not in f:
            out.append(f)
            continue
        t = tuple(x for x in f if x != w)
        if tri_side[t] == 0:
            out.append(f)
        else:
            out.append(tuple(sorted(list(t) + [w2])))
    for e in gedges:
        out.append(tuple(sorted([w, w2, e[0], e[1]])))
    return out


def literal_contract(facets, u, v):
    """Actually contract: relabel v -> u, drop degenerate facets."""
    out = []
    for f in facets:
        if u in f and v in f:
            continue
        g = tuple(sorted(u if x == v else x for x in f))
        out.append(g)
    return out


def state_probe(path):
    """Damage spectrum of the contraction channel on an arbitrary state.

    Damage is a DELTA (after - before) over the affected edges/vertices, so a
    move with (0, 0) leaves the illegal-edge count and the non-FK vertex count
    exactly where they were: a free move on the constraint surface.
    """
    m = ddg.Manifold.load(path, 3)
    facets = [tuple(sorted(int(x) for x in f)) for f in m.facets()]
    deg = edge_degrees(facets)
    star = star_map(facets)
    adj = defaultdict(set)
    for (a, b) in deg:
        adj[a].add(b)
        adj[b].add(a)

    f0, f1, f3 = len(adj), len(deg), len(facets)
    n_ill_e = sum(1 for d in deg.values() if d not in LEGAL)
    n_bad_v = 0
    for x in adj:
        n6, mm = n6_and_m(deg, [deg[key(x, y)] for y in adj[x]])
        if mm > 0 or n6 not in FK_N6:
            n_bad_v += 1

    print(f"\n=== {os.path.basename(path)}")
    print(f"    f0={f0} f1={f1} f3={f3}  e_mean={6.0 * f3 / f1:.6f} "
          f"(e*={E_FLAT:.6f})  Zbar={2.0 * f1 / f0:.4f}")
    print(f"    illegal edges {n_ill_e} ({100.0 * (1 - n_ill_e / f1):.3f}% legal)"
          f"   non-FK vertices {n_bad_v} "
          f"({100.0 * (1 - n_bad_v / f0):.3f}% FK)")

    res = []
    for (u, v) in sorted(deg):
        d = contraction_damage(facets, star, deg, adj, u, v)
        if d is not None:
            res.append(d)
    if not res:
        print("    no contractible edges")
        return
    ill = Counter(t[0] for t in res)
    joint = Counter((t[0], t[1]) for t in res)
    free = sum(n for (a, b), n in joint.items() if a <= 0 and b <= 0)
    print(f"    contractible {len(res)}/{f1}; d(illegal edges) histogram "
          f"{sorted(ill.items())[:8]}")
    print(f"    moves with d_ill<=0 AND d_nonFK<=0: {free} "
          f"({100.0 * free / len(res):.2f}% of contractible)")
    print(f"    best joint (d_ill, d_nonFK): "
          f"{sorted(joint.items())[:6]}")


def worm_probe(args):
    """contract (0 illegal edges, 1 hub) -> split back: how many ways?"""
    names = [args.crystal] if args.crystal else sorted(LIBRARY)
    for name in names:
        fname, _ = LIBRARY[name]
        m = ddg.Manifold.load(os.path.join(REF, fname), 3)
        facets = [tuple(sorted(int(x) for x in f)) for f in m.facets()]
        deg = edge_degrees(facets)
        star = star_map(facets)
        adj = defaultdict(set)
        for (a, b) in deg:
            adj[a].add(b)
            adj[b].add(a)

        clean = []
        for (u, v) in sorted(deg):
            d = contraction_damage(facets, star, deg, adj, u, v)
            if d is not None and d[0] == 0:
                clean.append(((u, v), d))
        print(f"\n=== {name}: {len(clean)} / {len(deg)} edges contract with "
              f"ZERO illegal edges ({100.0 * len(clean) / len(deg):.2f}%)")
        if not clean:
            continue
        print(f"    hub damage histogram (n_nonfk, r): "
              f"{Counter((d[1], d[2]) for _, d in clean).most_common()}")

        for (u, v), d in clean[:args.worm]:
            hf = literal_contract(facets, u, v)          # u survives as the hub
            hd = edge_degrees(hf)
            hadj = defaultdict(set)
            for (a, b) in hd:
                hadj[a].add(b)
                hadj[b].add(a)
            n6w = sum(1 for y in hadj[u] if hd[key(u, y)] >= 6)
            splits = hub_splits(hf, u)
            # is the resulting state the one we came from?
            orig = set(facets)
            novel = 0
            for gamma, nf in splits:
                w2 = max(max(f) for f in nf)
                # which side of the cut kept the label u is arbitrary, so test
                # both relabellings back onto the original {u, v}
                relab = [{w2: v}, {w2: u, u: v}]
                if not any({tuple(sorted(rl.get(x, x) for x in f)) for f in nf}
                           == orig for rl in relab):
                    novel += 1
            print(f"    contract ({u},{v}) r={d[2]} -> hub n6={n6w} (Z{12 + n6w}): "
                  f"{len(splits)} FK-restoring splits, {novel} of them NEW")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--crystal", default=None, choices=sorted(LIBRARY))
    ap.add_argument("--verify", action="store_true",
                    help="check the local degree rules against a literal "
                         "contraction on a sample of sites")
    ap.add_argument("--max-sites", type=int, default=0,
                    help="cap the number of edges censused (0 = all)")
    ap.add_argument("--worm", type=int, default=0, metavar="K",
                    help="on each crystal, take the first K zero-illegal "
                         "contractions and enumerate every FK-restoring split "
                         "of the resulting one-hub state. >1 restoring split "
                         "= a genuine FK-legal -> FK-legal composite move.")
    ap.add_argument("--state", nargs="*", default=None,
                    help="census these .mfd states instead of the library")
    args = ap.parse_args()

    if args.state:
        for p in args.state:
            state_probe(p)
        return

    if args.worm:
        worm_probe(args)
        return

    names = [args.crystal] if args.crystal else sorted(LIBRARY)

    print(f"flat mean edge degree e* = {E_FLAT:.7f}"
          f"   <->  Zbar* = 12/(6 - e*) = {12.0 / (6.0 - E_FLAT):.5f}"
          f"   <->  n6bar* = {12.0 / (6.0 - E_FLAT) - 12.0:.5f}\n")

    hdr = (f"{'crystal':7s} {'f0':>6s} {'Zbar':>8s} {'e_nat':>9s} "
           f"{'e-e*':>9s} {'contractible':>12s} {'min n_ill_e':>12s} "
           f"{'min n_nonfk':>12s} {'best (ill,fk,r)':>16s}")
    print(hdr)
    print("-" * len(hdr))

    for name in names:
        fname, zbar_ideal = LIBRARY[name]
        m = ddg.Manifold.load(os.path.join(REF, fname), 3)
        facets = [tuple(sorted(int(x) for x in f)) for f in m.facets()]
        deg = edge_degrees(facets)
        star = star_map(facets)
        adj = defaultdict(set)
        for (a, b) in deg:
            adj[a].add(b)
            adj[b].add(a)

        f0 = len(adj)
        f1 = len(deg)
        f3 = len(facets)
        zbar = 2.0 * f1 / f0
        e_nat = 6.0 * f3 / f1

        assert all(d in LEGAL for d in deg.values()), f"{name} not edge-legal"

        edges = sorted(deg)
        if args.max_sites:
            edges = edges[:args.max_sites]

        results = []
        for (u, v) in edges:
            d = contraction_damage(facets, star, deg, adj, u, v)
            if d is not None:
                results.append(d)

        if results:
            best = min(results, key=lambda t: (t[0], t[1]))
            min_ill = min(t[0] for t in results)
            min_fk = min(t[1] for t in results)
            bstr = f"({best[0]},{best[1]},r={best[2]})"
        else:
            min_ill = min_fk = -1
            bstr = "none"

        print(f"{name:7s} {f0:6d} {zbar:8.4f} {e_nat:9.6f} "
              f"{e_nat - E_FLAT:+9.6f} {len(results):12d} {min_ill:12d} "
              f"{min_fk:12d} {bstr:>16s}")

        if args.verify and results:
            # literal check on the first few contractible edges
            checked = 0
            for (u, v) in edges:
                d = contraction_damage(facets, star, deg, adj, u, v)
                if d is None:
                    continue
                nf = literal_contract(facets, u, v)
                nd = edge_degrees(nf)
                nadj = defaultdict(set)
                for (a, b) in nd:
                    nadj[a].add(b)
                    nadj[b].add(a)
                ill = sum(1 for x in nd.values() if x not in LEGAL)
                bad = 0
                for x in nadj:
                    inc = [nd[key(x, y)] for y in nadj[x]]
                    n6, mm = n6_and_m(nd, inc)
                    if mm > 0 or n6 not in FK_N6:
                        bad += 1
                status = "OK " if (ill, bad) == (d[0], d[1]) else "MISMATCH"
                print(f"    verify uv=({u},{v}) r={d[2]}: local=({d[0]},{d[1]}) "
                      f"literal=({ill},{bad})  {status}")
                checked += 1
                if checked >= 5:
                    break

    print("\nDamage columns are counts in the CONTRACTED complex; the inputs "
          "are all exactly FK-legal (0, 0).")


if __name__ == "__main__":
    main()
