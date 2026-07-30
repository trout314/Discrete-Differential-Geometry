"""Template census at all deg-4 anchors of lam35r_snap15000.

For every anchor:
  1. enumerate the worm class candidates (oracle);
  2. per candidate, build the MOVE-CLASS key: canonical form of
     (pre-state support patch, edge-degree decoration, anchor mark, move
     pair), minimized jointly so key equality == decorated isomorphism
     carrying the move list;
  3. build the LINK-PAIR key: canonical form of star(a) u star(b) with
     edge-degree decoration + anchor mark (the "pair of vertex links").

Then:
  A. count distinct move classes (= templates-with-context), coverage of
     the top-N (raw and e^-dS weighted);
  B. template-property test: does the link-pair class determine the
     candidate-set fingerprint (multiset of move-class keys)?
  C. fraction of candidates whose second move reaches OUTSIDE the
     link-pair patch (why links alone cannot be the template domain).
"""
import os, sys, json, math, time
from collections import Counter, defaultdict
from itertools import combinations, permutations, product

_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, os.path.join(_ROOT, "python"))
sys.path.insert(0, os.path.join(_ROOT, "scripts/defect_dynamics"))
import numpy as np
import discrete_differential_geometry as ddg
import worm_deg4_slide as W
from defect_state import canonical_key

SNAP = os.environ.get("SNAP",
                      os.path.join(_ROOT, "data/mgas/lam35r_snap15000.mfd"))
OUT = os.environ.get("OUT",
                     os.path.join(_ROOT, "data/template_census_results.json"))
PERM_CAP = 200000
ANCHOR_MARK = 1000                     # added to the anchor edge's ecolor


def mv_verts(mv):
    if mv[0] == "23":
        return set(mv[1]) | {mv[2], mv[3]}
    return set(mv[1]) | {mv[2], mv[3], mv[4]}


def relabel_mv(mv, idx):
    if mv[0] == "23":
        return ("23", tuple(sorted(idx[v] for v in mv[1])),
                tuple(sorted((idx[mv[2]], idx[mv[3]]))))
    return ("32", tuple(sorted(idx[v] for v in mv[1])),
            tuple(sorted((idx[mv[2]], idx[mv[3]], idx[mv[4]]))))


def refine(verts, ecolor):
    """1-WL colour refinement on the complete decorated graph over verts."""
    color = {v: 0 for v in verts}
    for _ in range(len(verts)):
        sig = {}
        for v in verts:
            ring = tuple(sorted(
                (ecolor.get(tuple(sorted((v, u))), 0), color[u])
                for u in verts if u != v))
            sig[v] = (color[v], ring)
        vals = {s: i for i, s in enumerate(sorted(set(sig.values())))}
        new = {v: vals[sig[v]] for v in verts}
        if new == color:
            break
        color = new
    return color


def canon_joint(verts, facets, ecolor, moves):
    """Canonical form of (decorated patch, move pair).  Returns (key, exact).

    Minimizes (relabeled facets, relabeled full ecolor, relabeled moves)
    over all vertex orders consistent with colour refinement."""
    verts = sorted(verts)
    color = refine(verts, ecolor)
    classes = defaultdict(list)
    for v in verts:
        classes[color[v]].append(v)
    groups = [sorted(classes[c]) for c in sorted(classes)]
    total = 1
    for g in groups:
        total *= math.factorial(len(g))
    exact = total <= PERM_CAP

    pairs = [tuple(sorted(p)) for p in combinations(verts, 2)]
    best = None
    seen = 0
    for combo in product(*(permutations(g) for g in groups)):
        order = [v for g in combo for v in g]
        idx = {v: i for i, v in enumerate(order)}
        rel_f = tuple(sorted(tuple(sorted(idx[v] for v in f))
                             for f in facets))
        rel_e = tuple(sorted(
            (tuple(sorted((idx[u], idx[w]))), ecolor.get((u, w), 0))
            for (u, w) in pairs))
        rel_m = tuple(sorted(relabel_mv(mv, idx) for mv in moves))
        rel = (rel_f, rel_e, rel_m)
        if best is None or rel < best:
            best = rel
        seen += 1
        if seen >= PERM_CAP:
            break
    return best, exact


def main():
    m = ddg.Manifold.load(SNAP, 3)
    L = W.Live(m)
    anchors = sorted(L.deg4())
    print(f"{os.path.basename(SNAP)}: {len(anchors)} deg-4 anchors")

    move_classes = Counter()           # key-hash -> count
    class_info = {}                    # key-hash -> dict
    anchor_rows = []
    nbhd_groups = defaultdict(list)    # nbhd key-hash -> [anchor index]
    inexact_joint = inexact_nbhd = 0
    n_cands_total = 0
    n_m2_outside = 0
    boltz_total = 0.0
    t0 = time.time()

    for ai, anchor in enumerate(anchors):
        a, b = anchor
        # ---- link-pair patch: star(a) u star(b), decorated ----
        star_tets = sorted(L.v2t[a] | L.v2t[b])
        star_verts = set().union(*star_tets)
        ec_star = {}
        for t in star_tets:
            for e in combinations(sorted(t), 2):
                ec_star[e] = L.edeg[e]
        ec_star[anchor] = L.edeg[anchor] + ANCHOR_MARK
        nkey, nexact = canonical_key(star_tets, ecolor=ec_star)
        if not nexact:
            inexact_nbhd += 1
        nh = hash(nkey)

        cands = W.candidates(L, anchor)
        n_cands_total += len(cands)
        fp = []
        for (m1, m2, ds, land, netk, gone4, new4) in cands:
            V = mv_verts(m1) | mv_verts(m2) | set(anchor)
            facets = set()
            for v in V:
                for t in L.v2t[v]:
                    if set(t) <= V:
                        facets.add(t)
            ec = {}
            for e in (tuple(sorted(p)) for p in combinations(sorted(V), 2)):
                d = L.edeg.get(e, 0)
                if d:
                    ec[e] = d
            ec[anchor] = L.edeg[anchor] + ANCHOR_MARK
            key, exact = canon_joint(V, facets, ec, (m1, m2))
            if not exact:
                inexact_joint += 1
            kh = hash(key)
            move_classes[kh] += 1
            k = len(gone4)
            info = class_info.setdefault(kh, {
                "kinds": m1[0] + "+" + m2[0], "k": k, "dS": set(),
                "anchors": set(), "n_verts": len(V),
                "n_facets": len(facets)})
            info["dS"].add(round(ds, 6))
            info["anchors"].add(ai)
            boltz_total += math.exp(-min(ds, 50.0))
            if not mv_verts(m2) <= star_verts:
                n_m2_outside += 1
            fp.append(kh)

        fp_key = tuple(sorted(Counter(fp).items()))
        nbhd_groups[nh].append(ai)
        anchor_rows.append({
            "anchor": list(anchor), "n_cands": len(cands),
            "nbhd_hash": nh, "fingerprint": fp_key,
            "star_nv": len(star_verts)})
        if (ai + 1) % 10 == 0 or ai == len(anchors) - 1:
            print(f"  [{ai+1}/{len(anchors)}] cands so far {n_cands_total}, "
                  f"move classes {len(move_classes)}, "
                  f"{time.time()-t0:.0f}s", flush=True)

    # ---- report A: move classes ----
    print(f"\n=== A. MOVE CLASSES (templates-with-context) ===")
    print(f"total candidates: {n_cands_total} at "
          f"{sum(1 for r in anchor_rows if r['n_cands'])} live anchors "
          f"({len(anchors)} anchors)")
    print(f"distinct move classes: {len(move_classes)} "
          f"(joint canon inexact: {inexact_joint}, "
          f"nbhd canon inexact: {inexact_nbhd})")
    ranked = move_classes.most_common()
    cum = cum_b = 0.0
    print(f"{'rank':>4} {'count':>6} {'cum%':>6} {'kinds':>6} {'k':>2} "
          f"{'#anch':>5} {'nV':>3} {'nF':>3}  dS values")
    for r, (kh, c) in enumerate(ranked[:25]):
        cum += c
        info = class_info[kh]
        ds_str = ",".join(f"{x:+.3f}" for x in sorted(info["dS"]))
        print(f"{r+1:>4} {c:>6} {100*cum/max(n_cands_total,1):>5.1f}% "
              f"{info['kinds']:>6} {info['k']:>2} "
              f"{len(info['anchors']):>5} {info['n_verts']:>3} "
              f"{info['n_facets']:>3}  {ds_str}")
    if len(ranked) > 25:
        print(f"  ... and {len(ranked)-25} more classes "
              f"(tail sizes: {[c for _, c in ranked[25:35]]}...)")
    for topn in (1, 5, 10, 20, len(ranked)):
        got = sum(c for _, c in ranked[:topn])
        gb = sum(c * math.exp(-min(min(class_info[kh]['dS']), 50.0))
                 for kh, c in ranked[:topn])
        print(f"  top-{topn:<3} coverage: {100*got/max(n_cands_total,1):5.1f}% raw, "
              f"{100*gb/max(boltz_total,1e-300):5.1f}% e^-dS-weighted")

    # ---- report B: template property of the link pair ----
    print(f"\n=== B. LINK-PAIR -> MOVE-LIST test ===")
    print(f"distinct link-pair classes among {len(anchors)} anchors: "
          f"{len(nbhd_groups)}")
    multi = {h: idxs for h, idxs in nbhd_groups.items() if len(idxs) > 1}
    n_ok = n_bad = 0
    for h, idxs in sorted(multi.items(), key=lambda kv: -len(kv[1])):
        fps = {anchor_rows[i]["fingerprint"] for i in idxs}
        ok = len(fps) == 1
        n_ok += ok
        n_bad += (not ok)
        if not ok or len(idxs) > 2:
            counts = [anchor_rows[i]["n_cands"] for i in idxs]
            print(f"  class ({len(idxs)} anchors) cand-counts {counts} "
                  f"-> {'SAME move list' if ok else 'DIFFERENT move lists'}")
    print(f"multi-anchor link-pair classes: {len(multi)}  "
          f"(determine moves: {n_ok}, fail to: {n_bad})")
    singles = len(nbhd_groups) - len(multi)
    print(f"singleton link-pair classes: {singles} (no test possible)")
    print(f"candidates whose 2nd move leaves the link-pair patch: "
          f"{n_m2_outside}/{n_cands_total}")

    # ---- persist ----
    out = {
        "snap": SNAP, "n_anchors": len(anchors),
        "n_cands": n_cands_total,
        "n_move_classes": len(move_classes),
        "inexact_joint": inexact_joint, "inexact_nbhd": inexact_nbhd,
        "n_m2_outside_linkpair": n_m2_outside,
        "classes": [
            {"count": c, **{k2: (sorted(v) if isinstance(v, set) else v)
                            for k2, v in class_info[kh].items()}}
            for kh, c in ranked],
        "anchors": anchor_rows,
        "n_linkpair_classes": len(nbhd_groups),
        "linkpair_multi_ok": n_ok, "linkpair_multi_bad": n_bad,
    }
    with open(OUT, "w") as f:
        json.dump(out, f, default=list)
    print(f"\nresults saved: {OUT}")


if __name__ == "__main__":
    main()
