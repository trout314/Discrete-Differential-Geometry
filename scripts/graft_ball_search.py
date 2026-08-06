#!/usr/bin/env python3
"""Search for vertex-count-changing BALL grafts (dV != 0 FK->FK block moves).

A ball graft cuts out a lump ~ B^3 (here: the union of closed stars of a small
vertex set, e.g. the 3 vertices of a face) and glues back a DIFFERENT FK-legal
filling of the SAME framed boundary -- boundary 2-sphere + per-boundary-edge
interior degree, i.e. the L2/L3 decoration of graft_signature.py restricted to
S^2.  A filling with fewer interior vertices is a defect-free local vertex
deletion.

Method: best-first enumeration in the spirit of defect_dynamics/
enumerate_fillings.py, with two extensions:
  * the interior move set gains 1->4 / 4->1 (vertex insertion/removal) --
    the vertex-preserving flips (2->3, 3->2, 4-4) conserve V_int, so the
    original searcher could never see a dV != 0 refilling;
  * exact framing accounting: interior 2->3 flips CAN change a boundary
    edge's interior degree (when a face edge lies on the boundary), so
    framing deviation is tracked move-by-move and acceptance requires exact
    restoration (frame_dev == 0), with deviation allowed mid-path.
Children are scored by exact delta bookkeeping at push time and analyzed only
when popped, so the state cap bounds the expensive work.  Priority =
(illegality + framing deviation, interior vertices, tets).

For every accepted filling with V_int != seed V_int the graft is performed on
the host crystal and fully validated (all {5,6}, no broken disclination
lines, chi = 0, orientable).  If none is found, the minimum
(illegality + framing deviation) reached at each smaller V_int quantifies how
non-FK such a graft would have to be ("at best").

Usage:
    caffeinate -i python scripts/graft_ball_search.py \
        [--mfd scripts/defect_dynamics/T3_C15_m3_N3672.mfd] \
        [--seeds vertex,edge,face] [--cap 50000] [--max-bad 6] [--grow 12]
"""
import argparse
import heapq
import os
import sys
import time
from collections import Counter
from itertools import combinations

import numpy as np

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_ROOT, "python"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_ROOT, "scripts", "defect_dynamics"))

import fk_moves as fk
from enumerate_fillings import neighbors as flip_neighbors
from graft_signature import validate_facets
from discrete_differential_geometry import Manifold
from discrete_differential_geometry.symmetry import CrystalSymmetry


def vertex_moves(tets, bverts, spare):
    """Interior 1->4 (with fresh label `spare`) and 4->1 moves.

    Returns [(rem, add, dV)].  4->1: an interior vertex in exactly 4 tets has
    link = boundary of a tetrahedron; replace its star by the outer tet.
    """
    out = []
    vt = {}
    tetset = set(tets)
    for t in tets:
        for v in t:
            vt.setdefault(v, []).append(t)
    for v, ts in vt.items():
        if v in bverts or len(ts) != 4:
            continue
        outer = set().union(*ts) - {v}
        if len(outer) != 4:
            continue
        newt = tuple(sorted(outer))
        if newt in tetset:
            continue
        out.append((frozenset(ts), frozenset([newt]), -1))
    if spare is not None:
        for t in tets:
            add = frozenset(tuple(sorted((*f, spare)))
                            for f in combinations(t, 3))
            out.append((frozenset([t]), add, +1))
    return out


def child_cost(a, rem, add, frame, bedges, n_ill, frame_dev):
    """Exact (illegality, framing deviation) of the child, from edge deltas."""
    delta = Counter()
    for t in rem:
        for e in combinations(t, 2):
            delta[e] -= 1
    for t in add:
        for e in combinations(t, 2):
            delta[e] += 1
    ill, fdev = n_ill, frame_dev
    for e, d in delta.items():
        if d == 0:
            continue
        old = a["edeg"].get(e, 0)
        new = old + d
        if e in bedges:
            fdev += abs(new - frame[e]) - abs(old - frame[e])
        else:
            ill += (new > 0 and new not in (5, 6)) - \
                   (old > 0 and old not in (5, 6))
    return ill, fdev


def violation_profile(a, frame, bedges):
    """Human-readable list of the state's violations."""
    out = []
    for e in a["iedges"]:
        if a["edeg"][e] not in (5, 6):
            out.append(f"iedge{e}=d{a['edeg'][e]}")
    for e in bedges:
        d = a["edeg"].get(e, 0)
        if d != frame[e]:
            out.append(f"bedge{e}=d{d}(want{frame[e]})")
    return out


def search_seed(host_F, seed_verts, label, cap, max_bad, grow, spares,
                out_dir, host_name, vfirst=False):
    """Best-first refilling search from the closed-star lump of seed_verts."""
    t0 = time.time()
    seed_set = set(seed_verts)
    lump_idx = np.nonzero(np.isin(host_F, list(seed_set)).any(axis=1))[0]
    seed = frozenset(tuple(sorted(int(x) for x in t)) for t in host_F[lump_idx])
    a0 = fk.analyze(list(seed))
    bfaces0, bedges = a0["bfaces"], a0["bedges"]
    bverts = a0["bverts"]
    frame = {e: a0["edeg"][e] for e in bedges}
    V0 = len(a0["iverts"])
    nb_v = len(bverts)
    nb_e = len(bedges)
    chi = nb_v - nb_e + len(bfaces0)
    print(f"\n[{label}] lump {len(seed)} tets, boundary {len(bfaces0)} faces "
          f"(chi={chi}), V_int={V0} {sorted(a0['iverts'])}")
    if chi != 2:
        print(f"  [{label}] boundary is not S^2 -- skipped")
        return {}
    if set(seed_set) != a0["iverts"]:
        print(f"  [{label}] note: interior = {sorted(a0['iverts'])} "
              f"differs from seed set")

    spare_base = int(host_F.max()) + 1
    spare_labels = list(range(spare_base, spare_base + spares))
    max_tets = len(seed) + grow

    popped = [seed]                        # index -> state
    seen = {hash(seed)}
    heap = []
    tie = 0
    # root has exact stats from a0
    root_ill = sum(1 for e in a0["iedges"] if a0["edeg"][e] not in (5, 6))
    stack0 = (root_ill, 0)
    heap.append(((root_ill, V0, len(seed)), tie, 0, None, None,
                 root_ill, 0, V0))
    heapq.heapify(heap)

    found = {}                             # filling_canon -> (V_int, tets)
    best_bad = {}                          # V_int -> min(ill+fdev) seen
    best_prof = {}                         # V_int -> violation profile
    hits = []
    pops = 0
    while heap and pops < cap:
        (_, _, pidx, rem, add, ill, fdev, vint) = heapq.heappop(heap)
        cur = popped[pidx] if rem is None else \
            frozenset((popped[pidx] - rem) | add)
        pops += 1
        a = fk.analyze(list(cur))
        fc = Counter()
        for t in a["tets"]:
            for f in combinations(t, 3):
                fc[f] += 1
        if (a["bfaces"] != bfaces0 or frozenset(a["tets"]) != cur
                or max(fc.values()) > 2):
            continue          # move damaged the boundary / dup tet / pinched
        my_idx = len(popped)
        popped.append(cur)
        bad = ill + fdev
        if vint != V0 and bad < best_bad.get(vint, 99):
            best_bad[vint] = bad
            best_prof[vint] = violation_profile(a, frame, bedges)
        if bad == 0 and all(fk._vfk(a, v) for v in a["iverts"]):
            key = fk.filling_canon(a)
            if key not in found:
                found[key] = (vint, sorted(cur))
                if vint != V0:
                    hits.append((vint, sorted(cur)))
                    print(f"  [{label}] *** HIT: FK filling with V_int={vint} "
                          f"(seed {V0}), {len(cur)} tets, pop {pops}")

        # spare label: lowest not used in this state
        used = {v for t in cur for v in t}
        spare = next((s for s in spare_labels if s not in used), None)
        moves = [(frozenset(r), frozenset(ad), 0)
                 for r, ad in flip_neighbors(cur, bfaces0)]
        moves += vertex_moves(cur, bverts, spare)
        for rem2, add2, dV in moves:
            child = frozenset((cur - rem2) | add2)
            if len(child) > max_tets:
                continue
            h = hash(child)
            if h in seen:
                continue
            c_ill, c_fdev = child_cost(a, rem2, add2, frame, bedges, ill, fdev)
            if c_ill + c_fdev > max_bad:
                continue
            seen.add(h)
            tie += 1
            prio = ((vint + dV, c_ill + c_fdev, len(child)) if vfirst
                    else (c_ill + c_fdev, vint + dV, len(child)))
            heapq.heappush(heap, (prio, tie, my_idx, rem2, add2,
                                  c_ill, c_fdev, vint + dV))

    n_by_v = Counter(v for v, _ in found.values())
    print(f"  [{label}] {pops} pops, {len(seen)} states, "
          f"{len(found)} FK fillings by V_int: {dict(n_by_v)}; "
          f"best (ill+fdev) at other V_int: "
          f"{ {v: b for v, b in sorted(best_bad.items())} } "
          f"({time.time()-t0:.0f}s)")
    for v in sorted(best_prof):
        print(f"    best V_int={v} profile: {best_prof[v]}")

    results = dict(label=label, seed_Vint=V0, n_tets=len(seed),
                   pops=pops, fillings={int(v): int(n) for v, n
                                        in n_by_v.items()},
                   best_bad={int(v): int(b) for v, b in best_bad.items()})
    # perform + validate any dV != 0 graft
    for vint, tets in hits:
        keep = np.ones(len(host_F), bool)
        keep[lump_idx] = False
        newF = np.vstack([host_F[keep],
                          np.array(tets, np.int64)])
        lab, inv = np.unique(newF, return_inverse=True)
        newF = inv.reshape(newF.shape).astype(np.int64)
        try:
            rep = validate_facets(newF)
        except AssertionError as exc:
            print(f"  [{label}] graft dV={vint - V0:+d}: INVALID ({exc})")
            continue
        good = (rep["all_56"] and rep["n_broken_disclination"] == 0
                and rep["euler_characteristic"] == 0 and rep["orientable"])
        print(f"  [{label}] graft dV={vint - V0:+d}: all56={rep['all_56']} "
              f"broken={rep['n_broken_disclination']} "
              f"chi={rep['euler_characteristic']} orient={rep['orientable']} "
              f"V={rep['n_vertices']} Z={rep['z_census']}")
        if good:
            path = os.path.join(out_dir,
                                f"T3_{host_name}_ballgraft_{label}"
                                f"_dV{vint - V0:+d}_f3{len(newF)}.mfd")
            Manifold(3, newF.tolist()).save(path, comments=[
                f"ball graft at seed {label} ({sorted(seed_set)}): "
                f"V_int {V0} -> {vint}",
                f"validation: {rep}"])
            print(f"  [{label}] SAVED -> {path}")
            results.setdefault("grafts", []).append(path)
    return results


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--mfd", default=os.path.join(
        _ROOT, "data", "tcp_reference", "T3_C15_m3_N3672.mfd"))
    ap.add_argument("--seeds", default="vertex,edge,face")
    ap.add_argument("--cap", type=int, default=50000)
    ap.add_argument("--max-bad", type=int, default=6)
    ap.add_argument("--grow", type=int, default=12)
    ap.add_argument("--spares", type=int, default=1)
    ap.add_argument("--vfirst", action="store_true",
                    help="goal-directed: prioritize low V_int over low "
                         "violation (tunnels toward dV<0 fillings)")
    ap.add_argument("--only", default=None,
                    help="comma list of seed labels to run, e.g. face2,edge0")
    ap.add_argument("--out", default=os.path.join(_ROOT, "data", "grafts"))
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    host_name = os.path.basename(args.mfd).split(".")[0]
    mfd = Manifold.load(args.mfd, 3)
    F = np.asarray(mfd.facets(), np.int64)
    print(f"host {host_name}: V={len(np.unique(F))} tets={len(F)}")
    sym = CrystalSymmetry.for_manifold_path(args.mfd)
    print(f"|Aut| = {sym.order}; orbits: " +
          ", ".join(f"{k} {sym.n_orbits(k)}"
                    for k in ("vertex", "edge", "face")))

    all_res = []
    for kind in args.seeds.split(","):
        reps = sym.orbit_representatives(kind)
        only = set(args.only.split(",")) if args.only else None
        for i, rep in enumerate(reps):
            if only and f"{kind}{i}" not in only:
                continue
            verts = [rep] if kind == "vertex" else list(rep)
            res = search_seed(F, verts, f"{kind}{i}", args.cap, args.max_bad,
                              args.grow, args.spares, args.out, host_name,
                              vfirst=args.vfirst)
            all_res.append(res)

    n_hits = sum(len(r.get("grafts", [])) for r in all_res if r)
    print(f"\n=== total dV!=0 validated grafts: {n_hits} ===")


if __name__ == "__main__":
    main()
