#!/usr/bin/env python3
"""Two-label dock census: LOCAL collision geometry x LOOP HOLONOMY.

Local label: both chain stacks developed into the frame of the WITNESS
tet -- a pristine tet containing a vertex of each complex's 5-vertex
set, which exists for every dock BY the dock criterion (the one-tet-
layer neighborhood test) -- via shortest canonical legs A-stack ->
witness -> B-stack. This is the relative axis geometry AT the
interaction site: cos^2 exact rational; legs recorded.

Holonomy label: the through-witness route and the direct shortest
A->B route generally differ by a nontrivial loop; the loop element is
extracted exactly as the rigid motion relating B's two developed
placements, and reported as cos(phi) = (tr R - 1)/2 (exact rational).

Outputs: local-class table with dock counts by s0 and P(recombine)
attached (prodA + prodB freed docks), and the loop-holonomy spectrum.
Same machinery applied to the phase-2 collider crossings (m4 host) so
contact channels get local labels too.
"""
import glob
import json
import os
import sys
from collections import Counter, defaultdict
from fractions import Fraction

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
for _p in ("../../python", "../../scripts", "../../tools", "."):
    sys.path.insert(0, os.path.join(_HERE, _p))
import discrete_differential_geometry as ddg
from discrete_differential_geometry.fpkmc import face_apex_maps
from discrete_differential_geometry import development as dev
from worm_helix import bc_orbit
from fp_dock_census_intrinsic import chain_seq, stack_tets, connector

ROOT = os.path.join(_HERE, "..", "..")
DATA = os.path.join(ROOT, "data", "fpkmc")


class Host:
    def __init__(self, ref):
        self.m = ddg.Manifold.load(ref, 3)
        F = np.asarray(self.m.facets())
        _, self.face_of = face_apex_maps(self.m)
        self.dual = {}
        self.v2t = defaultdict(list)
        f2t = defaultdict(list)
        for f4 in F:
            t = frozenset(int(x) for x in f4)
            for v in t:
                f2t[t - {v}].append(t)
                self.v2t[v].append(t)
        for fc, ts in f2t.items():
            for t in ts:
                for u in ts:
                    if u != t:
                        self.dual.setdefault(t, []).append(u)
        self.dual = {t: sorted(ns, key=sorted)
                     for t, ns in self.dual.items()}
        self.F = F

    def complex_verts(self, chord):
        c = tuple(sorted(int(x) for x in chord))
        return set(c) | set(self.face_of[c])

    def witness(self, chordA, chordB):
        """Canonical pristine tet touching both complexes' 5-vert sets."""
        va, vb = self.complex_verts(chordA), self.complex_verts(chordB)
        cands = [t for v in va for t in self.v2t[v] if t & vb]
        if not cands:
            return None
        return min(cands, key=sorted)

    def two_label(self, chordA, chordB):
        cA = tuple(sorted(int(x) for x in chordA))
        cB = tuple(sorted(int(x) for x in chordB))
        seqA = chain_seq(self.m, cA, list(self.face_of[cA]))
        seqB = chain_seq(self.m, cB, list(self.face_of[cB]))
        posA = dev.develop_chain(seqA)
        axA = dev.chain_axis(seqA, posA)
        tetsA, tetsB = stack_tets(seqA), stack_tets(seqB)

        def b_positions(path):
            start = tuple(path[0])
            placements = dev.develop_path(start,
                                          {v: posA[v] for v in start},
                                          [tuple(t) for t in path])
            arrive_tet, arrive_pos = tuple(path[-1]), placements[-1]
            i0 = next(k for k, t in enumerate(tetsB)
                      if t == frozenset(arrive_tet))
            posB = dict(arrive_pos)
            for k in range(i0, len(seqB) - 4):
                posB[seqB[k + 4]] = dev.reflect(
                    posB[seqB[k]], *(posB[v] for v in seqB[k + 1:k + 4]))
            for k in range(i0, 0, -1):
                posB[seqB[k - 1]] = dev.reflect(
                    posB[seqB[k + 3]], *(posB[v] for v in seqB[k:k + 3]))
            return posB

        # direct route
        p_dir = connector(self.dual, tetsA, tetsB)
        posB_dir = b_positions(p_dir)
        # witness route
        w = self.witness(cA, cB)
        if w is None:
            return None
        leg1 = connector(self.dual, tetsA, {w})
        leg2 = connector(self.dual, {w}, tetsB)
        p_wit = leg1 + leg2[1:]
        posB_wit = b_positions(p_wit)

        axB_wit = dev.chain_axis(seqB, posB_wit)
        local = dev.cos2_between(axA, axB_wit)
        axB_dir = dev.chain_axis(seqB, posB_dir)
        direct = dev.cos2_between(axA, axB_dir)

        # loop holonomy: rigid motion between B's two placements
        t0 = seqB[0:4]
        P = [posB_dir[v] for v in t0]
        Q = [posB_wit[v] for v in t0]
        M = [dev._sub(P[i + 1], P[0]) for i in range(3)]
        N = [dev._sub(Q[i + 1], Q[0]) for i in range(3)]
        MT = [tuple(M[i][j] for i in range(3)) for j in range(3)]
        tr = Fraction(0)
        for j in range(3):
            rhs = tuple(N[i][j] for i in range(3))
            tr += dev._solve3(MT, rhs)[j]
        cos_phi = (tr - 1) / 2
        return {"local": local, "direct": direct, "cos_phi": cos_phi,
                "legs": (len(leg1) - 1, len(leg2) - 1),
                "path_dir": len(p_dir)}


def ang(c2):
    return float(np.degrees(np.arccos(np.sqrt(float(c2)))))


def main():
    m2 = Host(os.path.join(ROOT, "data/tcp_reference/T3_R_m2_N7248.mfd"))
    F = m2.F
    orb = [int(x) for x in bc_orbit(m2.m, [int(x) for x in F[0]])]
    L = len(orb)

    # ---- FP docks: census (prodA) + freed recombination docks (prodB) ----
    items = []
    for pat, kind in (("prodA_s*_p*.json", "census"),
                      ("prodB_p*.json", "recombine"),
                      ("prodB2_p*.json", "recombine")):
        for f in sorted(glob.glob(os.path.join(DATA, pat))):
            d = json.load(open(f))
            for e in d["episodes"]:
                if not e.get("dock_chord"):
                    continue
                if kind == "recombine" and not (e.get("freed")
                                                and "outcome" in e):
                    continue
                items.append({"kind": kind, "sep": d["sep_sites"],
                              "chordA": e["dock_chord"],
                              "chordB": (orb[e["jB"] % L],
                                         orb[(e["jB"] + 4) % L]),
                              "outcome": e.get("outcome")})
    out_rows = []
    fails = 0
    for it in items:
        try:
            lab = m2.two_label(it["chordA"], it["chordB"])
        except Exception as ex:
            fails += 1
            continue
        if lab is None:
            fails += 1
            continue
        out_rows.append({**it, **lab})
    print(f"{len(out_rows)}/{len(items)} FP docks labeled ({fails} fails)")

    # ---- local-class table ----
    by_local = defaultdict(lambda: {"n": 0, "rec": 0, "esc": 0,
                                    "legs": Counter()})
    for r in out_rows:
        g = by_local[r["local"]]
        g["n"] += 1
        g["legs"][r["legs"]] += 1
        if r["outcome"] == "recombine":
            g["rec"] += 1
        elif r["outcome"] == "escape":
            g["esc"] += 1
    print(f"\n{len(by_local)} LOCAL classes (witness frame):")
    print(f"{'local cos^2':>16} {'angle':>8} {'docks':>6} "
          f"{'P(rec)':>7} {'common legs':>14}")
    for c2, g in sorted(by_local.items(), key=lambda kv: -kv[1]["n"]):
        tot = g["rec"] + g["esc"]
        p = f"{g['rec'] / tot:.2f}({tot})" if tot else "--"
        legs = ",".join(f"{a}+{b}" for (a, b), _ in
                        g["legs"].most_common(2))
        print(f"{str(c2):>16} {ang(c2):7.2f}d {g['n']:6d} {p:>7} "
              f"{legs:>14}")

    # ---- holonomy spectrum ----
    hol = Counter(r["cos_phi"] for r in out_rows)
    print(f"\nloop holonomy (witness vs direct route), "
          f"{len(hol)} classes:")
    for cp, n in hol.most_common(12):
        phi = float(np.degrees(np.arccos(max(-1.0, min(1.0, float(cp))))))
        print(f"  cos(phi) = {str(cp):>12} ({phi:6.2f} deg): {n}")

    # ---- collider crossings (m4) ----
    print("\ncollider crossings (m4 host), local labels:")
    m4 = Host(os.path.join(ROOT,
                           "data/tcp_reference/T3_R_m4_N57984.mfd"))
    for name in ("crossing_default", "crossing_steep"):
        p = os.path.join(DATA, name + ".json")
        if not os.path.exists(p):
            continue
        d = json.load(open(p))
        for r in d["results"]:
            try:
                lab = m4.two_label(r["chordA"], r["chordB"])
            except Exception as ex:
                print(f"  {name}#{r['crossing']}: fail {ex}")
                continue
            if lab is None:
                print(f"  {name}#{r['crossing']}: no witness "
                      f"(dmin {r['dmin']:.2f}) -- direct-route label only")
                continue
            o = r["outcome"]
            v = f"V={o['V']:+.3f}" if o.get("V") is not None else ""
            fstr = f" f={tuple(o['species_f'])}" if o.get("species_f") \
                else ""
            print(f"  {name}#{r['crossing']}: local {lab['local']} "
                  f"({ang(lab['local']):.2f}d) legs {lab['legs']} "
                  f"[direct {lab['direct']}] {o['type']} {v}{fstr}")

    with open(os.path.join(DATA, "fp_dock_local_census.json"), "w") as fh:
        json.dump([{**r, "local": str(r["local"]),
                    "direct": str(r["direct"]),
                    "cos_phi": str(r["cos_phi"]),
                    "chordA": [int(x) for x in r["chordA"]],
                    "chordB": [int(x) for x in r["chordB"]]}
                   for r in out_rows], fh)
    print(f"\nwrote {os.path.join(DATA, 'fp_dock_local_census.json')}")


if __name__ == "__main__":
    main()
