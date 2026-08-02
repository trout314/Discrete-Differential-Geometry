#!/usr/bin/env python3
"""Intrinsic (PL) dock-angle census -- the development-module test case.

Recomputes the FP production dock census (prodA_*.json) with INTRINSIC
contact angles per CONVENTIONS.md section 6: both chains' tet stacks are
rigidly developed into R^3 in a COMMON frame -- A's stack first, then a
shortest dual-graph tet path (canonical: BFS with sorted neighbors) from
A's stack to B's, then B's stack -- and the angle between the two exact
helix axes is an exact rational cos^2(theta), a PL invariant of the
local complex + connecting path. Exactness is enforced throughout
(congruence of every placed tet; step trace == -1/3; axis constancy
along each stack).

PREDICTION under test: the intrinsic census should QUANTIZE into
discrete cos^2 classes where the registry census (fp_dock_angles.py) was
a sin-theta-like quasi-continuum.

Sanity anchors run first: (i) two sites on the SAME chain develop to
cos^2 == 1 exactly; (ii) the head-on V5 geometry.
"""
import glob
import json
import os
import sys
from collections import Counter, deque
from fractions import Fraction

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

_HERE = os.path.dirname(os.path.abspath(__file__))
for _p in ("../../python", "../../scripts", "../../tools", "."):
    sys.path.insert(0, os.path.join(_HERE, _p))
import discrete_differential_geometry as ddg
from discrete_differential_geometry.fpkmc import face_apex_maps
from discrete_differential_geometry import development as dev
from chain_select import chain_for_run, add_chain_args
from worm_helix import bc_orbit

ROOT = os.path.join(_HERE, "..", "..")
DATA = os.path.join(ROOT, "data", "fpkmc")
REF = os.path.join(ROOT, "data", "tcp_reference", "T3_R_m2_N7248.mfd")
HALF = 8          # stack half-length in windows (2 precession periods)


def chain_seq(m, chord, face, half=HALF):
    """Self-consistent chain vertex sequence of 2*half windows through
    the chord's site: walk backward from the site, reverse the far
    window, walk forward through."""
    bwd = [int(x) for x in m.chain_walk([chord[1]] + list(face), half)]
    far = list(reversed(bwd[-4:]))
    seq = [int(x) for x in m.chain_walk(far, 2 * half + 1)]
    return seq


def stack_tets(seq):
    return [frozenset(seq[k:k + 4]) for k in range(len(seq) - 3)]


def connector(dual, a_tets, b_tets, cap=40):
    """Canonical shortest dual path from A's stack to B's stack:
    multi-source BFS, neighbors visited in sorted order."""
    b_set = set(b_tets)
    q = deque()
    prev = {}
    for t in a_tets:
        prev[t] = None
        q.append((t, 0))
    while q:
        t, d = q.popleft()
        if t in b_set:
            path = []
            while t is not None:
                path.append(t)
                t = prev[t]
            return list(reversed(path))
        if d >= cap:
            continue
        for nb in dual[t]:
            if nb not in prev:
                prev[nb] = t
                q.append((nb, d + 1))
    return None


def intrinsic_cos2(m, dual, chordA, faceA, chordB, faceB):
    seqA = chain_seq(m, chordA, faceA)
    seqB = chain_seq(m, chordB, faceB)
    posA = dev.develop_chain(seqA)
    axA = dev.chain_axis(seqA, posA)

    tetsA, tetsB = stack_tets(seqA), stack_tets(seqB)
    path = connector(dual, tetsA, tetsB)
    if path is None:
        raise RuntimeError("no dual path within cap")
    start = tuple(path[0])
    start_pos = {v: posA[v] for v in start}
    placements = dev.develop_path(start, start_pos, [tuple(t) for t in path])
    arrive_tet, arrive_pos = tuple(path[-1]), placements[-1]

    i0 = next(k for k, t in enumerate(tetsB) if t == frozenset(arrive_tet))
    posB = dict(arrive_pos)
    for k in range(i0, len(seqB) - 4):          # forward along the stack
        posB[seqB[k + 4]] = dev.reflect(
            posB[seqB[k]], *(posB[v] for v in seqB[k + 1:k + 4]))
    for k in range(i0, 0, -1):                  # backward
        posB[seqB[k - 1]] = dev.reflect(
            posB[seqB[k + 3]], *(posB[v] for v in seqB[k:k + 3]))
    axB = dev.chain_axis(seqB, posB)
    return dev.cos2_between(axA, axB), len(path)


def main():
    m = ddg.Manifold.load(REF, 3)
    F = np.asarray(m.facets())
    _cc, _kcls, _seq, chain_prov = chain_for_run(
        REF, F, None, seed_tet=0)
    orb = [int(x) for x in _seq]
    L = len(orb)
    _, face_of = face_apex_maps(m)
    dual = {}
    face2tets = {}
    for f4 in F:
        t = frozenset(int(x) for x in f4)
        for v in t:
            face2tets.setdefault(t - {v}, []).append(t)
    for fc, ts in face2tets.items():
        for t in ts:
            for u in ts:
                if u != t:
                    dual.setdefault(t, []).append(u)
    dual = {t: sorted(ns, key=sorted) for t, ns in dual.items()}

    def chord_face(chord):
        c = tuple(sorted(int(x) for x in chord))
        return c, list(face_of[c])

    # ---- sanity anchors ----
    cA, fA = chord_face((orb[40], orb[44]))
    cB, fB = chord_face((orb[64], orb[68]))
    c2, plen = intrinsic_cos2(m, dual, cA, fA, cB, fB)
    print(f"sanity same-chain (sites 10 apart): cos^2 = {c2} "
          f"(path {plen}) -> {'OK' if c2 == 1 else 'FAIL'}")

    # ---- production docks ----
    rows = []
    fails = 0
    for f in sorted(glob.glob(os.path.join(DATA, "prodA_s*_p*.json"))):
        d = json.load(open(f))
        for e in d["episodes"]:
            if not e.get("dock_chord"):
                continue
            try:
                cA, fA = chord_face(e["dock_chord"])
                cB, fB = chord_face((orb[e["jB"] % L],
                                     orb[(e["jB"] + 4) % L]))
                c2, plen = intrinsic_cos2(m, dual, cA, fA, cB, fB)
                rows.append({"cos2": c2, "path": plen,
                             "sep": d["sep_sites"],
                             "V": e.get("V_dock")})
            except Exception as ex:
                fails += 1
                if fails <= 5:
                    print(f"  fail: {type(ex).__name__}: {ex}")
    print(f"\n{len(rows)} docks developed, {fails} failures")

    # ---- exact class table ----
    classes = Counter(r["cos2"] for r in rows)
    print(f"\n{len(classes)} exact cos^2 classes:")
    print(f"{'cos^2 (exact)':>28} {'angle':>8} {'count':>6}")
    for c2, n in sorted(classes.items(), key=lambda kv: -kv[1])[:25]:
        ang = float(np.degrees(np.arccos(np.sqrt(float(c2)))))
        print(f"{str(c2):>28} {ang:7.2f}d {n:6d}")
    n_top = sum(n for _, n in classes.most_common(10))
    print(f"top 10 classes hold {n_top}/{len(rows)} docks")

    # ---- figure: intrinsic spectrum vs registry smear ----
    angs = [float(np.degrees(np.arccos(np.sqrt(float(r['cos2'])))))
            for r in rows]
    fig, ax = plt.subplots(figsize=(9, 4.5))
    ax.hist(angs, bins=np.arange(-0.25, 90.5, 0.5), color="tab:blue")
    ax.set_xlabel("intrinsic contact angle (deg), exact classes")
    ax.set_ylabel("dock count")
    ax.set_title(f"INTRINSIC dock-angle census ({len(rows)} docks, "
                 f"{len(classes)} exact cos^2 classes) -- "
                 "R m2 N7248, lam=0.40, development frame, canonical "
                 "shortest connecting path", fontsize=9)
    fig.tight_layout()
    out = os.path.join(DATA, "fp_dock_census_intrinsic.png")
    fig.savefig(out, dpi=140)
    print(f"\nSaved to: {os.path.abspath(out)}")
    with open(os.path.join(DATA, "fp_dock_census_intrinsic.json"),
              "w") as fh:
        json.dump([{"cos2": str(r["cos2"]), "path": r["path"],
                    "sep": r["sep"], "V": r["V"]} for r in rows], fh)


if __name__ == "__main__":
    main()
