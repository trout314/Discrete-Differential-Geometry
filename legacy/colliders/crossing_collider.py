#!/usr/bin/env python3
"""Phase-2 collider: knots on CROSSING BC chains -- angle-resolved contact.

Phase 1 (knot_collider.py) settled the same-helix channel: V(s) = 0 exactly
until chord-sharing contact at s = 4, which is repulsive (+4.4..+12.8 lam=1)
and reversible. This phase collides knots approaching on DIFFERENT chains:

  1. knot A sits at a chain-A window (static);
  2. candidate crossing chains are found by walking bounded BC stretches
     from every tet near A's chord (4 dropped-vertex choices per tet),
     keeping stretches whose closest approach to A is small and whose local
     tangent makes a distinct angle with A's;
  3. for each crossing: B's washboard S_single_B(j) is measured by walking
     B along the stretch on a CLEAN manifold (no A) -- which also validates
     the path template-wise; then, on a fresh manifold with A present, B is
     created far out and slid through the crossing, recording

         V(j) = S_total - S_single_A - S_single_B(j)

     exactly (washboard-corrected, the Phase-1 lesson).

Outcomes per crossing: V profile, closest-approach distance and shared
vertices, merge species / jam / transparent pass-through, and the unfuse
test on merges. Angle 0 (same chain) is Phase 1's channel; this script
reports whatever angles the R-phase chain geometry actually offers.
"""
import argparse
import json
import os
import sys
from collections import Counter

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
for _p in ("../../python", "../../scripts", "../../tools",
           "../../scripts/defect_dynamics", "."):   # archived: siblings stayed in dd/
    sys.path.insert(0, os.path.join(_HERE, _p))
import discrete_differential_geometry as ddg
import worm_slide as ws
from chain_select import chain_for_run, add_chain_args
from worm_helix import bc_orbit
from knot_collider import make_knot, slide_along, census_str, ESTAR
from cocycle_check import reference_frac_positions
import defect_state as ds


def minimg(d, box):
    return d - np.round(d / box) * box


def walk_stretch(m, window, n):
    """Bounded sliding-window walk: n steps forward from `window`.
    Returns the vertex list v[0..n+3] (window k = v[k..k+3])."""
    w = list(window)
    seq = [w[0]]
    for _ in range(n):
        face = w[1:]
        a1, a2 = m.face_apexes(*face)
        w = face + [a2 if a1 == w[0] else a1]
        seq.append(w[0])
    seq.extend(w[1:])
    return seq


def chord_mid(rp, chain, j, box):
    a, b = rp[chain[j]], rp[chain[j + 4]]
    return (a + minimg(b - a, box) / 2) % box


def tangent(rp, chain, j, box):
    d = minimg(rp[chain[j + 8]] - rp[chain[j]], box)
    n = np.linalg.norm(d)
    return d / n if n > 0 else d


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--ref",
                    default="data/tcp_reference/T3_R_m4_N57984.mfd")
    ap.add_argument("--mcell", type=int, default=4)
    ap.add_argument("--windowA", type=int, default=1624,
                    help="chain-A window for knot A (1624 = a rung-70 site)")
    ap.add_argument("--seed-tet", type=int, default=0)
    ap.add_argument("--kmax", type=int, default=6,
                    help="B starts kmax slides from the crossing")
    ap.add_argument("--crossings", type=int, default=6)
    ap.add_argument("--search-r", type=float, default=1.0,
                    help="cells: tet search radius around A's chord")
    ap.add_argument("--min-angle", type=float, default=0.0,
                    help="deg: only test crossings at least this steep. "
                         "Near-90 crossings pair a LEFT- and RIGHT-handed "
                         "spiral, a geometry the default closest-approach "
                         "picker never reaches")
    ap.add_argument("--out", default=None)
    add_chain_args(ap, default=None)
    args = ap.parse_args()

    box = float(args.mcell)
    rp = np.asarray(reference_frac_positions("r", args.mcell))

    ref = ddg.Manifold.load(args.ref, 3)
    F = np.asarray(ref.facets())
    _cc, _kcls, _seq, chain_prov = chain_for_run(
        args.ref, F, args.chain_class, seed_tet=args.seed_tet)
    chainA = [int(x) for x in _seq]
    LA = len(chainA)
    bA = args.windowA % LA
    em0 = ws.edeg_dict(ref)
    midA = chord_mid(rp, chainA, bA, box)
    tanA = tangent(rp, chainA, bA, box)
    chordA_verts = {chainA[(bA + i) % LA] for i in range(5)}

    # ---- candidate crossing stretches
    dv = minimg(rp - midA, box)
    near = np.where((dv ** 2).sum(1) < args.search_r ** 2)[0]
    nearset = set(int(x) for x in near)
    cand_tets = [t for t in F if all(int(v) in nearset for v in t)]
    print(f"chain A: L={LA}; knot A window {bA}, chord mid {midA.round(2)}, "
          f"{len(cand_tets)} tets within {args.search_r} cells")

    need = 4 * (args.kmax + 3)
    seen = set()
    crossings = []
    for t in cand_tets:
        for drop in range(4):
            w0 = [int(t[drop])] + [int(x) for x in t if int(x) != int(t[drop])]
            stretch = walk_stretch(ref, w0, 2 * need)
            # closest approach of stretch chords to A's chord midpoint
            best = None
            for j in range(0, len(stretch) - 8):
                mid = chord_mid(rp, stretch, j, box)
                d = float(np.linalg.norm(minimg(mid - midA, box)))
                if best is None or d < best[1]:
                    best = (j, d)
            jX, dmin = best
            if not (4 <= jX <= len(stretch) - 12):
                continue
            tanB = tangent(rp, stretch, max(0, min(jX, len(stretch) - 9)),
                           box)
            ang = float(np.degrees(np.arccos(
                min(1.0, abs(float(np.dot(tanA, tanB)))))))
            key = frozenset(stretch[jX:jX + 5])
            if key in seen:
                continue
            seen.add(key)
            if ang < args.min_angle:
                continue
            shared = len(set(stretch[jX:jX + 5]) & chordA_verts)
            # handedness: torsion sign of successive chain-edge vectors
            def hand(ch, j):
                e = [minimg(rp[ch[j + k + 1]] - rp[ch[j + k]], box)
                     for k in range(3)]
                return float(np.sign(np.dot(np.cross(e[0], e[1]), e[2])))
            h = float(np.mean([hand(stretch, j) for j in
                               range(max(0, jX - 8),
                                     min(len(stretch) - 4, jX + 8), 4)]))
            crossings.append({"w0": w0, "stretch": stretch, "jX": jX,
                              "dmin": dmin, "angle": ang, "shared": shared,
                              "hand": h})
    crossings.sort(key=lambda c: (c["dmin"], -c["angle"]))
    # keep an angle-diverse selection
    picked, used_ang = [], []
    for c in crossings:
        if len(picked) >= args.crossings:
            break
        if any(abs(c["angle"] - a) < 5 and abs(c["dmin"] - d) < 0.05
               for a, d in used_ang):
            continue
        if c["jX"] < 4 * args.kmax + 4:
            continue
        picked.append(c)
        used_ang.append((c["angle"], c["dmin"]))
    print(f"{len(crossings)} distinct crossing stretches; testing "
          f"{len(picked)} (angle-diverse):")
    for c in picked:
        print(f"   angle {c['angle']:5.1f} deg  dmin {c['dmin']:.3f} cells  "
              f"shared verts {c['shared']}  handed {c.get('hand', 0):+.1f}  "
              f"jX {c['jX']}")
    del ref

    # ---- S_single_A
    m = ddg.Manifold.load(args.ref, 3)
    make_knot(m, chainA, bA)
    S1A = ws.dS_between(em0, ws.edeg_dict(m), estar=ESTAR)
    cenA = census_str(m)
    print(f"\nknot A: sig {cenA[0][1]}, sumZ {cenA[0][2]}, "
          f"S1A = {S1A:.4f}")
    del m

    results = []
    for ci, c in enumerate(picked):
        stretch, jX = c["stretch"], c["jX"]
        j0 = jX - 4 * args.kmax
        # -- pass 1: B alone -- washboard + path validation
        m = ddg.Manifold.load(args.ref, 3)
        try:
            make_knot(m, stretch, j0)
        except AssertionError as e:
            print(f"\n== crossing {ci}: B creation failed ({e}); skip")
            del m
            continue
        Ssingle = {j0: ws.dS_between(em0, ws.edeg_dict(m), estar=ESTAR)}
        jB, ok = j0, True
        while jB < jX:
            jB2, nv, nc = slide_along(m, stretch, jB, fwd=True)
            if jB2 is None:
                print(f"\n== crossing {ci}: B path jams at j={jB} even "
                      f"without A ({nv} template-valid); skip")
                ok = False
                break
            jB = jB2
            Ssingle[jB] = ws.dS_between(em0, ws.edeg_dict(m), estar=ESTAR)
        del m
        if not ok:
            continue
        # -- pass 2: A present, B slides through the crossing
        m = ddg.Manifold.load(args.ref, 3)
        make_knot(m, chainA, bA)
        make_knot(m, stretch, j0)
        jB = j0
        rows, outcome = [], None
        print(f"\n== crossing {ci}: angle {c['angle']:.1f} deg, "
              f"dmin {c['dmin']:.3f}, shared {c['shared']}")
        while True:
            S = ws.dS_between(em0, ws.edeg_dict(m), estar=ESTAR)
            cen = census_str(m)
            base = Ssingle.get(jB)
            V = S - S1A - base if base is not None else None
            mid = chord_mid(rp, stretch, jB, box)
            d = float(np.linalg.norm(minimg(mid - midA, box)))
            rows.append({"j_rel": jB - jX, "d_cells": d, "V": V,
                         "n_components": len(cen),
                         "census": [list(map(list, x[:2])) + [x[2]]
                                    for x in cen]})
            print(f"   j-jX={jB-jX:+4d}  d={d:5.3f}  V="
                  f"{(f'{V:+.6f}' if V is not None else '  n/a  ')}  "
                  f"comps={len(cen)} "
                  f"{[x[1] for x in cen] if len(cen) <= 3 else ''}")
            if len(cen) < 2:
                outcome = {"type": "merged", "at_jrel": jB - jX,
                           "d": d, "species_f": list(cen[0][0]),
                           "sig": list(cen[0][1]), "sumZ": cen[0][2],
                           "V": V}
                break
            if jB >= jX:
                outcome = {"type": "transparent" if all(
                    (r["V"] or 0) == 0 or abs(r["V"] or 0) < 1e-9
                    for r in rows) else "contact-no-merge",
                    "at_jrel": jB - jX, "d": d, "V": V}
                break
            jB2, nv, nc = slide_along(m, stretch, jB, fwd=True)
            if jB2 is None:
                outcome = {"type": "jammed", "at_jrel": jB - jX, "d": d,
                           "template_valid": nv, "V": V}
                print(f"   JAMMED ({nv} template-valid)")
                break
            jB = jB2
        # unfuse test on merges
        rev = None
        if outcome["type"] == "merged":
            e3, _ = ws.knot_chords(m)
            rev = {"unfuses": False}
            for ch in e3:
                for cs, mv in ws.candidate_slides(m, tuple(ch)):
                    recs = ws.apply_slide(m, mv)
                    if recs is None:
                        continue
                    cen2 = census_str(m)
                    if (len(cen2) == 2 and all(
                            list(x[1]).count(3) == 1 for x in cen2)):
                        rev["unfuses"] = True
                        break
                    ws.undo_slide(m, recs)
                if rev["unfuses"]:
                    break
        results.append({"crossing": ci, "angle": c["angle"],
                        "dmin": c["dmin"], "shared": c["shared"],
                        "hand": c.get("hand"),
                        "chordA": [int(chainA[bA]),
                                   int(chainA[(bA + 4) % LA])],
                        "chordB": [int(stretch[jX]),
                                   int(stretch[jX + 4])],
                        "rows": rows, "outcome": outcome, "reversal": rev})
        del m

    print(f"\n{'='*70}\nSUMMARY (angle deg, dmin cells, V in lam=1 units):")
    for r in results:
        far = [abs(x["V"]) for x in r["rows"]
               if x["V"] is not None and x["j_rel"] < -8]
        o = r["outcome"]
        print(f"  angle {r['angle']:5.1f}  dmin {r['dmin']:.3f}  "
              f"shared {r['shared']}  hand {r.get('hand', 0):+.1f}  "
              f"max|V(far)| "
              f"{max(far, default=0):.2e}  -> {o['type']}"
              + (f" f={tuple(o['species_f'])} sumZ={o['sumZ']} "
                 f"V={o['V']:+.3f} unfuses="
                 f"{r['reversal']['unfuses']}"
                 if o["type"] == "merged" else
                 f" (V at crossing = "
                 f"{o.get('V') if o.get('V') is not None else 'n/a'})"))

    if args.out:
        with open(args.out, "w") as fh:
            json.dump({"ref": args.ref, "windowA": bA,
                       "S1A": S1A, "results": results}, fh, indent=1)
        print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
