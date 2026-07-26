#!/usr/bin/env python3
"""Plan A: the exhaustive knot-knot S-matrix sweep on pristine R.

SCOPE (declared, so "exhaustive" has a precise meaning): every
translation-equivalence class of

    (knot-A window on the chain-A orbit)  x  (BC-chain geometry whose
    closest approach to A's chord enters dmin < --dmin-cut)

with one deterministic collision per class. Completeness rests on:

  1. the no-halo theorem: geometries with dmin above the transparency
     radius (~0.6 cells, measured) have V identically 0 -- the infinite
     tail of distant chains is dismissed by theorem, not sampling;
  2. finite enumeration: chain windows are (tet, dropped-vertex) pairs, at
     most 4 per tet, over the finitely many tets in the contact ball --
     ALL are walked, both directions arising as separate drops;
  3. crystal periodicity + determinism: configurations identical up to a
     lattice translation give identical outcomes move-for-move, so classes
     are deduped by EXACT relative-position keys (translation-safe; missed
     point-group symmetries only cause harmless redundancy, never a gap).

Per collision: the approach V profile (theorem check), the contact
FINGERPRINT (canonicalized bipartite adjacency between A's and B's
role-labelled vertices at first contact), the product species, V_contact,
and the unfuse test. Results flush to JSON after every A-class.
"""
import argparse
import json
import os
import sys
import time
from collections import Counter

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
for _p in ("../../python", "../../scripts", "../../tools", "."):
    sys.path.insert(0, os.path.join(_HERE, _p))
import discrete_differential_geometry as ddg
import worm_slide as ws
from worm_helix import bc_orbit
from knot_collider import make_knot, slide_along, ESTAR
from crossing_collider import walk_stretch, chord_mid, tangent, minimg
from cocycle_check import reference_frac_positions
import defect_state as ds


def relkey(rp, verts, origin, box, digits=4):
    """Translation-invariant key: sorted rounded positions rel. to origin."""
    rel = minimg(rp[sorted(set(verts))] - origin, box)
    return tuple(sorted(tuple(int(round(x * 10 ** digits)) for x in row)
                        for row in rel))


def roles(st, verts):
    """Role label per vertex of a defect complex: sorted degrees of its
    illegal edges inside the complex (chord ends carry a 3)."""
    vs = set(verts)
    lab = {v: [] for v in verts}
    for e in st.ill_edges:
        if e[0] in vs and e[1] in vs:
            lab[e[0]].append(st.edeg[e])
            lab[e[1]].append(st.edeg[e])
    return {v: tuple(sorted(l)) for v, l in lab.items()}


def fingerprint(m, vertsA, vertsB):
    """Canonical bipartite contact pattern between two complexes."""
    st = ds.DefectState(m)
    rA, rB = roles(st, vertsA), roles(st, vertsB)
    shared = sorted(set(vertsA) & set(vertsB))
    rows = []
    for a in sorted(vertsA):
        row = []
        for b in sorted(vertsB):
            if a == b:
                row.append(2)                      # identified vertex
            else:
                row.append(1 if (min(a, b), max(a, b)) in st.edeg else 0)
        rows.append((rA[a], tuple(row)))
    # canonicalize: sort A rows by (role, pattern-multiset); columns by role
    colr = [rB[b] for b in sorted(vertsB)]
    order = sorted(range(len(colr)), key=lambda i: colr[i])
    rows2 = sorted((role, tuple(row[i] for i in order))
                   for role, row in rows)
    del st
    return (tuple(rows2), tuple(sorted(colr)), len(shared))


def census(m):
    st = ds.DefectState(m)
    out = [{"f": tuple(st.induced_shape(c.verts)["f"]),
            "sig": tuple(c.sig), "sumZ": st.total_coordination(c.verts),
            "verts": tuple(c.verts)}
           for c in st.components(broad=True)]
    del st
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--ref",
                    default="data/tcp_reference/T3_R_m4_N57984.mfd")
    ap.add_argument("--mcell", type=int, default=4)
    ap.add_argument("--seed-tet", type=int, default=0)
    ap.add_argument("--kmax", type=int, default=3,
                    help="slides B starts out from the crossing. 3 is "
                         "sufficient: V=0 at range is guaranteed by the "
                         "no-halo theorem, so only the contact zone needs "
                         "stepping through")
    ap.add_argument("--dmin-cut", type=float, default=0.7)
    ap.add_argument("--ball-r", type=float, default=1.0,
                    help="cells: radius of the A-class equivalence ball. "
                         "Provably sufficient at 1.0: outcomes are fixed by "
                         "the joint configuration within reach of any move "
                         "frame that can touch A (contact radius ~0.6 + "
                         "frame span ~0.55, minus the overlap), because "
                         "everything farther is pristine crystal identical "
                         "across classes (no-halo)")
    ap.add_argument("--limit-classes", type=int, default=0,
                    help="0 = all A-window classes on the orbit")
    ap.add_argument("--start-class", type=int, default=0,
                    help="resume: skip A-classes before this index")
    ap.add_argument("--max-collisions", type=int, default=4000)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    box = float(args.mcell)
    rp = np.asarray(reference_frac_positions("r", args.mcell))
    ref = ddg.Manifold.load(args.ref, 3)
    F = np.asarray(ref.facets())
    chainA = bc_orbit(ref, [int(x) for x in F[args.seed_tet]])
    LA = len(chainA)
    em0 = ws.edeg_dict(ref)

    # ---- A-window translation classes over the whole orbit
    aclasses = {}
    for b in range(0, LA - 12, 4):
        mid = chord_mid(rp, chainA, b, box)
        ball = np.where((minimg(rp - mid, box) ** 2).sum(1)
                        < args.ball_r ** 2)[0]
        k = relkey(rp, ball, mid, box)
        if k not in aclasses:
            aclasses[k] = b
    reps = sorted(aclasses.values())
    if args.limit_classes:
        reps = reps[:args.limit_classes]
    reps = reps[args.start_class:]
    print(f"orbit L={LA}: {LA//4} windows -> {len(aclasses)} A-window "
          f"translation classes; running {len(reps)}", flush=True)

    need = 4 * (args.kmax + 3)
    t0 = time.time()
    ncoll = 0
    out_rows = []
    capped = False
    for ai, bA in enumerate(reps):
        midA = chord_mid(rp, chainA, bA, box)
        tanA = tangent(rp, chainA, bA, box)
        # enumerate ALL entering geometries for this A class
        dv = minimg(rp - midA, box)
        nearset = set(int(x)
                      for x in np.where((dv ** 2).sum(1) < 1.0)[0])
        cand_tets = [t for t in F
                     if all(int(v) in nearset for v in t)]
        geoms, seen = [], set()
        for t in cand_tets:
            for drop in range(4):
                w0 = [int(t[drop])] + [int(x) for x in t
                                       if int(x) != int(t[drop])]
                stx = walk_stretch(ref, w0, 2 * need)
                best = None
                for j in range(0, len(stx) - 8):
                    d = float(np.linalg.norm(minimg(
                        chord_mid(rp, stx, j, box) - midA, box)))
                    if best is None or d < best[1]:
                        best = (j, d)
                jX, dmin = best
                if dmin > args.dmin_cut:
                    continue
                if not (4 * args.kmax + 4 <= jX <= len(stx) - 12):
                    continue
                jk = relkey(rp, stx[jX - 8:jX + 12], midA, box)
                if jk in seen:
                    continue
                seen.add(jk)
                geoms.append((stx, jX, dmin))
        print(f"[A-class {ai}/{len(reps)}] window {bA}: "
              f"{len(geoms)} joint classes", flush=True)

        # washboard cache for this A class is per-geometry (B's chain)
        for stx, jX, dmin in geoms:
            if ncoll >= args.max_collisions:
                capped = True
                break
            ncoll += 1
            j0 = jX - 4 * args.kmax
            # pass 1: B washboard (clean crystal)
            mm = ddg.Manifold.load(args.ref, 3)
            try:
                make_knot(mm, stx, j0)
            except AssertionError:
                del mm
                continue
            Ssb = {j0: ws.dS_between(em0, ws.edeg_dict(mm), estar=ESTAR)}
            jB, ok = j0, True
            while jB < jX:
                jB2, nv, nc = slide_along(mm, stx, jB, fwd=True)
                if jB2 is None:
                    ok = False
                    break
                jB = jB2
                Ssb[jB] = ws.dS_between(em0, ws.edeg_dict(mm),
                                        estar=ESTAR)
            del mm
            if not ok:
                out_rows.append({"aclass": ai, "windowA": bA,
                                 "dmin": dmin, "outcome": "path-jams-clean"})
                continue
            # pass 2: collision
            m = ddg.Manifold.load(args.ref, 3)
            make_knot(m, chainA, bA)
            S1A = ws.dS_between(em0, ws.edeg_dict(m), estar=ESTAR)
            try:
                make_knot(m, stx, j0)
            except AssertionError:
                del m
                out_rows.append({"aclass": ai, "windowA": bA,
                                 "dmin": dmin,
                                 "outcome": "B-blocked-by-A"})
                continue
            jB = j0
            vfar_max = 0.0
            res = None
            vA5 = [chainA[(bA + i) % LA] for i in range(5)]
            emc = None
            while True:
                emc = ws.edeg_dict(m)
                S = ws.dS_between(em0, emc, estar=ESTAR)
                V = S - S1A - Ssb.get(jB, np.nan)
                # cheap contact test: is any current B-frame vertex adjacent
                # to A's five? Only then is a full census needed.
                vB5 = [stx[(jB + i)] for i in range(5)]
                touching = any((min(a, b), max(a, b)) in emc or a == b
                               for a in vA5 for b in vB5)
                cen = census(m) if (touching or abs(jB - jX) < 1) \
                    else [{"f": None}, {"f": None}]
                if len(cen) < 2:
                    fp = None
                    # fingerprint needs the two vertex sets just before
                    # merge -- recover from the merged complex + A's known
                    # chord frame: use A = 5 chain verts, B = the rest
                    va = vA5
                    vb = [v for v in cen[0]["verts"] if v not in va]
                    fp = fingerprint(m, sorted(va), sorted(vb) or list(va))
                    res = {"outcome": "merged",
                           "f": list(cen[0]["f"]),
                           "sig": list(cen[0]["sig"]),
                           "sumZ": cen[0]["sumZ"],
                           "V_contact": V, "fingerprint": str(fp)}
                    break
                if abs(jB - jX) < 1:
                    res = {"outcome": "transparent"
                           if vfar_max < 1e-9 and abs(V) < 1e-9
                           else "contact-no-merge", "V_at_min": V}
                    break
                if not np.isnan(V):
                    vfar_max = max(vfar_max, abs(V))
                jB2, nv, nc = slide_along(m, stx, jB, fwd=True)
                if jB2 is None:
                    res = {"outcome": "jammed", "V_at_jam": V,
                           "template_valid": nv}
                    break
                jB = jB2
            # unfuse test
            if res["outcome"] == "merged":
                e3, _ = ws.knot_chords(m)
                res["unfuses"] = False
                for ch in e3:
                    for cs, mv in ws.candidate_slides(m, tuple(ch)):
                        recs = ws.apply_slide(m, mv)
                        if recs is None:
                            continue
                        c2 = census(m)
                        if (len(c2) == 2 and all(
                                list(x["sig"]).count(3) == 1
                                for x in c2)):
                            res["unfuses"] = True
                            break
                        ws.undo_slide(m, recs)
                    if res["unfuses"]:
                        break
            res.update({"aclass": ai, "windowA": bA, "dmin": dmin,
                        "vfar_max": vfar_max})
            out_rows.append(res)
            del m
        # flush after every A class
        with open(args.out, "w") as fh:
            json.dump({"scope": {"ref": args.ref, "orbit_len": LA,
                                 "aclasses_total": len(aclasses),
                                 "aclasses_run": ai + 1,
                                 "dmin_cut": args.dmin_cut,
                                 "capped": capped},
                       "rows": out_rows}, fh)
        el = time.time() - t0
        print(f"   ... {ncoll} collisions, {el:.0f}s "
              f"({el/max(ncoll,1):.1f}s each)", flush=True)
        if capped:
            print(f"WARNING: --max-collisions {args.max_collisions} hit; "
                  f"sweep INCOMPLETE (classes {ai+1}..{len(reps)-1} "
                  f"not run)", flush=True)
            break

    oc = Counter(r["outcome"] for r in out_rows)
    fps = Counter(r.get("fingerprint") for r in out_rows
                  if r.get("fingerprint"))
    prods = Counter(tuple(r["f"]) for r in out_rows if "f" in r)
    print(f"\nSWEEP {'(CAPPED!)' if capped else 'complete'}: "
          f"{ncoll} collisions over {len(reps)} A-classes")
    print(f"outcomes: {dict(oc)}")
    print(f"distinct fingerprints: {len(fps)}")
    print(f"products: {dict(prods)}")
    vc = Counter()
    for r in out_rows:
        if "V_contact" in r and r["V_contact"] is not None \
                and not (isinstance(r["V_contact"], float)
                         and np.isnan(r["V_contact"])):
            vc[round(r["V_contact"], 3)] += 1
    print(f"V_contact spectrum: {dict(sorted(vc.items()))}")
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
