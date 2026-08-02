#!/usr/bin/env python3
"""B.1 + B.2: spectroscopy of a manufactured two-knot compound, then
three-body chemistry -- a third knot slid into it.

Builds the free-association compound the way Phase 2 did (knot A static,
knot B slid in on a crossing chain until the census fuses, V = 0), then:

B.1  SPECTROSCOPY of the compound:
     - both chords are still degree 3: enumerate ALL 12 slide slots on each
       chord, classify every outcome -- translate-as-unit / unfuse /
       rearrange to another species / invalid -- with exact dS for each.
     - if any outcome translates the compound whole, that is composite
       transport, and its washboard is measured over a few steps.

B.2  THREE-BODY: find crossing chains near the compound, slide knot C in
     (washboard-corrected against pristine single-knot energies), record
     V(j) on approach, the contact product (15 vertices? free again?), and
     the unfuse test. The growth-by-accretion question: is association
     still free at the second step, 5 -> 10 -> 15?
"""
import argparse
import json
import os
import sys
from collections import Counter

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
for _p in ("../../python", "../../scripts", "../../tools", "."):
    sys.path.insert(0, os.path.join(_HERE, _p))
import discrete_differential_geometry as ddg
import worm_slide as ws
from chain_select import chain_for_run, add_chain_args
from worm_helix import bc_orbit
from knot_collider import make_knot, slide_along, census_str, ESTAR
from crossing_collider import walk_stretch, chord_mid, tangent, minimg
from cocycle_check import reference_frac_positions
import defect_state as ds


def find_crossings(ref, rp, box, mid0, tan0, avoid, search_r, need):
    """All distinct crossing stretches near mid0 (same recipe as Phase 2)."""
    F = np.asarray(ref.facets())
    dv = minimg(rp - mid0, box)
    near = np.where((dv ** 2).sum(1) < search_r ** 2)[0]
    nearset = set(int(x) for x in near)
    cands = [t for t in F if all(int(v) in nearset for v in t)]
    seen, out = set(), []
    for t in cands:
        for drop in range(4):
            w0 = [int(t[drop])] + [int(x) for x in t
                                   if int(x) != int(t[drop])]
            st = walk_stretch(ref, w0, 2 * need)
            best = None
            for j in range(0, len(st) - 8):
                mid = chord_mid(rp, st, j, box)
                d = float(np.linalg.norm(minimg(mid - mid0, box)))
                if best is None or d < best[1]:
                    best = (j, d)
            jX, dmin = best
            if not (4 * 8 + 4 <= jX <= len(st) - 12):
                continue
            key = frozenset(st[jX:jX + 5])
            if key in seen:
                continue
            seen.add(key)
            if avoid and (set(st[jX:jX + 5]) & avoid):
                shared = len(set(st[jX:jX + 5]) & avoid)
            else:
                shared = 0
            tb = tangent(rp, st, min(jX, len(st) - 9), box)
            ang = float(np.degrees(np.arccos(
                min(1.0, abs(float(np.dot(tan0, tb)))))))
            out.append({"stretch": st, "jX": jX, "dmin": dmin,
                        "angle": ang, "shared": shared})
    out.sort(key=lambda c: c["dmin"])
    return out


def classify_state(m, before_verts=None):
    st = ds.DefectState(m)
    comps = st.components(broad=True)
    out = []
    for cx in comps:
        f = tuple(st.induced_shape(cx.verts)["f"])
        out.append({"f": f, "sig": tuple(cx.sig),
                    "sumZ": st.total_coordination(cx.verts),
                    "verts": tuple(cx.verts)})
    del st
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--ref",
                    default="data/tcp_reference/T3_R_m4_N57984.mfd")
    ap.add_argument("--mcell", type=int, default=4)
    ap.add_argument("--windowA", type=int, default=1624)
    ap.add_argument("--seed-tet", type=int, default=0)
    ap.add_argument("--kmax", type=int, default=6)
    ap.add_argument("--three-body", type=int, default=2,
                    help="number of knot-C approach geometries to run")
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
    need = 4 * (args.kmax + 3)

    crossings = find_crossings(ref, rp, box, midA, tanA, chordA_verts,
                               1.0, need)
    del ref

    # ---- build the compound: first crossing that fuses at V = 0
    print("building the compound...")
    m = None
    build = None
    for c in crossings:
        mm = ddg.Manifold.load(args.ref, 3)
        make_knot(mm, chainA, bA)
        S1A = ws.dS_between(em0, ws.edeg_dict(mm), estar=ESTAR)
        st_, jX = c["stretch"], c["jX"]
        j0 = jX - 4 * args.kmax
        try:
            make_knot(mm, st_, j0)
        except AssertionError:
            del mm
            continue
        jB, ok = j0, True
        fused = False
        while True:
            cen = census_str(mm)
            if len(cen) < 2:
                fused = True
                break
            if jB >= jX:
                break
            jB2, nv, nc = slide_along(mm, st_, jB, fwd=True)
            if jB2 is None:
                ok = False
                break
            jB = jB2
        if fused:
            S = ws.dS_between(em0, ws.edeg_dict(mm), estar=ESTAR)
            cen = census_str(mm)
            if cen[0][0][0] == 10:      # want the 10-vertex free compound
                m = mm
                build = {"angle": c["angle"], "dmin": c["dmin"],
                         "S1A": S1A, "S_compound": S,
                         "f": cen[0][0], "sig": cen[0][1],
                         "sumZ": cen[0][2]}
                break
        del mm
    assert m is not None, "no fusing crossing found"
    print(f"  compound built: f={build['f']} sig={build['sig']} "
          f"sumZ={build['sumZ']}  (crossing angle {build['angle']:.1f}, "
          f"S_compound = {build['S_compound']:.4f}, "
          f"binding V = {build['S_compound'] - 2 * build['S1A']:+.6f} "
          f"vs 2 knots... using per-site S1s from phase 2 conventions)")
    em_compound = ws.edeg_dict(m)
    S_comp = build["S_compound"]
    state0 = classify_state(m)
    comp_verts = set(state0[0]["verts"])

    # ======== B.1: slide spectroscopy of the compound ========
    print(f"\n{'='*70}\nB.1 SPECTROSCOPY: all slide slots on both chords")
    e3, _ = ws.knot_chords(m)
    spec = []
    for chord in e3:
        for cs, mv in ws.candidate_slides(m, tuple(chord)):
            recs = ws.apply_slide(m, mv)
            if recs is None:
                spec.append({"chord": list(chord), "outcome": "invalid"})
                continue
            S1 = ws.dS_between(em0, ws.edeg_dict(m), estar=ESTAR)
            cen = classify_state(m)
            nc = len(cen)
            vs = set(cen[0]["verts"]) if nc == 1 else set()
            if nc == 2 and all(list(x["sig"]).count(3) == 1 for x in cen):
                oc = "unfuse"
            elif nc == 1 and cen[0]["f"] == build["f"] \
                    and vs != comp_verts:
                oc = "TRANSLATE"
            elif nc == 1 and cen[0]["f"] == build["f"]:
                oc = "same-verts-shape"
            elif nc == 1:
                oc = f"rearrange->{cen[0]['f']}"
            else:
                oc = f"split->{[x['f'] for x in cen]}"
            spec.append({"chord": list(chord), "outcome": oc,
                         "dS": S1 - S_comp,
                         "products": [list(x["f"]) for x in cen]})
            ws.undo_slide(m, recs)
    oc_count = Counter(x["outcome"].split("->")[0] for x in spec)
    print(f"  {len(spec)} slot attempts on {len(e3)} chords: {dict(oc_count)}")
    for x in spec:
        if x["outcome"] != "invalid":
            print(f"   chord {x['chord']}: {x['outcome']:<28s} "
                  f"dS = {x['dS']:+8.4f}")

    # ======== B.2: three-body -- knot C into the compound ========
    print(f"\n{'='*70}\nB.2 THREE-BODY: knot C slid into the compound")
    # compound centroid
    pv = rp[sorted(comp_verts)]
    rel = minimg(pv - pv[0], box)
    midC = (pv[0] + rel.mean(0)) % box
    ref2 = ddg.Manifold.load(args.ref, 3)
    cands = find_crossings(ref2, rp, box, midC, tanA, comp_verts, 1.2, need)
    del ref2
    tb_results = []
    tried = 0
    for c in cands:
        if tried >= args.three_body:
            break
        st_, jX = c["stretch"], c["jX"]
        j0 = jX - 4 * args.kmax
        # washboard for C on PRISTINE crystal
        mm = ddg.Manifold.load(args.ref, 3)
        try:
            make_knot(mm, st_, j0)
        except AssertionError:
            del mm
            continue
        Ssingle = {j0: ws.dS_between(em0, ws.edeg_dict(mm), estar=ESTAR)}
        jC, ok = j0, True
        while jC < jX:
            jC2, nv, nc = slide_along(mm, st_, jC, fwd=True)
            if jC2 is None:
                ok = False
                break
            jC = jC2
            Ssingle[jC] = ws.dS_between(em0, ws.edeg_dict(mm), estar=ESTAR)
        del mm
        if not ok:
            continue
        tried += 1
        # approach on a COPY of the compound state
        mc = ddg.Manifold.load(args.ref, 3)
        # rebuild compound exactly: A + B walk (same recipe)
        make_knot(mc, chainA, bA)
        stB, jXB = None, None
        # replay the build crossing
        for cb in crossings:
            if abs(cb["angle"] - build["angle"]) < 1e-9 \
                    and abs(cb["dmin"] - build["dmin"]) < 1e-9:
                stB, jXB = cb["stretch"], cb["jX"]
                break
        jB = jXB - 4 * args.kmax
        make_knot(mc, stB, jB)
        while len(census_str(mc)) >= 2:
            jB2, nv, nc = slide_along(mc, stB, jB, fwd=True)
            if jB2 is None:
                break
            jB = jB2
        try:
            make_knot(mc, st_, j0)
        except AssertionError:
            print(f"  C creation blocked at this geometry; skip")
            del mc
            tried -= 1
            continue
        jC = j0
        rows, outcome = [], None
        print(f"\n  -- C approach: angle-to-A {c['angle']:.1f} deg, "
              f"dmin {c['dmin']:.3f} to compound centroid")
        while True:
            S = ws.dS_between(em0, ws.edeg_dict(mc), estar=ESTAR)
            cen = classify_state(mc)
            base = Ssingle.get(jC)
            V = S - S_comp - base if base is not None else None
            mid = chord_mid(rp, st_, jC, box)
            d = float(np.linalg.norm(minimg(mid - midC, box)))
            rows.append({"j_rel": jC - jX, "d": d, "V": V,
                         "ncomp": len(cen),
                         "fs": [list(x["f"]) for x in cen]})
            print(f"     j-jX={jC-jX:+4d} d={d:5.3f} V="
                  f"{(f'{V:+.6f}' if V is not None else '  n/a  ')} "
                  f"comps={len(cen)} {[x['f'] for x in cen]}")
            if len(cen) < 2:
                outcome = {"type": "merged", "f": list(cen[0]["f"]),
                           "sig": list(cen[0]["sig"]),
                           "sumZ": cen[0]["sumZ"], "V": V}
                break
            if jC >= jX:
                outcome = {"type": "no-merge", "V": V, "d": d}
                break
            jC2, nv, nc = slide_along(mc, st_, jC, fwd=True)
            if jC2 is None:
                outcome = {"type": "jammed", "d": d, "V": V,
                           "template_valid": nv}
                print(f"     JAMMED ({nv} template-valid)")
                break
            jC = jC2
        tb_results.append({"angle": c["angle"], "dmin": c["dmin"],
                           "rows": rows, "outcome": outcome})
        del mc

    print(f"\n{'='*70}\nSUMMARY")
    print(f"compound: f={build['f']} sumZ={build['sumZ']}")
    print(f"B.1: {dict(oc_count)}")
    for r in tb_results:
        o = r["outcome"]
        far = [abs(x["V"]) for x in r["rows"]
               if x["V"] is not None and x["j_rel"] < -8]
        print(f"B.2 angle {r['angle']:5.1f} dmin {r['dmin']:.3f}: "
              f"max|V(far)| {max(far, default=0):.2e} -> {o['type']}"
              + (f" f={tuple(o['f'])} sumZ={o['sumZ']} V={o['V']:+.4f}"
                 if o["type"] == "merged" else ""))

    if args.out:
        with open(args.out, "w") as fh:
            json.dump({"build": {k: (list(v) if isinstance(v, tuple) else v)
                                 for k, v in build.items()},
                       "spectroscopy": spec,
                       "three_body": tb_results}, fh, indent=1, default=str)
        print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
