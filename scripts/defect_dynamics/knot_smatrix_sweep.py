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
  3. crystal symmetry + determinism: configurations related by ANY
     automorphism of the triangulation give identical outcomes move-for-move,
     so classes are deduped by exact Aut keys (symmetry.config_key), anchored
     on knot A's frame.

     This replaced a rounded relative-position hash of a BALL of surrounding
     vertices, which was wrong in both directions. It keyed on the
     ENVIRONMENT rather than the object, so two chords sitting in congruent
     balls but oriented differently within them were MERGED -- measured on
     R m4 chainA: of the 130 A-windows it dropped, 87 were not Aut-equivalent
     to the representative kept, i.e. silently missing classes in a sweep
     declaring exhaustiveness. It also missed every point-group equivalence,
     giving 680 classes where there are 271. The two partitions cross-cut;
     the old one was neither a refinement nor a coarsening of the truth.

     One consequence to expect: with A's window PINNED there is no residual
     symmetry to exploit, because Aut acts freely on frames -- any g fixing
     A's frame is the identity. So for a fixed A, distinct B geometries are
     now all kept, where the ball hash used to fuse some of them. Expect more
     collisions per A-class and fewer A-classes (271 vs 680); the sharing that
     still pays is the washboard cache, whose key is B's start frame alone and
     so is shared across all A-classes.

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
from chain_select import chain_for_run, add_chain_args
from worm_helix import bc_orbit
from knot_collider import make_knot, slide_along, ESTAR
from crossing_collider import walk_stretch, chord_mid, tangent, minimg
from cocycle_check import reference_frac_positions
import defect_state as ds


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


def walk_bidir(m, window, n):
    """Stretch spliced from n steps BACKWARD and n steps forward of
    `window`, so the seed sits at index ~n and a crossing found anywhere
    along the chain leaves room for B to start 4*kmax upstream. Walking
    forward only (the first version) silently dropped every chain whose
    in-ball tets sit near its closest-approach point -- including the
    same-chain channel, whose upstream seeds are always outside the search
    ball (jX >= 4*kmax+4 needs seeds ~1.1+ cells away).

    The backward walk uses the reversed window; its vertex sequence is
    reversed and spliced so index arithmetic (window k = v[k..k+3]) is
    preserved across the seam."""
    fwd = walk_stretch(m, window, n)
    bwd = walk_stretch(m, list(reversed(window)), n)
    # bwd = [w3', w2', ...] walking the chain the other way; reversed(bwd)
    # ends with the seed window's vertices, which fwd starts with -- drop
    # the 4-vertex overlap.
    return list(reversed(bwd))[:-4] + fwd


class Ledger:
    """Undo ledger: every mutation is an invertible move, so state is
    restored by unwinding instead of reloading (Manifold reloads leak ~70MB
    each under the conservative D GC -- the original sweep died of OOM at
    148 collisions)."""

    def __init__(self, m):
        self.m = m
        self.ops = []

    def knot(self, chain, j):
        chord = make_knot(self.m, chain, j)
        L = len(chain)
        face = sorted([chain[(j + 1) % L], chain[(j + 2) % L],
                      chain[(j + 3) % L]])
        self.ops.append(("bis", chord, face))
        return chord

    def slide(self, chain, j, fwd):
        recs = []
        out = slide_along(self.m, chain, j, fwd, recs_out=recs)
        if out[0] is not None:
            self.ops.append(("slide", recs[0]))
        return out

    def unwind(self):
        for kind, *pay in reversed(self.ops):
            if kind == "slide":
                ws.undo_slide(self.m, pay[0])
            else:
                chord, face = pay
                self.m.do_bistellar_move(sorted(chord), face)
        self.ops = []


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
    ap.add_argument("--max-collisions", type=int, default=45000)
    ap.add_argument("--audit-every", type=int, default=25,
                    help="verify exact facet-set restoration every N "
                         "collisions")
    ap.add_argument("--out", required=True)
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
    em0 = ws.edeg_dict(ref)
    # pristine adjacency, for the transparent-by-theorem prefilter
    nbr = {}
    for t in F:
        for a in t:
            nbr.setdefault(int(a), set()).update(int(x) for x in t)

    # ---- A-window classes over the whole orbit, EXACT under Aut
    winA = lambda b: [chainA[(b + i) % LA] for i in range(4)]
    aclasses = {}
    for b in range(0, LA - 12, 4):
        k = _cc.config_key([winA(b)])
        if k not in aclasses:
            aclasses[k] = b
    reps = sorted(aclasses.values())
    # randomized (seeded) order: a capped run is then a UNIFORM sample of
    # A-classes, so the fingerprint-saturation curve speaks for the whole
    # population, not for one end of the orbit
    rng = np.random.default_rng(20260725)
    rng.shuffle(reps)
    if args.limit_classes:
        reps = reps[:args.limit_classes]
    reps = reps[args.start_class:]
    print(f"orbit L={LA}: {LA//4} windows -> {len(aclasses)} A-window "
          f"classes (exact, under the full Aut); running {len(reps)}",
          flush=True)

    need = 4 * (args.kmax + 3)
    t0 = time.time()
    ncoll = 0
    out_rows = []
    capped = False
    # two long-lived manifolds, restored by ledger unwind (never reloaded)
    mW = ddg.Manifold.load(args.ref, 3)      # washboard passes
    mC = ddg.Manifold.load(args.ref, 3)      # collision passes
    fs0_n = mW.num_facets
    wb_cache = {}                     # B start-frame Aut key -> Ssb dict
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
                stx = walk_bidir(ref, w0, 2 * need)
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
                jk = _cc.config_key([winA(bA),
                                     [int(stx[jX + i]) for i in range(4)]])
                if jk in seen:
                    continue
                seen.add(jk)
                geoms.append((stx, jX, dmin))
        # SAME-CHAIN channels, both head-on orientations, by pure chain
        # arithmetic (the ball search cannot find these: their upstream
        # seeds lie outside any local ball). B is created 4*kmax+8 chain
        # steps from A's chord and slides to jX = the chord ONE step short
        # of A's -- merges register during the approach exactly as in
        # Phase 1 (contact at chord-sharing, s=4).
        jXs = 4 * args.kmax + 8
        for orient in (+1, -1):
            if orient == +1:
                s0 = bA - 4 - jXs          # seg[i] = chainA[s0 + i]
                seg = [chainA[(s0 + i) % LA]
                       for i in range(jXs + 4 * args.kmax + 25)]
            else:
                r0 = bA + 8 + jXs          # seg[i] = chainA[r0 - i]
                seg = [chainA[(r0 - i) % LA]
                       for i in range(jXs + 4 * args.kmax + 25)]
            jk = _cc.config_key([winA(bA),
                                 [int(seg[jXs + i]) for i in range(4)]])
            if jk not in seen:
                seen.add(jk)
                geoms.append((seg, jXs, 0.0))
        # TRANSPARENT BY THEOREM: if no slide frame along B's whole path
        # (frames reach stx[j..j+8] for j in [j0, jX]) touches A's five
        # vertices or their pristine neighbours, the configurations
        # decompose at every step (no-halo additivity): V = 0 identically,
        # no merge, no jam-by-A. Recorded without running.
        A5 = [chainA[(bA + i) % LA] for i in range(5)]
        reach = set(A5)
        for a in A5:
            reach |= nbr.get(a, set())
        run_geoms = []
        n_theorem = 0
        for stx, jX, dmin in geoms:
            j0 = jX - 4 * args.kmax
            path = set(stx[max(0, j0):jX + 9])
            if path & reach:
                run_geoms.append((stx, jX, dmin))
            else:
                n_theorem += 1
                out_rows.append({"aclass": ai, "windowA": bA,
                                 "dmin": dmin,
                                 "outcome": "transparent-by-theorem"})
        run_geoms.sort(key=lambda g: g[2])
        geoms = run_geoms
        # knot A stays in place for the whole class (its own ledger)
        LA_ledger = Ledger(mC)
        LA_ledger.knot(chainA, bA)
        S1A = ws.dS_between(em0, ws.edeg_dict(mC), estar=ESTAR)
        print(f"[A-class {ai}/{len(reps)}] window {bA}: "
              f"{len(geoms)} to run + {n_theorem} transparent by theorem "
              f"(incl. same-chain)", flush=True)

        # washboard cache for this A class is per-geometry (B's chain)
        for stx, jX, dmin in geoms:
            if ncoll >= args.max_collisions:
                capped = True
                break
            ncoll += 1
            j0 = jX - 4 * args.kmax
            # pass 1: B washboard -- cached by the EXACT Aut key of B's own
            # start frame. The BC walk is deterministic, so that frame
            # determines the whole segment, and the segment length is fixed;
            # the same chain geometry recurs across A-classes.
            seg = stx[max(0, j0 - 1):jX + 9]
            wkey = _cc.config_key([[int(stx[j0 + i]) for i in range(4)]])
            if wkey in wb_cache:
                Ssb = wb_cache[wkey]
                if Ssb is None:
                    out_rows.append({"aclass": ai, "windowA": bA,
                                     "dmin": dmin,
                                     "outcome": "path-jams-clean"})
                    continue
                Ssb = {j0 + k: v for k, v in Ssb.items()}
            else:
                LW = Ledger(mW)
                try:
                    LW.knot(stx, j0)
                except AssertionError:
                    LW.unwind()
                    wb_cache[wkey] = None
                    continue
                Ssb = {j0: ws.dS_between(em0, ws.edeg_dict(mW),
                                         estar=ESTAR)}
                jB, ok = j0, True
                while jB < jX:
                    jB2, nv, nc = LW.slide(stx, jB, fwd=True)
                    if jB2 is None:
                        ok = False
                        break
                    jB = jB2
                    Ssb[jB] = ws.dS_between(em0, ws.edeg_dict(mW),
                                            estar=ESTAR)
                LW.unwind()
                if not ok:
                    wb_cache[wkey] = None
                    out_rows.append({"aclass": ai, "windowA": bA,
                                     "dmin": dmin,
                                     "outcome": "path-jams-clean"})
                    continue
                wb_cache[wkey] = {k - j0: v for k, v in Ssb.items()}
            # pass 2: collision (A already in place), ledger-unwound
            LC = Ledger(mC)
            m = mC
            try:
                LC.knot(stx, j0)
            except AssertionError:
                LC.unwind()
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
                jB2, nv, nc = LC.slide(stx, jB, fwd=True)
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
                        good = (len(c2) == 2 and all(
                            list(x["sig"]).count(3) == 1 for x in c2))
                        ws.undo_slide(m, recs)      # ALWAYS undo
                        if good:
                            res["unfuses"] = True
                            break
                    if res["unfuses"]:
                        break
            res.update({"aclass": ai, "windowA": bA, "dmin": dmin,
                        "vfar_max": vfar_max})
            out_rows.append(res)
            LC.unwind()
            if args.audit_every and ncoll % args.audit_every == 0:
                assert mW.num_facets == fs0_n, "unwind failed (mW)"
                assert mC.num_facets == fs0_n + 1, \
                    "unwind failed (mC: A + crystal expected)"
        LA_ledger.unwind()
        if ws.edeg_dict(mC) != em0:
            sys.exit("class-level unwind failed")
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
