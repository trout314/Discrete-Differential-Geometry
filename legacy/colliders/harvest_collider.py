#!/usr/bin/env python3
"""Plan C: collide a knot into HARVESTED long-lived defects, in the frozen
thermal background where they actually live.

Load a thermal snapshot, never call the sampler (the background is static;
all dynamics are deterministic directed slides), pick a target complex from
the census, find a BC chain whose stretch passes near it through otherwise
clear crystal, create knot C far out on the chain, slide it in.

THE WASHBOARD REFERENCE IN A THERMAL STATE: C's single-knot site energies
cannot be measured in place (the target cannot be removed), but the
background is EXACTLY pristine crystal away from its defects (no-halo,
100.00% site registration measured by defect_halo). So every stretch vertex
maps to a pristine-crystal vertex by position, the same stretch is rebuilt
on the pristine reference, and the washboard is measured there. Then

    V(j) = [S(bg + C at j) - S(bg)]  -  W_pristine(site of j)

is exact, and the no-halo prediction is V = 0 until C's frame touches the
target's skin -- now tested in a thermal background against nature's own
long-lived species.
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
from discrete_differential_geometry import cocycle as coc
import worm_slide as ws
from knot_collider import make_knot, slide_along, ESTAR
from crossing_collider import walk_stretch, chord_mid, minimg
from cocycle_check import reference_frac_positions
import defect_state as ds

CELL = 1.0e6


def census_full(m):
    st = ds.DefectState(m)
    comps = []
    for cx in st.components(broad=True):
        comps.append({"f": tuple(st.induced_shape(cx.verts)["f"]),
                      "sig": tuple(cx.sig),
                      "sumZ": st.total_coordination(cx.verts),
                      "verts": tuple(cx.verts)})
    return st, comps


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--snap", default="data/mgas/lam40x_snap15000.mfd")
    ap.add_argument("--ref",
                    default="data/tcp_reference/T3_R_m4_N57984.mfd")
    ap.add_argument("--mcell", type=int, default=4)
    ap.add_argument("--target", type=int, default=-1,
                    help="census index; -1 = list targets and run the "
                         "largest")
    ap.add_argument("--kmax", type=int, default=4)
    ap.add_argument("--geoms", type=int, default=2,
                    help="approach geometries per target")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    box = float(args.mcell)
    pbox = CELL * args.mcell
    rp = np.asarray(reference_frac_positions("r", args.mcell))

    # ---- pristine reference + site registry (defect_halo recipe)
    ref = ddg.Manifold.load(args.ref, 3)
    em0 = ws.edeg_dict(ref)
    raw0 = np.round(CELL * rp).astype(np.int64) % int(pbox)
    site_of = {tuple(raw0[v]): v for v in range(len(rp))}

    # ---- thermal background
    m = ddg.Manifold.load(args.snap, 3)
    st, comps = census_full(m)
    em_bg = ws.edeg_dict(m)
    e1, w1, _ = coc.load_cocycle(args.snap.replace(".mfd", ".cocycle.npz"))
    e1 = np.asarray(e1)
    pos = coc.tree_positions(e1, np.asarray(w1), int(e1.max()) + 1)
    pos = pos[0]
    verts_all = sorted(st.v2t)
    cand = [v for v in verts_all if v < len(raw0)]
    offs = (pos[cand] - raw0[cand]) % int(pbox)
    off = np.array(Counter(map(tuple, offs)).most_common(1)[0][0])

    def to_site(v):
        return site_of.get(tuple((pos[v] - off) % int(pbox)))

    # inverse: site -> thermal label (unique in the crystal region; in the
    # defect skins a site can host a recycled label -- first wins, and any
    # ambiguity there surfaces as a jam, i.e. as data)
    inv_site = {}

    def cellpos(v):
        return (pos[v] - off) % int(pbox) / CELL

    for v in verts_all:
        s_ = to_site(v)
        if s_ is not None and s_ not in inv_site:
            inv_site[s_] = v
    nmatch = sum(1 for v in verts_all if to_site(v) is not None)
    print(f"background {os.path.basename(args.snap)}: "
          f"{len(verts_all)} vertices, site match "
          f"{100 * nmatch / len(verts_all):.2f}%")
    comps.sort(key=lambda c: -len(c["verts"]))
    print("census (by size):")
    for i, c in enumerate(comps):
        cen = np.mean([cellpos(v) for v in c["verts"]], axis=0)
        print(f"  [{i}] f={c['f']}  sig={ws.__name__ and c['sig']}  "
              f"sumZ={c['sumZ']}  at ~{cen.round(2)}")
    tsel = [args.target] if args.target >= 0 else [0]
    defect_pos = np.array([cellpos(v) for v in st.defect])

    results = []
    for ti in tsel:
        tgt = comps[ti]
        tverts = set(tgt["verts"])
        tpos = np.array([cellpos(v) for v in tgt["verts"]])
        print(f"\n{'='*70}\nTARGET [{ti}] f={tgt['f']} sumZ={tgt['sumZ']}")

        # ---- find approach stretches: near target, path otherwise clear
        tcen = tpos[0] + minimg(tpos - tpos[0], box).mean(0)
        near = {v for v in verts_all
                if float(np.linalg.norm(minimg(
                    np.array(cellpos(v)) - tcen, box))) < 1.3}
        F = np.asarray(m.facets())
        cand_tets = [t for t in F if all(int(x) in near for x in t)
                     and not any(int(x) in tverts for x in t)]
        need = 4 * (args.kmax + 3)
        geoms, seen = [], set()
        for t in cand_tets:
            for drop in range(4):
                w0 = [int(t[drop])] + [int(x) for x in t
                                       if int(x) != int(t[drop])]
                # CANONICAL chain: mirror the seed to pristine sites and
                # walk THERE. Walking the thermal manifold near the target
                # routes through deformed tets, and that walk's mirror is
                # not a pristine BC chain (the first version's washboard
                # jams). The canonical chain is walkable on pristine by
                # construction; its thermal image deviates only inside the
                # target's skin, where deviations ARE the experiment.
                w0s = [to_site(v) for v in w0]
                if any(x is None for x in w0s):
                    continue
                try:
                    stp = walk_stretch(ref, w0s, 2 * need)
                except Exception:
                    continue
                best = None
                for j in range(0, len(stp) - 8):
                    a, b = rp[stp[j]], rp[stp[j + 4]]
                    mid = a + minimg(b - a, box) / 2
                    d = float(np.linalg.norm(
                        minimg(mid - tpos, box), axis=1).min())
                    if best is None or d < best[1]:
                        best = (j, d)
                jX, dmin = best
                if dmin > 0.6 or not (4 * args.kmax + 4 <= jX
                                      <= len(stp) - 12):
                    continue
                stx = [inv_site.get(s_) for s_ in stp]
                path = stx[jX - 4 * args.kmax:jX + 8]
                if any(v is None for v in path):
                    continue
                if any(v in st.defect and v not in tverts for v in path):
                    continue
                key = frozenset(stp[jX:jX + 5])
                if key in seen:
                    continue
                seen.add(key)
                geoms.append((stp, stx, jX, dmin))
        geoms.sort(key=lambda g: g[3])
        print(f"  {len(geoms)} clear approach stretches (dmin<0.6); "
              f"running {min(args.geoms, len(geoms))}")

        ran = 0
        for gi, (stp, stx, jX, dmin) in enumerate(geoms):
            if ran >= args.geoms:
                break
            j0 = jX - 4 * args.kmax
            # ---- washboard on the canonical pristine chain
            mp = ddg.Manifold.load(args.ref, 3)
            try:
                make_knot(mp, stp, j0)
            except AssertionError:
                print(f"  geom {gi}: pristine creation fails; skip")
                del mp
                continue
            W = {j0: ws.dS_between(em0, ws.edeg_dict(mp), estar=ESTAR)}
            jC, ok = j0, True
            while jC < jX:
                jC2, nv, nc = slide_along(mp, stp, jC, fwd=True)
                if jC2 is None:
                    ok = False
                    break
                jC = jC2
                W[jC] = ws.dS_between(em0, ws.edeg_dict(mp), estar=ESTAR)
            del mp
            if not ok:
                print(f"  geom {gi}: pristine washboard jams; skip")
                continue
            ran += 1
            # ---- the real approach, frozen thermal background
            mt = ddg.Manifold.load(args.snap, 3)
            try:
                make_knot(mt, stx, j0)
            except AssertionError:
                print(f"  geom {gi}: creation blocked in background; skip")
                del mt
                continue
            jC = j0
            rows, outcome = [], None
            print(f"\n  -- geom {gi}: dmin {dmin:.3f}")
            while True:
                S = ws.dS_between(em_bg, ws.edeg_dict(mt), estar=ESTAR)
                V = S - W.get(jC, np.nan)
                st2, cen2 = census_full(mt)
                ntot = len(cen2)
                # which complex is C? the one containing C's chord verts
                cverts = {stx[jC], stx[jC + 4]}
                mine = [c for c in cen2 if cverts & set(c["verts"])]
                fused = bool(mine and tverts & set(mine[0]["verts"]))
                a, b = rp[stp[jC]], rp[stp[jC + 4]]
                mid = a + minimg(b - a, box) / 2
                d = float(np.linalg.norm(minimg(mid - tpos, box),
                                         axis=1).min())
                rows.append({"j_rel": jC - jX, "d": d, "V": float(V),
                             "ncomp": ntot})
                print(f"     j-jX={jC-jX:+4d} d={d:5.3f} "
                      f"V={V:+.6f} comps={ntot}"
                      + (f"  FUSED f={mine[0]['f']} "
                         f"sumZ={mine[0]['sumZ']}" if fused else ""))
                del st2
                if fused:
                    outcome = {"type": "merged", "f": list(mine[0]["f"]),
                               "sig": list(mine[0]["sig"]),
                               "sumZ": mine[0]["sumZ"], "V": float(V),
                               "d": d}
                    break
                if jC >= jX:
                    outcome = {"type": "no-merge", "V": float(V), "d": d}
                    break
                jC2, nv, nc = slide_along(mt, stx, jC, fwd=True)
                if jC2 is None:
                    outcome = {"type": "jammed", "V": float(V), "d": d,
                               "template_valid": nv}
                    print(f"     JAMMED ({nv} template-valid)")
                    break
                jC = jC2
            # unfuse test
            if outcome["type"] == "merged":
                e3, _ = ws.knot_chords(mt)
                outcome["unfuses"] = False
                for ch in e3:
                    if not ({stx[jC], stx[jC + 4]} & set(ch)):
                        continue
                    for cs, mv in ws.candidate_slides(mt, tuple(ch)):
                        recs = ws.apply_slide(mt, mv)
                        if recs is None:
                            continue
                        _, c3 = census_full(mt)
                        back = any(set(c["verts"]) == tverts for c in c3)
                        if back:
                            outcome["unfuses"] = True
                            break
                        ws.undo_slide(mt, recs)
                    if outcome["unfuses"]:
                        break
            results.append({"target": ti, "target_f": list(tgt["f"]),
                            "target_sumZ": tgt["sumZ"], "geom": gi,
                            "dmin": dmin, "rows": rows,
                            "outcome": outcome})
            del mt

    print(f"\n{'='*70}\nSUMMARY")
    for r in results:
        o = r["outcome"]
        far = [abs(x["V"]) for x in r["rows"]
               if not np.isnan(x["V"]) and x["j_rel"] < -8]
        print(f"  target {tuple(r['target_f'])} sumZ={r['target_sumZ']} "
              f"geom {r['geom']}: max|V(far)| {max(far, default=0):.2e} "
              f"-> {o['type']}"
              + (f" f={tuple(o['f'])} sumZ={o['sumZ']} V={o['V']:+.4f} "
                 f"unfuses={o.get('unfuses')}"
                 if o["type"] == "merged" else f" (V={o.get('V'):+.4f})"))

    if args.out:
        with open(args.out, "w") as fh:
            json.dump({"snap": args.snap, "results": results}, fh,
                      indent=1, default=str)
        print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
