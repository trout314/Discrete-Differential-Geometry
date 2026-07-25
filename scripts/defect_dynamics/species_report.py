#!/usr/bin/env python3
"""What defect species exist in an ensemble of states? -- a reporting front end.

Pools the defect complexes over a set of snapshots and reports the species
distribution, following notes/CONVENTIONS.md: the standard state-report block
(provenance, action + couplings, legality, curvature, census, certification
status) followed by the species table.

A species is labelled by its ILLEGAL-EDGE SIGNATURE -- the sorted multiset of
illegal edge degrees inside the complex -- with any non-FK coordination
decorations reported separately, so "(3,4,4)" keeps exactly the meaning it has
always had. Also reported per species: complex size s, closed-star size, total
disclination charge sum w(e) = sum (6 - deg e), and curvature charge Q_c
against the cell-mean reference.

Usage:
  species_report.py --glob 'data/mgas/lam40*_snap*.mfd' --lam 0.40 \
      --zleg 0.6 --cimp 1.0 --etarget 5.105025 \
      --provenance 'm4 r-crystal -> lam40 chains'
"""
import argparse
import glob as globmod
import os
import re
import sys
from collections import Counter, defaultdict

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
for _p in ("../../python", "../../scripts", "."):
    sys.path.insert(0, os.path.join(_HERE, _p))
import discrete_differential_geometry as ddg
import defect_state as ds


def sig_str(sig, width=30):
    """Compact signature: (3,4,4) stays literal; long runs of one degree
    collapse to d x n so the table stays readable."""
    if not sig:
        return "()"
    c = Counter(sig)
    if len(c) == 1 and len(sig) > 4:
        d, n = next(iter(c.items()))
        return f"({d} x{n})"
    s = "(" + ",".join(str(x) for x in sig) + ")"
    return s if len(s) <= width else s[:width - 3] + "...)"


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--glob", required=True, action="append",
                    help="snapshot glob (repeatable)")
    ap.add_argument("--lam", type=float, default=None)
    ap.add_argument("--edq-only", action="store_true",
                    help="the action is the volume pin + EDQ only (no n6); "
                         "reports lambda as lam_EDQ, never as lam")
    ap.add_argument("--zleg", type=float, default=0.6)
    ap.add_argument("--cimp", type=float, default=1.0)
    ap.add_argument("--etarget", type=float, default=None)
    ap.add_argument("--provenance", default=None)
    ap.add_argument("--certified", default=None,
                    help="gate numbers if certified; omitted => PROVISIONAL")
    ap.add_argument("--group",
                    choices=("sig", "induced", "decorated", "full"),
                    default="sig",
                    help="sig: illegal-edge signature (historical). "
                         "induced: isomorphism class of the subcomplex induced "
                         "by the defect's vertices -- its bare shape. "
                         "decorated: that subcomplex with ambient EDGE DEGREES "
                         "on its edges, which is the shape plus the physics. "
                         "full: edge degrees AND n6 on the vertices, which "
                         "pins the coordination Z and hence the curvature at "
                         "every vertex.")
    ap.add_argument("--top", type=int, default=25)
    ap.add_argument("--min-count", type=int, default=1)
    ap.add_argument("--out", default=None, help="write the table as JSON too")
    args = ap.parse_args()

    files = sorted({f for g in args.glob for f in globmod.glob(g)})
    if not files:
        sys.exit(f"no snapshots matched {args.glob}")
    chains = sorted({re.sub(r"_snap\d+\.mfd$", "", os.path.basename(f))
                     for f in files})

    spec = Counter()
    size_of, star_of, chg_of = (defaultdict(list) for _ in range(3))
    shape_of = {}
    n6_of = defaultdict(Counter)
    sig_of = defaultdict(Counter)
    nodes_of = defaultdict(Counter)
    inexact = 0
    ncx_per, nill_per, ffk_per, fe_per, edeg_per = [], [], [], [], []
    n3_states = []
    ncomp = 0
    for f in files:
        m = ddg.Manifold.load(f, 3)
        st = ds.DefectState(m)
        c = ds.census(st)
        comps = st.components()
        q = st.vertex_charges()
        qbar = float(np.mean([q[v] for v in st.v2t]))
        ncomp += len(comps)
        ncx_per.append(len(comps))
        nill_per.append(c["n_illegal"])
        ffk_per.append(c["legalvert_fk"])
        fe_per.append(c["legaledge"])
        edeg_per.append(c["mean_edeg"])
        n3_states.append(m.num_facets)
        for cx in comps:
            if args.group in ("induced", "decorated", "full"):
                fac = st.induced_facets(cx.verts)
                if args.group in ("decorated", "full"):
                    vc, ec = st.decorations(cx.verts)
                    ck, exact = ds.canonical_key(
                        fac, ecolor=ec,
                        vcolor=(vc if args.group == "full" else None))
                    n6_of[ck][tuple(sorted(st.n6[v] for v in cx.verts))] += 1
                else:
                    ck, exact = ds.canonical_key(fac)
                if not exact:
                    inexact += 1
                key = ck
                shape_of[key] = st.induced_shape(cx.verts)
                sig_of[key][cx.sig] += 1
                nodes_of[key][cx.nodes] += 1
            else:
                key = cx.key
            spec[key] += 1
            size_of[key].append(len(cx.verts))
            star_of[key].append(len(st.star(cx.verts)))
            chg_of[key].append(st.complex_charge(cx.verts, q, qbar))
        del st, m

    def pm(x):
        return f"{np.mean(x):.4g} +- {np.std(x):.2g}"

    print("=" * 78)
    print("DEFECT SPECIES CENSUS")
    print("=" * 78)
    print(f"Provenance   : {args.provenance or '(not stated)'}")
    print(f"               {len(files)} snapshots over {len(chains)} chains")
    print(f"               chains: {', '.join(chains[:8])}"
          + (" ..." if len(chains) > 8 else ""))
    if args.lam is not None:
        if args.edq_only:
            print(f"Action       : volume pin + EDQ only, lam_EDQ = {args.lam}"
                  f"   (NOT comparable to full-action lambda)")
        else:
            print(f"Action       : full (EDQ + n6), lambda = {args.lam}, "
                  f"zleg = {args.zleg}, cimp = {args.cimp}")
    else:
        print("Action       : (not stated)")
    if args.etarget:
        print(f"Curvature    : e* = {args.etarget:.6f}   <e> = {pm(edeg_per)}")
    else:
        print(f"Curvature    : <e> = {pm(edeg_per)}")
    print(f"Legality     : f_FK = {pm(ffk_per)}   f_e = {pm(fe_per)}")
    print(f"Census       : n_ill = {pm(nill_per)}   N_cx = {pm(ncx_per)}   "
          f"N_3 = {int(np.mean(n3_states))}")
    print("Crystallinity: not computed (crystal_grains covering map is "
          "disproportionately expensive here; run defect_census.py for "
          "N_gr / f_G1 / f_reg)")
    print(f"Certification: {args.certified or 'PROVISIONAL (uncertified)'}")
    print(f"Totals       : {ncomp} complexes, {len(spec)} distinct species")
    print()

    if args.group in ("induced", "decorated", "full"):
        what = {"induced": "induced subcomplex",
                "decorated": "induced subcomplex DECORATED BY EDGE DEGREE",
                "full": "induced subcomplex DECORATED BY EDGE DEGREE + n6"
                }[args.group]
        print(f"Grouped by ISOMORPHISM CLASS of the {what} "
              f"({len(spec)} classes)")
        if inexact:
            print(f"  WARNING: {inexact} complexes hit the canonical-form "
                  f"search limit -- their keys are strong invariants, not "
                  f"certificates")
        hdr = (f"{'#':<4} {'f=(v,e,t,T)':<17} {'facets':<7} {'degseq (1-skel)':<22} "
               f"{'n':>5} {'%':>6} {'star':>6} {'Q_c (rad)':>15}")
        print(hdr)
        print("-" * len(hdr))
        rows = []
        for i, (k, n) in enumerate(spec.most_common()):
            if n < args.min_count or i >= args.top:
                break
            sh = shape_of[k]
            dq = str(sh["degseq"])
            nf = len(k[0]) if args.group in ("decorated", "full") else len(k)
            rows.append({"class": i, "f": list(sh["f"]), "n_facets": nf,
                         "degseq": list(sh["degseq"]), "n": n,
                         "frac": n / ncomp,
                         "sigs": {str(a): b for a, b in sig_of[k].items()},
                         "nodes": {str(a): b for a, b in nodes_of[k].items()},
                         "star": float(np.mean(star_of[k])),
                         "Q_c": float(np.mean(chg_of[k])),
                         "Q_c_sd": float(np.std(chg_of[k]))})
            print(f"{i:<4} {str(sh['f']):<17} {nf:<7} {dq[:21]:<22} "
                  f"{n:5d} {100*n/ncomp:5.1f}% {np.mean(star_of[k]):6.0f} "
                  f"{np.mean(chg_of[k]):7.2f}+-{np.std(chg_of[k]):<6.2f}")
        # how well does shape determine the historical label, and vice versa?
        print()
        print("class -> illegal-edge signatures it contains "
              "(is it finer or coarser than the old label?)")
        for i, (k, n) in enumerate(spec.most_common(min(args.top, 12))):
            sigs = sig_of[k]
            pretty = ", ".join(f"{sig_str(a)}x{b}" for a, b in sigs.most_common(4))
            print(f"   class {i:<3} (n={n:4d}): {pretty}"
                  + ("  ..." if len(sigs) > 4 else ""))
        if args.group in ("decorated", "full"):
            print()
            print("n6 multiset within each class"
                  + (" -- would adding n6 to the key split it?"
                     if args.group == "decorated" else " (now IN the key)"))
            nsplit = sum(1 for k in spec if len(n6_of[k]) > 1)
            for i, (k, n) in enumerate(spec.most_common(min(args.top, 10))):
                ms = n6_of[k]
                pretty = ", ".join(f"{a}x{b}" for a, b in ms.most_common(3))
                print(f"   class {i:<3} (n={n:4d}): {len(ms)} distinct -> "
                      f"{pretty}" + ("  ..." if len(ms) > 3 else ""))
            print(f"   classes whose n6 is NOT constant: {nsplit} of "
                  f"{len(spec)}  -> adding n6 would split those")
        rev = defaultdict(set)
        for k in spec:
            for a in sig_of[k]:
                rev[a].add(k)
        multi = {a: len(v) for a, v in rev.items() if len(v) > 1}
        print(f"\n   signatures realised by more than one shape: "
              f"{len(multi)} of {len(rev)}")
        for a, c in sorted(multi.items(), key=lambda x: -x[1])[:6]:
            print(f"      {sig_str(a)} appears as {c} distinct shapes")
    else:
        hdr = (f"{'illegal-edge signature':<20} {'nodes':<9} {'n':>5} {'%':>6} "
               f"{'s (verts)':>13} {'star':>7} {'sum w':>6} {'Q_c (rad)':>15}")
        print(hdr)
        print("-" * len(hdr))
        rows = []
        for k, n in spec.most_common():
            if n < args.min_count or len(rows) >= args.top:
                break
            sig, nodes = k
            w = sum(6 - d for d in sig)
            rows.append({"sig": list(sig), "nodes": list(nodes), "n": n,
                         "frac": n / ncomp, "size": float(np.mean(size_of[k])),
                         "star": float(np.mean(star_of[k])), "sum_w": w,
                         "Q_c": float(np.mean(chg_of[k])),
                         "Q_c_sd": float(np.std(chg_of[k]))})
            print(f"{sig_str(sig):<20} {str(nodes) if nodes else '':<9} {n:5d} "
                  f"{100*n/ncomp:5.1f}% "
                  f"{np.mean(size_of[k]):6.1f}+-{np.std(size_of[k]):<5.1f} "
                  f"{np.mean(star_of[k]):7.0f} {w:6d} "
                  f"{np.mean(chg_of[k]):7.2f}+-{np.std(chg_of[k]):<6.2f}")

    # --- collapsed views
    by3 = Counter()
    bydec = Counter()
    sigsrc = ([(a, b) for k in spec for a, b in
               ((s_, sig_of[k][s_] * 1) for s_ in sig_of[k])]
              if args.group in ("induced", "decorated", "full")
              else [(k[0], n) for k, n in spec.items()])
    nodesrc = ([(a, b) for k in spec for a, b in
                ((s_, nodes_of[k][s_]) for s_ in nodes_of[k])]
               if args.group in ("induced", "decorated", "full")
               else [(k[1], n) for k, n in spec.items()])
    for sig, n in sigsrc:
        by3[sig.count(3)] += n
    for nodes, n in nodesrc:
        bydec[bool(nodes)] += n
    for _unused in ():
        pass
    print()
    print("slide handles (a slide needs a degree-3 chord):")
    for k in sorted(by3):
        print(f"   complexes with {k} degree-3 edge(s): {by3[k]:5d}  "
              f"{100*by3[k]/ncomp:5.1f}%")
    print("non-FK coordination decoration:")
    print(f"   undecorated                    : {bydec[False]:5d}  "
          f"{100*bydec[False]/ncomp:5.1f}%")
    print(f"   carrying >=1 non-FK vertex     : {bydec[True]:5d}  "
          f"{100*bydec[True]/ncomp:5.1f}%")

    if args.out:
        import json
        with open(args.out, "w") as fh:
            json.dump({"files": len(files), "chains": chains,
                       "complexes": ncomp, "species": len(spec),
                       "lam": args.lam, "edq_only": args.edq_only,
                       "rows": rows}, fh, indent=1)
        print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
