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
            spec[cx.key] += 1
            size_of[cx.key].append(len(cx.verts))
            star_of[cx.key].append(len(st.star(cx.verts)))
            chg_of[cx.key].append(st.complex_charge(cx.verts, q, qbar))
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
        row = {"sig": list(sig), "nodes": list(nodes), "n": n,
               "frac": n / ncomp, "size": float(np.mean(size_of[k])),
               "star": float(np.mean(star_of[k])), "sum_w": w,
               "Q_c": float(np.mean(chg_of[k])),
               "Q_c_sd": float(np.std(chg_of[k]))}
        rows.append(row)
        print(f"{sig_str(sig):<20} {str(nodes) if nodes else '':<9} {n:5d} "
              f"{100*n/ncomp:5.1f}% "
              f"{np.mean(size_of[k]):6.1f}+-{np.std(size_of[k]):<5.1f} "
              f"{np.mean(star_of[k]):7.0f} {w:6d} "
              f"{np.mean(chg_of[k]):7.2f}+-{np.std(chg_of[k]):<6.2f}")

    # --- collapsed views
    by3 = Counter()
    bydec = Counter()
    for (sig, nodes), n in spec.items():
        by3[sig.count(3)] += n
        bydec[bool(nodes)] += n
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
