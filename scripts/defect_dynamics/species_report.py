#!/usr/bin/env python3
"""What defect species exist in an ensemble of states? -- a reporting front end.

Pools the defect complexes over a set of snapshots and reports the species
distribution, following notes/CONVENTIONS.md: the standard state-report block
(provenance, action + couplings, legality, curvature, census, certification
status) followed by the species table.

By default a species is the isomorphism class of the subcomplex INDUCED by a
defect's vertices, decorated with ambient edge degrees on its edges and n6 on
its vertices (--group full). That is the reporting default because it is the
only key that determines the curvature charge: it pins the coordination Z at
every vertex, and then

    Q_c = (sum Z) * (pi - 3 theta) + 6 * n_verts * theta

exactly (notes/CONVENTIONS.md 1.1). Coarser keys remain available: --group
decorated drops n6, --group induced drops all decoration, and --group sig is
the historical illegal-edge signature, kept for continuity with earlier runs.

Also reported per species: complex size s, closed-star size, total
coordination sum Z, and the RAW curvature charge Q_c = sum q_R, with the
background-subtracted value alongside. Raw is primary: it depends only on the
complex's own vertices, so identical complexes give identical values, whereas
subtracting a state mean ties every measurement to the whole configuration
(measured sd over 60 classes: 4e-16 raw vs 2e-3 subtracted).

Pass --flicker with a crystal_flicker.py JSON to mark which species a single
2->3 on the PRISTINE crystal can produce, and to print the residual census with
that background separated. This matters more than it sounds: every single-2->3
excitation of a pristine crystal has f = (5,10,9,3), so the (3,4,4) knot is the
shape of one Pachner move against a perfect crystal rather than an excitation
of the melt, and at lam=0.40 96% of defect BIRTHS are exactly that. The flag is
an upper bound, not an attribution -- see load_flicker's docstring.

Usage:
  species_report.py --glob 'data/mgas/lam40*_snap*.mfd' --lam 0.40 \
      --zleg 0.6 --cimp 1.0 --etarget 5.105025 \
      --flicker flicker_R_m4.json \
      --provenance 'm4 r-crystal -> lam40 chains'
  ... add --group sig to reproduce a pre-2026-07-25 species table.
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
import crystal_flicker as cf


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
                    choices=("full", "decorated", "induced", "sig"),
                    default="full",
                    help="full (default): the induced subcomplex decorated by "
                         "ambient EDGE DEGREE on edges and n6 on vertices. This "
                         "is the only key that determines Q_c, so it is the "
                         "reporting default (CONVENTIONS.md 1.1). "
                         "sig: illegal-edge signature (historical, for "
                         "continuity with earlier runs). "
                         "induced: isomorphism class of the subcomplex induced "
                         "by the defect's vertices -- its bare shape. "
                         "decorated: edge degrees only, without n6.")
    ap.add_argument("--top", type=int, default=25)
    ap.add_argument("--min-count", type=int, default=8,
                    help="hide species seen fewer than this many times. The "
                         "default suits a few hundred complexes; the tail is "
                         "mostly one-off configurations. Whatever is hidden is "
                         "always reported as a count, never dropped silently. "
                         "Use 1 to see everything.")
    ap.add_argument("--flicker", default=None,
                    help="crystal_flicker.py JSON for the SAME crystal. Marks "
                         "each species reachable by a single 2->3 on the "
                         "pristine crystal, and reports the census with that "
                         "background separated out. Note this is an UPPER "
                         "BOUND on flicker contamination: a snapshot has no "
                         "move history, so 'reachable' is not 'was produced "
                         "that way'. For real attribution use "
                         "flicker_fraction.py on an event stream.")
    ap.add_argument("--out", default=None, help="write the table as JSON too")
    args = ap.parse_args()

    flick, flick_meta = (set(), None)
    if args.flicker:
        flick, flick_meta = cf.load_flicker(args.flicker, args.group)
        if not flick:
            sys.exit(f"{args.flicker} has no keys for --group {args.group}; "
                     f"regenerate it with the current crystal_flicker.py")

    files = sorted({f for g in args.glob for f in globmod.glob(g)})
    if not files:
        sys.exit(f"no snapshots matched {args.glob}")
    chains = sorted({re.sub(r"_snap\d+\.mfd$", "", os.path.basename(f))
                     for f in files})

    spec = Counter()
    size_of, star_of, chg_of = (defaultdict(list) for _ in range(3))
    rel_of, sumz_of = defaultdict(list), defaultdict(list)
    shape_of = {}
    n6_of = defaultdict(Counter)
    sig_of = defaultdict(Counter)
    nodes_of = defaultdict(Counter)
    is_flick = {}
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
                vc = ec = None           # --group induced carries no colours
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
                if args.flicker and key not in is_flick:
                    is_flick[key] = cf.flicker_key(
                        args.group, facets=fac, vcolor=vc,
                        ecolor=ec) in flick
            else:
                key = cx.key
                if args.flicker and key not in is_flick:
                    is_flick[key] = cf.flicker_key("sig", cx=cx) in flick
            spec[key] += 1
            size_of[key].append(len(cx.verts))
            star_of[key].append(len(st.star(cx.verts)))
            chg_of[key].append(st.complex_charge(cx.verts, q))
            rel_of[key].append(st.complex_charge_rel(cx.verts, q, qbar))
            sumz_of[key].append(st.total_coordination(cx.verts))
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
        fcol = "F " if args.flicker else ""
        hdr = (f"{'#':<4} {fcol}{'f=(v,e,t,T)':<17} {'facets':<6} "
               f"{'degseq (1-skel)':<20} {'n':>5} {'%':>6} {'sumZ':>6} "
               f"{'Q_c raw':>16} {'Q_c-bg':>9}")
        print(hdr)
        print("-" * len(hdr))
        rows = []
        for i, (k, n) in enumerate(spec.most_common()):
            if n < args.min_count or i >= args.top:
                break
            sh = shape_of[k]
            dq = str(sh["degseq"])
            nf = len(k[0]) if args.group in ("decorated", "full") else len(k)
            rows.append({"class": i, "flicker": is_flick.get(k),
                         "f": list(sh["f"]), "n_facets": nf,
                         "degseq": list(sh["degseq"]), "n": n,
                         "frac": n / ncomp,
                         "sigs": {str(a): b for a, b in sig_of[k].items()},
                         "nodes": {str(a): b for a, b in nodes_of[k].items()},
                         "star": float(np.mean(star_of[k])),
                         "sum_Z": float(np.mean(sumz_of[k])),
                         "Q_c_raw": float(np.mean(chg_of[k])),
                         "Q_c_raw_sd": float(np.std(chg_of[k])),
                         "Q_c_rel": float(np.mean(rel_of[k])),
                         "Q_c_rel_sd": float(np.std(rel_of[k]))})
            sz = np.mean(sumz_of[k])
            fm = ("F " if is_flick.get(k) else ". ") if args.flicker else ""
            print(f"{i:<4} {fm}{str(sh['f']):<17} {nf:<6} {dq[:19]:<20} "
                  f"{n:5d} {100*n/ncomp:5.1f}% {sz:6.1f} "
                  f"{np.mean(chg_of[k]):8.3f}+-{np.std(chg_of[k]):<7.2g} "
                  f"{np.mean(rel_of[k]):9.3f}")
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
        fcol = "F " if args.flicker else ""
        hdr = (f"{fcol}{'illegal-edge signature':<20} {'nodes':<9} {'n':>5} "
               f"{'%':>6} {'s (verts)':>13} {'star':>7} {'sum w':>6} "
               f"{'Q_c (rad)':>15}")
        print(hdr)
        print("-" * len(hdr))
        rows = []
        for k, n in spec.most_common():
            if n < args.min_count or len(rows) >= args.top:
                break
            sig, nodes = k
            w = sum(6 - d for d in sig)
            rows.append({"sig": list(sig), "nodes": list(nodes), "n": n,
                         "flicker": is_flick.get(k),
                         "frac": n / ncomp, "size": float(np.mean(size_of[k])),
                         "star": float(np.mean(star_of[k])), "sum_w": w,
                         "sum_Z": float(np.mean(sumz_of[k])),
                         "Q_c_raw": float(np.mean(chg_of[k])),
                         "Q_c_raw_sd": float(np.std(chg_of[k])),
                         "Q_c_rel": float(np.mean(rel_of[k])),
                         "Q_c_rel_sd": float(np.std(rel_of[k]))})
            fm = ("F " if is_flick.get(k) else ". ") if args.flicker else ""
            print(f"{fm}{sig_str(sig):<20} {str(nodes) if nodes else '':<9} {n:5d} "
                  f"{100*n/ncomp:5.1f}% "
                  f"{np.mean(size_of[k]):6.1f}+-{np.std(size_of[k]):<5.1f} "
                  f"{np.mean(star_of[k]):7.0f} {w:6d} "
                  f"{np.mean(chg_of[k]):7.2f}+-{np.std(chg_of[k]):<6.2f}")

    if args.flicker:
        fl_n = sum(n for k, n in spec.items() if is_flick.get(k))
        fl_s = sum(1 for k in spec if is_flick.get(k))
        print()
        print("=" * 78)
        print("FLICKER BACKGROUND  (F in the table above)")
        print("=" * 78)
        print(f"reference    : {flick_meta['cell']}")
        print(f"               {flick_meta['sites']} single-2->3 sites, "
              f"{flick_meta['species']} distinct species, N_3 = "
              f"{flick_meta['n3']}")
        print(f"grouping     : --group {args.group}")
        print(f"reachable    : {fl_s} of {len(spec)} species, "
              f"{fl_n} of {ncomp} complexes ({100*fl_n/ncomp:.1f}%)")
        print(f"NOT reachable: {len(spec)-fl_s} species, {ncomp-fl_n} "
              f"complexes ({100*(ncomp-fl_n)/ncomp:.1f}%)   <- the residual "
              f"census")
        print()
        print("  READ THIS AS AN UPPER BOUND on flicker contamination. A")
        print("  species marked F is one a single 2->3 on the PRISTINE crystal")
        print("  can produce; a snapshot carries no move history, so this")
        print("  cannot tell you that any particular instance arose that way.")
        print("  A long-lived complex that happens to be flicker-shaped is")
        print("  counted here too -- measured at lam=0.40, 49 of the 57 tracks")
        print("  living past 20 sweeps were flicker-born, so F does NOT mean")
        print("  transient. For a real per-instance attribution, which needs")
        print("  the birth move, run flicker_fraction.py on an event stream.")
        if args.group != "full":
            print()
            print(f"  NOTE --group {args.group} is coarser than the key the")
            print("  flicker list is most discriminating at. A coarser key")
            print("  merges species, so it can only OVERSTATE reachability.")
        print()
        print("residual census after removing flicker-reachable species:")
        res = [(k, n) for k, n in spec.most_common() if not is_flick.get(k)]
        if not res:
            print("   (nothing -- every observed species is flicker-reachable)")
        for i, (k, n) in enumerate(res[:args.top]):
            if args.group in ("induced", "decorated", "full"):
                sh = shape_of[k]
                print(f"   {i:<3} f={str(sh['f']):<17} n={n:5d} "
                      f"{100*n/ncomp:5.1f}% of all   "
                      f"{100*n/(ncomp-fl_n):5.1f}% of residual   "
                      f"sumZ {np.mean(sumz_of[k]):6.1f}  "
                      f"Q_c {np.mean(chg_of[k]):7.3f}")
            else:
                print(f"   {i:<3} {sig_str(k[0]):<20} n={n:5d} "
                      f"{100*n/ncomp:5.1f}% of all   "
                      f"Q_c {np.mean(chg_of[k]):7.3f}")
        if len(res) > args.top:
            print(f"   ... {len(res)-args.top} further residual species")

    kept = [(k, n) for k, n in spec.items() if n >= args.min_count]
    hid = [(k, n) for k, n in spec.items() if n < args.min_count]
    if hid:
        print()
        print(f"below --min-count {args.min_count}: {len(hid)} species holding "
              f"{sum(n for _, n in hid)} complexes "
              f"({100 * sum(n for _, n in hid) / ncomp:.1f}% of the "
              f"population) -- not shown, not dropped from the totals above")
    if len(kept) > args.top:
        print(f"beyond --top {args.top}: {len(kept) - args.top} further "
              f"species above the count cut are not shown")

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
                       "flicker": (None if not args.flicker else {
                           "reference": flick_meta,
                           "group": args.group,
                           "species_reachable": sum(
                               1 for k in spec if is_flick.get(k)),
                           "complexes_reachable": sum(
                               n for k, n in spec.items() if is_flick.get(k)),
                           "note": "UPPER BOUND -- reachable by a single 2->3 "
                                   "on the pristine crystal, not a per-instance "
                                   "attribution"}),
                       "rows": rows}, fh, indent=1)
        print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
