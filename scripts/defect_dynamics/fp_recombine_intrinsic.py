#!/usr/bin/env python3
"""P(recombine) by INTRINSIC dock class (development frame).

Reruns the recombination-by-angle question with exact intrinsic labels:
each freed dock in prodB_*/prodB2_* gets its exact rational cos^2 class
(canonical shortest connecting path, as in fp_dock_census_intrinsic),
and outcomes are tabulated per class -- plus the coarse aligned
(cos^2 > 1/2) vs crossed split the registry data only hinted at.
"""
import glob
import json
import os
import sys
from collections import Counter
from fractions import Fraction

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
for _p in ("../../python", "../../scripts", "../../tools", "."):
    sys.path.insert(0, os.path.join(_HERE, _p))
import discrete_differential_geometry as ddg
from discrete_differential_geometry.fpkmc import face_apex_maps
from worm_helix import bc_orbit
from fp_dock_census_intrinsic import (chain_seq, stack_tets, connector,
                                      intrinsic_cos2)

ROOT = os.path.join(_HERE, "..", "..")
DATA = os.path.join(ROOT, "data", "fpkmc")
REF = os.path.join(ROOT, "data", "tcp_reference", "T3_R_m2_N7248.mfd")


def main():
    m = ddg.Manifold.load(REF, 3)
    F = np.asarray(m.facets())
    orb = [int(x) for x in bc_orbit(m, [int(x) for x in F[0]])]
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

    rows = []
    fails = 0
    files = sorted(glob.glob(os.path.join(DATA, "prodB_p*.json"))) + \
        sorted(glob.glob(os.path.join(DATA, "prodB2_p*.json")))
    for f in files:
        d = json.load(open(f))
        for e in d["episodes"]:
            if not (e.get("dock_chord") and e.get("freed")
                    and "outcome" in e):
                continue
            try:
                cA = tuple(sorted(int(x) for x in e["dock_chord"]))
                cB = tuple(sorted((orb[e["jB"] % L],
                                   orb[(e["jB"] + 4) % L])))
                c2, plen = intrinsic_cos2(m, dual, cA, list(face_of[cA]),
                                          cB, list(face_of[cB]))
                rows.append({"cos2": c2, "outcome": e["outcome"],
                             "sep": d["sep_sites"], "path": plen})
            except Exception as ex:
                fails += 1
                print(f"  fail: {type(ex).__name__}: {ex}")
    print(f"{len(rows)} freed docks labeled, {fails} failures\n")

    stats = {}
    for r in rows:
        stats.setdefault(r["cos2"], Counter())[r["outcome"]] += 1
    print(f"{'cos^2 (exact)':>20} {'angle':>8} {'rec':>4} {'esc':>4} "
          f"{'other':>5} {'P(rec)':>7}")
    for c2, cnt in sorted(stats.items(),
                          key=lambda kv: -sum(kv[1].values())):
        nr, ne = cnt.get("recombine", 0), cnt.get("escape", 0)
        no = sum(cnt.values()) - nr - ne
        tot = nr + ne
        ang = float(np.degrees(np.arccos(np.sqrt(float(c2)))))
        p = f"{nr / tot:.2f}" if tot else "--"
        print(f"{str(c2):>20} {ang:7.2f}d {nr:4d} {ne:4d} {no:5d} {p:>7}")

    for label, sel in (("aligned (cos^2 > 1/2)",
                        lambda c: c > Fraction(1, 2)),
                       ("crossed (cos^2 <= 1/2)",
                        lambda c: c <= Fraction(1, 2))):
        nr = sum(c.get("recombine", 0) for c2, c in stats.items()
                 if sel(c2))
        ne = sum(c.get("escape", 0) for c2, c in stats.items() if sel(c2))
        tot = nr + ne
        if tot:
            p = nr / tot
            print(f"\n{label}: P(rec) = {nr}/{tot} = {p:.2f} "
                  f"+- {np.sqrt(p * (1 - p) / tot):.2f}")
    with open(os.path.join(DATA, "fp_recombine_intrinsic.json"),
              "w") as fh:
        json.dump([{"cos2": str(r["cos2"]), "outcome": r["outcome"],
                    "sep": r["sep"], "path": r["path"]} for r in rows],
                  fh)


if __name__ == "__main__":
    main()
