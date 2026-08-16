#!/usr/bin/env python3
"""Re-label the phase-2 crossing-collider table with exact intrinsic
contact classes (CONVENTIONS.md section 6).

Loads the deterministic collider reruns (crossing_default.json,
crossing_steep.json -- results now record chordA/chordB), computes each
crossing's exact rational cos^2 class in the development frame
(canonical shortest connecting path), and prints the table: registry
angle vs intrinsic class, alongside dmin, handedness, contact outcome
and V. The physics columns (outcome, V, products) are unchanged --
this replaces only their geometric labels.
"""
import json
import os
import sys
from fractions import Fraction

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
for _p in ("../../python", "../../scripts", "../../tools",
           "../../scripts/defect_dynamics", "."):   # archived: siblings stayed in dd/
    sys.path.insert(0, os.path.join(_HERE, _p))
import discrete_differential_geometry as ddg
from discrete_differential_geometry.fpkmc import face_apex_maps
from fp_dock_census_intrinsic import intrinsic_cos2

ROOT = os.path.join(_HERE, "..", "..")
DATA = os.path.join(ROOT, "data", "fpkmc")
# the crossing collider runs on the m4 host
REF = os.path.join(ROOT, "data", "tcp_reference", "T3_R_m4_N57984.mfd")


def main():
    m = ddg.Manifold.load(REF, 3)
    F = np.asarray(m.facets())
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
    for name in ("crossing_default", "crossing_steep"):
        p = os.path.join(DATA, name + ".json")
        if not os.path.exists(p):
            continue
        d = json.load(open(p))
        for r in d["results"]:
            cA = tuple(sorted(r["chordA"]))
            cB = tuple(sorted(r["chordB"]))
            try:
                c2, plen = intrinsic_cos2(m, dual, cA, list(face_of[cA]),
                                          cB, list(face_of[cB]))
            except Exception as ex:
                print(f"  {name}#{r['crossing']}: FAILED "
                      f"{type(ex).__name__}: {ex}")
                continue
            o = r["outcome"]
            rows.append({
                "run": name.split("_")[1], "i": r["crossing"],
                "reg_angle": r["angle"], "cos2": c2, "path": plen,
                "dmin": r["dmin"], "hand": r.get("hand"),
                "shared": r["shared"], "type": o["type"],
                "V": o.get("V"),
                "f": tuple(o["species_f"]) if o.get("species_f") else None})

    print(f"{'run':>8} {'reg ang':>8} {'intrinsic cos^2':>18} {'ang':>7} "
          f"{'path':>4} {'hand':>5} {'shr':>3} {'outcome':>12} {'V':>7}")
    for r in sorted(rows, key=lambda x: (x["run"], x["i"])):
        ang = float(np.degrees(np.arccos(np.sqrt(float(r["cos2"])))))
        v = f"{r['V']:+.3f}" if r["V"] is not None else "--"
        print(f"{r['run']:>8} {r['reg_angle']:7.1f}d {str(r['cos2']):>18} "
              f"{ang:6.2f}d {r['path']:4d} {r['hand']:+5.1f} "
              f"{r['shared']:3d} {r['type']:>12} {v:>7}"
              + (f"  f={r['f']}" if r["f"] else ""))

    with open(os.path.join(DATA, "crossing_intrinsic_labels.json"),
              "w") as fh:
        json.dump([{**r, "cos2": str(r["cos2"])} for r in rows], fh)
    print(f"\nwrote {os.path.join(DATA, 'crossing_intrinsic_labels.json')}")


if __name__ == "__main__":
    main()
