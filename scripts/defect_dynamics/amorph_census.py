#!/usr/bin/env python3
"""Defect census for AMORPHOUS (non-dilute) states -- the fk_amorphous anneal.

The dilute-gas instruments do not transfer here. At f_FK ~ 60% roughly a third
of vertices touch an illegal edge, so the illegal-edge graph percolates and a
species catalog degenerates to "one complex". What is informative instead is
the composition:

  EDGE side   the degree histogram. Legality is deg in {5,6}; everything else
              is the defect population, and its SIGN matters -- deg 3/4 are
              positive-curvature (the deg-4 quantum of the FK gas), deg 7+
              negative.

  VERTEX side n6(v) = # incident deg-6 edges. On a fully legal vertex the link
              sum rule gives Z = 12 + n6 exactly, so n6 IS the Frank-Kasper
              class: 0 -> Z12, 2 -> Z14, 3 -> Z15, 4 -> Z16. n6 = 1 is
              impossible (no fullerene has exactly one hexagon) and n6 >= 5 is
              a HUB -- edge-legal but not Frank-Kasper, and the known endpoint
              of every previous anneal. The census therefore splits the legal
              vertices into the FK classes and the hubs, which is the quantity
              zleg would price if it were switched on.

  INVARIANT   n6 + 2 n7 - n4 - 2 n3 = 6 f3 - 5 f1 = f3 - 5 f0, fixed by the
              two pins alone. Defects can only REDISTRIBUTE the six-web
              excess, never change its total, so the census prints it as a
              consistency check on both the state and the pins.

Usage:
    python scripts/defect_dynamics/amorph_census.py data/fk_amorph/ratchet/*.final.mfd \
        --ref data/tcp_reference/T3_A15_m5_N5750.mfd
"""
import argparse
import os
import sys
from collections import Counter

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(os.path.dirname(_HERE))
sys.path.insert(0, os.path.join(_ROOT, "python"))
sys.path.insert(0, os.path.join(_ROOT, "scripts"))
sys.path.insert(0, _HERE)

import discrete_differential_geometry as ddg
import defect_state as dsm

FK = {0: "Z12", 2: "Z14", 3: "Z15", 4: "Z16"}


def census(path):
    m = ddg.Manifold.load(path, 3)
    st = dsm.DefectState(m)
    f0, f1, f2, f3 = (int(x) for x in m.f_vector)
    V = list(st.v2t)
    n = len(V)

    edeg = Counter(st.edeg.values())
    n_ill_e = sum(c for d, c in edeg.items() if d < 5 or d > 6)
    legal_v = [v for v in V if st.imp[v] == 0]
    hubs = [v for v in legal_v if st.n6[v] >= 5]
    fk_cls = Counter(st.n6[v] for v in legal_v if st.n6[v] in FK)
    n6_one = sum(1 for v in legal_v if st.n6[v] == 1)

    # illegal-edge graph components (the "narrow" definition)
    parent, e2i = {}, {}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    ill = sorted(st.ill_edges)
    for i, e in enumerate(ill):
        parent[i] = i
        for v in e:
            if v in e2i:
                a, b = find(i), find(e2i[v])
                if a != b:
                    parent[a] = b
            e2i[v] = i
    groups = {}
    for i, e in enumerate(ill):
        groups.setdefault(find(i), set()).update(e)
    sizes = sorted((len(g) for g in groups.values()), reverse=True)

    inv = (edeg[6] + 2 * sum(c for d, c in edeg.items() if d == 7)
           - edeg[4] - 2 * edeg[3]
           + sum((d - 5) * c for d, c in edeg.items() if d >= 8))
    return dict(
        path=path, f=(f0, f1, f2, f3), n=n,
        ebar=6.0 * f3 / f1, zbar=2.0 * f1 / f0,
        f_FK=(n - sum(1 for v in V if st.imp[v] > 0)) / n,
        f_e=(edeg[5] + edeg[6]) / sum(edeg.values()),
        edeg=edeg, n_ill_e=n_ill_e,
        fk_cls=fk_cls, hubs=len(hubs), n6_one=n6_one,
        hub_n6=Counter(st.n6[v] for v in hubs),
        ncomp=len(sizes), sizes=sizes,
        n_ill_v=sum(1 for v in V if st.imp[v] > 0),
        inv=inv, inv_pred=f3 - 5 * f0,
        maxm=max((st.imp[v] for v in V), default=0))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("states", nargs="+")
    ap.add_argument("--ref", nargs="*", default=[],
                    help="reference .mfd files to census alongside")
    args = ap.parse_args()

    rows = [census(p) for p in list(args.ref) + list(args.states)
            if os.path.exists(p)]

    print(f"{'state':<26}{'f0':>6}{'f3':>7}{'f_FK%':>8}{'f_e%':>7}{'ill_e':>7}"
          f"{'Zbar':>9}{'ebar':>11}{'hubs':>6}{'n6=1':>6}")
    for r in rows:
        print(f"{os.path.basename(r['path'])[:25]:<26}{r['f'][0]:>6}{r['f'][3]:>7}"
              f"{100*r['f_FK']:>8.2f}{100*r['f_e']:>7.2f}{r['n_ill_e']:>7}"
              f"{r['zbar']:>9.4f}{r['ebar']:>11.7f}{r['hubs']:>6}{r['n6_one']:>6}")

    print("\nEDGE-DEGREE HISTOGRAM (legality is {5,6}; 3/4 positive curvature, 7+ negative)")
    degs = sorted({d for r in rows for d in r["edeg"]})
    print(f"{'state':<26}" + "".join(f"{('deg'+str(d)):>9}" for d in degs))
    for r in rows:
        print(f"{os.path.basename(r['path'])[:25]:<26}"
              + "".join(f"{r['edeg'].get(d,0):>9}" for d in degs))

    print("\nFRANK-KASPER COMPOSITION of the fully-legal vertices "
          "(n6 -> Z; hubs are n6>=5, legal but NOT FK)")
    print(f"{'state':<26}{'legal v':>9}{'Z12':>8}{'Z14':>8}{'Z15':>8}{'Z16':>8}"
          f"{'hubs':>7}{'hub n6':>16}")
    for r in rows:
        legal = sum(r["fk_cls"].values()) + r["hubs"] + r["n6_one"]
        hub = ",".join(f"{k}:{v}" for k, v in sorted(r["hub_n6"].items())[:4])
        print(f"{os.path.basename(r['path'])[:25]:<26}{legal:>9}"
              f"{r['fk_cls'].get(0,0):>8}{r['fk_cls'].get(2,0):>8}"
              f"{r['fk_cls'].get(3,0):>8}{r['fk_cls'].get(4,0):>8}"
              f"{r['hubs']:>7}{hub:>16}")

    print("\nILLEGAL-EDGE GRAPH (percolates once f_FK is low -- read frac_top1, not species)")
    print(f"{'state':<26}{'defect v':>9}{'ncomp':>7}{'top1':>7}{'frac_top1':>11}"
          f"{'max m':>7}{'invariant':>11}{'= f3-5f0':>10}")
    for r in rows:
        t1 = r["sizes"][0] if r["sizes"] else 0
        print(f"{os.path.basename(r['path'])[:25]:<26}{r['n_ill_v']:>9}{r['ncomp']:>7}"
              f"{t1:>7}{(t1/max(r['n_ill_v'],1)):>11.3f}{r['maxm']:>7}"
              f"{r['inv']:>11}{r['inv_pred']:>10}")


if __name__ == "__main__":
    main()
