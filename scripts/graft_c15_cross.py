#!/usr/bin/env python3
"""Cross-crystal grafting with cubic C15: (111) slabs vs hex basal slabs.

C15 (MgCu2) stacks its Laves layers along the cubic [111]; commensuration of
the boundary tori requires the (111) wrap of an m-cell cubic box (a 2m x 2m
hex torus of the Laves net) to equal the basal wrap of an m_h-cell hex box
(m_h x m_h), i.e. m_h = 2 m.  So c15 m=2 is the epitaxial partner of
c14/c36 m=4.

Census c15 m2 along (1,1,1) + c14/c36 m4 basal, join the decorated boundary
certificates at levels 1/2/3, then perform and validate cross grafts
(c15 slab into the c36 host, and into the c14 host).  Reports the graft
vertex-count change dV = donor interior vertices - host-lump interior
vertices: level-3 seams are FK-perfect by construction, so any dV != 0
match is a defect-free vertex-count-changing graft.

Usage:
    caffeinate -i python scripts/graft_c15_cross.py [--out data/grafts]
"""
import argparse
import json
import os
import sys
import time

import numpy as np

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_ROOT, "python"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from graft_c36_control import slab_census, group_report, cert_hash
from graft_signature import match_phis, graft, validate_facets


def cross_join(slabs_a, slabs_b):
    """{level: [(key_a, key_b)]} exact-cert matches (all pairs)."""
    out = {}
    for lvl in (3, 2, 1):
        by_cert = {}
        for key, (sig, _) in slabs_a.items():
            by_cert.setdefault(sig["certs"][lvl], []).append(key)
        pairs = []
        for key_b, (sig_b, _) in slabs_b.items():
            for key_a in by_cert.get(sig_b["certs"][lvl], []):
                pairs.append((key_a, key_b))
        out[lvl] = pairs
    return out


def do_graft(tag, ctx_host, slabs_host, ctx_donor, slabs_donor, pairs, out):
    """Graft the most interior-different donor slab into the host; validate."""
    def pref(pair):
        kd, kh = pair
        fpd = slabs_donor[kd][0]["diag"]["interior_fp"]
        fph = slabs_host[kh][0]["diag"]["interior_fp"]
        dv = abs(slabs_donor[kd][0]["diag"]["n_int_vertices"]
                 - slabs_host[kh][0]["diag"]["n_int_vertices"])
        return (fpd != fph, dv)
    kd, kh = max(pairs, key=pref)
    sigd, lumpd = slabs_donor[kd]
    sigh, lumph = slabs_host[kh]
    dv = (sigd["diag"]["n_int_vertices"] - sigh["diag"]["n_int_vertices"])
    print(f"\n[{tag}] graft donor {kd} ({sigd['diag']['n_tets']} tets) "
          f"into host {kh} ({sigh['diag']['n_tets']} tets), dV={dv:+d}:")
    for i, phi in enumerate(match_phis(sigd, sigh, level=3)):
        newF, info = graft(ctx_host, lumph, ctx_donor, lumpd, phi)
        rep = validate_facets(newF)
        good = (rep["all_56"] and rep["n_broken_disclination"] == 0
                and rep["euler_characteristic"] == 0 and rep["orientable"])
        print(f"  phi#{i}: all56={rep['all_56']} "
              f"broken={rep['n_broken_disclination']} "
              f"chi={rep['euler_characteristic']} orient={rep['orientable']} "
              f"V={rep['n_vertices']} Z={rep['z_census']}")
        if good:
            from discrete_differential_geometry import Manifold
            path = os.path.join(out, f"T3_{tag}_f3{len(newF)}.mfd")
            Manifold(3, newF.tolist()).save(path, comments=[
                f"graft: donor {kd} into host {kh} (level-3 match, phi #{i})",
                f"dV = {dv:+d} interior vertices", f"validation: {rep}"])
            print(f"  SAVED -> {path}")
            return dict(donor=list(kd), host=list(kh), dv=dv, phi=i,
                        validation={k: str(v) for k, v in rep.items()},
                        path=path)
    print(f"  [{tag}] no orientation-compatible phi validated!")
    return None


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", default=os.path.join(_ROOT, "data", "grafts"))
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    ctx15, slabs15 = slab_census("c15", 2, 3, 6, normal=(1, 1, 1))
    ctx14, slabs14 = slab_census("c14", 4, 3, 12)
    ctx36, slabs36 = slab_census("c36", 4, 3, 12)

    group_report(slabs15, 3, "within-c15(111)")

    summary = {"n_slabs": {"c15m2": len(slabs15), "c14m4": len(slabs14),
                           "c36m4": len(slabs36)}}
    for tag, (ctxh, slabsh) in (("c15xc14", (ctx14, slabs14)),
                                ("c15xc36", (ctx36, slabs36))):
        cj = cross_join(slabs15, slabsh)
        print(f"\ncross join {tag}: " +
              ", ".join(f"L{l}: {len(p)} pairs" for l, p in cj.items()))
        summary[tag] = {f"L{l}": [[list(a), list(b)] for a, b in p]
                        for l, p in cj.items()}
        if cj[3]:
            res = do_graft(f"{tag.split('x')[1].upper()}m4_graftC15m2",
                           ctxh, slabsh, ctx15, slabs15, cj[3], args.out)
            if res:
                summary[f"graft_{tag}"] = res

    jpath = os.path.join(args.out, "c15_cross.json")
    with open(jpath, "w") as fh:
        json.dump(summary, fh, indent=1, default=str)
    print(f"\nsummary JSON -> {jpath}")


if __name__ == "__main__":
    main()
