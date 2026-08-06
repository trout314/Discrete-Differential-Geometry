#!/usr/bin/env python3
"""C36 positive control for crystal grafting (see graft_signature.py).

C36 (MgNi2) is the hc-stacked Laves polytype: it interleaves C14-like (h) and
C15-like (c) stacking of the same Laves layer, so basal slabs of C36 should
share decorated boundary surfaces with basal slabs of pure C14 -- a
guaranteed-positive test of the whole decorated-surface matching pipeline.

Pipeline:
  1. Build c36 and c14 references in-memory at the same in-plane box size
     (m x m x m cells) with their deterministic vertex positions.
  2. Enumerate basal slabs: every cut sits between consecutive atomic
     z-planes; a slab is the set of tets whose (unwrap-consistent) centroid
     lies between two cuts spanning kmin..kmax planes.  Base cut restricted
     to the first cell (z-translations are exact symmetries).
  3. Exact decorated-boundary certificates at levels 1/2/3.
  4. Join: within-crystal groups (boundary-isomorphic slabs with different
     interiors = stacking-swap candidates) and the cross-crystal match table.
  5. Perform the best cross-crystal graft (c14 slab into the c36 host),
     validate (all edges {5,6}, no broken disclination lines, chi=0,
     orientable), and save the grafted manifold.

Usage:
    caffeinate -i python scripts/graft_c36_control.py [--m 4]
        [--kmin 3] [--kmax 12] [--out data/grafts]
"""
import argparse
import hashlib
import json
import os
import sys
import time
from collections import defaultdict

import numpy as np

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_ROOT, "python"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from tcp_reference import build_t3_triangulation
from cocycle_check import reference_frac_positions
from graft_signature import (CrystalContext, SurfaceError, lump_signature,
                             match_phis, graft, validate_facets)

EPS = 3.3e-5  # nudge cuts off exact centroid ties


def z_cuts(uvals, m, gap=0.02):
    """Cut positions between consecutive atomic layer-planes (period m)."""
    z = np.sort(uvals % m)
    brk = np.nonzero(np.diff(z) > gap)[0]
    centers = []
    start = 0
    for b in list(brk) + [len(z) - 1]:
        centers.append(z[start:b + 1].mean())
        start = b + 1
    centers = np.array(centers)
    # wrap-merge: first and last plane may be the same plane split at 0
    if (centers[0] + m) - centers[-1] < gap:
        centers = centers[:-1]
    nxt = np.roll(centers, -1)
    nxt[-1] += m
    return ((centers + nxt) / 2 + EPS) % m


def tet_layer_coord(F, uvals, m):
    """Unwrap-consistent centroid layer coordinate (mod m) per tet."""
    z = uvals % m
    zt = z[F]                                   # (nt, 4)
    z0 = zt[:, :1]
    d = (zt - z0 + m / 2) % m - m / 2
    return (z0[:, 0] + d.mean(axis=1)) % m


def slab_census(name, m, kmin, kmax, normal=(0, 0, 1)):
    """Signatures of all layer slabs of one crystal, cut along `normal`.

    normal is an integer plane normal in fractional coordinates; the layer
    coordinate u = frac . normal is well-defined mod m on the torus, and
    unit-cell translations shift u by integers, so base cuts are deduped to
    u0 in [0,1).  (0,0,1) = basal for the hex phases; (1,1,1) = the Laves
    stacking direction of cubic c15.
    """
    t0 = time.time()
    fac, nv = build_t3_triangulation(name, m)
    ctx = CrystalContext(fac, f"{name}m{m}")
    pos = reference_frac_positions(name, m)
    uvals = pos @ np.asarray(normal, float)
    cuts = np.sort(z_cuts(uvals, m))
    zc = tet_layer_coord(ctx.F, uvals, m)
    print(f"[{name} m{m} n={normal}] V={nv} tets={len(fac)} "
          f"cuts/box={len(cuts)} ({time.time()-t0:.1f}s build)")

    slabs = {}
    skipped = 0
    for i0, u0 in enumerate(cuts):
        if u0 >= 1.0:      # base cut in the first cell only
            continue
        for k in range(kmin, kmax + 1):
            i1 = (i0 + k) % len(cuts)
            u1 = cuts[i1]
            th = (u1 - u0) % m
            lump = np.nonzero(((zc - u0) % m) < th)[0]
            if len(lump) == 0 or len(lump) == len(ctx.F):
                continue
            try:
                sig = lump_signature(ctx, lump)
            except SurfaceError:
                skipped += 1
                continue
            slabs[(name, round(float(u0), 4), round(float(u1), 4), k)] = \
                (sig, lump)
    print(f"[{name} m{m}] {len(slabs)} slab signatures "
          f"({skipped} non-manifold cuts skipped, {time.time()-t0:.1f}s)")
    return ctx, slabs


def cert_hash(cert):
    return hashlib.blake2b(repr(cert).encode(), digest_size=8).hexdigest()


def group_report(slabs, level, label):
    """Group slabs by exact level-cert; report groups with >1 interior."""
    groups = defaultdict(list)
    for key, (sig, lump) in slabs.items():
        groups[sig["certs"][level]].append(key)
    multi = {c: ks for c, ks in groups.items() if len(ks) > 1}
    print(f"\n{label}: {len(groups)} distinct L{level} boundary classes "
          f"among {len(slabs)} slabs; {len(multi)} classes with >1 member")
    interesting = []
    for c, ks in sorted(multi.items(), key=lambda kv: -len(kv[1])):
        fps = {slabs[k][0]["diag"]["interior_fp"] for k in ks}
        if len(fps) > 1:
            interesting.append((c, ks, fps))
    print(f"  classes whose members have DIFFERENT interiors: "
          f"{len(interesting)}")
    for c, ks, fps in interesting[:6]:
        print(f"    {cert_hash(c)}: {len(ks)} slabs, {len(fps)} interiors "
              f"-> {[k[1:] for k in ks[:6]]}")
    return groups


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--m", type=int, default=4)
    ap.add_argument("--kmin", type=int, default=3)
    ap.add_argument("--kmax", type=int, default=12)
    ap.add_argument("--out", default=os.path.join(_ROOT, "data", "grafts"))
    ap.add_argument("--no-graft", action="store_true")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    ctx36, slabs36 = slab_census("c36", args.m, args.kmin, args.kmax)
    ctx14, slabs14 = slab_census("c14", args.m, args.kmin, args.kmax)

    g36 = group_report(slabs36, 3, "within-c36")
    g14 = group_report(slabs14, 3, "within-c14")

    # ---- cross-crystal join ------------------------------------------------
    print("\ncross-crystal join (c36 x c14):")
    cross = {}
    for lvl in (3, 2, 1):
        c36certs = {sig["certs"][lvl]: key
                    for key, (sig, _) in slabs36.items()}
        pairs = []
        for key14, (sig14, _) in slabs14.items():
            c = sig14["certs"][lvl]
            if c in c36certs:
                pairs.append((c36certs[c], key14))
        cross[lvl] = pairs
        print(f"  L{lvl}: {len(pairs)} matched slab pairs")

    summary = {
        "m": args.m, "kmin": args.kmin, "kmax": args.kmax,
        "n_slabs": {"c36": len(slabs36), "c14": len(slabs14)},
        "n_classes_L3": {"c36": len(g36), "c14": len(g14)},
        "cross_matches": {lvl: [[list(a), list(b)] for a, b in ps]
                          for lvl, ps in cross.items()},
    }

    # ---- perform the best cross-crystal graft ------------------------------
    if cross[3] and not args.no_graft:
        # prefer a pair whose interiors differ most (a genuine heterograft):
        def pref(pair):
            k36, k14 = pair
            fp36 = slabs36[k36][0]["diag"]["interior_fp"]
            fp14 = slabs14[k14][0]["diag"]["interior_fp"]
            return (fp36 != fp14, abs(fp36[0] - fp14[0]))
        k36, k14 = max(cross[3], key=pref)
        sig36, lump36 = slabs36[k36]
        sig14, lump14 = slabs14[k14]
        print(f"\ngrafting c14 slab {k14[1:]} into c36 host at {k36[1:]}:")
        print(f"  host lump: {sig36['diag']['n_tets']} tets, "
              f"donor: {sig14['diag']['n_tets']} tets")
        ok = None
        for i, phi in enumerate(match_phis(sig14, sig36, level=3)):
            newF, info = graft(ctx36, lump36, ctx14, lump14, phi)
            rep = validate_facets(newF)
            good = (rep["all_56"] and rep["n_broken_disclination"] == 0
                    and rep["euler_characteristic"] == 0 and rep["orientable"])
            print(f"  phi#{i}: all56={rep['all_56']} "
                  f"broken={rep['n_broken_disclination']} "
                  f"chi={rep['euler_characteristic']} "
                  f"orient={rep['orientable']} Z={rep['z_census']}")
            if good:
                ok = (newF, info, rep, i)
                break
        if ok:
            newF, info, rep, i = ok
            from discrete_differential_geometry import Manifold
            path = os.path.join(
                args.out, f"T3_C36m{args.m}_graftC14_"
                          f"{k36[1]}-{k36[2]}_f3{len(newF)}.mfd")
            Manifold(3, newF.tolist()).save(path, comments=[
                f"graft: c14 m{args.m} slab z=({k14[1]},{k14[2]}) planes={k14[3]}",
                f"  into c36 m{args.m} host, removed slab z=({k36[1]},{k36[2]})"
                f" planes={k36[3]}",
                f"level-3 decorated-boundary match, phi #{i}",
                f"validation: {rep}",
            ])
            print(f"  SAVED graft -> {path}")
            summary["graft"] = {"host": list(k36), "donor": list(k14),
                                "info": info, "validation":
                                {k: v for k, v in rep.items()},
                                "path": path}
    elif not cross[3]:
        print("\nno level-3 cross matches -- no graft performed")

    jpath = os.path.join(args.out, f"c36_control_m{args.m}.json")
    with open(jpath, "w") as fh:
        json.dump(summary, fh, indent=1, default=str)
    print(f"\nsummary JSON -> {jpath}")


if __name__ == "__main__":
    main()
