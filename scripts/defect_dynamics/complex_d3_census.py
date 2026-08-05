#!/usr/bin/env python3
"""Do NATURALLY OCCURRING defect complexes hold a deg-3 edge?

This is the ceiling on any contact-handoff rule: the lift can only be handed
to an existing deg-3 chord, so if a fused complex has none, a flight that
collides with it has nothing to hand to.

Why this script exists (2026-08-05). A census run on the ecmc_ab BLOB said
62-94% of blocking complexes hold a deg-3 edge, rising with size. That was an
artifact: `build_blob` injects 12 fliers into a radius-4 ball, each flier IS a
deg-3 edge, and they agglomerate -- so a large "blocker" there is a clump of
fliers, not a fused defect. Aaron caught it from the defect catalog, where the
only deg-3 edges in a thermal snapshot were the elementary (3,4,4) flickers
and every larger complex was all-deg-4.

So: measure it on thermal gas snapshots instead, where the complexes formed on
their own.

Two adjacency conventions are reported because they can disagree:
  illegal   components of the ILLEGAL-edge graph (matches defect_catalog's
            `sig`, which is what a human reading the catalog sees)
  induced   components of the subgraph induced on defect vertices using ALL
            edges (what contact_census.py used) -- fuses complexes that are
            joined by a legal edge, so it reports fewer, larger complexes

Two averages, because they answer different questions:
  per complex   what fraction of complexes hold a deg-3 edge
  per vertex    size-weighted -- a flight is likelier to collide with a big
                complex, so this is closer to "what will a blocker look like"

Usage:
    python scripts/defect_dynamics/complex_d3_census.py \
        data/mgas/lam40_snap*.mfd data/mgas/lam35*_snap*.mfd \
        --out data/ecmc_ab/complex_d3_census.json
"""
import argparse
import collections
import json
import os
import sys

import numpy as np

_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(_ROOT, "python"))

import discrete_differential_geometry as ddg
from discrete_differential_geometry.vertex_fields import edges_and_degrees

BINS = [(0, 6, "<=6  (flicker-scale)"), (7, 9, "7-9"),
        (10, 14, "10-14"), (15, 10 ** 9, ">=15")]


def components(pairs, verts):
    """Connected components over `pairs`, restricted to `verts`."""
    adj = collections.defaultdict(set)
    for u, v in pairs:
        adj[u].add(v)
        adj[v].add(u)
    seen, out = set(), []
    for s0 in verts:
        if s0 in seen:
            continue
        stack, cur = [s0], set()
        while stack:
            x = stack.pop()
            if x in cur:
                continue
            cur.add(x)
            stack.extend(adj[x] - cur)
        seen |= cur
        out.append(cur)
    return out


def census_one(path):
    F = np.asarray(ddg.Manifold.load(path, 3).facets())
    eu, ecnt, deg, V = edges_and_degrees(F)
    ill = (ecnt < 5) | (ecnt > 6)
    ill_e = eu[ill]
    ill_d = ecnt[ill]
    defect = set(int(v) for e in ill_e for v in e)

    out = {}
    for conv, pairs in (("illegal", ill_e),
                        ("induced", eu[np.isin(eu[:, 0], list(defect))
                                       & np.isin(eu[:, 1], list(defect))])):
        comps = components([tuple(map(int, p)) for p in pairs], defect)
        rows = []
        for c in comps:
            # illegal edges with BOTH ends in this complex, and their degrees
            m = np.array([e[0] in c and e[1] in c for e in ill_e], bool) \
                if len(ill_e) else np.zeros(0, bool)
            degs = sorted(int(d) for d in ill_d[m])
            rows.append({"n": len(c), "sig": degs,
                         "n_d3": sum(1 for d in degs if d == 3),
                         "n_d4": sum(1 for d in degs if d == 4)})
        out[conv] = rows
    return out


def report(rows, label):
    if not rows:
        print(f"\n=== {label}: none")
        return {}
    n_tot = len(rows)
    sz = np.array([r["n"] for r in rows], float)
    has = np.array([r["n_d3"] > 0 for r in rows], float)
    print(f"\n=== {label}   ({n_tot:,} complexes, {int(sz.sum()):,} defect vertices)")
    print(f"    per complex: {has.mean():6.1%} hold a deg-3 edge")
    print(f"    per vertex (size-weighted): "
          f"{float((has * sz).sum() / sz.sum()):6.1%}")
    print(f"    {'size bin':<22} {'n_cplx':>7} {'has d3':>8} {'mean n_d3':>10} "
          f"{'mean n_d4':>10}")
    out = {"n_complexes": n_tot, "per_complex": float(has.mean()),
           "per_vertex": float((has * sz).sum() / sz.sum()), "bins": {}}
    for lo, hi, lab in BINS:
        g = [r for r in rows if lo <= r["n"] <= hi]
        if not g:
            continue
        h = np.mean([r["n_d3"] > 0 for r in g])
        out["bins"][lab] = {"n": len(g), "has_d3": float(h),
                            "mean_n_d3": float(np.mean([r["n_d3"] for r in g])),
                            "mean_n_d4": float(np.mean([r["n_d4"] for r in g]))}
        print(f"    {lab:<22} {len(g):>7,} {h:>7.1%} "
              f"{np.mean([r['n_d3'] for r in g]):>10.2f} "
              f"{np.mean([r['n_d4'] for r in g]):>10.2f}")
    # the elementary flicker: n=5, sig (3,4,4)
    fl = [r for r in rows if r["sig"] == [3, 4, 4]]
    out["flicker_344"] = {"n": len(fl), "frac": len(fl) / n_tot}
    print(f"    elementary flicker sig=(3,4,4): {len(fl):,} "
          f"= {len(fl)/n_tot:.1%} of complexes, "
          f"{sum(r['n_d3'] for r in fl)}/{sum(r['n_d3'] for r in rows)} "
          f"of all deg-3 edges")
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("snapshots", nargs="+")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    pool = {"illegal": [], "induced": []}
    for p in args.snapshots:
        try:
            r = census_one(p)
        except Exception as ex:                       # noqa: BLE001
            print(f"  SKIP {os.path.basename(p)}: {ex}")
            continue
        for k in pool:
            pool[k] += r[k]
        print(f"  {os.path.basename(p):<28} "
              f"{len(r['illegal']):>3} complexes (illegal-graph)", flush=True)

    summ = {}
    for conv in ("illegal", "induced"):
        summ[conv] = report(pool[conv],
                            f"{conv.upper()}-graph components  "
                            f"[{'catalog sig convention' if conv == 'illegal' else 'contact_census convention'}]")

    if args.out:
        os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
        with open(args.out, "w") as f:
            json.dump({"summary": summ, "snapshots": args.snapshots,
                       "complexes": pool}, f, indent=1, default=float)
        print(f"\nwrote {os.path.abspath(args.out)}")


if __name__ == "__main__":
    main()
