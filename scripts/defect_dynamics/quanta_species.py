#!/usr/bin/env python3
"""Ensemble species catalog for few-quantum-strained states (quanta_heal.py).

A single snapshot of a DILUTE state holds only a handful of complexes, so a
species list read off one `.mfd` is anecdote. This drives each cell for a
while and censuses EVERY snapshot, so a species gets a population (mean count
per snapshot) with an error bar rather than a tick mark.

Species naming follows defect_state.Complex: the primary name is the
historical illegal-edge-degree signature `sig` (so "(3,4,4)" still means
exactly what it always has), with non-FK coordination anomalies `nodes`
reported alongside rather than folded in, plus the exact decorated-
isomorphism hash (`canonical_key` with ecolor = ambient edge degree,
vcolor = n6) so accidental sig collisions are visible.

Complexes are components of the ILLEGAL-EDGE GRAPH (`--narrow` in
defect_catalog): on a15 any vertex-adjacency definition percolates through
the native disclination network and fuses every complex into one.

Per species it reports population, size in vertices, the integer curvature
charge Q = sum_{illegal e} (5 - deg) (deg-4 -> +1, deg-3 -> +2, deg-7 -> -2),
the total coordination sum Z that quantises the raw charge Q_c, and the
illegal-subgraph anatomy (edges / components / cycles) -- the curve-like
invariant, not the dense induced one.

Run several cells in one invocation so the table is a comparison: the point
of a below-native catalog is what it is NOT (the deg-4 multimer ladder that
the above-native and native gases show).

Usage:
    python scripts/defect_dynamics/quanta_species.py \
        --cell "below dq+3=data/quanta_heal/uncap/uncap_isostrain.A.final.mfd,3" \
        --cell "above dq-3=data/quanta_heal/main/above_dq-3_c0.9_s11.A.final.mfd,-3" \
        --cimp 0.6 --sweeps 4000 --out data/quanta_heal/species/s
"""
import argparse
import hashlib
import json
import os
import sys
from collections import Counter, defaultdict

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(os.path.dirname(_HERE))
sys.path.insert(0, os.path.join(_ROOT, "python"))
sys.path.insert(0, os.path.join(_ROOT, "scripts"))
sys.path.insert(0, _HERE)

import discrete_differential_geometry as ddg
from discrete_differential_geometry.recording import Recorder
import defect_state as dsm
from crystal_gas_scan import LIBRARY
from quanta_heal import edge_degree_target


def narrow_components(st):
    """Components of the illegal-edge graph, as defect_state.Complex."""
    parent, e2i = {}, {}
    ill = sorted(st.ill_edges)

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for i, e in enumerate(ill):
        parent[i] = i
        for v in e:
            if v in e2i:
                r1, r2 = find(i), find(e2i[v])
                if r1 != r2:
                    parent[r1] = r2
            e2i[v] = i
    groups = defaultdict(list)
    for i, e in enumerate(ill):
        groups[find(i)].append(e)
    out = []
    for edges in groups.values():
        cv = sorted({v for e in edges for v in e})
        out.append(dsm.Complex(cv, tuple(sorted(st.edeg[e] for e in edges)),
                               tuple(sorted(st.n6[v] for v in cv
                                            if st.imp[v] == 0
                                            and st.n6[v] not in dsm.FK_N6))))
    return sorted(out, key=lambda c: -len(c.verts))


def canon_hash(st, c, max_size=60):
    """Short hash of the decorated-isomorphism class of a defect complex.

    Same object defect_catalog hashes: the INDUCED subcomplex on the defect's
    own vertices, decorated by ambient edge degree and n6 -- not the closed
    star, whose ~60 vertices put the canonicalisation search past its limit
    (measured: 14 s per complex, and it returns exact=False anyway)."""
    if len(c.verts) > max_size:
        return "big", False
    fac_i = st.induced_facets(c.verts)
    vc, ec = st.decorations(c.verts)
    key, exact = dsm.canonical_key(fac_i, vc, ec)
    return hashlib.blake2b(repr(key).encode(), digest_size=4).hexdigest(), exact


def run_cell(name, start, dq, args, cell_path, mcell):
    ddg.set_random_seed(args.seed)
    ref = ddg.Manifold.load(cell_path, 3)
    f0n, _, _, f3n = (int(x) for x in ref.f_vector)
    f3t = f3n - dq
    mfd = ddg.Manifold.load(start, 3)
    params = ddg.SamplerParams(
        num_facets_target=f3t, num_facets_coef=args.vol_coef,
        hinge_degree_target=edge_degree_target(f0n, f3t),
        num_hinges_coef=args.hinge_coef,
        hinge_degree_variance_coef=0.0, codim3_degree_variance_coef=0.0,
        hinge_degree_target_coef=0.0, codim3_degree_target_coef=0.0)
    s = ddg.ManifoldSampler(mfd, params)
    s.set_n6_potential(0.0, args.cimp, tilt=[0.0] * 5)
    v = s.manifold

    out = f"{args.out}.{name.replace(' ', '_').replace('+', 'p')}"
    rec = Recorder(s, out, chunk=args.chunk, census_every=1, snap_mid=False,
                   extra_meta={"cell": name, "start": start, "dq": dq,
                               "cimp": args.cimp, "f3_target": f3t,
                               "vol_coef": args.vol_coef,
                               "hinge_coef": args.hinge_coef})
    # Population statistics use the cheap per-snapshot quantities on EVERY
    # snapshot; the exact decorated-isomorphism hash costs a canonicalisation
    # search per complex, so it is subsampled (--canon-every) and reported
    # only as a class count within each species.
    pop = defaultdict(list)          # (sig, nodes) -> per-snapshot count
    info, hashes, nsnap, best = {}, defaultdict(set), 0, (-1, None)
    ncanon = 0
    while rec.sw < args.sweeps:
        rec.step(min(args.chunk, args.sweeps - rec.sw))
        if rec.sw < args.burn:
            continue
        st = dsm.DefectState(v)
        comps = narrow_components(st)
        nsnap += 1
        do_canon = (nsnap % args.canon_every == 0)
        ncanon += do_canon
        here = Counter()
        for c in comps:
            k = (c.sig, c.nodes)
            here[k] += 1
            if k not in info:
                sh = st.illegal_shape(c.verts)
                info[k] = {
                    "size": len(c.verts), "Q": sum(5 - d for d in c.sig),
                    "sumZ": st.total_coordination(c.verts),
                    "Qc": round(st.complex_charge(c.verts), 4),
                    "ill_edges": sh["n_edges"], "ill_comps": sh["components"],
                    "cycles": sh["cycles"]}
            if do_canon:
                h, exact = canon_hash(st, c)
                hashes[k].add(h)
                info[k]["exact"] = bool(exact)
        for k in set(pop) | set(here):
            pop[k].append(here.get(k, 0))
            if len(pop[k]) < nsnap:      # first sighting: back-fill zeros
                pop[k] = [0] * (nsnap - 1) + [here.get(k, 0)]
        if len(here) > best[0]:
            best = (len(here), rec.sw)
            v.save(out + ".rep.mfd",
                   [f"representative snapshot, {name}, sweep {rec.sw}",
                    f"dq={dq} c_imp={args.cimp} species here={len(here)}"])
    rec.finish()
    rows = []
    for k, series in pop.items():
        a = np.asarray(series, float)
        rows.append({"sig": list(k[0]), "nodes": list(k[1]),
                     "classes": len(hashes[k]),
                     "hashes": sorted(hashes[k]),
                     "pop": float(a.mean()), "pop_sd": float(a.std()),
                     "seen": int((a > 0).sum()), "nsnap": nsnap, **info[k]})
    rows.sort(key=lambda r: -r["pop"])
    return {"cell": name, "dq": dq, "start": start, "nsnap": nsnap,
            "n_canon_snapshots": ncanon, "rep": out + ".rep.mfd",
            "rows": rows}


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--cell", action="append", required=True,
                    help='"NAME=PATH.mfd,DQ" -- repeatable')
    ap.add_argument("--crystal", default="a15", choices=sorted(LIBRARY))
    ap.add_argument("--cimp", type=float, default=0.6)
    ap.add_argument("--vol-coef", type=float, default=1.0)
    ap.add_argument("--hinge-coef", type=float, default=30.0)
    ap.add_argument("--sweeps", type=int, default=4000)
    ap.add_argument("--burn", type=int, default=500)
    ap.add_argument("--chunk", type=int, default=50)
    ap.add_argument("--canon-every", type=int, default=10,
                    help="canonicalise (exact decorated-isomorphism hash) "
                         "every Nth snapshot only -- the search dominates "
                         "the census cost, and the class count it feeds is "
                         "a secondary column")
    ap.add_argument("--seed", type=int, default=31)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    fname, mcell = LIBRARY[args.crystal]
    cell_path = os.path.join(_ROOT, "data", "tcp_reference", fname)
    os.makedirs(os.path.dirname(os.path.abspath(args.out)) or ".",
                exist_ok=True)

    results = []
    for spec in args.cell:
        name, rest = spec.split("=", 1)
        path, dq = rest.rsplit(",", 1)
        print(f"[{name}] dq={dq} from {os.path.basename(path)} ...",
              flush=True)
        results.append(run_cell(name, path, int(dq), args, cell_path, mcell))

    with open(args.out + ".species.json", "w") as fh:
        json.dump({"cimp": args.cimp, "sweeps": args.sweeps,
                   "burn": args.burn, "seed": args.seed,
                   "vol_coef": args.vol_coef, "hinge_coef": args.hinge_coef,
                   "cells": results}, fh, indent=1)

    for r in results:
        tot = sum(x["pop"] for x in r["rows"])
        # The primary species name is (sig, nodes); the canonical hashes
        # inside it are distinct EMBEDDINGS of the same illegal signature in
        # the crystal, so they are reported as a class count, not as rows.
        ncls = sum(x["classes"] for x in r["rows"])
        print(f"\n=== {r['cell']}   (dq = {r['dq']:+d}, c_imp = {args.cimp}, "
              f"{r['nsnap']} snapshots, {len(r['rows'])} species, "
              f">= {ncls} decorated classes from "
              f"{r['n_canon_snapshots']} canonicalised snapshots, "
              f"{tot:.2f} complexes per snapshot)")
        print(f"{'sig (illegal degrees)':<24}{'nodes':>8}{'pop':>8}{'sd':>7}"
              f"{'%':>6}{'cls':>5}{'size':>6}{'Q':>4}{'sumZ':>6}"
              f"{'ill e/c/cyc':>13}")
        for x in r["rows"]:
            if x["pop"] < 0.005:
                continue
            anat = f"{x['ill_edges']}/{x['ill_comps']}/{x['cycles']}"
            print(f"{str(tuple(x['sig'])):<24}{str(tuple(x['nodes'])):>8}"
                  f"{x['pop']:>8.2f}{x['pop_sd']:>7.2f}"
                  f"{100*x['pop']/tot:>6.1f}{x['classes']:>5}"
                  f"{x['size']:>6}{x['Q']:>4}{x['sumZ']:>6}{anat:>13}")
        deg = Counter()
        for x in r["rows"]:
            for d in x["sig"]:
                deg[d] += x["pop"]
        s = sum(deg.values()) or 1
        print("  illegal-edge budget: " + ", ".join(
            f"deg-{d}: {deg[d]:.2f} ({100*deg[d]/s:.0f}%)"
            for d in sorted(deg)) +
            f"   |   net Q = {sum(x['Q']*x['pop'] for x in r['rows']):+.2f}"
            f" per snapshot")
        print(f"  representative snapshot: {r['rep']}")
    print(f"\nJSON: {os.path.abspath(args.out)}.species.json")


if __name__ == "__main__":
    main()
