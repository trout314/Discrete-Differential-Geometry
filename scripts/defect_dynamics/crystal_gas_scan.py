#!/usr/bin/env python3
"""Is there a stable dilute defect gas on THIS crystal, under the minimal action?

One run = one (crystal, c_imp) cell. The action is the user's minimal three-term
one -- nothing else is switched on:

    S = 0.1 * (f3 - f3_ref)^2                     volume pin
      + 1.0 * (f1 - 6 f3 / e*)^2                  mean-edge-degree pin
      + c_imp * sum_v m(v)^2                      anti-clustering

with e* = 2*pi/arccos(1/3) = 5.1042993 the FLAT edge-degree target (zero mean
Regge deficit) and m(v) = # edges at v whose degree is not in {5, 6} (the
"impurity valence"). No per-edge EDQ term, no U(n6) legality term, no degree
variance penalties -- the variance couplings default to ON in SamplerParams and
fight an FK/TCP pin, so they are explicitly zeroed here (see
notes/memory/vdv-hdv-conflicts-with-tcp.md).

WHY THE GAS EXISTS AT ALL

The two pins are not simultaneously satisfiable by the pristine crystal: at
fixed f3, the flat target demands f1 = 6 f3 / e*, but the crystal's native mean
edge degree e_nat = 6 - 12/CN is generally != e*. The mismatch

    gap = f1_native - 6 f3_native / e*

is an EXTENSIVE debt the sampler must pay in defects. Each 2->3 move changes it
by (1 - 6/e*) = -0.1755, each 1->4 by (4 - 18/e*) = +0.4734, so the crystal is
forced to carry ~|gap| / 0.1755 defect-creating moves. The SIGN of the gap sets
which move type is forced -- crystals with e_nat < e* (c15, c14, c36, mu, r)
pay in 2->3 flickers, crystals with e_nat > e* (a15, sigma, z, p, delta) pay in
the reverse. That is a per-crystal prediction this scan tests directly.

c_imp is then the only free dial, and it decides ARRANGEMENT, not number: too
small and the forced defects condense into one percolating clump (the melt);
large enough and they disperse into a gas. This scan finds that crossover per
crystal.

OBSERVABLES (per chunk; all center-free)

    n_ill     illegal-edge-incident vertices -- the defect population
    ncomp     connected defect complexes; = n_ill/typ.size when dispersed, 1 when melted
    top1      largest complex, in vertices -- the concrete failure mode
    frac_top1 top1 / n_ill: 1.0 = single clump, ~0 = gas
    max_m     largest impurity valence anywhere (deep burial)
    sum_m2    the clustering energy itself
    turnover  1 - Jaccard(defect vertex set now, one chunk ago) -- gas vs glass
    e_mean    achieved mean edge degree (pin satisfaction)

Usage:
    python scripts/defect_dynamics/crystal_gas_scan.py --crystal c15 --mcell 4 \
        --cimp 0.5 --burn 1500 --span 3000 --out data/crystal_gas/c15_c0.50
"""
import argparse
import json
import math
import os
import sys
import time
from collections import Counter

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(os.path.dirname(_HERE))
sys.path.insert(0, os.path.join(_ROOT, "python"))
sys.path.insert(0, os.path.join(_ROOT, "scripts"))
sys.path.insert(0, _HERE)

import discrete_differential_geometry as ddg
from discrete_differential_geometry import cocycle as coc
from discrete_differential_geometry.recording import Recorder
from cocycle_check import reference_frac_positions
from tcp_reference import STRUCTURES

E_FLAT = 2.0 * math.pi / math.acos(1.0 / 3.0)      # 5.1042993...

# The reference cell used per crystal. Chosen so every host lands in
# f0 = 1000..2500 vertices, which keeps the forced defect debt comparable
# across the library and every run inside a couple of minutes.
LIBRARY = {
    "a15":   ("T3_A15_m5_N5750.mfd",     5),
    "c15":   ("T3_C15_m4_N8704.mfd",     4),
    "c14":   ("T3_C14_m5_N8500.mfd",     5),
    "c36":   ("T3_C36_m4_N8704.mfd",     4),
    "sigma": ("T3_SIGMA_m4_N11008.mfd",  4),
    "z":     ("T3_Z_m6_N8640.mfd",       6),
    "mu":    ("T3_MU_m4_N14208.mfd",     4),
    "r":     ("T3_R_m2_N7248.mfd",       2),
    "p":     ("T3_P_m3_N8640.mfd",       3),
    "delta": ("T3_DELTA_m3_N8640.mfd",   3),
}


def pin_gap(name, mcell):
    """(e_native, gap in f1 units, forced 2->3-equivalent move count)."""
    _, sites, cn, _ = STRUCTURES[name]
    f0 = len(sites) * mcell ** 3
    f3, f1 = f0 * (cn / 2 - 1), f0 * cn / 2
    gap = f1 - 6.0 * f3 / E_FLAT
    return 6.0 - 12.0 / cn, gap, abs(gap) / abs(1.0 - 6.0 / E_FLAT)


def vertex_components(pairs):
    """Vertex sets of the connected components of the illegal-edge graph.

    recording._components returns component sizes in EDGES; every defect
    census in the tree (defect_state.census, mobile_gas) counts VERTICES, so
    sizes must be built over vertices to stay comparable with those series."""
    parent = {}

    def find(x):
        parent.setdefault(x, x)
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for a, b in map(tuple, pairs):
        ra, rb = find(int(a)), find(int(b))
        if ra != rb:
            parent[ra] = rb
    groups = {}
    for x in parent:
        groups.setdefault(find(x), []).append(x)
    return sorted(groups.values(), key=len, reverse=True)


def defect_stats(view):
    """Dispersal anatomy from the illegal-edge set alone (cheap, exact).

    m(v) is by definition the count of incident illegal edges, so the whole
    impurity ledger falls out of `illegal_edges()` with no O(N) rebuild."""
    pairs, degs = view.illegal_edges()
    pairs = np.asarray(pairs).reshape(-1, 2)
    degs = np.asarray(degs)
    if len(pairs) == 0:
        return {"n_ill": 0, "ncomp": 0, "top1": 0, "frac_top1": 0.0,
                "max_m": 0, "sum_m2": 0, "sizes": [], "deg_hist": {},
                "verts": frozenset()}
    m = Counter(pairs.ravel().tolist())
    comps = vertex_components(pairs)
    sizes = [len(c) for c in comps]
    n_ill = sum(sizes)
    return {
        "n_ill": n_ill,
        "ncomp": len(comps),
        "top1": sizes[0],
        "frac_top1": sizes[0] / n_ill,
        "max_m": max(m.values()),
        "sum_m2": int(sum(v * v for v in m.values())),
        "sizes": sizes[:12],
        "deg_hist": {int(k): int(c) for k, c in sorted(Counter(degs.tolist()).items())},
        "verts": frozenset(m),
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--crystal", required=True, choices=sorted(LIBRARY))
    ap.add_argument("--cimp", type=float, required=True)
    ap.add_argument("--burn", type=int, default=1500)
    ap.add_argument("--span", type=int, default=3000)
    ap.add_argument("--chunk", type=int, default=50)
    ap.add_argument("--seed", type=int, default=20260806)
    ap.add_argument("--contract-split", type=float, default=0.0,
                    help="probability of the contract/split channel per "
                         "MCMC step (0 = off). Unpins f0: the flat pin's "
                         "gap can then be paid in vertices instead of "
                         "defects (f1 = f0 + f3).")
    ap.add_argument("--vol-scale", type=float, default=1.0,
                    help="volume-pin target as a multiple of the native "
                         "f3_ref (1.0 = native). With the channel on, the "
                         "flat-pin line f0*(f3) = f3(6/e* - 1) is reachable "
                         "at any volume -- the question is what STRUCTURE "
                         "a non-native volume forces in the fixed box.")
    ap.add_argument("--out", required=True)
    ap.add_argument("--cocycle", action="store_true",
                    help="attach the harmonic cocycle at sweep 0 (needed for "
                         "the 3D viewer; post-hoc regen FAILS on these gases)")
    args = ap.parse_args()

    fname, mcell = LIBRARY[args.crystal]
    cell = os.path.join(_ROOT, "data", "tcp_reference", fname)
    e_nat, gap, n_forced = pin_gap(args.crystal, mcell)

    ddg.set_random_seed(args.seed)
    ref = ddg.Manifold.load(cell, 3)
    f3_ref = ref.num_facets
    f3_target = int(round(f3_ref * args.vol_scale))
    mfd = ddg.Manifold.load(cell, 3)

    params = ddg.SamplerParams(
        num_facets_target=f3_target, num_facets_coef=0.1,
        hinge_degree_target=E_FLAT, num_hinges_coef=1.0,
        # everything else OFF: the defaults are ON and fight an FK pin
        hinge_degree_variance_coef=0.0, codim3_degree_variance_coef=0.0,
        hinge_degree_target_coef=0.0, codim3_degree_target_coef=0.0)
    s = ddg.ManifoldSampler(mfd, params)
    s.set_n6_potential(0.0, args.cimp, tilt=[0.0] * 5)   # zleg = 0: m^2 only
    if args.contract_split > 0:
        s.set_contract_split(args.contract_split, 6)
    v = s.manifold

    if args.cocycle:
        edges = np.asarray(v.simplices(1))
        s.enable_cocycle(edges, coc.build_from_positions(
            edges, reference_frac_positions(args.crystal, mcell), mcell))

    os.makedirs(os.path.dirname(os.path.abspath(args.out)) or ".",
                exist_ok=True)
    meta = {"crystal": args.crystal, "cell": os.path.basename(cell),
            "mcell": mcell, "cimp": args.cimp, "e_flat": E_FLAT,
            "e_native": e_nat, "pin_gap_f1": gap, "n_forced_moves": n_forced,
            "burn": args.burn, "span": args.span, "f3_ref": f3_ref,
            "contract_split": args.contract_split,
            "vol_scale": args.vol_scale, "f3_target": f3_target}
    rec = Recorder(s, args.out, chunk=args.chunk, census_every=1,
                   snap_mid=False,
                   cocycle_box=mcell if args.cocycle else None,
                   extra_meta=meta)
    obs = open(args.out + ".obs.jsonl", "w")
    obs.write(json.dumps({"kind": "header", **meta}) + "\n")

    t0 = time.time()
    prev = None
    total = args.burn + args.span
    while rec.sw < total:
        row = rec.step(min(args.chunk, total - rec.sw))
        d = defect_stats(v)
        cur = d.pop("verts")
        d["turnover"] = (None if prev is None else
                         (1.0 - len(cur & prev) / len(cur | prev)
                          if (cur | prev) else 0.0))
        prev = cur
        fv = [int(x) for x in v.f_vector]
        d.update(sw=rec.sw, t=round(time.time() - t0, 1),
                 burned=rec.sw > args.burn, f=fv,
                 e_mean=6.0 * fv[3] / fv[1], obj=float(s.current_objective))
        if args.contract_split > 0:
            d["cs"] = list(s.contract_split_stats())
        obs.write(json.dumps(d) + "\n")
        obs.flush()

    path = rec.finish()
    obs.write(json.dumps({"kind": "final", "path": path,
                          "t": round(time.time() - t0, 1)}) + "\n")
    obs.close()
    print(f"[{args.crystal} c={args.cimp}] DONE {rec.sw} sw "
          f"({time.time()-t0:.0f}s) n_ill={d['n_ill']} ncomp={d['ncomp']} "
          f"top1={d['top1']} max_m={d['max_m']} "
          f"<edeg>={d['e_mean']:.6f} -> {path}", flush=True)


if __name__ == "__main__":
    main()
