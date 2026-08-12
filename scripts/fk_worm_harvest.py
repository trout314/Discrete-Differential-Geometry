#!/usr/bin/env python3
"""Harvest FK-legal triangulations from a chain that carries an illegality budget.

THE IDEA

No single move -- Pachner, hinge, contraction or split -- maps an FK-legal
state to another FK-legal state (proof and library census in
scripts/defect_dynamics/fk_channel_census.py), so the legal set is dust and no
strictly-legal sampler exists. The way around it is the standard trick for a
constrained system: let the chain carry a small, BOUNDED reservoir of defects
and read off the configurations at the moments the reservoir empties.

Cap the illegal-edge count at B (``ManifoldSampler.set_illegal_budget``). The
cap is an infinite energy, so the chain is an ordinary reversible Metropolis
chain on the capped set. On the legal manifold at fixed (f0, f3) the action is
EXACTLY degenerate -- volume pin, flat pin, U(n6) and V(m) are all constant
there -- so the sub-collection of visits with n_ill = 0 is an exact sample of
the UNIFORM measure over FK-legal triangulations with those counts. Crystals
are measure-zero in that measure; amorphous states are essentially all of it.

The defects are transport, not physics. The only question this script exists to
answer is whether the transport WORKS:

    do successive returns to n_ill = 0 land on DIFFERENT legal states,
    and how far apart are they?

If yes, the program is a matter of running longer. If the chain keeps returning
to the same configuration, the reservoir is too small or too expensive and the
budget needs raising (note the flicker quantum is 3-5 illegal edges, so B below
~10 leaves no room for even one excursion to travel).

WHAT IS REPORTED

    returns        chunks observed at n_ill = 0
    distinct       distinct legal states among them (exact facet-set hash)
    d_facets       facets changed since the PREVIOUS legal state, as a
                   fraction of f3 -- 0 means the chain came straight back
    d_six          symmetric difference of the degree-6 edge sets (the
                   six-web) between successive legal states, normalised.
                   This is the physically meaningful motion: the web IS the
                   state's curvature content.
    blocked        proposals rejected by the cap alone (is the cap binding?)

Usage:
    python scripts/fk_worm_harvest.py \
        --start data/tcp_reference/T3_R_m2_N7248.mfd \
        --budget 24 --sweeps 4000 --cs 0.05 --out data/fk_worm/r_B24
"""
import argparse
import hashlib
import json
import math
import os
import sys
import time

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_ROOT, "python"))
import numpy as np
import discrete_differential_geometry as ddg
from discrete_differential_geometry.recording import Recorder

E_FLAT = 2.0 * math.pi / math.acos(1.0 / 3.0)
FK_N6 = (0, 2, 3, 4)


def facet_hash(v):
    """Exact fingerprint of the labelled triangulation."""
    f = np.asarray(v.facets())
    f = np.sort(f, axis=1)
    f = f[np.lexsort(f.T[::-1])]
    return hashlib.blake2b(f.astype(np.int32).tobytes(), digest_size=16).hexdigest()


def six_web(v):
    """Set of degree-6 edges -- the disclination network."""
    edges = np.asarray(v.simplices(1))
    return {(int(a), int(b)) for a, b in edges
            if v.degree((int(a), int(b))) >= 6}


def fk_counts(v, cap=10):
    """(n_fk, n_hub, n_imp, f0) from the joint (n6, m) census."""
    c = np.asarray(v.valence_census(cap, cap))
    f0 = int(c.sum())
    n_imp = int(c[:, 1:].sum())
    n_fk = int(sum(c[k, 0] for k in FK_N6))
    return n_fk, f0 - n_imp - n_fk, n_imp, f0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--start", required=True)
    ap.add_argument("--budget", type=int, required=True,
                    help="illegal-edge cap B")
    ap.add_argument("--f3-target", type=int, default=0)
    ap.add_argument("--e-target", type=float, default=E_FLAT)
    ap.add_argument("--cn", type=float, default=0.1)
    ap.add_argument("--nh", type=float, default=30.0)
    ap.add_argument("--zleg", type=float, default=1.0,
                    help="hub penalty. Must be nonzero somewhere in the "
                         "potential for the budget's counter to be maintained")
    ap.add_argument("--mu-ill", type=float, default=0.0,
                    help="illegal-edge fugacity. 0 = a PURE budget (flat "
                         "inside the reservoir), which is the intended "
                         "setting: the cap does the confining, and a flat "
                         "interior maximises how far excursions travel")
    ap.add_argument("--slide", type=float, default=0.0,
                    help="CLEAN knot-slide probability. This is the worm's "
                         "transport move: a clean slide preserves the "
                         "illegal-degree multiset exactly, so it cannot "
                         "breach the cap, and without it the reservoir only "
                         "does birth/death and every excursion retraces.")
    ap.add_argument("--cs", type=float, default=0.05)
    ap.add_argument("--max-ring", type=int, default=6)
    ap.add_argument("--sweeps", type=int, default=4000)
    ap.add_argument("--chunk", type=float, default=5.0,
                    help="sweeps between legality checks; small enough to "
                         "catch returns, large enough to be cheap")
    ap.add_argument("--save-legal", type=int, default=0,
                    help="save the first K distinct legal states as .mfd")
    ap.add_argument("--seed", type=int, default=20260810)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    ddg.set_random_seed(args.seed)
    mfd = ddg.Manifold.load(args.start, 3)
    f3_target = args.f3_target or mfd.num_facets

    params = ddg.SamplerParams(
        num_facets_target=f3_target, num_facets_coef=args.cn,
        hinge_degree_target=args.e_target, num_hinges_coef=args.nh,
        hinge_degree_variance_coef=0.0, codim3_degree_variance_coef=0.0,
        hinge_degree_target_coef=0.0, codim3_degree_target_coef=0.0)
    s = ddg.ManifoldSampler(mfd, params)
    s.set_n6_potential(args.zleg, 0.0, tilt=[0.0] * 5, imp_lin=args.mu_ill)
    if args.slide > 0:
        s.set_slide_clean_only(True)      # required alongside the budget
    s.set_illegal_budget(args.budget)
    if args.slide > 0:
        s.set_slide_prob(args.slide)
    if args.cs > 0:
        s.set_contract_split(args.cs, args.max_ring)
    v = s.manifold

    os.makedirs(os.path.dirname(os.path.abspath(args.out)) or ".", exist_ok=True)
    meta = {"start": os.path.basename(args.start), "budget": args.budget,
            "f3_target": f3_target, "e_target": args.e_target, "cn": args.cn,
            "nh": args.nh, "zleg": args.zleg, "mu_ill": args.mu_ill,
            "cs": args.cs, "max_ring": args.max_ring, "sweeps": args.sweeps,
            "chunk": args.chunk, "seed": args.seed, "e_flat": E_FLAT}
    rec = Recorder(s, args.out, chunk=args.chunk, census_every=1,
                   snap_mid=False, extra_meta=meta)
    obs = open(args.out + ".obs.jsonl", "w")
    obs.write(json.dumps({"kind": "header", **meta}) + "\n")

    seen = {}
    prev_hash = None
    prev_web = None
    prev_facets = None
    returns = 0
    saved = 0
    t0 = time.time()

    while rec.sw < args.sweeps:
        rec.step()
        n_ill = len(v.illegal_edges()[0])
        cap, maintained, blocked = s.illegal_budget_stats()
        assert maintained == n_ill, (
            f"budget counter {maintained} != scan {n_ill}")
        # Defensive: the clean-slide channel is illegal-count-preserving by
        # construction, so nothing may tunnel past the cap. Checked rather
        # than trusted -- a silent breach would poison the harvest.
        assert n_ill <= args.budget, (
            f"budget breached: {n_ill} > {args.budget}")

        row = {"sw": rec.sw, "n_ill": n_ill, "blocked": blocked,
               "slide": s.slide_stats() if args.slide > 0 else None}
        if n_ill == 0:
            returns += 1
            h = facet_hash(v)
            web = six_web(v)
            facets = {tuple(sorted(int(x) for x in f)) for f in v.facets()}
            n_fk, n_hub, n_imp, f0 = fk_counts(v)
            row.update(novel=h not in seen, n_hub=n_hub, f0=f0,
                       f_fk=100.0 * n_fk / f0,
                       e_mean=6.0 * int(v.f_vector[3]) / int(v.f_vector[1]))
            if prev_hash is not None:
                row["d_facets"] = len(facets ^ prev_facets) / (2 * len(facets))
                row["d_six"] = (len(web ^ prev_web) / (2 * len(web))
                                if web else 0.0)
                row["same_as_prev"] = (h == prev_hash)
            if h not in seen:
                seen[h] = rec.sw
                if saved < args.save_legal:
                    v.save(f"{args.out}.legal{saved:03d}.mfd")
                    saved += 1
            prev_hash, prev_web, prev_facets = h, web, facets

        obs.write(json.dumps(row) + "\n")
        obs.flush()
        print(f"sw {rec.sw:7.0f} n_ill {n_ill:4d} | returns {returns:4d} "
              f"distinct {len(seen):4d} | blocked {blocked:11d} "
              f"[{time.time() - t0:.0f}s]", flush=True)

    v.save(args.out + ".final.mfd")
    obs.close()

    print(f"\n=== budget B = {args.budget}, {args.sweeps} sweeps")
    print(f"returns to n_ill = 0 : {returns}")
    print(f"distinct legal states: {len(seen)}")
    if returns and len(seen) <= 1:
        print("VERDICT: the chain returns to the SAME legal state -- the "
              "reservoir is too small or too costly to transport anything. "
              "Raise --budget (the flicker quantum is 3-5 illegal edges) or "
              "lower --mu-ill.")
    elif len(seen) > 1:
        print("VERDICT: the budgeted chain MOVES on the legal manifold.")
    print(f"final state: {os.path.abspath(args.out)}.final.mfd")


if __name__ == "__main__":
    main()
