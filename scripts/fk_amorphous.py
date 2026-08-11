#!/usr/bin/env python3
"""Drive a state toward an AMORPHOUS Frank-Kasper triangulation: n_ill -> 0,
all vertices Z12/14/15/16, mean edge degree pinned near the flat value.

WHY THIS IS A DIFFERENT EXPERIMENT FROM THE OLD VDV ANNEAL

Two exact facts set the strategy (both verified in
scripts/defect_dynamics/fk_channel_census.py):

1. With every edge degree in {5,6}, the link sum rule gives Z(v) = 12 + n6(v)
   and hence

       ebar = 6 - 12/Zbar,    Zbar = 12 + n6bar,

   so the mean edge degree is NOT an independent knob -- it is exactly the
   six-web line density.  ebar = e* = 2*pi/arccos(1/3) means Zbar = 13.39733,
   which sits inside the FK window [13.333 (C15), 13.5 (A15)].  Scanning ebar
   near flat = scanning the Z12:Z14:Z15:Z16 composition.

2. No single move -- Pachner, hinge, slide, contraction or split -- maps an
   FK-legal state to another FK-legal state.  For contraction the proof is
   two lines: edge-legality of the result forces every merged spoke at the
   surviving vertex w to have degree 6, so n6(w) >= deg(uv) >= 5 and w is a
   Z17+ hub.  Measured on the reference library: the R crystal has ZERO
   contractions with even non-positive damage, minimum (d_ill, d_nonFK) =
   (1, 3).  Perfect FK crystals are isolated points.

So the FK-legal set is dust, and any sampler must travel through a controlled
amount of illegality.  What is new is HOW MUCH.  On disordered near-legal
states the contraction channel carries moves with d_ill = -4 or -5: it
annihilates whole defect clusters in one step, which no Pachner move can do
(a 2->3 only ever nibbles).  That is the lever the old VDV/quench route did
not have, and it is why those runs arrested at f_e ~ 0.52.

ACTION

    S = c_N (f3 - f3*)^2                     volume pin
      + nh (f1 - 6 f3 / e*)^2                flat mean-edge-degree pin
      + zleg * sum_v dist^2(n6(v), {0,2,3,4})   hub penalty  (FK-ness)
      + mu_ill * sum_v m(v)                  illegal-edge chemical potential

f0 is left FREE and equilibrates through the contract/split channel; that is
what makes the flat pin satisfiable without forced defects (on a crystal the
composition is frozen and the pin becomes an extensive debt -- see
notes/memory/crystal-library-gas-campaign.md -- but an amorphous state can
just choose its composition).

zleg and mu_ill RAMP linearly over the run: low at the start so the state can
rearrange, high at the end so illegality and hubs are expelled.  Report the
final legality, not the average.

Usage:
    python scripts/fk_amorphous.py --start data/crystal_gas/a15_debt0_nh30.final.mfd \
        --f3-target 5700 --sweeps 8000 --cs 0.25 \
        --zleg 0.0 2.0 --mu-ill 0.5 4.0 --out data/fk_amorph/pilot
"""
import argparse
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


def census(v, cap=10):
    """(n_ill_edges, n_fk, n_hub, n_imp, f0) from the joint (n6, m) census."""
    c = np.asarray(v.valence_census(cap, cap))
    f0 = int(c.sum())
    n_imp = int(c[:, 1:].sum())                 # vertices with an illegal edge
    legal = c[:, 0]
    n_fk = int(sum(legal[k] for k in FK_N6))
    n_hub = f0 - n_imp - n_fk
    n_ill_e = len(v.illegal_edges()[0])
    return n_ill_e, n_fk, n_hub, n_imp, f0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--start", required=True, help="starting .mfd")
    ap.add_argument("--f3-target", type=int, default=0,
                    help="volume pin target (0 = the starting f3)")
    ap.add_argument("--e-target", type=float, default=E_FLAT,
                    help="mean-edge-degree target of the flat pin")
    ap.add_argument("--cn", type=float, default=0.1, help="volume-pin coef")
    ap.add_argument("--nh", type=float, default=30.0, help="flat-pin coef")
    ap.add_argument("--zleg", type=float, nargs=2, default=(0.0, 2.0),
                    metavar=("LO", "HI"), help="hub penalty ramp")
    ap.add_argument("--mu-ill", type=float, nargs=2, default=(0.5, 4.0),
                    metavar=("LO", "HI"), help="illegal-edge fugacity ramp")
    ap.add_argument("--cimp", type=float, default=0.0,
                    help="quadratic clustering term (0 = pure fugacity)")
    ap.add_argument("--ratchet-slack", type=int, default=-1,
                    help="ADAPTIVE budget (-1 = off, overrides --budget): the "
                         "cap is held at (running-minimum illegal-edge count) "
                         "+ slack, so lateral moves are always legal (it "
                         "cannot deadlock) while uphill excursions beyond the "
                         "record are forbidden -- n_ill is monotone "
                         "non-increasing on the coarse scale with NO fugacity "
                         "doing the work. A FIXED cap cannot be used from a "
                         "melt at all: the gate rejects any move whose "
                         "POST-move count exceeds the cap, so at n_ill 1400 "
                         "with cap 20 every move is blocked (sampler.d:632). "
                         "slack is the real dial -- 0 arrests in a glass, "
                         "large is just an uncapped melt.")
    ap.add_argument("--budget", type=int, default=-1,
                    help="hard cap on the illegal-edge count (-1 = off). A "
                         "budget, not a price: keeps a fixed small reservoir "
                         "circulating instead of pricing defects out, which "
                         "freezes the chain. Needs a nonzero n6 coupling.")
    ap.add_argument("--cs", type=float, default=0.25,
                    help="contract/split probability per MCMC step")
    ap.add_argument("--max-ring", type=int, default=6)
    ap.add_argument("--sweeps", type=int, default=8000)
    ap.add_argument("--chunk", type=int, default=100)
    ap.add_argument("--hold", type=float, default=0.25,
                    help="fraction of the run held at the final couplings")
    ap.add_argument("--snap-every", type=int, default=0,
                    help="write <out>.snap.mfd every N sweeps (0 = off). "
                         "Cheap insurance on multi-hour runs: the driver only "
                         "saves .final.mfd at the end, so an interrupt would "
                         "otherwise leave the annealed state unrecoverable "
                         "(the .obs.jsonl trajectory survives either way).")
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
    if args.cs > 0:
        s.set_contract_split(args.cs, args.max_ring)
    v = s.manifold

    os.makedirs(os.path.dirname(os.path.abspath(args.out)) or ".", exist_ok=True)
    meta = {"start": os.path.basename(args.start), "f3_target": f3_target,
            "e_target": args.e_target, "cn": args.cn, "nh": args.nh,
            "zleg": list(args.zleg), "mu_ill": list(args.mu_ill),
            "cimp": args.cimp, "cs": args.cs, "max_ring": args.max_ring,
            "ratchet_slack": args.ratchet_slack, "budget": args.budget,
            "sweeps": args.sweeps, "hold": args.hold, "seed": args.seed,
            "e_flat": E_FLAT}
    rec = Recorder(s, args.out, chunk=args.chunk, census_every=1,
                   snap_mid=False, extra_meta=meta)
    obs = open(args.out + ".obs.jsonl", "w")
    obs.write(json.dumps({"kind": "header", **meta}) + "\n")

    ramp_end = args.sweeps * (1.0 - args.hold)
    t0 = time.time()
    record, cap_now = None, None
    while rec.sw < args.sweeps:
        u = min(1.0, rec.sw / ramp_end) if ramp_end > 0 else 1.0
        zleg = args.zleg[0] + u * (args.zleg[1] - args.zleg[0])
        mu = args.mu_ill[0] + u * (args.mu_ill[1] - args.mu_ill[0])
        s.set_n6_potential(zleg, args.cimp, tilt=[0.0] * 5, imp_lin=mu)
        # the gate needs the vertex-potential state, so (re)apply it AFTER
        # every set_n6_potential call, not once at sweep 0
        if args.ratchet_slack >= 0:
            n_now = len(v.illegal_edges()[0])
            record = n_now if record is None else min(record, n_now)
            cap_now = record + args.ratchet_slack
            s.set_illegal_budget(cap_now)
        elif args.budget >= 0:
            s.set_illegal_budget(args.budget)
        rec.step()

        n_ill_e, n_fk, n_hub, n_imp, f0 = census(v)
        f1, f3 = int(v.f_vector[1]), int(v.f_vector[3])
        row = {"sw": rec.sw, "zleg": zleg, "mu_ill": mu,
               "f0": f0, "f1": f1, "f3": f3,
               "e_mean": 6.0 * f3 / f1, "zbar": 2.0 * f1 / f0,
               "n_ill_e": n_ill_e, "n_fk": n_fk, "n_hub": n_hub,
               "n_imp": n_imp,
               "gap": f1 - 6.0 * f3 / args.e_target,
               "cs": s.contract_split_stats(),
               "budget": s.illegal_budget_stats(),
               "ratchet_record": record, "cap": cap_now,
               # the user-facing f_FK: EVERY incident edge deg 5 or 6, hubs
               # allowed (CONVENTIONS `legalvert`); n_fk is the stricter
               # `legalvert_fk` that also demands n6 in {0,2,3,4}
               "f_legalvert": (f0 - n_imp) / f0}
        obs.write(json.dumps(row) + "\n")
        obs.flush()
        if args.snap_every and rec.sw % args.snap_every == 0:
            v.save(args.out + ".snap.mfd",
                   [f"periodic snapshot at sweep {rec.sw}",
                    f"n_ill_e={n_ill_e} f_legalvert={(f0 - n_imp) / f0:.4f}"])
        print(f"sw {rec.sw:6d} zleg {zleg:4.2f} mu {mu:4.2f} | "
              f"f0 {f0:5d} f3 {f3:6d} e {row['e_mean']:.5f} "
              f"Zbar {row['zbar']:.4f} | n_ill_e {n_ill_e:5d} "
              f"hubs {n_hub:4d} imp {n_imp:5d} "
              f"f_legalvert {100.0 * (f0 - n_imp) / f0:6.2f}% "
              f"(fk-coord {100.0 * n_fk / f0:6.2f}%) "
              f"cap {cap_now}  [{time.time() - t0:.0f}s]",
              flush=True)

    v.save(args.out + ".final.mfd")
    obs.close()
    print(f"\nfinal state: {os.path.abspath(args.out)}.final.mfd")


if __name__ == "__main__":
    main()
