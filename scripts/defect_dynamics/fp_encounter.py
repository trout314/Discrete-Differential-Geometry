#!/usr/bin/env python3
"""FP production: encounter kinetics + contact-pair recombination (M3).

Frozen-background (v1 FROZEN) slide-channel kinetics of a mobile (3,4,4)
knot A against a frozen knot B on R m2, on the validated FP driver
(V4/V5 PASS). Time is in attempted-move units (1 sweep = N3 attempts).

--mode encounter   independent episodes from separation --sep-sites:
                   FP-fly until dock (cap --max-flights, censored
                   episodes recorded as such). Per episode: flights,
                   t_dock, dock chord, V(dock) additivity audit, and the
                   dock geometry (on-chain site index if the dock chord
                   is a chain chord of the shared orbit, else off-chain)
                   -- the dock-geometry census.

--mode recombine   fly to dock, then unbind: explicit Metropolis hits
                   (attempted-move clock via geometric skips) until the
                   state leaves the dock zone (cap --unbind-hits), then
                   FP-fly until re-dock (RECOMBINE, t_return recorded)
                   or --escape-flights consecutive dock-free flights
                   (ESCAPE; operational free-pair definition ~ that many
                   ball radii of clear wandering).

Every episode is walked back through exact inverse slides and the
objective audited. One D scan per flight leaks MB-scale (design R6):
keep --episodes modest per process and split seeds across processes.
"""
import argparse
import json
import os
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
for _p in ("../../python", "../../scripts", "../../tools", "."):
    sys.path.insert(0, os.path.join(_HERE, _p))
import discrete_differential_geometry as ddg  # noqa: F401
from discrete_differential_geometry import fpkmc
from discrete_differential_geometry._sampler import SLIDE_SLOTS
from fpkmc_v4_fp import walk_back
from fpkmc_v5_contact import build, R_of_B, root_dock


def fly_to_dock(drv, max_flights):
    for _ in range(max_flights):
        reason, info = drv.flight()
        if reason in (fpkmc.ABSORB_DOCK, fpkmc.ABSORB_MULTI):
            return reason
    return "censored"


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--mode", choices=("encounter", "recombine"),
                    required=True)
    ap.add_argument("--ref", default="data/tcp_reference/T3_R_m2_N7248.mfd")
    ap.add_argument("--lam", type=float, default=0.40)
    ap.add_argument("--estar", type=float, default=5.105025)
    ap.add_argument("--window", type=int, default=40)
    ap.add_argument("--sep-sites", type=int, default=3)
    ap.add_argument("--dS-max", type=float, default=5.0)
    ap.add_argument("--depth", type=int, default=3)
    ap.add_argument("--episodes", type=int, default=4)
    ap.add_argument("--max-flights", type=int, default=60)
    ap.add_argument("--unbind-hits", type=int, default=500)
    ap.add_argument("--escape-flights", type=int, default=8)
    ap.add_argument("--seed", type=int, required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    s, orb, jA, jB, apxA, apxB, faceB, blocked = build(args)
    S_start = s.current_objective
    n3 = s.manifold.num_facets
    nu = fpkmc.nu_per_attempt(1.0, n3)
    Rinf = R_of_B(s, apxB, faceB)
    L = len(orb)
    chain_site = {tuple(sorted((orb[k], orb[(k + 4) % L]))): k
                  for k in range(L)}
    rng = np.random.default_rng(args.seed)
    p_hit = 3.0 / (6.0 * n3)
    episodes = []

    for ep in range(args.episodes):
        drv = fpkmc.FPDriver(s, apxA, blocked, nu, dS_max=args.dS_max,
                             max_depth=args.depth, rng=rng)
        reason = fly_to_dock(drv, args.max_flights)
        rec = {"flights": drv.nflight, "t": drv.t, "reason": str(reason)}
        if reason == fpkmc.ABSORB_DOCK:
            R = R_of_B(s, apxB, faceB)
            rec["V_dock"] = None if R is None else R - Rinf
            rec["dock_chord"] = list(drv.chord)
            k = chain_site.get(drv.chord)
            rec["on_chain_site"] = k
            rec["jB"] = jB

            if args.mode == "recombine":
                # unbind: explicit hits until the state leaves dock zone
                chord = drv.chord
                burst = []
                t_ub = 0
                freed = False
                for _ in range(args.unbind_hits):
                    t_ub += int(rng.geometric(p_hit))
                    slot = int(rng.integers(SLIDE_SLOTS))
                    try:
                        dS, arr = s.slide_at2(chord[0], chord[1], slot,
                                              commit=False)
                    except RuntimeError:
                        continue
                    if dS is None:
                        continue
                    if not (dS <= 0 or rng.random() < np.exp(-dS)):
                        continue
                    dS, arr = s.slide_at2(chord[0], chord[1], slot,
                                          commit=True)
                    new_chord = tuple(sorted(int(x) for x in arr))
                    burst.append((chord, new_chord, float(dS)))
                    chord = new_chord
                    d, nch = root_dock(s, chord, blocked)
                    if nch != 1:
                        break                     # reaction: stop episode
                    if d == 0:
                        freed = True
                        break
                rec.update({"t_unbind": t_ub, "unbind_commits": len(burst),
                            "freed": freed})
                if freed:
                    drv2 = fpkmc.FPDriver(s, chord, blocked, nu,
                                          dS_max=args.dS_max,
                                          max_depth=args.depth, rng=rng)
                    outcome = "escape"
                    clear = 0
                    while clear < args.escape_flights:
                        r2, _ = drv2.flight()
                        if r2 == fpkmc.ABSORB_DOCK:
                            outcome = "recombine"
                            break
                        if r2 == fpkmc.ABSORB_MULTI:
                            outcome = "reaction"
                            break
                        clear += 1
                    rec.update({"outcome": outcome, "t_return": drv2.t,
                                "flights_out": drv2.nflight})
                    walk_back(s, drv2.recs)
                walk_back(s, burst)
        walk_back(s, drv.recs)
        if abs(s.current_objective - S_start) > 1e-6:
            raise RuntimeError(f"episode {ep} restore failed")
        episodes.append(rec)
        print(f"  ep {ep}: {rec}", flush=True)

    with open(args.out, "w") as fh:
        json.dump({"mode": args.mode, "sep_sites": args.sep_sites,
                   "lam": args.lam, "n3": n3, "nu": nu, "seed": args.seed,
                   "max_flights": args.max_flights,
                   "escape_flights": args.escape_flights,
                   "episodes": episodes}, fh)
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
