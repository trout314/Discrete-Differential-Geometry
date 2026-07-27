#!/usr/bin/env python3
"""FP production: pure slide-channel transport of a lone knot (M3).

Frozen-background (v1 FROZEN) transport constant: a single (3,4,4) knot
on pristine R m2, no other defects, chained FP flights (no blocked set,
so absorption is only the dS/depth frontier -- protective-domain
transport). Per flight we record the attempted-move clock and the knot
position: the minimal-image-unwrapped midpoint of the chord's two
vertices in crystal coordinates (reference fractional positions x mcell,
units of cells; vertex labels are stable, no 1-4/4-1 in the slide
channel). Trajectories from many seeds/processes pool into MSD(dt) ->
D_slide. Transport channel: SLIDE ONLY, frozen background -- the clean
reference against defect_travel's thermal (flicker-dressed) caging.

A multichord flight exit ends the trajectory (recorded; the lone-knot
sector was left). One D scan per flight leaks MB-scale (R6): keep
--flights per process modest; split across processes.

FRAME: registry lift (extrinsic; CONVENTIONS.md sec 6). Scaling
exponents are frame-robust for near-pristine states
(quasi-isometry); lengths, displacements, D's and radii are
frame-gauge. Gauge-free transport constants: Cartan-development
MSD (fp_intrinsic_msd.py).
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
from cocycle_check import reference_frac_positions
from fpkmc_v3_hb import build


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--ref", default="data/tcp_reference/T3_R_m2_N7248.mfd")
    ap.add_argument("--mcell", type=int, default=2)
    ap.add_argument("--lam", type=float, default=0.40)
    ap.add_argument("--estar", type=float, default=5.105025)
    ap.add_argument("--window", type=int, default=40)
    ap.add_argument("--dS-max", type=float, default=5.0)
    ap.add_argument("--depth", type=int, default=3)
    ap.add_argument("--flights", type=int, default=200)
    ap.add_argument("--seed", type=int, required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    print("[frame] registry lift -- exponents frame-robust, lengths/D gauge; see fp_intrinsic_msd.py")

    m, s, apx, window = build(args.ref, args.estar, args.lam,
                              window=args.window)
    n3 = s.manifold.num_facets
    nu = fpkmc.nu_per_attempt(1.0, n3)
    frac = np.asarray(reference_frac_positions("r", args.mcell))
    box = float(args.mcell)

    def pos_of(chord):
        pa = frac[chord[0]] * box
        pb = frac[chord[1]] * box
        d = pb - pa
        d -= box * np.round(d / box)      # b's image nearest a
        return pa + 0.5 * d

    rng = np.random.default_rng(args.seed)
    drv = fpkmc.FPDriver(s, apx, [], nu, dS_max=args.dS_max,
                         max_depth=args.depth, rng=rng)
    pos = pos_of(drv.chord)
    traj = [{"t": 0, "x": list(pos), "dS": 0.0}]
    ended = None
    for i in range(args.flights):
        reason, info = drv.flight()
        if reason == fpkmc.ABSORB_MULTI:
            ended = "multichord"
            break
        raw = pos_of(drv.chord)
        prev_wrapped = np.mod(pos, box)
        step = raw - prev_wrapped
        step -= box * np.round(step / box)   # minimal-image flight step
        pos = pos + step
        traj.append({"t": int(drv.t), "x": list(pos),
                     "dS": float(drv.S_rel)})
        if (i + 1) % 50 == 0:
            print(f"  flight {i + 1}/{args.flights}: t={drv.t:.3e}, "
                  f"|x-x0|={np.linalg.norm(pos - np.asarray(traj[0]['x'])):.2f} cells",
                  flush=True)
    with open(args.out, "w") as fh:
        json.dump({"lam": args.lam, "n3": n3, "nu": nu, "seed": args.seed,
                   "mcell": args.mcell, "window": args.window,
                   "dS_max": args.dS_max, "depth": args.depth,
                   "ended": ended, "traj": traj}, fh)
    print(f"wrote {args.out} ({len(traj) - 1} flights, ended={ended})")


if __name__ == "__main__":
    main()
