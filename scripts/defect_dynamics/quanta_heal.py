#!/usr/bin/env python3
"""Strain a crystal by a FEW edge-degree quanta, then release it and watch it heal.

Two phases in ONE chain, on one crystal from the gas-campaign LIBRARY:

  A (STRAIN)  pins at (f3_ref, e_tar) = (f3_nat - dq, e(f3_nat - dq))
  B (HEAL)    pins snapped back to native (f3_nat, e_nat); same chain, same
              couplings -- only the two targets move.

WHAT ONE QUANTUM IS

At fixed vertex count a closed 3-manifold has f1 = f0 + f3, so the mean edge
degree is a function of f3 alone,

    e(f3) = 6 f3 / (f0 + f3),

and its achievable values form a LATTICE with spacing 6 f0 / f1^2 -- one tet.
That lattice spacing IS the quantum of mean edge degree, and it is the same
quantum the gas campaign counts as one "forced 2->3-equivalent move" (for a15
m5 the flat target e* sits 51.7 quanta below native, vs the campaign's 51.3
forced moves per 1000 vertices).

Both pins are therefore zeroed SIMULTANEOUSLY at (f0_nat, f3_nat - dq) --
that is the point of moving the volume target along with the degree target:
the strained state needs no change in vertex count, so the flat pin acts as a
pure restoring force on f0 (each 1<->4 costs num_hinges_coef * 0.478^2) while
the two targets jointly select f3.

WHAT THE STRAINED TARGET ASKS FOR, STRUCTURALLY

Solve n5 + n6 = f1 and 5 n5 + 6 n6 = 6 f3 at f0 fixed: an EDGE-LEGAL state
at dq quanta below native has

    n5 = f0 * 6      (unchanged!)          n6 = f3_nat - f0*5 - dq

i.e. exactly dq fewer sixfold edges and the SAME 6000 fivefold disclination
lines. Straining by dq quanta = asking the six-web to shed dq units of line
length while the fivefold web is untouched. A perfectly edge-legal answer
exists arithmetically at every dq; whether the sampler can find one -- and
whether it can find its way BACK when the target is restored -- is the
experiment. n_ill = 0 at the end of phase B with the native f-vector restored
and the crystal_grains registry intact means it healed to the pristine
crystal; n_ill = 0 with a broken registry means it healed to a DIFFERENT
edge-legal state.

ACTION (the campaign's minimal three-term one, nothing else on)

    S = A (f3 - f3_ref)^2 + B (f1 - 6 f3 / e_tar)^2 + c_imp sum_v m(v)^2

At fixed f0 the first two terms are the same quadratic in f3, with effective
stiffness c_eff = A + B (f0/f3_ref)^2, so a pristine crystal at dq quanta off
target pays c_eff*dq^2 and the first (necessarily uphill) 2->3 costs
c_eff*(2 dq + 1) on top of its m^2 price. B is thus doing double duty:
enforcement of the strain AND the pin on f0. Below B ~ 10 the vertex count is
not held and the premise of the run dissolves; far above B ~ 30 nothing can
nucleate. Both are printed at startup so a cell can be judged before it runs.

Usage:
    python scripts/defect_dynamics/quanta_heal.py --crystal a15 --dq 3 \
        --cimp 1.0 --hinge-coef 30 --strain 8000 --heal 8000 \
        --out data/quanta_heal/a15_dq3_nh30_c1.0
"""
import argparse
import json
import math
import os
import sys
import time

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
from crystal_gas_scan import LIBRARY, E_FLAT, defect_stats


def edge_degree_target(f0: int, f3: int) -> float:
    """Mean edge degree of a closed 3-manifold with this f0 and f3."""
    return 6.0 * f3 / (f0 + f3)


def edge_hist(view, lo: int = 3, hi: int = 9) -> dict:
    """Full edge-degree histogram (not just the illegal ones).

    n6 is the six-web line density -- the quantity the strain actually acts
    on -- so it is worth the O(E) call every chunk."""
    h = np.asarray(view.degree_histogram(1))     # h[i] = count of degree i+1
    out = {d: int(h[d - 1]) for d in range(lo, hi + 1)
           if d - 1 < len(h) and h[d - 1]}
    tail = int(h[hi:].sum()) if len(h) > hi else 0
    if tail:
        out[f">{hi}"] = tail
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--crystal", default="a15", choices=sorted(LIBRARY))
    ap.add_argument("--dq", type=int, default=3,
                    help="quanta BELOW native for the phase-A targets "
                         "(1 quantum = 1 tet = 1 sixfold edge at fixed f0); "
                         "negative strains upward. 0 = pure control run.")
    ap.add_argument("--cimp", type=float, default=1.0)
    ap.add_argument("--heal-cimp", type=float, default=None,
                    help="c_imp for phase B (default: same as phase A). "
                         "Raising it quenches: defects get more expensive at "
                         "the same time as the targets return to native, so "
                         "'heal' means anneal-by-cooling rather than pure "
                         "release. Reported separately -- it moves two knobs.")
    ap.add_argument("--hinge-coef", type=float, default=30.0,
                    help="num_hinges_coef B: flat-pin stiffness, and the only "
                         "restoring force on f0")
    ap.add_argument("--vol-coef", type=float, default=0.1,
                    help="num_facets_coef A (campaign standard 0.1)")
    ap.add_argument("--strain", type=int, default=8000, help="phase-A sweeps")
    ap.add_argument("--heal", type=int, default=8000, help="phase-B sweeps")
    ap.add_argument("--heal-dq", type=int, default=0,
                    help="quanta below native for the PHASE-B targets "
                         "(default 0 = release to native, the healing test). "
                         "Set it equal to --dq for an ISOSTRAIN quench: the "
                         "strain is held while --heal-cimp raises the defect "
                         "price, which asks whether an edge-legal state "
                         "exists at dq quanta below native and can be "
                         "annealed into, rather than whether the crystal can "
                         "get back to where it started.")
    ap.add_argument("--ramp-hold", type=int, default=0,
                    help="sweeps to hold at each intermediate quantum while "
                         "walking the phase-A targets down (0 = snap "
                         "straight to dq). Snapping puts the whole strain "
                         "into the FIRST uphill move -- on pristine a15 the "
                         "cheapest 2->3 costs c_eff*(2dq+1) + 8*c_imp, so at "
                         "dq = 3 the chain freezes before it can nucleate. "
                         "Ramping pays c_eff*(2k+1) one quantum at a time.")
    ap.add_argument("--contract-split", type=float, default=0.0,
                    help="probability per step of the contract/split channel "
                         "(0 = off). The ONE known composite move that "
                         "destroys edges in a single accepted step, so it is "
                         "the only candidate route to a below-native strain "
                         "that does not go through the melt -- but it moves "
                         "f0 by one vertex per use, which the fixed-vertex "
                         "premise of this script otherwise forbids.")
    ap.add_argument("--illegal-budget", type=int, default=-1,
                    help="hard cap on the number of illegal edges (-1 = off). "
                         "A BUDGET, not a price: it keeps a fixed small "
                         "reservoir of defects circulating instead of pricing "
                         "them out, which is the only way to be dilute AND "
                         "mobile at the same time. This is the lever for the "
                         "below-native direction, where a price either "
                         "freezes the crystal (c_imp >= 0.5) or lets it melt "
                         "(c_imp <= 0.4). Needs a nonzero c_imp.")
    ap.add_argument("--chunk", type=int, default=25)
    ap.add_argument("--seed", type=int, default=20260810)
    ap.add_argument("--start", default=None,
                    help="start from this .mfd instead of the pristine "
                         "reference (the quanta are still defined by the "
                         "REFERENCE f-vector, so dq keeps its meaning). Use "
                         "it to add an anneal stage on top of a finished "
                         "two-phase run -- e.g. lift the illegal-edge cap "
                         "that made a below-native strain reachable but then "
                         "blocks the last step to legality.")
    ap.add_argument("--cocycle", action="store_true",
                    help="attach the harmonic cocycle at sweep 0 (viewer)")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    fname, mcell = LIBRARY[args.crystal]
    cell = os.path.join(_ROOT, "data", "tcp_reference", fname)

    ddg.set_random_seed(args.seed)
    ref = ddg.Manifold.load(cell, 3)
    f0n, f1n, _, f3n = (int(x) for x in ref.f_vector)
    mfd = ddg.Manifold.load(args.start, 3) if args.start else ref
    if f1n != f0n + f3n:
        raise SystemExit(f"f1 != f0 + f3 ({f1n} != {f0n}+{f3n}): not a closed "
                         "3-manifold with chi = 0; the quantum arithmetic "
                         "this script is built on does not apply")
    e_nat = edge_degree_target(f0n, f3n)
    f3_str = f3n - args.dq
    e_str = edge_degree_target(f0n, f3_str)
    quantum = e_nat - edge_degree_target(f0n, f3n - 1)

    # c_eff and the nucleation cost, both at fixed f0 (see module docstring).
    grad = f0n / f3_str                     # |d gap / d f3|
    c_eff = args.vol_coef + args.hinge_coef * grad ** 2
    n5_leg, n6_leg = 6 * f0n, f3_str - 5 * f0n

    params = ddg.SamplerParams(
        num_facets_target=f3_str, num_facets_coef=args.vol_coef,
        hinge_degree_target=e_str, num_hinges_coef=args.hinge_coef,
        # everything else OFF: the defaults are ON and fight an FK pin
        hinge_degree_variance_coef=0.0, codim3_degree_variance_coef=0.0,
        hinge_degree_target_coef=0.0, codim3_degree_target_coef=0.0)
    s = ddg.ManifoldSampler(mfd, params)
    s.set_n6_potential(0.0, args.cimp, tilt=[0.0] * 5)      # zleg = 0: m^2 only
    if args.contract_split > 0:
        s.set_contract_split(args.contract_split, 6)
    if args.illegal_budget >= 0:
        s.set_illegal_budget(args.illegal_budget)
    v = s.manifold

    if args.cocycle:
        edges = np.asarray(v.simplices(1))
        s.enable_cocycle(edges, coc.build_from_positions(
            edges, reference_frac_positions(args.crystal, mcell), mcell))

    meta = {"crystal": args.crystal, "cell": os.path.basename(cell),
            "mcell": mcell, "dq": args.dq, "cimp": args.cimp,
            "hinge_coef": args.hinge_coef, "vol_coef": args.vol_coef,
            "heal_cimp": args.heal_cimp, "heal_dq": args.heal_dq, "start": args.start, "contract_split": args.contract_split, "illegal_budget": args.illegal_budget, "seed": args.seed, "f0_native": f0n, "f3_native": f3n,
            "f3_strain": f3_str, "e_native": e_nat, "e_strain": e_str,
            "e_flat": E_FLAT, "quantum": quantum,
            "quanta_native_to_eflat": (e_nat - E_FLAT) / quantum,
            "c_eff": c_eff, "S_pristine_strained": c_eff * args.dq ** 2,
            "dS_first_uphill": c_eff * (2 * args.dq + 1),
            "legal_target_n5": n5_leg, "legal_target_n6": n6_leg,
            "strain_sweeps": args.strain, "heal_sweeps": args.heal}

    print(f"[{args.crystal} m{mcell}] native f = ({f0n}, {f1n}, {2*f3n}, "
          f"{f3n}), e_nat = {e_nat:.9f}", flush=True)
    print(f"  quantum (1 tet at fixed f0) = {quantum:.4e}; e* is "
          f"{(e_nat - E_FLAT)/quantum:.2f} quanta below native", flush=True)
    print(f"  phase A: dq = {args.dq} -> f3_ref = {f3_str}, e_tar = "
          f"{e_str:.9f}; edge-legal answer is n5 = {n5_leg}, n6 = {n6_leg} "
          f"(native n6 = {f3n - 5*f0n})", flush=True)
    print(f"  c_eff = {c_eff:.4f} -> pristine pays {c_eff*args.dq**2:.2f}, "
          f"first uphill 2->3 costs {c_eff*(2*args.dq+1):.2f} + m^2; "
          f"1<->4 costs {args.hinge_coef*(1-3*grad)**2:.2f}", flush=True)

    os.makedirs(os.path.dirname(os.path.abspath(args.out)) or ".",
                exist_ok=True)
    obs = open(args.out + ".obs.jsonl", "w")
    obs.write(json.dumps({"kind": "header", **meta}) + "\n")

    t0 = time.time()
    events = []
    prev = None

    def phase(tag, sweeps, f3_ref, e_tar, ramp_hold=0):
        """Run one phase to completion, recording to <out>.<tag>.

        With ramp_hold > 0 the targets walk from native to (f3_ref, e_tar)
        one quantum every ramp_hold sweeps, then hold for the remainder."""
        nonlocal prev
        rec = Recorder(s, f"{args.out}.{tag}", chunk=args.chunk,
                       census_every=1, snap_mid=False,
                       cocycle_box=mcell if args.cocycle else None,
                       extra_meta={**meta, "phase": tag, "f3_ref": f3_ref,
                                   "e_tar": e_tar, "ramp_hold": ramp_hold})
        d = None
        dq_tot = f3n - f3_ref
        k_set, cur_ref, cur_tar = None, f3_ref, e_tar
        while rec.sw < sweeps:
            if ramp_hold and dq_tot:
                k = int(np.sign(dq_tot)) * min(abs(dq_tot),
                                               rec.sw // ramp_hold)
                if k != k_set:
                    cur_ref = f3n - k
                    cur_tar = edge_degree_target(f0n, cur_ref)
                    s.set_num_facets_target(cur_ref)
                    s.set_hinge_degree_target(cur_tar)
                    obs.write(json.dumps({"kind": "ramp", "phase": tag,
                                          "sw": rec.sw, "dq": int(k)}) + "\n")
                    k_set = k
            rec.step(min(args.chunk, sweeps - rec.sw))
            d = defect_stats(v)
            cur = d.pop("verts")
            d["turnover"] = (None if prev is None else
                             (1.0 - len(cur & prev) / len(cur | prev)
                              if (cur | prev) else 0.0))
            prev = cur
            fv = [int(x) for x in v.f_vector]
            d.update(phase=tag, sw=rec.sw, t=round(time.time() - t0, 1),
                     f=fv, e_mean=6.0 * fv[3] / fv[1],
                     dq_now=(f3n - fv[3]), edeg=edge_hist(v),
                     dq_ref=(f3n - cur_ref),
                     cs=(list(s.contract_split_stats())
                         if args.contract_split > 0 else None),
                     budget=(list(s.illegal_budget_stats())
                             if args.illegal_budget >= 0 else None),
                     gap=fv[1] - 6.0 * fv[3] / cur_tar,
                     obj=float(s.current_objective))
            obs.write(json.dumps(d) + "\n")
            obs.flush()
            # first passages: on-target volume, and full edge legality
            if fv[3] == f3_ref and not any(e["what"] == f"{tag}:on_target"
                                           for e in events):
                events.append({"what": f"{tag}:on_target", "sw": rec.sw})
            if d["n_ill"] == 0 and fv[3] == f3_ref and not any(
                    e["what"] == f"{tag}:legal" for e in events):
                events.append({"what": f"{tag}:legal", "sw": rec.sw})
        path = rec.finish()
        print(f"  [{tag}] {rec.sw} sw ({time.time()-t0:.0f}s) "
              f"f3={d['f'][3]} (dq={d['dq_now']}) f0={d['f'][0]} "
              f"n_ill={d['n_ill']} ncomp={d['ncomp']} top1={d['top1']} "
              f"<e>={d['e_mean']:.7f} edeg={d['edeg']} -> {path}", flush=True)
        return d, path

    dA, pathA = phase("A", args.strain, f3_str, e_str, args.ramp_hold)

    f3_heal = f3n - args.heal_dq
    e_heal = edge_degree_target(f0n, f3_heal)
    if args.heal > 0:
        # -- move the targets to the phase-B point; SAME chain --
        s.set_num_facets_target(f3_heal)
        s.set_hinge_degree_target(e_heal)
        if args.heal_cimp is not None and args.heal_cimp != args.cimp:
            s.set_n6_potential(0.0, args.heal_cimp, tilt=[0.0] * 5)
        obs.write(json.dumps({"kind": "release", "sw_A": args.strain,
                              "f3_ref": f3_heal, "e_tar": e_heal,
                              "heal_dq": args.heal_dq, "start": args.start,
                              "cimp": args.heal_cimp
                              if args.heal_cimp is not None
                              else args.cimp}) + "\n")
        dB, pathB = phase("B", args.heal, f3_heal, e_heal)
    else:
        dB, pathB = dA, pathA

    obs.write(json.dumps({"kind": "final", "events": events,
                          "A": {k: dA[k] for k in
                                ("n_ill", "ncomp", "top1", "f", "edeg")},
                          "B": {k: dB[k] for k in
                                ("n_ill", "ncomp", "top1", "f", "edeg")},
                          "pathA": pathA, "pathB": pathB,
                          "t": round(time.time() - t0, 1)}) + "\n")
    obs.close()
    healed = (dB["n_ill"] == 0
              and dB["f"] == [f0n, f0n + f3_heal, 2 * f3_heal, f3_heal])
    print(f"[{args.crystal} dq={args.dq}->{args.heal_dq} nh={args.hinge_coef} "
          f"c={args.cimp}] "
          f"{'EDGE-LEGAL at the phase-B f-vector' if healed else 'not edge-legal'}"
          f"  events={events}", flush=True)


if __name__ == "__main__":
    main()
