#!/usr/bin/env python3
"""Do larger, longer-lived defect complexes travel farther?

Per-WORLDLINE transport: each defect complex is tracked at accepted-move
resolution (defect_state.Worldlines) and its centroid is sampled through the
D-side incremental cocycle lift (O(1) per query, gauge fixed for the run, so
displacements need no glitch filtering). Channel: complex-centroid tracer.

Per track this records lifetime, size, and the centroid trajectory, and
reports net displacement, maximum excursion, and an effective diffusivity.

TWO CONFOUNDS, stated up front:

1. Lifetime itself. Under any diffusive motion net^2 ~ t, so long-lived
   complexes travel farther trivially. "Travels farther" only means something
   at fixed time: the comparison that matters is D_eff = net^2 / (6 * life)
   across size classes, plus the scaling exponent of net^2 vs life.

2. Membership churn. The centroid moves when a vertex joins or leaves even if
   nothing is transported (the run5h centroid-protocol lesson). Bigger
   complexes rearrange more (mode_frac 0.38 vs 0.87), so their centroids
   jitter more per unit time regardless of real motion. Net displacement over
   a long life is robust against this (jitter averages, transport
   accumulates); per-step path length is NOT and is reported only as the
   jitter scale, never as "distance travelled".
"""
import argparse
import json
import os
import sys
from collections import Counter

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
for _p in ("../../python", "../../scripts", "../../tools", "."):
    sys.path.insert(0, os.path.join(_HERE, _p))
import discrete_differential_geometry as ddg
from discrete_differential_geometry import cocycle as coc
import defect_state as ds
from fk_skeleton import edges_from_facets

CELL = 1.0e6


class Trav:
    __slots__ = ("born", "last", "states", "samples", "nv_sum", "nv_n")

    def __init__(self, clock):
        self.born = clock
        self.last = clock
        self.states = Counter()      # (shape f, sumZ, n_verts) -> steps
        self.samples = []            # (clock, centroid in cells)
        self.nv_sum = 0
        self.nv_n = 0


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--cell", required=True)
    ap.add_argument("--start", default=None)
    ap.add_argument("--startcoc", default=None)
    ap.add_argument("--mcell", type=int, default=4)
    ap.add_argument("--lam", type=float, required=True)
    ap.add_argument("--etarget", type=float, default=None)
    ap.add_argument("--zleg", type=float, default=0.6)
    ap.add_argument("--cimp", type=float, default=1.0)
    ap.add_argument("--slide-prob", type=float, default=0.0)
    ap.add_argument("--sweeps", type=int, default=2500)
    ap.add_argument("--chunk", type=int, default=25)
    ap.add_argument("--stride", type=int, default=2,
                    help="sample centroids every Nth accepted move")
    ap.add_argument("--logmb", type=float, default=256.0)
    ap.add_argument("--audit-every", type=int, default=8)
    ap.add_argument("--seed", type=int, default=4242)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    ddg.set_random_seed(args.seed)
    ref = ddg.Manifold.load(args.cell, 3)
    native = float(edges_from_facets(ref.facets())[1].mean())
    et = args.etarget if args.etarget is not None else native
    m = ddg.Manifold.load(args.start or args.cell, 3)
    params = ddg.SamplerParams(
        num_facets_target=ref.num_facets, num_facets_coef=0.1,
        hinge_degree_target=et, num_hinges_coef=0.0,
        hinge_degree_variance_coef=0.0, codim3_degree_variance_coef=0.0,
        hinge_degree_target_coef=args.lam * et / 6.0)
    s = ddg.ManifoldSampler(m, params)
    if args.zleg or args.cimp:
        s.set_n6_potential(args.zleg * args.lam, args.cimp * args.lam,
                           tilt=[0.0] * 5)
    if args.slide_prob:
        s.set_slide_prob(args.slide_prob)
    v = s.manifold
    if not args.startcoc:
        sys.exit("--startcoc is required (positions need a cocycle)")
    e0, w0, _ = coc.load_cocycle(args.startcoc)
    s.enable_cocycle(np.asarray(e0), np.asarray(w0))
    s.check_cocycle()
    pbox = CELL * args.mcell
    s.enable_cocycle_positions([int(pbox)] * 3)
    s.check_cocycle_positions()
    n3 = v.num_facets

    st = ds.DefectState(v)
    wl = ds.Worldlines(st)
    s.enable_event_log(args.logmb)
    s.drain_event_log()

    def centroid(verts):
        q = s.vertex_positions(verts).astype(float)
        d = (q - q[0] + pbox / 2) % pbox - pbox / 2
        return (q[0] + d.mean(0)) / CELL

    tracks = {}
    nev = 0

    def record(clock):
        for cx, tid in wl.step(clock):
            T = tracks.get(tid)
            if T is None:
                T = tracks[tid] = Trav(clock)
                T.samples.append((clock, centroid(cx.verts)))
            T.last = clock
            T.nv_sum += len(cx.verts)
            T.nv_n += 1
            T.states[(tuple(st.induced_shape(cx.verts)["f"]),
                      st.total_coordination(cx.verts), len(cx.verts))] += 1
            if nev % args.stride == 0:
                T.samples.append((clock, centroid(cx.verts)))

    record(0)
    done = 0
    audit_fail = []
    while done < args.sweeps:
        s.run(sweeps=args.chunk)
        done += args.chunk
        ev = s.drain_event_log()
        if s.event_log_overflowed():
            print("  WARN event-log overflow", flush=True)
        for e in ev:
            st.apply(e)
            nev += 1
            record(int(e["clock"]))
        if args.audit_every and (done // args.chunk) % args.audit_every == 0:
            p = st.audit(v)
            if p:
                audit_fail.append((done, p[:2]))
                print(f"  AUDIT FAILED sw{done}: {p[:2]}", flush=True)
            s.check_cocycle_positions()
        print(f"  sw{done}: {len(tracks)} worldlines, "
              f"{len(st.defect)} defect vertices", flush=True)

    rows = []
    for tid, T in tracks.items():
        if len(T.samples) < 2:
            continue
        cl = np.array([x[0] for x in T.samples], float) / n3
        raw = np.array([x[1] for x in T.samples])
        un = np.zeros_like(raw)
        for i in range(1, len(raw)):
            step = raw[i] - raw[i - 1]
            step -= np.round(step / args.mcell) * args.mcell
            un[i] = un[i - 1] + step
        disp = un - un[0]
        r2 = np.sum(disp ** 2, axis=1)
        # jitter scale: median per-sample step, NOT a distance travelled
        steps = np.sqrt(np.sum(np.diff(un, axis=0) ** 2, axis=1))
        (fvec, sumz, nv), _ = T.states.most_common(1)[0]
        life = (T.last - T.born) / n3
        rows.append({
            "tid": tid, "life_sweeps": life,
            "censored": tid not in wl.died,
            "n_verts_mean": T.nv_sum / T.nv_n,
            "n_verts_mode": nv, "sumZ": sumz, "f": list(fvec),
            "pure_shape": len({k[0] for k in T.states}) == 1,
            "n_samples": len(T.samples),
            "net2": float(r2[-1]),                  # cells^2, birth -> death
            "maxexc2": float(r2.max()),
            "jitter_step": float(np.median(steps)) if len(steps) else 0.0,
        })

    with open(args.out + ".json", "w") as fh:
        json.dump({"sweeps": done, "n3": n3, "mcell": args.mcell,
                   "lam": args.lam, "slide_prob": args.slide_prob,
                   "stride": args.stride, "seed": args.seed,
                   "audit_failures": audit_fail, "tracks": rows}, fh)

    obs = [r for r in rows if not r["censored"]]
    print(f"\n{len(rows)} worldlines with >=2 samples, {len(obs)} died")
    if audit_fail:
        print(f"AUDIT FAILURES: {audit_fail}")

    print(f"\nby size class (modal n_verts), died tracks only "
          f"[cells, sweeps]:")
    print(f"   {'nv':>3} {'n':>5} {'med life':>9} {'med net':>8} "
          f"{'med maxexc':>10} {'med D_eff':>10} {'jitter':>7}")
    for nv in sorted({r["n_verts_mode"] for r in obs}):
        g = [r for r in obs if r["n_verts_mode"] == nv]
        if len(g) < 5:
            continue
        net = np.sqrt([r["net2"] for r in g])
        exc = np.sqrt([r["maxexc2"] for r in g])
        lf = np.array([r["life_sweeps"] for r in g])
        deff = np.array([r["net2"] / (6 * max(r["life_sweeps"], 1e-9))
                         for r in g])
        jit = np.median([r["jitter_step"] for r in g])
        print(f"   {nv:3d} {len(g):5d} {np.median(lf):9.2f} "
              f"{np.median(net):8.3f} {np.median(exc):10.3f} "
              f"{np.median(deff):10.4f} {jit:7.4f}")

    # scaling: net2 vs life, pooled and per coarse size class
    print(f"\nnet^2 vs lifetime scaling (log-log OLS on died tracks with "
          f"life > 0.5 sw):")
    for lab, sel in (("all", obs),
                     ("small (nv <= 6)",
                      [r for r in obs if r["n_verts_mode"] <= 6]),
                     ("large (nv >= 7)",
                      [r for r in obs if r["n_verts_mode"] >= 7])):
        g = [r for r in sel if r["life_sweeps"] > 0.5 and r["net2"] > 0]
        if len(g) < 10:
            print(f"   {lab:18s} n={len(g)} (too few)")
            continue
        x = np.log([r["life_sweeps"] for r in g])
        y = np.log([r["net2"] for r in g])
        A = np.vstack([x, np.ones_like(x)]).T
        (a, b), res, *_ = np.linalg.lstsq(A, y, rcond=None)
        # block bootstrap over tracks for the slope error
        rng = np.random.default_rng(5)
        bs = []
        for _ in range(2000):
            i = rng.integers(0, len(g), len(g))
            bs.append(np.linalg.lstsq(A[i], y[i], rcond=None)[0][0])
        print(f"   {lab:18s} n={len(g):5d}  exponent {a:+.3f} "
              f"+- {np.std(bs):.3f}   (1 = diffusive, 0 = pinned)")
    print(f"\nwrote {args.out}.json")


if __name__ == "__main__":
    main()
