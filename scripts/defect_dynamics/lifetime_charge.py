#!/usr/bin/env python3
"""Is a defect complex's lifetime correlated with its curvature charge?

Tracks complexes at full accepted-move resolution (defect_state.Worldlines)
and records, per worldline, its lifetime and its charge.

THE CONFOUND, which is why the answer is reported two ways. Charge is very
nearly proportional to size: Q_c = (sum Z)(pi - 3 theta) + 6 n theta, and
sum Z ~ 14 n in this host, so Q_c ~ -0.33 n. Size is already known to predict
lifetime (big complexes persist, small ones blink). So a raw lifetime-vs-charge
correlation mostly re-measures lifetime-vs-size and says nothing new.

The sharp question is whether charge predicts lifetime AT FIXED SHAPE. The
(3,4,4) knot is ideal for it: its induced subcomplex is always the same object
(three tets of the 4-simplex boundary, f = (5,10,9,3)), so size and shape are
held exactly constant while sum Z still varies over a ladder of rungs
(67..72 at lam = 0.40). If lifetime differs across those rungs, charge matters
on its own; if not, the raw correlation was size all along.

CENSORING is handled explicitly rather than ignored: worldlines still alive at
the end of the run have unknown lifetimes and are counted separately. Ignoring
them biases every mean downward, and biases it MORE for long-lived groups --
exactly the groups a lifetime comparison is about.
"""
import argparse
import json
import os
import sys
from collections import Counter, defaultdict

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
for _p in ("../../python", "../../scripts", "../../tools", "."):
    sys.path.insert(0, os.path.join(_HERE, _p))
import discrete_differential_geometry as ddg
from discrete_differential_geometry import cocycle as coc
import defect_state as ds
from cocycle_check import reference_frac_positions
from fk_skeleton import edges_from_facets

THETA = np.arccos(1.0 / 3.0)
UNIT = np.pi - 3 * THETA          # -0.5512856 rad per unit of total Z


def charge(sumz, nverts):
    """Q_c from the closed form -- no need to touch vertex charges."""
    return sumz * UNIT + 6 * nverts * THETA


class Life:
    __slots__ = ("born", "last", "states", "shapes", "steps")

    def __init__(self, clock):
        self.born = clock
        self.last = clock
        # JOINT (shape, sumZ, n_verts) -> steps held. Recording shape and
        # coordination separately and taking each mode independently is wrong:
        # a track can be binned by a shape it held at one time and a
        # coordination it held at another, which produced 5-vertex "knots"
        # with sum Z = 150.
        self.states = Counter()
        self.steps = 0


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
    ap.add_argument("--sweeps", type=int, default=600)
    ap.add_argument("--chunk", type=int, default=25)
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
    if args.startcoc:
        e0, w0, _ = coc.load_cocycle(args.startcoc)
        s.enable_cocycle(np.asarray(e0), np.asarray(w0))
    n3 = v.num_facets

    st = ds.DefectState(v)
    wl = ds.Worldlines(st)
    s.enable_event_log(args.logmb)
    s.drain_event_log()

    lives = {}

    def record(clock):
        for cx, tid in wl.step(clock):
            L = lives.get(tid)
            if L is None:
                L = lives[tid] = Life(clock)
            L.last = clock
            L.steps += 1
            L.states[(tuple(st.induced_shape(cx.verts)["f"]),
                      st.total_coordination(cx.verts), len(cx.verts))] += 1

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
            record(int(e["clock"]))
        if args.audit_every and (done // args.chunk) % args.audit_every == 0:
            p = st.audit(v)
            if p:
                audit_fail.append((done, p[:2]))
                print(f"  AUDIT FAILED sw{done}: {p[:2]}", flush=True)
        print(f"  sw{done}: {len(lives)} worldlines so far, "
              f"{len(st.defect)} defect vertices", flush=True)

    end_clock = max((L.last for L in lives.values()), default=0)
    rows = []
    for tid, L in lives.items():
        (fvec, sumz, nv), nmode = L.states.most_common(1)[0]
        censored = tid not in wl.died          # still alive at the end
        shapes = {k[0] for k in L.states}
        rows.append({
            "tid": tid, "life_sweeps": (L.last - L.born) / n3,
            "censored": bool(censored), "sumZ": sumz, "n_verts": nv,
            "Q_c": charge(sumz, nv), "f": list(fvec), "steps": L.steps,
            # was the track this one shape for its whole life?
            "pure_shape": len(shapes) == 1,
            "mode_frac": nmode / L.steps})

    with open(args.out + ".json", "w") as fh:
        json.dump({"sweeps": done, "n3": n3, "slide_prob": args.slide_prob,
                   "lam": args.lam, "audit_failures": audit_fail,
                   "tracks": rows}, fh)

    obs = [r for r in rows if not r["censored"]]
    cen = [r for r in rows if r["censored"]]
    print(f"\n{len(rows)} worldlines: {len(obs)} died (lifetime observed), "
          f"{len(cen)} still alive at the end (right-censored)")
    if audit_fail:
        print(f"AUDIT FAILURES: {audit_fail}")

    def summarise(sel, label, keyfn):
        groups = defaultdict(list)
        for r in sel:
            groups[keyfn(r)].append(r["life_sweeps"])
        print(f"\n{label}")
        print(f"   {'group':>10} {'n':>5} {'median life':>12} "
              f"{'mean life':>11}  (sweeps)")
        for k in sorted(groups):
            L = groups[k]
            print(f"   {str(k):>10} {len(L):5d} {np.median(L):12.3f} "
                  f"{np.mean(L):11.3f}")

    # 1. raw -- confounded by size
    summarise(obs, "RAW lifetime by charge bin (CONFOUNDED by size: "
                   "Q_c ~ -0.33 * n_verts)",
              lambda r: f"{round(r['Q_c'])}")
    summarise(obs, "lifetime by size (the confound itself)",
              lambda r: r["n_verts"])

    # 2. the clean test: fixed shape, charge varying over the ladder
    def isknot(r):
        # PURE knots only: the same shape for the whole worldline, so size and
        # shape really are held constant and only sum Z varies.
        return tuple(r["f"]) == (5, 10, 9, 3) and r["pure_shape"]
    knots = [r for r in obs if isknot(r)]
    kcen = [r for r in cen if isknot(r)]
    mixed = sum(1 for r in obs if tuple(r["f"]) == (5, 10, 9, 3)
                and not r["pure_shape"])
    print(f"\n=== FIXED SHAPE: the (3,4,4) knot, f = (5,10,9,3) ===")
    print(f"{len(knots)} died, {len(kcen)} censored, {mixed} excluded for "
          f"changing shape during life. Size and shape are held exactly "
          f"constant in what remains; only sum Z (hence Q_c) varies.")
    if knots:
        summarise(knots, "lifetime by sum Z (each rung is a fixed Q_c)",
                  lambda r: r["sumZ"])
        z = np.array([r["sumZ"] for r in knots], float)
        t = np.array([r["life_sweeps"] for r in knots], float)
        if len(set(z)) > 1:
            def rank(a):
                o = np.argsort(a, kind="mergesort")
                r = np.empty(len(a), float)
                r[o] = np.arange(len(a), dtype=float)
                # average ties
                for val in np.unique(a):
                    m_ = a == val
                    r[m_] = r[m_].mean()
                return r
            rho = float(np.corrcoef(z, t)[0, 1])
            srho = float(np.corrcoef(rank(z), rank(t))[0, 1])
            se = 1.0 / np.sqrt(max(len(z) - 3, 1))
            print(f"\n   lifetimes are heavily skewed, so the rank "
                  f"correlation is the one to read:")
            print(f"     Spearman rho(sum Z, lifetime) = {srho:+.3f} "
                  f"({abs(srho)/se:.1f} sigma, n = {len(z)})")
            print(f"     Pearson  r(sum Z, lifetime)   = {rho:+.3f} "
                  f"({abs(rho)/se:.1f} sigma)  -- outlier-sensitive")
            print(f"   sign convention: Q_c falls as sum Z rises, so a "
                  f"NEGATIVE rho here means more-negative charge lives LONGER")
            cfrac = len(kcen) / max(len(kcen) + len(knots), 1)
            print(f"   censored fraction {100*cfrac:.1f}% -- these are the "
                  f"LONGEST-lived and are excluded, so every mean above is a "
                  f"lower bound, more so for longer-lived rungs")
    print(f"\nwrote {args.out}.json")


if __name__ == "__main__":
    main()
