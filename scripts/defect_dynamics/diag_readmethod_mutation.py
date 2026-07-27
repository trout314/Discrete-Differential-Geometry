#!/usr/bin/env python3
"""Prove the DefectState divergence is caused by a read-side query method
mutating state (2026-07-27).

By elimination: a clean apply-only replay stayed exact for 9000 sweeps on
the same seed/stream where reaction_census diverged at 8000. The only
difference is that reaction_census calls components() / induced_shape() /
total_coordination() every step. So a "read-only" method mutates state.

This runs TWO DefectState instances on ONE sampler (identical stream,
same process, fully deterministic):
    A  -- apply() only
    B  -- apply() + the exact read calls reaction_census makes per step
After each chunk, every maintained dict is compared A vs B, and B.edeg is
also checked against a fresh recompute over B.tets (the known-good set).
The first divergence is dumped with the offending key. This localises the
mutating method and the corrupted structure in a few hundred sweeps.
"""
import argparse
import os
import sys
from collections import Counter
from itertools import combinations

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
for _p in ("../../python", "../../scripts", "../../tools", "."):
    sys.path.insert(0, os.path.join(_HERE, _p))
import discrete_differential_geometry as ddg
import defect_state as ds


def snapshot(st):
    return {
        "edeg": {k: v for k, v in st.edeg.items() if v != 0},
        "v2t_keys": set(st.v2t),
        "n6_nonzero": {k: v for k, v in st.n6.items() if v != 0},
        "imp_nonzero": {k: v for k, v in st.imp.items() if v != 0},
        "defect": set(st.defect),
        "ntets": len(st.tets),
    }


def first_diff(a, b, name):
    if a == b:
        return None
    ka, kb = set(a) if isinstance(a, (set, dict)) else a, \
        set(b) if isinstance(b, (set, dict)) else b
    if isinstance(a, dict):
        only_a = {k: a[k] for k in list(set(a) - set(b))[:5]}
        only_b = {k: b[k] for k in list(set(b) - set(a))[:5]}
        val = [(k, a[k], b[k]) for k in (set(a) & set(b))
               if a[k] != b[k]][:5]
        return (f"{name}: A-only={only_a} B-only={only_b} "
                f"val-diff={val}")
    if isinstance(a, set):
        return (f"{name}: A-only={sorted(a - b)[:5]} "
                f"B-only={sorted(b - a)[:5]}")
    return f"{name}: A={a} B={b}"


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--start", default="data/mgas/lam35_snap17000.mfd")
    ap.add_argument("--cell", default="data/tcp_reference/T3_R_m4_N57984.mfd")
    ap.add_argument("--lam", type=float, default=0.35)
    ap.add_argument("--slide-prob", type=float, default=0.10)
    ap.add_argument("--etarget", type=float, default=5.105025)
    ap.add_argument("--seed", type=int, default=1005)
    ap.add_argument("--sweeps", type=int, default=3000)
    ap.add_argument("--chunk", type=int, default=25)
    args = ap.parse_args()

    ddg.set_random_seed(args.seed)
    ref = ddg.Manifold.load(args.cell, 3)
    m = ddg.Manifold.load(args.start, 3)
    et = args.etarget
    params = ddg.SamplerParams(
        num_facets_target=ref.num_facets, num_facets_coef=0.1,
        hinge_degree_target=et, num_hinges_coef=0.0,
        hinge_degree_variance_coef=0.0, codim3_degree_variance_coef=0.0,
        hinge_degree_target_coef=args.lam * et / 6.0)
    s = ddg.ManifoldSampler(m, params)
    s.set_n6_potential(0.6 * args.lam, 1.0 * args.lam, tilt=[0.0] * 5)
    if args.slide_prob:
        s.set_slide_prob(args.slide_prob)
    v = s.manifold

    A = ds.DefectState(v)
    B = ds.DefectState(v)
    s.enable_event_log(64.0)
    s.drain_event_log()

    done = 0
    while done < args.sweeps:
        s.run(sweeps=args.chunk)
        done += args.chunk
        ev = s.drain_event_log()
        for e in ev:
            A.apply(e)
            B.apply(e)
        # B does exactly what reaction_census does per step
        for cx, _ in _fake_worldline_step(B):
            _ = B.induced_shape(cx)
            _ = B.total_coordination(cx)

        sa, sb = snapshot(A), snapshot(B)
        diffs = []
        for key in ("edeg", "v2t_keys", "n6_nonzero", "imp_nonzero",
                    "defect", "ntets"):
            d = first_diff(sa[key], sb[key], key)
            if d:
                diffs.append(d)
        # independent: is each tracker's edeg consistent with its OWN tets?
        # (recompute from the tet set, which the audit proved is always right)
        def edeg_bad(st):
            recomp = Counter()
            for t in st.tets:
                for p in combinations(t, 2):
                    recomp[p] += 1
            return {k: (st.edeg.get(k, 0), recomp.get(k, 0))
                    for k in (set(st.edeg) | set(recomp))
                    if st.edeg.get(k, 0) != recomp.get(k, 0)}
        a_bad = edeg_bad(A)
        b_bad = edeg_bad(B)
        if diffs or a_bad or b_bad:
            print(f"\n=== DIVERGENCE at sweep {done} ===")
            print(f"  A (apply only)      edeg-vs-tets bad edges: {len(a_bad)}")
            print(f"  B (apply + reads)   edeg-vs-tets bad edges: {len(b_bad)}")
            verdict = ("READ METHODS (B corrupts, A clean)"
                       if a_bad == {} and b_bad
                       else "APPLY ITSELF (A corrupts too)"
                       if a_bad else "state A/B mismatch only")
            print(f"  => CULPRIT: {verdict}")
            for d in diffs:
                print("  " + d)
            if b_bad:
                items = list(b_bad.items())[:6]
                print(f"  B.edeg inconsistent with B.tets on "
                      f"{len(b_bad)} edges, e.g. {items}")
            if a_bad:
                items = list(a_bad.items())[:6]
                print(f"  A.edeg inconsistent with A.tets on "
                      f"{len(a_bad)} edges, e.g. {items}")
            # how many v2t keys did B gain that A lacks?
            extra = sb["v2t_keys"] - sa["v2t_keys"]
            print(f"  B has {len(extra)} v2t keys A lacks; sample: "
                  f"{sorted(extra)[:8]}")
            for w in sorted(extra)[:4]:
                print(f"    v2t[{w}] in B = {B.v2t[w]} (empty? "
                      f"{len(B.v2t[w]) == 0})")
            return
        if done % 250 == 0:
            print(f"  sweep {done}: A==B, both consistent "
                  f"(v2t A={len(sa['v2t_keys'])} B={len(sb['v2t_keys'])})",
                  flush=True)
    print(f"\nno divergence in {args.sweeps} sweeps")


def _fake_worldline_step(st):
    """Mimic reaction_census's per-step read pattern: components(), and for
    each, induced_shape + total_coordination on its vertices."""
    return [(c.verts, None) for c in st.components()]


if __name__ == "__main__":
    main()
