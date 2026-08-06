#!/usr/bin/env python3
"""Edge-contraction move: delete two adjacent stars, cone the boundary.

Aaron's donor-free graft: for adjacent vertices (u, v), remove the union of
their closed stars and glue back the cone over the boundary sphere from ONE
new vertex w.  This is exactly EDGE CONTRACTION (u,v) -> w: tets containing
both u and v contribute no boundary faces, every one-vertex tet survives with
u/v -> w.  Validity = the classical link condition, which is precisely the
surface-pinch check in graft_signature.lump_boundary.

Cost structure (Frank's rule collects a toll, exactly computable):
  * CN(w) = CN(u) + CN(v) - 2 - deg(uv)   (17 for Z12-Z12 along a 5-edge,
    up to 24 for Z16-Z16 along a 6-edge) -- w is always non-FK;
  * every edge of the link(uv) ring loses one tet (its d_in drops 3 -> 2),
    so deg-5 ring edges become illegal deg-4;
  * dV = -1 (more if extra vertices were swallowed), df3 = |bd faces| - |lump|.

The move digs through the wall that isolates FK states under single Pachner
moves; the sampler's job is then to dissolve the localized defect cluster at
the new, lower vertex count.

Modes:
    --demo    contract one Z12-Z12 edge of pristine c15 m3, report the exact
              defect cluster and validation.
    --relax   A/B experiment: inflate c15 m3 by k 1->4 insertions,
              pre-thermalize under the minimal gas action, then relax with
              (A) the plain sampler vs (B) sampler + Metropolis-screened
              contraction attempts.  Tracks f0/f1/f3/n_ill/sum m^2.

NOTE (prototype): branch B is a composite chain WITHOUT the reverse (vertex
split) channel, so it is an ANNEALING scheme, not an equilibrium sampler.
The custom .jsonl trajectory log here deviates from the Recorder convention
deliberately: samplers are re-created per cycle around the Python surgery.
A D-native reversible contract/split channel is the production follow-up.

Usage:
    caffeinate -i python scripts/contract_relax.py --demo
    caffeinate -i python scripts/contract_relax.py --relax --inflate 12 \
        --cimp 0.5 --cycles 40 --sweeps 20 --attempts 200
"""
import argparse
import json
import math
import os
import sys
import time

import numpy as np

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_ROOT, "python"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import discrete_differential_geometry as ddg
from graft_signature import (CrystalContext, SurfaceError, lump_boundary,
                             validate_facets)

E_FLAT = 2.0 * math.pi / math.acos(1.0 / 3.0)


# ------------------------------------------------------------- the move

def contraction_delta(ctx, u, v):
    """Exact consequences of contracting edge (u,v); no mutation.

    Returns None if the contraction is invalid (pinched boundary, extra
    swallowed interior vertices are allowed and reported).  Otherwise a dict
    with the lump, boundary, per-vertex new impurity valences, and the
    action-relevant deltas (df0, df1, df3, dsum_m2).
    """
    F = ctx.F
    lump_idx = np.nonzero(((F == u) | (F == v)).any(axis=1))[0]
    try:
        bd = lump_boundary(ctx, lump_idx)
    except SurfaceError:
        return None
    faces = bd["faces"]
    sverts = sorted({x for f in faces for x in f})
    n_s = len(sverts)
    # boundary must be a single S^2
    n_se = len(bd["edeco"])
    if n_s - n_se + len(faces) != 2:
        return None

    # new degrees: (w,x) = # boundary faces at x; surface edge -> d_out + 2
    wdeg = {}
    for f in faces:
        for x in f:
            wdeg[x] = wdeg.get(x, 0) + 1
    sedge_new = {e: (deg - din + 2) for e, (deg, din) in bd["edeco"].items()}

    ill = lambda d: d not in (5, 6)
    dsum_m2 = 0.0
    m_new = {}
    for x in sverts:
        m_out = int(ctx.m[x])
        for e, (deg, din) in bd["edeco"].items():
            if x in e and ill(deg):
                m_out -= 1                       # old surface edge, changes
        for e, din in bd["d_in"].items():
            if x in e and din == ctx.edge_deg[e] and ill(ctx.edge_deg[e]):
                m_out -= 1                       # interior edge, vanishes
        mx = m_out
        for e, dnew in sedge_new.items():
            if x in e and ill(dnew):
                mx += 1
        if ill(wdeg[x]):
            mx += 1
        m_new[x] = mx
        dsum_m2 += mx * mx - int(ctx.m[x]) ** 2
    mw = sum(1 for x in sverts if ill(wdeg[x]))
    dsum_m2 += mw * mw
    for x in bd["interior_vids"]:
        dsum_m2 -= int(ctx.m[x]) ** 2

    interior = bd["interior_vids"]
    df0 = 1 - len(interior)
    df3 = len(faces) - len(lump_idx)
    old_edges = sum(1 for e, din in bd["d_in"].items()
                    if din == ctx.edge_deg[e])   # interior edges vanish
    df1 = n_s - old_edges
    return dict(lump_idx=lump_idx, faces=faces, sverts=sverts, wdeg=wdeg,
                m_new=m_new, m_w=mw, df0=df0, df1=df1, df3=df3,
                dsum_m2=dsum_m2, interior=sorted(interior),
                cn_w=n_s, sedge_new=sedge_new)


def apply_contraction(F, delta, w):
    """Facet array after the contraction (w = fresh vertex label)."""
    keep = np.ones(len(F), bool)
    keep[delta["lump_idx"]] = False
    cone = np.array([sorted((w,) + f) for f in delta["faces"]], np.int64)
    newF = np.vstack([F[keep], cone])
    lab, inv = np.unique(newF, return_inverse=True)
    return inv.reshape(newF.shape).astype(np.int64)


# ------------------------------------------------------------------ modes

def action(f3, f1, sum_m2, f3_ref, cimp):
    return (0.1 * (f3 - f3_ref) ** 2
            + 1.0 * (f1 - 6.0 * f3 / E_FLAT) ** 2
            + cimp * sum_m2)


def state_terms(ctx):
    f3 = len(ctx.F)
    f1 = len(ctx.edge_deg)
    sum_m2 = float((ctx.m.astype(float) ** 2).sum())
    n_ill = int((ctx.m > 0).sum())
    return f3, f1, sum_m2, n_ill


def demo():
    path = os.path.join(_ROOT, "data", "tcp_reference", "T3_C15_m3_N3672.mfd")
    F = np.asarray(ddg.Manifold.load(path, 3).facets(), np.int64)
    ctx = CrystalContext(F, "c15m3")
    # a deg-5 edge between two Z12s
    edge = next(e for e, d in ctx.edge_deg.items()
                if d == 5 and ctx.cn[e[0]] == 12 and ctx.cn[e[1]] == 12)
    u, v = edge
    d = contraction_delta(ctx, u, v)
    assert d is not None, "contraction invalid on pristine crystal?"
    ring = sorted(d["sedge_new"].values())
    print(f"contracting Z12-Z12 edge {edge} (deg 5):")
    print(f"  lump {len(d['lump_idx'])} tets, boundary {len(d['faces'])} "
          f"faces, interior {d['interior']}")
    print(f"  CN(w) = {d['cn_w']}, m(w) = {d['m_w']}, "
          f"w-edge degrees: {sorted(d['wdeg'].values())}")
    print(f"  df0={d['df0']} df1={d['df1']} df3={d['df3']} "
          f"dsum_m2={d['dsum_m2']:.0f}")
    newF = apply_contraction(F, d, int(F.max()) + 1)
    rep = validate_facets(newF)
    print(f"  validation: V={rep['n_vertices']} degs={rep['edge_degrees']} "
          f"Z={rep['z_census']} broken={rep['n_broken_disclination']} "
          f"chi={rep['euler_characteristic']} orient={rep['orientable']}")
    ctx2 = CrystalContext(newF)
    f3, f1, sm2, n_ill = state_terms(ctx2)
    f30, f10, sm20, _ = state_terms(ctx)
    print(f"  exact-vs-predicted: df1 {f1-f10} vs {d['df1']}, "
          f"df3 {f3-f30} vs {d['df3']}, dsum_m2 {sm2-sm20:.0f} "
          f"vs {d['dsum_m2']:.0f}")


def make_sampler(mfd, f3_target, cimp):
    params = ddg.SamplerParams(
        num_facets_target=f3_target, num_facets_coef=0.1,
        hinge_degree_target=E_FLAT, num_hinges_coef=1.0,
        hinge_degree_variance_coef=0.0, codim3_degree_variance_coef=0.0,
        hinge_degree_target_coef=0.0, codim3_degree_target_coef=0.0)
    s = ddg.ManifoldSampler(mfd, params)
    s.set_n6_potential(0.0, cimp, tilt=[0.0] * 5)
    return s


def relax(args):
    rng = np.random.default_rng(args.seed)
    path = os.path.join(_ROOT, "data", "tcp_reference", "T3_C15_m3_N3672.mfd")
    mfd = ddg.Manifold.load(path, 3)
    f3_ref = mfd.num_facets
    f0_ref = int(mfd.f_vector[0])

    # preparation: equilibrate against an INFLATED volume target at cheap
    # c_imp -- the equilibrium there carries extra vertices.  Restoring the
    # true target at the quench leaves a state with too many vertices by
    # construction, glassy at high c_imp.
    f3_prep = f3_ref + args.prep_extra
    s = make_sampler(mfd, f3_prep, args.cimp_pre)
    for i in range(0, args.prethermal, 10):
        s.run(sweeps=10)
        v = s.manifold
        print(f"  prep sweep {i+10} (f3_target={f3_prep}, "
              f"cimp={args.cimp_pre}): f0={int(v.f_vector[0])} "
              f"f3={v.num_facets}")
    F0 = np.asarray(s.manifold.dup().facets(), np.int64)
    ctx0 = CrystalContext(F0)
    f0_start = len(np.unique(F0))
    print(f"prepared: f0={f0_start} (crystal {f0_ref}) f3={len(F0)} "
          f"n_ill={int((ctx0.m > 0).sum())}; quench to f3_target={f3_ref}, "
          f"cimp={args.cimp}")
    if f0_start <= f0_ref:
        print("WARNING: nothing left to relax -- raise --prep-extra")

    logs = {}
    for branch in ("A", "B"):
        F = F0.copy()
        log = []
        t0 = time.time()
        n_contract = 0
        for cyc in range(args.cycles):
            mfd = ddg.Manifold(3, F.tolist())
            s = make_sampler(mfd, f3_ref, args.cimp)
            s.run(sweeps=args.sweeps)
            F = np.asarray(s.manifold.dup().facets(), np.int64)
            ctx = CrystalContext(F)
            acc = 0
            if branch == "B":
                f3, f1, sm2, _ = state_terms(ctx)
                edges = list(ctx.edge_deg.keys())
                for _ in range(args.attempts):
                    if len(np.unique(F)) <= f0_ref:
                        break              # overshoot guard: at crystal f0
                    u, v = edges[rng.integers(len(edges))]
                    d = contraction_delta(ctx, int(u), int(v))
                    if d is None:
                        continue
                    dS = (action(f3 + d["df3"], f1 + d["df1"],
                                 sm2 + d["dsum_m2"], f3_ref, args.cimp)
                          - action(f3, f1, sm2, f3_ref, args.cimp))
                    if dS <= 0 or rng.random() < math.exp(-dS / args.temp):
                        F = apply_contraction(F, d, int(F.max()) + 1)
                        ctx = CrystalContext(F)
                        f3, f1, sm2, _ = state_terms(ctx)
                        edges = list(ctx.edge_deg.keys())
                        acc += 1
                        n_contract += 1
            f3, f1, sm2, n_ill = state_terms(ctx)
            f0 = len(np.unique(F))
            log.append(dict(cycle=cyc, f0=f0, f1=f1, f3=f3,
                            sum_m2=sm2, n_ill=n_ill,
                            S=action(f3, f1, sm2, f3_ref, args.cimp),
                            contracted=acc))
            if cyc % 5 == 0 or cyc == args.cycles - 1:
                print(f"  [{branch}] cyc {cyc:3d}: f0={f0} f3={f3} "
                      f"n_ill={n_ill} S={log[-1]['S']:.1f} "
                      f"(+{acc} contractions)")
        logs[branch] = log
        print(f"[{branch}] done: f0 {log[0]['f0']} -> {log[-1]['f0']} "
              f"(crystal {f0_ref}), S -> {log[-1]['S']:.1f}, "
              f"{n_contract} contractions, {time.time()-t0:.0f}s")

    out = os.path.join(_ROOT, "data", "grafts",
                       f"contract_relax_k{args.inflate}_c{args.cimp}.json")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    with open(out, "w") as fh:
        json.dump(dict(f0_ref=f0_ref, f3_ref=f3_ref,
                       args=vars(args), logs=logs), fh, indent=1)
    print(f"log -> {out}")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--demo", action="store_true")
    ap.add_argument("--relax", action="store_true")
    ap.add_argument("--inflate", type=int, default=12)   # (unused in v3 prep)
    ap.add_argument("--prep-extra", type=int, default=150)
    ap.add_argument("--cimp", type=float, default=0.5)
    ap.add_argument("--cimp-pre", type=float, default=0.2)
    ap.add_argument("--prethermal", type=int, default=50)
    ap.add_argument("--cycles", type=int, default=40)
    ap.add_argument("--sweeps", type=int, default=20)
    ap.add_argument("--attempts", type=int, default=200)
    ap.add_argument("--temp", type=float, default=1.0)
    ap.add_argument("--seed", type=int, default=7)
    args = ap.parse_args()
    if args.demo:
        demo()
    if args.relax:
        relax(args)


if __name__ == "__main__":
    main()
