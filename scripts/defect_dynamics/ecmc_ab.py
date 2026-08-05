#!/usr/bin/env python3
"""Blob-relaxation A/B: does the lifted ECMC flight channel beat the basic
dynamics on a transport-limited out-of-equilibrium state?

THE STATE. C15 pristine + k fliers (2->3 excitations) packed vertex-adjacent
inside one graph ball -- an m^2-loaded "hot blob". The pins are chosen
self-consistently AFTER insertion (num_facets_target = f3, e* = 6*f3/f1), so
both global pins are EXACTLY zero at t = 0 and at any redistributed k-flier
state: the entire excess action is c_imp * sum_v m(v)^2 from crowding, and
equilibrium is the dilute spread gas. Relaxation therefore requires MOVING
excitations, not making them -- unlike the e-target quench, which relaxes by
local birth/death and cannot be helped by a transport channel (consistent
with the strict-chord-channel null on percolation).

ARMS.
  A  plain sampler (local bistellar + hinge moves).
  B  plain + interleaved lifted ECMC flight episodes (ecmc_flight.Flight,
     expect_free=False since the action is pins + m^2, not EDQ).
  N  plain + the D-side diffusive nonlocal-slide channel (unlifted teleport
     baseline, set_nonlocal_slide_prob).

Honest competition: the soft pins leave the kill-and-rebirth channel open to
every arm, so arm A can "teleport" excitations at its own rate.

WORK ACCOUNTING (the percolation-A/B lesson): totalTried never sees the
flight kernel. Arm B's true work = d_tried + 2 * (nonlocal_slide_at calls)
(each call is one annihilate+recreate composite, trials included). The
sidecar records both raw counters; compare observables vs WORK, not sweeps.

OUTPUT per chain:  <out>/arm<X>_s<NN>.rec.jsonl  (Recorder, unconditional)
                   <out>/arm<X>_s<NN>.ab.jsonl   (sidecar, incremental:
                       work counters, sum m^2, components, spread from the
                       blob center, episode event totals)
                   <out>/manifest.json

    python scripts/defect_dynamics/ecmc_ab.py --arms AB --seeds 4 \
        --k-fliers 12 --cimp 0.5 --sweeps 4000 --chunk 20 \
        --eps-per-chunk 2 --ep-steps 40 --out-dir data/ecmc_ab/pilot1
"""
import argparse
import json
import os
import sys
import time
from collections import Counter, deque
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np

_ROOT = os.path.normpath(os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "..", ".."))
sys.path.insert(0, os.path.join(_ROOT, "python"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

REF = os.path.join(_ROOT, "data/tcp_reference/T3_C15_m3_N3672.mfd")


# ---------------------------------------------------------------------------
# state preparation
# ---------------------------------------------------------------------------

def _face_apex_map(m):
    f2a = {}
    for t in np.asarray(m.facets()):
        t = tuple(sorted(int(x) for x in t))
        for i in range(4):
            f2a.setdefault(t[:i] + t[i + 1:], []).append(t[i])
    return f2a


def _edge_set(m):
    return {tuple(sorted(map(int, e))) for e in m.simplices(1)}


def build_blob(m, k, rng, radius=4):
    """Insert k vertex-adjacent 2->3 fliers inside a graph ball; return
    (list of (face, chord), center vertex)."""
    # ball on the pristine crystal
    edges = _edge_set(m)
    V = int(np.asarray(m.facets()).max()) + 1
    adj = [[] for _ in range(V)]
    for a, b in edges:
        adj[a].append(b)
        adj[b].append(a)
    center = int(rng.integers(V))
    dist = np.full(V, -1, np.int32)
    dist[center] = 0
    dq = deque([center])
    while dq:
        u = dq.popleft()
        if dist[u] >= radius:
            continue
        for w in adj[u]:
            if dist[w] < 0:
                dist[w] = dist[u] + 1
                dq.append(w)
    ball = {v for v in range(V) if 0 <= dist[v] <= radius}

    placed, support = [], set()
    for _ in range(50 * k):                       # retry budget
        if len(placed) >= k:
            break
        f2a = _face_apex_map(m)
        edges = _edge_set(m)
        faces = [(f, ap) for f, ap in f2a.items() if len(ap) == 2]
        rng.shuffle(faces)
        for f, ap in faces:
            sup = set(f) | set(map(int, ap))
            if not sup <= ball:
                continue
            if placed and not (sup & support):
                continue                          # must touch the blob
            chord = tuple(sorted(map(int, ap)))
            if chord in edges:
                continue                          # 2->3 invalid
            m.do_bistellar_move(list(f), list(ap))
            placed.append((f, chord))
            support |= sup
            break
        else:
            break                                 # no candidate face at all
    if len(placed) < k:
        raise RuntimeError(f"only placed {len(placed)}/{k} fliers")
    return placed, center


# ---------------------------------------------------------------------------
# observables
# ---------------------------------------------------------------------------

def measure(view, center):
    """Clustering + spread observables of the current defect configuration."""
    from discrete_differential_geometry.recording import _components
    pairs, degs = view.illegal_edges()
    pairs = np.asarray(pairs).reshape(-1, 2)
    mv = Counter()
    for a, b in pairs:
        mv[int(a)] += 1
        mv[int(b)] += 1
    comps = _components(pairs)
    row = {
        "n_ill": int(len(degs)),
        "n3": int(np.sum(np.asarray(degs) == 3)),
        "m2": int(sum(c * c for c in mv.values())),
        "ncomp": len(comps),
        "top_comps": comps[:3],
    }
    if mv:
        # BFS on the CURRENT edge graph from the (fixed) blob center
        E = np.asarray(view.simplices(1))
        V = int(E.max()) + 1
        adj = [[] for _ in range(V)]
        for a, b in E:
            adj[int(a)].append(int(b))
            adj[int(b)].append(int(a))
        dist = np.full(V, -1, np.int32)
        dist[center] = 0
        dq = deque([center])
        while dq:
            u = dq.popleft()
            for w in adj[u]:
                if dist[w] < 0:
                    dist[w] = dist[u] + 1
                    dq.append(w)
        dv = [int(dist[v]) for v in mv if dist[v] >= 0]
        row["spread_mean"] = round(float(np.mean(dv)), 3) if dv else None
        row["spread_max"] = int(max(dv)) if dv else None
    else:
        row["spread_mean"] = row["spread_max"] = None
    return row


# ---------------------------------------------------------------------------
# the lifted episode (arm B)
# ---------------------------------------------------------------------------

def run_episode(s, rng, nstep, kscan, beta, p_hand, hand_rule, audit,
                refresh_every=10):
    """One flight episode from a uniformly chosen deg-3 chord. Returns a
    summary dict; 'd_calls' counts nonlocal_slide_at invocations (work).

    refresh_every caps the lift's persistence: the direction is resampled
    every that many steps (0 = never within an episode -- the episode
    boundary is then the only refresh, since each episode starts from a fresh
    chord and a uniform frame). It was hardcoded to 10 through the
    2026-08-05 A/B, which put a HARD 10-step ceiling on persistence and so on
    any ballistic (t^2) signal -- with contacts reflecting ~42% of steps the
    realised persistence was ~2.4. Keep 10 to reproduce those runs; raise it
    (or 0) for the MSD diagnostic."""
    from ecmc_flight import Flight
    pairs, degs = s.manifold.illegal_edges()
    pairs = np.asarray(pairs).reshape(-1, 2)
    d3 = [tuple(sorted(map(int, p)))
          for p, d in zip(pairs, np.asarray(degs)) if d == 3]
    if not d3:
        return {"skipped": "no deg-3 chord", "d_calls": 0}
    C = d3[int(rng.integers(len(d3)))]

    calls = [0]
    orig = s.nonlocal_slide_at

    def counted(*a, **kw):
        calls[0] += 1
        return orig(*a, **kw)

    s.nonlocal_slide_at = counted
    try:
        # link triangle of the chord, for the initial frame
        lk = set()
        for t in np.asarray(s.manifold.facets()):
            t = [int(x) for x in t]
            if C[0] in t and C[1] in t:
                lk |= {v for v in t if v not in C}
        fl = Flight(s, C, (C[0],) + tuple(sorted(lk)), kscan=kscan,
                    audit=audit, beta=beta, p_hand=p_hand,
                    hand_rule=hand_rule, expect_free=False)
        fl.refresh_frame(rng)          # initial direction draw
        for i in range(nstep):
            if refresh_every and i and i % refresh_every == 0:
                fl.refresh_frame(rng)
            fl.step(rng)
        return {"chord0": list(C), "chord1": list(fl.C),
                "dS": round(fl.dS_sum, 6), "d_calls": calls[0],
                "events": dict(fl.events)}
    except Exception as ex:                        # noqa: BLE001
        return {"chord0": list(C), "failed": repr(ex)[:200],
                "d_calls": calls[0]}
    finally:
        del s.nonlocal_slide_at                   # restore the class method


# ---------------------------------------------------------------------------
# one chain
# ---------------------------------------------------------------------------

def run_chain(arm, sidx, seed, cfg, out_dir):
    import discrete_differential_geometry as ddg
    from discrete_differential_geometry import (Manifold, ManifoldSampler,
                                                SamplerParams)
    from discrete_differential_geometry.recording import Recorder

    rng = np.random.default_rng(seed)
    ddg.set_random_seed(int(rng.integers(2 ** 31)))

    m = Manifold.load(cfg["ref"], 3)
    blob, center = build_blob(m, cfg["k_fliers"], rng, radius=cfg["radius"])
    fv = [int(x) for x in m.f_vector]
    estar = 6.0 * fv[3] / fv[1]                   # pins exactly zero now
    params = SamplerParams(
        num_facets_target=fv[3], num_facets_coef=cfg["nfc"],
        hinge_degree_target=estar, num_hinges_coef=cfg["k1"],
        hinge_degree_variance_coef=0.0, codim3_degree_variance_coef=0.0,
        hinge_degree_target_coef=0.0)
    s = ManifoldSampler(m, params)
    s.set_n6_potential(0.0, cfg["cimp"])
    if arm == "N":
        s.set_nonlocal_slide_prob(cfg["nslide_prob"], cfg["nslide_max_step"])

    tag = f"arm{arm}_s{sidx:02d}"
    out = os.path.join(out_dir, tag)
    rec = Recorder(s, out, chunk=cfg["chunk"], census_every=1,
                   snap_mid=False,
                   extra_meta={"arm": arm, "seed": int(seed),
                               "blob": [list(map(int, f)) for f, _ in blob],
                               "center": center, "estar": estar,
                               "ab_cfg": {k: v for k, v in cfg.items()
                                          if k != "ref"}})
    side = open(out + ".ab.jsonl", "w")

    def swrite(row):
        side.write(json.dumps(row) + "\n")
        side.flush()

    work_plain, work_ep = 0, 0
    m0 = measure(s.manifold, center)
    swrite({"sw": 0, "work_plain": 0, "work_ep": 0, **m0})

    nchunks = cfg["sweeps"] // cfg["chunk"]
    for ci in range(nchunks):
        row = rec.step()
        work_plain += int(row["d_tried"])
        ep_ev, ep_summ = Counter(), []
        if arm == "B":
            for _ in range(cfg["eps_per_chunk"]):
                r = run_episode(s, rng, cfg["ep_steps"], cfg["kscan"],
                                cfg["beta"], cfg["p_hand"], cfg["hand_rule"],
                                cfg["audit"], cfg.get("refresh_every", 10))
                work_ep += 2 * r["d_calls"]
                ep_ev.update(r.get("events", {}))
                if "failed" in r:
                    ep_ev["EPISODE_FAILED"] += 1
                ep_summ.append({k: r[k] for k in ("chord0", "chord1", "dS",
                                                  "failed", "skipped")
                                if k in r})
        mr = measure(s.manifold, center)
        out_row = {"sw": rec.sw, "work_plain": work_plain,
                   "work_ep": work_ep, "obj": row["obj"], **mr}
        if ep_ev:
            out_row["ep_events"] = dict(ep_ev)
        if ep_summ:
            out_row["episodes"] = ep_summ
        swrite(out_row)
    rec.finish()
    side.close()
    return tag, work_plain, work_ep


def _worker(job):
    arm, sidx, seed, cfg, out_dir = job
    t0 = time.time()
    tag, wp, we = run_chain(arm, sidx, seed, cfg, out_dir)
    return tag, wp, we, time.time() - t0


# ---------------------------------------------------------------------------

def main():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--ref", default=REF)
    p.add_argument("--arms", default="AB", help="subset of 'ABN'")
    p.add_argument("--seeds", type=int, default=4, help="chains per arm")
    p.add_argument("--seed0", type=int, default=1000)
    p.add_argument("--k-fliers", type=int, default=12)
    p.add_argument("--radius", type=int, default=4, help="blob ball radius")
    p.add_argument("--cimp", type=float, default=0.5)
    p.add_argument("--nfc", type=float, default=0.1,
                   help="volume-pin coefficient")
    p.add_argument("--k1", type=float, default=1.0,
                   help="f1-pin (num_hinges) coefficient")
    p.add_argument("--sweeps", type=int, default=4000)
    p.add_argument("--chunk", type=int, default=20)
    p.add_argument("--eps-per-chunk", type=int, default=2)
    p.add_argument("--ep-steps", type=int, default=40)
    p.add_argument("--kscan", type=int, default=60)
    p.add_argument("--refresh-every", type=int, default=10,
                   help="resample the lift direction every N kernel steps "
                        "(0 = never within an episode). 10 reproduces the "
                        "2026-08-05 A/B, where it was hardcoded and capped "
                        "persistence at 10 steps.")
    p.add_argument("--beta", type=float, default=1.0)
    p.add_argument("--p-hand", type=float, default=0.0)
    p.add_argument("--hand-rule", choices=("alg", "chir"), default="alg")
    p.add_argument("--no-audit", action="store_true")
    p.add_argument("--nslide-prob", type=float, default=0.1)
    p.add_argument("--nslide-max-step", type=int, default=8)
    p.add_argument("--out-dir", required=True)
    p.add_argument("--max-workers", type=int, default=8)
    args = p.parse_args()

    cfg = dict(ref=args.ref, k_fliers=args.k_fliers, radius=args.radius,
               cimp=args.cimp, nfc=args.nfc, k1=args.k1, sweeps=args.sweeps,
               chunk=args.chunk, eps_per_chunk=args.eps_per_chunk,
               ep_steps=args.ep_steps, kscan=args.kscan, beta=args.beta,
               refresh_every=args.refresh_every,
               p_hand=args.p_hand, hand_rule=args.hand_rule,
               audit=not args.no_audit, nslide_prob=args.nslide_prob,
               nslide_max_step=args.nslide_max_step)
    os.makedirs(args.out_dir, exist_ok=True)
    with open(os.path.join(args.out_dir, "manifest.json"), "w") as f:
        json.dump({"cfg": cfg, "arms": args.arms, "seeds": args.seeds,
                   "seed0": args.seed0, "argv": sys.argv}, f, indent=1)

    # same seed for the SAME index across arms => identical initial blobs
    jobs = [(arm, i, args.seed0 + i, cfg, args.out_dir)
            for arm in args.arms for i in range(args.seeds)]
    print(f"ecmc_ab: {len(jobs)} chains ({args.arms} x {args.seeds}), "
          f"k={args.k_fliers} cimp={args.cimp} pins {args.nfc}/{args.k1}, "
          f"{args.sweeps} sw chunk {args.chunk}, "
          f"B: {args.eps_per_chunk} x {args.ep_steps} kernel steps/chunk",
          flush=True)
    with ProcessPoolExecutor(max_workers=args.max_workers) as ex:
        futs = [ex.submit(_worker, j) for j in jobs]
        for fu in as_completed(futs):
            tag, wp, we, dt = fu.result()
            print(f"  {tag}: plain {wp:,}  episode {we:,}  ({dt:.0f}s)",
                  flush=True)
    print(f"done -> {args.out_dir}", flush=True)


if __name__ == "__main__":
    main()
