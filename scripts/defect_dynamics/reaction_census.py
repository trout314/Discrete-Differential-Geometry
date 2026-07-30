#!/usr/bin/env python3
"""Reaction census: merger/split chemistry of defect complexes under the
FULL thermal dynamics.

Everything we know about how two defect complexes react was measured in
one of two artificial settings: the STATIC colliders (knot_collider,
crossing_collider -- a frozen background and deterministic directed
slides, giving V(s) = 0 exactly beyond contact) and the SLIDE-ONLY FP
driver (fp_encounter -- one mobile knot, frozen everything else, giving
P(recombine | freed) = 0.70). Neither has a thermal vacuum. Under the
full dynamics the flicker background is present, complexes change
species as they move, and the actual reactions are MERGES and SPLITS of
worldlines -- which defect_state.Worldlines has always counted (two
integers) and never characterised.

This driver characterises them. Complexes are tracked at accepted-move
resolution; every identity change is streamed as an event with the
participants' sizes, ages and induced shapes. The observables:

  * k_merge, k_split as rates per sweep, and their scaling with the
    complex number density n (merge should be bimolecular ~ n^2, split
    unimolecular ~ n -- the test that the "reaction" picture is even
    the right language);
  * the association constant K = k_merge / k_split at steady state,
    which is a FREE ENERGY of association. The static two-body potential
    is exactly zero beyond contact, so any K != ideal is entropic /
    kinetic in origin -- a sharp, falsifiable number;
  * compound lifetimes: run5h observed "fused multimers are immortal"
    from snapshots; here it becomes a survival curve;
  * the species table: which induced shapes merge into which.

SLIDE ON/OFF is the clean control. The slide channel is certified to
satisfy detailed balance (V3b), so switching it on changes the KINETICS
and not the stationary ensemble: if reaction rates move, the chemistry
is transport-limited; if they do not, it is reaction-limited.

FRAME: none -- no positions are used. Events, sizes, shapes and graph
identity are all frame-free (CONVENTIONS.md sec 6). Separations and
intrinsic collision classes are deliberately left to an offline pass.

Streams (unbounded event rate, so nothing is accumulated in memory):
  <out>.events.jsonl   merge / split records
  <out>.tracks.jsonl   finalised worldline records (lifetime, max size)
  <out>.json           params, counters, census time series
"""
import argparse
import json
import os
import sys
import time
from collections import Counter

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
for _p in ("../../python", "../../scripts", "../../tools", "."):
    sys.path.insert(0, os.path.join(_HERE, _p))
import discrete_differential_geometry as ddg
from fk_skeleton import edges_from_facets
import defect_state as ds


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--cell", required=True)
    ap.add_argument("--start", default=None)
    ap.add_argument("--lam", type=float, required=True)
    ap.add_argument("--etarget", type=float, default=None)
    ap.add_argument("--zleg", type=float, default=0.6)
    ap.add_argument("--cimp", type=float, default=1.0)
    ap.add_argument("--slide-prob", type=float, default=0.0)
    ap.add_argument("--worm-prob", type=float, default=0.0,
                    help="deg-4 worm channel probability per step (D-side "
                         "catalysed transport move; 0 = off)")
    ap.add_argument("--sweeps", type=int, default=10 ** 9)
    ap.add_argument("--max-seconds", type=float, default=1e9,
                    help="wall-clock deadline; runs stop cleanly at it")
    ap.add_argument("--chunk", type=int, default=25)
    ap.add_argument("--logmb", type=float, default=256.0)
    ap.add_argument("--audit-every", type=int, default=40,
                    help="chunks between incremental-state audits")
    ap.add_argument("--keep-life", type=float, default=1.0,
                    help="stream a track record if it lived >= this many "
                         "sweeps (shorter ones are pure flicker; they "
                         "still enter the lifetime histogram)")
    ap.add_argument("--seed", type=int, default=4242)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    print("[frame] none -- graph identity only; no positions, no gauge")

    t_start = time.time()
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
    s.set_n6_potential(args.zleg * args.lam, args.cimp * args.lam,
                       tilt=[0.0] * 5)
    if args.slide_prob:
        s.set_slide_prob(args.slide_prob)
    if args.worm_prob:
        s.set_worm_prob(args.worm_prob)
    v = s.manifold
    n3 = v.num_facets

    # maximalist flight recorder: unconditional per-chunk instrumentation
    # (f-vector, <e>, acceptance, channel stats, GC heap, defect census)
    # alongside this script's own event pipeline; also saves mid + final
    # .mfd states. See discrete_differential_geometry.recording.
    from discrete_differential_geometry.recording import Recorder
    rec = Recorder(s, args.out, chunk=args.chunk,
                   extra_meta={"cell": args.cell, "start": args.start,
                               "lam": args.lam, "zleg": args.zleg,
                               "cimp": args.cimp, "seed": args.seed})

    st = ds.DefectState(v)
    wl = ds.Worldlines(st)
    s.enable_event_log(args.logmb)
    s.drain_event_log()

    ev_fh = open(args.out + ".events.jsonl", "w")
    tr_fh = open(args.out + ".tracks.jsonl", "w")

    # live per-track bookkeeping: tid -> [born_clock, max_size, n_merges_in,
    #                                     from_split, modal shapes]
    live = {}
    life_hist = Counter()          # lifetime in whole sweeps -> count
    counters = Counter()
    ts = []

    def shape_of(verts):
        try:
            return (tuple(int(x) for x in st.induced_shape(verts)["f"]),
                    int(st.total_coordination(verts)))
        except Exception:
            return None

    def finalize(tid, clock, how="death"):
        rec = live.pop(tid, None)
        if rec is None:
            return
        life_sw = (clock - rec[0]) / n3
        life_hist[int(min(life_sw, 999))] += 1
        counters["tracks_total"] += 1
        # always keep a track that ENDED BY MERGING, however short-lived:
        # its shape is one half of the reaction A + B -> C, and most
        # merge participants are short-lived flicker knots.
        if (how == "merge" or life_sw >= args.keep_life
                or rec[1] >= 8 or rec[2]):
            tr_fh.write(json.dumps({
                "tid": tid, "born_sw": rec[0] / n3, "life_sw": life_sw,
                "max_size": rec[1], "n_merge_in": rec[2],
                "ended": how, "from_split": rec[3],
                "shape": list(rec[4].most_common(1)[0][0][0])
                if rec[4] else None,
                "coord": rec[4].most_common(1)[0][0][1] if rec[4] else None,
            }) + "\n")

    def record(clock):
        out = wl.step(clock)
        for cx, tid in out:
            rec = live.get(tid)
            if rec is None:
                rec = live[tid] = [wl.born.get(tid, clock), 0, 0, False,
                                   Counter()]
            n = len(cx.verts)
            if n > rec[1]:
                rec[1] = n
            sh = shape_of(cx.verts)
            if sh is not None:
                rec[4][sh] += 1
        for e in wl.drain_events():
            k = e["k"]
            counters[k] += 1
            if k == "merge":
                into = live.get(e["into"])
                e["age_into_sw"] = ((e["clock"] - into[0]) / n3
                                    if into else None)
                e["from_ages_sw"] = [
                    (e["clock"] - live[f[0]][0]) / n3
                    if f[0] in live else None for f in e["from"]]
                e["from_max_sizes"] = [
                    live[f[0]][1] if f[0] in live else None
                    for f in e["from"]]
                e["clock_sw"] = e["clock"] / n3
                if into:
                    into[2] += 1
                ev_fh.write(json.dumps(e) + "\n")
                for f in e["from"]:
                    finalize(f[0], e["clock"], how="merge")
            elif k == "split":
                par = live.get(e["parent"])
                e["parent_age_sw"] = ((e["clock"] - par[0]) / n3
                                      if par else None)
                e["clock_sw"] = e["clock"] / n3
                ev_fh.write(json.dumps(e) + "\n")
                ch = live.get(e["child"])
                if ch is not None:
                    ch[3] = True
            elif k == "death":
                finalize(e["tid"], e["clock"])

    record(0)
    done = 0
    nchunk = 0
    audit_fail = []
    while done < args.sweeps and time.time() - t_start < args.max_seconds:
        rec.step(args.chunk)
        done += args.chunk
        nchunk += 1
        ev = s.drain_event_log()
        if s.event_log_overflowed():
            counters["event_log_overflow"] += 1
            print("  WARN event-log overflow", flush=True)
        for e in ev:
            st.apply(e)
            record(int(e["clock"]))
        c = ds.census(st)
        ts.append({"sw": done, "n_illegal": c["n_illegal"],
                   "n_components": c["n_components"],
                   "sizes": c["sizes"]})
        if args.audit_every and nchunk % args.audit_every == 0:
            p = st.audit(v)
            if p:
                audit_fail.append((done, p[:2]))
                print(f"  AUDIT FAILED sw{done}: {p[:2]}", flush=True)
            el = time.time() - t_start
            print(f"  sw{done} ({el:.0f}s): ncomp={c['n_components']} "
                  f"nill={c['n_illegal']} merges={counters['merge']} "
                  f"splits={counters['split']} births={counters['birth']}",
                  flush=True)
            ev_fh.flush()
            tr_fh.flush()

    for tid in list(live):
        finalize(tid, done * n3)
    ev_fh.close()
    tr_fh.close()
    rec.finish()                    # final .mfd snapshot + closing record

    nc = np.array([r["n_components"] for r in ts], float)
    summary = {
        "cell": args.cell, "start": args.start, "lam": args.lam,
        "zleg": args.zleg, "cimp": args.cimp,
        "slide_prob": args.slide_prob, "worm_prob": args.worm_prob,
        "worm_stats": list(s.worm_stats()), "etarget": et, "n3": n3,
        "seed": args.seed, "sweeps": done,
        "wall_s": time.time() - t_start,
        "counters": dict(counters),
        "life_hist": {str(k): v for k, v in sorted(life_hist.items())},
        "mean_ncomp": float(nc.mean()) if len(nc) else None,
        "mean_ncomp_sq": float((nc ** 2).mean()) if len(nc) else None,
        "k_merge_per_sweep": counters["merge"] / max(done, 1),
        "k_split_per_sweep": counters["split"] / max(done, 1),
        "audit_failures": audit_fail[:5],
        "ts": ts,
    }
    with open(args.out + ".json", "w") as fh:
        json.dump(summary, fh)
    print(f"\n{done} sweeps, {counters['merge']} merges, "
          f"{counters['split']} splits, {counters['birth']} births, "
          f"{counters['tracks_total']} tracks")
    print(f"wrote {args.out}.json / .events.jsonl / .tracks.jsonl")


if __name__ == "__main__":
    main()
