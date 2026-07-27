#!/usr/bin/env python3
"""Pinpoint the DefectState incremental-divergence bug (2026-07-27).

Every long reaction_census run eventually fails its self-audit with
`edge degrees differ` / `illegal-edge set differs` but NEVER `tet set
differs`. That asymmetry is the whole clue: the maintained tet set
(add/discard -- idempotent) stays exactly right while edge degrees
(edeg[p] += / -=  -- NOT idempotent) corrupt and never recover. So the
event stream must carry an anomaly that is a no-op on the tet SET but
not on the per-edge COUNTS -- a duplicated event is the textbook case
(replaying a 2->3 twice: the second discard/add no-op the set but shift
every incident edge by an extra +-1).

This driver reproduces a known-failing run and, per accepted event,
checks against a strict ground-truth tet set BEFORE applying:
  * every `rem` tet must currently be present   (else: phantom removal)
  * every `add` tet must currently be absent     (else: phantom addition)
The first violation is printed with the event, whether it duplicates a
recent event, and the running edeg-vs-groundtruth delta -- which is the
exact culprit, found in O(events) with no from-scratch rebuilds.

No fix is applied; this only diagnoses.
"""
import argparse
import os
import sys
from collections import deque, Counter
from itertools import combinations

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
for _p in ("../../python", "../../scripts", "../../tools", "."):
    sys.path.insert(0, os.path.join(_HERE, _p))
import discrete_differential_geometry as ddg
from fk_skeleton import edges_from_facets
import event_replay as er
import defect_state as ds


def tet_key(ev):
    return (int(ev["type"]), tuple(int(x) for x in ev["labels"]))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--cell", default="data/tcp_reference/T3_R_m4_N57984.mfd")
    ap.add_argument("--start", default="data/mgas/lam35_snap17000.mfd")
    ap.add_argument("--lam", type=float, default=0.35)
    ap.add_argument("--slide-prob", type=float, default=0.10)
    ap.add_argument("--etarget", type=float, default=5.105025)
    ap.add_argument("--zleg", type=float, default=0.6)
    ap.add_argument("--cimp", type=float, default=1.0)
    ap.add_argument("--seed", type=int, default=1005)
    ap.add_argument("--sweeps", type=int, default=9000)
    ap.add_argument("--chunk", type=int, default=25)
    ap.add_argument("--logmb", type=float, default=64.0)
    args = ap.parse_args()

    ddg.set_random_seed(args.seed)
    ref = ddg.Manifold.load(args.cell, 3)
    et = args.etarget
    m = ddg.Manifold.load(args.start, 3)
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
    v = s.manifold

    # ground-truth tet set, updated per event with STRICT presence checks
    G = {tuple(sorted(int(x) for x in t)) for t in np.asarray(v.facets())}
    Gdeg = Counter()
    for t in G:
        for e in combinations(t, 2):
            Gdeg[e] += 1

    st = ds.DefectState(v)          # the suspect incremental tracker
    s.enable_event_log(args.logmb)
    s.drain_event_log()

    recent = deque(maxlen=64)       # (idx, key) for duplicate detection
    seen_all = Counter()            # key -> times observed (whole run)
    idx = 0
    done = 0
    n_phantom_rem = n_phantom_add = 0

    def report(kind, ev, extra=""):
        k = tet_key(ev)
        dup = [i for i, kk in recent if kk == k]
        print(f"\n*** {kind} at event #{idx} (sweep ~{done}) ***")
        print(f"    type={k[0]} labels={k[1]}")
        print(f"    seen this key {seen_all[k]}x before; "
              f"recent-duplicate at events {dup if dup else 'none'}")
        rem, add = er.event_changes(ev)
        print(f"    rem={rem}")
        print(f"    add={add}")
        if extra:
            print(f"    {extra}")

    while done < args.sweeps:
        s.run(sweeps=args.chunk)
        done += args.chunk
        ev = s.drain_event_log()
        if s.event_log_overflowed():
            print(f"  [!] event-log overflow at sweep {done}")
        for e in ev:
            rem, add = er.event_changes(e)
            k = tet_key(e)
            # STRICT checks against ground truth, before applying
            bad = None
            for r in rem:
                if r not in G:
                    bad = ("PHANTOM REMOVE", r)
                    break
            if bad is None:
                for a in add:
                    if a in G:
                        bad = ("PHANTOM ADD", a)
                        break
            if bad is not None:
                n_phantom_rem += bad[0] == "PHANTOM REMOVE"
                n_phantom_add += bad[0] == "PHANTOM ADD"
                if n_phantom_rem + n_phantom_add <= 5:
                    report(bad[0], e, f"offending tet {bad[1]}")
            # advance ground truth (idempotent set ops, mirrors live)
            for r in rem:
                G.discard(r)
                for p in combinations(r, 2):
                    Gdeg[p] -= 1
                    if Gdeg[p] == 0:
                        del Gdeg[p]
            for a in add:
                G.add(a)
                for p in combinations(a, 2):
                    Gdeg[p] += 1
            st.apply(e)
            recent.append((idx, k))
            seen_all[k] += 1
            idx += 1

        # cross-check incremental edeg vs ground truth at chunk boundary
        mine = {e: d for e, d in st.edeg.items() if d > 0}
        gt = dict(Gdeg)
        if mine != gt:
            diff = set(mine.items()) ^ set(gt.items())
            print(f"\n=== edeg DIVERGED by sweep {done} "
                  f"({len(diff)} entries); "
                  f"phantom_rem={n_phantom_rem} phantom_add={n_phantom_add} ===")
            # is the ground truth itself still exact vs the live manifold?
            live = {tuple(sorted(int(x) for x in t))
                    for t in np.asarray(v.facets())}
            print(f"    ground-truth tets vs live: "
                  f"{'MATCH' if live == G else f'DIFFER by {len(live ^ G)}'}")
            print(f"    DefectState tets vs live:  "
                  f"{'MATCH' if live == st.tets else f'DIFFER by {len(live ^ st.tets)}'}")
            # show a few corrupted edges
            for e, d in list(diff)[:6]:
                print(f"    edge {e}: mine={mine.get(e)} truth={gt.get(e)}")
            print(f"    total accepted events so far: {idx}")
            print(f"    distinct keys seen >1x: "
                  f"{sum(1 for c in seen_all.values() if c > 1)}")
            return
        if done % 500 == 0:
            print(f"  sweep {done}: {idx} events, clean "
                  f"(phantom_rem={n_phantom_rem} phantom_add={n_phantom_add})",
                  flush=True)

    print(f"\nno divergence in {args.sweeps} sweeps ({idx} events)")


if __name__ == "__main__":
    main()
