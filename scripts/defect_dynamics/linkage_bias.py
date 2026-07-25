#!/usr/bin/env python3
"""Does frame-based overlap linkage invent deaths that are really motion?

pass1_kinematics identifies a defect across 150-sweep frames by member
overlap: >= 2 shared members, or >= 1 shared with Jaccard >= 0.3. Over 150
sweeps a defect's membership can turn over entirely, so a defect that MOVED
looks like a death plus a nearby birth. pass1 reports exactly that pattern --
1599 death->rebirth events within 0.7 cell, a recurrence fraction of
1599/1948 = 0.82 -- and interprets it as "blinking". The competing explanation
is that the linkage simply lost the object.

Those two explanations make opposite predictions, and per-move identity can
tell them apart. This script runs one chain and tracks it BOTH ways on exactly
the same data:

  * exact   -- defect_state.Worldlines, relabelled after every accepted move
               (a move touches <= 5 vertices, so untouched complexes keep
               their identity by construction);
  * overlap -- pass1's rule applied to the same components sampled at the same
               150-sweep frames.

Then it cross-tabulates: of the deaths the overlap rule reports, how many were
the SAME OBJECT still alive according to exact identity? A high fraction means
"blinking" and "94% frozen" are linkage artifacts, and the frame-based tracker
should not be used for kinematics.
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

CELL = 1.0e6


def overlap_link(prev, cur):
    """pass1's rule: greedy by shared-member count, requiring >= 2 shared or
    Jaccard >= 0.3. `prev`/`cur` are lists of vertex sets. Returns a dict
    cur_index -> prev_index for the matched ones."""
    cand = []
    for pi, a in enumerate(prev):
        for ci, b in enumerate(cur):
            sh = len(a & b)
            if sh >= 2 or (sh and sh / len(a | b) >= 0.3):
                cand.append((sh, pi, ci))
    used_p, used_c, out = set(), set(), {}
    for sh, pi, ci in sorted(cand, reverse=True):
        if pi in used_p or ci in used_c:
            continue
        used_p.add(pi)
        used_c.add(ci)
        out[ci] = pi
    return out


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
    ap.add_argument("--frame", type=int, default=150,
                    help="frame length in sweeps (pass1 used 150)")
    ap.add_argument("--frames", type=int, default=20)
    ap.add_argument("--logmb", type=float, default=256.0)
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
    else:
        edges = np.asarray(v.simplices(1))
        s.enable_cocycle(edges, coc.build_from_positions(
            edges, reference_frac_positions("r", args.mcell), args.mcell))
    pbox = CELL * args.mcell
    s.enable_cocycle_positions([int(pbox)] * 3)

    st = ds.DefectState(v)
    wl = ds.Worldlines(st)
    s.enable_event_log(args.logmb)
    s.drain_event_log()

    def centroid(verts):
        q = s.vertex_positions(verts).astype(float)
        d = (q - q[0] + pbox / 2) % pbox - pbox / 2
        return (q[0] + d.mean(0)) / CELL

    wl.step(0)
    frames = []          # per frame: list of (vertex set, exact tid, centroid)
    def snap():
        cur = []
        for c in st.components():
            tid = Counter(wl.label[x] for x in c.verts
                          if x in wl.label).most_common(1)
            cur.append((set(c.verts), tid[0][0] if tid else -1,
                        centroid(c.verts), c.sig))
        frames.append(cur)
    snap()

    n_ev = 0
    for f in range(args.frames):
        s.run(sweeps=args.frame)
        ev = s.drain_event_log()
        if s.event_log_overflowed():
            print("  WARN event-log overflow", flush=True)
        # advance the exact tracker move by move
        clk = 0
        for e in ev:
            st.apply(e)
            n_ev += 1
            clk = int(e["clock"])
            wl.step(clk)
        snap()
        print(f"  frame {f+1}/{args.frames}: {n_ev} events, "
              f"{len(frames[-1])} complexes", flush=True)

    # --- compare the two linkages frame to frame
    tot_deaths = tot_alive = tot_matched = 0
    dists = []
    for i in range(len(frames) - 1):
        prev, cur = frames[i], frames[i + 1]
        match = overlap_link([p[0] for p in prev], [c[0] for c in cur])
        matched_prev = set(match.values())
        cur_exact = {c[1] for c in cur}
        for pi, (pverts, ptid, pcen, psig) in enumerate(prev):
            if pi in matched_prev:
                tot_matched += 1
                continue
            tot_deaths += 1                      # overlap rule says: died
            if ptid in cur_exact:                # exact rule says: still alive
                tot_alive += 1
                for cverts, ctid, ccen, csig in cur:
                    if ctid == ptid:
                        d = ccen - pcen
                        d -= np.round(d / args.mcell) * args.mcell
                        dists.append(float(np.linalg.norm(d)))
                        break

    frac = tot_alive / tot_deaths if tot_deaths else 0.0
    summary = {
        "frame_sweeps": args.frame, "frames": len(frames), "events": n_ev,
        "slide_prob": args.slide_prob,
        "overlap_matched": tot_matched,
        "overlap_deaths": tot_deaths,
        "of_which_still_alive_exactly": tot_alive,
        "false_death_fraction": frac,
        "median_displacement_cells":
            float(np.median(dists)) if dists else None,
        "exact_merges": wl.merges, "exact_splits": wl.splits,
    }
    with open(args.out + ".json", "w") as f:
        json.dump(summary, f, indent=1)
    print(f"\nframes {len(frames)} of {args.frame} sweeps, {n_ev} events")
    print(f"  overlap linkage: {tot_matched} continuations, "
          f"{tot_deaths} deaths")
    print(f"  of those deaths, still alive under EXACT identity: {tot_alive}"
          f"  ({100*frac:.1f}%)")
    if dists:
        print(f"  their displacement over one frame: median "
              f"{np.median(dists):.3f} cells, p90 "
              f"{np.percentile(dists, 90):.3f}")
    print(f"  exact tracker: {wl.merges} merges, {wl.splits} splits")
    print(f"wrote {args.out}.json")


if __name__ == "__main__":
    main()
