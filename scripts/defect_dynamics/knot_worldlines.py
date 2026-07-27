#!/usr/bin/env python3
"""Defect worldlines and MSD at full accepted-move resolution.

Frame-based trackers (pass1_kinematics) link complexes between 150-sweep
snapshots by member overlap (>=2 shared members, or Jaccard >= 0.3). Over a
150-sweep gap a defect's membership can turn over completely, so a move reads
as a death plus a birth -- and the bias is ONE-SIDED, suppressing displacement
exactly where transport is fastest. With the slide channel on that is fatal:
it manufactures a "slides don't move defects" null.

Here identity is established per ACCEPTED MOVE instead, which removes almost
all of the inference. One Pachner move touches at most 5 vertices (center u
coCenter), so every complex not meeting that set is untouched -- same
vertices, same structure, identity trivially preserved. Only complexes meeting
those <=5 vertices can change, and they are relabelled by majority
inheritance, which at one-move granularity is unambiguous in a way 150-sweep
overlap linkage is not. Merges and splits are recorded rather than silently
becoming deaths and births.

Slides are identified exactly on top of that. logEvent stamps the ledger clock
and clock++ fires once per mcmcStep, so a slide's four moves all carry ONE
clock value while every ordinary move is one event per tick: a slide is a
same-clock run of four events typed (3->2, 2->3, 2->3, 3->2). So the tracker
can report how much observed transport came through the slide channel rather
than thermally -- which the previous version could not, because it scored
thermal chord relocation as zero displacement by construction.

Defect bookkeeping comes from defect_state, so this shares one audited census
with the rest of the tree and uses the BROADENED defect definition (illegal
edge incidence OR non-FK coordination). Positions come from the D-side
incremental cocycle lift: O(1) per query, gauge fixed for the run, so
displacements need no glitch filtering.

Time is the ledger clock (attempted moves) divided by N3 -- finer and more
honest than chunk boundaries.

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
T_32, T_23 = 2, 1
SLIDE_PATTERN = [T_32, T_23, T_23, T_32]


def group_by_clock(ev):
    """Consecutive events sharing a clock value = one mcmcStep."""
    out, cur, clk = [], [], None
    for e in ev:
        c = int(e["clock"])
        if clk is not None and c != clk:
            out.append(cur)
            cur = []
        cur.append(e)
        clk = c
    if cur:
        out.append(cur)
    return out


def is_slide(grp):
    return (len(grp) == 4
            and [int(e["type"]) for e in grp] == SLIDE_PATTERN)


class Track:
    __slots__ = ("tid", "born", "samples", "sig_seen", "nodes_seen",
                 "slide_steps", "thermal_steps", "merged_from", "died")

    def __init__(self, tid, clock):
        self.tid = tid
        self.born = clock
        self.samples = []          # (clock, centroid in cells)
        self.sig_seen = Counter()
        self.nodes_seen = Counter()
        self.slide_steps = 0
        self.thermal_steps = 0
        self.merged_from = 0
        self.died = None


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
    ap.add_argument("--sweeps", type=int, default=2000)
    ap.add_argument("--chunk", type=int, default=25,
                    help="sweeps between event drains -- an I/O cadence only; "
                         "identity is resolved per move regardless")
    ap.add_argument("--stride", type=int, default=1,
                    help="record a position sample every Nth accepted move")
    ap.add_argument("--species", default="3,4,4",
                    help="edge-degree signature to report MSD for, or 'all'")
    ap.add_argument("--logmb", type=float, default=256.0)
    ap.add_argument("--audit-every", type=int, default=4,
                    help="audit the incremental census every Nth chunk (0=off)")
    ap.add_argument("--seed", type=int, default=4242)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    print("[frame] registry lift -- exponents frame-robust, lengths/D gauge")

    want = None if args.species == "all" else tuple(
        int(x) for x in args.species.split(","))

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
    s.check_cocycle()
    pbox = CELL * args.mcell
    s.enable_cocycle_positions([int(pbox)] * 3)
    s.check_cocycle_positions()

    n3 = v.num_facets
    st = ds.DefectState(v)
    s.enable_event_log(args.logmb)
    s.drain_event_log()

    def centroid(verts):
        q = s.vertex_positions(verts).astype(float)
        d = (q - q[0] + pbox / 2) % pbox - pbox / 2
        return (q[0] + d.mean(0)) / CELL

    tracks = {}                 # tid -> Track
    label = {}                  # vertex -> tid
    counters = {"tid": 0, "ev": 0}
    n_slidegrp = 0
    overflow = False
    audit_fail = []

    def relabel(clock, was_slide):
        """Match the current components to the previous labelling."""
        comps = st.components()
        claimed, newlabel = set(), {}
        for c in comps:
            prev = Counter(label[x] for x in c.verts if x in label)
            t = None
            if prev:
                cand = prev.most_common(1)[0][0]
                if cand not in claimed:      # else this piece is a split-off
                    t = cand
                    if len(prev) > 1:
                        tracks[t].merged_from += len(prev) - 1
            if t is None:
                t = counters["tid"]
                counters["tid"] += 1
                tracks[t] = Track(t, clock)
            claimed.add(t)
            tr = tracks[t]
            tr.sig_seen[c.sig] += 1
            tr.nodes_seen[c.nodes] += 1
            if was_slide:
                tr.slide_steps += 1
            else:
                tr.thermal_steps += 1
            if counters["ev"] % args.stride == 0:
                tr.samples.append((clock, centroid(c.verts)))
            for x in c.verts:
                newlabel[x] = t
        for t in set(label.values()) - claimed:
            if tracks[t].died is None:
                tracks[t].died = clock
        label.clear()
        label.update(newlabel)

    relabel(0, False)
    done = 0
    while done < args.sweeps:
        s.run(sweeps=args.chunk)
        done += args.chunk
        ev = s.drain_event_log()
        if s.event_log_overflowed():
            overflow = True
            print(f"  WARN event-log overflow at sweep {done}", flush=True)
        for grp in group_by_clock(ev):
            sl = is_slide(grp)
            if sl:
                n_slidegrp += 1
            for e in grp:
                st.apply(e)
                counters["ev"] += 1
            relabel(int(grp[-1]["clock"]), sl)
        if args.audit_every and (done // args.chunk) % args.audit_every == 0:
            probs = st.audit(v)
            if probs:
                audit_fail.append((done, probs[:3]))
                print(f"  AUDIT FAILED at sweep {done}: {probs[:3]}", flush=True)
        c = ds.census(st)
        print(f"  sw{done}: events {counters['ev']} | defect {c['n_defect']} "
              f"(nonFK {c['n_nonfk_all_legal']}) | comps {c['n_components']} "
              f"| tracks {len(tracks)}", flush=True)

    # --- MSD over worldlines of the requested species
    sel = [t for t in tracks.values()
           if len(t.samples) >= 2
           and (want is None or (t.sig_seen and
                                 t.sig_seen.most_common(1)[0][0] == want))]
    bins = defaultdict(list)
    for t in sel:
        cl = np.array([x[0] for x in t.samples], float) / n3   # -> sweeps
        raw = np.array([x[1] for x in t.samples])
        un = np.zeros_like(raw)
        for i in range(1, len(raw)):
            step = raw[i] - raw[i - 1]
            step -= np.round(step / args.mcell) * args.mcell   # min image
            un[i] = un[i - 1] + step
        for i in range(len(un)):
            for j in range(i + 1, len(un)):
                bins[max(1, int(round(cl[j] - cl[i])))].append(
                    float(np.sum((un[j] - un[i]) ** 2)))
    rows = [{"dt_sweeps": k, "msd": float(np.mean(x)), "n": len(x)}
            for k, x in sorted(bins.items()) if len(x) >= 5]

    summary = {
        "slide_prob": args.slide_prob, "sweeps": done, "n3": n3,
        "events": counters["ev"], "slide_groups": n_slidegrp,
        "tracks_total": len(tracks), "tracks_used": len(sel),
        "species": args.species,
        "slide_steps": sum(t.slide_steps for t in sel),
        "thermal_steps": sum(t.thermal_steps for t in sel),
        "merges": sum(t.merged_from for t in tracks.values()),
        "overflow": overflow, "audit_failures": audit_fail,
        "msd": rows,
    }
    with open(args.out + ".msd.json", "w") as f:
        json.dump(summary, f, indent=1)
    life = Counter(len(t.samples) for t in tracks.values())
    print(f"\nevents {counters['ev']}; slide groups {n_slidegrp}; "
          f"tracks {len(tracks)} ({len(sel)} matching species {args.species})")
    print(f"worldline length histogram: "
          f"{dict(sorted(life.items())[:12])}")
    print(f"steps on selected tracks: slide {summary['slide_steps']}, "
          f"thermal {summary['thermal_steps']}; merges {summary['merges']}")
    for r in rows[:10]:
        print(f"   dt={r['dt_sweeps']:5d} sw   MSD={r['msd']:9.4f} cells^2  "
              f"(n={r['n']})")
    print(f"wrote {args.out}.msd.json")
    if audit_fail:
        print(f"AUDIT FAILURES: {audit_fail}")
        sys.exit(1)


if __name__ == "__main__":
    main()
