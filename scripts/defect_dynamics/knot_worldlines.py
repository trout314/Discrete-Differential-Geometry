#!/usr/bin/env python3
"""(3,4,4) knot worldlines and MSD from the accepted-move event stream.

Frame-based trackers (pass1_kinematics) link complexes between 150-sweep
snapshots by member overlap (>=2 shared members, or Jaccard >= 0.3). That
linkage is unsafe once the slide channel is on: a slide TRANSLATES a knot four
steps along its Boerdijk-Coxeter chain, turning over most of its membership, so
a hop can read as a death plus a birth. The bias is one-sided -- it suppresses
measured displacement exactly in the slides-on arm -- so it would manufacture a
"slides don't move defects" null.

This tracker avoids inference entirely by working at full accepted-move
resolution, following deg4_moves.py's protocol:

  * replay every accepted move (event_replay.event_changes) and maintain edge
    degrees incrementally -- no rescans;
  * identify a defect by its EDGE KEY, the sorted vertex pair. Vertex labels
    are stable under Pachner moves except at 1->4 / 4->1, so an edge keeps its
    identity for as long as it exists.

The one thing edge-key identity cannot see by itself is a slide, which
destroys chord (c0,c4) and creates (c4,c8) -- different keys, so it would read
as death + birth. Resolution rescues it: logEvent stamps the ledger clock and
clock++ fires once per mcmcStep, so a slide's four moves all carry ONE clock
value while every ordinary move is one event per tick. A slide is therefore
exactly identifiable as a same-clock run of four events typed
(3->2, 2->3, 2->3, 3->2), and within that window

    destroyed chord = labelsA of event 0      (3->2: A is the removed edge)
    created   chord = labelsB of event 2      (2->3: B is the created edge)

are the same defect before and after BY CONSTRUCTION. Worldlines are therefore
exact, not inferred.

Positions come from the cocycle tree lift, sampled at chunk boundaries (vertex
positions are a lift of the crystal, so per-event resolution buys nothing for
displacement). Steps use the validated protocol: min-image in cell units, with
single-step jumps above GLITCH cells dropped as spanning-tree gauge jumps.
"""
import argparse
import json
import os
import sys
from collections import Counter, defaultdict
from itertools import combinations

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
for p in ("../../python", "../../scripts", "../../tools", "."):
    sys.path.insert(0, os.path.join(_HERE, p))
import discrete_differential_geometry as ddg
from discrete_differential_geometry import cocycle as coc
import event_replay as er
from cocycle_check import reference_frac_positions
from fk_skeleton import edges_from_facets
from dopant_pairs import vertex_classes

CELL = 1.0e6                      # raw cocycle units per cell
GLITCH = 1.0                      # cells; larger single step = gauge jump
T_32, T_23 = 2, 1                 # event type codes
SLIDE_PATTERN = [T_32, T_23, T_23, T_32]


def ek(u, w):
    return (u, w) if u < w else (w, u)


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


def slide_chords(grp):
    """(destroyed_chord, created_chord) if this step is a slide, else None."""
    if len(grp) != 4:
        return None
    if [int(e["type"]) for e in grp] != SLIDE_PATTERN:
        return None
    a0, _ = er._support(grp[0])            # 3->2: A is the removed edge
    _, b2 = er._support(grp[2])            # 2->3: B is the created edge
    if len(a0) != 2 or len(b2) != 2:
        return None
    return ek(*a0), ek(*b2)


def knot_chords_now(v):
    """chord -> component, for every complex whose illegal-degree signature is
    exactly (3,4,4). This is the SPECIES definition; it is deliberately not
    "any degree-3 edge", because a 2->3 move creates its pole-pole edge AT
    degree 3, so ordinary MCMC churns transient degree-3 edges constantly."""
    fac = np.asarray(v.facets())
    eu, edg, V = edges_from_facets(fac)
    n6, imp, adj = vertex_classes(fac)
    lab = np.unique(fac)
    idx = {int(x): i for i, x in enumerate(lab)}
    illdeg = {frozenset(int(x) for x in e): int(d)
              for e, d in zip(eu, edg) if d not in (5, 6)}
    seen, out = set(), {}
    for i in range(len(lab)):
        if imp[i] <= 0 or i in seen:
            continue
        st, comp = [i], []
        seen.add(i)
        while st:
            u = st.pop()
            comp.append(u)
            for w in adj[u]:
                if imp[w] > 0 and w not in seen:
                    seen.add(w)
                    st.append(w)
        cv = {int(lab[x]) for x in comp}
        inside = {e: d for e, d in illdeg.items() if e <= cv}
        if sorted(inside.values()) != [3, 4, 4]:
            continue
        for e, d in inside.items():
            if d == 3:
                out[ek(*sorted(e))] = cv
    return out


class Track:
    __slots__ = ("tid", "born", "chord", "hops", "samples", "alive", "died",
                 "glitches", "sig_seen")

    def __init__(self, tid, sweep, chord):
        self.tid = tid
        self.born = sweep
        self.chord = chord
        self.hops = 0
        self.samples = []          # (sweep, unwrapped position in cells)
        self.alive = True
        self.died = None
        self.glitches = 0
        self.sig_seen = Counter()


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
    ap.add_argument("--sweeps", type=int, default=3000)
    ap.add_argument("--chunk", type=int, default=25,
                    help="sweeps between event drains / position samples")
    ap.add_argument("--logmb", type=float, default=256.0)
    ap.add_argument("--seed", type=int, default=4242)
    ap.add_argument("--verify", action="store_true",
                    help="also maintain the tet set and check the replay "
                         "against the live manifold each chunk (slow)")
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
    s.check_cocycle()
    # Incrementally-maintained vertex lift: O(1) per move and per query, with
    # the gauge FIXED for the run. The old path re-integrated the whole
    # cochain per sample (O(V+E) plus a full marshal) and rebuilt the spanning
    # tree each time, which can reassign a persisting vertex -- the gauge
    # glitch that had to be filtered. That cannot happen here.
    s.enable_cocycle_positions([int(CELL * args.mcell)] * 3)
    s.check_cocycle_positions()

    # --- incremental edge degrees (deg4_moves protocol)
    tets = {tuple(sorted(int(x) for x in t)) for t in v.facets()}
    edeg = Counter()
    for t in tets:
        for e in combinations(t, 2):
            edeg[e] += 1
    if not args.verify:
        tets = None

    pbox = CELL * args.mcell
    s.enable_event_log(args.logmb)
    s.drain_event_log()

    def chord_pos(chord, pos=None):
        """Midpoint of a chord in cells, min-imaged about its first endpoint.
        Two O(1) lookups -- no cochain marshal, no tree integration."""
        q = s.vertex_positions(chord).astype(float)
        d = (q[1] - q[0] + pbox / 2) % pbox - pbox / 2
        return (q[0] + d / 2.0) / CELL

    tracks = {}          # chord key -> Track
    finished = []
    tid = 0
    n_ev = n_slide = n_hop = 0
    overflow = False

    for c in knot_chords_now(v):
        tracks[c] = Track(tid, 0, c)
        tracks[c].samples.append((0, chord_pos(c)))
        tid += 1

    done = 0
    while done < args.sweeps:
        s.run(sweeps=args.chunk)
        done += args.chunk
        ev = s.drain_event_log()
        if s.event_log_overflowed():
            overflow = True
            print(f"  WARN event-log overflow at sweep {done} -- MSD is "
                  f"unreliable; raise --logmb or lower --chunk", flush=True)
        n_ev += len(ev)
        # --- full-resolution pass: maintain edge degrees, and follow every
        # tracked chord across any slide that moves it. This is the only thing
        # the event stream is needed for; species membership is settled at the
        # chunk boundary below, where the real state is available.
        for grp in group_by_clock(ev):
            sl = slide_chords(grp)
            for e in grp:
                rem, add = er.event_changes(e)
                if tets is not None:
                    for r in rem:
                        tets.discard(r)
                    tets.update(add)
                for t in rem:
                    for pp in combinations(t, 2):
                        edeg[pp] -= 1
                for t in add:
                    for pp in combinations(t, 2):
                        edeg[pp] += 1
            if sl is None:
                continue
            n_slide += 1
            old, new = sl
            tr = tracks.pop(old, None)
            if tr is None:
                continue                       # slid a chord we do not track
            if edeg.get(new, 0) == 3:
                tr.chord = new
                tr.hops += 1
                tracks[new] = tr
                n_hop += 1
            else:                              # landed without a handle
                tr.alive = False
                tr.died = done
                finished.append(tr)

        # --- chunk boundary: reconcile against the VERIFIED species set
        live = knot_chords_now(v)
        for c in list(tracks):
            if c in live:
                tracks[c].samples.append((done, chord_pos(c)))
                tracks[c].sig_seen[(3, 4, 4)] += 1
            else:
                tr = tracks.pop(c)
                tr.alive = False
                tr.died = done
                finished.append(tr)
        for c in live:
            if c not in tracks:
                tracks[c] = Track(tid, done, c)
                tracks[c].samples.append((done, chord_pos(c)))
                tracks[c].sig_seen[(3, 4, 4)] += 1
                tid += 1
        if tets is not None:
            actual = {tuple(sorted(int(x) for x in t)) for t in v.facets()}
            if tets != actual:
                print(f"  REPLAY DIVERGED at sweep {done}: "
                      f"symdiff {len(tets ^ actual)}", flush=True)
                sys.exit(1)

    finished.extend(tracks.values())
    from collections import Counter as _C
    print("worldline length histogram (samples/track):",
          dict(sorted(_C(len(t.samples) for t in finished).items())))
    print(f"replayed {n_ev} accepted moves over {done} sweeps; "
          f"{n_slide} slides seen, {n_hop} linked as hops"
          + ("  [OVERFLOW]" if overflow else ""), flush=True)

    # --- MSD over worldlines: min-image steps in cells, glitch-filtered
    disp = defaultdict(list)
    ntr = 0
    for tr in finished:
        if len(tr.samples) < 2:
            continue
        ntr += 1
        sw = [x[0] for x in tr.samples]
        raw = [x[1] for x in tr.samples]
        un = [np.zeros(3)]
        for i in range(1, len(raw)):
            step = raw[i] - raw[i - 1]
            step -= np.round(step / args.mcell) * args.mcell   # min image
            if np.linalg.norm(step) > GLITCH:
                tr.glitches += 1
                step = np.zeros(3)
            un.append(un[-1] + step)
        un = np.array(un)
        for i in range(len(un)):
            for j in range(i + 1, len(un)):
                disp[sw[j] - sw[i]].append(float(np.sum((un[j] - un[i]) ** 2)))

    rows = [{"dt": k, "msd": float(np.mean(vv)), "n": len(vv)}
            for k, vv in sorted(disp.items())]
    knots = sum(1 for tr in finished
                if tr.sig_seen and tr.sig_seen.most_common(1)[0][0] == (3, 4, 4))
    summary = {
        "slide_prob": args.slide_prob, "sweeps": done, "chunk": args.chunk,
        "events": n_ev, "slides": n_slide, "hops": n_hop,
        "tracks": ntr, "tracks_mostly_344": knots,
        "hops_per_track": (n_hop / ntr) if ntr else 0.0,
        "glitches": sum(tr.glitches for tr in finished),
        "overflow": overflow, "msd": rows,
    }
    with open(args.out + ".msd.json", "w") as f:
        json.dump(summary, f, indent=1)
    print(f"tracks {ntr} ({knots} predominantly (3,4,4)); hops/track "
          f"{summary['hops_per_track']:.2f}; glitches {summary['glitches']}")
    for r in rows[:8]:
        print(f"   dt={r['dt']:5d} sw   MSD={r['msd']:9.4f} cells^2  "
              f"(n={r['n']})")
    print(f"wrote {args.out}.msd.json")


if __name__ == "__main__":
    main()
