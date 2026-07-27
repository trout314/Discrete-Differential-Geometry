#!/usr/bin/env python3
"""Interactions between defect species, string <-> vacuum and string <-> string.

The statics problem: persistent complexes are rare (a handful per box) and
sessile (caged, defect_travel), so snapshot pair statistics cannot resolve
species-species interactions. But the FLICKER VACUUM refreshes every ~2
sweeps, so one long run yields thousands of independent samples of where
births happen RELATIVE to the standing strings. That is a species-resolved
interaction measurement with real statistics:

  g_bs(r)     flicker-birth rate vs distance to each string species,
              normalised by a MATCHED null (uniform points scored against the
              same instantaneous string set -- exact same geometry, no
              interaction by construction).
  P2 profile  births resolved by angle to the string's principal axis
              (gyration tensor of its member positions): the directional
              part of the string-vacuum coupling.
  absorption  which deaths are really merges into a string, at what birth
              distance -- simultaneously the string GROWTH PATHWAY: do
              strings grow by eating flickers or by internal conversion?
  lifetime    flicker lifetime vs birth distance to the nearest string.

String := worldline alive with age > --age-min sweeps OR modal n_verts >=
--big-nv. Species label = modal induced f-vector over the track's life.

String-string pair separations are also recorded (every --ss-every sweeps)
but the strings barely move, so a run contributes ~one independent
configuration; those numbers are reported with that caveat attached.

FRAME: registry lift (extrinsic; CONVENTIONS.md sec 6).
Exponents/class verdicts frame-robust; lengths gauge. P2-type
correlators compare DIRECTIONS AT SEPARATED POINTS: by sec 6
there is no intrinsic global direction field (holonomy), so
they are registry-gauge in an ESSENTIAL way, not merely in
their axis labels (the nonzero P2 null along crystal axes is
one symptom). Honest upgrade: transported-frame correlators
(development.TransportContext), with the transport path
stated.
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
from fk_skeleton import edges_from_facets

CELL = 1.0e6


def minimg(d, box):
    return d - np.round(d / box) * box


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--cell", required=True)
    ap.add_argument("--start", default=None)
    ap.add_argument("--startcoc", required=True)
    ap.add_argument("--mcell", type=int, default=4)
    ap.add_argument("--lam", type=float, required=True)
    ap.add_argument("--etarget", type=float, default=None)
    ap.add_argument("--zleg", type=float, default=0.6)
    ap.add_argument("--cimp", type=float, default=1.0)
    ap.add_argument("--slide-prob", type=float, default=0.0)
    ap.add_argument("--sweeps", type=int, default=2500)
    ap.add_argument("--chunk", type=int, default=25)
    ap.add_argument("--age-min", type=float, default=15.0)
    ap.add_argument("--big-nv", type=int, default=7)
    ap.add_argument("--null-k", type=int, default=8)
    ap.add_argument("--ss-every", type=int, default=5,
                    help="sweeps between string-string separation samples")
    ap.add_argument("--logmb", type=float, default=256.0)
    ap.add_argument("--seed", type=int, default=4242)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    print("[frame] registry lift -- P2 profile essentially gauge (sec 6)")

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
    e0, w0, _ = coc.load_cocycle(args.startcoc)
    s.enable_cocycle(np.asarray(e0), np.asarray(w0))
    s.check_cocycle()
    box = float(args.mcell)
    pbox = CELL * args.mcell
    s.enable_cocycle_positions([int(pbox)] * 3)
    s.check_cocycle_positions()
    n3 = v.num_facets

    st = ds.DefectState(v)
    wl = ds.Worldlines(st)
    s.enable_event_log(args.logmb)
    s.drain_event_log()
    rng = np.random.default_rng(args.seed + 1)

    def centroid_axis(verts):
        q = s.vertex_positions(verts).astype(float)
        rel = minimg(q - q[0], pbox)
        cen = ((q[0] + rel.mean(0)) / CELL) % box
        rc = (rel - rel.mean(0)) / CELL
        if len(verts) >= 3:
            w_, vec = np.linalg.eigh(rc.T @ rc)
            axis = vec[:, -1]
            elong = float(w_[-1] / max(w_.sum(), 1e-12))
        else:
            axis, elong = None, 0.0
        return cen, axis, elong

    # per-track bookkeeping
    born_clock, verts_of, shapes_of, nvmax = {}, {}, {}, {}
    birth_pos = {}
    pair_rec = []          # flicker-birth vs string records
    null_rec = []
    absorb = []            # (dead_tid, absorber_tid, r_at_birth or nan)
    death_clock = {}
    ss_rec = []
    last_ss = -1

    def strings_now(clock, skip_tid=None):
        out = []
        for t2, vv in verts_of.items():
            if t2 == skip_tid or t2 in death_clock:
                continue
            age = (clock - born_clock[t2]) / n3
            if age > args.age_min or nvmax.get(t2, 0) >= args.big_nv:
                cen, axis, elong = centroid_axis(vv)
                out.append((t2, cen, axis, elong))
        return out

    def record(clock):
        alive_now = set()
        for cx, tid in wl.step(clock):
            alive_now.add(tid)
            new = tid not in born_clock
            if new:
                born_clock[tid] = clock
            verts_of[tid] = list(cx.verts)
            f = tuple(st.induced_shape(cx.verts)["f"])
            shapes_of.setdefault(tid, Counter())[f] += 1
            nvmax[tid] = max(nvmax.get(tid, 0), len(cx.verts))
            if new:
                cen, _, _ = centroid_axis(cx.verts)
                birth_pos[tid] = cen
                for (t2, c2, ax, el) in strings_now(clock, skip_tid=tid):
                    d = minimg(np.array(cen) - c2, box)
                    r = float(np.linalg.norm(d))
                    c2ang = (float(np.dot(d / max(r, 1e-12), ax)) ** 2
                             if ax is not None and r > 1e-9 else np.nan)
                    pair_rec.append((tid, t2, clock, r, c2ang, el))
                    for u in rng.uniform(0, box, (args.null_k, 3)):
                        dn = minimg(u - c2, box)
                        rn = float(np.linalg.norm(dn))
                        cn = (float(np.dot(dn / max(rn, 1e-12), ax)) ** 2
                              if ax is not None and rn > 1e-9 else np.nan)
                        null_rec.append((t2, clock, rn, cn, el))
        # deaths + absorption
        for t2 in [t for t in verts_of if t not in alive_now
                   and t not in death_clock]:
            death_clock[t2] = clock
            heir = Counter(wl.label.get(x) for x in verts_of[t2]
                           if wl.label.get(x) is not None)
            if heir:
                ab = heir.most_common(1)[0][0]
                if ab in alive_now:
                    d = minimg(np.array(birth_pos[t2])
                               - np.array(birth_pos.get(ab, birth_pos[t2])),
                               box)
                    absorb.append((t2, ab, float(np.linalg.norm(d))))

    record(0)
    done = 0
    while done < args.sweeps:
        s.run(sweeps=args.chunk)
        done += args.chunk
        for e in s.drain_event_log():
            st.apply(e)
            clock = int(e["clock"])
            record(clock)
            if clock / n3 - last_ss >= args.ss_every:
                last_ss = clock / n3
                ss = strings_now(clock)
                for i in range(len(ss)):
                    for j in range(i + 1, len(ss)):
                        d = minimg(ss[i][1] - ss[j][1], box)
                        ss_rec.append((ss[i][0], ss[j][0],
                                       float(np.linalg.norm(d))))
        if s.event_log_overflowed():
            print("  WARN event-log overflow", flush=True)
        print(f"  sw{done}: {len(born_clock)} worldlines, "
              f"{len(pair_rec)} birth-string pairs", flush=True)

    end_clock = done * n3

    def species(tid):
        return tuple(shapes_of[tid].most_common(1)[0][0])

    def lifef(tid):
        return (death_clock.get(tid, end_clock) - born_clock[tid]) / n3

    out = {
        "sweeps": done, "n3": n3, "lam": args.lam, "box": box,
        "age_min": args.age_min, "big_nv": args.big_nv,
        "null_k": args.null_k, "seed": args.seed,
        "merges": wl.merges, "splits": wl.splits,
        "tracks": {str(t): {"born": born_clock[t] / n3,
                            "life": lifef(t), "species": list(species(t)),
                            "nvmax": nvmax[t],
                            "censored": t not in death_clock}
                   for t in born_clock},
        "pairs": [{"tid": a, "string": b, "sw": c / n3, "r": r,
                   "cos2": None if np.isnan(c2) else c2, "elong": el}
                  for (a, b, c, r, c2, el) in pair_rec],
        "null": [{"string": b, "sw": c / n3, "r": r,
                  "cos2": None if np.isnan(c2) else c2, "elong": el}
                 for (b, c, r, c2, el) in null_rec],
        "absorb": [{"dead": a, "into": b, "r_birth": r}
                   for (a, b, r) in absorb],
        "ss": [{"a": a, "b": b, "r": r} for (a, b, r) in ss_rec],
    }
    with open(args.out, "w") as fh:
        json.dump(out, fh)
    print(f"\n{len(born_clock)} worldlines; {len(pair_rec)} birth-string "
          f"pairs; {len(absorb)} absorptions; {len(ss_rec)} string-string "
          f"samples")
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
