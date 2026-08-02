#!/usr/bin/env python3
"""How much of an observed defect census is 2->3 flicker against clean crystal?

`crystal_flicker.py` says what a single 2->3 on a PRISTINE crystal produces:
always the f = (5,10,9,3) knot shape, in a short list of decorated species.
This script asks the converse question of a real run -- of the defect
worldlines actually born, how many are that, and how many are something the
crystal cannot make in one move?

Each worldline is attributed by three independent marks:

  birth move type   a flicker complex is born by a 2->3 and nothing else;
  death move type   a pure flicker is undone by the reciprocal 3->2;
  species key       the full decorated key (induced subcomplex + ambient edge
                    degrees + n6) must be one the pristine crystal can make.

The species key is the strict mark, and deliberately so: it carries the
AMBIENT degrees, so a 2->3 fired into a region that is already disturbed
produces a key the pristine enumeration never contains. Splitting births into
"2->3 on clean crystal" and "2->3 on disturbed background" is the point --
the first is pure background and carries no information about the melt, the
second is the melt reacting.
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

MOVE = {0: "1->4", 1: "2->3", 2: "3->2", 3: "4->1", 4: "4-4"}


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--cell", required=True)
    ap.add_argument("--start", default=None)
    ap.add_argument("--startcoc", default=None)
    ap.add_argument("--lam", type=float, required=True)
    ap.add_argument("--etarget", type=float, default=None)
    ap.add_argument("--zleg", type=float, default=0.6)
    ap.add_argument("--cimp", type=float, default=1.0)
    ap.add_argument("--slide-prob", type=float, default=0.0)
    ap.add_argument("--sweeps", type=int, default=400)
    ap.add_argument("--chunk", type=int, default=25)
    ap.add_argument("--logmb", type=float, default=256.0)
    ap.add_argument("--flicker", required=True,
                    help="crystal_flicker.py JSON for the same crystal")
    ap.add_argument("--seed", type=int, default=99)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    # the pristine-crystal species list, as a set of decorated keys
    fl = json.load(open(args.flicker))
    flick_keys = set()
    for sp in fl["species"]:
        for c in sp["components"]:
            flick_keys.add(json.dumps(c["_key"]) if "_key" in c else None)
    flick_keys.discard(None)
    if not flick_keys:
        print("flicker JSON has no keys -- rerun crystal_flicker.py "
              "(it must be the version that stores '_key')")
        return 1
    print(f"{len(flick_keys)} pristine-crystal species keys loaded")

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

    born = {}          # tid -> dict
    alive = set()
    cur_type = None

    def record(clock):
        seen = set()
        for cx, tid in wl.step(clock):
            seen.add(tid)
            if tid not in born:
                fac = st.induced_facets(cx.verts)
                vc, ec = st.decorations(cx.verts)
                ck = ds.canonical_key_exact(fac, ecolor=ec, vcolor=vc)
                born[tid] = {
                    "birth_move": cur_type, "clock": clock,
                    "nv": len(cx.verts), "f": list(st.induced_shape(
                        cx.verts)["f"]),
                    "sig": list(cx.sig),
                    "sumZ": st.total_coordination(cx.verts),
                    "flicker_key": json.dumps(ck) in flick_keys,
                    "death_move": None, "last": clock}
                alive.add(tid)
            else:
                born[tid]["last"] = clock
        for tid in list(alive):
            if tid not in seen:
                born[tid]["death_move"] = cur_type
                alive.discard(tid)

    record(0)
    done = 0
    while done < args.sweeps:
        s.run(sweeps=args.chunk)
        done += args.chunk
        for e in s.drain_event_log():
            cur_type = MOVE.get(int(e["type"]), "?")
            st.apply(e)
            record(int(e["clock"]))
        print(f"  sw{done}: {len(born)} worldlines, "
              f"{len(st.defect)} defect vertices", flush=True)

    rows = list(born.values())
    for r in rows:
        r["life_sweeps"] = (r["last"] - r["clock"]) / n3
    tot = len(rows)
    print(f"\n{tot} worldlines born in {done} sweeps\n")

    bt = Counter(r["birth_move"] for r in rows)
    print("birth move type:")
    for k, c in bt.most_common():
        print(f"   {str(k):>6} {c:7d} {100*c/tot:6.2f}%")
    dt = Counter(r["death_move"] for r in rows)
    print("death move type (None = still alive at the end):")
    for k, c in dt.most_common():
        print(f"   {str(k):>6} {c:7d} {100*c/tot:6.2f}%")

    knot = [r for r in rows if tuple(r["f"]) == (5, 10, 9, 3)]
    print(f"\nshape f = (5,10,9,3): {len(knot)} ({100*len(knot)/tot:.1f}%) "
          f"-- the only shape a single 2->3 on pristine crystal can make")

    nfk = sum(r["flicker_key"] for r in rows)
    print(f"\nATTRIBUTION")
    print(f"   born by 2->3 AND key in the pristine list "
          f"(flicker on clean crystal): "
          f"{sum(1 for r in rows if r['birth_move'] == '2->3' and r['flicker_key']):7d}"
          f" {100*sum(1 for r in rows if r['birth_move']=='2->3' and r['flicker_key'])/tot:6.2f}%")
    print(f"   born by 2->3, key NOT in the list "
          f"(2->3 into disturbed background):   "
          f"{sum(1 for r in rows if r['birth_move'] == '2->3' and not r['flicker_key']):7d}"
          f" {100*sum(1 for r in rows if r['birth_move']=='2->3' and not r['flicker_key'])/tot:6.2f}%")
    print(f"   born by any other move:                                       "
          f"     {sum(1 for r in rows if r['birth_move'] != '2->3'):7d}"
          f" {100*sum(1 for r in rows if r['birth_move']!='2->3')/tot:6.2f}%")
    print(f"   (total key matches to the pristine list: {nfk}, "
          f"{100*nfk/tot:.2f}%)")

    pure = [r for r in rows if r["birth_move"] == "2->3"
            and r["death_move"] == "3->2" and r["flicker_key"]]
    print(f"\n   strictest: 2->3 birth, 3->2 death, pristine key: "
          f"{len(pure)} ({100*len(pure)/tot:.2f}%)  "
          f"median life {np.median([r['life_sweeps'] for r in pure]):.3f} sweeps"
          if pure else "\n   strictest: none")

    print(f"\nlifetime by attribution class:")
    for lab, sel in (("flicker on clean crystal",
                      [r for r in rows if r["birth_move"] == "2->3"
                       and r["flicker_key"]]),
                     ("2->3 into disturbed background",
                      [r for r in rows if r["birth_move"] == "2->3"
                       and not r["flicker_key"]]),
                     ("other birth move",
                      [r for r in rows if r["birth_move"] != "2->3"])):
        if sel:
            L = [r["life_sweeps"] for r in sel]
            print(f"   {lab:35s} n={len(sel):6d}  median {np.median(L):8.3f}"
                  f"  mean {np.mean(L):9.3f}  max {max(L):9.2f}")

    if args.out:
        with open(args.out, "w") as fh:
            json.dump({"sweeps": done, "n3": n3, "lam": args.lam,
                       "tracks": rows}, fh)
        print(f"\nwrote {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
