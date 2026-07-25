#!/usr/bin/env python3
"""What defect species does a SINGLE 2->3 move on a pristine crystal make?

At fixed N the sampler's cheapest excitation is a 2->3 immediately undone by
the reciprocal 3->2. Every such pair injects a defect complex for a few moves
and then removes it. That is the "flicker background": species that are not
excitations of the melt at all, just the move set breathing against a perfect
crystal. Any census taken at per-move resolution is dominated by it, so we
need to know exactly what it looks like before reading anything into a
species table.

WHY THE ANSWER IS A SHORT LIST. A 2->3 on triangle abc with apexes x,y changes
the degree of exactly nine edges -- ab,ac,bc (-1), the six side edges (+1) --
and creates xy at degree 3. Every one of those has both endpoints inside
{a,b,c,x,y}, so NO vertex outside the five can change its edge degrees, its
n6, or its legality. The resulting defect is therefore always an induced
subcomplex on a subset of five vertices, and the species is fixed by the local
combinatorics of (abc, x, y) in the crystal. Inequivalent triangle orbits give
inequivalent species; there are only a handful.

METHOD. Seed a DefectState from the pristine crystal once, then for each
candidate triangle apply the synthetic 2->3 event, read the species off, and
apply the reciprocal 3->2 to restore. Both are O(1) incremental updates, so
the whole enumeration is linear in the number of triangles rather than
quadratic. The restoration is verified, not assumed.
"""
import argparse
import json
import os
import sys
from collections import Counter, defaultdict
from itertools import combinations

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
for _p in ("../../python", "../../scripts", "../../tools", "."):
    sys.path.insert(0, os.path.join(_HERE, _p))
import discrete_differential_geometry as ddg
from discrete_differential_geometry.move_geometry import EVENT_DTYPE
import defect_state as ds

THETA = np.arccos(1.0 / 3.0)
UNIT = np.pi - 3 * THETA


def load_flicker(path, group="full"):
    """The pristine-crystal species list as a set of keys at one grouping.

    Returns (keyset, meta). Membership is tested with json.dumps so tuples and
    the lists that survive a JSON round trip compare equal.

    WHAT MEMBERSHIP MEANS, and does not: a species in this set is one the
    pristine crystal CAN make with a single 2->3. It does not follow that any
    particular observed instance was made that way -- a static snapshot has no
    move history. So the flicker column on a snapshot census is an UPPER BOUND
    on flicker contamination. For a real attribution, which needs the birth
    move, use flicker_fraction.py on an event stream.
    """
    fl = json.load(open(path))
    keys = set()
    for sp in fl["species"]:
        for c in sp["components"]:
            k = c.get("_keys", {}).get(group)
            if k is None and group == "full":
                k = c.get("_key")
            if k is not None:
                keys.add(json.dumps(k))
    return keys, {"cell": fl.get("cell"), "sites": fl.get("sites"),
                  "n3": fl.get("n3"), "species": len(fl["species"])}


def flicker_key(group, cx=None, facets=None, vcolor=None, ecolor=None):
    """The key of an observed complex, built exactly as load_flicker's are."""
    if group == "sig":
        return json.dumps([list(cx.sig), list(cx.nodes)])
    if group == "full":
        k, _ = ds.canonical_key(facets, ecolor=ecolor, vcolor=vcolor)
    elif group == "decorated":
        k, _ = ds.canonical_key(facets, ecolor=ecolor)
    else:
        k, _ = ds.canonical_key(facets)
    return json.dumps(k)


def _ev(mtype, labels):
    """A synthetic event record in the sampler's own EVENT_DTYPE."""
    e = np.zeros(1, dtype=EVENT_DTYPE)[0]
    e["type"] = mtype
    lab = list(labels) + [0] * (6 - len(labels))
    for i, v in enumerate(lab):
        e["labels"][i] = v
    return e


def triangle_apexes(state):
    """Every triangle of the complex -> its two apex vertices."""
    tri = defaultdict(list)
    for t in state.tets:
        for f in combinations(t, 3):
            tri[f].append([v for v in t if v not in f][0])
    return tri


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--cell", required=True)
    ap.add_argument("--limit", type=int, default=0,
                    help="stop after this many triangles (0 = all)")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    m = ddg.Manifold.load(args.cell, 3)
    st = ds.DefectState(m)
    print(f"{args.cell}: N3 = {m.num_facets}, "
          f"{len(st.defect)} defect vertices before any move")
    if st.defect:
        print("  NOT PRISTINE -- the species below would be contaminated")
        return 1
    edeg0 = Counter(st.edeg)
    Z = {v: len(st.neighbours(v)) for v in st.v2t}
    zc = Counter(Z.values())
    print(f"  coordination census: "
          f"{', '.join(f'Z{k}x{v}' for k, v in sorted(zc.items()))}")
    print(f"  edge degrees: "
          f"{', '.join(f'{k}:{v}' for k, v in sorted(Counter(edeg0.values()).items()))}")

    tri = triangle_apexes(st)
    print(f"  {len(tri)} triangles")

    species = defaultdict(lambda: {"n": 0, "rep": None})
    skipped_existing_edge = 0
    ntried = 0
    for f, ap_ in tri.items():
        if len(ap_) != 2:
            continue
        x, y = sorted(ap_)
        if ds._ek(x, y) in st.edeg:
            # the 2->3 would duplicate an existing edge: not a valid move
            skipped_existing_edge += 1
            continue
        ntried += 1
        st.apply(_ev(1, list(f) + [x, y]))
        comps = st.components(broad=True)
        rec = []
        for cx in comps:
            fac = st.induced_facets(cx.verts)
            vc, ec = st.decorations(cx.verts)
            ck, exact = ds.canonical_key(fac, ecolor=ec, vcolor=vc)
            # Keys at EVERY grouping level species_report.py offers, because a
            # coarser key cannot be projected out of the full one -- dropping a
            # decoration changes which vertex ordering minimises the form. A
            # consumer grouping by --group decorated must compare against keys
            # built the same way, or the membership test silently misses.
            keys = {
                "full": ck,
                "decorated": ds.canonical_key(fac, ecolor=ec)[0],
                "induced": ds.canonical_key(fac)[0],
                "sig": [list(cx.sig), list(cx.nodes)],
            }
            sumz = st.total_coordination(cx.verts)
            rec.append({
                "key": ck, "keys": keys, "exact": exact, "nv": len(cx.verts),
                "sig": cx.sig, "nodes": cx.nodes, "sumZ": sumz,
                "Q_c": sumz * UNIT + 6 * len(cx.verts) * THETA,
                "f": st.induced_shape(cx.verts)["f"],
                "n6": tuple(sorted(st.n6[v] for v in cx.verts)),
                "Zs": tuple(sorted(len(st.neighbours(v)) for v in cx.verts)),
            })
        # species label for the WHOLE excitation (a move may split the
        # defect into more than one component, so key on the multiset)
        label = tuple(sorted(r["key"] for r in rec))
        s = species[label]
        s["n"] += 1
        if s["rep"] is None:
            s["rep"] = rec
        st.apply(_ev(2, [x, y] + list(f)))       # reciprocal 3->2 restores
        if args.limit and ntried >= args.limit:
            break

    bad = [e for e, d in st.edeg.items() if d != edeg0.get(e, 0)]
    print(f"\nrestoration check: {len(bad)} edges differ from pristine, "
          f"{len(st.defect)} defect vertices remain "
          f"-> {'OK' if not bad and not st.defect else 'FAILED'}")
    print(f"{ntried} valid 2->3 sites, {skipped_existing_edge} triangles "
          f"skipped (apexes already joined by an edge)")

    print(f"\n{len(species)} distinct species from a single 2->3:")
    rows = sorted(species.items(), key=lambda kv: -kv[1]["n"])
    print(f"  {'count':>7} {'share':>7}  {'ncomp':>5} {'f-vector':>16} "
          f"{'sig':>14} {'sumZ':>5} {'Q_c':>8}  {'n6':>16} {'Z':>18}")
    out = []
    for label, s in rows:
        rec = s["rep"]
        r0 = rec[0]
        print(f"  {s['n']:7d} {100*s['n']/ntried:6.2f}%  {len(rec):5d} "
              f"{str(tuple(r0['f'])):>16} {str(tuple(r0['sig'])):>14} "
              f"{r0['sumZ']:5d} {r0['Q_c']:8.3f}  {str(r0['n6']):>16} "
              f"{str(r0['Zs']):>18}")
        if len(rec) > 1:
            for r in rec[1:]:
                print(f"  {'':7} {'':7}  {'':5} {str(tuple(r['f'])):>16} "
                      f"{str(tuple(r['sig'])):>14} {r['sumZ']:5d} "
                      f"{r['Q_c']:8.3f}  {str(r['n6']):>16} "
                      f"{str(r['Zs']):>18}")
        out.append({"count": s["n"], "share": s["n"] / ntried,
                    "components": [
                        dict({k: (list(v) if isinstance(v, tuple) else v)
                              for k, v in r.items()
                              if k not in ("key", "keys")},
                             # the keys themselves, so downstream scripts can
                             # test observed defects for membership in this
                             # pristine list at whatever grouping they use
                             _key=r["key"], _keys=r["keys"])
                        for r in rec]})

    print("\nsumZ distribution over 2->3 sites (the flicker's charge ladder):")
    zdist = Counter()
    for label, s in rows:
        for r in s["rep"]:
            zdist[r["sumZ"]] += s["n"]
    tot = sum(zdist.values())
    for k in sorted(zdist):
        print(f"   sumZ {k:3d}  Q_c {k*UNIT + 30*THETA:+7.3f}  "
              f"{zdist[k]:7d}  {100*zdist[k]/tot:6.2f}%")

    if args.out:
        with open(args.out, "w") as fh:
            json.dump({"cell": args.cell, "n3": m.num_facets,
                       "sites": ntried, "species": out}, fh, indent=1)
        print(f"\nwrote {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
