#!/usr/bin/env python3
"""Census of the BC (Boerdijk-Coxeter) chain classes of a crystal triangulation.

A BC chain is the closed sliding-window walk of ``worm_helix.bc_orbit``:
window k = v[k..k+3], stepping to (v1,v2,v3, other apex of face v1v2v3). The
step map is a BIJECTION on ordered tets ("frames"), so the 24*nT frames
partition into disjoint cycles -- the COMPLETE set of directed BC chains,
enumerated in one linear pass by ``symmetry.CrystalSymmetry``. Aut(K) acts on
that set; the orbits are the chain CLASSES, and they are what a run should
name. Almost every script in ``defect_dynamics/`` instead takes the single
chain through ``F[0]``, whose class is neither chosen nor recorded.

Per class this reports:

  LENGTH / ORBIT / |Stab|
      L steps; how many chains are in the class; the stabilizer order.

  SCREW s
      The chain's own symmetry. An automorphism commutes with the sliding-
      window map, so any g fixing the chain setwise sends window k to window
      k+s for a fixed shift s -- the realized shifts form a subgroup of Z_L.
      The generator s is the crystallographic screw: the chain maps to itself
      advanced by s windows, and L/s = |Stab| such elements.

  HOLONOMY phi_L
      The EXACT geometric holonomy of the closed chain, from
      ``development.py``'s rational embedding. In the developed picture every
      step is the same screw motion S (this is what makes it a tetrahelix), so
      developing once around a chain that closes combinatorially returns the
      window rotated by R^L, R = rotation part of S with cos phi = -2/3
      exactly. By Niven's theorem cos phi = -2/3 is rational while phi/pi is
      NOT, so R^L is never the identity: a BC chain closes combinatorially on
      the torus and never closes geometrically. The holonomy is the exact
      measure of that frustration. cos phi_L is a rational with denominator
      3^L, so it is reported as a float with the exact form noted.

  WINDING
      The chain's class in H_1(T^3) = Z^3, in supercell lattice-vector units,
      summed from minimum-image steps of the reference fractional positions.
      Reported as the SET of windings realized across the class: translations
      preserve winding but the point group rotates it, so a class has a
      winding orbit, not a winding. The minimum-image reconstruction is
      certified by cochain closure before any of it is reported (see
      ``validate_positions``); if it fails, winding is suppressed rather than
      guessed.

  COVERAGE
      How many tet / vertex orbits the chain visits, out of all of them. A
      chain hitting every tet orbit is a single closed curve that samples
      every symmetry-distinct local environment in the crystal.

  REV (direction reversal, NOT chirality)
      Reversal (v0,v1,v2,v3) -> (v3,v2,v1,v0) conjugates the step map to its
      inverse, so it permutes chains. "rev" says whether a class is closed
      under traversing its chains backwards. This is NOT handedness: a helix
      walked from the other end has the SAME handedness, so reversal is not a
      mirror.

  MIRROR (this one is chirality)
      Handedness flips under ORIENTATION-REVERSING automorphisms, so the
      relevant group is Aut+ = the orientation-preserving subgroup (index 1 or
      2). A class is CHIRAL iff its Aut-orbit splits into two Aut+-orbits --
      a left- and a right-handed enantiomer that the full Aut identifies -- and
      ACHIRAL iff some orientation-reversing automorphism fixes it. If Aut has
      no orientation-reversing element at all the crystal itself is chiral
      (a Sohncke space group) and no enantiomer is present to compare with;
      that is reported separately. Computed by
      ``CrystalSymmetry.is_chiral``/``.orientation_preserving``.

Usage:
  python scripts/bc_chain_census.py data/tcp_reference/T3_R_m3_N24462.mfd
  python scripts/bc_chain_census.py <mfd> --json chains.json
"""
import argparse
import collections
import json
import os
import re
import sys
from fractions import Fraction

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
sys.path.insert(0, os.path.join(_ROOT, "python"))
sys.path.insert(0, _HERE)
sys.path.insert(0, os.path.join(_ROOT, "tools"))

from discrete_differential_geometry import CrystalSymmetry
from discrete_differential_geometry.development import (
    develop_chain, chain_step_quat)
# position / certification / winding helpers live in tools/chain_select.py so
# the drivers and this census share ONE implementation
from chain_select import (parse_ref_name as parse_name, reference_positions,
                          certify_positions as validate_positions, minimg,
                          chain_winding)


def class_windings(sym, chains, member_ids, rp, m):
    """The DISTINCT windings realized across a whole chain class.

    Lattice translations preserve winding but the point-group part of Aut
    ROTATES it, so a class does not have "a" winding -- reporting one
    representative's vector would misrepresent the class. Returns the sorted
    set of integer winding vectors and the worst distance from integrality
    (a diagnostic: it must be ~0 or the position labelling is wrong).
    """
    seen, resid = set(), 0.0
    for i in member_ids:
        w = chain_winding(sym, chains[i], rp, m)
        resid = max(resid, float(np.abs(w - np.round(w)).max()))
        seen.add(tuple(int(x) for x in np.round(w)))
    return sorted(seen), resid


def chain_holonomy(sym, chain):
    """(cos phi_L as an exact Fraction, phi_L in degrees, step axis).

    Develops only the first five windows -- the step screw is constant along
    the chain -- and raises R to the chain length. `chain_step_quat` asserts
    cos phi = -2/3, which is itself a check that these really are BC chains.
    """
    view = sym.view
    lab = view.labels
    seq = [int(lab[view.frame_window(int(chain[k % len(chain)]))[0]])
           for k in range(5)]
    if len(set(seq)) != 5:
        return None
    pos = develop_chain(seq)
    q = chain_step_quat(seq, pos)
    hol = q ** len(chain)
    c = hol.cos_phi()
    return c, float(np.degrees(np.arccos(max(-1.0, min(1.0, float(c)))))), \
        q.axis()


def reverse_frame(view, fid):
    w = view.frame_window(fid)
    return view.frame_id((w[3], w[2], w[1], w[0]))


def main(path, jsonout, no_cache):
    sym = CrystalSymmetry.for_manifold_path(path, cache=not no_cache)
    view = sym.view
    lab = view.labels
    _, chain_of, chains = sym._chain_tables()
    oid, members, reps = sym.chain_orbits()

    struct, m = parse_name(path)
    print(f"{os.path.basename(path)}: V={view.V} tets={view.nT}  "
          f"|Aut| = {sym.order}")
    print(f"  orbits: vertices {sym.n_orbits('vertex')} edges "
          f"{sym.n_orbits('edge')} faces {sym.n_orbits('face')} "
          f"tets {sym.n_orbits('tet')}")
    print(f"  BC chains: {len(chains)} directed cycles in {len(members)} "
          f"classes")

    rp = reference_positions(struct, m, view.V) if struct else None
    if rp is not None:
        worst, closure, alt, ok = validate_positions(rp, sym, m)
        if not ok:
            print(f"  WARNING: displacement reconstruction NOT certified "
                  f"(triangle-closure residual {closure:.2e}, alt-wrap ratio "
                  f"{alt:.2f}x); winding numbers SUPPRESSED rather than "
                  f"reported wrong.")
            rp = None
        else:
            print(f"  positions: {struct} m={m} -- certified: triangle "
                  f"closure {closure:.1e}; a mis-wrapped edge would be "
                  f"{alt:.1f}x the longest real one ({worst:.4f})")
    elif struct:
        print(f"  positions: unavailable for {struct} m={m}; winding suppressed")
    else:
        print(f"  positions: filename is not a tcp_reference crystal; "
              f"winding suppressed")

    # per-tet / per-vertex orbit id, in internal ids
    tmap = sym.orbit_id_map("tet")
    vmap = sym.orbit_id_map("vertex")
    tet_orb = np.array([tmap[tuple(sorted(int(lab[x]) for x in view.tets[t]))]
                        for t in range(view.nT)])
    vert_orb = np.array([vmap[int(lab[v])] for v in range(view.V)])
    n_tet_orb, n_vert_orb = sym.n_orbits("tet"), sym.n_orbits("vertex")

    # which class holds the default seed used across defect_dynamics/
    f0 = view.frame_id(tuple(int(view.index[int(x)])
                             for x in np.asarray(view.tets[0])))
    default_class = int(oid[int(chain_of[f0])])

    n_rev = (sym.order - sym.orientation_preserving.order
             if sym.has_orientation_reversing else 0)
    if n_rev:
        print(f"  orientation: Aut+ has index 2 "
              f"({sym.orientation_preserving.order} preserving, {n_rev} "
              f"reversing) -- mirror images are comparable")
    else:
        print(f"  orientation: every automorphism preserves orientation, so "
              f"the CRYSTAL is chiral; no chain has its mirror present to "
              f"compare with (mirror column = 'n/a')")

    rows = []
    for k, mem in enumerate(members):
        c = chains[reps[k]]
        L = len(c)
        pos_in = {int(f): i for i, f in enumerate(c)}
        w0 = view.frame_window(int(c[0]))

        shifts = set()
        for g in sym.elements:
            img = view.frame_id(tuple(int(g[v]) for v in w0))
            if int(chain_of[img]) == reps[k]:
                shifts.add(pos_in[img])
        screw = min((s for s in shifts if s), default=L)
        chiral = sym.is_chiral("chain", reps[k])
        mirror = None if chiral is None else not chiral

        tets = np.array([int(f) // 24 for f in c])
        verts = np.unique(view.tets[tets].ravel())
        tcov = len(np.unique(tet_orb[tets]))
        vcov = len(np.unique(vert_orb[verts]))

        rev_class = int(oid[int(chain_of[reverse_frame(view, int(c[0]))])])
        hol = chain_holonomy(sym, c)
        wind, resid = (class_windings(sym, chains, mem, rp, m)
                       if rp is not None else (None, None))

        rows.append(dict(
            cls=k, L=L, n_chains=len(mem), stab=sym.order // len(mem),
            screw=screw, shifts=len(shifts),
            n_tets=len(set(tets.tolist())), n_verts=len(verts),
            tet_orbits=tcov, vert_orbits=vcov,
            covers_all_tet_orbits=bool(tcov == n_tet_orb),
            covers_all_vert_orbits=bool(vcov == n_vert_orb),
            reverse_class=rev_class, reverse_closed=bool(rev_class == k),
            achiral=mirror,
            holonomy_deg=(hol[1] if hol else None),
            holonomy_cos_num_digits=(len(str(hol[0].numerator)) if hol else None),
            step_axis=(list(hol[2]) if hol else None),
            windings=([list(w) for w in wind] if wind is not None else None),
            n_distinct_windings=(len(wind) if wind is not None else None),
            winding_residual=resid,
            is_default_seed_class=bool(k == default_class)))

    print(f"\n  {'#':>3} {'L':>6} {'chains':>7} {'|stab|':>7} {'screw':>6} "
          f"{'tets':>6} {'tetorb':>7} {'vorb':>6} {'rev':>7} {'mirror':>7} "
          f"{'holonomy':>9}  windings in class")
    for r in rows:
        if r["windings"] is None:
            wd = "-"
        else:
            ws = ["(" + ",".join(f"{x:+d}" for x in w) + ")"
                  for w in r["windings"]]
            wd = (f"{len(ws)}x " if len(ws) > 1 else "") + " ".join(ws[:3]) \
                + (" ..." if len(ws) > 3 else "")
        mark = "  <- bc_orbit(F[0]) lands here" if r["is_default_seed_class"] else ""
        rv = "self" if r["reverse_closed"] else f"->#{r['reverse_class']}"
        mi = ("n/a" if r["achiral"] is None
              else ("achiral" if r["achiral"] else "CHIRAL"))
        print(f"  {r['cls']:>3} {r['L']:>6} {r['n_chains']:>7} {r['stab']:>7} "
              f"{r['screw']:>6} {r['n_tets']:>6} "
              f"{r['tet_orbits']}/{n_tet_orb:<5} "
              f"{r['vert_orbits']}/{n_vert_orb:<4} "
              f"{rv:>7} {mi:>7} {r['holonomy_deg']:>8.3f}d  {wd}{mark}")

    full = [r for r in rows if r["covers_all_tet_orbits"]]
    print(f"\n  chain classes visiting EVERY tet orbit ({n_tet_orb}): "
          + (", ".join(f"#{r['cls']} (L={r['L']})" for r in full)
             if full else "none"))
    fullv = [r for r in rows if r["covers_all_vert_orbits"]]
    print(f"  chain classes visiting EVERY vertex orbit ({n_vert_orb}): "
          + (", ".join(f"#{r['cls']} (L={r['L']})" for r in fullv)
             if fullv else "none"))
    nrev = sum(1 for r in rows if r["reverse_closed"])
    print(f"  direction reversal: {nrev} class(es) closed under reversal, "
          f"{len(rows) - nrev} swapped in pairs (NOT a handedness statement)")
    if n_rev:
        nach = sum(1 for r in rows if r["achiral"])
        nplus = sum(1 if r["achiral"] else 2 for r in rows)
        print(f"  chirality (mirror): {nach} achiral class(es), "
              f"{len(rows) - nach} chiral (Aut-orbit splits into an "
              f"enantiomer pair under Aut+)")
        print(f"  => {len(rows)} classes under Aut, but {nplus} under Aut+ "
              f"(handedness-resolved). Use the Aut count for anything the "
              f"sampler sees;\n     use the Aut+ count for handedness-"
              f"sensitive quantities (Wilson-line signs, the BC-helix Dirac "
              f"walk).")
    else:
        print(f"  chirality (mirror): undefined here -- the crystal has no "
              f"orientation-reversing automorphism at all")
    print(f"  holonomy: nonzero for every class (cos phi = -2/3 per step is "
          f"rational with phi/pi irrational, so R^L != I for all L)")

    if jsonout:
        with open(jsonout, "w") as fh:
            json.dump({"crystal": os.path.basename(path), "aut_order": sym.order,
                       "n_tet_orbits": n_tet_orb, "n_vertex_orbits": n_vert_orb,
                       "n_chains": len(chains), "classes": rows}, fh, indent=1)
        print(f"\n  wrote {os.path.abspath(jsonout)}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("crystal")
    ap.add_argument("--json", dest="jsonout", default=None)
    ap.add_argument("--no-cache", action="store_true")
    a = ap.parse_args()
    p = a.crystal if os.path.isabs(a.crystal) else os.path.join(_ROOT, a.crystal)
    main(p, a.jsonout, a.no_cache)
