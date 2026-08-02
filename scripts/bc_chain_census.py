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
      Requires the tcp_reference labelling v = cell*ns + site; the script
      VALIDATES that (every edge must be shorter than half a cell) and refuses
      to report winding rather than print a wrong number.

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
      Handedness flips under ORIENTATION-REVERSING automorphisms. Aut(K) is
      combinatorial and may contain them; Aut+ (orientation-preserving) has
      index 1 or 2. A chain class is ACHIRAL iff its stabilizer contains an
      orientation-reversing element -- then the mirror image of the chain is
      the chain itself, up to Aut+. If the stabilizer is entirely inside Aut+
      the class is CHIRAL: its Aut-orbit splits into two Aut+-orbits, a
      left- and a right-handed enantiomer. (Equivalently: orbit size under
      Aut+ is halved exactly when the stabilizer loses no element.) If Aut has
      no orientation-reversing element at all, the crystal itself is chiral
      and no enantiomer is present to compare with -- reported separately.

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

from discrete_differential_geometry import CrystalSymmetry
from discrete_differential_geometry.development import (
    develop_chain, chain_step_quat)


def parse_name(path):
    """(structure, m) from a tcp_reference filename, or (None, None)."""
    mm = re.match(r"T3_([A-Z0-9]+)_m(\d+)_N\d+\.mfd", os.path.basename(path))
    return (mm.group(1).lower(), int(mm.group(2))) if mm else (None, None)


def reference_positions(struct, m, nvert):
    """Fractional torus coordinates (period m) per vertex id, or None.

    Reproduces tcp_reference's site perturbation and id scheme exactly; see
    cocycle_check.reference_frac_positions.
    """
    from tcp_reference import STRUCTURES
    if struct not in STRUCTURES:
        return None
    L, sites, _, _ = STRUCTURES[struct]
    ns = len(sites)
    if ns * m ** 3 != nvert:
        return None
    rng = np.random.default_rng(12345)
    sites = sites + 1e-6 * rng.standard_normal(sites.shape)
    v = np.arange(nvert)
    c = v // ns
    return (sites[v % ns] + np.stack([c // (m * m), (c // m) % m, c % m],
                                     axis=1)) % m


def minimg(d, m):
    return d - m * np.round(d / m)


def validate_positions(rp, sym, m, tol=1e-9):
    """Certify that minimum-image reconstructs the TRUE edge displacements.

    Stored positions give each vertex's spot in the fundamental domain, not
    which periodic image an edge runs to, so every displacement is
    RECONSTRUCTED as the minimum image. That is correct only if the true
    displacement already lies within (-m/2, m/2) componentwise; otherwise the
    reconstruction is off by exactly m on that edge and the winding is off by
    exactly +-1 -- still a tidy integer vector, just the wrong one.

    Note that "the minimum image is under m/2" is NOT a test of this: minimg
    returns components in that range by construction, so it would only ever
    catch an exact tie. The real certificate is CLOSURE: the reconstructed
    displacement is a 1-cochain, and the true one is closed, so if the
    reconstruction sums to zero around every triangle it agrees with the truth
    up to a cocycle -- and a single mis-wrapped edge breaks closure on every
    triangle containing it. Returns (worst edge component, closure residual,
    margin = (m/2)/worst, ok).

    Closure alone leaves the theoretical loophole of a coherent mis-wrapped
    CUT SURFACE (a nonzero cocycle), which is why the margin is reported too:
    at margin >> 1 the alternative wrap would make an edge (m - worst) long,
    many times the longest edge that actually occurs, so no such surface
    exists.
    """
    view, lab = sym.view, sym.view.labels
    ext = lab.astype(np.int64)

    eu = np.array([[ext[u], ext[w]] for (u, w) in view.edges])
    de = minimg(rp[eu[:, 1]] - rp[eu[:, 0]], m)
    worst = float(np.abs(de).max())

    fa = np.array([[ext[a], ext[b], ext[c]] for (a, b, c) in view.faces])
    cyc = (minimg(rp[fa[:, 1]] - rp[fa[:, 0]], m)
           + minimg(rp[fa[:, 2]] - rp[fa[:, 1]], m)
           + minimg(rp[fa[:, 0]] - rp[fa[:, 2]], m))
    closure = float(np.abs(cyc).max())

    margin = (m / 2.0) / worst if worst > 0 else float("inf")
    return worst, closure, margin, (closure < tol and margin > 2.0)


def chain_winding(sym, chain, rp, m):
    """Winding of ONE chain, in supercell periods (its class in H_1(T^3))."""
    view, lab = sym.view, sym.view.labels
    seq = [int(lab[view.frame_window(int(f))[0]]) for f in chain]
    p = rp[seq]
    d = minimg(np.roll(p, -1, axis=0) - p, m)
    return d.sum(axis=0) / m


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


_PARITY4 = {}


def _perm_sign(p):
    """Sign of a permutation given as a 4-tuple of distinct indices."""
    if p not in _PARITY4:
        s, seen = 1, list(p)
        for i in range(4):
            for j in range(i + 1, 4):
                if seen[i] > seen[j]:
                    s = -s
        _PARITY4[p] = s
    return _PARITY4[p]


def orientation_parity(view):
    """parity[t] in {0,1}: the tet's coherently oriented vertex order is its
    sorted order, or sorted with the last two swapped when parity is 1.

    Adjacent tets are consistent iff they induce OPPOSITE orientations on the
    shared face; dropping vertex i from a sorted 4-tuple leaves a sorted face,
    so the induced sign is just (-1)^i. Same rule as
    development.TransportContext._orient, on TriView's arrays.
    """
    st = np.sort(np.asarray(view.tets), axis=1)
    keys = [tuple(int(x) for x in st[t]) for t in range(view.nT)]
    parity = np.full(view.nT, -1, np.int8)
    parity[0] = 0
    stack = [0]
    while stack:
        t = stack.pop()
        kt = keys[t]
        for i in range(4):
            face = kt[:i] + kt[i + 1:]
            u = [x for x in view.face2[face] if x[0] != t][0][0]
            ku = keys[u]
            iu = ku.index(next(v for v in ku if v not in face))
            want = parity[t] ^ ((i ^ iu ^ 1) & 1)
            if parity[u] < 0:
                parity[u] = want
                stack.append(u)
            elif parity[u] != want:
                raise ValueError("non-orientable complex")
    if (parity < 0).any():
        raise ValueError("disconnected complex")
    return st, keys, parity


def _canon_order(key, par):
    return (key[0], key[1], key[3], key[2]) if par else key


def orientation_sign(view, keys, parity, g):
    """+1 if automorphism `g` preserves the global orientation, -1 if it
    reverses it. Well defined: the sign is the same at every tet."""
    o = _canon_order(keys[0], int(parity[0]))
    img = tuple(int(g[v]) for v in o)
    t2 = view.tetid[tuple(sorted(img))]
    o2 = _canon_order(keys[t2], int(parity[t2]))
    return _perm_sign(tuple(o2.index(v) for v in img))


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
        worst, closure, margin, ok = validate_positions(rp, sym, m)
        if not ok:
            print(f"  WARNING: displacement reconstruction NOT certified "
                  f"(triangle-closure residual {closure:.2e}, wrap margin "
                  f"{margin:.2f}x); winding numbers SUPPRESSED rather than "
                  f"reported wrong.")
            rp = None
        else:
            print(f"  positions: {struct} m={m} -- reconstruction certified: "
                  f"triangle closure exact to {closure:.1e}, worst edge "
                  f"{worst:.4f} vs wrap scale m/2 = {m/2:.1f} "
                  f"({margin:.1f}x margin)")
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

    _, keys, parity = orientation_parity(view)
    signs = [orientation_sign(view, keys, parity, g) for g in sym.elements]
    n_rev = sum(1 for s in signs if s < 0)
    if n_rev:
        print(f"  orientation: Aut+ has index 2 ({sym.order - n_rev} preserving, "
              f"{n_rev} reversing) -- mirror images are comparable")
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

        shifts, stab_signs = set(), []
        for g, sg in zip(sym.elements, signs):
            img = view.frame_id(tuple(int(g[v]) for v in w0))
            if int(chain_of[img]) == reps[k]:
                shifts.add(pos_in[img])
                stab_signs.append(sg)
        screw = min((s for s in shifts if s), default=L)
        # achiral iff some orientation-REVERSING automorphism fixes the chain
        mirror = (None if not n_rev else any(s < 0 for s in stab_signs))

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
