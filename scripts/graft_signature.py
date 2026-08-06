#!/usr/bin/env python3
"""Decorated-boundary signatures of tet lumps, and graft surgery.

A GRAFT replaces the interior of a lump L_B (a set of tets) in crystal B with
a lump L_A from crystal A, glued along a simplicial isomorphism
phi: bd(L_A) -> bd(L_B) of the boundary surfaces. Matching levels:

  L1  simplicial iso of the boundary 2-complexes.  The glued object is
      automatically a closed 3-manifold (each interface vertex link becomes
      disc + disc = S^2), but the seam can be arbitrarily defective.
  L2  L1 + phi preserves each boundary edge's (total degree, interior tet
      count d_in).  Then every interface edge keeps its host degree, so ALL
      edges of the graft stay in {5,6} (interior edges of either piece keep
      their native legal degrees).
  L3  L2 + phi preserves each boundary vertex's interior decoration
      (cn, n6, n_int_edges, n_int_6edges).  Then every vertex keeps its host
      coordination and 6-edge count: the disclination web reconnects across
      the seam and every vertex is a legal Z class.  A level-3 graft has
      exactly zero local energy cost under any degree-only action.

Certificates are EXACT canonical forms computed by combinatorial-map
traversal of the (decorated) boundary surface -- minimum over all starting
darts of a deterministic BFS labelling.  No WL/hash bucketing is ever the
arbiter (notes/memory/crystal-symmetry-group.md), and the canonical vertex
order realizing the minimum doubles as the isomorphism phi for surgery.

Library module; `python scripts/graft_signature.py --selftest` runs a
star-swap smoke test on the in-memory c15 m3 reference.
"""
import argparse
import os
import sys
from collections import Counter, deque

import numpy as np

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_ROOT, "python"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

TET_EDGES = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
TET_FACES = [(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)]


# ----------------------------------------------------------------- context

class CrystalContext:
    """Per-triangulation precomputation: edge degrees, vertex cn / n6."""

    def __init__(self, facets, name=""):
        F = np.asarray(facets, np.int64)
        self.F = F
        self.name = name
        E = np.vstack([F[:, [i, j]] for i, j in TET_EDGES])
        E.sort(axis=1)
        eu, cnt = np.unique(E, axis=0, return_counts=True)
        self.edge_deg = {(int(a), int(b)): int(c)
                         for (a, b), c in zip(eu, cnt)}
        nv = int(F.max()) + 1
        cn = np.zeros(nv, np.int64)
        n6 = np.zeros(nv, np.int64)
        np.add.at(cn, eu[:, 0], 1)
        np.add.at(cn, eu[:, 1], 1)
        six = eu[cnt == 6]
        np.add.at(n6, six[:, 0], 1)
        np.add.at(n6, six[:, 1], 1)
        self.cn = cn
        self.n6 = n6

    def star_of_vertex(self, v):
        """Tet indices containing vertex v (a closed-star lump)."""
        return np.nonzero((self.F == v).any(axis=1))[0]


# ------------------------------------------------------- boundary extraction

class SurfaceError(ValueError):
    """Lump boundary is not a closed 2-manifold (pinched edge/vertex)."""


def lump_boundary(ctx, lump_idx):
    """Extract the decorated boundary surface of a lump.

    Returns dict with:
      faces      list of 3-tuples (global vids), the boundary triangles
      edeco      {(a,b): (deg, d_in)} for surface edges
      vdeco      {v: (cn, n6, n_int, n_int6)} for surface vertices
      interior_vids  set of lump vertices not on the surface
      d_in       {(a,b): lump tet count} for every edge of a lump tet
    Raises SurfaceError if the boundary is not a closed surface.
    """
    F = ctx.F
    lump_idx = np.asarray(lump_idx, np.int64)
    L = F[lump_idx]

    # boundary faces: triangles in exactly one lump tet
    faces_all = np.vstack([L[:, list(f)] for f in TET_FACES])
    faces_all.sort(axis=1)
    fu, fcnt = np.unique(faces_all, axis=0, return_counts=True)
    if fcnt.max() > 2:
        raise SurfaceError("face in >2 lump tets -- not a manifold lump")
    bfaces = [tuple(int(x) for x in f) for f in fu[fcnt == 1]]
    if not bfaces:
        raise SurfaceError("lump has empty boundary (whole manifold?)")

    # d_in for every edge of a lump tet
    E = np.vstack([L[:, [i, j]] for i, j in TET_EDGES])
    E.sort(axis=1)
    eu, cnt = np.unique(E, axis=0, return_counts=True)
    d_in = {(int(a), int(b)): int(c) for (a, b), c in zip(eu, cnt)}

    # surface edges & their decoration; consistency check vs d_in split
    surf_edges = set()
    for a, b, c in bfaces:
        surf_edges |= {(a, b), (a, c), (b, c)}
    edeco = {}
    for e in surf_edges:
        deg = ctx.edge_deg[e]
        di = d_in.get(e, 0)
        if not (0 < di < deg):
            raise SurfaceError(f"surface edge {e} has d_in={di}, deg={deg}")
        edeco[e] = (deg, di)
    mixed = {e for e, di in d_in.items() if di < ctx.edge_deg[e]}
    if mixed != surf_edges:
        raise SurfaceError("mixed-edge set != surface-edge set (pinched edge)")

    # surface vertices: interior-edge decorations
    surf_vids = {v for f in bfaces for v in f}
    n_int = Counter()
    n_int6 = Counter()
    for (a, b), di in d_in.items():
        deg = ctx.edge_deg[(a, b)]
        if di == deg:  # interior edge (all its tets are in the lump)
            for v in (a, b):
                if v in surf_vids:
                    n_int[v] += 1
                    if deg == 6:
                        n_int6[v] += 1
    vdeco = {v: (int(ctx.cn[v]), int(ctx.n6[v]), n_int[v], n_int6[v])
             for v in surf_vids}

    lump_vids = set(int(x) for x in np.unique(L))
    interior_vids = lump_vids - surf_vids

    # surface manifoldness: every surface edge in exactly 2 boundary faces,
    # and the boundary faces around each vertex form a single cycle
    e2f = {}
    for fi, (a, b, c) in enumerate(bfaces):
        for e in ((a, b), (a, c), (b, c)):
            e2f.setdefault(e, []).append(fi)
    for e, fl in e2f.items():
        if len(fl) != 2:
            raise SurfaceError(f"surface edge {e} lies in {len(fl)} faces")
    v2f = {}
    for fi, f in enumerate(bfaces):
        for v in f:
            v2f.setdefault(v, []).append(fi)
    for v, fl in v2f.items():
        n_e = sum(1 for e in surf_edges if v in e)
        if n_e != len(fl):
            raise SurfaceError(f"vertex {v}: {n_e} surface edges, {len(fl)} faces")
        # single cycle <=> the faces at v are connected through edges at v
        seen = {fl[0]}
        stack = [fl[0]]
        while stack:
            fi = stack.pop()
            a, b, c = bfaces[fi]
            for e in ((a, b), (a, c), (b, c)):
                if v in e:
                    for g in e2f[e]:
                        if g not in seen:
                            seen.add(g)
                            stack.append(g)
        if len(seen) != len(fl):
            raise SurfaceError(f"vertex {v}: pinched surface (link not a cycle)")

    return dict(faces=bfaces, edeco=edeco, vdeco=vdeco, e2f=e2f,
                interior_vids=interior_vids, d_in=d_in)


# ------------------------------------------------------------ canonical form

def _components(faces, e2f):
    """Partition face indices into connected components."""
    n = len(faces)
    comp = [-1] * n
    ncomp = 0
    for s in range(n):
        if comp[s] >= 0:
            continue
        comp[s] = ncomp
        stack = [s]
        while stack:
            fi = stack.pop()
            a, b, c = faces[fi]
            for e in ((a, b), (a, c), (b, c)):
                for g in e2f[e]:
                    if comp[g] < 0:
                        comp[g] = ncomp
                        stack.append(g)
        ncomp += 1
    return comp, ncomp


def _traverse(faces, e2f, start_fi, a, b, vdec, edec):
    """Deterministic BFS labelling from dart (start_fi, a->b).

    Returns (cert, order): cert is a tuple of per-face records
    (lab_u, lab_v, lab_w, vdec u, vdec v, vdec w, edec uv, edec vw, edec wu)
    in visit order; order[i] = global vid with canonical label i.
    """
    lab = {a: 0, b: 1}
    order = [a, b]
    visited = {start_fi}
    rec = []
    q = deque([(start_fi, a, b)])
    while q:
        fi, u, v = q.popleft()
        x, y, z = faces[fi]
        w = x if x != u and x != v else (y if y != u and y != v else z)
        if w not in lab:
            lab[w] = len(order)
            order.append(w)
        rec.append((lab[u], lab[v], lab[w], vdec(u), vdec(v), vdec(w),
                    edec(u, v), edec(v, w), edec(w, u)))
        for p, r in ((u, v), (v, w), (w, u)):
            e = (p, r) if p < r else (r, p)
            f1, f2 = e2f[e]
            g = f2 if f1 == fi else f1
            if g not in visited:
                visited.add(g)
                q.append((g, r, p))
    return tuple(rec), order


def _canon_component(faces, e2f, face_ids, vdec, edec):
    """Exact canonical form of one surface component.

    Returns (cert, orders) where orders is the list of ALL canonical vertex
    orders achieving the minimal cert (one per symmetry of the decorated
    component -- alternative phis for surgery).
    """
    darts = []
    for fi in face_ids:
        x, y, z = faces[fi]
        for a, b in ((x, y), (y, x), (x, z), (z, x), (y, z), (z, y)):
            darts.append((fi, a, b))

    def dart_key(d):
        fi, a, b = d
        x, y, z = faces[fi]
        c = x if x not in (a, b) else (y if y not in (a, b) else z)
        return (vdec(a), vdec(b), vdec(c), edec(a, b), edec(b, c), edec(c, a))

    kmin = min(dart_key(d) for d in darts)
    best = None
    orders = []
    for d in darts:
        if dart_key(d) != kmin:
            continue
        cert, order = _traverse(faces, e2f, d[0], d[1], d[2], vdec, edec)
        if best is None or cert < best:
            best, orders = cert, [order]
        elif cert == best:
            orders.append(order)
    return best, orders


_LEVELS = (1, 2, 3)


def lump_signature(ctx, lump_idx, levels=_LEVELS):
    """Signature of a lump: per-level sorted component certificates.

    Returns dict:
      certs   {level: tuple of component certs, sorted}
      orders  {level: [orders-list per component, in cert-sorted order]}
      diag    dict of diagnostics (face/component counts, interior fingerprint)
    """
    bd = lump_boundary(ctx, lump_idx)
    faces, e2f = bd["faces"], bd["e2f"]
    comp, ncomp = _components(faces, e2f)
    comp_faces = [[i for i in range(len(faces)) if comp[i] == k]
                  for k in range(ncomp)]

    null_v = lambda v: ()
    null_e = lambda a, b: ()
    edec_l = lambda a, b: bd["edeco"][(a, b) if a < b else (b, a)]
    vdec_l = lambda v: bd["vdeco"][v]
    deco = {1: (null_v, null_e), 2: (null_v, edec_l), 3: (vdec_l, edec_l)}

    certs, orders = {}, {}
    for lvl in levels:
        vdec, edec = deco[lvl]
        pairs = [_canon_component(faces, e2f, cf, vdec, edec)
                 for cf in comp_faces]
        pairs.sort(key=lambda p: p[0])
        certs[lvl] = tuple(p[0] for p in pairs)
        orders[lvl] = [p[1] for p in pairs]

    # interior fingerprint: distinguishes interiors of boundary-isomorphic lumps
    lump_idx = np.asarray(lump_idx, np.int64)
    int_edges = sorted(ctx.edge_deg[e] for e, di in bd["d_in"].items()
                       if di == ctx.edge_deg[e])
    int_zcls = sorted((int(ctx.cn[v]), int(ctx.n6[v]))
                      for v in bd["interior_vids"])
    diag = dict(n_tets=len(lump_idx), n_faces=len(faces), n_comps=ncomp,
                comp_sizes=sorted(len(cf) for cf in comp_faces),
                n_int_vertices=len(bd["interior_vids"]),
                interior_fp=(len(lump_idx), tuple(Counter(int_edges).items()
                             if int_edges else ()),
                             tuple(sorted(Counter(int_zcls).items()))))
    return dict(certs=certs, orders=orders, diag=diag, boundary=bd)


# ------------------------------------------------------------------- matching

def match_phis(sig_donor, sig_host, level=3, max_phis=64):
    """Candidate boundary isomorphisms phi: donor bd vids -> host bd vids.

    Requires equal cert tuples at `level`.  Yields dicts; multiple candidates
    arise from decorated-surface symmetries (and are tried in turn when a
    graft needs an orientation-compatible phi).
    """
    ca, cb = sig_donor["certs"][level], sig_host["certs"][level]
    if ca != cb:
        return
    oa, ob = sig_donor["orders"][level], sig_host["orders"][level]
    # per component: donor's first canonical order against each host order
    from itertools import product
    combos = product(*[range(len(obs)) for obs in ob])
    n = 0
    for combo in combos:
        phi = {}
        for k, (oas, obs) in enumerate(zip(oa, ob)):
            for va, vb in zip(oas[0], obs[combo[k]]):
                phi[va] = vb
        yield phi
        n += 1
        if n >= max_phis:
            return


# -------------------------------------------------------------------- surgery

def graft(ctx_host, lump_host_idx, ctx_donor, lump_donor_idx, phi):
    """Replace the lump in the host with the donor lump, glued along phi.

    Returns (facets, info): the new facet array (vertex ids compacted) and a
    dict describing the surgery.  Asserts that phi maps the donor boundary
    faces exactly onto the host boundary faces.
    """
    bd_h = lump_boundary(ctx_host, lump_host_idx)
    bd_d = lump_boundary(ctx_donor, lump_donor_idx)
    host_bfaces = {tuple(sorted(f)) for f in bd_h["faces"]}
    img = {tuple(sorted(phi[v] for v in f)) for f in bd_d["faces"]}
    if img != host_bfaces:
        raise ValueError("phi does not map donor boundary onto host boundary")

    F_h = ctx_host.F
    keep = np.ones(len(F_h), bool)
    keep[np.asarray(lump_host_idx, np.int64)] = False
    host_part = F_h[keep]

    fresh = int(F_h.max()) + 1
    relab = dict(phi)
    for v in sorted(bd_d["interior_vids"]):
        relab[v] = fresh
        fresh += 1
    donor_part = np.array([[relab[int(v)] for v in tet]
                           for tet in ctx_donor.F[np.asarray(lump_donor_idx)]],
                          np.int64)
    newF = np.vstack([host_part, donor_part])
    newF.sort(axis=1)
    uniq = np.unique(newF, axis=0)
    if len(uniq) != len(newF):
        raise ValueError("duplicate facets after graft")
    # compact vertex labels
    lab, inv = np.unique(newF, return_inverse=True)
    newF = inv.reshape(newF.shape).astype(np.int64)
    info = dict(n_host_kept=int(keep.sum()), n_donor=len(donor_part),
                n_interior_new=len(bd_d["interior_vids"]),
                n_boundary=len(phi), n_vertices=len(lab))
    return newF, info


# ----------------------------------------------------------------- validation

def validate_facets(F, expect_z=None):
    """Full FK validation of a facet array.

    Checks: closed 3-manifold (link sum rule, via vertex_class_census's
    assert), edge-degree census, Z-class census + broken disclination count,
    Euler characteristic, orientability.  Returns a report dict.
    """
    from fk_skeleton import edges_from_facets, vertex_class_census
    from discrete_differential_geometry import Manifold

    eu, edeg, V = edges_from_facets(F)
    fz, n_broken = vertex_class_census(eu, edeg, V)  # asserts sum rule
    deg_census = dict(sorted(Counter(edeg.tolist()).items()))
    mfd = Manifold(3, np.asarray(F).tolist())
    rep = dict(n_vertices=V, n_tets=len(F),
               edge_degrees=deg_census,
               all_56=set(deg_census) <= {5, 6},
               z_census={k: round(v * V) for k, v in fz.items()},
               n_broken_disclination=n_broken,
               euler_characteristic=int(mfd.euler_characteristic),
               orientable=bool(mfd.is_orientable))
    return rep


# ------------------------------------------------------------------ self-test

def _selftest():
    from tcp_reference import build_t3_triangulation
    print("building c15 m3 in-memory ...")
    fac, nv = build_t3_triangulation("c15", 3)
    ctx = CrystalContext(fac, "c15m3")
    base = validate_facets(fac)
    print(f"  base: V={base['n_vertices']} tets={base['n_tets']} "
          f"degs={base['edge_degrees']} Z={base['z_census']} "
          f"chi={base['euler_characteristic']} orient={base['orientable']}")

    # two Z12 vertices (cn=12): star-swap
    z12 = np.nonzero(ctx.cn == 12)[0]
    v0, v1 = int(z12[0]), int(z12[-1])
    s0, s1 = ctx.star_of_vertex(v0), ctx.star_of_vertex(v1)
    sig0 = lump_signature(ctx, s0)
    sig1 = lump_signature(ctx, s1)
    print(f"  star({v0}): {sig0['diag']}")
    for lvl in _LEVELS:
        same = sig0["certs"][lvl] == sig1["certs"][lvl]
        print(f"  L{lvl} certs match star({v0}) vs star({v1}): {same}")
        assert same, "Z12 stars in the same crystal must match at every level"
    n_sym = len(sig0["orders"][3][0])
    print(f"  decorated-link symmetries of Z12 star: {n_sym}")

    ok = False
    for phi in match_phis(sig1, sig0, level=3):
        newF, info = graft(ctx, s0, ctx, s1, phi)
        rep = validate_facets(newF)
        good = (rep["all_56"] and rep["n_broken_disclination"] == 0
                and rep["euler_characteristic"] == 0 and rep["orientable"]
                and rep["z_census"] == base["z_census"])
        print(f"  graft: {info} -> all56={rep['all_56']} "
              f"broken={rep['n_broken_disclination']} chi="
              f"{rep['euler_characteristic']} orient={rep['orientable']} "
              f"Z={rep['z_census']}")
        if good:
            ok = True
            break
    assert ok, "star-swap graft failed validation"
    print("SELFTEST PASS")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--selftest", action="store_true")
    args = ap.parse_args()
    if args.selftest:
        _selftest()
    else:
        ap.print_help()
