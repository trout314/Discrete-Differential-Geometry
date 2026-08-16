#!/usr/bin/env python3
"""The halo: how far from a defect does the crystal stop being the crystal?

Per-vertex differencing against the pristine reference, made possible by the
cocycle lift: legal melt vertices sit ON reference sites (build_from_positions
uses p = round(scale*frac), and moves preserve the lattice registration up to
one global gauge offset per snapshot, which is estimated and removed). Each
melt vertex is matched to its reference site by exact integer position, so
"deviation from crystal" is measured site-for-site, immune to the site-class
texture that a naive shell average would read as halo:

  mismatch(v)  legal vertex whose (Z, n6) class differs from the reference
               class of its site -- FK-legal but chemically wrong, the
               silent field invisible to the defect census
  dq(v)        q_R(melt vertex) - q_R(reference site): the local screening
               charge carried by legal vertices

Profiles reported:
  global    P(mismatch | d) and <dq>(d) vs distance to the NEAREST defect
            vertex of any complex -- the crystal-recovery length, no
            exclusion rules needed
  species   tube profiles (distance to nearest member vertex) around each
            complex species, with contaminated shells dropped (a shell that
            contains another complex's vertex is discarded for that complex)
  cumulative screening charge vs r per species, against the complex's own
            excess charge dQ_cx -- where the image charge saturates IS the
            screening length
  axis      charge/mismatch split along vs perpendicular to the string's
            principal axis (elongation-gated): the directional-halo test

Matching quality is printed per snapshot (fraction of vertices matched
exactly); a poor match invalidates the profiles and aborts.
"""
import argparse
import glob as globmod
import json
import os
import sys
from collections import Counter, defaultdict

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
for _p in ("../../python", "../../scripts", "../../tools",
           "../../scripts/defect_dynamics", "."):   # archived: siblings stayed in dd/
    sys.path.insert(0, os.path.join(_HERE, _p))
import discrete_differential_geometry as ddg
from discrete_differential_geometry import cocycle as coc
from cocycle_check import reference_frac_positions
import defect_state as ds

CELL = 1.0e6


def minimg(d, box):
    return d - np.round(d / box) * box


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--glob", required=True, action="append")
    ap.add_argument("--cell", required=True)
    ap.add_argument("--mcell", type=int, default=4)
    ap.add_argument("--match-tol", type=float, default=2.0,
                    help="raw units; exact match is 0-1 after rounding")
    ap.add_argument("--steiner", type=int, default=0,
                    help="also bin the GLOBAL recovery profile by the "
                         "INTRINSIC PL geodesic distance (unit-edge "
                         "metric, Steiner order N on the actual "
                         "snapshot triangulation; CONVENTIONS.md sec "
                         "6). 0 = registry-only (legacy)")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    box = float(args.mcell)
    pbox = CELL * args.mcell

    # ---- pristine reference: site -> (Z, n6, q) by integer raw position
    ref = ddg.Manifold.load(args.cell, 3)
    st0 = ds.DefectState(ref)
    q0 = st0.vertex_charges()
    frac = np.asarray(reference_frac_positions("r", args.mcell))
    raw0 = np.round(CELL * frac).astype(np.int64) % int(pbox)
    site = {}
    for v in st0.v2t:
        Z = len(st0.neighbours(v))
        site[tuple(raw0[v])] = (Z, st0.n6[v], q0[v])
    assert len(site) == len(st0.v2t), "reference sites not unique"
    qbar0 = float(np.mean([q0[v] for v in st0.v2t]))
    print(f"reference: {len(site)} sites, qbar = {qbar0:.4f}")

    files = sorted({f for g in args.glob for f in globmod.glob(g)})
    dg_edges = np.linspace(0.0, 2.0, 11)
    dg_mid = 0.5 * (dg_edges[1:] + dg_edges[:-1])

    glob_mm = defaultdict(lambda: [0, 0])       # bin -> [mismatch, count]
    glob_dq = defaultdict(list)                 # bin -> dq values
    # intrinsic (PL geodesic) global profile, unit-edge bins
    pl_edges = np.linspace(0.0, 4.0, 17)
    pl_mid = 0.5 * (pl_edges[1:] + pl_edges[:-1])
    pl_mm = defaultdict(lambda: [0, 0])
    pl_dq = defaultdict(list)
    tube = {}      # species-group -> bin -> dict(sums)
    cum = {}       # species-group -> list of (r_sorted dq arrays) per complex
    axis_prof = defaultdict(lambda: defaultdict(lambda: [0.0, 0, 0]))
    dQcx_of = defaultdict(list)
    nsnap = 0

    def group(f):
        if f == (5, 10, 9, 3):
            return "knot"
        nv = f[0]
        return "small(nv<=6)" if nv <= 6 else "string(nv>=7)"

    for f in files:
        cocf = f.replace(".mfd", ".cocycle.npz")
        if not os.path.exists(cocf):
            continue
        m = ddg.Manifold.load(f, 3)
        st = ds.DefectState(m)
        q = st.vertex_charges()
        e1, w1, _ = coc.load_cocycle(cocf)
        e1 = np.asarray(e1)
        pos = coc.tree_positions(e1, np.asarray(w1), int(e1.max()) + 1)[0]
        verts = sorted(st.v2t)
        # ---- gauge offset: modal (pos - raw0[label]) over persisting labels
        cand = [v for v in verts if v < len(raw0)]
        offs = (pos[cand] - raw0[cand]) % int(pbox)
        off = np.array([Counter(map(tuple, offs)).most_common(1)[0][0]])
        matched, klass = {}, {}
        nmiss = 0
        for v in verts:
            key = tuple((pos[v] - off[0]) % int(pbox))
            hit = site.get(key)
            if hit is None:
                nmiss += 1
            matched[v] = hit
        mfrac = 1 - nmiss / len(verts)
        print(f"  {os.path.basename(f)}: matched {100*mfrac:.2f}% "
              f"({nmiss} miss)", flush=True)
        if mfrac < 0.99:
            sys.exit("site matching failed -- gauge offset wrong?")
        nsnap += 1

        pv = pos[verts] / CELL % box
        vidx = {v: i for i, v in enumerate(verts)}
        defect = sorted(st.defect)
        dpos = pos[defect] / CELL % box
        comps = st.components(broad=True)
        legal = [v for v in verts if v not in st.defect]

        # per-vertex fields on legal vertices
        mm, dq = {}, {}
        for v in legal:
            hit = matched[v]
            if hit is None:
                mm[v], dq[v] = 1, q[v] - qbar0
                continue
            Z = len(st.neighbours(v))
            mm[v] = int((Z, st.n6[v]) != (hit[0], hit[1]))
            dq[v] = q[v] - hit[2]

        # ---- global recovery profile
        lp = pv[[vidx[v] for v in legal]]
        d = minimg(lp[:, None, :] - dpos[None, :, :], box)
        dnear = np.sqrt((d ** 2).sum(-1)).min(1)
        for v, dn in zip(legal, dnear):
            b = np.searchsorted(dg_edges, dn) - 1
            if 0 <= b < len(dg_mid):
                glob_mm[b][0] += mm[v]
                glob_mm[b][1] += 1
                glob_dq[b].append(dq[v])

        # ---- intrinsic global profile (PL geodesic, unit-edge metric)
        if args.steiner:
            import steiner_geodesic as sg
            from scipy.sparse import coo_matrix
            from scipy.sparse.csgraph import dijkstra
            Fm = np.asarray(m.facets())
            relab = vidx
            Tarr = np.array([[relab[int(x)] for x in f4] for f4 in Fm],
                            dtype=np.int64)
            G, Ntot = sg.build_steiner_graph(Tarr, len(verts),
                                             args.steiner)
            # super-source at index Ntot linked to every defect vertex
            Gc = G.tocoo()
            eps = 1e-9
            src = np.array([relab[v] for v in defect])
            r = np.concatenate([Gc.row, np.full(len(src), Ntot), src])
            c = np.concatenate([Gc.col, src, np.full(len(src), Ntot)])
            w = np.concatenate([Gc.data, np.full(2 * len(src), eps)])
            Ga = coo_matrix((w, (r, c)),
                            shape=(Ntot + 1, Ntot + 1)).tocsr()
            dall = dijkstra(Ga, directed=False, indices=Ntot)
            for v in legal:
                dn = float(dall[relab[v]])
                b = np.searchsorted(pl_edges, dn) - 1
                if 0 <= b < len(pl_mid):
                    pl_mm[b][0] += mm[v]
                    pl_mm[b][1] += 1
                    pl_dq[b].append(dq[v])

        # ---- per-species tube profiles
        for cx in comps:
            gk = group(tuple(st.induced_shape(cx.verts)["f"]))
            mem = pv[[vidx[v] for v in cx.verts]]
            other = [v for v in defect if v not in set(cx.verts)]
            opos = pv[[vidx[v] for v in other]] if other else None
            dm = minimg(lp[:, None, :] - mem[None, :, :], box)
            dtube = np.sqrt((dm ** 2).sum(-1)).min(1)
            # contaminated shells: bins containing another complex's vertex
            bad = set()
            if opos is not None and len(opos):
                do = minimg(opos[:, None, :] - mem[None, :, :], box)
                dob = np.sqrt((do ** 2).sum(-1)).min(1)
                for x in dob:
                    b = np.searchsorted(dg_edges, x) - 1
                    if 0 <= b < len(dg_mid):
                        bad.add(b)
            T = tube.setdefault(gk, defaultdict(lambda: [0.0, 0, 0, 0]))
            rows = []
            for v, dt in zip(legal, dtube):
                b = np.searchsorted(dg_edges, dt) - 1
                if 0 <= b < len(dg_mid) and b not in bad:
                    T[b][0] += dq[v]
                    T[b][1] += mm[v]
                    T[b][2] += 1
                    T[b][3] = T[b][3]
                if dt <= 2.0:
                    rows.append((dt, dq[v]))
            # cumulative screening charge for this complex (uncontaminated
            # complexes only -- any bad shell poisons the cumulative sum)
            dQcx = sum(q[v] for v in cx.verts) \
                - sum((matched[v][2] if matched[v] else qbar0)
                      for v in cx.verts)
            dQcx_of[gk].append(dQcx)
            if not bad:
                rows.sort()
                cum.setdefault(gk, []).append(
                    (np.array([r_ for r_, _ in rows]),
                     np.cumsum([x for _, x in rows]), dQcx))
            # axis-resolved (elongated complexes only)
            if len(cx.verts) >= 3:
                rel = minimg(mem - mem[0], box)
                rc = rel - rel.mean(0)
                w_, vec = np.linalg.eigh(rc.T @ rc)
                if w_[-1] / max(w_.sum(), 1e-12) >= 0.5:
                    axv = vec[:, -1]
                    cen = (mem[0] + rel.mean(0)) % box
                    dc = minimg(lp - cen, box)
                    rr = np.sqrt((dc ** 2).sum(-1))
                    with np.errstate(invalid="ignore"):
                        c2 = (dc @ axv) ** 2 / np.maximum(rr ** 2, 1e-12)
                    for v, r_, cc in zip(legal, rr, c2):
                        if r_ > 2.0:
                            continue
                        b = np.searchsorted(dg_edges, r_) - 1
                        if not (0 <= b < len(dg_mid)) or b in bad:
                            continue
                        sect = "axis" if cc > 0.5 else "equator"
                        A = axis_prof[(gk, sect)][b]
                        A[0] += dq[v]
                        A[1] += mm[v]
                        A[2] += 1
        del st, m

    # ---------------- report ----------------
    print(f"\n{nsnap} snapshots")
    print("\nGLOBAL RECOVERY: distance to nearest defect vertex (cells)")
    print(f"   {'d':>5} {'P(mismatch)':>12} {'<dq>':>9} {'sum dq':>9} "
          f"{'n':>7}")
    grows = []
    for b in range(len(dg_mid)):
        mmc = glob_mm[b]
        dqs = np.array(glob_dq[b]) if glob_dq[b] else np.array([0.0])
        grows.append({"d": float(dg_mid[b]),
                      "p_mismatch": mmc[0] / max(mmc[1], 1),
                      "mean_dq": float(dqs.mean()),
                      "sum_dq": float(dqs.sum()) / max(nsnap, 1),
                      "n": mmc[1]})
        print(f"   {dg_mid[b]:5.2f} {grows[-1]['p_mismatch']:12.5f} "
              f"{grows[-1]['mean_dq']:9.5f} {grows[-1]['sum_dq']:9.4f} "
              f"{mmc[1]:7d}")

    pl_rows = []
    if args.steiner:
        print(f"\nGLOBAL RECOVERY, INTRINSIC: PL geodesic distance "
              f"(unit-edge lengths; Steiner order {args.steiner}, "
              f"actual triangulation)")
        print(f"   {'d':>5} {'P(mismatch)':>12} {'<dq>':>9} "
              f"{'sum dq':>9} {'n':>7}")
        for b in range(len(pl_mid)):
            mmc = pl_mm[b]
            dqs = np.array(pl_dq[b]) if pl_dq[b] else np.array([0.0])
            pl_rows.append({"d": float(pl_mid[b]),
                            "p_mismatch": mmc[0] / max(mmc[1], 1),
                            "mean_dq": float(dqs.mean()),
                            "sum_dq": float(dqs.sum()) / max(nsnap, 1),
                            "n": mmc[1]})
            print(f"   {pl_mid[b]:5.2f} "
                  f"{pl_rows[-1]['p_mismatch']:12.5f} "
                  f"{pl_rows[-1]['mean_dq']:9.5f} "
                  f"{pl_rows[-1]['sum_dq']:9.4f} {mmc[1]:7d}")

    trows = {}
    for gk, T in tube.items():
        trows[gk] = []
        print(f"\nTUBE PROFILE around {gk}  (uncontaminated shells)")
        print(f"   {'d':>5} {'P(mismatch)':>12} {'<dq>':>9} {'n':>7}")
        for b in range(len(dg_mid)):
            s = T.get(b)
            if not s or s[2] == 0:
                continue
            trows[gk].append({"d": float(dg_mid[b]),
                              "p_mismatch": s[1] / s[2],
                              "mean_dq": s[0] / s[2], "n": s[2]})
            print(f"   {dg_mid[b]:5.2f} {s[1]/s[2]:12.5f} "
                  f"{s[0]/s[2]:9.5f} {s[2]:7d}")

    crows = {}
    print("\nCUMULATIVE screening charge vs r (mean over uncontaminated "
          "complexes) against -dQ_cx:")
    for gk, lst in cum.items():
        rg = np.linspace(0.2, 2.0, 10)
        vals = []
        for rr, cs, dQcx in lst:
            vals.append([float(cs[np.searchsorted(rr, r_) - 1])
                         if np.searchsorted(rr, r_) > 0 else 0.0
                         for r_ in rg])
        v = np.array(vals)
        tgt = -np.mean([x[2] for x in lst])
        crows[gk] = {"r": rg.tolist(), "cum": v.mean(0).tolist(),
                     "cum_se": (v.std(0) / np.sqrt(len(v))).tolist(),
                     "target": float(tgt), "n_complexes": len(lst),
                     "dQcx_mean": float(-tgt)}
        print(f"   {gk} (n={len(lst)}, -<dQ_cx> = {tgt:+.3f}):")
        for i, r_ in enumerate(rg):
            print(f"      r={r_:4.2f}  cum dq = {v.mean(0)[i]:+8.4f} "
                  f"+- {v.std(0)[i]/np.sqrt(len(v)):6.4f}")

    arows = {}
    print("\nAXIS vs EQUATOR (elongated complexes):")
    for (gk, sect), P in sorted(axis_prof.items()):
        key = f"{gk}:{sect}"
        arows[key] = []
        for b in range(len(dg_mid)):
            s = P.get(b)
            if not s or s[2] == 0:
                continue
            arows[key].append({"d": float(dg_mid[b]),
                               "p_mismatch": s[1] / s[2],
                               "mean_dq": s[0] / s[2], "n": s[2]})
        line = " ".join(f"{r['d']:.1f}:{r['p_mismatch']:.4f}"
                        for r in arows[key][:8])
        print(f"   {key:24s} {line}")

    with open(args.out, "w") as fh:
        json.dump({"snapshots": nsnap, "box": box,
                   "global": grows, "global_pl": pl_rows,
                   "tube": trows, "cumulative": crows,
                   "axis": arows,
                   "dQcx": {k: [float(x) for x in v]
                            for k, v in dQcx_of.items()}}, fh, indent=1)
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
