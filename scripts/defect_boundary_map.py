#!/usr/bin/env python3
"""Exact boundary distance map of the elementary 2->3 defect ball.

Computes d_B^R and d_B^D on a shared order-n Steiner grid of dB (six unit
equilateral triangles, identical in both configurations), validates them,
and reports the difference map -- the complete, compactly supported
metric signature of one 2->3 move. See
``discrete_differential_geometry/ball_boundary.py`` for the exactness
argument; every distance here is sqrt(exact Fraction) in unit-edge units.

The map is UNIVERSAL (all tets are congruent regular unit simplices), so
this table applies to every 2->3 in every triangulation -- it is not
sampled data and carries no action, couplings, or chains.

Usage:
  python scripts/defect_boundary_map.py --order 6
  python scripts/defect_boundary_map.py --order 8 --out data/figs
"""
import os, sys, json, argparse, itertools
from fractions import Fraction

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
sys.path.insert(0, os.path.join(_ROOT, "python"))

from discrete_differential_geometry.ball_boundary import (
    BallBoundaryMap, BOUNDARY, VERTEX_NAMES, boundary_surface_lengths,
    node_key, node_name, node_faces, face_name, steiner_interior,
    CONFIG_SYMMETRY, A, B, C, P, Q)

SQ83 = float(Fraction(8, 3)) ** 0.5          # apex-apex chord in R


def build(order):
    mR = BallBoundaryMap("R", order)
    mD = BallBoundaryMap("D", order)
    assert mR.keys == mD.keys, "boundary grids must coincide"
    MR, missR = mR.matrix()
    MD, missD = mD.matrix()
    return mR, mD, MR, MD, missR, missD


def validate(mR, MR, MD, missR, missD, order):
    ok = True
    def check(name, cond, detail=""):
        nonlocal ok
        ok &= bool(cond)
        print(f"  [{'PASS' if cond else 'FAIL'}] {name}" +
              (f"   {detail}" if detail else ""))

    check("every strip chord resolved (no minimizer touches PQ)",
          not missR and not missD,
          f"unresolved: R {len(missR)}, D {len(missD)}")
    check("symmetric", np.allclose(MR, MR.T) and np.allclose(MD, MD.T))
    check("zero diagonal", not MR.diagonal().any() and not MD.diagonal().any())

    idx = {k: i for i, k in enumerate(mR.keys)}
    def v(x):
        return idx[node_key({x: Fraction(1)})]
    edges = [(A, B), (A, C), (B, C), (A, P), (B, P), (C, P), (A, Q), (B, Q), (C, Q)]
    check("all 9 boundary edges have length exactly 1 in both",
          all(abs(MR[v(x), v(y)] - 1) < 1e-15 and abs(MD[v(x), v(y)] - 1) < 1e-15
              for x, y in edges))
    check("apex chord d^R(P,Q) = 2*sqrt(6)/3", abs(MR[v(P), v(Q)] - SQ83) < 1e-15,
          f"{MR[v(P), v(Q)]:.15f} vs {SQ83:.15f}")
    check("apex chord d^D(P,Q) = 1 (the new edge)", abs(MD[v(P), v(Q)] - 1) < 1e-15)

    for nm, M in (("R", MR), ("D", MD)):
        tri = (M[:, :, None] + M[None, :, :] - M[:, None, :]).min()
        check(f"triangle inequality in {nm}", tri > -1e-12, f"min slack {tri:+.2e}")

    # EQUIVARIANCE under the configuration's own symmetry group (S3 on the
    # face x S2 on the apexes, order 12). The ball is the whole world here, so
    # this group is an exact symmetry of the object and the distance map must
    # commute with it. Independent of every other check: a wrong strip
    # placement, chord test or node index would break it.
    # one permutation table serves both configs: build() asserts the two
    # boundary grids coincide (mR.keys == mD.keys), which is what lets d^D-d^R
    # be taken entrywise at all.
    perms = [mR.node_permutation(g) for g in CONFIG_SYMMETRY]
    for nm, M in (("R", MR), ("D", MD)):
        worst = 0.0
        for perm in perms:
            worst = max(worst, float(np.abs(M - M[np.ix_(perm, perm)]).max()))
        check(f"{nm}: equivariant under all {len(CONFIG_SYMMETRY)} "
              f"S3xS2 relabellings", worst == 0.0,
              f"max |d(x,y) - d(gx,gy)| = {worst:.2e} (must be EXACTLY 0)")

    S, _ = boundary_surface_lengths(order)
    for nm, M in (("R", MR), ("D", MD)):
        check(f"{nm}: solid never longer than its own boundary surface",
              (M <= S + 1e-12).all(),
              f"max mean shortcut {np.mean(S - M):.4f}")

    # independent bound: a Steiner graph over B's INTERIOR can represent
    # paths that touch PQ or hug dB, which the strip enumeration cannot.
    for nm, M in (("R", MR), ("D", MD)):
        gaps = []
        for g in (order, 2 * order):
            T, _ = steiner_interior(nm, order, grid=g)
            gaps.append((g, (M - T).max(), np.mean(T - M)))
        check(f"{nm}: interior Steiner graph never beats the strip minimum",
              all(mx < 1e-12 for _, mx, _ in gaps),
              "  ".join(f"grid {g}: max(exact-steiner) {mx:+.2e}, "
                        f"mean excess {ex:.5f}" for g, mx, ex in gaps))
    return ok


def report(mR, MR, MD):
    dd = MD - MR
    iu = np.triu_indices(len(MR), 1)
    d = dd[iu]
    print(f"\n  pairs: {len(d)}   Delta = d_B^D - d_B^R (unit-edge units)")
    print(f"    shortened  (Delta < -1e-12): {(d < -1e-12).mean()*100:6.2f}%")
    print(f"    unchanged                  : {(np.abs(d) <= 1e-12).mean()*100:6.2f}%")
    print(f"    LENGTHENED (Delta > +1e-12): {(d > 1e-12).mean()*100:6.2f}%")
    print(f"    min {d.min():+.6f}   mean {d.mean():+.6f}   max {d.max():+.6f}")
    print(f"    max shortcut vs 2*sqrt(6)/3 - 1 = {SQ83 - 1:.6f}: "
          f"{abs(d.min() + (SQ83 - 1)) < 1e-12}")

    keys = mR.keys
    for label, sel in (("shortened", np.argsort(d)[:4]),
                       ("LENGTHENED", np.argsort(d)[-4:][::-1])):
        print(f"\n  extreme {label} pairs:")
        for t in sel:
            i, j = iu[0][t], iu[1][t]
            wx = dict(keys[i]); wy = dict(keys[j])
            print(f"    {node_name(wx):>18} -> {node_name(wy):<18} "
                  f"d_R {MR[i,j]:.6f}  d_D {MD[i,j]:.6f}  Delta {dd[i,j]:+.6f}")

    # 6x6 face-block summary
    fidx = {t: [i for i, k in enumerate(keys) if set(dict(k)) <= set(t)]
            for t in BOUNDARY}
    print("\n  mean Delta by boundary-face pair (rows/cols = the six triangles):")
    hdr = "        " + "".join(f"{face_name(t):>9}" for t in BOUNDARY)
    print(hdr)
    blocks = np.zeros((6, 6))
    for a_, ta in enumerate(BOUNDARY):
        row = f"    {face_name(ta):<4}"
        for b_, tb in enumerate(BOUNDARY):
            sub = dd[np.ix_(fidx[ta], fidx[tb])]
            blocks[a_, b_] = sub.mean()
            row += f"{sub.mean():>9.4f}"
        print(row)
    return dd, blocks


def figure(mR, MR, MD, dd, blocks, order, outdir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Polygon

    keys = mR.keys
    iu = np.triu_indices(len(MR), 1)
    fig = plt.figure(figsize=(15, 9))
    gs = fig.add_gridspec(2, 3, height_ratios=[1, 1.15], hspace=0.32, wspace=0.28)

    ax = fig.add_subplot(gs[0, 0])
    ax.plot([0, MR.max()], [0, MR.max()], "k--", lw=1, zorder=3, label="no change")
    sgn = dd[iu]
    for m, col, lab in ((sgn < -1e-12, "#2471a3", "through-paths: shorter"),
                        (np.abs(sgn) <= 1e-12, "0.7", "unchanged"),
                        (sgn > 1e-12, "#b03a2e", "tangential: longer")):
        ax.scatter(MR[iu][m], MD[iu][m], s=7, alpha=0.5, color=col,
                   edgecolors="none", label=lab)
    ax.set_xlabel(r"$d_B^{R}$  (2-tet ball)")
    ax.set_ylabel(r"$d_B^{D}$  (3-tet ball)")
    ax.set_title("the ball acts as a lens, not a shortcut")
    ax.legend(fontsize=8, loc="upper left")

    ax = fig.add_subplot(gs[0, 1])
    ax.hist(dd[iu], bins=60, color="#b03a2e")
    ax.axvline(-(SQ83 - 1), color="k", ls="--", lw=1)
    ax.text(-(SQ83 - 1), ax.get_ylim()[1] * 0.92,
            r"  $-(2\sqrt{6}/3-1)$", fontsize=9, va="top")
    ax.set_xlabel(r"$\Delta = d_B^{D} - d_B^{R}$")
    ax.set_ylabel("boundary pairs")
    ax.set_title(r"$\Delta$ distribution")

    ax = fig.add_subplot(gs[0, 2])
    bl = np.abs(blocks).max()
    im = ax.imshow(blocks, cmap="coolwarm", vmin=-bl, vmax=bl)
    ax.set_xticks(range(6)); ax.set_yticks(range(6))
    names = [face_name(t) for t in BOUNDARY]
    ax.set_xticklabels(names, fontsize=8); ax.set_yticklabels(names, fontsize=8)
    for i in range(6):
        for j in range(6):
            ax.text(j, i, f"{blocks[i,j]:.3f}", ha="center", va="center",
                    fontsize=7, color="k")
    ax.set_title(r"mean $\Delta$ by boundary-face pair")
    fig.colorbar(im, ax=ax, fraction=0.046)

    # --- field on the unfolded boundary, from two sources ---
    h = np.sqrt(3) / 2
    def layout(face_i):
        col, row = face_i % 3, face_i // 3
        return np.array([col * 1.35, -row * 1.15])
    def nodexy(w, tri, off):
        loc = {tri[0]: np.array([0.0, 0.0]), tri[1]: np.array([1.0, 0.0]),
               tri[2]: np.array([0.5, h])}
        return off + sum(float(wt) * loc[v] for v, wt in w.items())

    def most_central(tri):
        """Grid node of `tri` closest to its centroid (exists at any order)."""
        best, bk = None, None
        for i, kk in enumerate(mR.keys):
            w = dict(kk)
            if set(w) != set(tri):
                continue
            s = sum((float(w[v]) - 1 / 3) ** 2 for v in tri)
            if best is None or s < best:
                best, bk = s, i
        if bk is None:                       # order < 3: no interior nodes
            bk = next(i for i, kk in enumerate(mR.keys)
                      if set(dict(kk)) <= set(tri))
        return bk

    for k, (si, label) in enumerate([
            (mR.keys.index(node_key({P: Fraction(1)})), "source = apex P"),
            (most_central((A, B, Q)), "source = interior of face ABQ")]):
        ax = fig.add_subplot(gs[1, k])
        field = dd[si]
        xs, ys, cs = [], [], []
        for fi, tri in enumerate(BOUNDARY):
            off = layout(fi)
            corners = np.array([nodexy({tri[m]: Fraction(1)}, tri, off)
                                for m in range(3)])
            ax.add_patch(Polygon(corners, closed=True, fc="none", ec="0.6", lw=0.8))
            for m in range(3):
                ax.text(*corners[m], VERTEX_NAMES[tri[m]], fontsize=7,
                        color="0.35", ha="center", va="center")
            for i, kk in enumerate(mR.keys):
                w = dict(kk)
                if set(w) <= set(tri):
                    p = nodexy(w, tri, off)
                    xs.append(p[0]); ys.append(p[1]); cs.append(field[i])
        lim = np.abs(dd).max()
        sc = ax.scatter(xs, ys, c=cs, s=34, cmap="coolwarm",
                        vmin=-lim, vmax=lim)
        ax.set_aspect("equal"); ax.axis("off")
        ax.set_title(f"{label}\n" + r"$\Delta$ over the six boundary triangles",
                     fontsize=10)
        fig.colorbar(sc, ax=ax, fraction=0.04)

    ax = fig.add_subplot(gs[1, 2]); ax.axis("off")
    ax.text(0, 1.0, "Elementary 2$\\to$3 defect ball\nexact boundary distance map",
            fontsize=13, va="top", weight="bold")
    ax.text(0, 0.78,
            "R: tets ABCP, ABCQ        (internal face ABC)\n"
            "D: tets ABPQ, ACPQ, BCPQ  (internal edge PQ, degree 3)\n"
            "$\\partial B$: ABP ACP BCP ABQ ACQ BCQ  (identical in both)\n\n"
            f"Steiner order n = {order}, {len(mR.nodes)} boundary nodes,\n"
            f"{len(iu[0])} pairs. Unit-edge units.\n\n"
            "EXACT: rigid development in the edge-$\\sqrt{2}$ rational\n"
            "embedding; every chord certified contained in its tet\n"
            "strip by exact rational sign tests. Not sampled data --\n"
            "no action, couplings, or chains. UNIVERSAL: applies to\n"
            "every 2$\\to$3 in every triangulation.\n\n"
            f"max shortcut  $2\\sqrt{{6}}/3 - 1$ = {-(SQ83-1):+.6f}\n"
            f"max lengthening          = {dd[iu].max():+.6f}\n"
            f"mean $\\Delta$ = {dd[iu].mean():+.6f}\n"
            f"shortened {(dd[iu] < -1e-12).mean()*100:.1f}%   "
            f"lengthened {(dd[iu] > 1e-12).mean()*100:.1f}%\n\n"
            "Blue = shorter through the defect, red = longer.\n"
            "Global distances take the min with routes AROUND B,\n"
            "so lengthening here need not lengthen d(x,y).",
            fontsize=8.5, va="top", family="monospace")

    os.makedirs(outdir, exist_ok=True)
    path = os.path.join(outdir, f"defect_boundary_map_n{order}.png")
    fig.savefig(path, dpi=140, bbox_inches="tight")
    return path


def dump(mR, order, outdir):
    """Exact table: squared distances as Fraction strings."""
    recs = {"order": order, "units": "unit edge length",
            "nodes": [node_name(dict(k)) for k in mR.keys],
            "faces": [[face_name(t) for t in node_faces(dict(k))]
                      for k in mR.keys]}
    for cfg in ("R", "D"):
        m = BallBoundaryMap(cfg, order)
        vals = {}
        for i in range(len(m.nodes)):
            for j in range(i + 1, len(m.nodes)):
                L = m.dist_sq(m.nodes[i], m.nodes[j])
                vals[f"{i},{j}"] = None if L is None else str(L)
        recs[f"dist_sq_{cfg}"] = vals
    os.makedirs(outdir, exist_ok=True)
    path = os.path.join(outdir, f"defect_boundary_map_n{order}.json")
    with open(path, "w") as fh:
        json.dump(recs, fh)
    return path


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--order", type=int, default=6,
                    help="Steiner subdivision order on the boundary triangles")
    ap.add_argument("--out", default=os.path.join(_ROOT, "data", "figs"))
    ap.add_argument("--no-dump", action="store_true")
    args = ap.parse_args()

    print(f"exact boundary distance map of the 2->3 ball, order n={args.order}")
    mR, mD, MR, MD, missR, missD = build(args.order)
    print(f"  boundary grid: {len(mR.nodes)} nodes "
          f"({len(mR.strips)} strips in R, {len(mD.strips)} in D)\n")
    ok = validate(mR, MR, MD, missR, missD, args.order)
    dd, blocks = report(mR, MR, MD)
    p = figure(mR, MR, MD, dd, blocks, args.order, args.out)
    print(f"\n  figure: {p}")
    if not args.no_dump:
        print(f"  exact table: {dump(mR, args.order, args.out)}")
    print(f"\n  validation: {'ALL PASS' if ok else 'FAILURES ABOVE'}")
    sys.exit(0 if ok else 1)
