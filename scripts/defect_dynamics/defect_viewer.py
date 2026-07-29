#!/usr/bin/env python3
"""Interactive 3D viewer for crystal defect complexes.

For a chosen defect complex of a snapshot, writes a self-contained
zoomable/rotatable HTML scene (plotly) containing:

  * the defect complex itself, decorated:
      - illegal edges THICK: orange = degree 4, red = degree 3
        (degree >= 7, if ever present, dark maroon);
      - legal edges between defect vertices: medium gray;
      - defect vertices incident to an illegal edge: dark markers;
      - defect vertices NOT incident to any illegal edge (the anomalous-
        coordination nodes): PURPLE, annotated with n6 - 4 (their count of
        degree-6 edges in excess of 4);
  * context: the closed star of the defect's vertices -- the one-shell
    region of legal crystal the defect occupies -- as thin translucent
    gray edges and small pale vertices.

Positions are the harmonic torus coordinates (<snap>.cocycle.npz),
min-imaged about the complex so the scene never wraps the box.

Usage:
  python defect_viewer.py data/mgas/lam35r_snap15000.mfd            # largest
  python defect_viewer.py SNAP.mfd --complex 3                      # by index
  python defect_viewer.py SNAP.mfd --all                            # every one
Complex indices are by descending size (0 = largest); use --list to see them.
FRAME: harmonic lift chart (CONVENTIONS sec 6); geometry is smoothed near
defects by the harmonic embedding.
"""
import argparse
import os
import sys
from collections import defaultdict
from itertools import combinations

import numpy as np

_ROOT = os.path.normpath(os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "..", ".."))
sys.path.insert(0, os.path.join(_ROOT, "python"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import discrete_differential_geometry as ddg
from discrete_differential_geometry import cocycle as coc
import defect_state as dsm

EDGE_COLORS = {3: "#d62728", 4: "#ff7f0e"}          # red / orange
EDGE_COLOR_HIGH = "#7f1010"                          # degree >= 7 (dark maroon)


def positions(snap, fac):
    edges, omega, _ = coc.load_cocycle(snap[:-4] + ".cocycle.npz")
    frac, basis = coc.torus_positions(fac, edges, omega)
    X = frac * np.diag(basis)
    P = np.abs(np.diag(basis))
    return X, P


def unwrap(X, P, verts, ref, adj=None):
    """Unwrapped positions of `verts`: BFS from `ref`, each vertex min-imaged
    RELATIVE TO an already-placed neighbour (robust for regions extending
    beyond half the box, where a single-reference min-image would fold)."""
    if adj is None:                     # fall back: min-image about ref
        d = X[sorted(verts)] - X[ref]
        d -= np.round(d / P) * P
        return {v: X[ref] + d[i] for i, v in enumerate(sorted(verts))}
    pos = {ref: X[ref].astype(float)}
    frontier = [ref]
    while frontier:
        nxt = []
        for u in frontier:
            for w in adj[u]:
                if w in pos or w not in verts:
                    continue
                d = X[w] - X[u]
                d -= np.round(d / P) * P
                pos[w] = pos[u] + d
                nxt.append(w)
        frontier = nxt
    for v in verts:                     # isolated fall-back (shouldn't happen)
        if v not in pos:
            d = X[v] - X[ref]
            d -= np.round(d / P) * P
            pos[v] = X[ref] + d
    return pos


def edge_trace(segs, pos, name, color, width, opacity=1.0, hover=None):
    import plotly.graph_objects as go
    xs, ys, zs = [], [], []
    for a, b in segs:
        pa, pb = pos[a], pos[b]
        xs += [pa[0], pb[0], None]
        ys += [pa[1], pb[1], None]
        zs += [pa[2], pb[2], None]
    return go.Scatter3d(
        x=xs, y=ys, z=zs, mode="lines", name=name,
        line=dict(color=color, width=width), opacity=opacity,
        hoverinfo="skip" if hover is None else "text", hovertext=hover)


def render(st, cx, X, P, out, title, plotlyjs=True):
    import plotly.graph_objects as go
    V = set(cx.verts)
    star_tets = st.star(V)
    region = {v for t in star_tets for v in t}

    # edge classification
    star_edges = set()
    for t in star_tets:
        for e in combinations(sorted(t), 2):
            star_edges.add(e)
    adj = defaultdict(set)
    for a, b in star_edges:
        adj[a].add(b)
        adj[b].add(a)
    pos = unwrap(X, P, region, cx.verts[0], adj)
    ill_in = {e: st.edeg[e] for e in st.ill_edges
              if e[0] in V and e[1] in V}
    induced_legal = {e for e in star_edges
                     if e[0] in V and e[1] in V and e not in ill_in}
    context = [e for e in star_edges
               if not (e[0] in V and e[1] in V)]

    illv = {v for e in ill_in for v in e}
    nodes = [v for v in cx.verts if v not in illv]     # purple, no illegal edge
    cores = sorted(illv)
    shell = sorted(region - V)

    def Z(v):
        return (len(st.v2t[v]) + 4) // 2

    fig = go.Figure()
    # context: closed star of the defect's vertices
    fig.add_trace(edge_trace(context, pos, "closed star (context)",
                             "#9a9a9a", 1.2, opacity=0.25))
    fig.add_trace(edge_trace(sorted(induced_legal), pos, "defect legal edges",
                             "#555555", 4, opacity=0.8))
    # illegal edges by degree, thick
    bydeg = defaultdict(list)
    for e, d in ill_in.items():
        bydeg[d].append(e)
    for d in sorted(bydeg):
        col = EDGE_COLORS.get(d, EDGE_COLOR_HIGH)
        fig.add_trace(edge_trace(
            bydeg[d], pos, f"degree-{d} edges (x{len(bydeg[d])})", col, 10))

    def scatter(vs, name, color, size, text=None, textpos="top center"):
        return go.Scatter3d(
            x=[pos[v][0] for v in vs], y=[pos[v][1] for v in vs],
            z=[pos[v][2] for v in vs],
            mode="markers+text" if text else "markers", name=name,
            marker=dict(color=color, size=size),
            text=text or None, textposition=textpos,
            textfont=dict(color="#7d26cd", size=14),
            hovertext=[f"v{v}  Z={Z(v)}  n6={st.n6[v]}  imp={st.imp[v]}"
                       for v in vs], hoverinfo="text")

    if shell:
        fig.add_trace(scatter(shell, "star shell", "#c9c9c9", 2.5))
    if cores:
        fig.add_trace(scatter(cores, "defect vertices", "#2b2b2b", 5))
    if nodes:
        fig.add_trace(scatter(
            nodes, "n6 nodes (excess over 4)", "#7d26cd", 8,
            text=[f"+{st.n6[v] - 4}" for v in nodes]))

    fig.update_layout(
        title=title, scene=dict(aspectmode="data",
                                xaxis_visible=False, yaxis_visible=False,
                                zaxis_visible=False),
        legend=dict(itemsizing="constant"),
        margin=dict(l=0, r=0, t=40, b=0))
    fig.write_html(out, include_plotlyjs=plotlyjs)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("snap")
    ap.add_argument("--complex", type=int, default=0,
                    help="complex index, by descending size (default 0)")
    ap.add_argument("--all", action="store_true", help="render every complex")
    ap.add_argument("--list", action="store_true", help="list complexes only")
    ap.add_argument("--out-dir", default=os.path.join(_ROOT, "data", "viz"))
    args = ap.parse_args()

    m = ddg.Manifold.load(args.snap, 3)
    fac = np.asarray(m.facets())
    st = dsm.DefectState(m)
    comps = sorted(st.components(), key=lambda c: -len(c.verts))
    print(f"{os.path.basename(args.snap)}: {len(comps)} complexes")
    if args.list:
        for i, c in enumerate(comps):
            print(f"  [{i}] n={len(c.verts)} sig={c.sig} nodes={c.nodes}")
        return
    X, P = positions(args.snap, fac)
    os.makedirs(args.out_dir, exist_ok=True)
    base = os.path.basename(args.snap)[:-4]
    todo = range(len(comps)) if args.all else [args.complex]
    for i in todo:
        cx = comps[i]
        out = os.path.join(args.out_dir, f"{base}_cx{i}.html")
        title = (f"{base} complex {i}: n={len(cx.verts)} "
                 f"sig={cx.sig} nodes={cx.nodes} "
                 f"[harmonic lift chart]")
        render(st, cx, X, P, out, title)
        print(f"  [{i}] n={len(cx.verts)} sig={cx.sig} -> "
              f"Saved to: {os.path.abspath(out)}")


if __name__ == "__main__":
    main()
