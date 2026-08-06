#!/usr/bin/env python3
"""Defect catalog: one HTML page per snapshot with a sortable table of every
defect complex's properties, each row linked to its interactive 3D model
(defect_viewer scenes, rendered into the same directory).

Per complex: size, illegal-edge composition (deg-4 / deg-3 / deg-7+),
anomalous-n6 nodes, induced f-vector, total coordination, illegal-graph
anatomy (fragments, tips, branch points, independent cycles = rings), net
curvature charge Q vs crystal background, gyration radius in the harmonic
chart, closed-star size, and the decorated-isomorphism class (short hash of
canonical_key; equal hash = same decorated species).

Output directory (default data/viz/<snap-base>/):
  index.html + cx<i>.html scenes + plotly.min.js (shared, so scenes stay
  small).  Open index.html locally; links are relative.

FRAME: harmonic lift chart (CONVENTIONS sec 6).
"""
import argparse
import hashlib
import html as html_mod
import os
import sys
from collections import Counter, defaultdict
from itertools import combinations

import numpy as np

_ROOT = os.path.normpath(os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "..", ".."))
sys.path.insert(0, os.path.join(_ROOT, "python"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import discrete_differential_geometry as ddg
import defect_state as dsm
import defect_viewer as dv

CELL = 1e6


def illegal_graph_anatomy(ill_edges):
    """(#fragments, #tips, #branches, #cycles, largest fragment edge count)
    of the illegal-edge graph."""
    adj = defaultdict(set)
    for a, b in ill_edges:
        adj[a].add(b)
        adj[b].add(a)
    seen, frags = set(), []
    for v0 in adj:
        if v0 in seen:
            continue
        stack, comp = [v0], set()
        seen.add(v0)
        while stack:
            u = stack.pop()
            comp.add(u)
            for w in adj[u]:
                if w not in seen:
                    seen.add(w)
                    stack.append(w)
        ne = sum(1 for a, b in ill_edges if a in comp and b in comp)
        frags.append((len(comp), ne))
    tips = sum(1 for v, nb in adj.items() if len(nb) == 1)
    branches = sum(1 for v, nb in adj.items() if len(nb) >= 3)
    nv = len(adj)
    ne = len(ill_edges)
    cycles = ne - nv + len(frags)          # first Betti number
    biggest = max((ne for _, ne in frags), default=0)
    return len(frags), tips, branches, cycles, biggest


TABLE_JS = """
<script>
function sortTable(k) {
  const t = document.getElementById('cat');
  const rows = Array.from(t.rows).slice(1);
  const dir = t.dataset['d' + k] === 'a' ? -1 : 1;
  t.dataset['d' + k] = dir === 1 ? 'a' : 'b';
  rows.sort((r1, r2) => {
    const a = r1.cells[k].dataset.v ?? r1.cells[k].innerText;
    const b = r2.cells[k].dataset.v ?? r2.cells[k].innerText;
    const na = parseFloat(a), nb = parseFloat(b);
    if (!isNaN(na) && !isNaN(nb)) return dir * (na - nb);
    return dir * a.localeCompare(b);
  });
  rows.forEach(r => t.appendChild(r));
}
</script>
"""

CSS = """
<style>
body { font-family: -apple-system, Helvetica, sans-serif; margin: 24px; }
table { border-collapse: collapse; font-size: 13px; }
th, td { border: 1px solid #ccc; padding: 4px 8px; text-align: right; }
th { background: #f0f0f0; cursor: pointer; position: sticky; top: 0; }
td.l, th.l { text-align: left; }
tr:hover { background: #fff7e6; }
a { color: #d97706; font-weight: 600; text-decoration: none; }
.meta { color: #555; margin-bottom: 14px; }
</style>
"""


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("snap")
    ap.add_argument("--out-dir", default=None,
                    help="default data/viz/<snap-base>/")
    ap.add_argument("--note", default="",
                    help="provenance line for the page header (action, "
                         "couplings, host -- CONVENTIONS sec 2)")
    ap.add_argument("--max-canon", type=int, default=60,
                    help="skip canonical-class hashing above this size")
    ap.add_argument("--layout", choices=("harmonic", "mds"),
                    default="harmonic",
                    help="vertex positions for the 3D scenes: 'harmonic' "
                         "needs a .cocycle.npz sidecar; 'mds' computes a "
                         "per-complex classical-MDS layout from graph "
                         "distances of the closed-star region -- use for "
                         "states without a maintained cocycle (e.g. "
                         "contract/split-channel runs, where the lift "
                         "cannot follow vertex creation/deletion). "
                         "COMBINATORIAL layout: local shape faithful, no "
                         "global chart, gyration radii in graph units.")
    ap.add_argument("--narrow", action="store_true",
                    help="use the historical complex definition (illegal-"
                         "edge incidence only, components(broad=False)). "
                         "On hosts with a dense native disclination network "
                         "(a15, sigma) the broad definition percolates "
                         "through n6-anomalous vertices and fuses every "
                         "complex into one.")
    args = ap.parse_args()

    base = os.path.basename(args.snap)[:-4]
    out_dir = args.out_dir or os.path.join(_ROOT, "data", "viz", base)
    os.makedirs(out_dir, exist_ok=True)

    m = ddg.Manifold.load(args.snap, 3)
    fac = np.asarray(m.facets())
    st = dsm.DefectState(m)
    if args.narrow:
        # Components of the ILLEGAL-EDGE GRAPH (edges linked by shared
        # vertices) -- the gas-scan census definition. At high defect-vertex
        # density (a15-like hosts) any vertex-adjacency definition
        # percolates; edge-graph clustering is what "dilute" refers to.
        e2i = {}
        parent = list(range(len(st.ill_edges)))

        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x

        ill = sorted(st.ill_edges)
        for i, e in enumerate(ill):
            for v in e:
                if v in e2i:
                    r1, r2 = find(i), find(e2i[v])
                    if r1 != r2:
                        parent[r1] = r2
                else:
                    e2i[v] = i
        groups = defaultdict(list)
        for i, e in enumerate(ill):
            groups[find(i)].append(e)
        comps = []
        for edges in groups.values():
            cv = sorted({v for e in edges for v in e})
            cvs = set(cv)
            sig = sorted(st.edeg[e] for e in edges)
            nodes = sorted(st.n6[v] for v in cv
                           if st.imp[v] == 0 and st.n6[v] not in dsm.FK_N6)
            comps.append(dsm.Complex(cv, tuple(sig), tuple(nodes)))
        comps.sort(key=lambda c: -len(c.verts))
    else:
        comps = sorted(st.components(), key=lambda c: -len(c.verts))
    if args.layout == "harmonic":
        X, P = dv.positions(args.snap, fac)
    else:
        # per-complex MDS layout, filled lazily before each render; the huge
        # period makes the viewer's min-imaging inert
        X = np.zeros((int(fac.max()) + 1, 3))
        P = np.array([1e9, 1e9, 1e9])

    def fill_mds(region_verts, edges_in_region):
        """Classical MDS of the region's graph metric into X (in place)."""
        verts = sorted(region_verts)
        idx = {v: i for i, v in enumerate(verts)}
        n = len(verts)
        adj = [[] for _ in range(n)]
        for a, b in edges_in_region:
            adj[idx[a]].append(idx[b])
            adj[idx[b]].append(idx[a])
        D = np.full((n, n), np.inf)
        for s in range(n):
            D[s, s] = 0
            frontier = [s]
            d = 0
            while frontier:
                d += 1
                nxt = []
                for u in frontier:
                    for w in adj[u]:
                        if D[s, w] == np.inf:
                            D[s, w] = d
                            nxt.append(w)
                frontier = nxt
        D[~np.isfinite(D)] = D[np.isfinite(D)].max() + 1
        D2 = D ** 2
        J = np.eye(n) - np.ones((n, n)) / n
        B = -0.5 * J @ D2 @ J
        w_, V_ = np.linalg.eigh(B)
        order = np.argsort(w_)[::-1][:3]
        coords = V_[:, order] * np.sqrt(np.maximum(w_[order], 1e-12))
        for v, i in idx.items():
            X[v] = coords[i]

    q = st.vertex_charges()
    qb = st.qbar(q)
    n_ill_tot = len(st.ill_edges)
    print(f"{base}: {len(comps)} complexes, {n_ill_tot} illegal edges")

    rows = []
    for i, cx in enumerate(comps):
        V = set(cx.verts)
        scene = f"cx{i}.html"
        chart = ("harmonic lift chart" if args.layout == "harmonic"
                 else "combinatorial MDS layout")
        title = (f"{base} complex {i}: n={len(cx.verts)} sig={cx.sig} "
                 f"nodes={cx.nodes} [{chart}]")
        if args.layout == "mds":
            star_tets = st.star(V)
            region = {v for t in star_tets for v in t}
            eir = set()
            for t in star_tets:
                for a, b in combinations(sorted(t), 2):
                    eir.add((a, b))
            fill_mds(region, eir)
        dv.render(st, cx, X, P, os.path.join(out_dir, scene), title,
                  plotlyjs="directory")

        ill_in = {e: st.edeg[e] for e in st.ill_edges
                  if e[0] in V and e[1] in V}
        degs = Counter(ill_in.values())
        nfrag, tips, branches, cycles, bigfrag = \
            illegal_graph_anatomy(list(ill_in))
        sh = st.induced_shape(cx.verts)
        f = tuple(int(x) for x in sh["f"])
        Z = int(st.total_coordination(cx.verts))
        Q = float(st.complex_charge_rel(cx.verts, q, qb))
        star_tets = st.star(V)
        # gyration radius in the harmonic chart (BFS-unwrapped), in cells
        star_edges = {e for t in star_tets
                      for e in combinations(sorted(t), 2)}
        adj = defaultdict(set)
        for a, b in star_edges:
            adj[a].add(b)
            adj[b].add(a)
        pos = dv.unwrap(X, P, V, cx.verts[0], adj)
        A = np.array([pos[v] for v in cx.verts])
        rg = float(np.sqrt(((A - A.mean(0)) ** 2).sum(1).mean())) / CELL
        # decorated species hash
        if len(cx.verts) <= args.max_canon:
            fac_i = st.induced_facets(cx.verts)
            vc, ec = st.decorations(cx.verts)
            key, exact = dsm.canonical_key(fac_i, vc, ec)
            species = hashlib.sha1(repr(key).encode()).hexdigest()[:8]
            if not exact:
                species += "~"
        else:
            species = "-"
        rows.append(dict(
            i=i, link=scene, n=len(cx.verts), nill=len(ill_in),
            d4=degs.get(4, 0), d3=degs.get(3, 0),
            d7=sum(v for k, v in degs.items() if k >= 7),
            nodes=" ".join(f"+{x - 4}" for x in cx.nodes) or "-",
            f=f, Z=Z, nfrag=nfrag, tips=tips, branches=branches,
            cycles=cycles, bigfrag=bigfrag, Q=Q, rg=rg,
            star=len(star_tets), species=species))
        print(f"  [{i}] n={len(cx.verts)} ill={len(ill_in)} -> {scene}")

    # ---- index page --------------------------------------------------
    hdr = ["#", "3D", "n vert", "n ill", "deg-4", "deg-3", "deg-7+",
           "n6 nodes", "f-vector", "Z", "frags", "tips", "branch",
           "rings", "max frag", "Q rel", "Rg (cells)", "star tets",
           "species"]
    th = "".join(f'<th class="{"l" if j in (1, 7, 8, 18) else ""}" '
                 f'onclick="sortTable({j})">{h}</th>'
                 for j, h in enumerate(hdr))
    body = []
    for r in rows:
        cells = [
            f'<td data-v="{r["i"]}">{r["i"]}</td>',
            f'<td class="l"><a href="{r["link"]}">view</a></td>',
            f'<td>{r["n"]}</td>', f'<td>{r["nill"]}</td>',
            f'<td>{r["d4"]}</td>', f'<td>{r["d3"]}</td>',
            f'<td>{r["d7"]}</td>',
            f'<td class="l">{html_mod.escape(r["nodes"])}</td>',
            f'<td class="l" data-v="{r["f"][0]}">{r["f"]}</td>',
            f'<td>{r["Z"]}</td>', f'<td>{r["nfrag"]}</td>',
            f'<td>{r["tips"]}</td>', f'<td>{r["branches"]}</td>',
            f'<td>{r["cycles"]}</td>', f'<td>{r["bigfrag"]}</td>',
            f'<td data-v="{r["Q"]:.3f}">{r["Q"]:+.2f}</td>',
            f'<td data-v="{r["rg"]:.4f}">{r["rg"]:.3f}</td>',
            f'<td>{r["star"]}</td>',
            f'<td class="l"><code>{r["species"]}</code></td>',
        ]
        body.append("<tr>" + "".join(cells) + "</tr>")
    note = html_mod.escape(args.note) if args.note else ""
    page = f"""<!DOCTYPE html>
<html><head><meta charset="utf-8"><title>Defect catalog: {base}</title>
{CSS}</head><body>
<h2>Defect catalog &mdash; {base}</h2>
<div class="meta">N3 = {m.num_facets} &middot; {len(comps)} complexes &middot;
{n_ill_tot} illegal edges ({dict(Counter(st.edeg[e] for e in st.ill_edges))})
&middot; frame: harmonic lift chart &middot; complexes: broad definition
(illegal-edge incidence or non-FK coordination)<br>
{note}<br>
Click a column header to sort. "rings" = independent cycles of the
illegal-edge graph. "species" = decorated-isomorphism class hash
(~ suffix: refinement-truncated, strong invariant only).
Q rel = curvature charge minus n&times;crystal-median.</div>
<table id="cat"><tr>{th}</tr>
{os.linesep.join(body)}
</table>{TABLE_JS}</body></html>
"""
    out = os.path.join(out_dir, "index.html")
    with open(out, "w") as fh:
        fh.write(page)
    print(f"\nSaved to: {os.path.abspath(out)}")


if __name__ == "__main__":
    main()
