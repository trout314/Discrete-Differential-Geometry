#!/usr/bin/env python3
"""Graph-distance distributions for seed triangulations.

For one or more seed ``.mfd`` files, compute the distribution of shortest-path
distances d(x, y) in two graphs:

  * ``edge``  -- the 1-skeleton (vertices, joined by edges of the triangulation).
  * ``dual``  -- the dual graph (one node per facet; facets joined when they
                 share a ridge / codimension-1 face).
  * ``steiner`` -- the TRUE PL geodesic vertex-to-vertex distance via Steiner-point
                 Dijkstra (steiner_geodesic.py): a rigorous upper bound converging
                 down to the exact flat-cone-metric distance as the subdivision
                 order rises. Unlike ``edge`` this removes the direction-dependent
                 lattice distortion, so it reports the metric SHAPE (not just the
                 scaling). Continuous-valued, so it is binned rather than shelled.

Distances for ``edge``/``dual`` come from unweighted BFS. Computing *all* pairs is
O(V*E); instead we sample a set of source nodes uniformly and histogram the
distances from each to every other node. By homogeneity the single-source distance
distribution equals the pairwise one, so sampling sources is unbiased (and we
average over the provided seeds/replicas to reduce noise). Pass ``--sources all``
for exact all-pairs on small graphs. ``steiner`` samples sources the same way but
builds the Steiner graph once per seed, then runs weighted Dijkstra.

Output: a tidy CSV (graph, distance, count, pdf, cdf) and a printed summary.
Optionally a PNG overlay of the two PDFs with ``--plot``.

Examples
--------
    python scripts/distance_distribution.py seeds/S3_N1e4_1e-1_ED5p0043_1e-1_s000.mfd
    python scripts/distance_distribution.py 'seeds/S3_N1e4_1e-1_ED5p0043_1e-1_s0*.mfd' \
        --sources 800 --plot out/dist_N1e4.png --csv out/dist_N1e4.csv
"""

import argparse
import glob
import os
import sys
import tempfile

import numpy as np

try:
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import shortest_path
except ImportError:
    sys.exit("This script needs scipy (pip install scipy).")

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "python"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # sibling: steiner_geodesic
from discrete_differential_geometry import Manifold


def read_edge_list(path):
    """Read a whitespace 'a b' edge-list file -> (E,2) int array of raw labels."""
    if os.path.getsize(path) == 0:
        return np.empty((0, 2), dtype=np.int64)
    return np.loadtxt(path, dtype=np.int64, ndmin=2)


def build_csr(edges, n_hint=None):
    """Symmetric unweighted CSR from a raw edge list, remapping node labels to
    a contiguous 0..N-1 range. Returns (csr, n_nodes, label_of_index)."""
    if edges.size == 0:
        raise ValueError("empty edge list")
    labels, inv = np.unique(edges, return_inverse=True)
    inv = inv.reshape(edges.shape)
    n = len(labels)
    if n_hint is not None and n != n_hint:
        print(f"  note: {n} distinct nodes in graph vs f-vector hint {n_hint}",
              file=sys.stderr)
    data = np.ones(len(inv), dtype=np.int8)
    # Symmetrize: add both directions so BFS is undirected regardless of order.
    rows = np.concatenate([inv[:, 0], inv[:, 1]])
    cols = np.concatenate([inv[:, 1], inv[:, 0]])
    g = coo_matrix((np.concatenate([data, data]), (rows, cols)), shape=(n, n))
    return g.tocsr(), n, labels


def distance_histogram(csr, n, n_sources, rng, chunk=32):
    """Accumulate a shell histogram of BFS distances from sampled sources.

    Returns (hist, n_sources_used) where hist[r] = number of (source, target)
    ordered pairs at distance r (r >= 1; self-distance 0 excluded)."""
    n_sources = n if n_sources is None or n_sources >= n else n_sources
    sources = (np.arange(n) if n_sources == n
               else rng.choice(n, size=n_sources, replace=False))
    hist = np.zeros(1, dtype=np.int64)
    for i in range(0, len(sources), chunk):
        batch = sources[i:i + chunk]
        dmat = shortest_path(csr, method="D", unweighted=True,
                             directed=False, indices=batch)
        d = dmat[np.isfinite(dmat)].astype(np.int64)
        d = d[d > 0]  # drop self-distances
        if d.size:
            counts = np.bincount(d)
            if len(counts) > len(hist):
                hist = np.pad(hist, (0, len(counts) - len(hist)))
            hist[:len(counts)] += counts
    return hist, len(sources)


def steiner_distance_histogram(mfd, order, width, n_sources, rng, chunk=48):
    """Continuous PL-geodesic (Steiner-Dijkstra) vertex-to-vertex distance histogram.

    Builds the Steiner graph once (a rigorous upper bound at the given subdivision
    order), then runs multi-source weighted Dijkstra and bins vertex-to-vertex
    distances at ``width`` (edge-length units). Returns (hist, n_vertices,
    n_sources_used) with hist[r] = # pairs with r*width <= d < (r+1)*width."""
    from steiner_geodesic import build_steiner_graph
    F = np.asarray(mfd.facets(), np.int64)
    labels, inv = np.unique(F, return_inverse=True)
    T = inv.reshape(F.shape)
    V = len(labels)
    G, _ = build_steiner_graph(T, V, order)
    ns = V if n_sources is None or n_sources >= V else n_sources
    src = np.arange(V) if ns == V else rng.choice(V, size=ns, replace=False)
    hist = np.zeros(1, dtype=np.int64)
    for i in range(0, len(src), chunk):
        batch = src[i:i + chunk]
        dmat = shortest_path(G, method="D", directed=False, indices=batch)
        d = dmat[:, :V]                       # vertex-to-vertex only (drop Steiner nodes)
        d = d[np.isfinite(d)]
        d = d[d > 1e-9]                       # drop self-distances
        if d.size:
            counts = np.bincount(np.floor(d / width).astype(np.int64))
            if len(counts) > len(hist):
                hist = np.pad(hist, (0, len(counts) - len(hist)))
            hist[:len(counts)] += counts
    return hist, V, len(src)


def summarize(hist, width=1.0):
    """Mean, median, and max (diameter proxy) from a shell/bin histogram.

    ``width`` is the physical size of one bin (1.0 for the integer edge/dual
    shells; the Steiner bin width for the continuous metric)."""
    r = np.arange(len(hist))
    total = hist.sum()
    if total == 0:
        return dict(mean=0.0, median=0.0, diameter=0.0, pairs=0)
    mean = (r * hist).sum() / total * width
    cdf = np.cumsum(hist) / total
    median = float(np.searchsorted(cdf, 0.5)) * width
    diameter = float(r[hist > 0][-1]) * width
    return dict(mean=float(mean), median=median, diameter=diameter, pairs=int(total))


def analyze(seed_paths, dim, n_sources, rng, metrics=("edge", "dual"),
            steiner_order=3, steiner_bin=0.25):
    """Aggregate distance histograms over the given seeds for the chosen metrics.

    Returns (agg, width, n_nodes, n_srcs): dicts keyed by metric. ``width`` is the
    per-metric bin size (1.0 for edge/dual integer shells, ``steiner_bin`` for the
    continuous Steiner metric)."""
    agg = {m: np.zeros(1, dtype=np.int64) for m in metrics}
    width = {m: (steiner_bin if m == "steiner" else 1.0) for m in metrics}
    n_nodes = {m: 0 for m in metrics}
    n_srcs = {m: 0 for m in metrics}
    with tempfile.TemporaryDirectory() as tmp:
        for path in seed_paths:
            mfd = Manifold.load(path, dim)
            fv = list(mfd.f_vector)
            for kind in metrics:
                if kind == "steiner":
                    hist, n, used = steiner_distance_histogram(
                        mfd, steiner_order, steiner_bin, n_sources, rng)
                else:
                    saver, n_hint = {
                        "edge": (mfd.save_edge_graph, fv[0]),
                        "dual": (mfd.save_dual_graph, fv[dim]),
                    }[kind]
                    gfile = os.path.join(tmp, f"{kind}.txt")
                    saver(gfile)
                    csr, n, _ = build_csr(read_edge_list(gfile), n_hint)
                    hist, used = distance_histogram(csr, n, n_sources, rng)
                if len(hist) > len(agg[kind]):
                    agg[kind] = np.pad(agg[kind], (0, len(hist) - len(agg[kind])))
                agg[kind][:len(hist)] += hist
                n_nodes[kind] = n          # nodes per seed (seeds share params)
                n_srcs[kind] += used
    return agg, width, n_nodes, n_srcs


def write_csv(path, agg, width):
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w") as f:
        f.write("graph,distance,count,pdf,cdf\n")
        for kind, hist in agg.items():
            w = width[kind]
            total = hist.sum() or 1
            cum = 0
            for r, c in enumerate(hist):
                cum += c
                f.write(f"{kind},{r*w:.6g},{int(c)},{c/total:.8g},{cum/total:.8g}\n")


def maybe_plot(path, agg, width):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("  (matplotlib not available; skipping plot)", file=sys.stderr)
        return
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    fig, ax = plt.subplots(figsize=(7, 4.5))
    for kind, hist in agg.items():
        w = width[kind]
        total = hist.sum() or 1
        r = np.arange(len(hist)) * w
        dens = hist / (total * w)          # probability density (per unit distance)
        style = dict(lw=1.8) if kind == "steiner" else dict(marker=".", ms=4, lw=1)
        ax.plot(r, dens, label=f"{kind}-distance", **style)
    ax.set_xlabel("distance  (edge-length units)"); ax.set_ylabel("density  P(d)/Δ")
    ax.set_title("Distance distribution"); ax.legend()
    fig.tight_layout(); fig.savefig(path, dpi=130)
    print(f"  wrote plot {path}")


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("seeds", nargs="+", help="Seed .mfd file(s) or glob(s).")
    p.add_argument("--dim", type=int, default=3)
    p.add_argument("--sources", default="1000",
                   help="Number of source nodes, or 'all' for exact all-pairs "
                        "(default 1000).")
    p.add_argument("--metric", default="edge,dual",
                   help="Comma list of metrics: edge, dual, steiner (default edge,dual). "
                        "'steiner' is the true PL geodesic (Steiner-Dijkstra).")
    p.add_argument("--steiner-order", type=int, default=3,
                   help="Steiner subdivision order n (1 recovers edge distance; "
                        "n=3 essentially converges the bound; default 3).")
    p.add_argument("--steiner-bin", type=float, default=0.25,
                   help="Bin width for the continuous Steiner histogram, in "
                        "edge-length units (default 0.25).")
    p.add_argument("--rng-seed", type=int, default=0, help="RNG seed for source sampling.")
    p.add_argument("--csv", default=None, help="Write the distribution to this CSV.")
    p.add_argument("--plot", default=None, help="Write a PNG overlay of the metric PDFs.")
    args = p.parse_args()

    paths = sorted({q for pat in args.seeds for q in glob.glob(pat)})
    if not paths:
        sys.exit(f"No files matched: {args.seeds}")
    metrics = tuple(m.strip() for m in args.metric.split(",") if m.strip())
    bad = set(metrics) - {"edge", "dual", "steiner"}
    if bad:
        sys.exit(f"unknown metric(s): {sorted(bad)} (choose from edge, dual, steiner)")
    n_sources = None if args.sources == "all" else int(args.sources)
    rng = np.random.default_rng(args.rng_seed)

    extra = f", steiner n={args.steiner_order} bin={args.steiner_bin}" if "steiner" in metrics else ""
    print(f"Analyzing {len(paths)} seed(s), dim={args.dim}, metrics={','.join(metrics)}, "
          f"sources={'all' if n_sources is None else n_sources}{extra}", flush=True)
    agg, width, n_nodes, n_srcs = analyze(paths, args.dim, n_sources, rng, metrics,
                                          args.steiner_order, args.steiner_bin)

    print(f"\n{'metric':>8} {'nodes':>9} {'sources':>8} {'pairs':>12} "
          f"{'mean':>8} {'median':>8} {'diameter':>9}")
    for kind in metrics:
        s = summarize(agg[kind], width[kind])
        print(f"{kind:>8} {n_nodes[kind]:>9} {n_srcs[kind]:>8} {s['pairs']:>12} "
              f"{s['mean']:>8.3f} {s['median']:>8.2f} {s['diameter']:>9.2f}")

    if args.csv:
        write_csv(args.csv, agg, width); print(f"\nwrote {args.csv}")
    if args.plot:
        maybe_plot(args.plot, agg, width)


if __name__ == "__main__":
    main()
