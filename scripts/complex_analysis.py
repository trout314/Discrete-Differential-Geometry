#!/usr/bin/env python3
"""Complex-resolved dopant interaction analysis in doped TCP crystals.

Upgrades dopant_pairs.py from vertex-level to COMPLEX-level statistics: the
physical particles are the bound defect complexes (joint clusters of all
non-host legal classes), not individual dopant vertices. Per state:

  * complex census: size, composition (per Z class), NET CURVATURE CHARGE
    Q = sum_{v in complex} qR(v) - size * qbar_host (the multipole test:
    Q ~ 0 => neutral dielectric gas, the hyperuniformity mechanism)
  * complex-complex pair correlation g_cc(r) (min BFS distance between
    vertex sets) vs a random-placement null -> u_eff(r) = -ln g_cc(r)
  * mass-action spectrum: n_k populations and binding free energies
    dF_k = -ln( n_k / n_1^k ) + (k-1) ln V   (ideal law of mass action)

Usage:
    python scripts/complex_analysis.py data/dope_hold/c15big_z14_mu3.0_final.mfd \
        --host-classes 0 4 --json out/complex_analysis_c15.json
"""
import argparse
import collections
import json
import os
import sys

import numpy as np

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_ROOT, "python"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import discrete_differential_geometry as ddg
from dopant_pairs import vertex_classes
from curvature_hyperuniformity_g import load_fields


def joint_complexes(n6, imp, adj, host_n6):
    """Connected clusters of legal-but-non-host-class vertices."""
    legal_nonhost = set(np.where((imp == 0) & ~np.isin(n6, host_n6)
                                 & np.isin(n6, [0, 2, 3, 4]))[0])
    seen, comps = set(), []
    for s in legal_nonhost:
        if s in seen:
            continue
        stack, comp = [s], []
        seen.add(s)
        while stack:
            u = stack.pop()
            comp.append(u)
            for w in adj[u]:
                if w in legal_nonhost and w not in seen:
                    seen.add(w)
                    stack.append(w)
        comps.append(comp)
    return comps


def multi_bfs(verts, adj, V):
    """Distance from a vertex SET to all vertices."""
    d = np.full(V, -1, int)
    q = collections.deque()
    for v in verts:
        d[v] = 0
        q.append(v)
    while q:
        u = q.popleft()
        for w in adj[u]:
            if d[w] < 0:
                d[w] = d[u] + 1
                q.append(w)
    return d


def cc_hist(comps, adj, V, rmax):
    hist = np.zeros(rmax + 1, int)
    for i, c in enumerate(comps):
        d = multi_bfs(c, adj, V)
        for c2 in comps[i + 1:]:
            r = int(min(d[v] for v in c2))
            if 0 < r <= rmax:
                hist[r] += 1
    return hist


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("samples", nargs="+")
    ap.add_argument("--host-classes", type=int, nargs="+", required=True,
                    help="native n6 values of the host (C15: 0 4; A15: 0 2)")
    ap.add_argument("--rmax", type=int, default=12)
    ap.add_argument("--shuffles", type=int, default=100)
    ap.add_argument("--json", default=None)
    args = ap.parse_args()

    rng = np.random.default_rng(11)
    out = {}
    for path in args.samples:
        facets = np.asarray(ddg.Manifold.load(path, 3).facets())
        n6, imp, adj = vertex_classes(facets)
        V = len(n6)
        qR = load_fields(path)[0]
        comps = joint_complexes(n6, imp, adj, args.host_classes)
        name = os.path.basename(path)
        if len(comps) < 3:
            print(f"[{name}] only {len(comps)} complexes — skipped")
            out[name] = {"n_complexes": len(comps)}
            continue

        # charge census: host baseline from vertices far from any complex
        allc = set(v for c in comps for v in c)
        dall = multi_bfs(list(allc), adj, V)
        host_far = np.where(dall >= 3)[0]
        qbar = qR[host_far].mean()
        charges = [float(qR[c].sum() - len(c) * qbar) for c in comps]
        sizes = [len(c) for c in comps]

        # mass-action spectrum
        nk = collections.Counter(sizes)
        n1 = nk.get(1, 0)
        spectrum = {}
        for k in sorted(nk):
            if n1 > 0 and k > 1:
                dF = -np.log(nk[k] / V) + k * np.log(n1 / V)
                spectrum[k] = dict(count=nk[k], dF=float(dF))
            else:
                spectrum[k] = dict(count=nk[k], dF=None)

        # g_cc vs null (random single-vertex placements, same count)
        hist = cc_hist(comps, adj, V, args.rmax)
        null = np.zeros((args.shuffles, args.rmax + 1), int)
        for s in range(args.shuffles):
            rand = [[int(v)] for v in rng.choice(V, size=len(comps), replace=False)]
            null[s] = cc_hist(rand, adj, V, args.rmax)
        nm, ns = null.mean(0), null.std(0)
        g = np.where(nm > 0, hist / np.maximum(nm, 1e-12), np.nan)

        print(f"[{name}] V={V} complexes={len(comps)} sizes={sorted(sizes, reverse=True)}")
        print(f"  charges Q: mean={np.mean(charges):+.3f} rms={np.std(charges):.3f} "
              f"(units: rad of deficit; per-vertex host qbar={qbar:+.4f})")
        print(f"  Q per complex: " + " ".join(f"{q:+.2f}" for q in sorted(charges)))
        print(f"  mass-action: " + " ".join(
            f"n{k}={v['count']}" + (f"(dF={v['dF']:.1f})" if v["dF"] else "")
            for k, v in spectrum.items()))
        print(f"  r:     " + " ".join(f"{r:>6d}" for r in range(1, args.rmax + 1)))
        print(f"  g_cc:  " + " ".join(f"{g[r]:>6.2f}" for r in range(1, args.rmax + 1)))
        ueff = [None if (np.isnan(g[r]) or g[r] <= 0) else float(-np.log(g[r]))
                for r in range(args.rmax + 1)]
        print(f"  u_eff: " + " ".join(
            f"{'  --  ' if ueff[r] is None else f'{ueff[r]:>+6.2f}'}"
            for r in range(1, args.rmax + 1)))
        out[name] = dict(V=V, n_complexes=len(comps), sizes=sizes,
                         charges=charges, qbar_host=float(qbar),
                         spectrum={str(k): v for k, v in spectrum.items()},
                         g_cc=[None if np.isnan(x) else float(x) for x in g],
                         pair_hist=hist.tolist(), null_mean=nm.tolist())

    if args.json:
        os.makedirs(os.path.dirname(args.json) or ".", exist_ok=True)
        with open(args.json, "w") as f:
            json.dump(out, f, indent=1)
        print(f"wrote {args.json}")


if __name__ == "__main__":
    main()
