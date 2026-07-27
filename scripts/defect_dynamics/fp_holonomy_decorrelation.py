#!/usr/bin/env python3
"""Worldline holonomy decorrelation: how far can a knot carry a frame?

--run: lone knot on pristine R m2, chained FP flights with the ACTUAL
jump sequences recorded (FPFlight.sample(record_path=True) -- not the
materialization tree path). The worldline is the visited chord
sequence; saved per process.

--analyze: each worldline hop (consecutive chords, one slide apart)
gets an exact Wilson line: the transport between the two sites'
canonical window tets along the canonical shortest dual path
(TransportContext; pristine development -- exact in clear regions by
no-halo). Hop quats are cached by tet pair. The decorrelation curve is
<tr R / 3> and <cos phi> of the composed line vs segment length l
(hops), with the spin-sign fraction as the double-cover diagnostic.
<tr R/3> decays 1 -> 0 (Haar) with decay length l_h = the frame-
memory length of the medium; <cos phi> decays 1 -> -1/2 (Haar).
Float prefix quaternions for the curve; exact spot-checks asserted.

Transport claims: slide channel, frozen background (v1); worldline
transport uses the pristine proxy (contact-free by construction --
lone knot).
"""
import argparse
import glob
import json
import os
import sys
from collections import deque

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
for _p in ("../../python", "../../scripts", "../../tools", "."):
    sys.path.insert(0, os.path.join(_HERE, _p))
import discrete_differential_geometry as ddg
from discrete_differential_geometry import fpkmc
from discrete_differential_geometry import development as dev
from discrete_differential_geometry.fpkmc import face_apex_maps
from fpkmc_v3_hb import build

REF = "data/tcp_reference/T3_R_m2_N7248.mfd"


def run(args):
    m, s, apx, window = build(args.ref, args.estar, args.lam,
                              window=args.window)
    rng = np.random.default_rng(args.seed)
    nu = fpkmc.nu_per_attempt(1.0, s.manifold.num_facets)
    chord = tuple(sorted(apx))
    chords = [list(chord)]
    for fl_i in range(args.flights):
        g = s.slide_graph_scan(chord, dS_max=args.dS_max,
                               max_depth=args.depth)
        fl = fpkmc.FPFlight(g, args.dS_max, args.depth, nu)
        j, t, jumps, nodes = fl.sample(rng, record_path=True)
        # worldline chords: visited single-chord nodes (skip node 0 --
        # it's the previous exit)
        stop = False
        for n in nodes[1:]:
            if int(g["n_chords"][n]) != 1:
                stop = True
                break
            chords.append(sorted(int(x) for x in g["chord"][n]))
        # materialize the exit to continue (tree path; endpoint-exact)
        par = fpkmc.HBDriver._parents(g)
        path = fpkmc.HBDriver._path(g, par, j)
        for i in path:
            ch = (int(g["edge_chord"][i][0]), int(g["edge_chord"][i][1]))
            dS = s.slide_at(ch[0], ch[1], int(g["edge_slot"][i]),
                            commit=True)
            if dS is None:
                raise RuntimeError("materialize failed")
        if stop:
            print(f"  flight {fl_i}: multichord node on worldline; "
                  f"ending run", flush=True)
            break
        chord = tuple(sorted((int(g["chord"][j][0]),
                              int(g["chord"][j][1]))))
        if (fl_i + 1) % 20 == 0:
            print(f"  flight {fl_i + 1}/{args.flights}: "
                  f"{len(chords)} worldline chords", flush=True)
    with open(args.out, "w") as fh:
        json.dump({"seed": args.seed, "lam": args.lam,
                   "chords": chords}, fh)
    print(f"wrote {args.out} ({len(chords)} chords)")


def analyze(args):
    m = ddg.Manifold.load(args.ref, 3)
    F = np.asarray(m.facets())
    _, face_of = face_apex_maps(m)
    ctx = dev.TransportContext([tuple(int(x) for x in f) for f in F])

    def site_tet(chord):
        c = tuple(sorted(chord))
        f = face_of[c]
        t1 = tuple(sorted((c[0],) + tuple(f)))
        t2 = tuple(sorted((c[1],) + tuple(f)))
        return frozenset(min(t1, t2))

    def bfs_path(a, b):
        if a == b:
            return [a]
        prev = {a: None}
        q = deque([a])
        while q:
            t = q.popleft()
            for nb in ctx.dual[t]:
                if nb in prev:
                    continue
                prev[nb] = t
                if nb == b:
                    out = []
                    while nb is not None:
                        out.append(nb)
                        nb = prev[nb]
                    return list(reversed(out))
                q.append(nb)
        raise RuntimeError("no dual path")

    hop_cache = {}

    def hop_quat(c1, c2):
        t1, t2 = site_tet(c1), site_tet(c2)
        key = (t1, t2)
        if key not in hop_cache:
            hop_cache[key] = ctx.wilson_line(bfs_path(t1, t2))
        return hop_cache[key]

    # ---- collect worldlines, compute hop quats ----
    seqs = []
    for f in sorted(glob.glob(args.glob)):
        d = json.load(open(f))
        ch = [tuple(c) for c in d["chords"]]
        hops = []
        for k in range(len(ch) - 1):
            hops.append(hop_quat(ch[k], ch[k + 1]))
        seqs.append(hops)
        print(f"  {os.path.basename(f)}: {len(hops)} hops "
              f"({len(hop_cache)} cached tet pairs)", flush=True)

    # exact spot-check on one short window per sequence
    for hops in seqs[:2]:
        W = hops[0]
        for q in hops[1:8]:
            W = q * W
        fq = _to_float(hops[0])
        for q in hops[1:8]:
            fq = _fmul(_to_float(q), fq)
        assert abs(float(W.cos_phi()) - (2 * fq[0] ** 2 - 1)) < 1e-9

    # ---- decorrelation curves (float prefixes) ----
    # RAW Wilson line W is dominated by the deterministic local twist
    # (a pure chain walk rotates hugely yet has PERFECT memory). The
    # memory observable is the fluctuation loop H = D^{-1} W: worldline
    # transport against the canonical shortest-path transport between
    # the same endpoint sites -- the absorber's frame-reconstruction
    # error given only the endpoints. <tr H>/3: 1 -> 0 (Haar).
    Ls = [1, 2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64, 96, 128, 192, 256]
    acc = {L: [] for L in Ls}       # loop H trace
    raw = {L: [] for L in Ls}       # raw W trace (for the record)
    sign = {L: [] for L in Ls}      # spin sign of H
    seq_chords = []
    for f in sorted(glob.glob(args.glob)):
        seq_chords.append([tuple(c) for c in json.load(open(f))["chords"]])
    for hops, chs in zip(seqs, seq_chords):
        n = len(hops)
        pre = [np.array([1.0, 0, 0, 0])]
        for q in hops:
            p = _fmul(_to_float(q), pre[-1])
            pre.append(p / np.linalg.norm(p))
        for L in Ls:
            if L >= n:
                continue
            stride = max(1, (n - L) // 300, L // 2)
            for i in range(0, n - L, stride):
                W = _fmul(pre[i + L], _fconj(pre[i]))
                raw[L].append(1.0 + 2.0 * (2 * W[0] ** 2 - 1))
                t1 = site_tet(chs[i])
                t2 = site_tet(chs[i + L])
                key = (t1, t2)
                if key not in hop_cache:
                    hop_cache[key] = ctx.wilson_line(bfs_path(t1, t2))
                D = _to_float(hop_cache[key])
                H = _fmul(_fconj(D), W)
                cphi = 2 * H[0] ** 2 - 1
                acc[L].append(1.0 + 2.0 * cphi)
                sign[L].append(1.0 if H[0] >= 0 else 0.0)

    print(f"\n{'l (hops)':>9} {'<tr H>/3':>10} {'SE':>7} "
          f"{'raw <tr W>/3':>13} {'P(+lift H)':>10} {'n':>6}")
    xs, ys = [], []
    for L in Ls:
        a = np.array(acc[L])
        if len(a) < 10:
            continue
        tr3 = a.mean() / 3.0
        se = a.std() / np.sqrt(len(a)) / 3.0
        rw = np.mean(raw[L]) / 3.0
        ps = np.mean(sign[L])
        print(f"{L:9d} {tr3:10.4f} {se:7.4f} {rw:13.4f} {ps:10.3f} "
              f"{len(a):6d}")
        if tr3 > 0.02:
            xs.append(L)
            ys.append(tr3)
    xs, ys = np.array(xs, float), np.array(ys)
    mfit = ys > 0
    lam_h = None
    if mfit.sum() >= 3:
        k, b = np.polyfit(xs[mfit], np.log(ys[mfit]), 1)
        lam_h = -1.0 / k
        print(f"\nframe-memory length l_h = {lam_h:.1f} hops "
              f"(exp fit of <tr R>/3)")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(8, 5))
    L_all = [L for L in Ls if len(acc[L]) >= 10]
    ax.errorbar(L_all, [np.mean(acc[L]) / 3 for L in L_all],
                yerr=[np.std(acc[L]) / np.sqrt(len(acc[L])) / 3
                      for L in L_all],
                fmt="o-", label="<tr H>/3 (frame memory)")
    ax.plot(L_all, [np.mean(raw[L]) / 3 for L in L_all], "s--",
            alpha=0.5, label="raw <tr W>/3")
    ax.plot(L_all, [np.mean(sign[L]) for L in L_all], "^:",
            label="P(+ spin lift of H)")
    ax.axhline(0, color="k", lw=0.5)
    ax.axhline(0.5, color="gray", lw=0.5, ls="-.")
    if lam_h:
        ax.plot(xs, np.exp(b) * np.exp(-xs / lam_h), "r-", lw=1,
                label=f"exp fit: l_h = {lam_h:.1f} hops")
    ax.set_xscale("log")
    ax.set_xlabel("worldline segment length l (slide hops)")
    ax.set_ylabel("holonomy observables")
    ax.set_title("Worldline frame memory: H = D^-1 W (worldline vs "
                 "direct transport) -- lone knot, R m2,\nlam=0.40, "
                 "slide channel, frozen bg; pristine development, "
                 "canonical shortest paths", fontsize=9)
    ax.legend(fontsize=8)
    fig.tight_layout()
    out = "data/fpkmc/fp_holonomy_decorrelation.png"
    fig.savefig(out, dpi=140)
    print(f"Saved to: {os.path.abspath(out)}")


def _to_float(q):
    v = np.array([q.a, q.b, q.c, q.d], float)
    return v / np.linalg.norm(v)


def _fmul(p, q):
    a1, b1, c1, d1 = p
    a2, b2, c2, d2 = q
    return np.array([a1 * a2 - b1 * b2 - c1 * c2 - d1 * d2,
                     a1 * b2 + b1 * a2 + c1 * d2 - d1 * c2,
                     a1 * c2 - b1 * d2 + c1 * a2 + d1 * b2,
                     a1 * d2 + b1 * c2 - c1 * b2 + d1 * a2])


def _fconj(p):
    return np.array([p[0], -p[1], -p[2], -p[3]])


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--run", action="store_true")
    ap.add_argument("--analyze", action="store_true")
    ap.add_argument("--ref", default=REF)
    ap.add_argument("--lam", type=float, default=0.40)
    ap.add_argument("--estar", type=float, default=5.105025)
    ap.add_argument("--window", type=int, default=40)
    ap.add_argument("--dS-max", type=float, default=5.0)
    ap.add_argument("--depth", type=int, default=3)
    ap.add_argument("--flights", type=int, default=100)
    ap.add_argument("--seed", type=int, default=61)
    ap.add_argument("--out", default=None)
    ap.add_argument("--glob", default="data/fpkmc/holo_w*.json")
    args = ap.parse_args()
    if args.run:
        run(args)
    elif args.analyze:
        analyze(args)
    else:
        raise SystemExit("need --run or --analyze")


if __name__ == "__main__":
    main()
