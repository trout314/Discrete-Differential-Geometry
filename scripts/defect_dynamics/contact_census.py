#!/usr/bin/env python3
"""Census of REJECTED CONTACTS: what could a contact handoff actually hand to?

Motivation (2026-08-05). In the blob A/B, 41.8% of flight kernel steps end in
`reject_contact_flip` -- a proposed slide whose support touches a defect,
priced by Metropolis (docking is repulsive, +4.4..+12.8 from the collider
work), rejected, and the lift REVERSED. In factorized/soft ECMC a rejected
move does not reverse the lift; it hands it to the partner that caused the
rejection. Those ~400k reversals are the momentum transfers we are throwing
away.

Before building that rule, measure what is actually reachable at a contact.
The existing transmit/absorb census (14.5% at lam=0.35, 28.1% at lam=0.40)
was taken on the R m4 melt, NOT on this C15 m3 blob, which is far denser --
so it does not transfer.

Per rejected contact this records:

  blocker      which background-illegal complex the support touches, its
               size, and whether it holds a deg-3 edge (= can it carry the
               lift at all)
  carrier      a handoff target must be an EXISTING deg-3 chord (the active
               object is one, and a slide TARGET only becomes deg-3 after
               the move -- so "walk on and find a deg-3 arrival" is
               vacuous, it is always a plain crystal edge). Recorded: does
               the blocking complex hold one, and if not, the graph distance
               from the contact support to the nearest background deg-3 edge
               -- i.e. how far a handoff would have to reach.
  resume       the first clean same-rung (dS = 0) site past the blocker --
               the pass-through option, which maximises displacement but
               transfers no momentum
  dS           how repulsive the rejected move was (is the rejection
               marginal or hard?)

Answers, in order, design questions (a) does the lift stop at the struck
defect or pass beyond, (b) what to do in the absorb case, (c) must the
handoff be free, (d) how far may hscan reach.

Usage:
    python scripts/defect_dynamics/contact_census.py --seeds 2 --sweeps 1000 \
        --out data/ecmc_ab/contact_census.json
"""
import argparse
import collections
import json
import os
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(os.path.dirname(_HERE))
sys.path.insert(0, os.path.join(_ROOT, "python"))
sys.path.insert(0, _HERE)

import discrete_differential_geometry as ddg
from discrete_differential_geometry import Manifold, ManifoldSampler, SamplerParams
from ecmc_flight import Flight, ed_of
from ecmc_ab import build_blob

REF = os.path.join(_ROOT, "data/tcp_reference/T3_C15_m3_N3672.mfd")


class CensusFlight(Flight):
    """Flight that records the anatomy of every rejected contact."""

    def __init__(self, *a, hscan=40, sink=None, sw=0, cimp=0.5,
                 imp_offset=0, imp_lin=0.0, **kw):
        super().__init__(*a, **kw)
        self.cimp = float(cimp)
        self.imp_offset = int(imp_offset)
        self.imp_lin = float(imp_lin)
        self.hscan = hscan
        self.sink = sink if sink is not None else []
        self.sw = sw
        self._comp_cache = (None, None)

    # -- background-illegal complexes, cached per background -----------------
    def complexes(self):
        """{vertex -> component id}, [component vertex sets], from the
        BACKGROUND illegal edges (the state minus the active flier)."""
        key = self.ill_snapshot
        if self._comp_cache[0] == key:
            return self._comp_cache[1]
        adj = collections.defaultdict(set)
        for e in self.ill_snapshot:
            adj[e[0]].add(e[1])
            adj[e[1]].add(e[0])
        seen, comp_of, comps = set(), {}, []
        for v in adj:
            if v in seen:
                continue
            stack, cur = [v], set()
            while stack:
                x = stack.pop()
                if x in cur:
                    continue
                cur.add(x)
                stack.extend(adj[x] - cur)
            for x in cur:
                comp_of[x] = len(comps)
            comps.append(cur)
            seen |= cur
        out = (comp_of, comps)
        self._comp_cache = (key, out)
        return out

    def deg3_edges_of(self, verts):
        """Background deg-3 edges with BOTH ends inside `verts`."""
        return [e for e, d in self.edB.items()
                if d == 3 and e[0] in verts and e[1] in verts]

    def adjacency(self):
        """1-skeleton adjacency of the background, cached per background."""
        if getattr(self, "_adj_key", None) == self.ill_snapshot:
            return self._adj
        adj = collections.defaultdict(set)
        for (u, v) in self.edB:
            adj[u].add(v)
            adj[v].add(u)
        self._adj_key, self._adj = self.ill_snapshot, adj
        return adj

    def dist_to_nearest_d3(self, src, exclude, cap=8):
        """Graph distance from any vertex of `src` to the nearest endpoint of
        a background deg-3 edge not inside `exclude`. None if beyond `cap`."""
        targets = {v for e, d in self.edB.items() if d == 3
                   and not (e[0] in exclude and e[1] in exclude)
                   for v in e}
        if not targets:
            return None
        if src & targets:
            return 0
        adj, seen, frontier = self.adjacency(), set(src), set(src)
        for r in range(1, cap + 1):
            frontier = {w for v in frontier for w in adj[v]} - seen
            if not frontier:
                return None
            if frontier & targets:
                return r
            seen |= frontier
        return None

    def on_reject(self, kind, k, arr, face, dS, rng):
        if kind != "contact":
            return False
        comp_of, comps = self.complexes()
        sup = set(face) | set(arr)
        hits = sorted(sup & self.BV)
        cids = sorted({comp_of[v] for v in hits if v in comp_of})
        blk = set().union(*[comps[c] for c in cids]) if cids else set()
        blk_d3 = self.deg3_edges_of(blk)

        # walk to the contact, then continue THROUGH it
        ww = self.w
        for _ in range(k):
            ww = self.stepw(ww)
            if ww is None:
                break
        resume_k = None
        if ww is not None:
            for j in range(1, self.hscan + 1):
                ww = self.stepw(ww)
                if ww is None:
                    break
                a2 = self.arrival(ww)
                if a2 is None:
                    break
                f2 = tuple(sorted(ww[1:4]))
                if not ((set(f2) | set(a2)) & self.BV) \
                        and self._rung_of_face(f2, a2) == self.Q:
                    resume_k = j
                    break

        # Decompose the barrier: how much of dS is the m^2 anti-clustering
        # term we CHOSE (V(m) = cimp * m^2, m = #incident edges of degree
        # outside {5,6}), versus the pins + geometry? Re-price the same
        # proposal with cimp temporarily zeroed. Setting cimp=0 for the whole
        # RUN is not a control -- the state degenerates (rungs run negative
        # and the walk dies), because the m^2 term is what holds the gas.
        dS_pins = None
        if self.cimp or self.imp_lin:
            slot = self._slot_for(self.w)
            if slot is not None:
                self.s.set_n6_potential(0.0, 0.0, None, self.imp_offset, 0.0)
                try:
                    r0 = self.s.nonlocal_slide_at(self.C[0], self.C[1],
                                                  slot, k, commit=False)
                    if r0 is not None:
                        dS_pins = float(r0[0])
                finally:
                    self.s.set_n6_potential(0.0, self.cimp, None,
                                            self.imp_offset, self.imp_lin)

        n_d3_bg = sum(1 for d in self.edB.values() if d == 3)
        self.sink.append({
            "sw": self.sw, "k": int(k), "dS": float(dS), "Q": int(self.Q),
            "n_hit": len(hits), "n_comp": len(cids),
            "blk_size": len(blk),
            "blk_has_d3": bool(blk_d3), "blk_n_d3": len(blk_d3),
            "n_d3_bg": n_d3_bg,
            "d3_dist": self.dist_to_nearest_d3(sup, blk),
            "dS_pins": dS_pins,
            "dS_m2": None if dS_pins is None else float(dS) - dS_pins,
            "resume_k": resume_k,
        })
        return False        # census only: still reverse, as production does


def run_chain(seed, cfg):
    rng = np.random.default_rng(seed)
    ddg.set_random_seed(int(rng.integers(2 ** 31)))
    m = Manifold.load(REF, 3)
    blob, center = build_blob(m, cfg["k_fliers"], rng, radius=cfg["radius"])
    fv = [int(x) for x in m.f_vector]
    estar = 6.0 * fv[3] / fv[1]
    s = ManifoldSampler(m, SamplerParams(
        num_facets_target=fv[3], num_facets_coef=cfg["nfc"],
        hinge_degree_target=estar, num_hinges_coef=cfg["k1"],
        hinge_degree_variance_coef=0.0, codim3_degree_variance_coef=0.0,
        hinge_degree_target_coef=0.0))
    s.set_n6_potential(0.0, cfg["cimp"], None, cfg.get("imp_offset", 0),
                       cfg.get("imp_lin", 0.0))

    sink, eps = [], []
    nchunks = cfg["sweeps"] // cfg["chunk"]
    for ci in range(nchunks):
        s.run(sweeps=cfg["chunk"])
        sw = (ci + 1) * cfg["chunk"]
        for _ in range(cfg["eps_per_chunk"]):
            pairs, degs = s.manifold.illegal_edges()
            pairs = np.asarray(pairs).reshape(-1, 2)
            d3 = [tuple(sorted(map(int, p)))
                  for p, d in zip(pairs, np.asarray(degs)) if d == 3]
            if not d3:
                continue
            C = d3[int(rng.integers(len(d3)))]
            lk = set()
            for t in np.asarray(s.manifold.facets()):
                t = [int(x) for x in t]
                if C[0] in t and C[1] in t:
                    lk |= {v for v in t if v not in C}
            if len(lk) != 3:
                continue
            try:
                fl = CensusFlight(s, C, (C[0],) + tuple(sorted(lk)),
                                  kscan=cfg["kscan"], audit=False,
                                  beta=1.0, p_hand=0.0, expect_free=False,
                                  hscan=cfg["hscan"], sink=sink, sw=sw,
                                  cimp=cfg["cimp"],
                                  imp_offset=cfg.get("imp_offset", 0),
                                  imp_lin=cfg.get("imp_lin", 0.0))
                fl.refresh_frame(rng)
                for i in range(cfg["ep_steps"]):
                    if cfg["refresh_every"] and i and i % cfg["refresh_every"] == 0:
                        fl.refresh_frame(rng)
                    fl.step(rng)
                eps.append({"sw": sw, "Q": int(fl.Q),
                            "events": dict(fl.events)})
            except Exception as ex:                      # noqa: BLE001
                print(f"  seed {seed} sw {sw}: episode failed: {ex}")
    return sink, eps


def report(rows, label):
    n = len(rows)
    if not n:
        print(f"\n=== {label}: no rejected contacts")
        return {}
    A = lambda k: np.array([r[k] for r in rows], float)          # noqa: E731
    res = np.array([r["resume_k"] if r["resume_k"] is not None else -1
                    for r in rows])
    dd = np.array([r["d3_dist"] if r["d3_dist"] is not None else -1
                   for r in rows])
    dS = A("dS")
    out = {
        "n": n,
        "blk_has_d3": float(A("blk_has_d3").mean()),
        "blk_size_mean": float(A("blk_size").mean()),
        "blk_n_d3_mean": float(A("blk_n_d3").mean()),
        "n_d3_bg_mean": float(A("n_d3_bg").mean()),
        "d3_within_cap": float((dd >= 0).mean()),
        "d3_dist_median": float(np.median(dd[dd >= 0])) if (dd >= 0).any() else None,
        "resume_found": float((res > 0).mean()),
        "resume_k_median": float(np.median(res[res > 0])) if (res > 0).any() else None,
        "dS_median": float(np.median(dS)),
        "dS_q10": float(np.quantile(dS, 0.10)),
        "dS_frac_lt_1": float((dS < 1.0).mean()),
        "dS_frac_lt_3": float((dS < 3.0).mean()),
    }
    print(f"\n=== {label}   ({n:,} rejected contacts)")
    print(f"  blocking complex size            {out['blk_size_mean']:.1f} vertices"
          f"   ({out['n_d3_bg_mean']:.1f} deg-3 chords in the whole background)")
    print(f"  (a) blocker HOLDS a deg-3 edge   {100*out['blk_has_d3']:5.1f}%"
          f"   (mean {out['blk_n_d3_mean']:.2f})  <- ceiling for 'struck defect carries'")
    print(f"      else nearest deg-3 elsewhere {100*out['d3_within_cap']:5.1f}% within cap,"
          f" median {out['d3_dist_median']} graph steps")
    print(f"  (b) RESUME clean same-rung past  {100*out['resume_found']:5.1f}%"
          f"   median {out['resume_k_median']} walk steps  <- ceiling for pass-through")
    print(f"  (c) rejected dS: median {out['dS_median']:.2f}, q10 {out['dS_q10']:.2f}, "
          f"frac<1 = {100*out['dS_frac_lt_1']:.1f}%, frac<3 = {100*out['dS_frac_lt_3']:.1f}%")
    dec = [r for r in rows if r.get("dS_pins") is not None]
    if dec:
        pin = np.array([r["dS_pins"] for r in dec])
        m2 = np.array([r["dS_m2"] for r in dec])
        tot = pin + m2
        out["dS_pins_median"] = float(np.median(pin))
        out["dS_m2_median"] = float(np.median(m2))
        out["m2_share"] = float(np.median(m2 / np.where(tot != 0, tot, np.nan)))
        out["frac_pins_lt_3"] = float((pin < 3.0).mean())
        print(f"      decomposed (n={len(dec):,}): pins+geometry median "
              f"{out['dS_pins_median']:.2f}, m^2 term median "
              f"{out['dS_m2_median']:.2f}  -> m^2 is "
              f"{100*out['m2_share']:.0f}% of the barrier")
        print(f"      WITHOUT the m^2 term, {100*out['frac_pins_lt_3']:.1f}% "
              f"of these contacts would cost < 3")
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--seeds", type=int, default=2)
    ap.add_argument("--seed0", type=int, default=1000)
    ap.add_argument("--sweeps", type=int, default=1000)
    ap.add_argument("--chunk", type=int, default=20)
    ap.add_argument("--eps-per-chunk", type=int, default=4)
    ap.add_argument("--ep-steps", type=int, default=40)
    ap.add_argument("--refresh-every", type=int, default=10)
    ap.add_argument("--kscan", type=int, default=60)
    ap.add_argument("--hscan", type=int, default=40)
    ap.add_argument("--k-fliers", type=int, default=12)
    ap.add_argument("--radius", type=int, default=4)
    ap.add_argument("--cimp", type=float, default=0.5)
    ap.add_argument("--imp-lin", type=float, default=0.0,
                    help="pure chemical potential on impure edges (imp_lin*m)")
    ap.add_argument("--imp-offset", type=int, default=0,
                    help="flat foot: V(m) = cimp * max(0, m - offset)^2")
    ap.add_argument("--nfc", type=float, default=30.0)
    ap.add_argument("--k1", type=float, default=1.0)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    cfg = vars(args)

    rows, episodes = [], []
    for i in range(args.seeds):
        print(f"[chain {i}] seed {args.seed0 + i} ...", flush=True)
        r, e = run_chain(args.seed0 + i, cfg)
        rows += r
        episodes += e

    summ = {"all": report(rows, "ALL rejected contacts")}
    half = args.sweeps // 2
    summ["early"] = report([r for r in rows if r["sw"] <= half],
                           f"EARLY (sw <= {half}, blob still packed)")
    summ["late"] = report([r for r in rows if r["sw"] > half],
                          f"LATE  (sw >  {half}, dispersing)")

    # --- stratify by RUNG. Q is conserved for a whole episode (no Q-refresh),
    # and the per-rung free network is NOT uniform -- some rungs carry a
    # crystal-spanning free web, others fragment -- so pooling rungs hides
    # exactly the sector dependence that decides whether a handoff is
    # available. See notes/memory/bc-washboard-not-free-spirals.md.
    qs = collections.Counter(r["Q"] for r in rows)
    print(f"\n--- rung occupancy of the sampled flights: "
          + ", ".join(f"Q={q}: {n:,} ({100*n/len(rows):.1f}%)"
                      for q, n in sorted(qs.items())))
    summ["by_Q"] = {}
    for q, n in sorted(qs.items()):
        if n >= 50:
            summ["by_Q"][str(q)] = report([r for r in rows if r["Q"] == q],
                                          f"RUNG Q={q}")

    # per-rung kernel-step budget: where does a flight in this sector spend
    # its steps? A sector whose free web does not span shows up here as a
    # collapsed accept_free share and a runaway reject_contact share.
    per_q = collections.defaultdict(collections.Counter)
    nep = collections.Counter()
    for e in episodes:
        per_q[e["Q"]].update(e["events"])
        nep[e["Q"]] += 1
    print("\n--- per-rung kernel-step budget (share of that rung's steps)")
    print(f"    {'Q':>4} {'eps':>6} {'steps':>8} {'free':>8} {'rejctc':>8} "
          f"{'accctc':>8} {'noslot':>8}")
    for q in sorted(per_q):
        c = per_q[q]
        st = sum(v for k2, v in c.items()
                 if k2.startswith(("accept_", "reject_")) or k2 in
                 ("noslot_flip", "stuck_flip", "illegal_flip"))
        if st == 0:
            continue
        print(f"    {q:>4} {nep[q]:>6} {st:>8,} "
              f"{c['accept_free']/st:>7.1%} {c['reject_contact_flip']/st:>7.1%} "
              f"{c['accept_contact']/st:>7.1%} {c['noslot_flip']/st:>7.1%}")
    summ["per_rung_steps"] = {str(q): dict(per_q[q]) for q in per_q}

    if args.out:
        os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
        with open(args.out, "w") as f:
            json.dump({"summary": summ, "cfg": cfg, "rows": rows,
                       "episodes": episodes}, f,
                      indent=1, default=float)
        print(f"\nwrote {os.path.abspath(args.out)}")


if __name__ == "__main__":
    main()
