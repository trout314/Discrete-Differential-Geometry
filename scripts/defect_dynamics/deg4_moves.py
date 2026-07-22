#!/usr/bin/env python3
"""Catch realized minimal disclination moves at full accepted-move resolution.

Run a constrained R chain from a DEFECTED start (so defects are active), replay
each accepted move in lockstep, maintain edge degrees + adjacency + vertex->tet
incidence incrementally, and:
  (1) DIAGNOSE how often a single ISOLATED deg-4 edge (the minimal disclination
      unit) is even present -- our size-2-never-exists result predicts ~never;
  (2) capture the closed star of the affected vertices BEFORE/AFTER every move
      that changes the local ILLEGAL-EDGE set; keep framing-preserving pairs as
      realized block moves, classified by block defect content.

args: ref_cell start_state bump cimp zleg sweeps chunk seed
"""
import os
import sys
from collections import Counter, defaultdict
from itertools import combinations

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
for p in ("../../python", "../../scripts", "../../tools", "."):
    sys.path.insert(0, os.path.join(_HERE, p))
import discrete_differential_geometry as ddg
import fk_moves as fk
import event_replay as er
from fk_skeleton import edges_from_facets

a = sys.argv
REF, START = a[1], a[2]
BUMP = float(a[3]) if len(a) > 3 else 8e-4
CIMP = float(a[4]) if len(a) > 4 else 1.0
ZLEG = float(a[5]) if len(a) > 5 else 0.6
SWEEPS = int(a[6]) if len(a) > 6 else 400
CHUNK = int(a[7]) if len(a) > 7 else 5
SEED = int(a[8]) if len(a) > 8 else 1


def ek(u, w):
    return (u, w) if u < w else (w, u)


def isolated_deg4(e, edeg, nbr):
    if edeg.get(e, 0) != 4:
        return False
    u, w = e
    for x in nbr[u]:
        if x != w and edeg[ek(u, x)] not in (5, 6):
            return False
    for x in nbr[w]:
        if x != u and edeg[ek(w, x)] not in (5, 6):
            return False
    return True


ddg.set_random_seed(SEED)
ref = ddg.Manifold.load(REF, 3)
native = float(edges_from_facets(ref.facets())[1].mean())
m = ddg.Manifold.load(START, 3)
params = ddg.SamplerParams(
    num_facets_target=ref.num_facets, num_facets_coef=0.1,
    hinge_degree_target=native + BUMP, num_hinges_coef=0.0,
    hinge_degree_variance_coef=0.0, codim3_degree_variance_coef=0.0,
    hinge_degree_target_coef=1.0 * (native + BUMP) / 6.0)
s = ddg.ManifoldSampler(m, params)
s.set_n6_potential(ZLEG, CIMP, tilt=[0.0] * 5)
v = s.manifold
s.enable_event_log(64.0)
s.drain_event_log()

tets = {tuple(sorted(int(x) for x in t)) for t in v.facets()}
edeg = Counter()
nbr = defaultdict(set)
v2t = defaultdict(set)
for t in tets:
    for e in combinations(t, 2):
        edeg[e] += 1
        nbr[e[0]].add(e[1]); nbr[e[1]].add(e[0])
    for w in t:
        v2t[w].add(t)


def all_isolated():
    iso = 0
    for e, d in edeg.items():
        if d == 4 and isolated_deg4(e, edeg, nbr):
            iso += 1
    return iso


captures = []
n_events = n_trans = 0
max_iso = all_isolated()
n_deg4_edges_seen = Counter()

for c in range(SWEEPS // CHUNK):
    s.run(sweeps=CHUNK)
    ev = s.drain_event_log()
    if s.event_log_overflowed():
        print("  WARN event overflow", flush=True)
    for e in ev:
        n_events += 1
        rem, add = er.event_changes(e)
        AV = set().union(*[set(t) for t in rem + add])
        before_edges = {ek(*p) for t in rem for p in combinations(t, 2)}
        before_iso = {x for x in before_edges if isolated_deg4(x, edeg, nbr)}
        before_illegal = {x for x in before_edges if edeg.get(x, 0) not in (5, 6, 0)}
        before_tets = [t for w in AV for t in v2t[w]]
        before_tets = list({t for t in before_tets})
        for r in rem:
            tets.discard(r)
            for w in r:
                v2t[w].discard(r)
            for p in combinations(r, 2):
                edeg[p] -= 1
                if edeg[p] == 0:
                    nbr[p[0]].discard(p[1]); nbr[p[1]].discard(p[0])
        for aa in add:
            tets.add(aa)
            for w in aa:
                v2t[w].add(aa)
            for p in combinations(aa, 2):
                if edeg[p] == 0:
                    nbr[p[0]].add(p[1]); nbr[p[1]].add(p[0])
                edeg[p] += 1
        after_edges = {ek(*p) for t in add for p in combinations(t, 2)}
        after_iso = {x for x in (after_edges | before_iso)
                     if isolated_deg4(x, edeg, nbr)}
        changed = before_iso ^ after_iso            # isolated deg-4 edges flipped
        if not changed:
            continue
        n_trans += 1
        # MINIMAL block: closed star of the changed edge's 2 endpoints only.
        for e in changed:
            u, w = e
            before_block = [t for t in before_tets if u in t or w in t]
            after_block = list(v2t[u] | v2t[w])
            captures.append((e in after_iso, before_block, after_block))
    if c % 40 == 0:
        max_iso = max(max_iso, all_isolated())

print(f"replayed {n_events} accepted moves; {n_trans} local illegal-edge changes",
      flush=True)
print(f"MAX simultaneous ISOLATED deg-4 edges ever seen: {max_iso}", flush=True)

# classify framing-preserving block moves
kinds = Counter()
moves = {}
for created, bt, at in captures:
    ab, aa = fk.analyze(bt), fk.analyze(at)
    okb, spb = fk.interior_ok(ab, "deg4")
    oka, spa = fk.interior_ok(aa, "deg4")
    if not (okb and oka):
        kinds["not_minimal_block"] += 1
        continue
    if fk.framed_boundary_canon(ab) != fk.framed_boundary_canon(aa):
        kinds["boundary_changed"] += 1
        continue
    bcb, bca = fk.ball_canon(ab), fk.ball_canon(aa)
    if bcb == bca:
        kinds["no_op"] += 1
        continue
    label = {("vac", "deg4"): "CREATION", ("deg4", "vac"): "CREATION",
             ("deg4", "deg4"): "HOP", ("vac", "vac"): "vacuum"}.get(
        tuple(sorted((spb, spa))), f"{spb}->{spa}")
    kinds[label] += 1
    moves.setdefault((fk.framed_boundary_canon(ab), frozenset((bcb, bca))),
                     (label, spb, spa, len(bt), len(at)))

print("\ncapture classification:")
for k, n in kinds.most_common():
    print(f"  {k:22s}: {n}")
print(f"\nDISTINCT minimal (single-deg4) framing-preserving block moves: {len(moves)}")
for i, (_, mv) in enumerate(moves.items()):
    label, spb, spa, nb, na = mv
    print(f"  move {i}: {label}  {spb}->{spa}  ({nb}->{na} tets)")
