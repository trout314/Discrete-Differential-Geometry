#!/usr/bin/env python3
"""Does a same-rung non-local slide really cost exactly zero?  (It does.)

Predicts dS = c*(Q_target - Q_source) from the PRISTINE rung table and checks
it against the SAMPLER's exact arithmetic via nonlocal_slide_at -- so this
validates ecmc.face_rung against the D core rather than against itself.

Why it matters: the non-local slide is a 3->2 at the source and a 2->3 at the
target, so the manifold changes only at the endpoints and dS depends on the
endpoints ALONE. Nothing is traversed, so there is no barrier between them: a
hop to any same-rung site is free at any distance. Proposing a uniform STEP
lands off-rung most of the time (mean acceptance 0.307); proposing a uniform
same-rung TARGET is accepted with probability exactly 1.

Measured 2026-08-02, R m2, pure edge action:

    468 legal (slot, steps) targets     max |dS - c*dQ| = 0.000e+00  (EXACT)
    124 same-rung targets               max |dS|        = 0.000e+00
                                        acceptance      = 1.000000000000
                                        displacements   = 1..15 chain steps
    mean acceptance, uniform proposal   0.307
    mean acceptance, same-rung proposal 1.000

GOTCHA: the vdv coefficients must be zeroed explicitly (num_hinges_coef,
hinge_degree_variance_coef, codim3_degree_variance_coef) or the residual
reappears and same-rung hops stop being exactly free -- SamplerParams zeroes
only what you pass it.
"""
import sys, os
sys.path.insert(0, "python"); sys.path.insert(0, "scripts")
import numpy as np
import discrete_differential_geometry as ddg
from discrete_differential_geometry import Manifold, ManifoldSampler, SamplerParams
from discrete_differential_geometry.ecmc import face_rung

P = "data/tcp_reference/T3_R_m2_N7248.mfd"
ESTAR = 5.105025

m0 = Manifold.load(P, 3)
edeg0 = {tuple(sorted(map(int, e))): m0.degree(e) for e in m0.simplices(1)}
F = np.asarray(m0.facets())
f2a = {}
for t in F:
    t = tuple(sorted(int(x) for x in t))
    for i in range(4):
        f2a.setdefault(t[:i] + t[i + 1:], []).append(t[i])
RUNG = {f: face_rung(f, tuple(ap), edeg0) for f, ap in f2a.items() if len(ap) == 2}
# a 2->3 on face f with apexes (x,y) creates the chord (x,y), so an arrival
# chord names its face -- provided the apex pair is unique per face, asserted:
A2F = {}
for f, ap in f2a.items():
    if len(ap) == 2:
        k = tuple(sorted(int(x) for x in ap))
        assert k not in A2F, f"apex pair {k} is shared by two faces"
        A2F[k] = f
print(f"apex-pair -> face map is injective over {len(A2F)} faces")

# PURE EDGE ACTION -- the contamination gotcha: zero the vdv coefficients
p = SamplerParams(num_facets_target=len(F), num_facets_coef=0.0,
                  hinge_degree_target=ESTAR, hinge_degree_target_coef=1.0,
                  num_hinges_coef=0.0, hinge_degree_variance_coef=0.0,
                  codim3_degree_variance_coef=0.0)

# create one 2->3 defect at a chosen face
seed_face = sorted(RUNG, key=lambda f: (RUNG[f], f))[0]
apex = tuple(f2a[seed_face])
m = Manifold.load(P, 3)
m.do_bistellar_move(list(seed_face), list(apex))
s = ManifoldSampler(m, p)
chord = tuple(sorted(apex))
print(f"source face {seed_face} rung Q={RUNG[seed_face]}  -> chord {chord}")

c = None
rows = []
for slot in range(12):
    for steps in range(1, 40):
        r = s.nonlocal_slide_at(chord[0], chord[1], slot, steps, commit=False)
        if r is None:
            continue
        dS, dn3, arr = r
        if arr is None:
            continue
        tgt = A2F.get(tuple(sorted(int(x) for x in arr)))
        if tgt is None:
            continue
        rows.append((slot, steps, dS, RUNG[seed_face], RUNG[tgt], dn3))

dq = np.array([r[4] - r[3] for r in rows], float)
ds = np.array([r[2] for r in rows], float)
nz = dq != 0
if nz.any():
    c = float(np.median(ds[nz] / dq[nz]))
print(f"\n{len(rows)} legal (slot, steps) targets;  fitted c = {c:.6f}  "
      f"(recorded 0.34034)")
resid = np.abs(ds - c * dq)
print(f"max |dS - c*dQ| over all targets: {resid.max():.3e}")

free = [r for r in rows if r[4] == r[3]]
print(f"\nSAME-RUNG targets: {len(free)} of {len(rows)}  "
      f"({100*len(free)/len(rows):.0f}%)")
if free:
    fds = np.abs([r[2] for r in free])
    print(f"  max |dS| on same-rung hops: {fds.max():.3e}   mean {fds.mean():.3e}")
    print(f"  -> acceptance min(1,e^-dS) = {np.exp(-fds.max()):.12f}")
    print(f"  displacement (steps) achieved: {sorted({r[1] for r in free})[:14]}")
off = [r for r in rows if r[4] != r[3]]
if off:
    print(f"\nOFF-rung targets: {len(off)};  |dS| range "
          f"{min(abs(r[2]) for r in off):.4f} .. {max(abs(r[2]) for r in off):.4f}")
    print(f"  mean acceptance if proposed uniformly: "
          f"{np.mean([min(1.0, np.exp(-r[2])) for r in rows]):.3f}")
    print(f"  mean acceptance if proposed SAME-RUNG only: "
          f"{np.mean([min(1.0, np.exp(-r[2])) for r in free]):.3f}")
