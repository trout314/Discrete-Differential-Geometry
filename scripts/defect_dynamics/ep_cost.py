"""How much wall time does one chord episode cost, and how much WORK
(nH+nG attempted moves) does it deliver? Determines whether an A/B at a
meaningful episode work-fraction is feasible at all.

Compares against the plain sampler's cost per attempted move at the SAME
density (cimp=0.30 mid-relaxation), which is the honest baseline -- the
plain move rate itself changes as the crystal melts.
"""
import os, sys, time
import numpy as np
_R = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, os.path.join(_R, "python"))
sys.path.insert(0, os.path.join(_R, "scripts", "defect_dynamics"))
import discrete_differential_geometry as ddg

START = os.environ["START"]
CIMP = float(os.environ.get("CIMP", "0.30"))
ET, Z0 = 5.1, 12.0
PRE = float(os.environ.get("PRE", "50"))     # sweeps of pre-melt
NEP = int(os.environ.get("NEP", "300"))
CMAXSTEP = int(os.environ.get("CMAXSTEP", "800"))

ddg.set_random_seed(5)
m = ddg.Manifold.load(START, 3)
p = ddg.SamplerParams(
    num_facets_target=m.num_facets, num_facets_coef=0.1,
    hinge_degree_target=ET, num_hinges_coef=1.0,
    hinge_degree_variance_coef=0.0, codim3_degree_variance_coef=0.0,
    hinge_degree_target_coef=0.0, codim3_degree_target_coef=0.0)
s = ddg.ManifoldSampler(m, p)
s.set_n6_potential(0.0, CIMP, tilt=[0.0] * 5)
s.set_worm_f0({(3, 3, 3, 3): 0.0}, [0.0] * 6, 0.0, Z0, lmax=4096,
              zeta=1.0, aof=0.5, ph=0.5, pg=0.3, bcf=1e-4,
              bc4=1.0 - 0.5 - 0.3 - 1e-4, maxstep=CMAXSTEP,
              ucap_hi=50.0, ucap_lo=-50.0, mu=3.0)
s.set_worm_pair(zeta2=float("nan"), bcp=0.05, chain_k=20)

print(f"cimp={CIMP}  maxstep={CMAXSTEP}   pre-melting {PRE} sweeps...")
s.run(sweeps=PRE)
pairs, degs = s.manifold.illegal_edges()
print(f"  n_ill now {len(np.asarray(degs).reshape(-1))}\n")

# --- plain sampler cost per attempted move, at THIS density
t0 = time.time()
st0 = s.get_stats().total_tried
s.run(sweeps=25)
dt_plain = time.time() - t0
dn_plain = s.get_stats().total_tried - st0
print(f"plain: {dn_plain:,} attempted moves in {dt_plain:.2f}s "
      f"= {1e6*dt_plain/dn_plain:.3f} us/move")

# --- episode cost
for rmax, beta in ((0, 0.0), (2, 2.0)):
    s.set_worm_chord_agg(region_max=rmax, agg_beta=beta)
    t0 = time.time()
    w, c, steps = 0, 0, 0
    for _ in range(NEP):
        r = s.worm_chord_strict_episode()
        w += r["nH"] + r["nG"]
        steps += r["steps"]
        c += 1 if r["changed"] else 0
    dt = time.time() - t0
    print(f"\nrmax={rmax} beta={beta}:  {NEP} episodes in {dt:.2f}s")
    print(f"  {dt/NEP*1e3:.2f} ms/episode   work {w:,} "
          f"({w/NEP:.1f} attempts/ep)   commits {c}")
    print(f"  {1e6*dt/max(w,1):.3f} us per attempted move "
          f"= {dt/max(w,1)/(dt_plain/dn_plain):.1f}x the plain sampler")
    # what it takes to hit a 50% episode work fraction
    per_chunk_plain = dn_plain
    need_ep = per_chunk_plain / max(w / NEP, 1)
    print(f"  for 50% work-fraction vs a {dn_plain:,}-move chunk: "
          f"~{need_ep:,.0f} episodes = {need_ep*dt/NEP:,.0f}s/chunk")
