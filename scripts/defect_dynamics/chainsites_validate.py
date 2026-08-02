"""Validate + benchmark wf0ChainSites phase 1.

The enumerated site COUNT is the balance-relevant quantity: `pick` is
uniform over [0, n), so if n matches for every state and kmax, the
proposal distribution is unchanged. Run this on the OLD code and the NEW
code and diff the printed counts -- the melting uses only the plain
sampler at a fixed seed, so the states are bit-identical across runs.
"""
import os, sys, time, json, ctypes
import numpy as np
_R = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, os.path.join(_R, "python"))
import discrete_differential_geometry as ddg
from discrete_differential_geometry._dlang import _lib

_lib.ddg_sampler_chain_sites.argtypes = [ctypes.c_void_p, ctypes.c_int]
_lib.ddg_sampler_chain_sites.restype = ctypes.c_long

START = os.environ["START"]
CIMP = float(os.environ.get("CIMP", "0.30"))
TAG = os.environ.get("TAG", "new")
OUT = os.environ.get("OUT", "/tmp/chainsites")
ET = 5.1

ddg.set_random_seed(4242)
m = ddg.Manifold.load(START, 3)
p = ddg.SamplerParams(
    num_facets_target=m.num_facets, num_facets_coef=0.1,
    hinge_degree_target=ET, num_hinges_coef=1.0,
    hinge_degree_variance_coef=0.0, codim3_degree_variance_coef=0.0,
    hinge_degree_target_coef=0.0, codim3_degree_target_coef=0.0)
s = ddg.ManifoldSampler(m, p)
s.set_n6_potential(0.0, CIMP, tilt=[0.0] * 5)

rows = []
print(f"[{TAG}] start={os.path.basename(START)} cimp={CIMP}")
print(f"{'sweeps':>8} {'n_ill':>7} {'n3':>6} " +
      " ".join(f"{'k=' + str(k):>9}" for k in (5, 10, 20, 40)) +
      f" {'us/call':>9}")
for sw in (0, 25, 50, 100, 200, 400):
    if sw:
        s.run(sweeps=sw - rows[-1]["sw"] if rows else sw)
    pairs, degs = s.manifold.illegal_edges()
    degs = np.asarray(degs).reshape(-1)
    counts, t_tot, ncall = {}, 0.0, 0
    for k in (5, 10, 20, 40):
        t0 = time.perf_counter()
        for _ in range(5):
            c = _lib.ddg_sampler_chain_sites(s._handle, k)
        t_tot += time.perf_counter() - t0
        ncall += 5
        counts[k] = int(c)
    us = 1e6 * t_tot / ncall
    row = {"sw": sw, "n_ill": int(len(degs)),
           "n3": int((degs == 3).sum()), "counts": counts, "us": us}
    rows.append(row)
    print(f"{sw:>8} {row['n_ill']:>7} {row['n3']:>6} " +
          " ".join(f"{counts[k]:>9}" for k in (5, 10, 20, 40)) +
          f" {us:>9.1f}")

with open(f"{OUT}_{TAG}.json", "w") as fh:
    json.dump(rows, fh, indent=1)
print(f"\nwrote {OUT}_{TAG}.json")
