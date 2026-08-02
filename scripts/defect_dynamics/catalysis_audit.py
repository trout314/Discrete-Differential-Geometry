"""Does the strict chord channel cross barriers the thermal background
cannot? Structurally it must not: U == 0 and every walk move is priced
alone, so an accepted uphill move of size X costs e^-X exactly as in an
ordinary sweep. The observable test: pool the accepted UPHILL moves and
check that the largest sits near log N, i.e. that the tail is the bare
Metropolis e^-X and nothing rides above it.

Reported per arm: N accepted uphill moves, the largest accepted single-
move dS, log N for comparison, and the max running excursion
max_k (S_k - S_0) -- the barrier an episode actually surmounts.
"""
import json, math, os, sys
from collections import Counter
import numpy as np
_R = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, os.path.join(_R, "python"))
sys.path.insert(0, os.path.join(_R, "scripts", "defect_dynamics"))
sys.argv = ["t", "u", "u"]
import f0_worm as F
import discrete_differential_geometry as ddg

START = os.environ["START"]
NEP = int(os.environ.get("NEP", "3000"))
SEED = int(os.environ.get("SEED", "1"))
PG, PH = 0.3, 0.5

ddg.set_random_seed(SEED)
m = ddg.Manifold.load(START, 3)
s, L = F.fresh(m)
s.set_worm_f0({(3, 3, 3, 3): 0.0}, [0.0] * 6, 0.0, F.Z0, lmax=F.LMAX,
              zeta=1.0, aof=0.5, ph=PH, pg=PG, bcf=1e-4,
              bc4=1.0 - PG - PH - 1e-4, maxstep=800,
              ucap_hi=50.0, ucap_lo=-50.0, mu=F.MU)
s.set_worm_pair(zeta2=float("nan"), bcp=0.05, chain_k=20)

dsmax, nup, exc, census = [], 0, [], Counter()
arm = [[], []]        # per-episode max accepted dS, [head, global]
narm = [0, 0]         # accepted uphill counts, [head, global]
commits = 0
for i in range(NEP):
    r = s.worm_chord_strict_episode()
    if r["nup"]:
        dsmax.append(r["dsmax"])
        nup += r["nup"]
        exc.append(r["exc"])
    for k in (0, 1):
        if r["nuparm"][k]:
            arm[k].append(r["dsarm"][k])
            narm[k] += r["nuparm"][k]
    if r["changed"]:
        commits += 1
        census[r["df"]] += 1

a = np.array(dsmax) if dsmax else np.zeros(1)
e = np.array(exc) if exc else np.zeros(1)
print(f"episodes {NEP}, commits {commits}")
print(f"accepted UPHILL moves pooled: N = {nup}")
print(f"  log N                       = {math.log(max(nup,1)):8.2f}")
print(f"  largest accepted dS         = {a.max():8.2f}")
print(f"  99th pct of per-ep max dS   = {np.percentile(a, 99):8.2f}")
print(f"  median per-episode max dS   = {np.median(a):8.2f}")
print(f"max running excursion max_k(S_k - S_0):")
print(f"  largest                     = {e.max():8.2f}")
print(f"  median                      = {np.median(e):8.2f}")
print("\nARM SPLIT  (global = plain thermal, forbidden to touch marks)")
print(f"{'arm':>8} {'N up':>7} {'max dS':>8} {'p99':>7} {'median':>7}")
for k, nm in ((0, "head"), (1, "global")):
    if arm[k]:
        v = np.array(arm[k])
        print(f"{nm:>8} {narm[k]:>7} {v.max():>8.2f} "
              f"{np.percentile(v, 99):>7.2f} {np.median(v):>7.2f}")
    else:
        print(f"{nm:>8} {narm[k]:>7}    (none)")
rx = sum(k for df, k in census.items() if df != (0, 0, 0, 0))
print(f"census: {commits} commits, {rx} with df != 0 "
      f"({100.0*rx/max(commits,1):.0f}%)")
json.dump({"nup": nup, "dsmax": list(map(float, a)),
           "exc": list(map(float, e)), "commits": commits},
          open(os.environ.get("OUT", "/dev/null"), "w"))
