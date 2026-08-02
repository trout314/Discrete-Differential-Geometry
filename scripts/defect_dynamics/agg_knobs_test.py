"""Regression + knob test for the chord channel's aggregation knobs.

(regionMax, aggBeta) = (0, 0.0) must reproduce the certified channel
bit-for-bit -- the earlier 3000-episode SEED=1 run gave 69 commits.
Then regionMax >= 1 should admit commits whose new chord has a
support-overlapping partner, which at regionMax = 0 is impossible.

n_pairs = flickers with at least one partner whose support shares a
vertex (the exact factorization-failure condition).
"""
import os, sys
from collections import Counter
_R = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, os.path.join(_R, "python"))
sys.path.insert(0, os.path.join(_R, "scripts", "defect_dynamics"))
sys.argv = ["t", "u", "u"]
import f0_worm as F
import discrete_differential_geometry as ddg
import worm_deg4_slide as W

NEP = int(os.environ.get("NEP", "3000"))
RMAX = int(os.environ.get("RMAX", "0"))
BETA = float(os.environ.get("BETA", "0.0"))
PG, PH = 0.3, 0.5

ddg.set_random_seed(int(os.environ.get("SEED", "1")))
m = ddg.Manifold.load(os.environ["START"], 3)
s, L = F.fresh(m)
s.set_worm_f0({(3, 3, 3, 3): 0.0}, [0.0] * 6, 0.0, F.Z0, lmax=F.LMAX,
              zeta=1.0, aof=0.5, ph=PH, pg=PG, bcf=1e-4,
              bc4=1.0 - PG - PH - 1e-4, maxstep=800,
              ucap_hi=50.0, ucap_lo=-50.0, mu=F.MU)
s.set_worm_pair(zeta2=float("nan"), bcp=0.05, chain_k=20)
s.set_worm_chord_agg(region_max=RMAX, agg_beta=BETA)


def npairs():
    Lv = W.Live(s.manifold)
    d3 = [e for e, d in Lv.edeg.items() if d == 3]
    sup = {}
    for e in d3:
        ts = [q for q in Lv.v2t[e[0]] if e[0] in q and e[1] in q]
        sup[e] = set(e) | ({x for q in ts for x in q} - set(e))
    return sum(1 for a in d3 if any(a != b and sup[a] & sup[b]
                                    for b in d3))


commits, census, paired = 0, Counter(), 0
for i in range(NEP):
    r = s.worm_chord_strict_episode()
    if r["changed"]:
        commits += 1
        census[r["df"]] += 1
        if commits % 50 == 0 and npairs() > 0:
            paired += 1
rx = sum(k for df, k in census.items() if df != (0, 0, 0, 0))
print(f"regionMax={RMAX} aggBeta={BETA}  episodes={NEP}")
print(f"  commits {commits}   df!=0 {rx} ({100.0*rx/max(commits,1):.0f}%)")
print(f"  n_pairs>0 at {paired}/{max(commits//50,1)} sampled checkpoints")
print(f"  final n_pairs = {npairs()}")
