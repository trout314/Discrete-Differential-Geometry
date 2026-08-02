"""Two-sided balance test for the strict chord channel.

A balanced sampler shows NO drift where it belongs and a RESTORING drift
where it does not, so the test is: run one chain from an equilibrated
state and one from a deliberately displaced state, and compare their
late-window slopes and their up/down reaction censuses.

BOTH STARTS MUST LIE IN THE SAME f0 SECTOR. This channel conserves f0
EXACTLY -- its marks are pure flags, its walk uses only 2<->3 moves, and
its global kernel skips 1<->4 (every entry of the df census has first
component 0). An earlier run of this test compared the f0 campaign's
quenched state (f0 = 1536) against its harvested one (f0 = 1522, the
same crystal after 14 forced vertex removals). Those sit in DISJOINT
sectors, so they could never converge, and every conclusion drawn from
their non-convergence was void. Prepare the displaced start instead with
NPREP forced 2->3 flips, which conserve f0 and raise n3 and f3 directly.

Result on quench_down5q, 5000 episodes each (LO unprepared, HI +60
flips), which is the channel's balance certificate:

    LO  n3 17  -> late 17.8 +- 0.9   slope +1.37/1000ep
        census 48 up / 47 down  (ratio 1.02 -- no bias, at equilibrium)
    HI  n3 119 -> late 26.0 +- 1.7   slope -3.24/1000ep
        census 48 up / 69 down  (ratio 0.70 -- restoring bias)

The up-rate is IDENTICAL (48) in both; only the removal channel responds
to the displacement, which is what a correct restoring force looks like.
HI relaxed n3 119 -> 26 and S 4671 -> 601, closing 92% of the initial
gap, and was still descending at the end -- a residual ~8 in n3 is
incompleteness, not a barrier. Read slopes, never consecutive
checkpoints: three identical ones here were a plateau inside the noise,
not the stall they looked like.

Env: START OUT NEP SEED NPREP
"""
import json, os, sys
from collections import Counter
_R = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, os.path.join(_R, "python"))
sys.path.insert(0, os.path.join(_R, "scripts", "defect_dynamics"))
sys.argv = ["t", "u", "u"]
import f0_worm as F
import discrete_differential_geometry as ddg

START = sys.argv[1] if len(sys.argv) > 1 else None
START = os.environ["START"]
OUT = os.environ["OUT"]
NEP = int(os.environ.get("NEP", "4000"))
SEED = int(os.environ.get("SEED", "1"))
PG, PH = 0.3, 0.5

ddg.set_random_seed(SEED)
m = ddg.Manifold.load(START, 3)
s, L = F.fresh(m)

# PREP: forced 2->3 flips. They conserve f0 exactly (which this channel
# also does), so the prepared state stays in the SAME f0 sector as the
# unprepared one -- that is what makes the two-sided comparison valid.
# Each flip raises n3 and f3 by one.
NPREP = int(os.environ.get("NPREP", "0"))
if NPREP:
    done = 0
    tets = [t for ts in L.v2t.values() for t in ts]
    for k, tt in enumerate(tets):
        if done >= NPREP:
            break
        for i in range(4):
            face = tuple(sorted(x for j, x in enumerate(tt) if j != i))
            fs = set(face)
            sh = [q for q in L.v2t[face[0]] if fs <= set(q)]
            if len(sh) != 2:
                continue
            ap = tuple(sorted({x for q in sh for x in q} - fs))
            if len(ap) != 2 or ap in L.edeg:
                continue
            try:
                L.do(face, ap)
            except Exception:
                continue
            done += 1
            break
    print(f"prep: {done} forced 2->3 flips", flush=True)
s.set_worm_f0({(3, 3, 3, 3): 0.0}, [0.0] * 6, 0.0, F.Z0, lmax=F.LMAX,
              zeta=1.0, aof=0.5, ph=PH, pg=PG, bcf=1e-4,
              bc4=1.0 - PG - PH - 1e-4, maxstep=800,
              ucap_hi=50.0, ucap_lo=-50.0, mu=F.MU)
s.set_worm_pair(zeta2=float("nan"), bcp=0.05, chain_k=20)

log = open(OUT, "w")
census = Counter()
commits = 0
for i in range(NEP):
    r = s.worm_chord_strict_episode()
    if r["changed"]:
        commits += 1
        census[r["df"]] += 1
    if i % 100 == 0 or i == NEP - 1:
        fv = [int(x) for x in s.manifold.f_vector]
        n3 = sum(1 for e, d in __import__("worm_deg4_slide").Live(
            s.manifold).edeg.items() if d == 3)
        log.write(json.dumps({"ep": i, "f": fv, "n3": n3,
                              "S": round(s.current_objective, 3),
                              "commits": commits}) + "\n")
        log.flush()
        print(f"ep {i:5d}  f3={fv[3]}  n3={n3:3d}  "
              f"S={s.current_objective:9.3f}  commits={commits}",
              flush=True)
rx = sum(k for df, k in census.items() if df != (0, 0, 0, 0))
print(f"done: commits {commits}, reactions {rx} "
      f"({100.0 * rx / max(commits, 1):.0f}%)")
log.write(json.dumps({"final_census":
                      {str(k): v for k, v in census.items()}}) + "\n")
log.close()
