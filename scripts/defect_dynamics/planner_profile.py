"""Where does plan_collapse spend its time, and how many nodes does it
touch? Determines the achievable D speedup.

Reports: nodes pushed, nodes popped, dedup hits, heap max size, and a
cProfile breakdown by function.
"""
import os, sys, time, cProfile, pstats, io, heapq
from itertools import combinations
_R = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, os.path.join(_R, "python"))
sys.path.insert(0, os.path.join(_R, "scripts", "defect_dynamics"))
sys.argv = ["t", "u", "u"]
import f0_worm as F
import discrete_differential_geometry as ddg
import link_planner as LP

V = int(os.environ.get("V", "657"))
m = ddg.Manifold.load(os.environ["START"], 3)
s, L = F.fresh(m)
S_ref = s.current_objective
LADDER = [((684, 749, 1056), (657, 687)),
          ((595, 651, 766), (657, 672)),
          ((595, 657, 684), (658, 766)),
          ((657, 749, 766), (684, 748))]

STATS = {}


def plan_instrumented(L, v, cimp, pin_coef, e_target, f30, f10,
                      f3_target, target_z=4, nodecap=200000,
                      optimize="z"):
    """verbatim copy of plan_collapse + counters"""
    s0 = LP.CollarState(L, v, cimp, pin_coef, e_target, f30, f10, f3_target)
    Z0 = len({x for f in s0.faces for x in f})
    cnt = 0
    npop = 0
    ndup = 0
    maxheap = 0

    def prio(z, ds, ln):
        return ((ds, z, ln) if optimize == "cost" else (z, ds, ln))
    heap = [(prio(Z0, 0.0, 0), cnt, s0, [])]
    seen = set()
    while heap and cnt < nodecap:
        maxheap = max(maxheap, len(heap))
        (_, _, _), _, st, ops = heapq.heappop(heap)
        npop += 1
        Z = len({x for f in st.faces for x in f})
        if Z <= target_z:
            STATS.update(pushed=cnt, popped=npop, dup=ndup,
                         maxheap=maxheap, seen=len(seen))
            return ops, st.dS
        k = st.key()
        if k in seen:
            ndup += 1
            continue
        seen.add(k)
        verts = {x for fc in st.faces for x in fc}
        for u in sorted(verts):
            if st.deg.get(LP._e(st.v, u)) == 3:
                nx = st.clone()
                nx.delete(u)
                cnt += 1
                heapq.heappush(heap, (
                    prio(len({x for fc in nx.faces for x in fc}),
                         round(nx.dS, 9), len(ops) + 1),
                    cnt, nx, ops + [("delete", u)]))
        for ab, xy in st.flip_candidates():
            for flavor in ("23", "32"):
                if flavor == "23" and xy in st.deg:
                    continue
                if flavor == "32" and st.deg.get(ab) != 3:
                    continue
                nx = st.clone()
                nx.flip(ab, xy, flavor)
                cnt += 1
                heapq.heappush(heap, (
                    prio(len({x for fc in nx.faces for x in fc}),
                         round(nx.dS, 9), len(ops) + 1),
                    cnt, nx, ops + [("flip", ab, xy, flavor)]))
    STATS.update(pushed=cnt, popped=npop, dup=ndup, maxheap=maxheap,
                 seen=len(seen))
    return None, None


def run(opt, nodecap=200000):
    fvn = [int(x) for x in s.manifold.f_vector]
    t0 = time.perf_counter()
    ops, c = plan_instrumented(L, V, F.CIMP, 1.0, F.ETARGET, fvn[3],
                               fvn[1], F.F3T, nodecap=nodecap,
                               optimize=opt)
    return ops, c, time.perf_counter() - t0


print(f"target v={V}  z={len(F.nb_of(L, V))}\n")
print("=== node counts ===")
print(f"{'k':>2} {'mode':>5} {'time_s':>8} {'pushed':>8} {'popped':>8} "
      f"{'dup':>7} {'maxheap':>8} {'us/node':>8}")
placed = []
for k in (0, 4):
    if k == 4:
        for f, a in LADDER:
            L.do(f, a)
            placed.append((f, a))
    for opt in ("z", "cost"):
        ops, c, dt = run(opt)
        us = 1e6 * dt / max(STATS['pushed'], 1)
        print(f"{k:>2} {opt:>5} {dt:>8.3f} {STATS['pushed']:>8} "
              f"{STATS['popped']:>8} {STATS['dup']:>7} "
              f"{STATS['maxheap']:>8} {us:>8.1f}")

for f, a in reversed(placed):
    L.do(a, f)
assert abs(s.current_objective - S_ref) < 1e-6

print("\n=== cProfile, cost mode, k=0 ===")
pr = cProfile.Profile()
pr.enable()
for _ in range(3):
    run("cost")
pr.disable()
sio = io.StringIO()
pstats.Stats(pr, stream=sio).sort_stats("tottime").print_stats(18)
print(sio.getvalue())
print("state restored exactly")
