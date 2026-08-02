"""Density scan v2. Two fixes over v1:

1. COMPOUND THRESHOLD. v1 assumed a lone flicker is one size-5 component.
   Wrong: its deg-3 apex chord shares no vertex with its three deg-4 face
   edges, so on pristine crystal one 2->3 gives components {1, 3, 1}
   (chord alone | the deg-4 triangle | the deg-7). A COMPOUND is therefore
   a component of size >= 4. v1's ">5" undercounted.

2. THE lam ROWS WERE NOT THE CENSUS ACTION. reaction_census.py defaults
   are --zleg 0.6 --cimp 1.0, BOTH scaled by lam, with num_hinges_coef=0.
   v1 set zleg=cimp=0, deleting the n6 potential, which is why every lam
   row percolated. Faithful lam=0.35: zleg=0.21, cimp=0.35.

Reports deg-4 clusters specifically: components of size >= 4 holding at
least one deg-4 edge, and how many deg-4 edges live in them.
"""
import os, sys, json
from collections import Counter, defaultdict
import numpy as np

_R = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, os.path.join(_R, "python"))
sys.path.insert(0, os.path.join(_R, "scripts", "defect_dynamics"))
import discrete_differential_geometry as ddg
from discrete_differential_geometry.recording import Recorder

START = os.environ["START"]
ET = float(os.environ.get("ET", "5.1"))
BURN = float(os.environ.get("BURN", "150"))
SAMP = int(os.environ.get("NSAMP", "8"))
CHUNK = float(os.environ.get("CHUNK", "25"))
SEED = int(os.environ.get("SEED", "11"))
OUT = os.environ.get("OUT", "/tmp/deg4_scan2")

# (label, lam, zleg, cimp, hinges_coef)
GRID = json.loads(os.environ.get("GRID", json.dumps([
    ["m2 cimp0.50", 0.0, 0.0, 0.50, 1.0],
    ["m2 cimp0.45", 0.0, 0.0, 0.45, 1.0],
    ["m2 cimp0.40", 0.0, 0.0, 0.40, 1.0],
    ["m2 cimp0.35", 0.0, 0.0, 0.35, 1.0],
    ["m2 cimp0.30", 0.0, 0.0, 0.30, 1.0],
    ["census l0.40", 0.40, 0.24, 0.40, 0.0],
    ["census l0.35", 0.35, 0.21, 0.35, 0.0],
])))


def build(lam, zleg, cimp, hcoef):
    m = ddg.Manifold.load(START, 3)
    p = ddg.SamplerParams(
        num_facets_target=m.num_facets, num_facets_coef=0.1,
        hinge_degree_target=ET, num_hinges_coef=hcoef,
        hinge_degree_variance_coef=0.0, codim3_degree_variance_coef=0.0,
        hinge_degree_target_coef=lam * ET / 6.0,
        codim3_degree_target_coef=0.0)
    s = ddg.ManifoldSampler(m, p)
    s.set_n6_potential(zleg, cimp, tilt=[0.0] * 5)
    return s


def analyse(mv):
    pairs, degs = mv.illegal_edges()
    pairs = np.asarray(pairs).reshape(-1, 2)
    degs = np.asarray(degs).reshape(-1)
    n = len(degs)
    if n == 0:
        return dict(n_ill=0, n3=0, n4=0, ncomp=0, ncmpd=0, mx=0,
                    nd4c=0, d4_in=0, sizes=[])
    byv = defaultdict(list)
    for i, (a, b) in enumerate(pairs):
        byv[int(a)].append(i)
        byv[int(b)].append(i)
    seen = [False] * n
    sizes, ncmpd, nd4c, d4_in, mx = [], 0, 0, 0, 0
    for i in range(n):
        if seen[i]:
            continue
        stack, members = [i], []
        seen[i] = True
        while stack:
            j = stack.pop()
            members.append(j)
            for v in (int(pairs[j][0]), int(pairs[j][1])):
                for k in byv[v]:
                    if not seen[k]:
                        seen[k] = True
                        stack.append(k)
        z = len(members)
        sizes.append(z)
        mx = max(mx, z)
        if z >= 4:
            ncmpd += 1
            nd4 = sum(1 for j in members if degs[j] == 4)
            if nd4:
                nd4c += 1
                d4_in += nd4
    return dict(n_ill=n, n3=int((degs == 3).sum()),
                n4=int((degs == 4).sum()), ncomp=len(sizes),
                ncmpd=ncmpd, mx=mx, nd4c=nd4c, d4_in=d4_in, sizes=sizes)


print(f"start {os.path.basename(START)}  e*={ET}  burn {BURN} sw, "
      f"{SAMP}x{CHUNK} sw")
print("compound = illegal-edge component of size >= 4 "
      "(lone flicker = {1,3,1})\n")
print(f"{'config':>13} {'obj':>9} {'n_ill':>7} {'n3':>6} {'n4':>7} "
      f"{'ncomp':>6} {'cmpd':>5} {'d4cl':>5} {'d4_in':>6} {'max':>5}")
out = []
for label, lam, zleg, cimp, hcoef in GRID:
    ddg.set_random_seed(SEED)
    s = build(lam, zleg, cimp, hcoef)
    tag = label.replace(" ", "_").replace(".", "")
    rec = Recorder(s, f"{OUT}_{tag}", chunk=CHUNK, extra_meta={
        "label": label, "lam": lam, "zleg": zleg, "cimp": cimp,
        "hinges_coef": hcoef})
    rec.run(BURN)
    acc, agg = Counter(), Counter()
    objs = 0.0
    for _ in range(SAMP):
        row = rec.step()
        a = analyse(s.manifold)
        objs += row["obj"]
        for k in ("n_ill", "n3", "n4", "ncomp", "ncmpd", "nd4c", "d4_in"):
            agg[k] += a[k]
        agg["mx"] = max(agg["mx"], a["mx"])
        acc.update(a["sizes"])
    rec.finish()
    k = float(SAMP)
    print(f"{label:>13} {objs/k:>9.1f} {agg['n_ill']/k:>7.1f} "
          f"{agg['n3']/k:>6.1f} {agg['n4']/k:>7.1f} {agg['ncomp']/k:>6.1f} "
          f"{agg['ncmpd']/k:>5.1f} {agg['nd4c']/k:>5.1f} "
          f"{agg['d4_in']/k:>6.1f} {agg['mx']:>5}")
    out.append({"label": label, "lam": lam, "zleg": zleg, "cimp": cimp,
                "hinges_coef": hcoef, "obj": objs / k,
                **{x: agg[x] / k for x in
                   ("n_ill", "n3", "n4", "ncomp", "ncmpd", "nd4c", "d4_in")},
                "max": agg["mx"], "size_hist": dict(sorted(acc.items()))})

with open(OUT + ".json", "w") as fh:
    json.dump(out, fh, indent=1)
print(f"\nsize histograms (pooled over {SAMP} samples) -> {OUT}.json")
for r in out:
    h = {a: b for a, b in r["size_hist"].items() if a >= 4}
    print(f"  {r['label']:>13} size>=4: {h}")
