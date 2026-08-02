"""DENSITY PILOT for the percolation A/B test.

From PRISTINE crystal, at each cimp, watch the illegal-edge graph
assemble and time the approach to percolation. Goal: find the cimp where
pristine -> percolated is resolvable (many samples before it happens,
but it does happen).

Observables per sample, vs ATTEMPTED MOVES:
  n_ill      total illegal edges
  S_max      largest illegal-edge component
  P = S_max/n_ill                       percolation order parameter
  N_c        number of components
  chi = sum' s^2 n_s / sum' s n_s       susceptibility (' = excl. largest);
                                        peaks AT the transition, so its
                                        peak time is a threshold-free
                                        timing signal

Reports first-passage (in attempted moves) to P >= 0.25 and P >= 0.50.
Single seed per cimp -- this is scoping, not the measurement.
"""
import os, sys, json
from collections import defaultdict
import numpy as np

_R = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, os.path.join(_R, "python"))
sys.path.insert(0, os.path.join(_R, "scripts", "defect_dynamics"))
import discrete_differential_geometry as ddg
from discrete_differential_geometry.recording import Recorder

START = os.environ["START"]
ET = float(os.environ.get("ET", "5.1"))
CHUNK = float(os.environ.get("CHUNK", "25"))
NCHUNK = int(os.environ.get("NCHUNK", "240"))
SEED = int(os.environ.get("SEED", "3"))
OUT = os.environ.get("OUT", "/tmp/perc_pilot")
CIMPS = [float(x) for x in
         os.environ.get("CIMPS", "0.35,0.30,0.25,0.20,0.15").split(",")]
STOP_P = float(os.environ.get("STOP_P", "0.60"))   # early exit
STOP_N = int(os.environ.get("STOP_N", "3"))        # consecutive samples


def build(cimp):
    m = ddg.Manifold.load(START, 3)
    p = ddg.SamplerParams(
        num_facets_target=m.num_facets, num_facets_coef=0.1,
        hinge_degree_target=ET, num_hinges_coef=1.0,
        hinge_degree_variance_coef=0.0, codim3_degree_variance_coef=0.0,
        hinge_degree_target_coef=0.0, codim3_degree_target_coef=0.0)
    s = ddg.ManifoldSampler(m, p)
    s.set_n6_potential(0.0, cimp, tilt=[0.0] * 5)
    return s


def analyse(mv):
    pairs, degs = mv.illegal_edges()
    pairs = np.asarray(pairs).reshape(-1, 2)
    degs = np.asarray(degs).reshape(-1)
    n = len(degs)
    if n == 0:
        return dict(n_ill=0, n3=0, n4=0, smax=0, ncomp=0, P=0.0, chi=0.0)
    byv = defaultdict(list)
    for i, (a, b) in enumerate(pairs):
        byv[int(a)].append(i)
        byv[int(b)].append(i)
    seen = [False] * n
    sizes = []
    for i in range(n):
        if seen[i]:
            continue
        stack, z = [i], 0
        seen[i] = True
        while stack:
            j = stack.pop()
            z += 1
            for v in (int(pairs[j][0]), int(pairs[j][1])):
                for k in byv[v]:
                    if not seen[k]:
                        seen[k] = True
                        stack.append(k)
        sizes.append(z)
    sizes.sort(reverse=True)
    smax = sizes[0]
    rest = sizes[1:]
    num = sum(z * z for z in rest)
    den = sum(rest)
    return dict(n_ill=n, n3=int((degs == 3).sum()),
                n4=int((degs == 4).sum()), smax=smax, ncomp=len(sizes),
                P=smax / n, chi=(num / den) if den else 0.0)


print(f"start {os.path.basename(START)}  e*={ET}  m2-only, edge pin ON")
print(f"chunk {CHUNK} sw, up to {NCHUNK} chunks, seed {SEED}\n")
summary = []
for cimp in CIMPS:
    ddg.set_random_seed(SEED)
    s = build(cimp)
    tag = f"cimp{cimp:g}".replace(".", "")
    rec = Recorder(s, f"{OUT}_{tag}", chunk=CHUNK,
                   extra_meta={"pilot": "percolation density",
                               "cimp": cimp, "et": ET, "seed": SEED})
    print(f"=== cimp = {cimp:g} ===")
    print(f"{'chunk':>5} {'moves':>12} {'n_ill':>7} {'n3':>6} {'n4':>7} "
          f"{'S_max':>7} {'N_c':>6} {'P':>6} {'chi':>8} {'obj':>10}")
    traj, tried, run_hi = [], 0, 0
    fp = {0.25: None, 0.50: None}
    for i in range(NCHUNK):
        row = rec.step()
        tried += row["d_tried"] or 0
        a = analyse(s.manifold)
        a.update(chunk=i + 1, moves=tried, obj=row["obj"], t=row["t"])
        traj.append(a)
        for th in fp:
            if fp[th] is None and a["P"] >= th and a["n_ill"] >= 20:
                fp[th] = tried
        if i < 6 or (i + 1) % 20 == 0:
            print(f"{i+1:>5} {tried:>12} {a['n_ill']:>7} {a['n3']:>6} "
                  f"{a['n4']:>7} {a['smax']:>7} {a['ncomp']:>6} "
                  f"{a['P']:>6.3f} {a['chi']:>8.2f} {row['obj']:>10.1f}")
        run_hi = run_hi + 1 if a["P"] >= STOP_P else 0
        if run_hi >= STOP_N:
            print(f"  -> percolated (P >= {STOP_P} for {STOP_N} samples), "
                  f"stopping at chunk {i+1}")
            break
    rec.finish()
    chi = [x["chi"] for x in traj]
    ipk = int(np.argmax(chi))
    print(f"  first passage P>=0.25: "
          f"{fp[0.25] if fp[0.25] else 'not reached':>12}")
    print(f"  first passage P>=0.50: "
          f"{fp[0.50] if fp[0.50] else 'not reached':>12}")
    print(f"  chi peak: {chi[ipk]:.2f} at chunk {traj[ipk]['chunk']} "
          f"({traj[ipk]['moves']:,} moves)   final P = {traj[-1]['P']:.3f}"
          f"   final n_ill = {traj[-1]['n_ill']}\n")
    summary.append({"cimp": cimp, "fp25": fp[0.25], "fp50": fp[0.50],
                    "chi_peak": chi[ipk], "chi_peak_moves":
                    traj[ipk]["moves"], "final_P": traj[-1]["P"],
                    "final_n_ill": traj[-1]["n_ill"],
                    "chunks": len(traj), "traj": traj})

with open(OUT + ".json", "w") as fh:
    json.dump(summary, fh, indent=1)
print("=" * 78)
print(f"{'cimp':>6} {'fp P>=.25':>12} {'fp P>=.50':>12} {'chi_pk':>9} "
      f"{'final_P':>8} {'final n_ill':>12} {'chunks':>7}")
for r in summary:
    print(f"{r['cimp']:>6.2f} "
          f"{(r['fp25'] if r['fp25'] else 0):>12} "
          f"{(r['fp50'] if r['fp50'] else 0):>12} "
          f"{r['chi_peak']:>9.2f} {r['final_P']:>8.3f} "
          f"{r['final_n_ill']:>12} {r['chunks']:>7}")
print(f"\nwrote {OUT}.json   (0 = not reached)")
