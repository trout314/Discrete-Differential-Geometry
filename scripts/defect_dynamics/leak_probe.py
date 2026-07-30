"""Measure sampler RSS+swap footprint growth with the worm channel on.

Amplified: worm prob 1e-3 (10x the campaign) on the lam35r m4 gas start,
so a leak shows in ~1 minute.  Prints footprint every 10 sweeps.
"""
import os, sys, time, subprocess

_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, os.path.join(_ROOT, "python"))
import discrete_differential_geometry as ddg

WORM = float(os.environ.get("WORM", "1e-3"))
SWEEPS = int(os.environ.get("SWEEPS", "60"))
CELL = os.path.join(_ROOT, "data/mgas/lam35r_snap15000.mfd")
LAM, ESTAR, ZLEG, CIMP = 0.35, 5.105025, 0.6, 1.0


def footprint_mb():
    out = subprocess.run(
        ["ps", "-o", "rss=", "-p", str(os.getpid())],
        capture_output=True, text=True).stdout
    rss = int(out.strip()) / 1024.0
    # phys_footprint includes compressed/swapped pages
    fp = subprocess.run(["footprint", "-p", str(os.getpid())],
                        capture_output=True, text=True).stdout
    tot = ""
    for ln in fp.splitlines():
        if "phys_footprint" in ln.lower() or ln.strip().startswith("TOTAL"):
            tot = ln.strip()
            break
    return rss, tot


m = ddg.Manifold.load(CELL, 3)
N = m.num_facets
p = ddg.SamplerParams(
    num_facets_target=N, hinge_degree_target=ESTAR,
    num_facets_coef=0.1, num_hinges_coef=0.0,
    hinge_degree_variance_coef=0.0, codim3_degree_variance_coef=0.0,
    hinge_degree_target_coef=LAM * ESTAR / 6.0)
s = ddg.ManifoldSampler(m, p)
s.set_n6_potential(ZLEG * LAM, CIMP * LAM, tilt=[0.0] * 5)
s.set_worm_prob(WORM)
print(f"N={N}  worm_prob={WORM}  sweeps={SWEEPS}", flush=True)

t0 = time.time()
r0 = None
for k in range(SWEEPS):
    s.run(sweeps=1)
    if (k + 1) % 10 == 0:
        rss, tot = footprint_mb()
        if r0 is None:
            r0 = rss
        tries, acc, noc = s.worm_stats()
        used, free = ddg.gc_stats()
        print(f"sweep {k+1:4d}  t={time.time()-t0:6.1f}s  "
              f"rss={rss:8.1f}MB (d={rss-r0:+8.1f})  "
              f"Dheap used={used/1e6:8.1f}MB free={free/1e6:8.1f}MB  "
              f"worm tries={tries} acc={acc} nocand={noc}",
              flush=True)
rss, tot = footprint_mb()
print(f"FINAL rss growth: {rss - r0:+.1f} MB over {SWEEPS} sweeps", flush=True)
ddg.gc_collect()
ddg.gc_minimize()
used, free = ddg.gc_stats()
rss2, _ = footprint_mb()
print(f"after explicit GC.collect+minimize: rss={rss2:.1f}MB "
      f"(reclaimed {rss - rss2:+.1f})  Dheap used={used/1e6:.1f}MB "
      f"free={free/1e6:.1f}MB", flush=True)
