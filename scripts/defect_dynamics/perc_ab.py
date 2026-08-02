"""Percolation A/B driver, with HONEST work accounting.

THE CONTROL: ddg_capi.d increments s.totalTried in exactly two places,
both inside the main run() loop. The worm episodes NEVER touch it. So
Recorder's d_tried counts arm A's moves only, and an arm-B run would get
its episode-internal moves for free -- any "speedup" would be an
artifact. Arm B's true work is

    work = d_tried  +  sum over episodes of (nH + nG)

where nH/nG are incremented once per PROPOSED move in the head and
global arms of wormChordStrictEpisode (sampler.d: res.nH++ / res.nG++
inside the walk loop). Both are exposed in the episode's return dict.
Every first-passage below is reported in these units.

Arms:
  thermal            plain sampling only
  chord[:rmax:beta]  plain sampling + CHORD episodes interleaved per chunk

Observables (illegal-edge graph): P = S_max/n_ill, N_c, and
chi = sum' s^2 n_s / sum' s n_s (' = excluding largest), which peaks AT
the transition.
"""
import os, sys, json, time
from collections import defaultdict
import numpy as np

_R = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, os.path.join(_R, "python"))
sys.path.insert(0, os.path.join(_R, "scripts", "defect_dynamics"))
import discrete_differential_geometry as ddg
from discrete_differential_geometry.recording import Recorder

START = os.environ["START"]
ET = float(os.environ.get("ET", "5.1"))
CIMP = float(os.environ.get("CIMP", "0.30"))
CHUNK = float(os.environ.get("CHUNK", "25"))
NCHUNK = int(os.environ.get("NCHUNK", "200"))
ARM = os.environ.get("ARM", "thermal")
SEEDS = [int(x) for x in os.environ.get("SEEDS", "1").split(",")]
NEP = int(os.environ.get("NEP", "200"))        # chord episodes per chunk
RMAX = int(os.environ.get("RMAX", "0"))
BETA = float(os.environ.get("BETA", "0.0"))
OUT = os.environ.get("OUT", "/tmp/perc_ab")
STOP_P = float(os.environ.get("STOP_P", "0.60"))
STOP_N = int(os.environ.get("STOP_N", "3"))
THRESH = [0.25, 0.50]

# chord-channel config (mirrors f0_worm's push_chord defaults)
CPH = float(os.environ.get("CPH", "0.5"))
CPG = float(os.environ.get("CPG", "0.3"))
CMAXSTEP = int(os.environ.get("CMAXSTEP", "800"))
CBCP = float(os.environ.get("CBCP", "0.05"))
CHAINK = int(os.environ.get("CHAINK", "20"))
MU = float(os.environ.get("MU", "3.0"))
LMAX = int(os.environ.get("LMAX", "4096"))
Z0 = 12.0


def build(seed):
    ddg.set_random_seed(seed)
    m = ddg.Manifold.load(START, 3)
    p = ddg.SamplerParams(
        num_facets_target=m.num_facets, num_facets_coef=0.1,
        hinge_degree_target=ET, num_hinges_coef=1.0,
        hinge_degree_variance_coef=0.0, codim3_degree_variance_coef=0.0,
        hinge_degree_target_coef=0.0, codim3_degree_target_coef=0.0)
    s = ddg.ManifoldSampler(m, p)
    s.set_n6_potential(0.0, CIMP, tilt=[0.0] * 5)
    if ARM.startswith("chord"):
        # U == 0 for the strict channel, so the tube is irrelevant; this
        # is the minimal valid upload.
        s.set_worm_f0({(3, 3, 3, 3): 0.0}, [0.0] * 6, 0.0, Z0, lmax=LMAX,
                      zeta=1.0, aof=0.5, ph=CPH, pg=CPG, bcf=1e-4,
                      bc4=1.0 - CPH - CPG - 1e-4, maxstep=CMAXSTEP,
                      ucap_hi=50.0, ucap_lo=-50.0, mu=MU)
        s.set_worm_pair(zeta2=float("nan"), bcp=CBCP, chain_k=CHAINK)
        s.set_worm_chord_agg(region_max=RMAX, agg_beta=BETA)
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
    rest = sizes[1:]
    den = sum(rest)
    return dict(n_ill=n, n3=int((degs == 3).sum()),
                n4=int((degs == 4).sum()), smax=sizes[0],
                ncomp=len(sizes), P=sizes[0] / n,
                chi=(sum(z * z for z in rest) / den) if den else 0.0)


def run_one(seed):
    s = build(seed)
    tag = f"{ARM.replace(':', '')}_s{seed}"
    rec = Recorder(s, f"{OUT}_{tag}", chunk=CHUNK, extra_meta={
        "ab": "percolation", "arm": ARM, "cimp": CIMP, "seed": seed,
        "nep": NEP if ARM.startswith("chord") else 0,
        "rmax": RMAX, "beta": BETA})
    work_plain, work_ep, commits, eps = 0, 0, 0, 0
    fp = {t: None for t in THRESH}
    fp_chi, chi_best = None, -1.0
    traj, run_hi, t0 = [], 0, time.time()
    for i in range(NCHUNK):
        row = rec.step()
        work_plain += row["d_tried"] or 0
        if ARM.startswith("chord"):
            for _ in range(NEP):
                r = s.worm_chord_strict_episode()
                work_ep += r["nH"] + r["nG"]
                eps += 1
                if r["changed"]:
                    commits += 1
        work = work_plain + work_ep
        a = analyse(s.manifold)
        a.update(chunk=i + 1, work=work, work_plain=work_plain,
                 work_ep=work_ep, obj=row["obj"], commits=commits)
        traj.append(a)
        for t in THRESH:
            if fp[t] is None and a["P"] >= t and a["n_ill"] >= 20:
                fp[t] = work
        if a["chi"] > chi_best:
            chi_best, fp_chi = a["chi"], work
        run_hi = run_hi + 1 if a["P"] >= STOP_P else 0
        if run_hi >= STOP_N:
            break
    rec.finish()
    epfrac = work_ep / max(work_plain + work_ep, 1)
    print(f"  seed {seed:>4}: chunks {len(traj):>4}  "
          f"work {work:>12,}  ep_frac {epfrac:>5.1%}  "
          f"commits {commits:>5}/{eps:<6}  "
          f"fp25 {fp[0.25] or 0:>11,}  fp50 {fp[0.50] or 0:>11,}  "
          f"chi_pk {chi_best:>7.2f} @ {fp_chi:>11,}  "
          f"({time.time()-t0:.0f}s)")
    return {"seed": seed, "arm": ARM, "chunks": len(traj), "work": work,
            "work_plain": work_plain, "work_ep": work_ep,
            "ep_frac": epfrac, "commits": commits, "episodes": eps,
            "fp25": fp[0.25], "fp50": fp[0.50], "chi_peak": chi_best,
            "chi_peak_work": fp_chi, "final_P": traj[-1]["P"],
            "traj": traj}


print(f"ARM={ARM}  cimp={CIMP}  e*={ET}  chunk={CHUNK}sw  "
      f"nchunk<={NCHUNK}  seeds={SEEDS}")
if ARM.startswith("chord"):
    print(f"  chord: {NEP} episodes/chunk  rmax={RMAX} beta={BETA} "
          f"ph={CPH} pg={CPG} maxstep={CMAXSTEP}")
print(f"  work = d_tried + sum(nH+nG)   [episodes do NOT touch "
      f"totalTried]\n")
res = [run_one(sd) for sd in SEEDS]
with open(OUT + f"_{ARM.replace(':', '')}.json", "w") as fh:
    json.dump(res, fh, indent=1)


def stats(key):
    v = np.array([r[key] for r in res if r[key]], float)
    if len(v) < 2:
        return None
    lv = np.log(v)
    return dict(n=len(v), mean=v.mean(), sd=v.std(ddof=1),
                cv=v.std(ddof=1) / v.mean(), med=np.median(v),
                gmean=np.exp(lv.mean()), sd_log=lv.std(ddof=1),
                sem_log=lv.std(ddof=1) / np.sqrt(len(v)))


print(f"\n{'observable':>12} {'n':>3} {'mean':>13} {'sd':>12} {'CV':>6} "
      f"{'gmean':>13} {'sd(log)':>8} {'sem(log)':>9}")
for k in ("fp25", "fp50", "chi_peak_work"):
    st = stats(k)
    if st:
        print(f"{k:>12} {st['n']:>3} {st['mean']:>13,.0f} "
              f"{st['sd']:>12,.0f} {st['cv']:>6.2f} {st['gmean']:>13,.0f} "
              f"{st['sd_log']:>8.3f} {st['sem_log']:>9.3f}")
        # chains needed to resolve a given ratio at 2 sigma, per arm
        for ratio in (1.2, 1.5, 2.0):
            need = 2 * (2 * st['sd_log'] / np.log(ratio)) ** 2
            print(f"{'':>12}   to resolve a {ratio:.1f}x speedup at 2sd: "
                  f"~{int(np.ceil(need)):>3} chains PER ARM")
print(f"\nwrote {OUT}_{ARM.replace(':', '')}.json")
