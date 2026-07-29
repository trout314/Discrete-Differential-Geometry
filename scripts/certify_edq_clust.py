"""Certify an equilibrated ensemble at lambda=0.35 with the EDQ + volume pin +
m^2 illegal-edge CLUSTERING penalty action (no n6 legality term).

8 over-dispersed chains (4 from pristine, 4 from an over-defected start) run
until the formal gate -- quantized_split_rhat < 1.01 AND ESS >= MIN_ESS on
EVERY observable -- passes. The observable mix includes the DEFECT-COMPLEX
CENSUS (number of complexes, largest-complex size, count of large complexes),
so convergence of the slow large-complex modes is judged, not just the count.

On pass: saves the 8 chains' states as the certified ensemble + a report.
"""
import os, sys, itertools, threading, time
from collections import Counter, defaultdict

import numpy as np

_ROOT = os.path.join(os.path.dirname(__file__), "..")
sys.path.insert(0, os.path.join(_ROOT, "python"))
import discrete_differential_geometry as ddg
from discrete_differential_geometry import _dlang
from discrete_differential_geometry._sampler import ManifoldSampler, SamplerParams
from discrete_differential_geometry.convergence import (
    quantized_split_rhat, integrated_autocorrelation_time)

CELL = os.path.join(_ROOT, "data/tcp_reference/T3_R_m3_N24462.mfd")
OUT = os.path.join(_ROOT, "data/edq_clust_lam035")
ESTAR, LAM, CIMP = 5.105025, 0.35, 1.0

N_CHAINS = 8
BURNIN = 3000            # sweeps discarded before recording
REC_EVERY = 40           # sweeps between recorded samples
MIN_REC = 60             # start gating after this many samples/chain
MAX_REC = 500            # give up (report anyway) after this many
MIN_ESS = 100.0
RHAT_MAX = 1.01

# observable -> quantum (min change with other integer counts fixed; 0=continuous)
QUANTA = dict(num_facets=1.0, n_ill=1.0, n_complexes=1.0, max_size=1.0,
              n_large=1.0, objective=0.0)
OBS = list(QUANTA)

def params(N3):
    return SamplerParams(num_facets_target=N3, num_facets_coef=0.1,
                         hinge_degree_target=ESTAR, num_hinges_coef=0.0,
                         hinge_degree_variance_coef=0.0,
                         codim3_degree_variance_coef=0.0,
                         hinge_degree_target_coef=LAM * ESTAR / 6.0)

def make_sampler(m):
    s = ManifoldSampler(m, params(m.num_facets))
    s.set_n6_potential(0.0, CIMP * LAM, tilt=[0.0] * 5)   # m^2 clustering only
    return s

def observables(view, objective):
    F = np.sort(np.asarray(view.facets()), axis=1)
    deg = Counter()
    for t in F.tolist():
        for i, j in itertools.combinations(range(4), 2):
            deg[(t[i], t[j])] += 1
    ill = [e for e, c in deg.items() if c not in (0, 5, 6)]
    # defect complexes = connected components of the illegal-edge graph
    parent = {}
    def find(x):
        parent.setdefault(x, x); r = x
        while parent[r] != r: r = parent[r]
        while parent[x] != r: parent[x], x = r, parent[x]
        return r
    for a, b in ill:
        ra, rb = find(a), find(b)
        if ra != rb: parent[ra] = rb
    comps = defaultdict(set)
    for a, b in ill:
        r = find(a); comps[r].add(a); comps[r].add(b)
    sizes = sorted((len(v) for v in comps.values()), reverse=True)
    return dict(num_facets=float(view.num_facets), n_ill=float(len(ill)),
                n_complexes=float(len(sizes)),
                max_size=float(sizes[0] if sizes else 0),
                n_large=float(sum(1 for s in sizes if s >= 4)),
                objective=float(objective))

# ---------------------------------------------------------------------------
series = [defaultdict(list) for _ in range(N_CHAINS)]
samplers = [None] * N_CHAINS
lock = threading.Lock()
stop = threading.Event()

def worker(idx, seed, start_path):
    _dlang._lib.ddg_set_random_seed(seed)          # seed THIS thread's rng
    m = ddg.Manifold.load(start_path, 3)
    s = make_sampler(m)
    N3 = m.num_facets
    s.run(moves=BURNIN * N3)                        # burn-in
    while not stop.is_set():
        s.run(moves=REC_EVERY * N3)
        obs = observables(s.manifold, s.current_objective)
        with lock:
            for k, v in obs.items():
                series[idx][k].append(v)
            samplers[idx] = s

def gate():
    with lock:
        data = {k: [list(series[i][k]) for i in range(N_CHAINS)] for k in OBS}
    nrec = min(len(data["n_ill"][i]) for i in range(N_CHAINS))
    rows, passed = [], (nrec >= 4)
    for k in OBS:
        chains = [np.asarray(c, float) for c in data[k] if len(c) >= 4]
        rh = quantized_split_rhat(chains, QUANTA[k]) if len(chains) >= 2 else float("nan")
        taus = [integrated_autocorrelation_time(c) for c in chains]
        tau = float(np.nanmedian(taus)) if taus else float("nan")
        nsamp = min((len(c) for c in chains), default=0)
        ess = (len(chains) * nsamp / tau) if (tau == tau and tau > 0) else float("nan")
        halfs = [integrated_autocorrelation_time(c[:len(c)//2])
                 for c in chains if len(c) >= 16]
        tauh = float(np.nanmedian(halfs)) if halfs else float("nan")
        growth = (tau / tauh) if (tauh == tauh and tauh > 0) else float("nan")
        pooled_mean = float(np.concatenate(chains).mean()) if chains else float("nan")
        ok = (rh == rh and rh < RHAT_MAX) and (ess == ess and ess >= MIN_ESS)
        passed = passed and ok
        rows.append((k, pooled_mean, rh, ess, tau, growth, ok))
    return passed, nrec, rows

def report(nrec, rows):
    print(f"\n  gate @ {nrec} recorded samples/chain "
          f"({BURNIN}+{nrec*REC_EVERY} sweeps/chain):", flush=True)
    print(f"    {'obs':11s} {'mean':>9s} {'qRhat':>7s} {'ESS':>8s} "
          f"{'tau':>6s} {'tau_grow':>8s}  ok", flush=True)
    for k, mu, rh, ess, tau, gr, ok in rows:
        print(f"    {k:11s} {mu:9.2f} {rh:7.4f} {ess:8.0f} {tau:6.1f} "
              f"{gr:8.2f}  {'Y' if ok else '.'}", flush=True)

def main():
    os.makedirs(OUT, exist_ok=True)
    print(f"certify EDQ + m^2 clustering, lam={LAM}, m3 R N3=24462, "
          f"{N_CHAINS} chains (4 pristine + 4 over-defected)", flush=True)
    # over-defected start: melt at lam=0.30 (EDQ+m^2) briefly
    print("  building over-defected start (lam=0.30, 700 sweeps)...", flush=True)
    _dlang._lib.ddg_set_random_seed(999)
    mo = ddg.Manifold.load(CELL, 3)
    so = ManifoldSampler(mo, SamplerParams(
        num_facets_target=mo.num_facets, num_facets_coef=0.1,
        hinge_degree_target=ESTAR, num_hinges_coef=0.0,
        hinge_degree_variance_coef=0.0, codim3_degree_variance_coef=0.0,
        hinge_degree_target_coef=0.30 * ESTAR / 6.0))
    so.set_n6_potential(0.0, CIMP * 0.30, tilt=[0.0] * 5)
    so.run(moves=700 * mo.num_facets)
    over_path = os.path.join(OUT, "_overdefected_start.mfd")
    so.manifold.save(over_path)
    od = observables(so.manifold, so.current_objective)
    print(f"    over-defected start: n_ill={od['n_ill']:.0f}, "
          f"max_size={od['max_size']:.0f}", flush=True)

    starts = [CELL] * 4 + [over_path] * 4
    threads = [threading.Thread(target=worker, args=(i, 100 + i, starts[i]),
                                daemon=True) for i in range(N_CHAINS)]
    for t in threads: t.start()

    certified = False
    while True:
        time.sleep(45)
        passed, nrec, rows = gate()
        if nrec >= 4:
            report(nrec, rows)
        if nrec >= MIN_REC and passed:
            certified = True; break
        if nrec >= MAX_REC:
            break
    stop.set()
    for t in threads: t.join(timeout=120)

    print(f"\n=== {'CERTIFIED' if certified else 'NOT CERTIFIED (cap hit)'} ===",
          flush=True)
    # save the ensemble (decorrelated final states)
    for i in range(N_CHAINS):
        s = samplers[i]
        if s is None: continue
        p = os.path.join(OUT, f"sample_{i}.mfd")
        s.manifold.save(p, comments=[
            f"EDQ + m^2 clustering (imp=1.0*lam) lam={LAM} e*={ESTAR}",
            f"m3 R certified ensemble, chain {i}, seed {100+i}",
            f"start={'pristine' if starts[i]==CELL else 'over-defected'}",
            f"gate: {'CERTIFIED' if certified else 'uncertified'}"])
    print(f"  saved {N_CHAINS} samples to {OUT}", flush=True)

if __name__ == "__main__":
    main()
