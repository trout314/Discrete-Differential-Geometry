#!/usr/bin/env python3
"""Fixed-coefficient VDV equilibrium sampling with two-sided convergence checks.

Goal: find, for a given N, the largest VDV coefficient beta at which we can still
*equilibrate* (the "edge of equilibratability"), and measure the equilibrium
vertex-degree variance <VDV>_beta there -- so we can compare it against what
(non-equilibrium) annealing achieves.

Method. At each beta we run chains from BOTH sides of the equilibrium VDV:
  * "above": start from a plain high-VDV seed and hold at beta -> VDV falls.
  * "below": start from the same seed, over-squeeze at a higher warmup
             coefficient, then hold at beta -> VDV rises.
If the two sides converge to a common value (split-R-hat -> 1, small gap) and the
integrated autocorrelation time tau is small enough for a usable ESS, beta is
equilibratable and <VDV>_beta is bracketed. Every other objective parameter
(facet-count pin, edge-degree pin below the flat value, HDV) is held fixed, so
only VDV is being equilibrated.

Two modes:
  --chain   run ONE chain, write its per-sample time series to --out (a worker).
  (default) driver: spawn chains over a beta-grid x {above,below} x replicas,
            then report R-hat / tau / ESS / <VDV> and the two-sided gap per beta.

Example
-------
    python scripts/equilibrium_vdv.py --n-target 1000 \
        --seed-file seeds/S3_N1e3_1e-1_ED5p0043_1e-1_s000.mfd \
        --coefs 1 5 20 50 100 200 400 --replicas 2 \
        --thin 5 --n-samples 500 --burnin 100 --output-dir data/equil_vdv/N1e3
"""

import argparse
import csv
import glob
import json
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "python"))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "tools"))
import discrete_differential_geometry as ddg
from discrete_differential_geometry import (
    Manifold, ManifoldSampler, SamplerParams,
    integrated_autocorrelation_time, split_rhat, rank_normalized_rhat,
    quantized_split_rhat,
)

# Coupled observables to check for joint convergence: the free/penalized ones
# (VDV, HDV) and the pinned ones (edge degree, facet count). Equilibrium requires
# ALL to converge -- a drift in a pinned coordinate masquerades as (non-)
# convergence in a penalized one, since they are coupled.
OBSERVABLES = ["vdv", "edge_deg", "hdv", "num_facets"]
from seed_utils import (build_metadata_comments, build_seed_filename,
                        get_total_memory_gb)

FIELDS = ["side", "replica", "coef", "sample", "sweeps",
          "num_facets", "vdv", "edge_deg", "hdv", "accept"]

# Per-chain resident memory ~ manifold footprint: ~2 MB per 1000 facets + base
# (rounded up from the ~2 KB/facet measured on reburn). At N=1e7 this is ~20 GB,
# so one chain per ~24 GB of RAM; the driver caps concurrency accordingly.
_MB_PER_1K_FACETS = 2.5
_BASE_MB = 250.0


def est_chain_gb(n_target):
    return (_BASE_MB + _MB_PER_1K_FACETS * (n_target / 1000)) / 1024


def effective_workers(args):
    """Cap concurrent chains so total estimated RAM fits the budget."""
    budget = args.max_memory_gb if args.max_memory_gb > 0 else max(
        1.0, get_total_memory_gb() - 4.0)
    per = est_chain_gb(args.n_target)
    cap = max(1, int(budget / per))
    return max(1, min(args.max_workers, cap)), per, budget


# ----------------------------------------------------------------------------
# Worker: one fixed-beta chain
# ----------------------------------------------------------------------------

def _write_rows(path, rows):
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    tmp = path + ".tmp"
    with open(tmp, "w", newline="") as f:
        w = csv.writer(f); w.writerow(FIELDS); w.writerows(rows)
    os.replace(tmp, path)


def _read_rows(path, n):
    if not path or not os.path.exists(path):
        return []
    with open(path) as f:
        r = csv.reader(f); next(r, None)
        rows = list(r)
    return rows[:n]


def _hist_path(out):
    return (out + ".hist.npz") if out else None


def _pad_stack(rows):
    """Stack ragged 1-D count arrays into a zero-padded int64 2-D array."""
    w = max((len(a) for a in rows), default=0)
    out = np.zeros((len(rows), w), dtype=np.int64)
    for i, a in enumerate(rows):
        out[i, :len(a)] = a
    return out


def _write_hists(path, vh, eh):
    """Atomically persist per-sample vertex/edge degree histograms (index i =
    degree i+1, matching the D accessor's 1-indexed convention)."""
    if not path:
        return
    tmp = path + ".tmp"
    with open(tmp, "wb") as f:
        np.savez(f, vhist=_pad_stack(vh), ehist=_pad_stack(eh))
    os.replace(tmp, path)


def _read_hists(path, n):
    """Reload the first n samples of each histogram as ragged lists (for resume)."""
    if not path or not os.path.exists(path):
        return [], []
    with np.load(path) as z:
        return [r for r in z["vhist"][:n]], [r for r in z["ehist"][:n]]


def run_chain(args):
    """One fixed-beta chain, resumable. Phases: warmup -> burnin -> measure.
    Resume is exact-enough for equilibrium: on restart we reload the checkpointed
    manifold and continue (a fresh RNG is fine -- the chain is Markov, and post-
    checkpoint sweeps only add relaxation). The JSON is written last, so it is the
    commit point; if the manifold on disk is slightly ahead of it, resume just
    re-runs a few harmless sweeps."""
    ck_mfd = ck_json = None
    do_ckpt = bool(args.save_config) and args.checkpoint_every > 0
    if do_ckpt:
        ck_mfd, ck_json = args.save_config + ".ckpt.mfd", args.save_config + ".ckpt.json"

    if args.save_config and os.path.exists(args.save_config):
        return  # already finished

    st = None
    if do_ckpt and os.path.exists(ck_json) and os.path.exists(ck_mfd):
        try:
            st = json.load(open(ck_json))
        except (OSError, ValueError):
            st = None

    rec_hist = bool(args.record_histograms) and bool(args.out)
    hp = _hist_path(args.out) if rec_hist else None
    if st:
        mfd = Manifold.load(ck_mfd, args.dim)
        phase, ph = st["phase"], st["ph_sweeps"]
        samples, meas_total = st["samples"], st["meas_total"]
        rows = _read_rows(args.out, samples)
        vhists, ehists = _read_hists(hp, samples) if rec_hist else ([], [])
    else:
        mfd = Manifold.load(args.init, args.dim)
        phase = "warmup" if args.warmup_sweeps > 0 else "burnin"
        ph = samples = meas_total = 0
        rows = []
        vhists, ehists = [], []

    params = SamplerParams(
        num_facets_target=args.n_target, num_facets_coef=args.num_facets_coef,
        hinge_degree_target=args.hinge_target, num_hinges_coef=args.num_hinges_coef,
        hinge_degree_variance_coef=args.hdv_coef,
        codim3_degree_variance_coef=args.coef)
    s = ManifoldSampler(mfd, params)
    v = s.manifold

    def checkpoint():
        if not do_ckpt:
            return
        tmp = ck_mfd + ".tmp"; v.save(tmp); os.replace(tmp, ck_mfd)   # manifold first
        if args.out:
            _write_rows(args.out, rows)
        if rec_hist:
            _write_hists(hp, vhists, ehists)
        tmp = ck_json + ".tmp"                                         # JSON last = commit
        with open(tmp, "w") as f:
            json.dump({"phase": phase, "ph_sweeps": ph, "samples": samples,
                       "meas_total": meas_total}, f)
        os.replace(tmp, ck_json)

    ce = args.checkpoint_every if do_ckpt else 10 ** 12  # chunk size (sweeps)

    if phase == "warmup":
        s.set_codim3_degree_variance_coef(args.warmup_coef)
        if args.warmup_hinge_coef > 0:      # hold edge degree while squeezing VDV
            s.set_num_hinges_coef(args.warmup_hinge_coef)
        while ph < args.warmup_sweeps:
            step = min(ce, args.warmup_sweeps - ph)
            s.run(sweeps=step); ph += step; checkpoint()
        s.set_codim3_degree_variance_coef(args.coef)
        if args.warmup_hinge_coef > 0:
            s.set_num_hinges_coef(args.num_hinges_coef)
        phase, ph = "burnin", 0; checkpoint()

    if phase == "burnin":
        while ph < args.burnin_sweeps:
            step = min(ce, args.burnin_sweeps - ph)
            s.run(sweeps=step); ph += step; checkpoint()
        phase, ph = "measure", 0; checkpoint()

    if phase == "measure":
        ck_samples = max(1, ce // max(args.thin, 1))
        while samples < args.n_samples:
            s.reset_stats(); s.run(sweeps=args.thin); meas_total += args.thin
            stt = s.get_stats()
            acc = stt.total_accepted / stt.total_tried if stt.total_tried else 0.0
            rows.append([args.side, args.replica, args.coef, samples, meas_total,
                         v.num_facets, v.degree_variance(0), v.mean_degree(1),
                         v.degree_variance(1), acc])
            if rec_hist:
                vhists.append(v.degree_histogram(0))
                ehists.append(v.degree_histogram(1))
            samples += 1
            if samples % ck_samples == 0:
                checkpoint()

    if args.out:
        _write_rows(args.out, rows)
    if rec_hist:
        _write_hists(hp, vhists, ehists)
    if args.save_config:
        comments = build_metadata_comments(
            topology=args.topology, dimension=args.dim,
            initial_triangulation=os.path.basename(args.init),
            num_facets_target=args.n_target, num_facets_coef=args.num_facets_coef,
            hinge_degree_target=args.hinge_target, num_hinges_coef=args.num_hinges_coef,
            hinge_degree_variance_coef=args.hdv_coef,
            codim3_degree_variance_coef=args.coef,
            growth_step_size=0, eq_sweeps_per_step=0,
            equilibration_sweeps=args.warmup_sweeps + args.burnin_sweeps
                                 + args.n_samples * args.thin,
            manifold_view=v, objective=s.current_objective,
            sampler_stats=s.get_stats())
        os.makedirs(os.path.dirname(args.save_config) or ".", exist_ok=True)
        v.save(args.save_config, comments=comments)
    if do_ckpt:                                     # clean up on success
        for p in (ck_mfd, ck_json):
            try:
                os.remove(p)
            except OSError:
                pass


# ----------------------------------------------------------------------------
# Driver: grid + diagnostics
# ----------------------------------------------------------------------------

def spawn(args, coef, side, rep, out):
    # "below" start: either a supplied low-VDV config (e.g. an annealed one) held
    # at beta so VDV rises to equilibrium, or the plain seed over-squeezed by a
    # warmup at a higher coefficient. "above" always starts from the plain seed.
    if side == "below" and args.below_init:
        init, warm_sweeps = args.below_init, 0
    elif side == "below":
        init, warm_sweeps = args.seed_file, args.warmup_sweeps
    else:
        init, warm_sweeps = args.seed_file, 0
    cmd = [sys.executable, os.path.abspath(__file__), "--chain",
           "--init", init, "--dim", str(args.dim),
           "--n-target", str(args.n_target),
           "--num-facets-coef", str(args.num_facets_coef),
           "--hinge-target", str(args.hinge_target),
           "--num-hinges-coef", str(args.num_hinges_coef),
           "--hdv-coef", str(args.hdv_coef),
           "--coef", str(coef),
           "--warmup-coef", str(coef * args.warmup_factor),
           "--warmup-sweeps", str(warm_sweeps),
           "--warmup-hinge-coef", str(args.warmup_hinge_coef),
           "--thin", str(args.thin), "--n-samples", str(args.n_samples),
           "--side", side, "--replica", str(rep), "--out", out]
    return subprocess.run(cmd).returncode


def load_series(path, burnin):
    """Return the post-burnin (vdv, edge_deg, hdv, accept) arrays from a chain CSV."""
    v, e, h, a = [], [], [], []
    with open(path) as f:
        for row in csv.DictReader(f):
            if int(row["sample"]) >= burnin:
                v.append(float(row["vdv"])); e.append(float(row["edge_deg"]))
                h.append(float(row["hdv"])); a.append(float(row["accept"]))
    return np.array(v), np.array(e), np.array(h), np.array(a)


def load_observables(path, burnin=0):
    """Return {observable -> post-burnin array} for OBSERVABLES + accept."""
    cols = {k: [] for k in OBSERVABLES + ["accept"]}
    with open(path) as f:
        for row in csv.DictReader(f):
            if int(row["sample"]) >= burnin:
                for k in cols:
                    cols[k].append(float(row[k]))
    return {k: np.array(v) for k, v in cols.items()}


def load_degree_fracs(out_path, dim, burnin=0):
    """One chain's degree distribution over time: return (fracs, n_simp) where
    fracs is an (n_samples x bins) array of per-degree fractions p_d (bin d ->
    degree d+1, the D accessor's 1-indexed convention) and n_simp is the mean
    number of dimension-`dim` simplices (the denominator whose reciprocal is the
    p_d discretization quantum). Returns (None, None) if the sidecar is absent."""
    hp = _hist_path(out_path)
    if not hp or not os.path.exists(hp):
        return None, None
    with np.load(hp) as z:
        counts = z["vhist" if dim == 0 else "ehist"].astype(float)
    if counts.size == 0:
        return None, None
    counts = counts[burnin:]
    totals = counts.sum(axis=1)
    denom = np.where(totals == 0, 1.0, totals)
    return counts / denom[:, None], float(totals.mean())


def analyze(coef, out_dir, burnin, thin):
    chains = {"above": [], "below": []}
    echains = []
    edge_all, acc_all = [], []
    for path in sorted(glob.glob(os.path.join(out_dir, f"coef{coef:g}_*.csv"))):
        side = "above" if "_above_" in path else "below"
        v, e, h, a = load_series(path, burnin)
        if len(v) >= 4:
            chains[side].append(v); echains.append(e)
            edge_all.append(e.mean()); acc_all.append(a.mean())
    allv = chains["above"] + chains["below"]
    if not allv:
        return None
    rhat = rank_normalized_rhat(allv) if len(allv) >= 2 else float("nan")
    rhat_e = rank_normalized_rhat(echains) if len(echains) >= 2 else float("nan")
    emeans = np.array([e.mean() for e in echains]) if echains else np.array([0.0])
    edge_spread = float(emeans.max() - emeans.min())
    taus = [integrated_autocorrelation_time(v) for v in allv]
    tau = float(np.nanmedian(taus))
    n_post = sum(len(v) for v in allv)
    ess = n_post / tau if tau and not np.isnan(tau) else float("nan")
    mean_above = np.mean([v.mean() for v in chains["above"]]) if chains["above"] else np.nan
    mean_below = np.mean([v.mean() for v in chains["below"]]) if chains["below"] else np.nan
    vdv_all = np.concatenate(allv)
    return dict(coef=coef, rhat=rhat, rhat_e=rhat_e, edge_spread=edge_spread,
                tau_samples=tau, tau_sweeps=tau * thin, ess=ess, vdv=vdv_all.mean(),
                vdv_std=vdv_all.std(), above=mean_above, below=mean_below,
                gap=abs(mean_above - mean_below), edge=np.mean(edge_all),
                accept=np.mean(acc_all), n_chains=len(allv))


def produce(args):
    """Run K equilibrium chains at args.beta; if split-R-hat passes, copy the
    configs into seeds/ with library names. Half the chains start from the plain
    seed (above) and half from --below-init (below) so R-hat tests convergence
    from both directions."""
    import shutil
    if not (args.seed_file and args.beta):
        sys.exit("--produce needs --seed-file and --beta")
    params = SamplerParams(
        num_facets_target=args.n_target, num_facets_coef=args.num_facets_coef,
        hinge_degree_target=args.hinge_target, num_hinges_coef=args.num_hinges_coef,
        hinge_degree_variance_coef=args.hdv_coef,
        codim3_degree_variance_coef=args.beta)
    stem = build_seed_filename(args.topology, params, seed_index=0).rsplit("_s", 1)[0]
    if glob.glob(os.path.join(args.seeds_dir, stem + "_s*.mfd")):
        print(f"seeds already exist for these params (stem {stem}); not duplicating.",
              file=sys.stderr)
        return

    K = args.replicas
    stage = os.path.join(args.output_dir, "staging")
    os.makedirs(stage, exist_ok=True)
    print(f"Producing {K} equilibrium chains: N={args.n_target}, beta={args.beta:g} "
          f"(beta/N={args.beta/args.n_target:.4g}), burnin={args.production_burnin} + "
          f"{args.n_samples}x{args.thin} measured", flush=True)

    def pspawn(i):
        below = (i % 2 == 1) and bool(args.below_init)
        init = args.below_init if below else args.seed_file
        cfg = os.path.join(stage, f"chain_{i}.mfd")
        out = os.path.join(stage, f"chain_{i}.csv")
        cmd = [sys.executable, os.path.abspath(__file__), "--chain",
               "--init", init, "--dim", str(args.dim), "--n-target", str(args.n_target),
               "--num-facets-coef", str(args.num_facets_coef),
               "--hinge-target", str(args.hinge_target),
               "--num-hinges-coef", str(args.num_hinges_coef),
               "--hdv-coef", str(args.hdv_coef), "--coef", str(args.beta),
               "--burnin-sweeps", str(args.production_burnin),
               "--thin", str(args.thin), "--n-samples", str(args.n_samples),
               "--side", "below" if below else "above", "--replica", str(i),
               "--topology", args.topology, "--out", out, "--save-config", cfg]
        if args.record_histograms:
            cmd.append("--record-histograms")
        return i, subprocess.run(cmd).returncode, cfg, out

    nw, per, budget = effective_workers(args)
    print(f"  concurrency: {nw} chains (~{per:.2f} GB/chain, {budget:.0f} GB budget)",
          flush=True)
    res = {}
    with ThreadPoolExecutor(max_workers=nw) as pool:
        for fut in as_completed([pool.submit(pspawn, i) for i in range(K)]):
            i, rc, cfg, out = fut.result()
            res[i] = (rc, cfg, out)
            if rc != 0:
                print(f"  chain {i} FAILED rc={rc}", file=sys.stderr)

    # Joint convergence over all coupled observables, judged by ONE statistic of
    # any observable's nature: quantization-floored split-R-hat. The within-chain
    # variance is floored at Delta^2/12, Delta being the minimum change in the
    # observable with the other underlying integer counts held fixed. For a
    # genuinely-fluctuating coordinate (VDV, HDV) the true variance dominates and
    # this is ordinary split-R-hat; for a near-deterministic one (edge degree =
    # 6 f3/f1, facet count) the floor sets the scale, so chains that agree to
    # within a few quanta pass instead of the statistic blowing up on a vanishing
    # denominator. It is also autocorrelation-aware (B): run length is judged
    # against the slowest coordinate's ESS.
    RH = args.rhat_max
    obs_series = {o: [] for o in OBSERVABLES}
    for i in range(K):
        rc, cfg, out = res.get(i, (1, None, None))
        if rc == 0 and out and os.path.exists(out):
            cols = load_observables(out)
            if len(cols["vdv"]) >= 4:
                for o in OBSERVABLES:
                    obs_series[o].append(cols[o])
    n_ok = len(obs_series["vdv"])
    n_samp = min((len(s) for s in obs_series["vdv"]), default=0)

    def observable_quantum(o, pooled):
        if o == "edge_deg":                       # 6 f3/f1; step f1 at fixed f3
            f3 = np.concatenate(obs_series["num_facets"]).mean()
            return (pooled.mean() ** 2) / (6.0 * f3) if f3 else 0.0
        if o == "num_facets":                     # a 2-3 move changes f3 by 1
            return 1.0
        return 0.0                                # variances: effectively continuous

    print(f"\nEnsemble: {n_ok} chains x {n_samp} samples", flush=True)
    print(f"  {'observable':>12} {'mean':>12} {'quantum':>10} {'qRhat':>8} "
          f"{'tau(smp)':>9} {'ESS':>7} {'gate':>14}")
    diag, bad = {}, []
    for o in OBSERVABLES:
        ser = obs_series[o]
        pooled = np.concatenate(ser) if ser else np.array([0.0])
        q = observable_quantum(o, pooled)
        rh = quantized_split_rhat(ser, q) if n_ok >= 2 else float("nan")
        taus = [integrated_autocorrelation_time(s) for s in ser]
        tau = float(np.nanmedian(taus)) if taus else float("nan")
        ess = (n_ok * n_samp / tau) if tau and not np.isnan(tau) else float("nan")
        ok = not np.isnan(rh) and rh < RH
        if not ok:
            bad.append(o)
        diag[o] = dict(rhat=rh, tau=tau, ess=ess, mean=pooled.mean())
        print(f"  {o:>12} {pooled.mean():>12.4f} {q:>10.5f} {rh:>8.3f} "
              f"{tau:>9.1f} {ess:>7.0f} {'qRhat<'+format(RH,'g'):>9}:"
              f"{'OK' if ok else 'FAIL':>4}")

    # Richer check: the full vertex- and edge-degree DISTRIBUTIONS, not just their
    # first two moments. Each degree bin's occupancy p_d = count_d / n_simplices is
    # another observable (quantum 1/n_simplices, count_d integer). "Distributions
    # statistically indistinguishable" == every material bin converges. Bins with
    # too little pooled mass are noise and would false-fail, so only bins whose mean
    # count >= --bin-mass-floor are gated; the sparse tail is reported, not gated.
    hist_bad = []
    if args.record_histograms:
        for dim, label in ((0, "vtx-dist"), (1, "edge-dist")):
            per_chain, nsimps = [], []
            for i in range(K):
                rc, cfg, out = res.get(i, (1, None, None))
                if rc == 0 and out:
                    fr, ns = load_degree_fracs(out, dim)
                    if fr is not None and len(fr) >= 4:
                        per_chain.append(fr); nsimps.append(ns)
            if len(per_chain) < 2:
                continue
            width = max(fr.shape[1] for fr in per_chain)
            per_chain = [np.pad(fr, ((0, 0), (0, width - fr.shape[1])))
                         for fr in per_chain]
            nsimp = float(np.mean(nsimps))
            q = 1.0 / nsimp if nsimp else 0.0
            mean_count = np.concatenate(per_chain, axis=0).mean(axis=0) * nsimp
            material = [d for d in range(width)
                        if mean_count[d] >= args.bin_mass_floor]
            worst, nfail = 0.0, 0
            for d in material:
                rh = quantized_split_rhat([fr[:, d] for fr in per_chain], q)
                if not np.isnan(rh):
                    worst = max(worst, rh)
                    nfail += rh >= RH
            if material:
                print(f"  {label:>12} {len(material):>3d} bins (deg "
                      f"{material[0] + 1}-{material[-1] + 1}) worst qRhat="
                      f"{worst:.3f}  {nfail} fail")
            if nfail:
                hist_bad.append(f"{label}({nfail})")

    min_ess = min((diag[o]["ess"] for o in OBSERVABLES
                   if not np.isnan(diag[o]["ess"])), default=float("nan"))
    all_bad = bad + hist_bad
    if all_bad or (not np.isnan(min_ess) and min_ess < args.min_ess):
        why = (f"not converged in {all_bad}" if all_bad
               else f"ESS {min_ess:.0f} < {args.min_ess} (run longer)")
        print(f"{why}; configs left in {stage}, not copied.", file=sys.stderr)
        return
    os.makedirs(args.seeds_dir, exist_ok=True)
    copied = 0
    for i in range(K):
        rc, cfg, _ = res.get(i, (1, None, None))
        if rc == 0 and cfg and os.path.exists(cfg):
            shutil.copy2(cfg, os.path.join(args.seeds_dir,
                         build_seed_filename(args.topology, params, seed_index=copied)))
            copied += 1
    print(f"CONVERGED (all observables, ESS>={args.min_ess:.0f}) -> copied {copied} "
          f"equilibrium seeds to {args.seeds_dir}/ ({stem}_s000..s{copied - 1:03d})",
          flush=True)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--chain", action="store_true", help="Worker mode: run one chain.")
    p.add_argument("--init", default=None, help="Initial .mfd config (worker).")
    p.add_argument("--dim", type=int, default=3)
    # objective template
    p.add_argument("--n-target", type=int, required=True)
    p.add_argument("--num-facets-coef", type=float, default=0.1)
    p.add_argument("--hinge-target", type=float, default=5.0043)
    p.add_argument("--num-hinges-coef", type=float, default=0.1)
    p.add_argument("--hdv-coef", type=float, default=0.0)
    # chain controls
    p.add_argument("--coef", type=float, help="Study VDV coefficient beta (worker).")
    p.add_argument("--thin", type=int, default=5, help="Sweeps per recorded sample.")
    p.add_argument("--n-samples", type=int, default=500)
    p.add_argument("--warmup-coef", type=float, default=0.0)
    p.add_argument("--warmup-sweeps", type=int, default=0)
    p.add_argument("--warmup-hinge-coef", type=float, default=0.0,
                   help="Stiffen the edge-degree pin to this value DURING the "
                        "over-squeeze warmup (0 = keep --num-hinges-coef), so the "
                        "'below' start stays near the study edge degree instead of "
                        "drifting. Restored to --num-hinges-coef afterward.")
    p.add_argument("--side", default="above")
    p.add_argument("--replica", type=int, default=0)
    p.add_argument("--out", default=None)
    p.add_argument("--burnin-sweeps", type=int, default=0,
                   help="Sweeps at beta before recording (worker).")
    p.add_argument("--save-config", default=None,
                   help="Save the final manifold here (worker).")
    p.add_argument("--topology", default="S3")
    p.add_argument("--checkpoint-every", type=int, default=200,
                   help="Checkpoint (resumable) every this many sweeps when "
                        "--save-config is set. 0 disables. For multi-day large-N "
                        "chains keep this modest so an interruption loses little.")
    # driver controls
    p.add_argument("--seed-file", help="Plain high-VDV init (driver).")
    p.add_argument("--below-init", default=None,
                   help="Low-VDV config to start 'below' chains from (e.g. an "
                        "annealed one), held at beta so VDV rises to equilibrium. "
                        "If unset, 'below' uses the seed over-squeezed by a warmup.")
    p.add_argument("--coefs", type=float, nargs="+", help="Beta grid (driver).")
    p.add_argument("--replicas", type=int, default=2)
    p.add_argument("--burnin", type=int, default=100, help="Samples discarded as burn-in.")
    p.add_argument("--warmup-factor", type=float, default=8.0,
                   help="'below' warmup coefficient = factor * beta.")
    p.add_argument("--warmup-sweeps-driver", type=int, default=400,
                   help="Sweeps of over-squeeze warmup for 'below' chains.")
    p.add_argument("--max-workers", type=int, default=max(1, (os.cpu_count() or 4) - 2))
    p.add_argument("--max-memory-gb", type=float, default=0.0,
                   help="RAM budget for concurrent chains (0 = auto: total - 4 GB). "
                        "Concurrency is capped so ~2 MB/1k-facet chains fit; e.g. "
                        "one chain at N=1e7 on 30 GB, ~6 on 128 GB.")
    p.add_argument("--output-dir", default="data/equil_vdv")
    # production mode
    p.add_argument("--produce", action="store_true",
                   help="Produce a validated equilibrium ensemble and, if R-hat<"
                        "--rhat-max, copy the configs into --seeds-dir.")
    p.add_argument("--beta", type=float, help="Study coefficient for --produce.")
    p.add_argument("--production-burnin", type=int, default=1500,
                   help="Burn-in sweeps per production chain.")
    p.add_argument("--rhat-max", type=float, default=1.05)
    p.add_argument("--edge-tol", type=float, default=0.03,
                   help="Reported spread tolerance for the frontier OK flag.")
    p.add_argument("--min-ess", type=float, default=100.0,
                   help="Minimum effective sample size (slowest observable) "
                        "required to accept a --produce ensemble.")
    p.add_argument("--hdv-tol", type=float, default=0.1,
                   help="Deprecated (unused): superseded by quantized split-R-hat.")
    p.add_argument("--facet-tol-frac", type=float, default=0.005,
                   help="Deprecated (unused): superseded by quantized split-R-hat.")
    p.add_argument("--record-histograms", action="store_true",
                   help="Record per-sample vertex/edge degree histograms to a "
                        ".hist.npz sidecar; --produce then gates the full degree "
                        "distributions (every material bin) for convergence.")
    p.add_argument("--bin-mass-floor", type=float, default=5.0,
                   help="Minimum pooled mean count for a degree bin to be gated "
                        "(sparser bins are reported but not gated).")
    p.add_argument("--seeds-dir", default="seeds")
    args = p.parse_args()

    if args.chain:
        run_chain(args)
        return
    if args.produce:
        produce(args)
        return

    # ---- driver ----
    if not (args.seed_file and args.coefs):
        p.error("driver mode needs --seed-file and --coefs")
    args.warmup_sweeps = args.warmup_sweeps_driver
    os.makedirs(args.output_dir, exist_ok=True)

    jobs = []
    for coef in args.coefs:
        for side in ("above", "below"):
            for rep in range(args.replicas):
                out = os.path.join(args.output_dir, f"coef{coef:g}_{side}_{rep}.csv")
                jobs.append((coef, side, rep, out))
    print(f"N={args.n_target}: {len(jobs)} chains "
          f"({len(args.coefs)} coefs x 2 sides x {args.replicas} reps), "
          f"{args.n_samples} samples x {args.thin} sweeps each", flush=True)

    done = 0
    nw, per, budget = effective_workers(args)
    with ThreadPoolExecutor(max_workers=nw) as pool:
        futs = {pool.submit(spawn, args, c, s, r, o): (c, s, r) for c, s, r, o in jobs}
        for fut in as_completed(futs):
            rc = fut.result(); done += 1
            if rc != 0:
                print(f"  FAILED {futs[fut]} rc={rc}", file=sys.stderr)
            if done % max(1, len(jobs) // 10) == 0 or done == len(jobs):
                print(f"  {done}/{len(jobs)} chains done", flush=True)

    print(f"\n{'coef':>8} {'accept':>7} {'rRhatV':>7} {'rRhatE':>7} {'edgeSpr':>8} "
          f"{'tau_swp':>8} {'ESS':>6} {'<VDV>':>8} {'above':>8} {'below':>8} {'edgeDeg':>7}")
    rows = []
    for coef in args.coefs:
        r = analyze(coef, args.output_dir, args.burnin, args.thin)
        if not r:
            continue
        rows.append(r)
        eq = "OK" if (r["rhat"] < args.rhat_max and r["rhat_e"] < args.rhat_max) else "  "
        print(f"{r['coef']:>8g} {r['accept']:>7.3f} {r['rhat']:>7.3f} {r['rhat_e']:>7.3f} "
              f"{r['edge_spread']:>8.4f} {r['tau_sweeps']:>8.1f} {r['ess']:>6.0f} "
              f"{r['vdv']:>8.3f} {r['above']:>8.3f} {r['below']:>8.3f} "
              f"{r['edge']:>7.4f} {eq}")

    summ = os.path.join(args.output_dir, "summary.csv")
    with open(summ, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys())) if rows else None
        if w:
            w.writeheader(); w.writerows(rows)
    print(f"\nwrote {summ}")
    print("edge-of-equilibratability = largest coef with R-hat<1.05, small gap, usable ESS")


if __name__ == "__main__":
    main()
