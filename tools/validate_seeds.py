#!/usr/bin/env python3
"""
Validate that existing seed triangulations are stationary (equilibrated).

Runs K=4 truly independent reference chains (each from a fresh standardSphere
-> ramped_grow) to estimate the stationary distribution. Then checks:

1. Cross-chain R-hat on the reference chains confirms they've converged
   (R-hat < 1.05), so the reference distribution is trustworthy.
2. Geweke z-test: the seed's initial observables are compared against the
   pooled reference distribution. |z| > 3 means the seed is an outlier.

Usage:
    python tools/validate_seeds.py standard_triangulations/3d_equilibrated_*.mfd --dim 3
    python tools/validate_seeds.py standard_triangulations/equilibrated_*.mfd --dim 3
"""

import argparse
import math
import sys

import numpy as np

sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parent.parent / 'python'))
import discrete_differential_geometry as ddg
from discrete_differential_geometry.convergence import split_rhat, effective_sample_size


NUM_REF_CHAINS = 4
RHAT_THRESHOLD = 1.05
GEWEKE_Z_THRESHOLD = 3.0
NUM_SAMPLES = 200

DIM_DEFAULTS = {
    2: {'thin': 2, 'warmup_factor': 10, 'hinge_degree_target': 6.0, 'codim3_dv': 0.0},
    3: {'thin': 5, 'warmup_factor': 20, 'hinge_degree_target': 5.1, 'codim3_dv': 0.1},
    4: {'thin': 10, 'warmup_factor': 30,
        'hinge_degree_target': 2 * math.pi / math.acos(0.25), 'codim3_dv': 0.1},
}


def collect_observables(sampler):
    """Snapshot observables from a sampler."""
    mfd = sampler.manifold
    dim = mfd.dimension
    obs = {
        'current_objective': sampler.current_objective,
        'num_facets': mfd.num_facets,
    }
    fv = mfd.f_vector
    for i, count in enumerate(fv):
        obs[f'f_{i}'] = count
    for d in range(max(0, dim - 2), dim):
        obs[f'degree_variance_{d}'] = mfd.degree_variance(d)
        obs[f'mean_degree_{d}'] = mfd.mean_degree(d)
    return obs


def make_params(size, defaults):
    return ddg.SamplerParams(
        num_facets_target=size,
        hinge_degree_target=defaults['hinge_degree_target'],
        num_facets_coef=0.1,
        num_hinges_coef=0.05,
        hinge_degree_variance_coef=0.2,
        codim3_degree_variance_coef=defaults['codim3_dv'],
    )


def validate_seed(path, dim):
    """Validate a single seed file. Returns (passed, max_rhat, max_z)."""
    defaults = DIM_DEFAULTS[dim]
    thin = defaults['thin']

    print(f"\nValidating: {path}")
    m = ddg.Manifold.load(path, dim)
    size = m.num_facets
    print(f"  Loaded: dim={dim}, num_facets={size}")

    params = make_params(size, defaults)
    warmup = max(1, defaults['warmup_factor'] * size // 1000)

    # Record seed's initial observables BEFORE any MCMC
    seed_sampler = ddg.ManifoldSampler(m.copy(), params)
    initial_obs = collect_observables(seed_sampler)
    del seed_sampler

    # Build K truly independent reference chains
    ref_chain_data = [[] for _ in range(NUM_REF_CHAINS)]
    ref_samplers = []
    step_size = max(100, size // 20)
    for k in range(NUM_REF_CHAINS):
        print(f"  Ref chain {k}: growing from standardSphere({dim})...")
        ref_m = ddg.Manifold.standard_sphere(dim)
        s = ddg.ManifoldSampler(ref_m, params)
        s.ramped_grow(size, step_size=step_size, eq_sweeps_per_step=3)
        print(f"  Ref chain {k}: grown to {s.manifold.num_facets} facets, "
              f"running {warmup} warmup sweeps...")
        s.run(sweeps=warmup)
        ref_samplers.append(s)
        ddg.gc_collect()

    # Collect samples from all reference chains
    print(f"  Collecting {NUM_SAMPLES} samples (thin={thin} sweeps) "
          f"from {NUM_REF_CHAINS} reference chains...")
    for sample_idx in range(NUM_SAMPLES):
        for k, s in enumerate(ref_samplers):
            s.run(sweeps=thin)
            obs = collect_observables(s)
            ref_chain_data[k].append(obs)

    del ref_samplers
    ddg.gc_collect()

    # Compute diagnostics
    obs_names = [k for k in ref_chain_data[0][0].keys()]
    print()
    header = (f"  {'Observable':<30s}  {'R-hat':>7s}  {'|z|':>6s}  "
              f"{'Init':>12s}  {'Ref Mean':>12s}  {'Ref Std':>12s}")
    print(header)
    print("  " + "-" * (len(header) - 2))

    all_pass = True
    max_rhat = 0.0
    max_z = 0.0
    ref_converged = True

    for obs_name in obs_names:
        # Cross-chain R-hat on reference chains
        chains = []
        for k in range(NUM_REF_CHAINS):
            values = np.array([s[obs_name] for s in ref_chain_data[k]])
            chains.append(values)
        rhat = split_rhat(chains)

        # Pooled reference distribution
        pooled = np.concatenate(chains)
        ref_mean = pooled.mean()
        ref_std = pooled.std()

        # Geweke z-score
        init_val = initial_obs[obs_name]
        if ref_std > 0:
            z = abs(init_val - ref_mean) / ref_std
        else:
            z = 0.0

        if not np.isnan(rhat):
            max_rhat = max(max_rhat, rhat)
        max_z = max(max_z, z)

        # Flags
        flag = ""
        rhat_bad = not np.isnan(rhat) and rhat >= RHAT_THRESHOLD
        z_bad = z >= GEWEKE_Z_THRESHOLD and ref_std > 0

        if rhat_bad:
            ref_converged = False

        if z_bad:
            flag = " *** FAIL (z)"
            all_pass = False
        elif rhat_bad:
            flag = " !  REF NOT CONVERGED"
        elif not np.isnan(rhat) and rhat >= RHAT_THRESHOLD - 0.02:
            flag = " *  ref marginal"
        elif z >= GEWEKE_Z_THRESHOLD - 0.5 and ref_std > 0:
            flag = " *  marginal (z)"

        print(f"  {obs_name:<30s}  {rhat:7.4f}  {z:6.2f}  "
              f"{init_val:12.4f}  {ref_mean:12.4f}  {ref_std:12.4f}{flag}")

    if not ref_converged:
        print("\n  WARNING: Reference chains did not fully converge (R-hat >= "
              f"{RHAT_THRESHOLD}). z-test results may be unreliable.")
        all_pass = False

    ddg.gc_minimize()
    return all_pass, max_rhat, max_z


def main():
    parser = argparse.ArgumentParser(description='Validate seed triangulations are stationary')
    parser.add_argument('files', nargs='+', help='Seed .mfd files to validate')
    parser.add_argument('--dim', type=int, required=True, choices=[2, 3, 4],
                        help='Manifold dimension')

    args = parser.parse_args()

    results = []
    for path in args.files:
        passed, max_rhat, max_z = validate_seed(path, args.dim)
        results.append((path, passed, max_rhat, max_z))

    # Summary
    print(f"\n{'='*60}")
    print("Validation Summary")
    print(f"{'='*60}")
    for path, passed, max_rhat, max_z in results:
        status = "PASS" if passed else "FAIL"
        print(f"  [{status}] {path} (ref R-hat = {max_rhat:.4f}, max |z| = {max_z:.2f})")

    n_fail = sum(1 for _, p, _, _ in results if not p)
    if n_fail > 0:
        print(f"\n{n_fail}/{len(results)} seeds FAILED validation")
        sys.exit(1)
    else:
        print(f"\nAll {len(results)} seeds PASSED validation.")


if __name__ == '__main__':
    main()
