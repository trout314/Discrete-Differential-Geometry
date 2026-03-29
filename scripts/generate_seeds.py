#!/usr/bin/env python3
"""
Generate equilibrated seed triangulations with multi-chain R-hat convergence gating.

For each (dim, size), runs K=4 independent chains from separate standardSphere -> ramped_grow
starting points. Samples observables every `thin` sweeps, checking split R-hat convergence
every 50 samples. A seed is accepted only when all observables show R-hat < 1.01 for two
consecutive checks.

Usage:
    python scripts/generate_seeds.py --dim 3 --sizes 1000,2000,5000
    python scripts/generate_seeds.py --dim 2 --all
    python scripts/generate_seeds.py --all-dims --all
"""

import argparse
import csv
import math
import os
import sys
import time

import numpy as np

sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parent.parent / 'python'))
import discrete_differential_geometry as ddg
from discrete_differential_geometry.convergence import split_rhat


NUM_CHAINS = 4
RHAT_THRESHOLD = 1.01
CONSECUTIVE_PASSES_REQUIRED = 2
CHECK_INTERVAL = 50  # samples between convergence checks

DIM_CONFIGS = {
    2: {
        'sizes': [200, 500, 1000, 2000, 5000, 10000, 20000, 50000, 100000,
                  200000, 500000, 1000000, 2000000],
        'thin': 2,
        'warmup_factor': 10,  # warmup = factor * N/1000 sweeps
        'max_samples': 500,
        'hinge_degree_target': 6.0,
        'codim3_degree_variance_coef': 0.0,
    },
    3: {
        'sizes': [200, 500, 1000, 2000, 5000, 10000, 20000, 50000, 100000,
                  200000, 500000, 1000000, 2000000],
        'thin': 5,
        'warmup_factor': 20,
        'max_samples': 500,
        'hinge_degree_target': 5.1,
        'codim3_degree_variance_coef': 0.1,
    },
    4: {
        'sizes': [200, 500, 1000, 2000, 5000, 10000, 20000, 50000],
        'thin': 10,
        'warmup_factor': 30,
        'max_samples': 300,
        'hinge_degree_target': 2 * math.pi / math.acos(0.25),
        'codim3_degree_variance_coef': 0.1,
    },
}


def make_sampler_params(dim, size, config):
    """Create SamplerParams for the given dimension and size."""
    return ddg.SamplerParams(
        num_facets_target=size,
        hinge_degree_target=config['hinge_degree_target'],
        num_facets_coef=0.1,
        num_hinges_coef=0.05,
        hinge_degree_variance_coef=0.2,
        codim3_degree_variance_coef=config['codim3_degree_variance_coef'],
    )


def collect_observables(sampler):
    """Snapshot all tracked observables from a sampler."""
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


def run_chain(dim, size, config, chain_idx):
    """Grow and sample one independent chain. Returns list of obs dicts."""
    print(f"  Chain {chain_idx}: growing from standardSphere({dim})...")
    m = ddg.Manifold.standard_sphere(dim)
    params = make_sampler_params(dim, size, config)
    sampler = ddg.ManifoldSampler(m, params)

    step_size = max(100, size // 20)
    sampler.ramped_grow(size, step_size=step_size, eq_sweeps_per_step=3)
    print(f"  Chain {chain_idx}: grown to {sampler.manifold.num_facets} facets")

    warmup = max(1, config['warmup_factor'] * size // 1000)
    print(f"  Chain {chain_idx}: running {warmup} warmup sweeps...")
    sampler.run(sweeps=warmup)

    return sampler


def generate_seed(dim, size, output_dir='standard_triangulations', data_dir='data/seed_gen'):
    """Orchestrate K chains with adaptive convergence stopping."""
    config = DIM_CONFIGS[dim]
    thin = config['thin']
    max_samples = config['max_samples']

    print(f"\n{'='*60}")
    print(f"Generating {dim}d seed with {size} facets")
    print(f"  {NUM_CHAINS} chains, thin={thin}, max_samples={max_samples}")
    print(f"{'='*60}")

    os.makedirs(data_dir, exist_ok=True)
    os.makedirs(output_dir, exist_ok=True)

    # Create all chains (truly independent: separate grow from standardSphere)
    samplers = []
    for k in range(NUM_CHAINS):
        sampler = run_chain(dim, size, config, k)
        samplers.append(sampler)
        ddg.gc_collect()

    # Collect samples with periodic convergence checks
    chain_data = [[] for _ in range(NUM_CHAINS)]
    consecutive_passes = 0
    converged = False
    t_start = time.monotonic()

    for sample_idx in range(max_samples):
        # Run thinning sweeps on all chains
        for k, sampler in enumerate(samplers):
            sampler.run(sweeps=thin)

        # Collect observables
        for k, sampler in enumerate(samplers):
            obs = collect_observables(sampler)
            obs['sample'] = sample_idx
            chain_data[k].append(obs)

        # Check convergence every CHECK_INTERVAL samples
        n_samples = sample_idx + 1
        if n_samples >= 2 * CHECK_INTERVAL and n_samples % CHECK_INTERVAL == 0:
            passed, report = check_convergence(chain_data)
            elapsed = time.monotonic() - t_start
            print(f"  Sample {n_samples}/{max_samples} ({elapsed:.0f}s): "
                  f"max R-hat = {report['max_rhat']:.4f} "
                  f"({'PASS' if passed else 'FAIL'})")

            if passed:
                consecutive_passes += 1
                if consecutive_passes >= CONSECUTIVE_PASSES_REQUIRED:
                    converged = True
                    print(f"  Converged after {n_samples} samples!")
                    break
            else:
                consecutive_passes = 0

    if not converged:
        print(f"  WARNING: Did not converge after {max_samples} samples. "
              f"Saving anyway for inspection.")

    # Save chain CSVs
    for k in range(NUM_CHAINS):
        csv_path = os.path.join(data_dir, f'{dim}d_{size}_chain_{k}.csv')
        save_chain_csv(chain_data[k], csv_path)
        print(f"  Saved chain data: {csv_path}")

    # Save chain 0's final manifold as the seed
    seed_path = os.path.join(output_dir, f'{dim}d_equilibrated_{size}.mfd')
    samplers[0].manifold.save(seed_path)
    print(f"  Saved seed: {seed_path}")

    # Clean up
    del samplers
    ddg.gc_collect()
    ddg.gc_minimize()

    return seed_path, converged


def check_convergence(chain_data):
    """Compute split R-hat across all chains for every observable."""
    obs_names = [k for k in chain_data[0][0].keys() if k != 'sample']
    report = {'observables': {}, 'max_rhat': 0.0}
    all_pass = True

    for obs_name in obs_names:
        chains = []
        for k in range(len(chain_data)):
            values = np.array([s[obs_name] for s in chain_data[k]])
            chains.append(values)

        rhat = split_rhat(chains)
        report['observables'][obs_name] = rhat

        if not np.isnan(rhat):
            report['max_rhat'] = max(report['max_rhat'], rhat)
            if rhat >= RHAT_THRESHOLD:
                all_pass = False

    return all_pass, report


def save_chain_csv(data, path):
    """Save chain observable data as CSV."""
    if not data:
        return
    keys = list(data[0].keys())
    with open(path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=keys)
        writer.writeheader()
        writer.writerows(data)


def main():
    parser = argparse.ArgumentParser(
        description='Generate equilibrated seed triangulations with convergence gating')
    parser.add_argument('--dim', type=int, choices=[2, 3, 4],
                        help='Manifold dimension')
    parser.add_argument('--sizes', type=str,
                        help='Comma-separated list of sizes')
    parser.add_argument('--all', action='store_true', dest='all_sizes',
                        help='Generate all sizes for the given dimension(s)')
    parser.add_argument('--all-dims', action='store_true',
                        help='Generate for all dimensions (2, 3, 4)')
    parser.add_argument('--output-dir', default='standard_triangulations',
                        help='Directory for seed .mfd files')
    parser.add_argument('--data-dir', default='data/seed_gen',
                        help='Directory for chain CSV data')

    args = parser.parse_args()

    if args.all_dims:
        dims = [2, 3, 4]
    elif args.dim is not None:
        dims = [args.dim]
    else:
        parser.error('Must specify --dim or --all-dims')

    results = []
    for dim in dims:
        config = DIM_CONFIGS[dim]
        if args.all_sizes:
            sizes = config['sizes']
        elif args.sizes:
            sizes = [int(s) for s in args.sizes.split(',')]
        else:
            parser.error('Must specify --sizes or --all')

        for size in sizes:
            seed_path, converged = generate_seed(
                dim, size,
                output_dir=args.output_dir,
                data_dir=args.data_dir,
            )
            results.append((dim, size, seed_path, converged))

    # Summary
    print(f"\n{'='*60}")
    print("Summary")
    print(f"{'='*60}")
    for dim, size, path, converged in results:
        status = "CONVERGED" if converged else "NOT CONVERGED"
        print(f"  {dim}d N={size}: {status} -> {path}")


if __name__ == '__main__':
    main()
