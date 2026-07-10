#!/usr/bin/env python3
"""
Compute Gelman-Rubin split-R-hat convergence diagnostics from multiple
MCMC chain CSV files produced by manifold_sampler.

Usage:
    python3 tools/convergence_analysis.py data/convergence_chain_*.dat

Prints R-hat for key observables. R-hat < 1.01 indicates convergence.
"""

import sys
import csv
import numpy as np

# Use the library implementations
sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parent.parent / 'python'))
from discrete_differential_geometry.convergence import split_rhat, effective_sample_size


def load_chain(path):
    """Load a CSV file into a dict of observable name -> numpy array."""
    with open(path) as f:
        reader = csv.DictReader(f)
        rows = list(reader)
    if not rows:
        raise ValueError(f"Empty chain file: {path}")
    data = {}
    for key in rows[0]:
        try:
            data[key] = np.array([float(r[key]) for r in rows])
        except (ValueError, KeyError):
            pass
    return data


def main():
    if len(sys.argv) < 2:
        print(__doc__.strip())
        sys.exit(1)

    paths = sys.argv[1:]
    print(f"Loading {len(paths)} chains...")
    chain_data = [load_chain(p) for p in paths]

    # Key observables to check
    key_observables = [
        "objective",
        "volume_penalty",
        "global_curvature_penalty",
        "local_curvature_variance_penalty",
        "local_solid_angle_curvature_variance_penalty",
        "num_0_simplices",
        "num_3_simplices",
        "deg_var_0_simplices",
        "deg_var_1_simplices",
        "deg_var_2_simplices",
    ]

    # Also include any degree histogram columns present
    all_keys = set(chain_data[0].keys())
    for obs in sorted(all_keys):
        if obs.startswith("codim") and obs not in key_observables:
            key_observables.append(obs)

    # Filter to observables present in all chains
    available = set.intersection(*(set(c.keys()) for c in chain_data))
    observables = [o for o in key_observables if o in available]

    min_len = min(len(c[observables[0]]) for c in chain_data)
    print(f"Samples per chain: {min_len}")
    print()

    # Discard first half as burn-in for R-hat computation
    burn = min_len // 2
    print(f"Using samples {burn}..{min_len} (discarding first half as burn-in)")
    print()

    header = f"{'Observable':<48s}  {'R-hat':>7s}  {'ESS':>8s}  {'Mean':>12s}  {'Std':>12s}"
    print(header)
    print("-" * len(header))

    any_bad = False
    for obs in observables:
        chains = [c[obs][burn:min_len] for c in chain_data]

        rhat = split_rhat(chains)
        ess = effective_sample_size(chains)
        pooled = np.concatenate(chains)
        mean = pooled.mean()
        std = pooled.std()

        flag = ""
        if rhat > 1.1:
            flag = " *** NOT CONVERGED"
            any_bad = True
        elif rhat > 1.01:
            flag = " *  MARGINAL"
            any_bad = True

        print(f"{obs:<48s}  {rhat:7.4f}  {ess:8.0f}  {mean:12.4f}  {std:12.4f}{flag}")

    print()
    if any_bad:
        print("WARNING: Some observables have R-hat > 1.01.")
        print("Consider increasing eqSweepsPerStep, growthStepSize, or maxSweeps.")
    else:
        print("All observables have R-hat < 1.01 — chains appear converged.")


if __name__ == "__main__":
    main()
