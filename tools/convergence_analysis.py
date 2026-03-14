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


def split_rhat(chains):
    """
    Compute split-R-hat (Gelman et al. 2013) from a list of 1-D arrays.

    Each chain is split in half, doubling the number of chains.
    This detects non-stationarity within individual chains.
    """
    # Split each chain in half
    split = []
    for c in chains:
        n = len(c)
        half = n // 2
        if half < 2:
            return float('nan')
        split.append(c[:half])
        split.append(c[half:2 * half])

    m = len(split)                       # number of split chains
    n = len(split[0])                    # length of each split chain
    chain_means = np.array([s.mean() for s in split])
    grand_mean = chain_means.mean()

    # Between-chain variance
    B = n / (m - 1) * np.sum((chain_means - grand_mean) ** 2)

    # Within-chain variance
    W = np.mean([s.var(ddof=1) for s in split])

    if W == 0:
        return 1.0 if B == 0 else float('inf')

    # Pooled variance estimate
    var_hat = (n - 1) / n * W + B / n

    return np.sqrt(var_hat / W)


def effective_sample_size(chains):
    """
    Estimate bulk effective sample size (ESS) across chains.
    Uses the initial positive sequence estimator for autocorrelation.
    """
    all_samples = np.concatenate(chains)
    n_total = len(all_samples)
    mean = all_samples.mean()
    var = all_samples.var()
    if var == 0:
        return float('inf')

    # Compute autocorrelation using FFT
    centered = all_samples - mean
    fft = np.fft.fft(centered, n=2 * n_total)
    acf = np.fft.ifft(fft * np.conj(fft)).real[:n_total] / (var * n_total)

    # Sum autocorrelations in pairs (initial positive sequence)
    tau = 1.0
    for i in range(1, n_total // 2):
        pair_sum = acf[2 * i - 1] + acf[2 * i]
        if pair_sum < 0:
            break
        tau += 2 * pair_sum

    return n_total / tau


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
