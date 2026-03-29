"""Convergence diagnostics for MCMC chains.

Provides split R-hat (Gelman et al. 2013) and effective sample size (ESS)
estimators for assessing convergence of multiple MCMC chains.
"""

import numpy as np


def split_rhat(chains):
    """
    Compute split-R-hat (Gelman et al. 2013) from a list of 1-D arrays.

    Each chain is split in half, doubling the number of chains.
    This detects non-stationarity within individual chains.

    Parameters
    ----------
    chains : list of array-like
        Each element is a 1-D array of samples from one chain.

    Returns
    -------
    float
        The split R-hat statistic. Values < 1.01 indicate convergence.
    """
    # Split each chain in half
    split = []
    for c in chains:
        c = np.asarray(c, dtype=float)
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

    return float(np.sqrt(var_hat / W))


def effective_sample_size(chains):
    """
    Estimate bulk effective sample size (ESS) across chains.

    Uses the initial positive sequence estimator for autocorrelation.

    Parameters
    ----------
    chains : list of array-like
        Each element is a 1-D array of samples from one chain.

    Returns
    -------
    float
        The estimated effective sample size.
    """
    all_samples = np.concatenate([np.asarray(c, dtype=float) for c in chains])
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
