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


def integrated_autocorrelation_time(x):
    """
    Estimate the integrated autocorrelation time tau of a single 1-D series.

    Uses the initial-positive-sequence (Geyer) estimator on the FFT
    autocorrelation. tau = 1 + 2 * sum_{k>=1} rho(k), truncated at the first
    non-positive consecutive pair. The number of effectively independent
    samples in the series is len(x) / tau.

    Parameters
    ----------
    x : array-like
        A 1-D array of samples from one chain (assumed post-warmup).

    Returns
    -------
    float
        The integrated autocorrelation time (>= 1). Returns 1.0 for a
        constant series (zero variance) and nan for series too short to
        estimate (< 4 samples).
    """
    x = np.asarray(x, dtype=float)
    n = len(x)
    if n < 4:
        return float('nan')
    mean = x.mean()
    var = x.var()
    if var == 0:
        return 1.0

    centered = x - mean
    fft = np.fft.fft(centered, n=2 * n)
    acf = np.fft.ifft(fft * np.conj(fft)).real[:n] / (var * n)

    tau = 1.0
    for i in range(1, n // 2):
        pair_sum = acf[2 * i - 1] + acf[2 * i]
        if pair_sum < 0:
            break
        tau += 2 * pair_sum

    return tau


def weighted_ess(weights):
    """
    Kish effective sample size of a set of importance weights.

    For a reweighted estimator with weights w_i, the effective number of
    independent samples is (sum w)^2 / sum(w^2). This measures the cost of
    non-uniform weights *only*; it does not account for autocorrelation
    between samples. Multiply by (n / tau) / n, i.e. combine with the
    autocorrelation ESS, to get the overall effective sample size.

    Parameters
    ----------
    weights : array-like
        Non-negative importance weights (e.g. 1/V(x) per sample).

    Returns
    -------
    float
        Kish ESS in [1, len(weights)]. Returns len(weights) for uniform
        weights and 0.0 for an empty or all-zero weight vector.
    """
    w = np.asarray(weights, dtype=float)
    if len(w) == 0:
        return 0.0
    s1 = w.sum()
    s2 = (w * w).sum()
    if s2 == 0:
        return 0.0
    return float(s1 * s1 / s2)
