"""Convergence diagnostics for MCMC chains.

Provides split R-hat (Gelman et al. 2013) and effective sample size (ESS)
estimators for assessing convergence of multiple MCMC chains.
"""

import numpy as np
from scipy.special import ndtri
from scipy.stats import rankdata


def _split(chains):
    """Split each chain in half, doubling the chain count (catches within-chain
    non-stationarity). Returns None if any chain is too short."""
    out = []
    for c in chains:
        c = np.asarray(c, dtype=float)
        half = len(c) // 2
        if half < 2:
            return None
        out.append(c[:half])
        out.append(c[half:2 * half])
    return out


def _gelman_rubin(chains, var_floor=0.0):
    """Gelman-Rubin R-hat on the given (already-split) chains.

    ``var_floor`` floors the within-chain variance W at a physical minimum (the
    variance Delta**2/12 of one discretization quantum). This keeps R-hat finite
    and meaningful for near-deterministic / tightly-pinned observables, whose true
    within-chain variance collapses toward zero and would otherwise make classic
    R-hat blow up on physically-negligible between-chain offsets. With var_floor=0
    this is the textbook estimator."""
    m = len(chains)
    if m < 2:
        return float('nan')
    n = min(len(s) for s in chains)
    chains = [s[:n] for s in chains]
    chain_means = np.array([s.mean() for s in chains])
    B = n / (m - 1) * np.sum((chain_means - chain_means.mean()) ** 2)
    W = max(np.mean([s.var(ddof=1) for s in chains]), var_floor)
    if W == 0:
        return 1.0 if B == 0 else float('inf')
    var_hat = (n - 1) / n * W + B / n
    return float(np.sqrt(var_hat / W))


def split_rhat(chains):
    """
    Compute split-R-hat (Gelman et al. 2013) from a list of 1-D arrays.

    Each chain is split in half, doubling the number of chains, which detects
    non-stationarity within individual chains.

    NOTE: classic R-hat divides by the within-chain variance, so it diverges for
    near-constant series (e.g. a tightly-pinned observable). Use
    ``rank_normalized_rhat`` for those; it is the recommended default.

    Parameters
    ----------
    chains : list of array-like
        Each element is a 1-D array of samples from one chain.

    Returns
    -------
    float
        The split R-hat statistic. Values < 1.01 indicate convergence.
    """
    split = _split(chains)
    return float('nan') if split is None else _gelman_rubin(split)


def _rank_normalize(chains):
    """Pool, average-rank (ties shared), and Blom-transform to normal scores;
    split back into per-chain arrays."""
    arrs = [np.asarray(c, dtype=float) for c in chains]
    lengths = [len(a) for a in arrs]
    pooled = np.concatenate(arrs)
    z = ndtri((rankdata(pooled) - 3.0 / 8.0) / (len(pooled) - 0.25))
    out, i = [], 0
    for L in lengths:
        out.append(z[i:i + L]); i += L
    return out


def rank_normalized_rhat(chains):
    """
    Rank-normalized split-R-hat (Vehtari et al. 2021).

    Ranks are pooled across chains and mapped to normal scores before computing
    split-R-hat, which makes it invariant to monotone transforms, robust to heavy
    tails, and -- crucially -- well-behaved for near-constant / tightly-pinned
    observables where classic R-hat divides by a vanishing within-chain variance
    and blows up. Interleaved (converged) chains give ~1; chains sitting at
    different values give a finite, meaningful value.

    Returns max(rank R-hat, folded rank R-hat), where the folded term is R-hat on
    |x - median(pooled)| and catches chains that agree in location but differ in
    scale/variance. Convergence: < 1.01.

    Parameters
    ----------
    chains : list of array-like
        Each element is a 1-D array of samples from one chain.
    """
    split = _split(chains)
    if split is None or len(split) < 2:
        return float('nan')
    bulk = _gelman_rubin(_rank_normalize(split))
    med = np.median(np.concatenate([np.asarray(c, dtype=float) for c in chains]))
    folded = _gelman_rubin(_rank_normalize([np.abs(s - med) for s in split]))
    return float(max(bulk, folded))


def quantized_split_rhat(chains, quantum=0.0):
    """
    Variance-floored split-R-hat: one convergence statistic for observables of
    any nature, fluctuating or near-deterministic.

    Standard split-R-hat divides the total-variance estimate by the within-chain
    variance W, so it diverges for a tightly-pinned / near-deterministic observable
    (W -> 0) even when the chains differ by a physically-negligible amount. Here W
    is floored at ``quantum**2 / 12`` -- the variance of a single discretization
    step, ``quantum`` being the minimum change in the observable when the other
    underlying integer counts are held fixed. For a genuinely fluctuating
    observable the true W dominates and the floor is inert, recovering ordinary
    split-R-hat; for a near-deterministic one the floor sets the scale, so the
    statistic reports "chains agree to within a few quanta" as convergence instead
    of blowing up. This makes a per-observable rank-normalization unnecessary: the
    floor cures the near-constant pathology directly and in physical units.

    The between-chain term uses every chain's mean (not a peak-to-peak range), so
    it does not inflate with the number of chains.

    Returns max(bulk, folded), where folded is the same statistic on
    ``|x - median(pooled)|`` (catches chains that agree in location but differ in
    spread). Convergence: < 1.01.

    Parameters
    ----------
    chains : list of array-like
        Each element is a 1-D array of samples from one chain.
    quantum : float, optional
        Minimum separation between attainable values with the other underlying
        counts held fixed. 0 (default) recovers ordinary split-R-hat.
    """
    split = _split(chains)
    if split is None or len(split) < 2:
        return float('nan')
    var_floor = (quantum * quantum) / 12.0
    bulk = _gelman_rubin(split, var_floor)
    med = np.median(np.concatenate([np.asarray(c, dtype=float) for c in chains]))
    folded = _gelman_rubin([np.abs(s - med) for s in split], var_floor)
    return float(max(bulk, folded))


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
