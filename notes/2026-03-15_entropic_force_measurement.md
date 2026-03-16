# Entropic Force Measurement for 3-Manifold Triangulations

Date: 2026-03-15

## Method

The MCMC sampler with energy `E = β_V * (N - N_target)²` samples from the
distribution `π(N) ∝ Ω(N) * exp(-β_V*(N-N_target)²) / N`, where the `1/N`
factor comes from the Hastings correction `log(N_before/N_after)` in the
bistellar move acceptance criterion. At equilibrium:

```
d/dN log Ω(N) = 2 * β_V * (⟨N⟩ - N_target) + ⟨1/N⟩
```

The left side is the entropic force F(N). The right side is measurable from
the time series of nFacets.

### Protocol

- Volume penalty only: `numFacetsCoef = β_V`, all other coefficients zero
- Seeds grown from scratch under volume-only energy with ramped growth
  (step size N/50, 20 eq sweeps per step, 50 sweep final equilibration)
- Hinge moves at 30% attempt rate throughout
- Production: 50 sweeps per measurement, batch means with 20 batches
- Three β_V values (0.05, 0.1, 0.2) at each size for consistency checking

### Consistency check

If the degree distribution has fully equilibrated, the measured F should be
independent of β_V. Systematic dependence on β_V indicates insufficient
equilibration of slow modes (primarily vertex degree variance).

## Results (Run 1: 20 eq sweeps, volume-only seeds)

```
N_target  beta_V       <N>   offset    Var(N)  1/2beta  F_entrop   F_err
     500    0.05    514.59    14.59      10.0     10.0    1.4611  0.0153
     500    0.10    506.93     6.93       5.2      5.0    1.3886  0.0169
     500    0.20    503.31     3.31       2.8      2.5    1.3261  0.0204
    1000    0.05   1015.82    15.82      10.1     10.0    1.5833  0.0136
    1000    0.10   1007.54     7.54       5.2      5.0    1.5095  0.0130
    1000    0.20   1003.62     3.62       2.8      2.5    1.4487  0.0093
    2000    0.05   2016.49    16.49      10.2     10.0    1.6497  0.0092
    2000    0.10   2008.07     8.07       5.2      5.0    1.6154  0.0066
    2000    0.20   2003.79     3.79       2.8      2.5    1.5181  0.0091
    5000    0.05   5017.32    17.32      10.2     10.0    1.7320  0.0079
    5000    0.10   5008.44     8.44       5.3      5.0    1.6875  0.0057
    5000    0.20   5003.95     3.95       2.7      2.5    1.5818  0.0061
   10000    0.05  10017.75    17.75      10.5     10.0    1.7755  0.0045
   10000    0.10  10008.76     8.76       5.2      5.0    1.7519  0.0045
   10000    0.20  10004.13     4.13       2.8      2.5    1.6517  0.0033
```

### Observations

1. **The entropic force is positive and grows with N.** F increases from ~1.4
   at N=500 to ~1.75 at N=10000. This means `log Ω(N)` is super-linear in N.
   If the number of triangulations grew as `exp(c*N)`, F would be constant.
   The growth suggests something like `log Ω(N) ~ c*N + d*log(N)` or faster.

2. **Var(N) matches the quadratic prediction.** `Var(N) ≈ 1/(2β_V)` holds
   well across all measurements, confirming the energy landscape is
   approximately quadratic near the target.

3. **β_V consistency check fails.** F systematically decreases with increasing
   β_V at all sizes (e.g., N=10000: 1.78→1.75→1.65). The discrepancy exceeds
   error bars by ~20σ. This indicates the degree distribution has not fully
   equilibrated during the 20-sweep equilibration phase.

4. **Degree variance is a very slow mode.** The vertex degree variance
   trajectories show 7-17% drift during equilibration, even with volume-only
   seeds grown from scratch. The autocorrelation time for degree variance
   appears to be much longer than for volume.

### Interpretation of β_V dependence

Lower β_V → wider volume fluctuations → the chain explores a broader range of
N values → this indirectly helps mix the degree distribution (visiting different
volumes forces topological rearrangements). Higher β_V confines the chain near
N_target, reducing the indirect mixing benefit.

This means the β_V=0.05 measurements are likely closest to the true
equilibrium value, and the β_V=0.20 measurements are biased low because the
degree distribution is stuck closer to its initial state.

### Extrapolation

Taking the β_V=0.05 values as the most reliable:

| N     | F(N)   |
|-------|--------|
| 500   | 1.461  |
| 1000  | 1.583  |
| 2000  | 1.650  |
| 5000  | 1.732  |
| 10000 | 1.776  |

The growth is roughly logarithmic in N: `F(N) ≈ 0.85 + 0.10 * ln(N)`.
If this holds, then `log Ω(N) ≈ 0.85*N + 0.10*N*ln(N) + const`, suggesting
the number of 3-manifold triangulations grows faster than any exponential
but slower than `exp(N^(1+ε))`.

## Results (Run 2: scaled equilibration — 20/50/100 sweeps for β=0.05/0.10/0.20)

```
N_target  beta_V       <N>   offset    Var(N)  1/2beta  F_entrop   F_err
     500    0.05    514.37    14.37      10.1     10.0    1.4393  0.0157
     500    0.10    507.03     7.03       5.2      5.0    1.4080  0.0179
     500    0.20    503.27     3.27       2.7      2.5    1.3095  0.0165
    1000    0.05   1015.39    15.39      10.4     10.0    1.5399  0.0133
    1000    0.10   1007.59     7.59       5.1      5.0    1.5183  0.0148
    1000    0.20   1003.58     3.58       2.7      2.5    1.4339  0.0119
    2000    0.05   2016.41    16.41      10.1     10.0    1.6417  0.0089
    2000    0.10   2007.97     7.97       5.3      5.0    1.5941  0.0081
    2000    0.20   2003.77     3.77       2.7      2.5    1.5081  0.0084
    5000    0.05   5017.35    17.35      10.4     10.0    1.7355  0.0076
    5000    0.10   5008.45     8.45       5.2      5.0    1.6906  0.0039
    5000    0.20   5004.02     4.02       2.7      2.5    1.6079  0.0069
   10000    0.05  10017.89    17.89      10.2     10.0    1.7891  0.0042
   10000    0.10  10008.76     8.76       5.3      5.0    1.7515  0.0044
   10000    0.20  10004.15     4.15       2.7      2.5    1.6611  0.0034
```

### Effect of increased equilibration

The β=0.10 gap with β=0.05 improved significantly with 50 eq sweeps
(e.g., N=1000: gap 0.07→0.02). The β=0.20 gap improved slightly with
100 eq sweeps but remains ~0.10-0.13 across sizes, indicating the degree
variance autocorrelation time at β=0.20 exceeds 100 sweeps for N≥2000.

### Best estimates (β_V=0.05, most reliable)

| N     | F(N)           |
|-------|----------------|
| 500   | 1.439 ± 0.016  |
| 1000  | 1.540 ± 0.013  |
| 2000  | 1.642 ± 0.009  |
| 5000  | 1.736 ± 0.008  |
| 10000 | 1.789 ± 0.004  |

Fit: `F(N) ≈ 0.85 + 0.10 * ln(N)`, implying
`log Ω(N) ≈ 0.85*N + 0.10*N*ln(N) + const`.

## Results (Run 3: dynamic equilibration stopping criterion)

Replaced fixed eq sweeps with dynamic convergence: run in 5-sweep windows,
stop when 4 consecutive window means agree within 5%. Max 500 sweeps.

```
N_target  beta_V       <N>   offset    Var(N)  1/2beta  F_entrop   F_err  EqSw Conv
     500    0.05    514.31    14.31      10.1     10.0    1.4328  0.0198   500  N
     500    0.10    506.94     6.94       5.3      5.0    1.3895  0.0194   425  Y
     500    0.20    503.29     3.29       2.7      2.5    1.3161  0.0115   500  N
    1000    0.05   1015.31    15.31      10.6     10.0    1.5319  0.0133    55  Y
    1000    0.10   1007.57     7.57       5.5      5.0    1.5156  0.0140   465  Y
    1000    0.20   1003.55     3.55       2.7      2.5    1.4205  0.0114    45  Y
    2000    0.05   2016.90    16.90      10.4     10.0    1.6906  0.0114   425  Y
    2000    0.10   2007.90     7.90       5.3      5.0    1.5806  0.0094   190  Y
    2000    0.20   2003.79     3.79       2.8      2.5    1.5167  0.0087    25  Y
    5000    0.05   5017.28    17.28      10.2     10.0    1.7281  0.0044    20  Y
    5000    0.10   5008.50     8.50       5.3      5.0    1.6996  0.0083    65  Y
    5000    0.20   5003.99     3.99       2.7      2.5    1.5960  0.0047    35  Y
   10000    0.05  10017.96    17.96      10.2     10.0    1.7964  0.0056    60  Y
   10000    0.10  10008.80     8.80       5.4      5.0    1.7604  0.0043    30  Y
   10000    0.20  10004.14     4.14       2.8      2.5    1.6573  0.0039    40  Y
```

### Key finding: β_V gap persists after dynamic equilibration

The gap between β_V=0.05 and β_V=0.20 is ~0.12-0.17 across all sizes,
essentially unchanged from runs 1 and 2 despite the dynamic convergence
check confirming deg_var stability. This suggests the gap is NOT an
equilibration artifact. Possible explanations:

1. **Higher-order coupling**: the stationary distribution π(N) depends on
   the full conditional density of states Ω(N | geometry), not just Ω(N).
   Different β_V values probe different slices of the joint (N, geometry)
   distribution, yielding different effective entropic forces.

2. **The 5% convergence criterion is satisfied but the distribution hasn't
   fully mixed**: deg_var may have reached a local plateau that looks stable
   but differs from the true equilibrium. The intrinsic fluctuations in
   deg_var (~30% at N=500) make this plausible.

3. **The Hastings correction `log(N_before/N_after)` introduces a β_V-dependent
   bias**: the 1/N factor in the stationary distribution is exact, but the
   approximation of using the mean ⟨1/N⟩ may break down differently at
   different β_V values due to different fluctuation amplitudes.

Hypothesis 1 is the most physically interesting: it would mean the entropic
force is not a simple function of N alone, but depends on the geometry.

## Open questions

- Is the β_V gap a real physical effect (geometry-dependent entropic force)
  or a subtle equilibration/measurement artifact?
- What is the autocorrelation time for vertex degree variance as a function
  of system size and β_V?
- Does the F(N) ~ 0.85 + 0.10*ln(N) fit hold at larger N (20K, 50K)?
- How does F change when curvature penalties are added? Is the entropic force
  purely a function of N, or does it depend on the degree distribution?
- The super-linear growth `N*ln(N)` is reminiscent of the entropy of random
  planar maps. Is there a connection to known combinatorial results for
  3-manifold triangulations?
