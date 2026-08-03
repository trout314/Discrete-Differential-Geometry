---
name: five-illegal-knot
description: Structure of the 5-illegal-vertex complex (the elementary illegal defect) in pushed TCP crystals
metadata: 
  node_type: memory
  type: project
  originSessionId: 1ed7f5d0-df5c-4552-b65e-5060180eafe2
  modified: 2026-07-21T17:20:55.462Z
---

2026-07-21: characterized the recurring 5-illegal-vertex complex `5(0·0·0·0|5ill)`
that forms when a TCP crystal (R, C15) is pushed off its native edge degree. See
[[c15-z14-dope-ensemble]], [[crystal-grains-tool]].

**STRUCTURE (exact, 11/11 instances across a two-sided equilibrated R ensemble):**
the illegal-edge core is ALWAYS `(3,4,4)` — exactly one degree-3 edge + two
degree-4 edges knotted over 5 vertices (rest of the ~50 spokes legal 5/6). The 5
vertices' coordinations vary around this fixed core (most common 12,13,13,14,15;
n6 usually 2,3,3,3,4). It is a COMPACT POSITIVE-DISCLINATION KNOT: degree-3/4 edges
are strongly under-coordinated (<<5 tets = strong +disclination), the opposite of
the degree-6 (-disclination) lines that carry a raised mean edge degree. So by
global Gauss-Bonnet it is the topological COUNTER-CHARGE to the disclination
network — raising qbar adds deg-6 lines, forcing this compensating +disclination
knot. Appears in BOTH push directions and in C15 & R = the universal elementary
illegal defect.

**CERTIFICATION infra (built this session, in scratchpad):** proper two-sided
equilibration gate reusing the repo's convergence primitives
(discrete_differential_geometry.quantized_split_rhat / integrated_autocorrelation_time),
same criteria as equilibrium_vdv (qRhat<1.05, min-ESS>=100) + two-sided gap +
tau_growth/min_cross blind-spot checks. Scripts: r_knot_chain.py (one chain, logs
illegal-component size time series; START arg for two-sided, native/facet-target
from perfect CELL, VDV off), cert_analyze.py (the gate). USE THIS PATTERN for any
"certified ensemble" ask — do NOT hand-roll R-hat. [[vdv-hdv-conflicts-with-tcp]]

**BEHAVIOR so far:** n_5knots is a FAST mode — CERTIFIED (qRhat 1.002, ESS 382,
two-sided gap 0.02), turns over on tau~250 sweeps. Total illegal count is a SLOW
metastable mode (larger 12-24-atom merged components). The illegal defects form a
SIZE SPECTRUM (5,7,10-26); the (3,4,4) 5-knot is the ELEMENTARY member. "Isolated
5-knots only" is NOT tunable — dilute=few, denser=they merge. cimp=1.0 is the
sweet spot (cimp>=1.5 pushes budget into spread-out low-m tangles AWAY from the
compact high-m 5-knot). 5-knots are TRANSIENT so counts need the TIME SERIES.

**5-HOUR PRODUCTION RUN: first launch (task b34vk8bwi, 04:04) was KILLED ~12min in
by machine idle-sleep. RELAUNCHED under `caffeinate -i` 2026-07-21 09:29 (done
~14:25), task bfynduql6 — PreventUserIdleSystemSleep assertion confirmed active.**
Partial-run preliminary result (kept): illegal-defect SIZE SPECTRUM has PREFERRED
sizes (peaks at 5, 11-12, 14); the (3,4,4) 5-knot is the single most common
component (20%). Not certified (only ~30 samples). Full run params below.
Regime: R m4 (V=10176, mcell=4), bump +8e-4, cimp=1.0, zleg=0.6, nhinge=6, VDV/HDV
OFF, no tilt. 8 TWO-SIDED chains: below0-3 (start perfect r_m4.mfd, illegal rises
from 0), above10-13 (start r_m4_over.mfd = 289 illegal, falls). burn 2500, sample
every 150 sw, snap every 5000, TLIMIT 17700s. Data in scratchpad/run5h/ (TEMP —
copy to data/ if keeping): {below,above}SEED.ts.jsonl (per-sample: sizes, n_5knots,
per-component members[manifold labels]+cocycle centroid), .events.bin (full move
stream, dtype move_geometry.EVENT_DTYPE), .ledger.npz (cumulative geometry ledger),
_snapN.mfd + _snapN.cocycle.npz, _final.*. All logging validated ~0% overhead,
cocycle closed+no-drift on R. Scripts: r_knot_run.py, run5h.sh, cert_analyze.py.

**DEEP ANALYSIS of the salvaged ~20k-sweep data (run killed twice by the HARNESS
terminating long background tasks — NOT sleep; lid open + caffeinate active. 8
chains x ~117 samples + 4 snapshots each + event logs). KEY RESULTS:**
- STRUCTURE: illegal defects are organized by degree-3-edge CORE COUNT (the deg-3
  edge = elementary +disclination charge). 5-knot=1 core (3,4,4) [32/32], size-12
  =2 cores (3,3,4...), size-14=3 cores (3,3,3,4...). Bound states share a deg-4
  screening cloud. SEPARATE deg-4-only family at size 7-8 (no core).
- DYNAMICS (2 populations): the 1-core 5-KNOT is FLEETING (99.6% flicker, gone in
  <150 sweeps, rarely accretes) — sampling too coarse to resolve its lifetime.
  The 2-3 core BOUND STATES are PERSISTENT (49% 2->2, 39% 3->3 per sample) and
  slowly interconvert by accreting/emitting single cores (2<->3 ~5%). High
  turnover (~2 births+2 deaths / 150 sweeps). Defects largely PINNED (median
  inter-sample centroid disp = 0), NOT freely diffusing. The bound states are the
  SLOW mode = why total-illegal didn't certify.
- CERTIFICATION (81 samples, trim<8000): n_5knots qRhat 1.043 ESS 561 (borderline,
  two-sided gap 0.56); n_illegal qRhat 1.547, n_components 1.293 = NOT converged
  (slow multi-core modes). Spectrum well-sampled: peaks 5,11,12,14; 5-knot 21%.
- Analysis scripts: cert_analyze.py + inline dissection. Structure/dynamics done
  from .ts.jsonl (comps: members+centroid) and _snap*.mfd. Cocycle centroids are
  tree-lift units (period ~few x1e6).
- CLEAN FOLLOW-UP (no new run): EVENT-LOG REPLAY (run5h/*.events.bin, dtype
  move_geometry.EVENT_DTYPE) to reconstruct core trajectories at MOVE resolution
  -> separate "knot died" vs "knot moved", real lifetimes + diffusion coeff.
- HYPERUNIFORMITY (via scripts/sk_torus.py = curvature structure factor S(k) on
  T3 through cocycle harmonic coords, field=n6 disclination density, vs exact
  permutation null; ratio<1 = hyperuniform). Ran on run5h snapshots (need .mfd +
  sibling .cocycle.npz, both saved). RESULT: perfect R crystal S(k)/null=0.000
  (perfectly hyperuniform); DOPED R equilibrium = 0.194±0.017 low-k (16 snaps) =
  STRONGLY hyperuniform, ~5x below random -- defects degrade it modestly, don't
  destroy it. But NEARLY not perfectly HU: ratio lowest at intermediate k
  (~0.07 @ k=5) and RISES toward longest wavelength (0.23 @ k=1) = defect-induced
  residual long-wavelength fluctuation, consistent with defects being PINNED
  correlated bound states (not a HU-dispersed gas). Ties to the project's
  Hamiltonian-constraint / curvature-HU theme. Also try --field deficit (only
  sees impurities). Built r_m4.cocycle.npz for the perfect baseline.
  VALIDATED with two controls (metric behaves correctly): perfect crystal =
  0.000 (S(k)=0 between Bragg peaks, Bragg peak at k=6); MELT (r_m3, FK potential
  OFF, fFK->0) = 2.14 (>1 at ALL k = clustered/hyper-fluctuating, NOT HU). So the
  3-way: crystal 0.00 / doped 0.19 / melt 2.14 — doped sits next to crystal, far
  from melt => doped HU is genuine. Melt script melt_r.py (no n6 potential melts
  instantly, fFK->0 in <500 sweeps); cocycle survives melting (check OK).
  BEST FIELD = q_R = 1/2 sum_e delta_e (Regge scalar curvature, delta=2pi-theta*deg,
  theta=arccos(1/3); flat edge degree 2pi/theta=5.1043 = native TCP, NOT 6) = the ADM
  Hamiltonian-constraint quantity (spatial scalar curvature). Already implemented as
  q_R in hyperuniformity.py + curvature_hyperuniformity_g.py (GRAPH estimators: BFS
  window-variance + Laplacian spectral vs shuffle null), and volume-normalized
  R=sum delta/(2*D_v/4) in adm_constraint.py (which also builds trK from volume_flux
  + tests R=a+lambda*K2). sk_torus (real-position S(k)) only has n6/deficit -- q_R
  NOT wired in (~2-line add). deficit=(6-deg) is DEGENERATE (Euler const 12); physical
  weight is (2pi/theta-deg)=(5.1043-deg). RESULT (graph q_R HU via qr_hu.py copying
  hyperuniformity.py estimators): perfect crystal window 0.18/spectral 0.000 (perfectly
  HU); DOPED window 0.19/spectral 0.001-0.11 (strongly HU, ~=crystal!); MELT window
  19.9 grows w/ m / spectral 1.7-28 (clustered, NOT HU). So spatial scalar curvature
  is hyperuniform in the doped samples, nearly at crystal level. CAVEAT: hyperuniformity.py
  runs its ensemble loop AT IMPORT (hardcoded seeds/ S3 files) = not cleanly importable
  = a HARMONIZATION target (Aaron wants to unify the HU scripts: hyperuniformity.py graph
  + sk_torus positions + curvature_hyperuniformity_g + adm_constraint, all share q_R/delta).

**(superseded) ANALYSIS TODO when run5h_done.flag appears:**
1. CERTIFY: two-sided gate (cert_analyze.py pattern, quantized_split_rhat<1.05,
   ESS>=100, gap, tau_growth/min_cross) on n_5knots/n_illegal/n_components from
   .ts.jsonl; trim initial transient (the above side falling from 289).
2. SIZE SPECTRUM: distribution of illegal-component sizes.
3. STRUCTURE: dissect knots from _snap*.mfd -> (3,4,4)-core stats (already 11/11).
4. MOBILITY: link 5-knots across samples by member-overlap (manifold labels stable
   across moves) -> trajectories -> centroid MSD -> DIFFUSION COEFF. Centroid is
   cocycle tree-lift (clean period jumps). Event log = move-by-move mechanism.
5. DISCLINATION context: disclination.py census (n_six, segments, loops, fray).
