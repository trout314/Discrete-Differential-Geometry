---
name: fpkmc-m1-status
description: "M0 DONE (all symbolic proofs pass), M1 core LANDED (chain_walk + site_survey in D, 1ms/site, cross-validated 1e-9); V1 census OVERTURNS chain picture — slide channel lives on a degree-12 graph, most chain transport is species-changing"
metadata: 
  node_type: memory
  type: project
  originSessionId: d562feea-0887-419b-9cd9-9899d180c7c2
  modified: 2026-07-26T23:12:12.451Z
---

2026-07-26. M0 DONE: `scripts/defect_dynamics/fpkmc_m0_derivations.py`
(SymPy exact, ALL PASS): telescoping identity, detailed balance |I1, cycle
condition, splitting closed form, discrete-time invariance (holding probs
drop out of splitting but NOT times — nu enters time bookkeeping only),
w-equation, phase-type tail. Design doc §3.2 updated.

M1 CORE LANDED (uncommitted at checkpoint time — commit soon):
- `manifold.d`: nextChainWindow + unittest (5-cycle on standardSphere!3;
  meson test 12/12).
- `ddg_capi.d`: ddg_manifold_chain_walk; ddg_sampler_site_survey — per
  window: create knot via speculative-delta path, trial all 12 slots
  (dS + destination + CLEAN flag via second cleanOnly trial), undo,
  restoration audited in-function. SURVEY_STRIDE = 1+4*12.
- Python: Manifold.chain_walk, ManifoldSampler.site_survey.
- Cross-validated: chain_walk == bc_orbit; dS_create == worm_slide oracle
  (0.1 volume-pin + lam*dS_between) to 1e-9; slot dS == slide_at trials.
- Performance: ~1 ms/site; full R m4 orbit (3252 sites) in 3.4 s.

**V1 CENSUS (the load-bearing result):** R m4 orbit, all sites:
- 12 legal slots per site everywhere: exactly 1 chain-fwd + 1 chain-bwd
  (I1 exact, 0 violations) + 10 off-chain.
- CLEAN subclass: only 900/3252 chain slides clean — most chain transport
  is SPECIES-CHANGING (rung ladder = multiset changes; "dirty" by the
  cleanliness test). Clean off-chain 0-9/site (median 3), dS sym about 0.
- CONSEQUENCE: physical slide channel (dirty-allowed, project default)
  lives on the full degree-12 slide GRAPH with species change along
  edges; sparse-graph formulation is the MAIN path; 1D tridiagonal
  survives only for the clean subclass. Design doc M1 section updated.

M1 COMPLETE (same day, later): edge-symmetry 24/24 verified;
fpkmc.py landed (guard holds, CleanKnotGraph with TWO exact zeros:
antisymmetry + S1-consistency => pure-knot graph is a potential
landscape); clean components are CLOSED ~24-node islands => transport is
inherently species-changing. D graph scan landed:
slideDecode/slideApplyKeep/slideRollback in sampler.d (SlideRec now
public) + ddg_sampler_slide_graph_scan (full legal graph from current
state, exact overlay keys, dS_max pruning, restoration audited) +
ManifoldSampler.slide_graph_scan. First scans: 185 nodes/0.33s at
(dS_max 4, depth 3); 7 species sigs; all single-chord; found states
0.48 BELOW the rung-70 knot (bare knot is NOT the single-defect ground
state at lam=0.4). GOTCHA: meson compile piped through grep eats the
real rc — check compiles directly.

V2 DONE: nu VERIFIED at 0.974+-0.075 of slide_prob*(3/N3)*(1/6)*(1/12)
via the slide_prob=1.0 trick (slide branch swallows every attempt ->
thermal never runs, knot immortal, zero contamination; 2e7 attempts).
GOTCHA logged: in-thermal measurements read 3.87x — flicker chords, and
even a 1.2-cell frozen patch holds ~1100 unfrozen vertices that brew
flicker beside the knot. Scan now emits per-node representative chords.

HB DRIVER LANDED (commit d771d21): fpkmc.HBDriver, Metropolized-ball
independence sampler; walk-back made DETERMINISTIC via slide_at2 (new C
API reporting each slot's arrival chord from the frame decode — the
first dS-guessing walk-back hit a same-dS impostor within 60 steps).
Smoke test: 80 steps, 31 accepts, ALL running-action audits pass at
1e-6 (~1.8 s/step, two ball scans).

V3 (scripts/defect_dynamics/fpkmc_v3_hb.py): HB 1500 steps vs
brute-force slide-channel emulation (uniform facet/pair/slot +
Metropolis via slide_at trial + force; no thermal channel in-process so
the defect is immortal) vs per-state Boltzmann. PASS = pairwise
ln-occupancy-ratio vs -dS slope 1 for both chains. States keyed
(chord, round(S_rel,6)); both chains carry audited running actions.

V3 OOM SAGA (2026-07-26): killed TWICE by jetsam. (1) unbatched run at
step ~800 (~25 MB/step D-GC scan leak, GC.collect ineffective, dominant
allocator = overlayKey per-visit baseline.keys). (2) exec-restart
batching at --hb-batch 250 STILL killed 100 steps into batch 3 (6 GB
peak/batch too much beside 3.1/4 GB swap used); worse, the .state files
saved in the session SCRATCHPAD were deleted by something external
(machine up 20 days — not a reboot wipe) => 500 steps lost. LESSONS:
batch <= 60 steps (~1.5 GB peak); resume state goes in repo data/, NEVER
/private/tmp. THIRD kill: brute phase's 4M raw slide_at2 calls -> ~4M D
exception throws (each allocates exception+trace under the GC) ->
jetsam; FIXED by Python-side deg-3 edge prefilter (same proposal dist,
same rng stream) + HB histogram persisted to .hb.json before brute
starts. PROPER FIX still queued (R6): allocation-free scratch in
ddg_sampler_slide_graph_scan.

V3 FIRST COMPLETE RUN (lam=0.40, data/fpkmc/v3_hb.json): infrastructure
clean (all audits pass) but occupancy test VACUOUS at lam=0.40 — knot
too mobile: HB 682 states/1500 steps, brute 62 states/91 moves, ~no
revisits => chains transient, per-state occupancy measures escape rates
+ translational entropy, NOT e^{-S} (slopes ~0, and S=+4 states still
got >=5 visits — impossible in equilibrium, expected for a wanderer).
Confinement rerun at lam=1.0 ALSO vacuous (757 states — the ball
sampler tunnels the washboard; chain glides on 3 near-degenerate levels
S in {0,0.6,1.2}): occupancy validation structurally dead here.

V3b PASSED — M2 VALIDATION CLOSED (2026-07-26):
scripts/defect_dynamics/fpkmc_v3b_db.py = numerical detailed-balance
certification of the HB kernel, NO equilibration needed. Both DB sides
reduce to min(e^{-d}/Zx, 1/Zy) => test the 4 implementation facts:
A1 dS antisymmetry (5.9e-14 max), A2 membership asymmetry safely
rejected (alpha=0 kills both fluxes; 22-32/60 uniform draws hit it),
A3 restoration determinism (node count exact; sorted dS to ~2e-14 —
scan dS has PATH-DEPENDENT float accumulation, DFS order differs after
restoration, so bitwise Z equality is unattainable; gate = count +
sorted-dS 1e-9), A4 scan-tree path-validity symmetry (the REAL DB risk:
fwd-valid/rev-invalid one-way flux; 0/237 pairs). ln-flux residual
<= 4.6e-14. 4 seeds x 60 pairs, lam=0.4, m2. Design doc section 9
updated (V3 lesson + V3b entry). Kernel obeys DB to ~5e-14/transition.
V3/V3b COMMITTED+PUSHED (b17b213).

M3 IN PROGRESS (2026-07-26): FP frozen driver on the GRAPH (protective-
domain construction, not the 1D chain — V1 census). D: slide_graph_scan
gained optional trailing params (blocked_verts, n_blocked, node_dock) —
dock = knot complex touches ONE-TET-LAYER nbhd of blocked set (static,
computed at root state; state-independent boundary => exact FP for any
dynamics absorbing on the same rule); dock nodes enumerated, never
expanded; null blocked => byte-identical HB behavior (V3b cert intact).
COMPLETE-INTERIOR guarantee: node with depth<max_depth & dS<=dS_max &
dock==0 was fully expanded (fresh children ALWAYS added as nodes; only
expansion is gated). Python: fpkmc.FPFlight (classify interior/absorb
by dock/dS_frontier/depth_frontier/multichord; exact jump-chain sample
with geometric holding in attempted-move units; dense exact splitting +
mean-time solves) + fpkmc.FPDriver (flights, materialize exit via scan-
tree replay — tree intermediates are always interior). V4 harness
scripts/defect_dynamics/fpkmc_v4_fp.py: two knots on one orbit, B
frozen (its chord's proposals excluded = v1 FROZEN), brute = exact
marginalized law (geometric attempts-to-hit 3/(6N3), uniform slot,
Metropolis via slide_at2), state matched to scan nodes by
(chord, round(S,6)) with collision guard, walk-back + restore audits.
V4 PASS BOTH CONFIGS (2026-07-26, data/fpkmc/v4_*.json): frontier
(sep5/depth3, 38 interior): chi2 p 0.879, KS p 0.703; dock (sep3/
depth4, 72 interior + 23 dock, P(dock)=0.0415): chi2 p 0.578, KS p
0.492. 0/600 unmatched landings — complete-interior guarantee exact.
FP speedup: ~6e4 attempts collapsed to ~10-100 jump draws per flight.
V4 batch COMMITTED+PUSHED (a12409d).

V5 PASS BOTH SEEDS (2026-07-26, fpkmc_v5_contact.py, data/fpkmc/
v5_*.json): Part 1 head-on march: V=0 to 1e-9 at every non-dock sep;
dock and interaction onset COINCIDE at sep 1 (boundary TIGHT not
conservative — still exact: interior additivity is the license, dock
nodes carry exact energies). Contact V=+1.76 = 0.40*4.4 = collider
phase-1 rung value via independent path. Part 2: 9/12 FP episodes
docked, |V(dock)| <= 1.1e-12 at arbitrary dock geometries, all bursts
BOUNCE (repulsive+reversible, matches collider), one burst found the
-0.48 sub-knot level; every episode restored exactly. FPDriver now
records recs for walk-back; slide_at2(commit=True) used in materialize.
M3 CORE DONE (driver + V4 + V5), V5 batch COMMITTED+PUSHED (3c135e8).

FP PRODUCTION RUN (2026-07-26, ~3h, task bl732mddg, log
data/fpkmc/prod.log): Phase A encounter kinetics tau(s0), s0 in
{2,3,4,5,6} x 3 procs x 4 episodes (fp_encounter.py --mode encounter;
censored-at-60-flights episodes recorded); Phase C lone-knot slide-only
transport, 8 procs x 200 flights (fp_transport.py; positions =
minimal-image-unwrapped chord midpoints in cells via
reference_frac_positions("r",2)); Phase B recombination from sep 3,
5 procs x 4 episodes (unbind burst then escape-vs-recombine, escape =
8 dock-free flights). Outputs data/fpkmc/prodA_s*_p*.json, prodC_p*,
prodB_p*.

3 waves DONE clean (21+74+24 min — flights ~0.2s, 8x faster than
estimated; waves 2-3 added s0=8,10, 120-flight de-censoring, lam=1.0
transport, sep-2 recombination). RESULTS (fp_prod_report.py; figure
data/fpkmc/fp_production_report.png, sent to user):
- tau(s0) KM medians: 7.3e5 (s0=3) -> ~5e6 attempts (s0=6-8), DROPS at
  s0=10 (2.7e6) = finite-box wrap (m2 = 2 cells across; chain-site sep
  != spatial sep at large s0). Real tau(r) needs m4 host.
- Dock census 253 docks: ~83% OFF-CHAIN (1D head-on is a minority
  channel); 251/253 additive (V < 1e-6), 2 contact docks (V = +4.5 at
  s0=2; flights can absorb one step INTO the wall — still exact).
- TRANSPORT SURPRISE: D_slide is lam-INDEPENDENT: 6.6e-7 (lam 0.4, 32
  traj) vs 6.3e-7 cells^2/attempt (lam 1.0, 20 traj); exponents
  0.88/0.90. Mechanism: transport rides the dS=0 flat directions
  (degenerate translates), lam-independent — washboard barely matters.
  (The earlier 2x suppression at 8 traj was noise.) = 4.8e-3
  cells^2/sweep; MSD short-lag plateau ~0.4 cells^2 (ball scale),
  ~linear beyond 3e5 attempts.
- Recombination: P(rec|freed) = 0.73 from sep 3 (44/60), 0.63 from
  sep 2 (17/27); t_unbind median 6.1e4, t_return median 8.5e4, return
  tail to 1e7 attempts.
DOCK-ANGLE ANALYSIS (fp_dock_angles.py, figure fp_dock_angles.png,
sent): rod axis at a chord = chain_walk +-8 windows (2 precession
periods) end-to-end displacement in reference registry. DISCOVERY: the
BC walk direction PRECESSES with period 8 windows (autocorr +0.90 at
lag 8, +0.73 at 16, decays to 0.16 by lag 48) — chains are coiled
helices with finite persistence ~tens of windows, NOT straight rods
(the closed 3252-window orbit in a 2-cell box already implied this;
"on-chain" does NOT mean parallel). Dock-angle spectrum ~ ISOTROPIC
sin-theta (mild 45-60 excess / 60-75 deficit); no strong angular
selectivity in encounters. P(rec|freed) by angle: aligned docks
(<30 deg) 13/15 = 0.87 vs crossed 48/72 = 0.67 — suggestive (~1.6
sigma), needs more episodes. Both contact docks at large angle
(70/75 deg).
NOTE manifold.link/star ALLOCATE full facet list per call (R6 target).
UNCOMMITTED: fp_encounter.py, fp_transport.py, fp_prod_report.py,
fp_dock_angles.py. See [[fpkmc-design]].
