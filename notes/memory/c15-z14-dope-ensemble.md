---
name: c15-z14-dope-ensemble
description: Equilibrated C15+Z14-tilt doping ensemble — working regime + the Z15-pair finding
metadata: 
  node_type: memory
  type: project
  originSessionId: 1ed7f5d0-df5c-4552-b65e-5060180eafe2
  modified: 2026-07-21T05:24:35.472Z
---

2026-07-21: revisited the "Z14-Z15 dopant complex in dynamic C15" picture with the
new crystal_grains defect identifier, on a freshly built **C15 m6 cell** (V=5184,
native qbar=5.100000 exact, Z12:3456 Z16:1728). See [[crystal-grains-tool]].

**Key theory (derived this session):** mean edge degree is quantized,
qbar = 6·f3/(f0+f3), and the finest step (one 2->3 Pachner move) is
Δqbar_min = (6-qbar)²/(6V) ≈ 0.135/V. For C15 m6 that's 2.6e-5; excess-disclination
DENSITY = 7.41·Δqbar_held (V-independent), so dilution needs a small *held* bump and
a big enough cell that the few complexes sit isolated instead of percolating (the m4
V=1536 cell percolated). Bigger cell => finer quantum, better isolation.

**WORKING REGIME (dilute, isolated, equilibrated):** start from perfect C15 m6;
SamplerParams(num_facets_coef=0.1, hinge_degree_target=native+8e-4, num_hinges_coef=6.0,
hinge_degree_target_coef=1.0*edge_target/6); set_n6_potential(zleg=0.6, cimp=1.0,
tilt=[0,0,-1.5,0,0]) i.e. mu=1.5 toward Z14. ~1500-2000 sweeps equilibrates
(flat from ~500). nhinge=6 only HELD ~10 of the 31 requested quanta (qbar reached
native+2.6e-4) — the FK potential resists; honest "held elevated edge degree".
Tuning lesson: strong cimp (>=~2) FREEZES nucleation (needs illegal intermediates);
weak cimp percolates illegals; cimp~1.0 is the sweet spot. Scripts in session
scratchpad: c15_dope_run.py (start seed mu bump total chunk zleg cimp nhinge out),
run_ensemble.sh (8 seeds parallel), c15_m6.mfd (the cell).

**RESULT — 8 replicas, tight equilibrium:** cryst 98.87±0.11%, 4.2±1.1 isolated
complexes, Z14 3.5±1.7, Z15 9.9±1.4, Z14-Z15 adjacencies 2.1±1.2 per replica.
The DOMINANT elementary excitation is a **Z15 PAIR** — recurring size-9 complex
`2 Z15 + 7 illegal-halo` — NOT a lone Z14. Z14 appears only in LARGER mixed
complexes (15(1·2·12), 24(2·4·18), 40(4·9·26)), where it pairs with Z15. So
Z15:Z14 ≈ 3:1 despite the Z14 tilt — raising C15's edge degree nucleates Z15
pairs first; the "Z14-Z15 complex" is the grown end of that family. Each complex =
a few FK dopants + a distortion (illegal-edge) halo; crystal_grains groups the
whole thing as one off-registry cluster (physically correct).

**NO-TILT control (mu=0, everything else identical) — the DEFINITIVE probe:**
8 replicas, cryst 99.37±0.10% (cleaner + more dilute than tilted), #cx 2.5±0.5,
Z14 1.12±0.93, Z15 4.25±1.20. **Z15/Z14 = 3.8** (vs 2.8 WITH the Z14 tilt) — so
the Z15 dominance is INTRINSIC and the tilt was fighting it (pumping extra Z14
lowered the ratio). Complex-motif census over 8 replicas, by (Z14,Z15):
(1,2)×7 [the size-15 1·2·12i unit, most common], (0,2)×5 [Z15 pair], (0,1)×4
[lone Z15], (1,3)×2, (0,0)×2 [pure distortion knots]. TWO robust facts:
(1) Z15 is in ~every complex; ZERO pure-Z14 complexes — Z14 only ever appears
FLANKED BY Z15 (canonical grown unit = 1 Z14 + 2 Z15). (2) the elementary carrier
is a lone Z15 or a Z15 PAIR; the Z14-Z15 object is the next member up the family,
not the ground-state excitation. So the original "Z14 paired with Z15" picture is
real but is the grown end; the intrinsic C15 edge-degree excitation is Z15-centric.

**RAISE vs LOWER ASYMMETRY (both no-tilt):** targeting a LOWER edge degree
(bump -1e-3, 25% bigger push) produces essentially NOTHING — all 8 replicas land
on the IDENTICAL state: 99.90% cryst, native Z16=1728, Z14=Z15=0, one 5-atom
pure-illegal knot (0·0·5i). And it's ~3x STIFFER: -1e-3 target realized only
-6.9e-5 deficit vs +8e-4 realizing +2.1e-4. So C15's edge degree is SOFT UPWARD
(nucleates Z15-centric complexes), effectively INCOMPRESSIBLE DOWNWARD. Reason:
C15 sits near the FLOOR of the TCP mean-edge-degree range (CN 13.33, among lowest
FK phases). Lowering coordination = REMOVING deg-6 edges from the major
disclination skeleton, but those lines are topologically continuous — excising a
segment leaves a dangling illegal disclination that cimp forbids. Raising = ADDING
closed skeleton loops, always available. The unique dopant-free downward
equilibrium (8/8 identical) underscores it.

**R CRYSTAL (m3, V=4293, native qbar=5.1042, census Z12:2187 Z14:972 Z15:486
Z16:648 — SATURATED, all 4 Kasper classes) same 2 experiments, no tilt:**
- LOWER -1e-3: PERFECTLY RIGID — 8/8 replicas 100.00% cryst, 0 complexes, 0
  illegal, exactly native. Even stiffer than C15 down.
- RAISE +8e-4: marginally soft but HETEROGENEOUS (cryst 88.3-100%, ±3.68 — near a
  nucleation threshold, not a clean equilibrium). Realized excess only +3e-5
  (~7x stiffer than C15's +2e-4). KEY: because R already has all species, raising
  makes OFF-REGISTRY FK-CLEAN patches, NOT new dopants — recurring motif
  46(26 Z12·0 Z14·12 Z15·8 Z16 | ZERO illegal): 46 legal shells in the WRONG
  registry. seed4 outlier: 501-atom off-registry region, 98% clean FK (10 ill).
  These are exactly what the OLD local tools call crystalline; only the
  covering-map registry check flags them — a clean demo of crystal_grains.
  Per-complex illegal: raise 0/5.2/10 (4/10 pure-illegal 5-atom knots); lower 0.

**R RAISE pushed harder (+2e-3, 2.5x): crosses a THRESHOLD into disorder.**
Response becomes UNIFORM (was near-threshold/heterogeneous at +8e-4): cryst
94.4±4.3%, 5.9 complexes/replica, realization jumps 4%->35%. CHARACTER CHANGES
from clean off-registry FK patches to ILLEGAL-LADEN disorder: illegal fraction of
defect mass 7% -> 33%; recurring motif 41(16·3·3·5|14ill) in 6/8 replicas; census
Z14 & Z16 both drop, extra coordination dumps into illegal degree-7 edges not
clean higher-FK. seed0 nucleated a 676-vertex partial-disorder region (83% cryst).
So R (saturated) has NO clean soft mode — forced to yield it frays into local MELT
(illegal proliferation), unlike C15 which makes clean new Z15-pair dopants. This
is a disordering transition, not a dilute equilibrium.

**FK-PUSH TURNED DOWN (zleg 0.6->0.3, cimp 1.0->0.4), R both directions:**
splits the two directions cleanly. RAISE +8e-4 still LARGE MOSTLY-FK (illegal-frac
only 20%, sizes to 709) — raising has a CLEAN FK CHANNEL (excess coordination ->
higher-Z off-registry patches, FK-legal, so soft cimp just lets R yield bigger).
LOWER -1e-3 becomes SMALL PURE-ILLEGAL: 99.7% cryst, every complex a
5(0·0·0·0|5ill) 5-atom pure-illegal knot — lowering has NO clean FK channel (must
break skeleton = illegal), and at strong cimp it was forbidden (rigid, nothing),
soft cimp UNLOCKS these. So 5(0·0·0·0|5ill) is the ELEMENTARY ILLEGAL DEFECT of
these triangulations (same knot appeared for C15 lower). To get "small illegal
defects not large FK ones", the LOWER direction + soft cimp is the regime; RAISE
inherently prefers FK-clean patches.

GENERAL: both C15 and R are SOFT to raising edge degree, STIFF to lowering
(add closed disclination loops = easy; excise continuous skeleton = illegal =
forbidden). C15 (binary Z12/Z16) has an easy soft mode -> new Z15-pair dopants
(illegal-halo-dressed). R (saturated) has no soft mode -> what it yields is
registry defects (0-illegal clean-FK patches). Illegal-per-complex separates the
two physics.

Ensemble states in session scratchpad (gitignored/temp): C15 rep*/nt*/lo*.mfd,
R rhi*/rlo*.mfd + .json; run scripts run_ensemble*.sh, run_R.sh, c15_dope_run.py
(now takes phase arg 11). Next: bigger negative bump, mobility, dissect motifs,
promote c15_dope_run.py to scripts/ if keeping this line.
