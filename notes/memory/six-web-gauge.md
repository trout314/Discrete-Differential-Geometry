---
name: six-web-gauge
description: "Emergent-gauge reading of the six-edge web MEASURED — nematic L/T=0.60±0.04 frozen; vector channel (web_vector.py) = exact Gauss-law identity, crystal L/T~1e-26, monopole charges anticorrelated"
metadata: 
  node_type: memory
  type: project
  originSessionId: 1ed7f5d0-df5c-4552-b65e-5060180eafe2
  modified: 2026-07-22T22:03:05.902Z
---

Reframe (2026-07-22, user-driven): the six-edge disclination web = emergent
gauge field; FK class condition n6 in {0,2,3,4} = lattice GAUSS LAW (closed
lines: 0/2/3/4 = none/passthrough/Y/X; n6=1 = line ENDPOINT = MONOPOLE);
zleg/cimp = monopole mass (penalty, NOT constraint); knots/blinkers = monopole
plasma. Direct mapping to spin-ice/dimer Coulomb phases: exact constraint =>
entropy-generated algebraic correlations + pinch points; soft constraint =>
monopole screening with width ~ sqrt(n_mono). On T3 the closed web carries a
Z^3 HOMOLOGY winding W — exactly conserved by local moves on the legal
manifold (the clean topological conservation law; monopoles let W drift).

MEASURED (scripts/defect_dynamics/web_s66.py: six-edge director-tensor
structure factor, TT/V/L helicity split vs occupation-shuffle null):
- crystal: stealthy (0.000 sub-Bragg).
- L/T at k=1-2, pooled: FROZEN lam1 = 0.603+/-0.040 (10 sigma below 1!);
  liquid lam0.4 = 0.844+/-0.064; lam0.35 (186 mono) ~1; lam0.30 web (3642
  mono) ~1 fully screened. Longitudinal suppression = smeared pinch point =
  emergent Gauss law CONFIRMED at low monopole density, k-trend right
  (strongest at smallest k). Stiffer constraint (lam1) suppresses better than
  the liquid at equal monopole COUNT — violation ACTIVITY (flicker), not just
  census, sets the screening.

EDQ-only ladder (aesthetic-action test, n6 OFF): lam_EDQ 0.4 -> crumple-melt
(<edeg> runs to 5.62!); 0.8 -> dense disclination fluid legalvert 0.06. R^2
term alone does NOT hold the constraint manifold => the six-web closure term
(decoded n6) is a NECESSARY third action term: S = volume pin + R^2(EDQ) +
line-closure(Bianchi). Stiff rungs 1.2/1.6/2.4 pending.

Bump ladder: carrier density is BUMP-INDEPENDENT (8e-4..4e-3 all give nill~
50-63): thermal, not charge-driven; the soft pin's demand goes unsatisfied
(15% at lam0.4, 66% at lam1 — certified two-sided cert states 5.104755 both
sides). Penalties don't make conservation laws — neutrality NOT enforced.

VECTOR CHANNEL (scripts/defect_dynamics/web_vector.py, 2026-07-22): the
nematic L/T=0.60 is a smeared inequality (Q quadratic in nhat; bending leaks
into k.Q.k). The identity-grade observable: orient the web by CYCLE-SPACE
PROJECTION (phi = sigma0 - B^T mu, BB^T mu = B sigma0, LSMR; div allowed only
at imp>0 monopole vertices — Y-junctions carry fractional flux, the abelian
shadow of the non-abelian disclination algebra), Fourier with the EXACT
SEGMENT INTEGRAL kernel e^{ikx0} e^{i th/2} sinc(th/2pi), th = k.d. Then k.J
TELESCOPES to (1/i) sum_v D_v e^{ikx_v} => L_vec(k) = S_mono(k)/k^2 exactly.
MEASURED: identity verified to 1e-10 abs on 1e10 powers; crystal L/T = 3e-26
(exact zero — upgrades nematic 0.60); monopole states L/T ~ 1e-2 purely
monopole-sourced; AND S_mono(k)/sum(D^2) = 0.10-0.13 at smallest shell in ALL
monopole states (rising to ~0.5 at large k): monopole charges themselves
anticorrelated at long wavelength (SL-type organization of the violation
plasma). Summary doc: notes/two_term_action_and_emergent_gauge.tex (also
covers fullerene theorem => two-term action verdict, EDQ ladder, S66,
penalties-vs-constraints, block moves, frozen sampler).

NEXT (agreed direction): drive monopole density -> 0 while keeping web mobile:
(1) budgeted-illegality sampler — reject moves pushing global illegal count
above a small cap (worm-like; trivial on the frozen-vertex rejection pattern,
n6 potState already counts); (2) watch L/T -> 0 and pinch sharpen as the
budget tightens; (3) W-sector physics on the legal manifold. lam0.4 ensemble
NOT yet certified (qR-hat 1.062, ESS 83, one-sided init) — two-sided cert
campaign queued (5 extends + 3 above-side chains from run5h states; mobile_gas
now has resume args start/startcoc).

See [[mobile-gas-liquid]], [[fk-move-search]], [[defect-kinetics-run5h]].
