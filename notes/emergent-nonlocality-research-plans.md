# Emergent Nonlocality from Curvature Condensates — Research Plans

*Speculative program connecting the hub/de-hubbing result on edge-pinned 3D DT
triangulations to emergent-locality and quantum-nonlocality ideas. Drafted
2026-07-11. See memory `project-vertex-degree-structure` for the underlying
empirical result.*

## Background result (what motivates all of this)

The k=2 flat-edge family (`1e-1_ED5p1043_2`, β=0) decomposes into two intertwined
geometric pieces on the same combinatorial S³:

- a **flat, extended, ~N^(1/3)-diameter bulk skeleton** (the de-hubbed 1-skeleton;
  diameter ∝ (giant size)^0.35 ≈ size^(1/3) after removing degree > ~5·flat), and
- a **connected, negatively-curved, small-diameter hub core** (a "rich club":
  vertices of degree > c·flat form a single connected component, each wired to
  ~16–18 other hubs, max degree ∝ √N, degree-tail exponent α_deg ≈ 3).

The hub core short-circuits the bulk into a small-world graph; removing it exposes
the extended 3D metric. Locality is therefore *metric-relative*: one object, two
geometries. The speculation is that this is a toy model of how quantum
nonlocality could be the shadow of hidden graph connectivity on an emergent
low-dimensional space (disordered locality / ER=EPR / Wolfram-hypergraph themes).

---

## Plan A — Rigorous two-metric characterization (do first; publishable alone)

Make "one complex, two geometries" theorem-grade.

1. **Spectral dimension** d_s via random-walk return probability p(t) ~ t^(−d_s/2),
   for (i) intact graph, (ii) de-hubbed bulk, (iii) hub core alone. Predict bulk
   d_s → 3; core d_s large/ill-defined; intact graph scale-dependent d_s (compare
   CDT dimensional-reduction curves).
2. **Volume growth** B(r) per piece: bulk ~ r³ with sin²-type shell (tie to the
   roundness pipeline); core ~ exponential (hyperbolic). Optionally Gromov
   δ-hyperbolicity of the core.
3. **Degree-corrected rich-club null:** Maslov–Sneppen degree-preserving rewiring
   to get the *enrichment* of hub–hub connectivity above chance (fixes the current
   raw "16–18 edges/hub" claim).
4. **Interface geometry:** dimension of the attachment set (bulk vertices adjacent
   to the core) — surface (d=2, horizon-like), space-filling (d=3), or punctures
   (d=0). Determines which physics analogy is available.

Success: scale-dependent d_s(t) with clean plateaus + degree-corrected rich-club
coefficient ≫ 1.

## Plan B — Disordered locality: field correlations on the two-metric graph (do second; cheap)

Markopoulou–Smolin "disordered locality" is literally this object; we can measure
their conjectured consequences from a controlled action.

1. Free scalar on the graph, covariance (L + m²)^(−1), L = graph Laplacian.
   Two-point G(x,y) **binned by bulk (de-hubbed) distance**.
2. Predict: G ~ flat 3D propagator e^(−m·d)/d at short bulk distance, then a
   **plateau/floor** from core leakage (correlation between bulk-far, core-near
   pairs). Measure floor vs hub-tail strength (tunable via β/N).
3. Vary m to find the crossover mass m* below which locality fails.
4. Cheaper variant: heat kernel on the same Laplacian.

Deliverable: "propagator anatomy of disordered locality in equilibrium 3D DT."

## Plan C — Correlation architecture: pairwise vs shared-resource

ER=EPR is pairwise; our core is one connected multipartite object. Which
architecture does it induce?

1. From Plan B's field, take k bulk-distant, core-near regions; compute full
   covariance + mutual-information graph.
2. Test **monogamy vs promiscuity**: pairwise-wormhole ⇒ monogamous; shared-register
   (rich club) ⇒ promiscuous, high multi-information, conditioning on the core
   kills all pairwise correlation.
3. Fit a **single-latent-node common-cause** Gaussian model; compare likelihood.

Note: classical, so Bell inequalities are *satisfied* — the point is the
architecture (GHZ-like common cause vs Bell-pair-like), à la Susskind ER=EPR+GHZ.

## Plan D — Dynamical hubs & the measurement-independence question

Static shortcuts = nonlocal hidden variables (inert); interesting regime is a
*responsive* core.

1. **Quench:** equilibrate, locally perturb the objective (pin a small ball of edge
   degrees off-target = crude "setting"), continue MCMC, track hub-core
   reorganization near the perturbation and its timescale (sweeps) vs bulk
   relaxation.
2. **Setting–core mutual information:** zero ⇒ frozen common cause; nonzero + fast ⇒
   responsive medium. Doubles as a curvature-transport measurement.

Warning: likely null (core reorganizes on the slow VDV timescale). A principled
null result is the honest way to kill the QM speculation.

## Plan E — Universality & continuum limit (prerequisite clincher)

Extend `1e-1_ED5p1043_2` to N ≥ 1e5 (currently ≤1e4). N-adjusted de-hub cutoff;
nail bulk exponent (→1/3?) and d_s → 3; test whether α_deg ≈ 3 and core
connectivity are **universal** across edge-degree targets (ED5p0043 / 5p1043 /
5p2043) and VDV couplings, or tuned. Re-run A–D at the largest N that certifies.

---

## Reading list (for a mathematician: discrete geometry / GR / QFT, non-expert in QM interpretation)

**Bell/nonlocality (minimum, skip interpretation lit):**
1. J. S. Bell (1964), "On the Einstein–Podolsky–Rosen paradox," *Physics* 1, 195.
2. A. Fine (1982), "Hidden variables, joint probability, and the Bell inequalities," *PRL* 48, 291. — LHV exists ⇔ a joint distribution exists; Bell as a marginal problem.
3. B. S. Cirel'son (1980), "Quantum generalizations of Bell's inequality," *Lett. Math. Phys.* 4, 93. — the 2√2 bound.
4. Brunner, Cavalcanti, Pironio, Scarani, Wehner (2014), "Bell nonlocality," *Rev. Mod. Phys.* 86, 419 (§II–III).
5. B. Toner (2009), "Monogamy of non-local quantum correlations," *Proc. R. Soc. A* 465, 59.

**Emergent space with nonlocal links (closest existing program):**
6. F. Markopoulou & L. Smolin (2007), "Disordered locality in loop quantum gravity states," *CQG* 24, 3813 (gr-qc/0702044). — the paper for Plan B.
7. Konopka, Markopoulou, Smolin (2006), "Quantum graphity" (hep-th/0611197); Konopka, Markopoulou, Severini (2008), *PRD* 77, 104029.
8. F. Markopoulou & L. Smolin (2004), "Quantum theory from quantum gravity," *PRD* 70, 124029.
9. S. Hossenfelder & T. Palmer (2020), "Rethinking superdeterminism," *Front. Phys.* 8, 139.
10. G. 't Hooft (2016), *The Cellular Automaton Interpretation of QM* (arXiv:1405.1548), Part I.

**Entanglement–geometry (gravity side):**
11. M. Van Raamsdonk (2010), "Building up spacetime with quantum entanglement," *GRG* 42, 2323.
12. J. Maldacena & L. Susskind (2013), "Cool horizons for entangled black holes," *Fortschr. Phys.* 61, 781. — ER=EPR.
13. L. Susskind (2014), "ER=EPR, GHZ, and the consistency of quantum mechanics" (arXiv:1412.8483). — Plan C.
14. Cao, Carroll, Michalakis (2017), "Space from Hilbert space," *PRD* 95, 024031.
15. (opt.) Swingle (2012), *PRD* 86, 065007; Pastawski, Yoshida, Harlow, Preskill (2015), *JHEP* 06, 149. — hyperbolic tensor networks.

**Dynamical triangulations (where our object sits):**
16. Catterall, Kogut, Renken, Thorleifsson (1996), "Singular vertices and the triangulation space of the D-sphere," *Nucl. Phys. B* 468, 263; Catterall et al. (1998), *Phys. Lett. B* 416, 274. — singular-vertex condensation = ancestor of our rich club.
17. Ambjørn, Jurkiewicz, Loll (2005), "Spectral dimension of the universe," *PRL* 95, 171301. — d_s method for Plan A.
18. Ambjørn, Durhuus, Jonsson (1997), *Quantum Geometry: A Statistical Field Theory Approach* (CUP).

**Network-geometry toolkit (methods for A/C):**
19. Krioukov et al. (2010), "Hyperbolic geometry of complex networks," *PRE* 82, 036106.
20. Colizza, Flammini, Serrano, Vespignani (2006), "Detecting rich-club ordering," *Nat. Phys.* 2, 110; Maslov & Sneppen (2002), *Science* 296, 910 (degree-preserving null).

**Suggested reading order:** 2 → 3 → 6 → 8 → 16 → 17 → 11 → 13, with 4 as reference.
