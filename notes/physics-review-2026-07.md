# Four lines of investigation: physics review

*2026-07-14. Written for a mathematician with a strong physics background (some GR/QFT, not expert).
Covers: (I) phase structure of the triangulation substance, (II) emergent roundness,
(III) vertex-degree structure and hubs, (IV) the |Ψ|² program.
Figures: `out/review_figs/fig1–fig4`, generated from data in this repository.*

---

## 0. Common setup and notation

All four investigations concern one statistical-mechanical object. A **configuration** is a
combinatorial triangulation T of the 3-sphere; we write its f-vector (f₀, f₁, f₂, f₃) for the
counts of vertices, edges, triangles, tetrahedra. Configuration space is the set of all such
triangulations, connected by **bistellar (Pachner) moves**, and we sample the Gibbs measure

  P(T) ∝ exp(−S(T)),  temperature fixed at 1,

by Metropolis–Hastings over Pachner moves. The action ("objective") actually implemented in
`source/sampler.d` is

  S(T) = c_f · (f₃ − N)² + k · (f₁ − 6f₃/d̄_t)² + γ · HDV(T) + β · VDV(T),

with the following meaning:

| symbol | meaning |
|---|---|
| N | target number of tetrahedra (the "volume"); c_f = 0.1 throughout |
| d_e | **edge (hinge) degree**: number of tetrahedra containing edge e; mean d̄ = 6f₃/f₁ |
| d̄_t, k | target mean edge degree and pin stiffness (k = `num_hinges_coef`; the pin is implemented on the equivalent hinge *count* f₁) |
| D_v | **vertex degree**: number of tetrahedra containing vertex v; mean D̄ = 4f₃/f₀ |
| VDV | vertex-degree variance Var_v(D_v), minus its quantization floor (below); coefficient β |
| HDV | edge-degree variance Var_e(d_e), minus its integer floor; coefficient γ |
| g | the **reduced VDV coupling** g = β/N (empirically the right variable — see §I) |

Geometry enters through Regge calculus. All edges have unit length, so every tetrahedron is a
regular Euclidean simplex with dihedral angle θ_tet = arccos(1/3) ≈ 70.53°, and curvature is
concentrated on edges with **deficit angle** δ_e = 2π − d_e·θ_tet. With unit lengths the total
scalar curvature is

  ∫R dV ≈ 2 Σ_e δ_e = 2 f₁ θ_tet (d̄* − d̄),  d̄* ≡ 2π/θ_tet ≈ **5.10430**,

so the *sign of the mean curvature is the sign of (d̄* − d̄)*: mean edge degree below flat ⇒
positively curved on average. There is a companion flat value for vertices: the solid angle of a
regular tetrahedron at a vertex is Ω_tet = 3 arccos(1/3) − π ≈ 0.5513 sr, giving flat vertex
degree D* = 4π/Ω_tet ≈ **22.795**.

Everything below refers to the **certified seed library**: families of 32 independent replicas
per parameter cell, passed through a multi-chain convergence gate (quantized split-R̂ over all
scalar observables *and* full degree distributions, ESS floors, τ-growth audits). "We measured X
in the ensemble" always means "over gate-certified equilibrium replicas."

---

## I. Phase structure of the triangulation substance

### I.1 Background and motivation

The dynamical-triangulations (DT) program interprets a sum over triangulations as a regularized
Euclidean path integral for gravity: Z = Σ_T e^{−S_Regge(T)}, with the continuum limit sought at
phase transitions of the statistical model [1–3]. The classic Euclidean results are sobering:
in 3D and 4D the unconstrained models have a **first-order** transition between a *crumpled*
phase (vanishing extent, singular structures) and a *branched-polymer* phase (Hausdorff
dimension 2) — with no smooth geometry on either side [3–5]. Much of the modern interest
(e.g. CDT, §II) comes from modifying the ensemble to evade this.

Our ensemble is a different slice through the same configuration space: we **fix** topology
(S³), edge lengths (all 1), volume (pinned f₃), and *mean* curvature (pinned d̄), and we pay an
energy for curvature *inhomogeneity* (β·VDV). In heat-bath language: a substance whose
microstates are PL 3-geometries, held at temperature 1, with a tunable "curvature-smoothing"
coupling. Before making any continuum claims one must answer the classical materials-science
questions: what is the equation of state? are there phase transitions along g? where does
ergodicity fail? That last question turned out to be governed not by critical phenomena but by
**glass physics** — the phenomenology of supercooled liquids: super-Arrhenius relaxation,
history dependence (fictive temperature), and kinetic arrest without a thermodynamic
singularity [6–8].

**Reading list.**
[1] J. Ambjørn, B. Durhuus, T. Jonsson, *Quantum Geometry: A Statistical Field Theory Approach*
(Cambridge UP, 1997) — the canonical DT text; Chs. 4–6 cover branched polymers and simplicial gravity.
[2] F. David, *Simplicial quantum gravity and random lattices*, Les Houches lectures 1992
(arXiv:hep-th/9303127) — the best "for theorists" introduction.
[3] J. Ambjørn, S. Varsted, *Three-dimensional simplicial quantum gravity*, Nucl. Phys. B 373 (1992) —
3D Euclidean DT phase structure.
[4] P. Białas, Z. Burda, A. Krzywicki, B. Petersson, *Focusing on the fixed point of 4d simplicial
gravity*, Nucl. Phys. B 472 (1996); B. de Bakker, Phys. Lett. B 389 (1996) — the first-order verdict in 4D.
[5] S. Catterall, J. Kogut, R. Renken, *Phase structure of four-dimensional simplicial quantum
gravity*, Phys. Lett. B 328 (1994).
[6] P. Debenedetti, F. Stillinger, *Supercooled liquids and the glass transition*, Nature 410 (2001) —
short authoritative review.
[7] A. Cavagna, *Supercooled liquids for pedestrians*, Phys. Rep. 476 (2009) — the friendliest
serious introduction to glass phenomenology.
[8] L. Berthier, G. Biroli, *Theoretical perspective on the glass transition*, Rev. Mod. Phys. 83 (2011).

### I.2 What we learned

**(a) The reduced coupling and the equation of state.** A single Pachner move changes any global
degree statistic by O(1/N), so the physically meaningful coupling is g = β/N, not β. In the
fluid window the ensemble obeys a clean equation of state with full N-collapse
(Fig. 1, left; N = 10³, 10⁴, 10⁵ on one curve):

  ⟨VDV⟩ ≈ 0.11 / g.

This is quantitative equipartition. Writing the penalty as E = β·VDV = (β/f₀) Σ_v δ_v² with
δ_v = D_v − D̄, each of the f₀ vertex-degree fluctuation modes is harmonic with stiffness β/f₀;
at temperature 1 each carries ⟨(β/f₀)δ_v²⟩ = ½, so

  ⟨VDV⟩ = f₀/(2β) ≈ (0.176 N)/(2β) ≈ 0.09/g,

within ~20% of the measured 0.11/g (the excess is consistent with anharmonicity and the
quantization floor). Correspondingly the energy density E/N = g·⟨VDV⟩ ≈ 0.10–0.12 is
*coupling-independent*, and the specific heat c = Var(E)/N ≈ 0.11–0.14 k_B per facet is flat in
both g and N. **There is no thermodynamic transition anywhere in the fluid window** — no peak,
no latent heat, no feature drifting with N. The medium behaves as a gas of ≈ N/8 harmonic
"degree phonons."

**(b) The kinetic glass at g\* ≈ 0.01–0.03.** Pushing g upward, the *dynamics* fails before the
*thermodynamics* does anything (Fig. 1, middle and right): Metropolis acceptance collapses
(0.33 → 0.01), the relaxation time τ is non-monotonic — first τ ∝ 1/g (an overdamped soft mode:
weaker spring, slower relaxation), then a super-Arrhenius blow-up at the same reduced g\* for
N = 10³ and 10⁵. Beyond g\*, chains freeze on distinct plateaus (four chains at N = 10³,
raw β = 200 froze at VDV ≈ 0.6, 9, 17, 21) and become history-dependent: *over-squeezed* initial
states trap at higher VDV than gently-cooled ones — a fictive-temperature effect, and the
concrete reason we abandoned annealing in favor of fixed-g sampling with a convergence gate.

The mechanism is a **quantization floor**. A vertex link is a triangulated S², which has an even
number of faces, so vertex degrees are *even integers*; the minimum achievable VDV at mean D̄ is

  VDV_floor = 4 y (1 − y),  y = frac(D̄/2),

(implemented in the sampler, subtracted from the penalty). Near-regular triangulations — those
close to saturating the floor — are combinatorially rare and connected by few Pachner moves, so
the landscape at large g is rugged: *a liquid that vitrifies before it can crystallize*.

**(c) The lower boundary is soft: entropic self-confinement.** Extrapolating ⟨VDV⟩ = 0.11/g and
τ ∝ 1/g to g → 0 predicts divergences. Both are wrong. Below the window the density of states
takes over: Ω(VDV) has its own restoring force with finite width set by N, so at β = 0 the VDV
saturates at a finite, N-dependent plateau (≈ 570 at N = 10³, ≈ 1500 at N = 10⁴) with finite τ.
Empirically, β = 0 ensembles **certify at standard production length** (VDV R̂ ≈ 1.02,
ESS ≈ 2000 at N = 10⁴). Equipartition holds only where the thermal well is narrower than the
entropic one.

**(d) A refinement from quench experiments.** The "glass wall" at the upper end resolves into
three distinct regimes: (i) at small N (≲ 600), genuine finite-size freezing (too few moves);
(ii) at intermediate N and g = 6–8×10⁻³, the equilibrium state *exists and is reachable* — a
quench from a nearby equilibrium relaxes ergodically (energy 412 → 86 in ~20 sweeps, betweeen-chain
spread collapsing to 4%) — but slow collective modes in the edge/facet sector (τ ~ 50–75 samples)
make over-dispersed multi-chain certification expensive; (iii) only g ≳ 10⁻² is true glass. The
distinction between "equilibratable" and "certifiable in standard runs" matters operationally
and was invisible to static tests alone.

**(e) A tradeoff, not a knob: the edge pin.** The pin stiffness k is part of the definition of
the ensemble. Stiffening k fixes the slow edge-degree mode (τ ∝ 1/k) but freezes the facet
count — the two are algebraically coupled through d̄ = 6f₃/f₁ — with the joint optimum near
k = 2. No choice of k opens the high-g frontier at small N: the boundary is a genuine joint
reachability limit of the (d̄, f₃) sector, not an artifact of one bad coordinate.

### I.3 Commentary

Two things here strike me as genuinely interesting physics rather than methodology.

First, **the glass is geometric**. Fragile-glass phenomenology — super-Arrhenius, fictive
temperature, plateau freezing — usually needs quenched disorder or frustrated interactions. Here
it emerges in a perfectly homogeneous combinatorial ensemble from a *parity constraint*: links
of vertices are 2-spheres, 2-spheres have even face counts, hence the even-degree lattice and
the 4y(1−y) floor. The obstruction to "crystallizing" (reaching degree-regular triangulations)
is the sparseness of Pachner paths into an exponentially thin corner of configuration space.
That is a purely topological-combinatorial origin for a canonical dynamical phenomenon, and it
sets a hard boundary — g < g\* — on which ensembles can honestly be called equilibrium at all.

Second, **equipartition works**, quantitatively, with the right mode count. It licenses thinking
of curvature inhomogeneity as a *phonon gas riding on the geometry*, with β¹ playing an inverse
temperature for that gas alone. The flat specific heat then says the substance has no internal
structure along g until dynamics kills you — the interesting axis for thermodynamic transitions
is d̄ (the curvature axis), which we have so far only pinned, not swept. The known 3D DT
first-order transition lives on that axis; finding (or excluding) it in *our* fixed-edge-length,
variance-controlled slice is an open question the library was built to answer.

---

## II. Emergent roundness: does the ensemble look like round S³?

### II.1 Background and motivation

The central question of any discrete-gravity program: does the statistical ensemble of discrete
geometries admit a scaling limit that is a *smooth space*, and if so which one? For Euclidean DT
the answer was famously "no" — crumpled or branched-polymer, never extended [1–3 of §I]. The
causal (CDT) modification produced the landmark counterexample: a 4D ensemble whose volume
profile matches Euclidean de Sitter space, i.e. a round S⁴ [9–11]. For our ensemble — Euclidean,
fixed S³ topology, curvature-variance-controlled — the natural analog target is the **round S³**.

The mathematically rigorous precedent cuts the other way and sharpens the question. Uniform
random 2D triangulations converge (Gromov–Hausdorff) to the **Brownian sphere**: topologically
S², but with Hausdorff dimension 4 and nowhere a manifold metrically [12, 13]. Entropy
overwhelmingly favors fractal geometry; *smoothness is atypical and must be paid for*. Our β·VDV
term is exactly such a payment, and the question "how much smoothing buys roundness?" is
quantitative.

Two more pieces of background frame the measurements. (i) **Curvature lives on edges** (Regge
[14]); with fixed unit edge lengths the deficit quantum is θ_tet ≈ 70.5°, so pointwise curvature
never refines — Cheeger–Müller–Schrader convergence [15] does not apply, and the correct target
is a *coarse-grained/GH* limit in which micro-roughness self-averages. (ii) On S³ the
Kazdan–Warner classification [16] puts every smooth function in the realizable-scalar-curvature
class, so sweeping d̄ across d̄* meets no obstruction — any breakdown is entropic/dynamical, not
topological. Finally, comparing metric *shape* requires an honest metric: graph distances are
quasi-isometric proxies; the heat method [17] and Steiner-type shortest-path bounds give PL
geodesics (only the latter with a guarantee).

**Reading list.**
[9] J. Ambjørn, J. Jurkiewicz, R. Loll, *Emergence of a 4D world from causal quantum gravity*,
PRL 93 (2004); [10] *Planckian birth of a quantum de Sitter universe*, PRL 100 (2008).
[11] R. Loll, *Quantum gravity from causal dynamical triangulations: a review*, CQG 37 (2020) —
start here; also J. Ambjørn, A. Görlich, J. Jurkiewicz, R. Loll, Phys. Rep. 519 (2012).
[12] J.-F. Le Gall, *Uniqueness and universality of the Brownian map*, Ann. Probab. 41 (2013);
G. Miermont, Acta Math. 210 (2013).
[13] J.-F. Le Gall, *Random geometry on the sphere*, ICM lecture (2014) — a survey written for
mathematicians; ideal entry point.
[14] T. Regge, *General relativity without coordinates*, Nuovo Cimento 19 (1961).
[15] J. Cheeger, W. Müller, R. Schrader, *On the curvature of piecewise flat spaces*,
Commun. Math. Phys. 92 (1984).
[16] J. Kazdan, F. Warner, *Scalar curvature and conformal deformation of Riemannian structure*,
J. Diff. Geom. 10 (1975).
[17] K. Crane, C. Weischedel, M. Wardetzky, *Geodesics in heat*, ACM Trans. Graph. 32 (2013).
[18] D. Burago, Y. Burago, S. Ivanov, *A Course in Metric Geometry* (AMS, 2001) — GH convergence.
[19] J. Ambjørn, J. Jurkiewicz, R. Loll, *Spectral dimension of the universe*, PRL 95 (2005) —
the complementary diffusion-based dimension probe (we haven't run it yet; it's on the list).

### II.2 What we learned

**Diagnostics.** On the unit round S³, the geodesic-distance distribution from any point has
shell density p(u) = (2/π) sin²(u) and CDF

  F(u) = (u − sin u cos u)/π,  u ∈ [0, π],  mean π/2.

We measure the empirical distance distribution (BFS on the 1-skeleton or dual graph; now also
true PL geodesics), rescale its mean to π/2 (comparing *shape* only), and report
**Δ_S3 = KS distance to F** (0 = perfectly round) and a **Hausdorff dimension d_H** from the
small-r volume-growth slope log N(≤r)/log r (round S³ ⇒ 3; branched polymer ⇒ 2; crumpled ⇒
effectively ∞, diameter ~ log N).

**(a) The low-VDV ensemble flows toward round S³ as N grows.** For the positive-curvature,
low-VDV families the finite-size trend is unambiguous and monotone. In the original ladder
(VDV ≈ 13): d_H = 1.41 → 1.79 → 2.56 and Δ_S3 = 0.147 → 0.074 → 0.022 for N = 10³ → 10⁴ → 10⁵.
In the current certified library (g = 2×10⁻³, to N = 5.6×10⁴): d_H → 1.75–1.79 and
Δ_S3 → ~0.10, all three distance metrics agreeing (Fig. 2). Both diagnostics are *still moving*
at our largest N — the N = 10⁵ tiers now queued are precisely what the limit needs.

**(b) High VDV is a different phase.** The same test on β = 0 ensembles (VDV ≈ 570–3100) gives
d_H = 2.7 → 4.2 → 5.9 (diverging through and past 3 — unphysical for a manifold, the
small-world signature) with Δ_S3 *rising* (0.13 → 0.21) and near-stagnant diameters. Low- and
high-VDV ensembles flow to **opposite fixed points** under N: extended-round vs
crumpled-small-world. This is our slice's version of the DT crumpled↔extended divide, with the
VDV coupling as the control parameter, and it proves curvature-variance control is *necessary*
for roundness.

**(c) Curvature sign is sub-leading; homogeneity is the driver.** At fixed low VDV, the three
edge-degree targets d̄_t = 5.0043 (positive), 5.1043 (≈flat), 5.2043 (negative) converge to
round nearly identically (Δ_S3 = 0.022/0.027/0.027 at N = 10⁵ in the older family). Curvature
shows up in the *scale* — effective radius, visible in the diameters (396/464/335) — but not in
the *shape*. Roundness, in this ensemble, is bought by homogeneity (low VDV), with the ±2%
mean-curvature differences acting as a radius knob. (A scale-aware fit that does *not* rescale
the mean — fitting R directly — is the natural follow-up and is on the queue-fed roadmap, along
with the new d̄_t = 4.8/5.4 families that triple the curvature lever arm.)

**(d) The verdict survives the honest metric (this week's result).** Graph distance is only a
quasi-isometric proxy for the PL metric, with *direction-dependent* distortion — a standing
worry for shape claims specifically. We implemented the true PL geodesic distance
(Steiner-point Dijkstra: a rigorous upper bound converging to the true flat-cone geodesic from
above; validated on the 600-cell against the analytic round metric, where it also agrees with
the independently-validated heat method). Findings: the true geodesic is a **near-constant
~0.89–0.91 contraction** of edge distance across N = 10³ → 5.6×10⁴ (mean ratio 0.887 → 0.907),
and the mean-rescaled distribution *shape is unchanged*. Rerunning the full roundness analysis
under the true metric moves every number slightly **toward** round (Δ_S3 lower by ~0.01–0.016;
d_H equal or higher; Fig. 2). So the edge-distance conclusions were not lattice artifacts — if
anything they *understated* roundness. (Cautionary tale logged on the way: the heat method,
which carries no error guarantee, silently violated even the crude edge upper bound on
disordered seeds while remaining essentially exact on the uniform 600-cell.)

### II.3 Commentary

The result I keep coming back to is **(c)**. A priori one might expect the mean curvature — the
thing that couples to the leading term of any effective action — to control the approach to
roundness. Instead the ensemble says: *first become homogeneous, then your mean curvature just
sets your radius.* That is exactly the structure of the smooth story (constant-curvature metrics
as the homogeneous representatives in a conformal class, curvature only fixing the scale), and
it emerged from a combinatorial measure with no metric input beyond "all edges equal." It also
gives the program a clean division of labor: β·VDV is a *shape* coupling, d̄_t is a *scale*
coupling.

Second, the Brownian-sphere contrast deserves emphasis for a mathematical audience: the 2D
theorem [12] shows the unbiased measure produces a fractal, so any claim of smooth emergent
geometry *must* be a statement about a biased measure — and (a)/(b) exhibit precisely the
bias-controlled transition between the fractal-like and smooth-like basins, in 3D, with a
one-parameter dial. Making claim (a) rigorous even at the level of "d_H → 3 in probability"
would be a genuinely new kind of theorem, and the finite-size data now behave as if such a
statement were true.

Third — a small methodological point with physical content — the ~0.9 edge-to-geodesic ratio
being *constant in N* is itself evidence of metric self-averaging: the lattice-level anisotropy
of the 1-skeleton does not accumulate at large scales. That is the concrete, measurable face of
the "coarse-graining kills the curvature quantum" hope in §II.1(i).

---

## III. Vertex-degree structure: hubs, condensates, and two kinds of crumpling

### III.1 Background and motivation

Once an ensemble is known to be small-world/crumpled, the mechanism matters. Complex-network
theory offers two archetypes: **hub-driven** shortcuts — scale-free degree distributions
P(D) ~ D^(−α) whose extreme members dominate connectivity (Barabási–Albert [20]; attack
tolerance [21]) — versus **diffuse** densification, where the whole graph thickens (small
worlds without hubs, Watts–Strogatz [22]). Dynamical triangulations has its own version of the
first: the crumpled phase of 4D DT is dominated by *singular vertices* whose stars contain a
finite fraction of the total volume [23]. Whether our crumpled phases are singular-vertex-like
or diffuse decides whether "surgery" (removing the hubs) can recover an extended geometry — and
connects to speculative "disordered locality" ideas in which an emergent geometry is spoiled by
a sparse network of nonlocal links. The probes are network-theoretic: degree distributions and
their N-scaling, percolation of degree-thresholded subgraphs (targeted attack [21, 24]), and
rich-club structure [25].

**Reading list.**
[20] A.-L. Barabási, R. Albert, *Emergence of scaling in random networks*, Science 286 (1999).
[21] R. Albert, H. Jeong, A.-L. Barabási, *Error and attack tolerance of complex networks*,
Nature 406 (2000).
[22] D. Watts, S. Strogatz, *Collective dynamics of "small-world" networks*, Nature 393 (1998).
[23] S. Catterall, G. Thorleifsson, J. Kogut, R. Renken, *Singular vertices and the triangulation
space of the D-sphere*, Nucl. Phys. B 468 (1996) — singular structures in the DT crumpled phase.
[24] D. Callaway, M. Newman, S. Strogatz, D. Watts, *Network robustness and fragility:
percolation on random graphs*, PRL 85 (2000).
[25] V. Colizza, A. Flammini, M. Serrano, A. Vespignani, *Detecting rich-club ordering in
complex networks*, Nature Physics 2 (2006).

### III.2 What we learned

Both studies concern β = 0 (VDV uncontrolled) ensembles, differing only in the edge pin. A
useful exact fact first: vertex degrees are always **even** (the link of a vertex is a
triangulated S²; 3F = 2E forces even F), and for a flat-pinned S³ triangulation the f-vector
identity (χ = 0, f₂ = 2f₃, d̄ = 6f₃/f₁ = d̄*) forces

  D̄ = 4f₃/f₀ = 4d̄*/(6 − d̄*) ≈ 22.795 = D*,

i.e. **pinning the mean edge degree at flat pins the mean vertex degree at flat, exactly,
independent of N.** All the freedom is in the shape of the distribution.

**(a) k = 0 (no edge pin): crumpling without hubs — hypothesis refuted.** The natural guess —
that the log-like distances are caused by a few singular hubs — is wrong here. The maximum
degree grows *at the same rate as the mean* (maxD ∝ N^0.41 vs D̄ ∝ N^0.38; ratio ~2, constant),
no vertex ever exceeds 10× the mean, and the whole degree PMF **collapses under D/⟨D⟩**
(Fig. 3, left): a self-similar distribution that simply inflates. The mechanism is
**densification**: f₀ ∝ N^0.61 (sublinear), so D̄ ∝ N^0.38 diverges and the 1-skeleton becomes
asymptotically dense — at N = 5.6×10⁴ a typical vertex is adjacent to ~⅕ of all vertices, so
diameter ~3 is trivial. Corroboration: the *dual* graph is exactly 4-regular — hubless by
construction — yet also small-world; the crumpling lives in expander-like gluing, not in special
vertices. The distribution is bimodal at large N: ~12–14% degree-4 "cone tips" (maximal positive
curvature) plus a broad negative-curvature "sea" around the mean; the tips form an (almost
exact) independent set — dust hanging off the connected sea. Nothing localizes at flat.

**(b) k = 2 flat pin, β = 0: the hub condensate.** With D̄ nailed to D* (measured exponent
N^0.00), the entropy goes instead into the *tail*: CV = std/mean grows ∝ N^0.19 without bound,
maxD ∝ N^0.50, and the distribution develops a **scale-free tail with exponent α ≈ 3**
(P(>D) ~ D^(−2); Fig. 3, right — the extreme-value scaling maxD ~ N^(1/(α−1)) = √N checks
independently). These hubs are real, and organized: the induced subgraph on high-degree vertices
is a single connected **rich-club core**, each hub wired to ~16–18 other hubs. Geometrically the
core is strange: by *count* it is codimension-1-like (|core| ∝ N^(2/3), 100% of its vertices on
its surface — no interior), but by *contact* it is space-filling (∂(bulk) ∝ N^1: ~94% of all
bulk vertices sit one edge from a hub). The right image is **vascular** — a capillary network
permeating the bulk, not an organ with a skin.

**(c) The surgery mirage.** Removing hubs (degree > c·D*) from the k = 2 ensemble at N ≤ 10⁴
produced an apparently extended remnant — diameter ∝ size^(1/3), i.e. genuinely 3-dimensional
looking. Extending the tiers to N = 5.6×10⁴ *killed* this: the connectivity-safe cut must climb
(θ ∝ N^0.4), the removable fraction shrinks, and the de-hubbed diameter stalls (exponent
α ≈ 0.09 — still small-world). The hub structure is a *hierarchy*, not a separable sliver: cut
the top tier and the next tier of near-hubs re-shortcuts the bulk. Fixed-cut versions of the
claim turn out to be percolation artifacts (the "giant" label migrating to stringy fragments).
Every quantitative trace of this hierarchy carries the same exponent: CV ∝ N^0.19, the
"both-pieces-connected" cut window shrinks as N^(−0.37), its θ-edges climb as N^0.35–0.41, and
family VDV ∝ N^0.37. One scale-free hierarchy, one exponent.

**(d) β > 0 has no hub core — and ordinary percolation instead.** In the extended
(VDV-controlled) phase, degree-thresholding behaves like textbook site percolation on a
homogeneous medium: there is **no** cut with both pieces connected; instead a double-percolation
gap opens with N, whose edges are set by targeted-attack asymmetry (removing high-degree
vertices fragments the remnant efficiently ⇒ θ_L ≈ 1.37 D*; the high-degree piece alone
percolates only while it is the majority ⇒ θ_H → D*). The gap narrows as β increases (the
degree distribution narrows, CV: 0.32 → 0.21 from g = 2→4×10⁻³), but the glass wall (§I) caps
how narrow the distribution can be made — the gap has a glass-limited minimum width. The
β = 0 → β > 0 transition is thus visible purely combinatorially: *delocalized permeating core →
homogeneous medium*.

### III.3 Commentary

The sharpest lesson is the **k = 0 / k = 2 dichotomy at identical β = 0**: the same entropic
pressure, given different conserved quantities, condenses in completely different ways. Unpinned,
it inflates the *mean* (densification — a condensation in the f-vector itself); pinned at flat,
it cannot touch the mean (the f-vector identity above is a hard constraint) and instead builds a
*scale-free tail* — a connected condensate of curvature-variance. This is a nice instance of a
general statistical-mechanics moral: constrain a mean, and fluctuations reorganize into extreme
structure. It is also, to my knowledge, a different beast from the DT singular-vertex phase
[23]: our core is an extensive *hierarchy* of hubs (α ≈ 3, max ~ √N) rather than O(1)
macroscopic singular vertices — the marginal case, since α = 3 is exactly the boundary where
the variance diverges but the mean does not.

Second, the surgery mirage is a cautionary tale I'd want any reader of finite-size claims to
internalize: at N = 10⁴ the de-hubbed graph sat *exactly* on the size^(1/3) line — as clean a
"hidden 3D skeleton" as one could ask for — and it was a finite-size window. The diagnosis was
only possible because the certified library let us extend the *same* family 5.6× further. The
physical conclusion is worth stating positively: **locality cannot be excavated from a β = 0
sample; it must be selected by the measure** (β > 0). The intact β > 0 ensembles are more
extended than any surgically-cleaned β = 0 remnant.

Third, the recurring 0.37–0.38 exponent (tail growth, window shrinkage, VDV growth — and, in
§IV, the spectral gap) is exactly the kind of coincidence that usually means one underlying
scaling field. We have not identified it analytically; α ≈ 3 plus extreme-value scaling gives
√N for the max, and the ~N^0.19 CV points to the variance-marginality of α = 3, but the 0.37 is
so far an empirical fixed point of several independent measurements. I'd flag it as the most
concrete open puzzle in the combinatorial sector.

---

## IV. The |Ψ|² program: the ensemble as a quantum ground state

### IV.1 Background and motivation

In canonical GR the configuration space is Wheeler's **superspace** — 3-geometries modulo
diffeomorphisms — and a quantum state of the universe is a functional Ψ[h] on it, obeying the
Wheeler–DeWitt constraint [26, 27]. Proposals like Hartle–Hawking [28] specify Ψ; |Ψ[h]|² is
then a probability measure on 3-geometries. Our sampler provides the discrete mirror image:
the certified ensembles *are* probability measures on a discrete superspace (triangulations of
S³ modulo isomorphism), so one can ask whether P(T) = e^{−S(T)} behaves like |Ψ₀[T]|² for an
interesting quantum state.

The correspondence is exact, not an analogy. For a reversible Markov chain with transition
matrix W and stationary distribution P, the operator

  H = 𝟙 − D^(1/2) W D^(−1/2),  D = diag(P),

is symmetric and positive semidefinite, with ground state ψ₀(T) = √P(T) = e^(−S(T)/2) at energy
0; the Metropolis dynamics is exactly *imaginary-time Schrödinger evolution* generated by H, and
every measured autocorrelation time is an inverse energy gap: τ_n = 1/E_n. This is the
**Rokhsar–Kivelson (RK) construction** [29]: classical Gibbs measures realized as ground-state
wavefunctions squared, with a stochastic parent Hamiltonian. The RK literature [30–32] supplies
the phenomenology to test against — e.g. dynamical exponent z ≈ 2 at RK criticality (gap closing
as L^(−2)) and equal-time correlations that are exactly classical-critical. Two important
honesty clauses: the *ground state* depends only on S (move-set independent), but the *excited
spectrum* — everything dynamical — depends on the chosen move set; and the construction is
Euclidean/timeless (relational time would need matter clocks). The flagship physics target is
the **curvature structure factor**: in a constrained (gravity-like) theory the long-wavelength
scalar mode is suppressed, S(k→0) → 0 — in condensed-matter language, *hyperuniformity* of the
curvature distribution [33]; a generic thermal field instead has S(k→0) finite. Spectroscopy of
H is done with the standard lattice-QCD tool: the generalized eigenvalue problem (GEVP) on a
matrix of correlators [34].

**Reading list.**
[26] B. DeWitt, *Quantum theory of gravity I: the canonical theory*, Phys. Rev. 160 (1967).
[27] J. Wheeler, *Superspace and the nature of quantum geometrodynamics*, in Battelle Rencontres (1968).
[28] J. Hartle, S. Hawking, *Wave function of the universe*, Phys. Rev. D 28 (1983).
[29] D. Rokhsar, S. Kivelson, *Superconductivity and the quantum hard-core dimer gas*, PRL 61 (1988).
[30] C. Henley, *From classical to quantum dynamics at Rokhsar–Kivelson points*, J. Phys.:
Condens. Matter 16 (2004) — the clearest exposition of the chain↔Hamiltonian dictionary.
[31] E. Ardonne, P. Fendley, E. Fradkin, *Topological order and conformal quantum critical
points*, Ann. Phys. 310 (2004).
[32] C. Castelnovo, C. Chamon, C. Mudry, P. Pujol, *From quantum mechanics to classical
statistical physics...*, Ann. Phys. 318 (2005) — the general stochastic-matrix-form construction.
[33] S. Torquato, *Hyperuniform states of matter*, Phys. Rep. 745 (2018) — S(k→0)→0 as a
physical organizing principle.
[34] M. Lüscher, U. Wolff, Nucl. Phys. B 339 (1990) — the GEVP method.

### IV.2 What we learned

**Method.** The production runs already record, per sample, full vertex- and edge-degree
histograms — a basis of 100–300 observables {O_i} at zero additional simulation cost, over
~6300 chains. We form connected correlators C_ij(t) = ⟨δO_i(s) δO_j(s+t)⟩ and solve the GEVP
C(t)v = λ(t,t₀)C(t₀)v, extracting E_n = −ln λ_n/(t−t₀) — variational upper bounds on the low
spectrum of H. (Numerical care: whiten with the top PCA modes of C(0), then restrict to the
positive eigenspace of the whitened C(t₀); jackknife over the 32 chains.)

**The three-phase gap taxonomy (Fig. 4).** Tracking the first gap E₁ along N-ladders:

| phase | E₁(N) | interpretation |
|---|---|---|
| β > 0 extended (all k, g = 2–4×10⁻³) | ∝ N⁰ (τ₁ ≈ 30–60 sweeps) | **gapped**; gap set by the coupling, not the size |
| β = 0, k = 2 flat pin (scale-free hub phase) | ∝ N^(−0.38) | gap closes with the **hub-tail exponent** |
| β = 0, k = 0 (crumpled/densifying) | ∝ N^(−0.19) | slow closure, expander/mean-field-like |

The gapped result **refutes the naive RK expectation z ≈ 2 for the extended phase**: with
L ~ N^(1/3), z = 2 would require E₁ to fall by ×3.2 across the k = 4 ladder — excluded (Fig. 4,
green). In RK language our extended states sit in a gapped (non-critical) region of their
phase diagram, consistent with the τ ∝ 1/g fluid EOS of §I. The physical mode content is
meaningful too: the slowest eigenvector is consistently the *curvature-variance relaxation*
channel (VDV plus the low-edge-degree histogram bins), and the second mode is an edge-histogram
shape mode — with the weight shifting between vertex- and edge-sector bins as g changes.

The scale-free β = 0 phase behaves oppositely: its gap closes as N^(−0.38) — numerically the
same exponent as every hub-hierarchy observable in §III. Gap closure tracks scale-free-ness,
not size: geometry with a growing condensate has soft collective modes; homogeneous geometry
does not.

**Consequence for the flagship test.** A gapped RK state generically has S(k→0) finite, whereas
gravitational constraint-suppression predicts S(k→0) → 0. These now *disagree* for our extended
phase, so the curvature structure factor — first pass computable from stored snapshots —
falsifies something either way: either the "ensemble as gravitational wavefunction" reading
loses its flagship signature, or we find hyperuniformity coexisting with a spectral gap, which
would be a strong and specific structural claim about the state. The quench machinery
(relaxation from prepared initial ensembles, two-tier recording of observables and full
snapshots) is validated and ready for the dynamical side (z from relaxation, aging/FDT).

### IV.3 Commentary

What I find most satisfying here is that the RK identity converts *bookkeeping we already had*
into spectroscopy. Every convergence diagnostic the certification gate computes — every τ — was
secretly an energy gap; the GEVP just organizes ~300 observables into variational eigenstates.
The taxonomy in Fig. 4 is, to my eyes, the single most physics-dense figure the project has
produced: it cleanly separates the three phases *dynamically*, assigns the soft mode a name
(curvature-variance relaxation), and quantitatively ties the spectral gap of the scale-free
phase to a combinatorial exponent measured by entirely static means (§III). That last link —
gap ∝ N^(−0.38) where CV, window width, and VDV all scale with the same 0.37–0.38 — is either a
coincidence to be dissolved or a small piece of genuine emergent scaling theory.

The gapped extended phase deserves a comment. One might have hoped the roundness-approaching
ensembles would be *critical* (they are, after all, supposed to approach a continuum). Instead
they are gapped at fixed g, which in RK phenomenology means: fine-tuning toward a critical
point is required for a relativistic-like continuum of soft modes, and our g-window may simply
not contain one. The honest possibilities are (i) the continuum limit here is of the
"coarse-grained smooth background" type without gapless excitations (perfectly adequate for the
§II roundness claims), or (ii) criticality lives elsewhere in the (g, d̄, k) space — e.g. on
the d̄ axis, where the DT-style transition should sit. The structure-factor measurement will
discriminate faster than more spectroscopy will.

Finally, the epistemic framing: this program is deliberately speculative, but it is *cheap* and
*falsifiable* — its first two measurements (GEVP, S(k)) recycle existing data, and its flagship
prediction now stands in explicit tension with its own spectral result. That is exactly the
position a speculative program should be in.

---

## Cross-cutting summary

One sentence per thread. The substance is a **harmonic fluid of curvature phonons with a
quantization-floor glass** bounding the accessible couplings (§I). Inside the accessible window,
**homogeneity buys roundness and mean curvature buys radius**, robustly under the true PL
metric (§II). Without homogeneity control, entropy condenses either into **densification (no
pin) or a scale-free vascular hub condensate (flat pin)** — and no surgery recovers locality;
the measure must do the work (§III). Viewed as a quantum state, the extended phase is **gapped**
with curvature-variance as its softest excitation, the scale-free phase closes its gap with the
hub exponent, and the next measurement (curvature hyperuniformity) is set up to falsify the
gravitational reading or the naive RK reading (§IV).

*Figures: `out/review_figs/fig1_phase_eos.png` (EOS collapse, acceptance, τ),
`fig2_roundness.png` (Δ_S3 and d_H vs N under three metrics; CDF overlay),
`fig3_hub_structure.png` (k=0 self-similar collapse vs k=2 scale-free tail),
`fig4_gevp_gaps.png` (spectral-gap taxonomy).*
