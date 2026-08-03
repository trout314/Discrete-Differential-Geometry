---
name: crystal-symmetry-group
description: "Exact Aut(K) for crystal triangulations is cheap (develop one frame); symmetry.py delivers orbits/stabilizers/BC-chain classes — and Aut(K) may legitimately EXCEED the physical crystal's space group"
metadata: 
  node_type: memory
  type: project
  originSessionId: c395a1df-09fb-4d44-bb80-ddc781dd7a5b
  modified: 2026-08-02T18:04:31.602Z
---

`python/discrete_differential_geometry/symmetry.py` (added 2026-08-02) computes
the **exact** combinatorial automorphism group of a 3-manifold triangulation,
replacing the brute-force-then-bucket pattern that was everywhere.

**THE CAVEAT TO REMEMBER.** A *triangulation* built from a real TCP crystal can
legitimately have a **larger** symmetry group than the physical crystal. Aut(K)
is guaranteed only to *contain* the space group's image on the supercell torus.
Three independent reasons: (1) K records adjacency, not distance — a relabelling
preserving every adjacency while distorting lengths is an automorphism but not
an isometry; (2) we triangulate the m×m×m torus, so lattice translations become
a *finite* group and wrap-around can manufacture coincidences the infinite
crystal lacks; (3) `tcp_reference.py` perturbs sites ~1e-6 before tiling, which
breaks the point group *metrically* while the connectivity may retain it. For
sampler physics Aut(K) is nevertheless the right group (the sampler sees K, so
Aut(K) is exactly what leaves the action/rates/species invariant) — but never
report |Aut| as a crystallographic order without checking the factorisation.
Empirically, at m=2,3,4 every TCP crystal tested has **no accidental
automorphisms** (|Aut| = m³ × centering × point group exactly) — a measurement,
not a theorem.

**Why it's cheap.** An automorphism of a connected 3-manifold triangulation is
determined by the image of ONE ordered tet (a "frame") — crossing a face forces
the new apex. So Aut = {developments of a base frame that close globally}, which
reuses the same propagation as `crystal_grains.develop`. Aut acts *freely* on
frames, so 24·nT / |Aut| is an integer (a free self-check) and every frame in
the known orbit needs no development at all. R m4 (57984 tets) takes 14 s;
results cache to a `.sym.npz` sidecar keyed by a facet-set hash.

**Measured (all factorise exactly, no accidentals):**

    A15 m3   |Aut| = 1296 = 27×48        orbits V/E/F/T = 2/3/4/3     3 chain classes
    C15 m3   |Aut| = 5184 = 27×48×4      2/3/4/3                      2 chain classes
    C15 m4   |Aut| = 12288 = 64×192      2/3/4/3                      2
    SIGMA m3 |Aut| = 432 = 27×16         5/20/29/17                   5
    C14 m3   |Aut| = 648 = 27×24         3/8/11/7                     4
    R m2/m3/m4 |Aut| = 144/486/1152      11/61/102/53 at ALL m       14 chain classes

Orbit counts are supercell-INDEPENDENT (properties of the crystal); chain
*lengths* scale with the box, the number of chain *classes* does not.

**Aut vs Aut+ — WHICH GROUP TO USE.** `sym.orientation_preserving` is itself a
`CrystalSymmetry`, so every orbit/stabilizer/chain method works on it unchanged;
`sym.is_chiral(kind, obj)` says whether an Aut-orbit splits into an enantiomer
pair (`None` if the group has no orientation-reversing element at all).

- **Aut** for anything the sampler sees. The sampler's state IS the abstract
  complex and its action/moves are combinatorial, so the Markov chain commutes
  with Aut and every sampler observable is constant on Aut orbits. Aut-orbit
  counts = the number of distinct cases to measure.
- **Aut+** for quantities that exist only after an orientation is chosen:
  Wilson-line spinor signs, helix handedness, anything with an ε tensor, the
  QRW Dirac coin. Aut-invariant = scalar; flips under Aut∖Aut+ = pseudoscalar,
  and averaging a pseudoscalar over an Aut orbit gives exactly 0.
- Rule of thumb: if computing it required fixing a sign, an orientation, or a
  seed parity, it's Aut+. `TransportContext`'s seed_tet gauge is exactly such a
  choice — its "parity flip mirrors every axis" warning IS this distinction.
- Would change if the action ever gains a parity-odd term (CS-like, or a
  chirality coupling off the six-web gauge field): then Aut+ governs sampling too.
- MEASURED: |Aut+| = |Aut|/2 for A15/C15/C14/SIGMA/R, and NO BC chain class is
  ever achiral (a tetrahelix must be handed), so every chain orbit splits —
  R: 14 under Aut, **28 under Aut+**. Crystals achiral, helices handed = racemic,
  so chirality-odd bulk averages MUST vanish; nonzero = bug or gauge leak.
- delta phase (P2₁2₁2₁, Sohncke/chiral group) has ZERO orientation-reversing
  automorphisms — recovered from combinatorics alone. Good validation test.

**NO-SILENT-DEGRADATION POLICY** (user directive, 2026-08-02: avoid silently
using incorrect data even at the cost of breaking existing data/scripts):

- WL is NEVER a fallback classification, only a labelled lower-bound
  diagnostic. A number that is sometimes the class count and sometimes an
  undercount is worse than no number — it can't be compared across inputs and
  corrupts anything aggregating it. On |Aut| = 1 the script states the exact
  (trivial) answer and says the count carries no information, rather than
  substituting WL, which measures a different quantity anyway.
- Disconnected input RAISES. Development from one frame can't see
  component-permuting automorphisms, so it would report a subgroup (usually
  trivial) — a silently wrong answer.
- The `.sym.npz` cache is VERIFIED, never trusted: digest match, every
  generator re-checked to permute the facet set, order re-derived from the
  closure (not read), free-action re-asserted. Any failure warns and recomputes.
- Non-manifold input (a face with ≠ 2 cofacets) raises in `TriView.nbr`.

**What this fixed / exposed:**

- `scripts/move_site_census.py` now classifies by exact face orbits. Its old WL
  bucketing reported **97** classes on R where there are **102** — WL fused 5
  inequivalent site classes. WL is kept as cross-check and as the fallback when
  |Aut| = 1 (melts), where orbits are correct but vacuous.
- `crystal_grains.develop` is now a thin adapter over `symmetry.develop_partial`;
  output verified byte-identical. Grain growth and Aut are the same traversal
  under two mismatch policies (boundary vs failure).
- **STILL OPEN — the BC-chain provenance hole.** ~10 scripts in
  `defect_dynamics/` do `bc_orbit(m, F[0])`, i.e. one arbitrary chain from tet 0.
  R has **14 inequivalent chain classes** with lengths spanning 99–2439 at m=3.
  Every pristine-crystal dock census / collider / S-matrix / FP-transport number
  was measured on whichever class tet 0 landed in, class unrecorded. Use
  `CrystalSymmetry.chain_orbits()` / `.chains` to pick and *record* a class.
- Also unconverted: `worm_catalog.canon_sig` (may over-merge — the "registry
  pass" its docstring promises is this), `knot_smatrix_sweep.relkey` (float
  fingerprint dedup of "A-window translation classes"), and the `% ns` registry
  in `crystal_grains.interior_vertices` / `defect_census` (works only while the
  tcp_reference labelling `v = cell·ns + site` survives relabelling).

- **`canonical_key`'s `exact` flag is discarded by some callers** —
  `crystal_flicker.py:78,80,82,157,158` and `flicker_fraction.py:113` do
  `k, _ = ds.canonical_key(...)`, using a possibly-truncated search result as a
  species certificate. Same class of bug (silently trusting a bound). NOT yet
  fixed: deciding what to do on `exact=False` changes the flicker JSON schema.

Terminology trap: "orbit" in `defect_dynamics/` means a **BC chain**, never a
group orbit. `symmetry.py` says "chain" for the walk, "orbit" only for groups.

Related: [[crystal-grains-tool]], [[flicker-background]],
[[bc-washboard-not-free-spirals]], [[tcp-provenance-wiring]],
[[fpkmc-m1-status]], [[cocycle-vertex-lift]]
