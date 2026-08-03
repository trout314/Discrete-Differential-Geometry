---
name: qrw-dirac-connection
description: "The user's SEPARATE project github.com/trout314/quantum-random-walks builds a Dirac quantum walk on the BC helix; it shares DDG's exact primitives (same twist Quat(1,0,1,2), same reflect(), same Swierczkowski uniqueness). The 3D coin problem = emergent-gravity coupling. How the two projects join."
metadata: 
  node_type: memory
  type: reference
  originSessionId: d562feea-0887-419b-9cd9-9899d180c7c2
  modified: 2026-07-27T15:51:18.704Z
---

2026-07-27. The user pointed to a second repo,
**github.com/trout314/quantum-random-walks** (cloned to scratch/qrw),
which constructs a discrete-time quantum walk on tetrahedral
triangulations that reduces to the **Dirac equation** at long distance
along a single BC chain. It is the "coherent spinor transport" idea I
raised — already built and working in 1D. Key facts and the precise
bridges to [[intrinsic-geometry]].

**QRW construction (theirs):**
- Coin space = 4D Dirac spinor. Operators tau_a = (sqrt7/4) beta +
  (3/4)(e_a . alpha), e_a the 4 regular-tetrahedron vertex directions.
  Unique family (in span{beta,alpha_i}) with Dirac correspondence
  Sum_a e_a^i tau_a = alpha_i, involutory, Hermitian, eigenvalues
  {-1,-1,1,1}. The (sqrt7/4, 3/4) split is forced by tetrahedral
  isotropy EE^T = (4/3)I.
- Walk W = S.C2.C1: shift S along the helix with per-site P+/P-
  projectors and unitary frame transport; coins C = cos(theta)I -
  i sin(theta)(f.alpha). Pure frame transport alone = "conveyor belt"
  (erases phase, no dispersion); a site-local coin breaks it and yields
  **massive 1D Dirac**: E^2 = m^2 + k^2, mass m ~ 0.59 sin(theta)
  (OPEN: why 0.59), exactly unitary. 3D walk (W = S_R.C.S_L.C) gives
  isotropic spreading but the coin/dispersion is UNTESTED in 3D (their
  frontier); lattice too small (radius ~22, saturates in ~2 steps).

**EXACT SHARED PRIMITIVES with DDG (verified, not analogy):**
- Twist per step = arccos(-2/3) = 131.81 deg = EXACTLY DDG chain-step
  Quat(1,0,1,2) (norm 6, cos = 2/6-1 = -2/3).
- Their reflection rule R_a(v) = v - 2(v.e_a)e_a IS development.py
  `reflect()`. Same seed tetrahedron.
- Their "Swierczkowski: unique walker paths" (listed OPEN re: two-slit)
  = DDG path-injectivity theorem, which DDG already PROVED in strong
  form (Wilson rotation alone AND endpoint alone each determine the
  path; ~18k paths, zero collisions). So DDG ANSWERS their two-slit
  question: distinct paths always carry distinguishable which-path data
  (arithmetic height), so genuine 2-path interference needs equal
  endpoint AND transport => same path.
- Their cut-and-project / 600-cell "decurving / geometric frustration"
  = DDG "curvature foam" (measured today: NO flat cycles, deficit
  spectrum 241/243 etc., e* = 2pi/arccos(1/3) flat-on-average).
- Their "12 helices, 6 perp L/R pairs through a point" = DDG dock
  geometry (12 BC helices per tet; local classes cos = k/5).
- Measured: their per-step frame transport rotates the 4D Dirac coin by
  +-25.659 deg; exact centroid step = sqrt(6)/6, axial advance =
  1/sqrt(10). DDG exact tools can pin these (and 0.59, c) symbolically.

**THE THREE JOIN POINTS (told the user, ranked):**
1. DDG already answers their two-slit open question (above).
2. **3D coin problem = emergent-gravity coupling.** Their 1D success is
   on a single chain = a FLAT direction. A 3D tetrahedral triangulation
   is a curvature foam, so a 3D walk cannot converge to FLAT Dirac — it
   converges to Dirac in a curved/torsionful background, and the DDG
   holonomy spectrum IS the connection it propagates in. Whether their
   walk_gen.c builds an idealized-flat lattice or the frustrated crystal
   is the fork between flat-3D-Dirac and Dirac-coupled-to-geometry.
3. **Do the 3D walk on the 600-cell first** (their data/600cell.mfd):
   S^3 tiled by 600 regular tets, 20 BC rings period 30, NO frustration
   => finite discrete Dirac spectrum, diagonalizable exactly, no
   boundary saturation; DDG development can verify ring holonomies close.

**Highest-value merge:** rebuild their frame transport on DDG's EXACT
integer-quaternion spinor lift (Quat/TransportContext) instead of numpy
polar decomposition => symbolic 1D walk, 0.59 and c fall out as exact
radicals. This is the concrete bridge between the emergent-gravity
(DDG) and matter-field (QRW) halves of the same tetrahedral geometry.
