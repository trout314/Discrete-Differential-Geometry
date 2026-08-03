---
name: bc-washboard-not-free-spirals
description: The 2->3 excitation has an exact 4-rung integer cost ladder Q; free slide <=> same rung; slide dS = c*dQ exactly; per-rung free network mapped
metadata: 
  node_type: memory
  type: project
  originSessionId: d562feea-0887-419b-9cd9-9899d180c7c2
  modified: 2026-07-27T22:37:00.648Z
---

**THE 2->3 PERSPECTIVE (R m2, 2026-07-27). Supersedes the earlier same-session
"washboard / ~1/8 clean" numbers, which were computed with a CONTAMINATED
action -- see the gotcha below.**

FUNDAMENTAL EXCITATION = a single 2->3 move on a triangle. Every triangle t has
an EXACT integer creation cost (change in sum_e deg_e^2):
  Q(t) = 2*(Scross - Sface) + 18,
  Sface = sum of the 3 face-edge degrees, Scross = sum of the 6 apex-shell
  edge degrees (the two apexes x,y of t's two tets, each to a,b,c).
On R m2 Q takes ONLY 4 values -- a 4-rung ladder:
  Q=46 (3456 tris, ALL [5,5,6]=pure (3,4,4), all-deg-5 apex shell -> lightest),
  Q=48 (3408), Q=50 (6480), Q=52 (1152, no (3,4,4)). Total 14496 = 2*N3.
  dS_create = c*Q + const, c = lam*e*/6 = 0.34034, rung gap 2c = 0.681.

SLIDE ENERGY IS EXACTLY THE RUNG DIFFERENCE (the 4-move mechanics are a RED
HERRING for energetics): a slide preserves N3 hence #edges and sum deg, so
  dS_slide(t1->t2) = c*(Q(t2) - Q(t1))   EXACT to 9e-16.
Verified: slide the (3,4,4) forward, measured dS vs c*dQ, max|.|=8.9e-16 with
pure edge action. FREE SLIDE <=> dQ=0 <=> same rung. ~37-45% of step-4 slides
are free (the earlier "~10-14% clean" was the contaminated action).

**ACTION-CONTAMINATION GOTCHA (this cost several wrong runs + a wrong "overlap"
story):** SamplerParams zeroes ONLY what you pass. Its constructor carries LIVE
DEFAULTS num_hinges_coef=0.05 AND codim3_degree_variance_coef=0.1 (vertex-degree
variance). The codim3 one perturbs on every slide (~1e-3, intensive) and it was
the mysterious residual on dQ=0 slides -- NOT the FK/n6 term. For FK/TCP
edge-only action you MUST pass num_hinges_coef=0, hinge_degree_variance_coef=0,
codim3_degree_variance_coef=0 explicitly (see [[vdv-hdv-conflicts-with-tcp]]).
With pure edge action the reframing is exact and dQ=0 slides are EXACTLY free.

FREE-SLIDE GRAPH per rung: nodes = triangles (defect sites); edge = a step-4 hop
along a BC spiral between two SAME-rung triangles (free, dS=0). Connected
component = free-reachable set (NOT just the flat same-Q set). "free-valence" of
a node = # free slides available there (= the ~4-5 clean slots/site: 1 fwd + 1
back + 2-3 off-chain). MEASURED (free_coverage.py, 432 window-cycles):
  Q=46: ONE connected branching web, covers ALL 1272 verts, <deg>~4.2,
        144/3456 sites ISOLATED (trapped defects).
  Q=50: ONE web, full coverage, <deg>~6.
  Q=48: FRAGMENTS -> 1 big (2352) + 288 tiny size-4 clusters, covers 1200.
  Q=52: 48 disconnected size-24 pockets (24 = small-helix period?), covers 912.
So the "free spirals" are really ONE connected same-rung WEB per energy level,
NOT 1D rails; the ground (3,4,4) roams the whole crystal (energetically), excited
species (48,52) localize. Thermal caging = rung-hopping barriers (dQ=+-2,+-4,+-6
at 2c/unit) on top of this free connectivity (see [[defect-travel]]).

Terminology: Q-rung = triangles with equal 2->3 dS. Component = free-reachable
piece. free-valence = #free slides from a triangle (NOT edge/vertex degree).
Scripts (scratchpad): reframe_2to3.py, clean_validate.py, free_coverage.py,
q_pattern.py. See [[worm-sampler-program]], [[fp-kinetics-findings]],
[[fpkmc-design]], [[flicker-background]].
