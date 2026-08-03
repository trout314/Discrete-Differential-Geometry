---
name: vdv-hdv-conflicts-with-tcp
description: Degree-VARIANCE penalties (VDV/HDV) conflict with the TCP n6 potential — leave them OFF for FK/TCP work
metadata: 
  node_type: memory
  type: feedback
  originSessionId: 1ed7f5d0-df5c-4552-b65e-5060180eafe2
  modified: 2026-07-21T06:05:36.516Z
---

The degree-VARIANCE objective terms conflict with the TCP/Frank-Kasper-encouraging
part of the objective and should generally be LEFT OFF whenever doing FK/TCP work
(anything using `set_n6_potential`: doping, melting, defect studies). Aaron
flagged this 2026-07-21.

**The two variance coefs** (SamplerParams): `hinge_degree_variance_coef` (HDV,
edge-degree variance, default **0.2**) and `codim3_degree_variance_coef` (VDV,
vertex-degree variance, default **0.1**). Both penalize spread of degrees toward
the mean.

**Why they conflict:** FK/TCP structure REQUIRES a spread of degrees — vertices at
12/14/15/16 (Z12/Z14/Z15/Z16), edges at both 5 and 6. A variance penalty pulls all
degrees toward the mean (~5.1 edges, ~13.4 verts), directly opposing the n6
potential that wants the discrete FK classes. So VDV/HDV fight the TCP objective.

**How to apply:** when building SamplerParams for any FK/TCP run, explicitly set
`hinge_degree_variance_coef=0.0` AND `codim3_degree_variance_coef=0.0` (do NOT rely
on defaults — the defaults are 0.2/0.1, ON). `run_point` in `tcp_melt.py` does this
correctly. NOTE: `dope_hold.py` and the scratchpad `c15_dope_run.py` did NOT — they
left the defaults on, so the earlier C15/R doping experiments in the
[[c15-z14-dope-ensemble]] line were run with VDV+HDV ON. RE-CHECKED 2026-07-21 with
variance OFF: the C15 raise/lower results are ROBUST (VDV had negligible effect —
strong FK potential dominates the dilute C15 regime: raise Z15/Z14 3.8->4.4 same
motifs, lower identical 5-atom knot), and R's cimp=0 result was also robust. Only
R's strong-FK (zleg0.6 cimp1.0) raise/lower + harder-push results remain un-rechecked
(R is stiffer / more finely balanced, so more likely to shift) — redo those VDV-off
before trusting the numbers. Script c15_dope_run.py now zeros both variance coefs.

Distinct context where VDV is deliberately ON: the seed-generation campaign
(`equilibrium_vdv.py`, `anneal_vdv.py`) STUDIES equilibrium vertex-degree variance
as its objective — there VDV IS the point. The "leave off" rule is specifically for
FK/TCP-structure experiments driven by the n6 potential.
