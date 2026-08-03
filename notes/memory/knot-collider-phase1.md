---
name: knot-collider-phase1
description: "Phase-1 collider DONE — V(s)=0 to 1e-14 until chord-sharing contact (s=4, the same-helix minimum); contact is repulsive (+4.4 to +12.8 lam=1) and reversible; WASHBOARD measured (17-51, explains caging); merged species (9,23,23,8) universal, NOT the decamer"
metadata: 
  node_type: memory
  type: project
  originSessionId: d562feea-0887-419b-9cd9-9899d180c7c2
  modified: 2026-07-25T23:38:33.840Z
---

Measured 2026-07-25. Tool: `scripts/defect_dynamics/knot_collider.py` —
deterministic two-knot collider on one BC helix, no sampling. Prep exploits
the no-halo theorem: create A on a chain window, slide +8 (directed clean
slides via worm_slide frames, steered by exact chain arithmetic from
worm_helix.bc_orbit), create B at A's origin (site exactly pristine again),
walk A back. Host pristine R m4; orbit L=3252 chain steps, winding (43,29,8),
chain REVISITS vertices (2928 distinct) so track j explicitly, never infer
from vertex->index. JSON scratchpad `collider_v2.json`; figure
`data/figs/knot_collider_phase1.png`.

**V(s) = 0 to 1e-14 for every s >= 8 chain steps, all 4 sites** — the
no-halo theorem verified dynamically to machine precision. CRITICAL
methodology: V(s) must be washboard-corrected — V = S(A,B) − S_single(A at
its CURRENT site) − S_single(B), with S_single recorded during the outbound
walk. Subtracting 2x creation cost instead gives spurious V ~ -30 (the
knot's site energy varies along the helix).

**The WASHBOARD:** single-knot energy S_single varies 17.1–51.0 (λ=1 units)
from site to site along the helix — rung-70 (3,4,4) sites are cheapest
(18.3), rung-74 most expensive (51.0). At λ=0.4 that is barriers of up to
~13 action units — this quantitatively explains the caging of thermal knots
([[defect-travel]]) and qualitatively the rung-population reweighting
([[flicker-background]]). Transport physics = washboard hopping.

**Contact (s=4, chords sharing one vertex — the minimum same-helix
separation):** merge into universal shape f=(9,23,23,8) (9 = 5+5−1), BOTH
deg-3 chords intact, sumZ 128–130. V_contact = +4.4 (rung70+70) to +12.8
(rung74 pair): REPULSIVE, rung-dependent — first entries of the S-matrix.
unfuses=True everywhere: one clean slide re-separates into two decorated
knots. So same-helix contact is a reversible repulsive doorstep, NOT a trap,
and NOT the decamer — the immortal (10,22,18,5) is not a same-helix
chord-sharing pair; it must form off-chain or thermally.

Gotcha fixed: the unfuse test must accept DECORATED knot sigs (exactly one
3 each), not bare (3,4,4).

Next (Phase 2): S-matrix over many sites and rung pairs; off-chain contact
geometries (does the decamer form?); Phase 3 thermal kinetics with slide
channel on. See [[no-halo-verdict]], [[worm-sampler-program]].
