---
name: m2-only-gas
description: "m^2-only action (global pins + c_imp*sum m^2, NO legality/EDQ terms) holds a crystal-preserving dilute mobile defect gas on R AND C15; threshold ~0.4-0.5 both hosts; simplifies to 3 species by c_imp~0.7; melt-start shows glassy bistability"
metadata: 
  node_type: memory
  type: project
  originSessionId: d562feea-0887-419b-9cd9-9899d180c7c2
  modified: 2026-07-30T16:08:40.848Z
---

**Program (2026-07-30): how much action is needed for a defect gas?**
Action: 0.1(f3−f3c)² + 1.0(f1 − 6f3/e*)² + c_imp·Σ_v m(v)², nothing else
(no per-edge EDQ term, no U(n6) legality term). Scan script:
scratchpad/m2_penalty_scan.py (CELL/ETARGET/CIMPS/SWEEPS env).

**Globals-only baseline (c_imp=0)**: crystal melts <100 sw into the
max-entropy fluid — ONE percolated component, ~76% illegal, geometric
degree distribution (mode deg-3, ratio ~0.7), ZERO imp=0 vertices
(= random-chance legality, no grains, no boundary — single-phase).
IDENTICAL for R and C15 (bin-by-bin after volume scaling).

**m² scan from crystal (C15 flat target, R native+3q)**: fluid thins
monotonically 0.02→0.35; threshold ~0.4–0.5 BOTH hosts; at 0.5: stable
dilute mobile gas (n_ill ~50–90, comps ≤8 edges, maxm ≤4, turnover
>0.9 per 150–250 sw, ~13–19 acc moves/sw). At 0.7–0.85 simplifies to
THREE species: shared-vertex deg-4 pairs (4,4) DOMINANT, lone deg-4,
lone deg-3; multimers gone; (3,4,4) knot only transient. (4,4) is
m²-DISFAVORED (6c vs 4c for two monomers) → kinetic species (stripped
flicker remnant), echoing the λ-gas Z15-pair motif.

**Bistability**: melt-start at c=0.5 does NOT evaporate — sticks at
n_ill≈4450, one component, maxm 34 (cost ~580/vertex!), turnover 0.15
= glassy arrest; m² widens coexistence (consistent with EDQ-context
memory). m² PREVENTS condensation but cannot reverse it.

**e-target = density dial**: gas density monotone in (e_target −
e_native) tension; each 2→3 closes the f1-pin gap by 1−6/e* when
target above native. FLAT target = 2π/arccos(1/3) = 5.1042993 (zero
MEAN Regge deficit; user's canonical choice — "flat edge degree
target"). R native is ~1 f1-quantum from flat; C15 native ~9.

**Chemistry census enabled**: reaction_census.py --hinges-coef
--abs-couplings (commit e235ba8); first run on C15 c=0.7 flat target
(scratchpad/m2gas_chem.*). Cocycle-attached twin for viewer/catalog:
scratchpad/m2gas_c15_cimp0.5_coc_sw1500.* + data/viz pages (post-hoc
regen_cocycle FAILS on these gases — monodromy; attach from sweep 0).

Related: [[edq-only-melting]], [[c15-z14-dope-ensemble]],
[[maximalist-recording]].
