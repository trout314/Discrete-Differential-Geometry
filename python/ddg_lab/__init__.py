"""ddg_lab — shared research-grade machinery for the DDG experiment scripts.

The second tier of the 2026-08 script-library cleanup (see
notes/CLEANUP_PLAN.md): importable modules that many research CLIs share but
that are not yet stable enough for the `discrete_differential_geometry`
package proper. The old locations under scripts/, scripts/defect_dynamics/,
and tools/ are thin shims, so both `from ddg_lab import worm_slide` and the
legacy `import worm_slide` (with the scripts dirs on sys.path) work.

Modules (no eager imports here -- several load the D core on import):

  defect core   : defect_state, event_replay, crystal_flicker
  worm machinery: worm_moves, worm_helix, worm_slide, worm_deg4_slide,
                  tip_retract_search, dressed_generators, fk_moves
  f0 sector     : f0_worm, link_planner
  crystals      : tcp_melt, dopant_pairs, graft_signature
  geometry      : steiner_geodesic, heat_geodesic, sixhundred_cell
  seed pipeline : grid_sweep

Data paths inside these modules assume the source-tree layout (repo root
three levels up); this package is research code, not a distributable.
"""
