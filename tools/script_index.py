#!/usr/bin/env python3
"""Generate scripts/INDEX.md: the classified manifest of every research script.

Mechanical columns (docstring, last commit, import fan-in/out) are recomputed
from the tree + git on every run; the curated columns (program, status) live in
the CLASSIFICATION table below.  Unclassified files are listed loudly at the
bottom of the manifest so drift is visible.

Status vocabulary (the cleanup-plan tiers, see notes/CLEANUP_PLAN.md):
  lib        promote into a package (ddg or ddg_lab) -- imported by other scripts
  shim       thin re-export / CLI shim; core already promoted into the ddg package
  active     current research surface, stays in scripts/
  tool       reusable reporting / viewing / construction front-end, stays
  validator  encodes exact invariants -> convert to pytest (tests/)
  dormant    program paused with results in; stays in place, tagged, untested
  closed     program scientifically closed -> archive to legacy/<program>/

Usage:  python3 tools/script_index.py           # writes scripts/INDEX.md
        python3 tools/script_index.py --check   # exit 1 if unclassified files
"""

import ast
import os
import subprocess
import sys
from collections import defaultdict

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

ROOTS = ["scripts", "scripts/defect_dynamics", "tools"]

# ---------------------------------------------------------------------------
# Curated classification: relkey -> (program, status)
# relkey is the bare filename for scripts/ and tools/, "dd/<name>" for
# scripts/defect_dynamics/.  Keep sorted by program for reviewability.
# ---------------------------------------------------------------------------
CLASSIFICATION = {
    # --- seed-pipeline: production seed generation & library maintenance ----
    "equilibrium_vdv.py":        ("seed-pipeline", "active"),
    "grow_seed.py":              ("seed-pipeline", "active"),
    "anneal_vdv.py":             ("seed-pipeline", "active"),   # kept for comparison vs equilibrium
    "explore_grid.py":           ("seed-pipeline", "active"),
    "produce_grid.py":           ("seed-pipeline", "active"),
    "extend_library.py":         ("seed-pipeline", "active"),
    "seed_queue_runner.py":      ("seed-pipeline", "active"),
    "overnight_large_queue.sh":  ("seed-pipeline", "active"),
    "compare_twins.py":          ("seed-pipeline", "active"),
    "quench.py":                 ("seed-pipeline", "dormant"),
    "grid_sweep.py":             ("seed-pipeline", "lib"),
    "seed_utils.py":             ("seed-pipeline", "shim"),     # -> ddg.seed_utils
    "reburn_batch.py":           ("seed-pipeline", "active"),
    "reburn_family.py":          ("seed-pipeline", "active"),
    "reburn_seed.py":            ("seed-pipeline", "active"),

    # --- analysis-chain: distance -> roundness -> curvature ----------------
    "distance_distribution.py":  ("analysis-chain", "active"),
    "roundness_analysis.py":     ("analysis-chain", "active"),
    "scale_curvature.py":        ("analysis-chain", "active"),

    # --- symmetry-census: exact Aut(K) orbit machinery ----------------------
    "crystal_grains.py":         ("symmetry-census", "shim"),   # -> ddg.crystal_grains
    "move_site_census.py":       ("symmetry-census", "active"),
    "bc_chain_census.py":        ("symmetry-census", "active"),
    "orbit_decoration_census.py":("symmetry-census", "active"),
    "chain_select.py":           ("symmetry-census", "shim"),   # -> ddg.chain_select

    # --- fk-census: FK / disclination-skeleton classification ---------------
    "fk_skeleton.py":            ("fk-census", "active"),       # core -> ddg.fk_skeleton
    "link_classes.py":           ("fk-census", "shim"),         # -> ddg.link_classes
    "defect_census.py":          ("fk-census", "tool"),

    # --- references: TCP crystal construction + melt spectroscopy -----------
    "tcp_reference.py":          ("references", "active"),      # core -> ddg.tcp_reference
    "tcp_melt.py":               ("references", "lib"),
    "reference_campaign.py":     ("references", "active"),
    "reference_summary.py":      ("references", "active"),

    # --- crystal-library: the cross-crystal catalog --------------------------
    "crystal_library_data.py":   ("crystal-library", "active"),
    "crystal_library_page.py":   ("crystal-library", "tool"),

    # --- cocycle: T^3 integer-cocycle instrumentation ------------------------
    # helpers are a de-facto library for 23 scripts; the validation ladder
    # itself becomes pytest.
    "cocycle_check.py":          ("cocycle", "active"),         # helper -> ddg.tcp_reference
    "dd/regen_cocycle.py":       ("cocycle", "tool"),

    # --- fk-amorphous: the ACTIVE amorphous-FK program -----------------------
    "fk_amorphous.py":           ("fk-amorphous", "active"),
    "fk_anneal.py":              ("fk-amorphous", "active"),
    "fk_worm_harvest.py":        ("fk-amorphous", "active"),
    "contract_relax.py":         ("fk-amorphous", "active"),
    "validate_contract_split.py":("fk-amorphous", "validator"),
    "dd/amorph_census.py":       ("fk-amorphous", "active"),
    "dd/fk_channel_census.py":   ("fk-amorphous", "active"),
    "dd/fk_amorphous_ab_figure.py": ("fk-amorphous", "active"),

    # --- fk-move-search: FK->FK block moves / glass producer -----------------
    "dd/fk_moves.py":            ("fk-move-search", "lib"),
    "dd/glass_mine.py":          ("fk-move-search", "dormant"),
    "dd/glass_quench.py":        ("fk-move-search", "dormant"),
    "dd/enumerate_fillings.py":  ("fk-move-search", "dormant"),

    # --- grafting: decorated-boundary lump surgery ---------------------------
    "graft_signature.py":        ("grafting", "lib"),
    "graft_ball_search.py":      ("grafting", "active"),
    "graft_c15_cross.py":        ("grafting", "active"),
    "graft_c36_control.py":      ("grafting", "active"),

    # --- defect-core: identification, bookkeeping, viewing -------------------
    "dd/defect_state.py":        ("defect-core", "lib"),
    "dd/event_replay.py":        ("defect-core", "lib"),
    "dd/crystal_flicker.py":     ("defect-core", "lib"),
    "dd/defect_viewer.py":       ("defect-core", "tool"),
    "dd/defect_catalog.py":      ("defect-core", "tool"),
    "dd/species_report.py":      ("defect-core", "tool"),
    "dd/knot_worldlines.py":     ("defect-core", "tool"),
    "dd/centroid_spatial.py":    ("defect-core", "tool"),
    "dd/complex_d3_census.py":   ("defect-core", "dormant"),

    # --- worm-core: exact move combinatorics + knot-slide machinery ----------
    "dd/worm_moves.py":          ("worm-core", "lib"),
    "dd/worm_helix.py":          ("worm-core", "lib"),
    "dd/worm_slide.py":          ("worm-core", "lib"),
    "dd/dressed_generators.py":  ("worm-core", "lib"),
    "dd/crossval_moves.py":      ("worm-core", "validator"),
    "dd/worm_crossval2.py":      ("worm-core", "validator"),
    "dd/samerung_validate.py":   ("worm-core", "validator"),
    "dd/fusion_verify.py":       ("worm-core", "validator"),
    "dd/worm_cycles.py":         ("worm-core", "validator"),
    "dd/worm_catalog.py":        ("worm-core", "dormant"),

    # --- deg4-transport: the deg-4 worm / translation program (delivered) ----
    "dd/worm_deg4_slide.py":     ("deg4-transport", "lib"),
    "dd/worm_deg4.py":           ("deg4-transport", "dormant"),
    "dd/worm_walk.py":           ("deg4-transport", "dormant"),
    "dd/worm_helix_launch": ("deg4-transport", "dormant"),  # placeholder, removed below if absent
    "dd/tip_retract_search.py":  ("deg4-transport", "dormant"),
    "dd/search_positive_control.py": ("deg4-transport", "dormant"),
    "dd/template_census.py":     ("deg4-transport", "dormant"),
    "dd/hinge_by_position.py":   ("deg4-transport", "dormant"),
    "dd/translation_episodes.py":("deg4-transport", "dormant"),
    "dd/deg4_moves.py":          ("deg4-transport", "dormant"),
    "dd/deg4_observable.py":     ("deg4-transport", "dormant"),
    "dd/deg4_pair_census.py":    ("deg4-transport", "dormant"),
    "dd/deg4_provenance.py":     ("deg4-transport", "dormant"),
    "dd/deg4_density_scan.py":   ("deg4-transport", "dormant"),
    "dd/launch_worm_campaign.sh":("deg4-transport", "dormant"),

    # --- f0-sector: vertex-number surgery (ACTIVE) ---------------------------
    "dd/f0_worm.py":             ("f0-sector", "lib"),
    "dd/f0_channel.py":          ("f0-sector", "active"),
    "dd/edge_removal.py":        ("f0-sector", "active"),
    "dd/harvest_f0.py":          ("f0-sector", "active"),
    "dd/harvest_f0_verdict.py":  ("f0-sector", "active"),
    "dd/link_planner.py":        ("f0-sector", "lib"),
    "dd/planner_profile.py":     ("f0-sector", "tool"),
    "dd/chainsites_validate.py": ("f0-sector", "validator"),

    # --- chord-bilocal: strict chord channel + bilocal carriers (CLOSED:
    #     5 impossibility results, transport blocked; see bilocal-program-saga)
    "dd/twosided_chord.py":      ("chord-bilocal", "validator"),
    "dd/agg_knobs_test.py":      ("chord-bilocal", "validator"),

    # --- fpkmc: paused pending D-side optimization (D-first policy) ----------
    "dd/fpkmc_m0_derivations.py":("fpkmc", "validator"),
    "dd/fpkmc_v2_nu.py":         ("fpkmc", "dormant"),
    "dd/fpkmc_v3_hb.py":         ("fpkmc", "dormant"),
    "dd/fpkmc_v3b_db.py":        ("fpkmc", "validator"),
    "dd/fpkmc_v4_fp.py":         ("fpkmc", "dormant"),
    "dd/fpkmc_v5_contact.py":    ("fpkmc", "dormant"),
    "dd/fp_transport.py":        ("fpkmc", "dormant"),
    "dd/fp_encounter.py":        ("fpkmc", "dormant"),
    "dd/fp_prod_report.py":      ("fpkmc", "dormant"),
    "dd/fp_dock_census_intrinsic.py": ("fpkmc", "dormant"),
    "dd/fp_dock_local_census.py":("fpkmc", "dormant"),
    "dd/fp_recombine_intrinsic.py": ("fpkmc", "dormant"),
    "dd/fp_holonomy_decorrelation.py": ("fpkmc", "dormant"),
    "dd/fp_intrinsic_msd.py":    ("fpkmc", "dormant"),
    "dd/run_recomb_campaign.sh": ("fpkmc", "dormant"),

    # --- intrinsic-geometry: development / parallel-transport toolkit tests --
    "dd/transport_tests.py":     ("intrinsic-geometry", "validator"),
    "dd/transport_injectivity_check.py": ("intrinsic-geometry", "validator"),

    # --- crystal-gas: library-wide dilute defect gas (ACTIVE) ----------------
    "dd/crystal_gas_scan.py":    ("crystal-gas", "active"),
    "dd/crystal_gas_sweep.py":   ("crystal-gas", "active"),
    "dd/crystal_gas_catalog.py": ("crystal-gas", "active"),
    "dd/crystal_gas_report.py":  ("crystal-gas", "active"),
    "dd/cimp_scan.py":           ("crystal-gas", "active"),
    "dd/cimp_scan_figure.py":    ("crystal-gas", "active"),
    "dd/m_histogram.py":         ("crystal-gas", "active"),

    # --- mobile-gas: constrained knot liquid (results in) --------------------
    "dd/mobile_gas.py":          ("mobile-gas", "dormant"),
    "dd/mgas_analyze.py":        ("mobile-gas", "dormant"),
    "dd/mobility_sweep.py":      ("mobile-gas", "dormant"),
    "dd/knot_migration.py":      ("mobile-gas", "dormant"),

    # --- quanta-strain: strain/heal quanta program (ACTIVE) ------------------
    "dd/quanta_heal.py":         ("quanta-strain", "active"),
    "dd/quanta_heal_report.py":  ("quanta-strain", "active"),
    "dd/quanta_species.py":      ("quanta-strain", "active"),
    "dd/certify_strain_gas.py":  ("quanta-strain", "active"),
    "dd/strain_kinematics.py":   ("quanta-strain", "active"),

    # --- hu-statics: hyperuniformity / structure factors / screening ---------
    "dd/defect_density_hu.py":   ("hu-statics", "active"),
    "dd/defect_density_hu_figure.py": ("hu-statics", "active"),
    "dd/carrier_gr.py":          ("hu-statics", "active"),
    "dd/sl_verdict.py":          ("hu-statics", "active"),
    "dd/cross_spectrum.py":      ("hu-statics", "active"),
    "dd/curvature_balls.py":     ("hu-statics", "active"),
    "dd/curv_scale_real.py":     ("hu-statics", "active"),
    "dd/sk_exponent.py":         ("hu-statics", "active"),
    "dd/tt_channel.py":          ("hu-statics", "active"),
    "dd/defect_statics.py":      ("hu-statics", "active"),
    "sk_torus.py":               ("hu-statics", "active"),
    "hyperuniformity.py":        ("hu-statics", "dormant"),
    "curvature_hyperuniformity_g.py": ("hu-statics", "dormant"),
    "w_liquid_ladder.py":        ("hu-statics", "dormant"),

    # --- six-web: emergent-gauge / disclination-web program ------------------
    "dd/web_s66.py":             ("six-web", "dormant"),
    "dd/web_vector.py":          ("six-web", "dormant"),
    "dd/web_transport.py":       ("six-web", "dormant"),
    "flip_analysis.py":          ("six-web", "dormant"),
    "line_tension.py":           ("six-web", "dormant"),

    # --- flicker: the (3,4,4)-knot flicker spectroscopy (results in) ---------
    "dd/flicker_fraction.py":    ("flicker", "dormant"),
    "dd/flicker_ladder.py":      ("flicker", "dormant"),
    "dd/flicker_quantum.py":     ("flicker", "dormant"),
    "dd/flicker_site.py":        ("flicker", "dormant"),
    "dd/flicker_spectrum.py":    ("flicker", "dormant"),
    "dd/flicker_subclass.py":    ("flicker", "dormant"),

    # --- run5h-passes: the run5h post-hoc analysis suite (CLOSED) ------------
    "dd/pass2_structure.py":     ("run5h-passes", "active"),  # recently reused for amorph census

    # --- reaction-census: thermal merge/split chemistry (results in) ---------
    "dd/reaction_census.py":     ("reaction-census", "dormant"),
    "dd/reaction_report.py":     ("reaction-census", "dormant"),
    "dd/species_interactions.py":("reaction-census", "dormant"),
    "dd/run_reaction_campaign.sh": ("reaction-census", "dormant"),

    # --- doping: R/C15 dopant solubility campaign (results in) ---------------
    "dopant_pairs.py":           ("doping", "lib"),          # 17 importers -- extract helpers
    "dope_hold.py":              ("doping", "dormant"),
    "complex_analysis.py":       ("doping", "dormant"),
    "defect_pairs.py":           ("doping", "dormant"),
    "r_solubility.py":           ("doping", "dormant"),
    "r_isotherm.py":             ("doping", "dormant"),
    "replica_hold.py":           ("doping", "dormant"),
    "replica_aggregate.py":      ("doping", "dormant"),

    # --- lapse-adm: combinatorial lapse / ADM constraint program -------------
    "adm_constraint.py":         ("lapse-adm", "dormant"),
    "local_lapse.py":            ("lapse-adm", "dormant"),
    "lapse_Nsweep.py":           ("lapse-adm", "dormant"),
    "lapse_g_sweep.py":          ("lapse-adm", "dormant"),

    # --- edq: EDQ-only melting / certification campaign (results in) ---------
    "certify_edq_clust.py":      ("edq", "dormant"),

    # --- geometry: geodesics + exact constructions ---------------------------
    "steiner_geodesic.py":       ("geometry", "lib"),
    "heat_geodesic.py":          ("geometry", "tool"),
    "sixhundred_cell.py":        ("geometry", "tool"),
    "defect_boundary_map.py":    ("geometry", "tool"),
    # TODO(generalize): currently hard-wired to the 2->3 lens study; make it a
    # general lens/optics front-end when we next touch it (see CLEANUP_PLAN.md).
    "defect_lens_field.py":      ("geometry", "tool"),

    "r_pinscan_queue.sh":        ("doping", "dormant"),

    # --- infra ---------------------------------------------------------------
    "link_memory.sh":            ("infra", "tool"),
    "script_index.py":           ("infra", "tool"),
}
# remove placeholder entries for files that don't exist
CLASSIFICATION = {k: v for k, v in CLASSIFICATION.items() if not k.endswith("worm_helix_launch")}

STATUS_ORDER = ["lib", "shim", "active", "tool", "validator", "dormant", "closed"]

STATUS_LEGEND = {
    "lib":       "promote into a package (imported by other scripts)",
    "shim":      "thin re-export/CLI shim; core promoted into the ddg package",
    "active":    "current research surface; stays in `scripts/`",
    "tool":      "reusable reporting/viewing front-end; stays",
    "validator": "encodes exact invariants; convert to pytest",
    "dormant":   "program paused, results recorded; stays in place, tagged",
    "closed":    "program scientifically closed; archive to `legacy/<program>/`",
}


def relkey(root, fname):
    return ("dd/" + fname) if root.endswith("defect_dynamics") else fname


def first_doc_line(path):
    if path.endswith(".sh"):
        with open(path) as fh:
            for line in fh:
                s = line.strip()
                if s.startswith("#") and not s.startswith("#!"):
                    return s.lstrip("# ").strip()
        return ""
    try:
        doc = ast.get_docstring(ast.parse(open(path).read()))
        return (doc or "").strip().split("\n")[0]
    except Exception:
        return "(parse error)"


def git_last_dates():
    out = subprocess.run(
        ["git", "log", "--format=%ad|%H", "--date=short", "--name-only"],
        capture_output=True, text=True, cwd=_ROOT).stdout
    last, cur = {}, None
    for line in out.splitlines():
        if "|" in line and len(line.split("|")[0]) == 10:
            cur = line.split("|")[0]
        elif line.strip() and cur:
            last.setdefault(line.strip(), cur)
    return last


def collect():
    files = {}   # relkey -> abs path
    for r in ROOTS:
        absr = os.path.join(_ROOT, r)
        for f in sorted(os.listdir(absr)):
            if f.endswith((".py", ".sh")) and not f.startswith("__"):
                files[relkey(r, f)] = os.path.join(absr, f)
    return files


def import_graph(files):
    modnames = {os.path.basename(p)[:-3]: k for k, p in files.items() if p.endswith(".py")}
    uses = defaultdict(set)      # relkey -> set of relkeys it imports
    for k, p in files.items():
        if not p.endswith(".py"):
            continue
        try:
            tree = ast.parse(open(p).read())
        except Exception:
            continue
        for node in ast.walk(tree):
            names = []
            if isinstance(node, ast.Import):
                names = [a.name.split(".")[0] for a in node.names]
            elif isinstance(node, ast.ImportFrom) and node.module:
                names = [node.module.split(".")[0]]
            for n in names:
                if n in modnames and modnames[n] != k:
                    uses[k].add(modnames[n])
    used_by = defaultdict(set)
    for k, deps in uses.items():
        for d in deps:
            used_by[d].add(k)
    return uses, used_by


def main():
    check = "--check" in sys.argv
    files = collect()
    uses, used_by = import_graph(files)
    dates = git_last_dates()

    unclassified = [k for k in files if k not in CLASSIFICATION]
    stale = [k for k in CLASSIFICATION if k not in files]

    by_program = defaultdict(list)
    for k in files:
        prog, status = CLASSIFICATION.get(k, ("UNCLASSIFIED", "?"))
        by_program[prog].append((k, status))

    lines = []
    lines.append("# Research-script manifest")
    lines.append("")
    lines.append("Auto-generated by `tools/script_index.py` -- do not edit the tables by")
    lines.append("hand; edit the CLASSIFICATION map in the generator and re-run.")
    lines.append("")
    lines.append("## Status legend")
    lines.append("")
    for s in STATUS_ORDER:
        lines.append(f"- **{s}** -- {STATUS_LEGEND[s]}")
    lines.append("")

    # summary
    counts = defaultdict(int)
    for k in files:
        counts[CLASSIFICATION.get(k, ("?", "?"))[1]] += 1
    lines.append("## Summary")
    lines.append("")
    lines.append("| status | files |")
    lines.append("|---|---|")
    for s in STATUS_ORDER:
        lines.append(f"| {s} | {counts.get(s, 0)} |")
    if unclassified:
        lines.append(f"| **UNCLASSIFIED** | {len(unclassified)} |")
    lines.append(f"| **total** | {len(files)} |")
    lines.append("")

    def relpath_of(k):
        return os.path.relpath(files[k], _ROOT)

    for prog in sorted(by_program):
        rows = sorted(by_program[prog],
                      key=lambda ks: (STATUS_ORDER.index(ks[1]) if ks[1] in STATUS_ORDER else 99, ks[0]))
        lines.append(f"## {prog}  ({len(rows)})")
        lines.append("")
        lines.append("| script | status | last commit | fan-in | what it is |")
        lines.append("|---|---|---|---|---|")
        for k, status in rows:
            p = relpath_of(k)
            doc = first_doc_line(files[k]).replace("|", "\\|")[:110]
            date = dates.get(p, "?")
            fan_in = len(used_by.get(k, ()))
            fi = str(fan_in) if fan_in else ""
            lines.append(f"| `{p}` | {status} | {date} | {fi} | {doc} |")
        lines.append("")

    if unclassified:
        lines.append("## UNCLASSIFIED (fix the generator!)")
        lines.append("")
        for k in sorted(unclassified):
            lines.append(f"- `{relpath_of(k)}`")
        lines.append("")
    if stale:
        lines.append("## Stale classification entries (file gone)")
        lines.append("")
        for k in sorted(stale):
            lines.append(f"- `{k}`")
        lines.append("")

    out = os.path.join(_ROOT, "scripts", "INDEX.md")
    with open(out, "w") as fh:
        fh.write("\n".join(lines))
    print(f"wrote {out}: {len(files)} scripts, "
          f"{len(unclassified)} unclassified, {len(stale)} stale entries")
    if check and (unclassified or stale):
        sys.exit(1)


if __name__ == "__main__":
    main()
