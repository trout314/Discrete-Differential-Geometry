"""The longest-lived defect complexes and their edge-degree decorated
structures.

Lifetimes come from the 8-chain reaction-census tracks (species = modal
induced-shape f-vector + coordination).  Structures come from the gas
snapshots via DefectState: per complex the EXACT decorated-isomorphism class
(canonical_key of the induced subcomplex, ecolor = ambient edge degree,
vcolor = n6), joined to the track species by (f-vector, coord).

CERTIFICATION.  canonical_key returns (key, exact).  Only exact=True makes the
key a canonical form, where key equality <=> isomorphism BOTH ways.  When the
search is truncated the minimum is taken over a prefix of orderings that
depends on the input vertex LABELS, so the key stops being a function of the
isomorphism class: equal keys still prove isomorphism, but two isomorphic
complexes can get DIFFERENT keys and be counted as two classes.  A table of
"decorated-isomorphism classes" built from truncated keys therefore
OVER-SPLITS and its class count is meaningless.

So uncertified complexes are excluded from the class table and reported
separately rather than silently inflating it, and the script exits non-zero
when any occur.  Raise the cap with DDG_CANON_LIMIT to try to certify them
(the cost is factorial in the size of the largest colour class).  Measured
2026-08-02 on both gas snapshots: 44 complexes, 0 uncertified.
"""
import glob, json, os, sys
from collections import Counter, defaultdict
import numpy as np

_ROOT = os.path.normpath(os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "..", ".."))
sys.path.insert(0, os.path.join(_ROOT, "python"))
sys.path.insert(0, os.path.join(_ROOT, "scripts/defect_dynamics"))
import discrete_differential_geometry as ddg
import defect_state as ds

SNAPS = [os.path.join(_ROOT, "data/mgas/lam35r_snap15000.mfd"),
         os.path.join(_ROOT, "data/mgas/lam35r_m4_over.mfd")]

# ---- Part A: lifetime by species from the campaign tracks ------------------
life = defaultdict(list)      # (f, coord) -> [life_sw...]
cens = Counter()              # (f, coord) -> censored (still alive at end)
for p in sorted(glob.glob(os.path.join(_ROOT, "data/rxn_lam035_m4/c?.tracks.jsonl"))):
    with open(p) as fh:
        for line in fh:
            if not line.strip():
                continue
            r = json.loads(line)
            sh = r.get("shape")
            if not sh:
                continue
            k = (tuple(sh), r.get("coord"))
            life[k].append(r["life_sw"])
            if r.get("ended") not in ("death", "merge"):
                cens[k] += 1

# ---- Part B: decorated structures from the snapshots -----------------------
KEY_LIMIT = int(os.environ.get("DDG_CANON_LIMIT", 200000))

classes = {}      # canonical key -> dict(count, sig, f, coord, nodes, key)
uncertified = []  # complexes whose key was truncated -- NOT isomorphism classes
for snap in SNAPS:
    m = ddg.Manifold.load(snap, 3)
    st = ds.DefectState(m)
    for cx in st.components():
        sh = st.induced_shape(cx.verts)
        f = tuple(int(x) for x in sh["f"])
        coord = int(st.total_coordination(cx.verts))
        fac = st.induced_facets(cx.verts)
        vc, ec = st.decorations(cx.verts)
        key, exact = ds.canonical_key(fac, vc, ec, limit=KEY_LIMIT)
        if not exact:
            uncertified.append(dict(nv=len(cx.verts), sig=cx.sig, f=f,
                                    coord=coord,
                                    snap=os.path.basename(snap)))
            continue
        c = classes.setdefault(key, dict(
            count=0, sig=cx.sig, nodes=cx.nodes, f=f, coord=coord, key=key))
        c["count"] += 1

# ---- join + report ---------------------------------------------------------
rows = []
for c in classes.values():
    k = (c["f"], c["coord"])
    lv = life.get(k, [])
    rows.append((c, lv))

def med(x): return float(np.median(x)) if x else float("nan")
def mx(x): return float(np.max(x)) if x else float("nan")

# order by species median lifetime (then count)
rows.sort(key=lambda r: (-(med(r[1]) if r[1] else -1), -r[0]["count"]))

print(f"{len(classes)} CERTIFIED decorated-isomorphism classes across "
      f"{sum(c['count'] for c in classes.values())} complexes "
      f"({', '.join(os.path.basename(s) for s in SNAPS)})")
if uncertified:
    print(f"\n  !! {len(uncertified)} complex(es) EXCLUDED: canonical_key "
          f"truncated at limit={KEY_LIMIT}, so their keys are not functions "
          f"of the\n     isomorphism class and would over-split the table "
          f"below. Sizes: "
          f"{sorted(u['nv'] for u in uncertified)}")
    print(f"     Re-run with a larger DDG_CANON_LIMIT to certify them "
          f"(cost is factorial in the largest colour class).")
else:
    print(f"  all complexes certified (exact isomorphism keys)")
print()
print(f"{'#':>3s} {'cnt':>4s} {'sig (illegal degs)':22s} {'f-vector':14s} "
      f"{'Z':>4s} {'n6 nodes':10s} {'tracks':>7s} {'med life':>9s} "
      f"{'max life':>9s} {'alive@end':>9s}")
for i, (c, lv) in enumerate(rows[:24]):
    k = (c["f"], c["coord"])
    print(f"{i:3d} {c['count']:4d} {str(c['sig']):22s} {str(c['f']):14s} "
          f"{c['coord']:4d} {str(c['nodes']):10s} {len(lv):7d} "
          f"{med(lv):9.1f} {mx(lv):9.1f} {cens.get(k, 0):9d}")

# decorated exemplars for the top classes: canonical relabeled edge list
print("\n--- decorated 1-skeletons (canonical labels; edges as i-j:deg, "
      "illegal degs UPPERCASE-marked with *) ---")
for i, (c, lv) in enumerate(rows[:10]):
    rel_f, rel_e, rel_v = c["key"]
    ill = [f"{u}-{w}:{d}*" for (u, w), d in rel_e if d not in (5, 6)]
    leg = [f"{u}-{w}:{d}" for (u, w), d in rel_e if d in (5, 6)]
    print(f"\n[{i}] x{c['count']}  sig={c['sig']}  f={c['f']}  "
          f"n6(vertices)={list(rel_v)}  [certified]")
    print(f"    illegal: {' '.join(ill) if ill else '(none)'}")
    print(f"    legal:   {' '.join(leg)}")

if uncertified:
    sys.exit(1)          # a truncated key must not pass silently in a pipeline
