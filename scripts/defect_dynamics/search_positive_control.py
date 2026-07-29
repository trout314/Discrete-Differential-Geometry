"""POSITIVE CONTROL for tip_retract_search: plant the validated 5-move CREATE
(worm_deg4.CREATE_SEQ) on the pristine m2 crystal -- a state where a 5-move
full retract PROVABLY exists (the exact reverse, entirely within a radius-2
core of tip edge (163,358)).  Run the same search with GOAL=empty; it MUST
find a depth-5 solution or the search logic is buggy."""
import os, sys
os.environ.setdefault("GOAL", "empty")
os.environ.setdefault("MAXD", "5")
os.environ.setdefault("ILLX", "3")

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
_ROOT = os.path.normpath(os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "..", ".."))
sys.path.insert(0, os.path.join(_ROOT, "scripts/defect_dynamics"))

import numpy as np
import tip_retract_search as T
from worm_deg4 import CREATE_SEQ, apply_seq, REF_CELL
import discrete_differential_geometry as ddg

m = ddg.Manifold.load(REF_CELL, 3)
apply_seq(m, CREATE_SEQ)
F = [tuple(sorted(int(x) for x in t)) for t in np.asarray(m.facets())]
tip = (163, 358)                      # one of the planted deg-4 edges
Vball = T.ball(F, tip, T.R_BALL)
Vcore = frozenset(T.ball(F, tip, T.R_CORE))
region = [t for t in F if set(t) <= Vball]
P = T.Patch([frozenset(t) for t in region])
sols = T.search(P, Vcore, tip, "POSITIVE CONTROL (planted 5-edge cluster)")
if sols:
    print(f"\nCONTROL PASS: search found {len(sols)} full-retract "
          f"solution(s) at depth {len(sols[0][0])}")
else:
    print("\nCONTROL FAIL: search missed a provably-existing retract -- "
          "the search logic is buggy")
    sys.exit(1)
