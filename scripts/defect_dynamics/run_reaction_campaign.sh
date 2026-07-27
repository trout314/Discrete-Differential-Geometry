#!/bin/bash
# Phase A: thermal reaction census fleet.
#
# Two levers, crossed:
#   DENSITY   lam=0.40 (~8 complexes in the m4 box) vs lam=0.35 (~22).
#             Merge is bimolecular, split unimolecular, so k_merge/k_split
#             must scale differently with n -- the test that "reaction"
#             is the right language at all.
#   TRANSPORT slide channel off vs on.  The slide move is certified to
#             satisfy detailed balance (V3b), so it changes the kinetics
#             and NOT the stationary ensemble: a rate that moves means
#             the chemistry is transport-limited, one that does not means
#             reaction-limited.
#
# Each worker starts from a distinct equilibrated mgas snapshot.  Note
# the runs begin recording immediately; the census time series carries
# n_illegal / n_components per chunk, so an initial re-equilibration
# window (if the snapshot's etarget differed) is visible and can be
# dropped offline.
set -u
ROOT="/Users/atrout/Desktop/Discrete-Differential-Geometry"
CELL="$ROOT/data/tcp_reference/T3_R_m4_N57984.mfd"
OUT="$ROOT/data/reactions"
DEADLINE=${1:-25200}
mkdir -p "$OUT"
cd "$ROOT"

launch () {   # name lam slideprob startsnap seed
  caffeinate -i python scripts/defect_dynamics/reaction_census.py \
      --cell "$CELL" --start "data/mgas/$4" \
      --lam "$2" --slide-prob "$3" --etarget 5.105025 \
      --max-seconds "$DEADLINE" --chunk 25 --logmb 64 \
      --audit-every 40 --seed "$5" --out "$OUT/$1" \
      > "$OUT/$1.log" 2>&1 &
  echo "  launched $1 (lam=$2 slide=$3 start=$4 seed=$5) pid $!"
}

launch l40a 0.40 0.0  ab1_snap14000.mfd    1001
launch l40b 0.40 0.0  ab2_snap14000.mfd    1002
launch l40s 0.40 0.10 ab3_snap14000.mfd    1003
launch l35a 0.35 0.0  lam35_snap20000.mfd  1004
launch l35s 0.35 0.10 lam35_snap17000.mfd  1005
wait
