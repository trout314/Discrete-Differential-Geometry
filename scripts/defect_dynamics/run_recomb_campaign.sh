#!/bin/bash
# Phase B: high-statistics FP recombination-by-intrinsic-class campaign.
#
# One episode per process: the D slide-graph scan leaks ~12 MB/scan (R6),
# so a process that runs many flights walks into GB-scale RSS.  Fresh
# process per episode caps peak memory at ~1 GB and costs ~5 s startup.
#
# Host is m2 deliberately: the exact intrinsic classes were shown to be
# HOST-INDEPENDENT (m4 collider crossings land in the same classes as m2
# FP docks), and m2 runs ~8x faster, so m2 buys ~8x the statistics for
# the same physics.  Frozen background, slide channel only.
#
# usage: run_recomb_campaign.sh <worker-id> <deadline-seconds>
set -u
W=${1:?worker id}
DEADLINE=${2:-25200}
ROOT="/Users/atrout/Desktop/Discrete-Differential-Geometry"
OUT="$ROOT/data/fpkmc/prodD"
mkdir -p "$OUT"
cd "$ROOT"

i=0
while [ $SECONDS -lt $DEADLINE ]; do
  i=$((i + 1))
  seed=$((20000 + W * 100000 + i))
  python scripts/defect_dynamics/fp_encounter.py \
      --mode recombine \
      --ref data/tcp_reference/T3_R_m2_N7248.mfd \
      --lam 0.40 --sep-sites 3 --episodes 1 --max-flights 100 \
      --seed "$seed" --out "$OUT/w${W}_e${i}.json" \
      >> "$OUT/w${W}.log" 2>&1
done
echo "[w$W] done: $i episodes in ${SECONDS}s" >> "$OUT/w${W}.log"
