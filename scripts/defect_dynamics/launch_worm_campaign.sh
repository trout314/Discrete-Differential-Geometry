#!/bin/bash
# worms-on lambda=0.35 m4 campaign: mirror of the baseline slide-off arm
# (c4-c7 in data/rxn_lam035_m4) with the D-side deg-4 worm channel enabled.
# 8 chains, over-dispersed starts (4 dilute + 4 over-defected), 18000s each.
# Wrap the invocation in `caffeinate -i` (see CLAUDE.md run rules).
cd "$(dirname "$0")/../.." || exit 1
export DDG_NO_AUTOBUILD=1
OUT=data/rxn_lam035_m4_worm
CELL=data/tcp_reference/T3_R_m4_N57984.mfd
OVER=data/mgas/lam35r_m4_over.mfd
DILUTE=data/mgas/lam35r_snap15000.mfd
echo "worm campaign start $(date); out=$OUT; 18000s/chain"
mkdir -p "$OUT"
i=0
for spec in \
  "$DILUTE 201" "$OVER 202" "$DILUTE 203" "$OVER 204" \
  "$DILUTE 205" "$OVER 206" "$DILUTE 207" "$OVER 208"; do
  set -- $spec
  START=$1; SEED=$2
  PYTHONPATH=python python scripts/defect_dynamics/reaction_census.py \
    --cell $CELL --lam 0.35 --start "$START" --seed $SEED \
    --slide-prob 0.0 --worm-prob 1e-4 --max-seconds 18000 \
    --out $OUT/w$i > $OUT/w$i.launch.log 2>&1 &
  echo "  launched w$i (start=$(basename $START) worm=1e-4 seed=$SEED) pid $!"
  i=$((i+1))
  sleep 3
done
wait
echo "worm campaign done $(date)"
