#!/bin/bash
# Overnight large-tier reference campaign: R glass (completes R large), then
# A15 / sigma / C15 large tiers, run SEQUENTIALLY (each family drains before
# the next) so the heavy V=10-17k chains never pile up. Launch detached:
#   setsid nohup bash scripts/overnight_large_queue.sh &
cd /home/aaron/Desktop/Discrete-Differential-Geometry

launch() {   # struct  seed  fleets...
  local struct=$1 seed=$2; shift 2
  echo "=== $(date '+%H:%M') launching $struct large: $* ==="
  python3 scripts/reference_campaign.py --struct "$struct" --tier large \
    --fleets "$@" --replicas 8 --base-seed "$seed" \
    --mem-floor-gb 12 --stagger 25
}

drain() {    # out-csv path substring identifying this family's large chains
  sleep 45
  while ps -eo args | grep '[d]ope_hold' | grep -q "reference/large/$1"; do
    sleep 90
  done
  echo "=== $(date '+%H:%M') $1 large drained ==="
}

echo "### overnight large queue started $(date) ###"

# wait for the small-tier stragglers to finish first
while ps -eo args | grep '[d]ope_hold' | grep -q 'reference/small'; do sleep 30; done
echo "small tier clear; starting large tiers"

launch r     20000 glass ; drain "r_glass_large"
launch a15   21000 crystal ownpin glass ; drain "a15_"
launch sigma 22000 crystal ownpin glass ; drain "sigma_"
launch c15   23000 crystal ownpin glass flat_vac ; drain "c15_"

echo "### overnight large queue COMPLETE $(date) ###"
