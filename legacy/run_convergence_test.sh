#!/bin/bash
# Run multiple independent MCMC chains in parallel for convergence diagnostics.
# Usage: ./tools/run_convergence_test.sh [num_chains] [params_file]
#
# Defaults: 4 chains, tools/convergence_test.params
# Each chain gets its own CSV output: data/convergence_chain_N.dat

set -euo pipefail

NUM_CHAINS=${1:-4}
PARAMS_FILE=${2:-tools/convergence_test.params}
SAMPLER=./builddir-release/manifold_sampler

if [ ! -x "$SAMPLER" ]; then
    echo "Error: $SAMPLER not found. Build with:"
    echo "  meson compile -C builddir-release manifold_sampler"
    exit 1
fi

if [ ! -f "$PARAMS_FILE" ]; then
    echo "Error: $PARAMS_FILE not found."
    exit 1
fi

echo "Launching $NUM_CHAINS independent chains..."

PIDS=()
for i in $(seq 1 "$NUM_CHAINS"); do
    # Create per-chain params file with unique prefix
    CHAIN_PARAMS="/tmp/convergence_chain_${i}.params"
    sed "s|sampleFilesPrefix = .*|sampleFilesPrefix = data/convergence_chain_${i}|" \
        "$PARAMS_FILE" > "$CHAIN_PARAMS"

    echo "  Chain $i -> data/convergence_chain_${i}.dat"
    $SAMPLER -p "$CHAIN_PARAMS" &
    PIDS+=($!)
done

echo "Waiting for all chains to finish..."

FAILED=0
for i in "${!PIDS[@]}"; do
    if wait "${PIDS[$i]}"; then
        echo "  Chain $((i+1)) finished."
    else
        echo "  Chain $((i+1)) FAILED (exit code $?)."
        FAILED=1
    fi
done

# Clean up temp param files
rm -f /tmp/convergence_chain_*.params

if [ "$FAILED" -eq 1 ]; then
    echo "Some chains failed. Check output above."
    exit 1
fi

echo ""
echo "All chains complete. Run convergence analysis:"
echo "  python3 tools/convergence_analysis.py data/convergence_chain_{1..$NUM_CHAINS}.dat"
