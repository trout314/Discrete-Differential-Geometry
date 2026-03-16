#!/usr/bin/env python3
"""Batch runner for equilibrate_seed.py.

Reads a JSON parameter table and runs equilibrate_seed.py for each combination,
optionally in parallel.

Example JSON input (pass via --params-file):
[
    {
        "topology": "S3",
        "num_facets_target": 1000,
        "num_facets_coef": 0.1,
        "hinge_degree_target": 5.1,
        "num_hinges_coef": 0.1,
        "hinge_degree_variance_coef": 0.0,
        "codim3_degree_variance_coef": 0.1,
        "seed_index": 0,
        "growth_step_size": 50
    }
]
"""

import argparse
import json
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

SCRIPT = str(Path(__file__).parent / "equilibrate_seed.py")


def run_one(combo: dict, shared_args: dict) -> int:
    """Run equilibrate_seed.py for one parameter combination.

    Subprocess inherits stdout/stderr directly so output streams live.
    Returns the process return code.
    """
    cmd = [sys.executable, SCRIPT, "--batch-mode"]

    # Merge shared args with per-combo overrides
    merged = {**shared_args, **combo}

    arg_map = {
        "topology": "--topology",
        "dimension": "--dimension",
        "num_facets_target": "--num-facets-target",
        "num_facets_coef": "--num-facets-coef",
        "hinge_degree_target": "--hinge-degree-target",
        "num_hinges_coef": "--num-hinges-coef",
        "hinge_degree_variance_coef": "--hinge-degree-variance-coef",
        "codim3_degree_variance_coef": "--codim3-degree-variance-coef",
        "equilibration_sweeps": "--equilibration-sweeps",
        "growth_step_size": "--growth-step-size",
        "eq_sweeps_per_step": "--eq-sweeps-per-step",
        "output_dir": "--output-dir",
        "seed_file": "--seed-file",
        "seed_index": "--seed-index",
    }

    for key, flag in arg_map.items():
        if key in merged:
            cmd.extend([flag, str(merged[key])])

    result = subprocess.run(cmd)
    return result.returncode


def main():
    parser = argparse.ArgumentParser(description="Batch equilibration runner.")
    parser.add_argument("--params-file", required=True, help="JSON file with parameter combos")
    parser.add_argument("--output-dir", default="seeds")
    parser.add_argument("--equilibration-sweeps", type=int, default=500)
    parser.add_argument("--eq-sweeps-per-step", type=int, default=5)
    parser.add_argument("--workers", type=int, default=1, help="Parallel workers (default: 1)")
    args = parser.parse_args()

    with open(args.params_file) as f:
        combos = json.load(f)

    shared = {
        "output_dir": args.output_dir,
        "equilibration_sweeps": args.equilibration_sweeps,
        "eq_sweeps_per_step": args.eq_sweeps_per_step,
    }

    total = len(combos)
    print(f"Running {total} jobs with {args.workers} worker(s)", flush=True)

    failed = 0
    completed = 0

    if args.workers <= 1:
        for combo in combos:
            rc = run_one(combo, shared)
            completed += 1
            if rc != 0:
                failed += 1
                print(f"FAILED (rc={rc}): {combo}", file=sys.stderr, flush=True)
    else:
        with ThreadPoolExecutor(max_workers=args.workers) as pool:
            futures = {pool.submit(run_one, combo, shared): combo for combo in combos}
            for future in as_completed(futures):
                completed += 1
                rc = future.result()
                if rc != 0:
                    failed += 1
                    combo = futures[future]
                    print(f"FAILED (rc={rc}): {combo}", file=sys.stderr, flush=True)
                if completed % 100 == 0 or completed == total:
                    print(f"Progress: {completed}/{total} done, {failed} failed", flush=True)

    print(f"\nBatch complete: {completed} jobs, {failed} failed.", flush=True)


if __name__ == "__main__":
    main()
