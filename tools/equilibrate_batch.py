#!/usr/bin/env python3
"""Batch runner for equilibrate_seed.py.

Reads a JSON parameter table and runs equilibrate_seed.py for each combination,
optionally in parallel.  By default, concurrency is automatically limited so that
the estimated total memory of running jobs stays within --max-memory-gb (24 GB).

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
import threading
from concurrent.futures import Future, ThreadPoolExecutor, as_completed
from pathlib import Path

from seed_utils import get_free_memory_gb

SCRIPT = str(Path(__file__).parent / "equilibrate_seed.py")

# Must match equilibrate_seed.py
EXIT_LOW_MEMORY = 42

# Linear model fitted from measurements:
#   1K facets -> 38 MB, 10K -> 42 MB, 100K -> 112 MB, 1M -> 2191 MB
# ~2.15 MB per 1K facets + ~36 MB base overhead.
_MB_PER_1K_FACETS = 2.15
_BASE_MB = 36.0


def estimate_memory_mb(num_facets_target: int) -> float:
    """Estimate peak RSS in MB for a job with the given facet target."""
    return _BASE_MB + _MB_PER_1K_FACETS * (num_facets_target / 1000)


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
        "min_free_memory_gb": "--min-free-memory-gb",
    }

    for key, flag in arg_map.items():
        if key in merged:
            cmd.extend([flag, str(merged[key])])

    if merged.get("skip_existing"):
        cmd.append("--skip-existing")

    result = subprocess.run(cmd)
    return result.returncode


def run_memory_aware(combos, shared, max_memory_mb, min_free_gb):
    """Schedule jobs largest-first, limiting concurrency by estimated memory."""
    # Sort largest jobs first so they don't pile up at the end
    indexed = list(enumerate(combos))
    indexed.sort(key=lambda ic: ic[1].get("num_facets_target", 0), reverse=True)

    lock = threading.Lock()
    memory_in_use = 0.0
    memory_available = threading.Event()
    memory_available.set()
    low_memory_abort = False

    total = len(combos)
    completed = 0
    failed = 0

    def wrapped_run(combo, est_mb):
        nonlocal memory_in_use
        try:
            return run_one(combo, shared), est_mb
        finally:
            with lock:
                memory_in_use -= est_mb
            memory_available.set()

    # Use a large pool — concurrency is controlled by memory, not pool size
    max_pool = max(1, int(max_memory_mb / (_BASE_MB + _MB_PER_1K_FACETS)))
    with ThreadPoolExecutor(max_workers=min(max_pool, total)) as pool:
        futures: dict[Future, tuple[int, dict]] = {}
        submit_idx = 0

        def submit_ready():
            """Submit as many jobs as fit in the memory budget."""
            nonlocal submit_idx, memory_in_use
            if low_memory_abort:
                return
            while submit_idx < len(indexed):
                # Check actual free memory before submitting
                free_gb = get_free_memory_gb()
                if free_gb < min_free_gb:
                    return
                orig_idx, combo = indexed[submit_idx]
                est = estimate_memory_mb(combo.get("num_facets_target", 0))
                with lock:
                    if memory_in_use + est > max_memory_mb and memory_in_use > 0:
                        # No room — wait for something to finish
                        return
                    memory_in_use += est
                fut = pool.submit(wrapped_run, combo, est)
                futures[fut] = (orig_idx, combo)
                submit_idx += 1

        submit_ready()

        while futures:
            # Wait for any one future to complete
            done_iter = as_completed(futures)
            future = next(done_iter)

            rc, _ = future.result()
            orig_idx, combo = futures.pop(future)
            completed += 1
            if rc == EXIT_LOW_MEMORY:
                low_memory_abort = True
                failed += 1
                print(
                    f"LOW MEMORY: child exited with code {EXIT_LOW_MEMORY}, "
                    f"stopping batch (waiting for {len(futures)} running job(s) to finish).",
                    file=sys.stderr, flush=True,
                )
            elif rc != 0:
                failed += 1
                print(f"FAILED (rc={rc}): {combo}", file=sys.stderr, flush=True)
            if completed % 100 == 0 or completed == total or low_memory_abort:
                with lock:
                    mem = memory_in_use
                free = get_free_memory_gb()
                print(
                    f"Progress: {completed}/{total} done, {failed} failed, "
                    f"~{mem / 1024:.1f} GB estimated in use, {free:.1f} GB free",
                    flush=True,
                )

            if not low_memory_abort:
                # Try to submit more now that memory freed up
                submit_ready()

    return completed, failed


def main():
    parser = argparse.ArgumentParser(description="Batch equilibration runner.")
    parser.add_argument("--params-file", required=True, help="JSON file with parameter combos")
    parser.add_argument("--output-dir", default="seeds")
    parser.add_argument("--equilibration-sweeps", type=int, default=500)
    parser.add_argument("--eq-sweeps-per-step", type=int, default=5)
    parser.add_argument(
        "--workers", type=int, default=None,
        help="Max parallel workers (default: auto, limited by --max-memory-gb)",
    )
    parser.add_argument(
        "--max-memory-gb", type=float, default=24.0,
        help="Max total estimated memory for running jobs in GB (default: 24)",
    )
    parser.add_argument("--skip-existing", action="store_true", help="Skip seeds whose output file already exists")
    parser.add_argument(
        "--min-free-memory-gb", type=float, default=4.0,
        help="Abort batch if system free memory drops below this threshold in GB (default: 4)."
    )
    args = parser.parse_args()

    with open(args.params_file) as f:
        combos = json.load(f)

    shared = {
        "output_dir": args.output_dir,
        "equilibration_sweeps": args.equilibration_sweeps,
        "eq_sweeps_per_step": args.eq_sweeps_per_step,
        "skip_existing": args.skip_existing,
        "min_free_memory_gb": args.min_free_memory_gb,
    }

    total = len(combos)
    max_memory_mb = args.max_memory_gb * 1024

    # If --workers is set, use fixed concurrency (legacy behavior)
    if args.workers is not None:
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
    else:
        # Memory-aware scheduling (default)
        # Show the plan
        size_counts: dict[int, int] = {}
        for c in combos:
            t = c.get("num_facets_target", 0)
            size_counts[t] = size_counts.get(t, 0) + 1
        print(f"Running {total} jobs, max memory budget: {args.max_memory_gb:.0f} GB", flush=True)
        for size in sorted(size_counts, reverse=True):
            est = estimate_memory_mb(size)
            max_conc = max(1, int(max_memory_mb / est))
            print(
                f"  {size:>12,} facets: {size_counts[size]:>4} jobs, "
                f"~{est / 1024:.1f} GB each, max {max_conc} concurrent",
                flush=True,
            )
        completed, failed = run_memory_aware(combos, shared, max_memory_mb, args.min_free_memory_gb)

    print(f"\nBatch complete: {completed} jobs, {failed} failed.", flush=True)


if __name__ == "__main__":
    main()
