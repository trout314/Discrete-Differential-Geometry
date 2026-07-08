#!/usr/bin/env python3
"""Re-burn the whole seed library under the current objective build.

For each family (a set of ``_s{NNN}`` replicas sharing one parameter set):

  1. CALIBRATE the per-family burn length N_burn with the ensemble-drift gate
     (see reburn_family.calibrate).  N_burn values are cached to a JSON file so
     a resumed run does not recalibrate.

  2. PRODUCE: re-burn every member for N_burn sweeps via memory-aware parallel
     subprocesses (reburn_seed.py), writing new seeds into --output-dir.

The old library in --seed-dir is never modified.  Use --skip-existing to resume.
"""

import argparse
import glob
import json
import os
import shutil
import subprocess
import sys
import threading
import time
from collections import defaultdict
from concurrent.futures import Future, ThreadPoolExecutor, as_completed

sys.path.insert(0, os.path.dirname(__file__))
from seed_utils import get_free_memory_gb, get_total_memory_gb, load_seed_metadata
from reburn_family import calibrate, family_stem, params_from_metadata

WORKER = os.path.join(os.path.dirname(__file__), "reburn_seed.py")
EXIT_LOW_MEMORY = 42

# Per-worker resident memory scales ~linearly with facet count.  Calibrated from
# the peak RSS of production workers in the July 2026 OOM incident: a chain of
# ~316k facets was resident at ~30 GB, i.e. ~95 MB per 1000 facets.  The old
# model (2.15 MB per 1000 facets, borrowed from equilibrate_batch.py) undercounted
# by ~45x, so the scheduler believed dozens of 30 GB workers fit in 24 GB and
# oversubscribed RAM until the machine thrashed to a hard freeze.  100 MB/1k is a
# deliberately round, slightly conservative fit.
_MB_PER_1K_FACETS = 100.0
_BASE_MB = 200.0


def estimate_memory_mb(num_facets_target: int) -> float:
    return _BASE_MB + _MB_PER_1K_FACETS * (num_facets_target / 1000)


def discover_families(seed_dir: str, filt: list[str] | None) -> dict[str, list[str]]:
    fams: dict[str, list[str]] = defaultdict(list)
    for path in glob.glob(os.path.join(seed_dir, "*.mfd")):
        fams[family_stem(path)].append(path)
    for stem in fams:
        fams[stem].sort()
    if filt:
        fams = {k: v for k, v in fams.items() if any(f in k for f in filt)}
    return dict(sorted(fams.items()))


def systemd_run_available() -> bool:
    """True if `systemd-run --user` can create memory-capped transient scopes."""
    if not shutil.which("systemd-run"):
        return False
    try:
        # `--scope true` is a no-op scope; succeeds only if the user manager and
        # a delegated memory controller are actually reachable.
        r = subprocess.run(
            ["systemd-run", "--user", "--scope", "--quiet",
             "-p", "MemoryMax=64M", "--", "true"],
            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=15)
        return r.returncode == 0
    except (OSError, subprocess.SubprocessError):
        return False


def run_worker(seed_file, params, n_burn, output_dir, min_free_gb,
               worker_mem_max_gb=None) -> int:
    cmd = [
        sys.executable, WORKER,
        "--seed-file", seed_file,
        "--output-dir", output_dir,
        "--n-burn", str(n_burn),
        "--num-facets-target", str(params.num_facets_target),
        "--num-facets-coef", str(params.num_facets_coef),
        "--hinge-degree-target", str(params.hinge_degree_target),
        "--num-hinges-coef", str(params.num_hinges_coef),
        "--hinge-degree-variance-coef", str(params.hinge_degree_variance_coef),
        "--codim3-degree-variance-coef", str(params.codim3_degree_variance_coef),
        "--min-free-memory-gb", str(min_free_gb),
        "--skip-existing",
    ]
    if worker_mem_max_gb:
        # Run the worker in its own memory-capped cgroup.  If it exceeds the cap
        # the kernel OOM-kills *only this scope* (MemorySwapMax=0 => it can't
        # thrash swap first), so one oversized chain fails cleanly with a nonzero
        # return code instead of dragging the whole machine into a freeze.
        mb = int(worker_mem_max_gb * 1024)
        cmd = [
            "systemd-run", "--user", "--scope", "--quiet",
            "-p", f"MemoryMax={mb}M", "-p", "MemorySwapMax=0",
            "--", *cmd,
        ]
    return subprocess.run(cmd).returncode


def schedule_memory_aware(jobs, output_dir, max_memory_mb, min_free_gb, max_workers,
                          worker_mem_max_gb=None):
    """jobs: list of (seed_file, params, n_burn, est_mb). Largest-first."""
    jobs = sorted(jobs, key=lambda j: j[3], reverse=True)
    lock = threading.Lock()
    mem_in_use = 0.0
    total = len(jobs)
    completed = failed = 0
    abort = False

    def wrapped(job):
        nonlocal mem_in_use
        sf, params, n_burn, est = job
        try:
            return run_worker(sf, params, n_burn, output_dir, min_free_gb,
                              worker_mem_max_gb)
        finally:
            with lock:
                mem_in_use -= est

    with ThreadPoolExecutor(max_workers=min(max_workers, total)) as pool:
        futures: dict[Future, tuple] = {}
        idx = 0

        def submit_ready():
            nonlocal idx, mem_in_use
            if abort:
                return
            while idx < len(jobs):
                if len(futures) >= max_workers:
                    return
                if get_free_memory_gb() < min_free_gb:
                    return
                est = jobs[idx][3]
                with lock:
                    if mem_in_use + est > max_memory_mb and mem_in_use > 0:
                        return
                    mem_in_use += est
                futures[pool.submit(wrapped, jobs[idx])] = jobs[idx]
                idx += 1

        submit_ready()
        while futures:
            fut = next(as_completed(futures))
            job = futures.pop(fut)
            rc = fut.result()
            completed += 1
            if rc == EXIT_LOW_MEMORY:
                abort = True
                failed += 1
                print(f"LOW MEMORY: worker exited {EXIT_LOW_MEMORY}; draining "
                      f"{len(futures)} running job(s).", file=sys.stderr, flush=True)
            elif rc != 0:
                failed += 1
                print(f"FAILED (rc={rc}): {os.path.basename(job[0])}", file=sys.stderr, flush=True)
            if completed % 50 == 0 or completed == total or abort:
                with lock:
                    mem = mem_in_use
                print(f"  produce: {completed}/{total} done, {failed} failed, "
                      f"~{mem/1024:.1f} GB in use, {get_free_memory_gb():.1f} GB free", flush=True)
            if not abort:
                submit_ready()
    return completed, failed


def load_cache(path):
    if path and os.path.exists(path):
        with open(path) as f:
            return json.load(f)
    return {}


def save_cache(path, cache):
    if not path:
        return
    tmp = path + ".tmp"
    with open(tmp, "w") as f:
        json.dump(cache, f, indent=2, sort_keys=True)
    os.replace(tmp, path)


def main():
    p = argparse.ArgumentParser(description="Re-burn the seed library under the current build.")
    p.add_argument("--seed-dir", default="seeds")
    p.add_argument("--output-dir", default="seeds_reburned")
    p.add_argument("--families", nargs="*", default=None,
                   help="Only families whose stem contains one of these substrings.")
    p.add_argument("--nburn-cache", default="seeds_reburned/nburn_cache.json")
    p.add_argument("--recalibrate", action="store_true", help="Ignore cached N_burn.")
    p.add_argument("--calibrate-only", action="store_true")
    p.add_argument("--nburn-margin", type=float, default=1.0,
                   help="Multiply calibrated N_burn by this safety factor (default 1.0).")
    p.add_argument("--nburn-floor", type=int, default=200,
                   help="Never burn fewer than this many sweeps (default 200).")
    p.add_argument("--calib-max-facets", type=int, default=20000,
                   help="Families larger than this skip the (slow) gate and use "
                        "--large-nburn; large-N drift is negligible and N_burn is "
                        "non-increasing in N, so the small-N max is conservative there.")
    p.add_argument("--large-nburn", type=int, default=500,
                   help="Fixed burn length applied to families above --calib-max-facets.")
    # calibration gate
    p.add_argument("--calib-subset", type=int, default=16)
    p.add_argument("--chunk-sweeps", type=int, default=100)
    p.add_argument("--window", type=int, default=5)
    p.add_argument("--patience", type=int, default=2)
    p.add_argument("--t-thresh", type=float, default=1.5)
    p.add_argument("--min-sweeps", type=int, default=200)
    p.add_argument("--max-sweeps", type=int, default=3000)
    # scheduler
    p.add_argument("--max-memory-gb", type=float, default=24.0)
    p.add_argument("--min-free-memory-gb", type=float, default=4.0)
    p.add_argument("--max-workers", type=int, default=max(1, (os.cpu_count() or 4) - 2),
                   help="Cap on concurrent worker processes (default: ncores-2).")
    p.add_argument("--worker-mem-max-gb", type=float, default=None,
                   help="Hard per-worker memory cap enforced via a systemd-run "
                        "cgroup (MemoryMax, swap disabled): a runaway worker is "
                        "OOM-killed inside its own scope instead of freezing the "
                        "machine.  Default: --max-memory-gb when systemd-run "
                        "--user is available.")
    p.add_argument("--no-worker-cgroup", action="store_true",
                   help="Disable the per-worker cgroup cap even if available.")
    args = p.parse_args()

    # Resolve the hard per-worker memory cap (fix: cgroup isolation so one
    # oversized chain can't take down the whole box).
    worker_mem_max_gb = None
    if not args.no_worker_cgroup:
        worker_mem_max_gb = args.worker_mem_max_gb or args.max_memory_gb
        if not systemd_run_available():
            print("WARNING: systemd-run --user unavailable; running workers WITHOUT "
                  "a hard memory cap. A single oversized chain can still freeze the "
                  "machine. Consider adding swap or running on a bigger host.",
                  file=sys.stderr, flush=True)
            worker_mem_max_gb = None
        else:
            print(f"Per-worker memory cap: {worker_mem_max_gb:.0f} GB "
                  f"(systemd-run cgroup, swap disabled).", flush=True)

    fams = discover_families(args.seed_dir, args.families)
    if not fams:
        print("No families matched.", file=sys.stderr)
        sys.exit(1)

    total_members = sum(len(v) for v in fams.values())
    print(f"Discovered {len(fams)} families, {total_members} seeds total.", flush=True)

    cache = {} if args.recalibrate else load_cache(args.nburn_cache)

    # ---- Phase A: calibrate every family (cached) ----
    dim_cache: dict[str, tuple] = {}
    for stem, members in fams.items():
        meta = load_seed_metadata(members[0])
        params, topology, dimension = params_from_metadata(meta)
        dim_cache[stem] = (params, topology, dimension)

        if stem in cache and not args.recalibrate:
            print(f"[calib] {stem}: cached N_burn={cache[stem]['n_burn']}", flush=True)
            continue

        # Large families: skip the slow gate, use a conservative fixed burn.
        if params.num_facets_target > args.calib_max_facets:
            n_burn = max(args.nburn_floor, args.large_nburn)
            cache[stem] = {"n_burn": n_burn, "num_facets_target": params.num_facets_target,
                           "calib": "fixed (large-N)"}
            os.makedirs(os.path.dirname(args.nburn_cache) or ".", exist_ok=True)
            save_cache(args.nburn_cache, cache)
            print(f"[calib] {stem}: N={params.num_facets_target} > "
                  f"{args.calib_max_facets}, fixed N_burn={n_burn}", flush=True)
            continue

        # shrink subset if memory would be tight.  Calibration loads `subset`
        # chains simultaneously *in this process* and then samples them, so it
        # is its own OOM vector.  Budget only half of max-memory for the resident
        # chains, leaving headroom for sampling transients and the OS/desktop.
        est = estimate_memory_mb(params.num_facets_target)
        calib_budget_mb = 0.5 * args.max_memory_gb * 1024
        subset = max(2, min(args.calib_subset,
                            int(calib_budget_mb / max(est, 1.0)),
                            len(members)))

        print(f"[calib] {stem}: N={params.num_facets_target} subset={subset}", flush=True)
        t0 = time.monotonic()
        n_burn, traj = calibrate(
            members, params, dimension,
            subset=subset, chunk_sweeps=args.chunk_sweeps,
            window=args.window, patience=args.patience, t_thresh=args.t_thresh,
            min_sweeps=args.min_sweeps, max_sweeps=args.max_sweeps,
            log=lambda m: None,  # quiet; keep the summary line below
        )
        n_burn = max(args.nburn_floor, int(round(n_burn * args.nburn_margin)))
        cache[stem] = {
            "n_burn": n_burn,
            "num_facets_target": params.num_facets_target,
            "edge_saved": traj[0][1], "edge_final": traj[-1][1],
            "calib_sweeps": traj[-1][0], "calib_subset": subset,
        }
        os.makedirs(os.path.dirname(args.nburn_cache) or ".", exist_ok=True)
        save_cache(args.nburn_cache, cache)
        print(f"[calib] {stem}: N_burn={n_burn} "
              f"(edge {traj[0][1]:.5f}->{traj[-1][1]:.5f}, {time.monotonic()-t0:.0f}s)", flush=True)

    if args.calibrate_only:
        print("DONE (calibrate-only).", flush=True)
        return

    # ---- Phase B: produce ----
    total_ram_gb = get_total_memory_gb()
    jobs = []
    oversized = []
    for stem, members in fams.items():
        params, _, _ = dim_cache[stem]
        n_burn = cache[stem]["n_burn"]
        est = estimate_memory_mb(params.num_facets_target)
        # A family whose single-worker estimate exceeds the per-worker cap (or,
        # absent a cap, physical RAM) cannot be re-burned here: every worker will
        # be OOM-killed (or thrash).  Flag it loudly rather than freeze silently.
        cap_gb = worker_mem_max_gb or total_ram_gb
        if est / 1024 > cap_gb:
            oversized.append((stem, params.num_facets_target, est / 1024))
        for sf in members:
            jobs.append((sf, params, n_burn, est))

    if oversized:
        print(f"WARNING: {len(oversized)} family(ies) exceed the per-worker memory "
              f"budget ({(worker_mem_max_gb or total_ram_gb):.0f} GB) and will fail:",
              file=sys.stderr, flush=True)
        for stem, n, gb in sorted(oversized, key=lambda x: -x[2]):
            print(f"  {stem}: N={n} ~{gb:.0f} GB/worker", file=sys.stderr, flush=True)
        print("  Re-run these on a larger host, or add swap and use "
              "--no-worker-cgroup to let them spill.", file=sys.stderr, flush=True)

    print(f"Producing {len(jobs)} seeds -> {args.output_dir} "
          f"(scheduler budget {args.max_memory_gb:.0f} GB, {total_ram_gb:.0f} GB RAM)",
          flush=True)
    t0 = time.monotonic()
    completed, failed = schedule_memory_aware(
        jobs, args.output_dir, args.max_memory_gb * 1024, args.min_free_memory_gb,
        args.max_workers, worker_mem_max_gb)
    print(f"\nDONE: {completed} seeds, {failed} failed, {time.monotonic()-t0:.0f}s.", flush=True)


if __name__ == "__main__":
    main()
