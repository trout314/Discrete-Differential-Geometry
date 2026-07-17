#!/usr/bin/env python3
"""Churn through a persistent, appendable seed-generation queue.

Reads a wishlist file (default data/seed_queue.txt) of desired seed families and,
each pass, produces the next cell whose certified seeds don't yet exist in seeds/:
found a starting triangulation -> grow_seed to (N, objective) -> equilibrium_vdv
--produce (multi-chain R-hat gate) -> copy the 32 certified replicas into seeds/.

Completion is derived purely from seeds/ (the exact library filename encoding), so
the queue file is never mutated -- it stays a clean human-edited to-do list you can
append to at any time; a --watch runner picks up appends on its next pass.

    python scripts/seed_queue_runner.py --watch           # daemon: churn forever
    python scripts/seed_queue_runner.py --once             # single pass then exit
    python scripts/seed_queue_runner.py --once --dry-run   # show the plan, build nothing
"""
import argparse
import fcntl
import glob
import os
import re
import shutil
import sys
import time
from types import SimpleNamespace

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_ROOT, "python"))
sys.path.insert(0, os.path.join(_ROOT, "tools"))
sys.path.insert(0, os.path.join(_ROOT, "scripts"))

from seed_utils import encode_float, build_seed_filename, load_seed_metadata  # noqa: E402
import grid_sweep as G                                                         # noqa: E402
from extend_library import grow, LADDER                                        # noqa: E402

DEFAULT_BRACKET = ["vdv", "edge_deg", "num_facets"]


def parse_queue(path):
    """Return an ordered list of cell dicts from the wishlist file.
    Each family line expands to one cell per N; cells sorted by (prio, N)."""
    cells = []
    with open(path) as f:
        for lineno, raw in enumerate(f, 1):
            line = raw.split("#", 1)[0].strip()
            if not line:
                continue
            kv = {}
            for tok in line.split():
                if "=" not in tok:
                    raise SystemExit(f"queue line {lineno}: bad token {tok!r}")
                k, v = tok.split("=", 1)
                kv[k] = v
            edge = float(kv.get("edge", 5.10430))
            k = float(kv.get("k", 2))
            bon = float(kv.get("beta_over_n", 0))
            hon = float(kv.get("hdv_over_n", 0))
            topology = kv.get("topology", "S3")
            # Fixed-target quadratics (VDQ/EDQ), specified as PER-TET scaled
            # couplings (same scale as beta_over_n/hdv_over_n; raw per-element
            # coefs are lambda/(6/edge-1) resp. lambda*edge/6).
            vdq = float(kv.get("vdq", 0))
            edq = float(kv.get("edq", 0))
            prio = int(kv.get("prio", 100))
            if "N" not in kv:
                raise SystemExit(f"queue line {lineno}: missing N=")
            for tokn in kv["N"].split(","):
                n = int(float(tokn))
                cells.append(dict(edge=edge, k=k, bon=bon, hon=hon,
                                  vdq=vdq, edq=edq, n=n, topology=topology,
                                  prio=prio, line=lineno))
    cells.sort(key=lambda c: (c["prio"], c["n"], c["edge"], c["bon"], c["hon"],
                              c["vdq"], c["edq"], c["topology"]))
    return cells


def vdq_raw(cell):
    """Raw per-vertex VDQ coupling c_v from the per-tet scaled queue value."""
    return cell["vdq"] / (6.0 / cell["edge"] - 1.0) if cell["vdq"] else 0.0


def edq_raw(cell):
    """Raw per-edge EDQ coupling c_e from the per-tet scaled queue value."""
    return cell["edq"] * cell["edge"] / 6.0 if cell["edq"] else 0.0


def params_at(cell, n):
    """A SamplerParams-like namespace for cell at facet target n (raw coefs)."""
    return SimpleNamespace(
        num_facets_target=n, num_facets_coef=0.1,
        hinge_degree_target=cell["edge"], num_hinges_coef=cell["k"],
        codim3_degree_variance_coef=cell["bon"] * n,
        hinge_degree_variance_coef=cell["hon"] * n,
        codim3_degree_target_coef=vdq_raw(cell),
        hinge_degree_target_coef=edq_raw(cell))


def seed_path(seeds_abs, cell, n, idx=0):
    return os.path.join(seeds_abs,
                        build_seed_filename(cell["topology"], params_at(cell, n), idx))


# Root triangulations for topologies with no library presence yet. T3 root is a
# validated TCP-crystal quotient (scripts/tcp_reference.py); it must be SMALLER
# than the smallest queued N: growing UP through the volume pin melts the
# crystal (volume-adding moves are pin-favored), whereas a root above target is
# pin-frozen (no volume-removing moves exist in a perfect {5,6} crystal).
TOPOLOGY_ROOTS = {
    "T3": os.path.join(_ROOT, "standard_triangulations", "T3_A15_m2_N368.mfd"),
}


def find_source(seeds_abs, cell):
    """A starting s000 triangulation to grow from: prefer this family's own largest
    tier below N (a stepping stone), else the nearest existing SAME-TOPOLOGY family
    (by edge then beta/N) at the largest N <= target, else the topology's root
    triangulation -- grow_seed relaxes to the objective."""
    n = cell["n"]
    # 1. same family, largest smaller tier
    best_n = -1; best = None
    for tok, nn in LADDER:
        if nn >= n:
            continue
        p = seed_path(seeds_abs, cell, nn)
        if os.path.exists(p) and nn > best_n:
            best, best_n = p, nn
    if best:
        return best, f"stepping-stone N={best_n}"
    # 2. nearest existing same-topology family at largest N <= target
    topo = cell["topology"]
    cand = None; cand_key = None
    for p in glob.glob(os.path.join(seeds_abs, f"{topo}_N*_1e-1_*_s000.mfd")):
        m = re.match(rf"{topo}_N([0-9e]+)_", os.path.basename(p))
        nn = next((v for t, v in LADDER if t == m.group(1)), None) if m else None
        if nn is None or nn > n:
            continue
        try:
            md = load_seed_metadata(p)
            e = float(md.get("hinge_degree_target", 5.1043))
            # Effective reduced vertex coupling: variance beta/N plus the VDQ
            # per-tet lambda (same scale by the exact coupling map), so old-mode
            # families rank as founders for VDQ cells at matched physics.
            b = float(md.get("codim3_degree_variance_coef", 0.0)) / max(
                float(md.get("num_facets_target", nn)), 1)
            b += float(md.get("codim3_degree_target_coef", 0.0)) * (6.0 / e - 1.0)
        except Exception:
            continue
        # rank: closest N (prefer largest<=n), then nearest edge, then nearest
        # effective reduced coupling
        key = (-nn, abs(e - cell["edge"]), abs(b - (cell["bon"] + cell["vdq"])))
        if cand_key is None or key < cand_key:
            cand, cand_key = p, key
    if cand:
        return cand, f"founded from {os.path.basename(cand)}"
    root = TOPOLOGY_ROOTS.get(topo)
    if root and os.path.exists(root):
        return root, f"topology root {os.path.basename(root)}"
    return None, "no source"


def log(msg):
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {msg}", flush=True)


def record(status_csv, cell, n, verdict, detail):
    new = not os.path.exists(status_csv)
    with open(status_csv, "a") as f:
        if new:
            f.write("time,edge,k,beta_over_n,hdv_over_n,N,verdict,detail\n")
        f.write(f"{time.strftime('%Y-%m-%dT%H:%M:%S')},{cell['edge']},{cell['k']},"
                f"{cell['bon']},{cell['hon']},{n},{verdict},\"{detail}\"\n")


def do_cell(args, seeds_abs, out_root, cell):
    """Grow + produce one cell. Returns verdict string."""
    n = cell["n"]
    edge, k = cell["edge"], cell["k"]
    beta, hdv = cell["bon"] * n, cell["hon"] * n
    cv, ce = vdq_raw(cell), edq_raw(cell)
    topo = cell["topology"]
    stem = build_seed_filename(topo, params_at(cell, n)).replace(
        f"{topo}_N{encode_float(n)}_", "").replace(".mfd", "")
    if topo != "S3":
        stem = f"{topo}_{stem}"
    tag = f"N={n} {stem}"
    src, how = find_source(seeds_abs, cell)
    if src is None:
        log(f"  {tag}: NO SOURCE — skipping"); return "NO-SOURCE"
    log(f"  {tag}: {how}")
    if args.dry_run:
        log(f"  {tag}: [dry-run] would grow -> produce (beta={beta:g}, hdv={hdv:g}, "
            f"vdq_raw={cv:g}, edq_raw={ce:g})")
        return "DRY"
    obj = dict(edge=edge, k=k, nfc=0.1)
    grown = os.path.join(args.scratch_dir, f"grow_N{encode_float(n)}_{stem}.mfd")
    if grow(_ROOT, src, n, obj, beta, hdv, grown, args.grow_min_free_gb,
            vdq_coef=cv, edq_coef=ce, topology=topo) is None:
        log(f"  {tag}: GROW-FAIL (transient? retried next pass)"); return "GROW-FAIL"
    out_dir = os.path.join(out_root, f"{stem}_N{encode_float(n)}")
    # Fresh staging every attempt: equilibrium_vdv --produce RESUMES chains whose
    # --save-config already exists, so a reused out_dir with prior-session staging
    # would re-gate stale (often frozen/short) data and instant-FAIL. Clear it so
    # each attempt runs clean. (Sacrifices mid-produce resume of a killed runner --
    # acceptable: the cell is the unit of work and re-doing one is cheap.)
    if os.path.isdir(out_dir):
        shutil.rmtree(out_dir, ignore_errors=True)
    verdict, detail = G.run_cell(
        _ROOT, grown, n, edge, beta, bracket=DEFAULT_BRACKET,
        replicas=args.replicas, burnin=args.burnin, nsamp=args.n_samples,
        thin=args.thin, dry_run=False, seeds_dir=args.seeds_dir, out_dir=out_dir,
        num_hinges_coef=k, hdv_coef=hdv, vdq_coef=cv, edq_coef=ce,
        max_workers=args.max_workers, max_memory_gb=args.max_memory_gb,
        topology=topo)
    if verdict == "FAIL" and G.retry_worthwhile(detail):
        log(f"  {tag}: FAIL ({detail}) — one longer retry")
        verdict, detail = G.run_cell(
            _ROOT, grown, n, edge, beta, bracket=DEFAULT_BRACKET,
            replicas=args.replicas, burnin=args.retry_burnin,
            nsamp=args.retry_n_samples, thin=args.thin, dry_run=False,
            seeds_dir=args.seeds_dir, out_dir=out_dir + "_retry",
            num_hinges_coef=k, hdv_coef=hdv, vdq_coef=cv, edq_coef=ce,
            max_workers=args.max_workers, max_memory_gb=args.max_memory_gb,
            topology=topo)
    log(f"  {tag}: {verdict} {detail}")
    return verdict


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--queue", default=os.path.join(_ROOT, "data", "seed_queue.txt"))
    p.add_argument("--status", default=os.path.join(_ROOT, "data", "seed_queue_status.csv"))
    p.add_argument("--out-dir", default=os.path.join(_ROOT, "data", "seed_queue"))
    p.add_argument("--scratch-dir",
                   default=os.path.join(_ROOT, "data", "seed_queue", "_grown"))
    p.add_argument("--seeds-dir", default="seeds")
    p.add_argument("--max-workers", type=int, default=8)
    p.add_argument("--max-memory-gb", type=float, default=18.0)
    p.add_argument("--replicas", type=int, default=32)
    p.add_argument("--burnin", type=int, default=6000)
    p.add_argument("--n-samples", type=int, default=2000)
    p.add_argument("--thin", type=int, default=5)
    p.add_argument("--retry-burnin", type=int, default=10000)
    p.add_argument("--retry-n-samples", type=int, default=3000)
    p.add_argument("--grow-min-free-gb", type=float, default=6.0)
    p.add_argument("--watch", action="store_true",
                   help="Loop forever, re-reading the queue each pass (daemon).")
    p.add_argument("--watch-interval", type=int, default=300,
                   help="Seconds to sleep when the queue is fully satisfied.")
    p.add_argument("--once", action="store_true", help="Single pass then exit.")
    p.add_argument("--dry-run", action="store_true",
                   help="Report the plan; build nothing.")
    args = p.parse_args()

    seeds_abs = args.seeds_dir if os.path.isabs(args.seeds_dir) \
        else os.path.join(_ROOT, args.seeds_dir)
    os.makedirs(args.out_dir, exist_ok=True)
    os.makedirs(args.scratch_dir, exist_ok=True)

    # Single-runner lock (skip for dry-run so it can inspect anytime).
    lock_fh = None
    if not args.dry_run:
        lock_fh = open(os.path.join(_ROOT, "data", ".seed_queue.lock"), "w")
        try:
            fcntl.flock(lock_fh, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except OSError:
            sys.exit("another seed_queue_runner holds the lock; exiting.")

    attempted_nonpass = set()   # cells tried this session that didn't PASS (avoid spin)

    def cell_key(c):
        return (c["edge"], c["k"], c["bon"], c["hon"],
                c["vdq"], c["edq"], c["n"], c["topology"])

    while True:
        cells = parse_queue(args.queue)
        pending = [c for c in cells
                   if not os.path.exists(seed_path(seeds_abs, c, c["n"]))
                   and cell_key(c) not in attempted_nonpass]
        done = sum(1 for c in cells if os.path.exists(seed_path(seeds_abs, c, c["n"])))
        log(f"queue: {len(cells)} cells, {done} already in seeds/, "
            f"{len(pending)} pending, {len(attempted_nonpass)} deferred this session")
        if args.dry_run:
            for c in pending:
                do_cell(args, seeds_abs, args.out_dir, c)
            return
        if not pending:
            if args.once:
                log("queue satisfied; exiting (--once)."); break
            time.sleep(args.watch_interval); continue
        c = pending[0]
        verdict = do_cell(args, seeds_abs, args.out_dir, c)
        qinfo = (f"vdq={c['vdq']:g} edq={c['edq']:g}"
                 if (c["vdq"] or c["edq"]) else "")
        record(args.status, c, c["n"], verdict, qinfo)
        if verdict not in ("PASS", "SKIP-EXISTS"):
            # GROW-FAIL is transient (memory) but we still defer for the session to
            # avoid spinning; a runner restart re-attempts everything deferred.
            attempted_nonpass.add(cell_key(c))
        if args.once and not any(
                not os.path.exists(seed_path(seeds_abs, c2, c2["n"]))
                and cell_key(c2) not in attempted_nonpass
                for c2 in cells):
            log("queue satisfied; exiting (--once)."); break


if __name__ == "__main__":
    main()
