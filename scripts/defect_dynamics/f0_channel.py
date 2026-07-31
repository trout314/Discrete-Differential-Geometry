"""Reversible MH channel for the f0 sector: CBMC vertex removal/insertion.

The bilocal design doc (2.4 scheme B, 2.5) applied to the frozen-f0
benchmark. Both directions are STOCHASTIC constructors over one reshaper
alphabet on star(v) -- link flips (both 3D flavors), spoke deletions,
attaches (reverse of deletion), and STOP -- with ONE shared stepwise
rule w = exp(-beta*dS - eta*dZ - gamma*dPhi), with Phi the concave
spoke potential sum sqrt(d(v,u)-3). Detailed balance holds PER PATH
(super-detailed balance): the reverse leg is priced by RETRACING the
reversed op sequence under the same rule, so every forward path has
computable, nonzero reverse support (the 0/24 deterministic-mirror
verdict forces stochastic constructors). Because the rule is shared,
the per-step normalizers W(s) cancel exactly between a path and its
reverse (same intermediate states), and the -eta*dZ term is a state
potential (telescopes), so acceptance carries only boundary terms --
this is what keeps log(alpha) O(1) instead of O(-path length).

  removal:   seed v ~ U(vertices) -> reshape until STOP; abort unless
             the walk stopped at Z=4 -> offline 4->1 (Manifold-level
             dup; capi forbids 1<->4 through the sampler).
  insertion: seed tet ~ U(tets), label ~ U(free labels < LMAX) ->
             offline 1->4 -> reshape until STOP (any Z).

  log alpha = -(S_y - S_x) + log q(y->x) - log q(x->y)
  q(x->y) = 1/2 * p_seed * prod_k w(op_k)/W_k   (STOP factors included)

Every candidate op is priced exactly by do-undo through the sampler.
The label pool is the fixed window {0..LMAX-1} minus used labels, so a
removed label is always re-proposable and the vacancy placement entropy
(f0, f3, pool-size factors) enters the acceptance ratio explicitly.

Lab driver (__main__): cycles of thermal sweeps + channel attempts on
the down-quenched m^2-only C15 gas; the sampled f0 trajectory IS the
free-energy comparison between the gap-open (f0=1536) and gap-closed
(f0=1522, defect-financed) sectors. Prototype channel: flight recording
is a local .jsonl here (channel moves happen outside the D event log).

Usage:
  python scripts/defect_dynamics/f0_channel.py START.mfd OUTBASE \
      [ETARGET [CIMP [F3T]]]
Env: CYCLES RELAX ATTEMPTS BETA ETA GAMMA WSTOP MAXSTEPS SEED LMAX
     VALIDATE
"""
import json
import math
import os
import random
import sys
import time
from itertools import combinations

_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, os.path.join(_ROOT, "python"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import discrete_differential_geometry as ddg
import worm_deg4_slide as W

ETARGET = float(sys.argv[3]) if len(sys.argv) > 3 else 5.1067907
CIMP = float(sys.argv[4]) if len(sys.argv) > 4 else 0.7
F3T = int(sys.argv[5]) if len(sys.argv) > 5 else 8704
BETA = float(os.environ.get("BETA", "0.3"))
ETA = float(os.environ.get("ETA", "1.0"))
GAMMA = float(os.environ.get("GAMMA", "1.5"))
WSTOP = float(os.environ.get("WSTOP", "0.05"))
MAXSTEPS = int(os.environ.get("MAXSTEPS", "60"))
LMAX = int(os.environ.get("LMAX", "4096"))
VALIDATE = bool(int(os.environ.get("VALIDATE", "0")))


def dup_of(x):
    """Manifold has .copy(), ManifoldView has .dup() -- unify."""
    return x.dup() if hasattr(x, "dup") else x.copy()


def _e(a, b):
    return (a, b) if a < b else (b, a)


def fresh(m):
    p = ddg.SamplerParams(
        num_facets_target=F3T, num_facets_coef=0.1,
        hinge_degree_target=ETARGET, num_hinges_coef=1.0,
        hinge_degree_variance_coef=0.0, codim3_degree_variance_coef=0.0,
        hinge_degree_target_coef=0.0, codim3_degree_target_coef=0.0)
    s = ddg.ManifoldSampler(m, p)
    s.set_n6_potential(0.0, CIMP, tilt=[0.0] * 5)
    return s, W.Live(s.manifold, driver=s.do_bistellar_move)


def facet_set(m):
    return {frozenset(int(x) for x in t) for t in np.asarray(m.facets())}


# ---------------------------------------------------------------- alphabet
def zdeg(L, v):
    return sum(1 for e in L.edeg if v in e)


def enum_ops(L, v):
    """Reshaper alphabet at star(v). Ops:
      ('flip', ab, xy, flavor)  -- link-edge flip, flavor '23'/'32'
      ('delete', u, face)       -- 3->2 on spoke (v,u); face = the link
                                   face left behind (reverse-op datum)
      ('attach', face, u)       -- 2->3 across link face; u = outside
                                   apex (determined by ambient state)
    Legality here is the cheap ambient screen; the sampler is the final
    arbiter (pricing failures drop the op)."""
    faces = [tuple(sorted(set(t) - {v})) for t in L.v2t[v]]
    wing = {}
    for f in faces:
        for pr in combinations(f, 2):
            wing.setdefault(_e(*pr), []).append(
                next(iter(set(f) - set(pr))))
    link_edges = set(wing)
    ops = []
    for ed, ws in wing.items():
        if len(ws) != 2 or _e(*ws) in link_edges:
            continue
        xy = _e(*ws)
        if xy not in L.edeg:                      # chord absent: 2->3 ok
            ops.append(("flip", ed, xy, "23"))
        if L.edeg.get(ed) == 3:                   # outside tet: 3->2 ok
            ops.append(("flip", ed, xy, "32"))
    for u in sorted({x for f in faces for x in f}):
        if L.edeg.get(_e(v, u)) == 3:
            tets = [t for t in L.v2t[v] if u in t]
            lk = sorted({x for t in tets for x in t} - {v, u})
            if len(tets) == 3 and len(lk) == 3:
                ops.append(("delete", u, tuple(lk)))
    for f in set(faces):
        ts = [t for t in L.v2t[f[0]] if set(f) <= set(t)]
        if len(ts) != 2:
            continue
        tout = ts[0] if v not in ts[0] else ts[1]
        if v in tout:
            continue
        u = next(iter(set(tout) - set(f)))
        if _e(v, u) not in L.edeg:
            ops.append(("attach", f, u))
    return ops


def op_moves(v, op):
    """3D (center, cocenter) of a reshaper op."""
    if op[0] == "flip":
        _, ab, xy, flavor = op
        if flavor == "23":
            return tuple(sorted((v,) + ab)), xy
        return ab, tuple(sorted((v,) + xy))
    if op[0] == "delete":
        return _e(v, op[1]), op[2]
    return op[1], _e(v, op[2])                    # attach


def rev_op(op):
    if op[0] == "flip":
        _, ab, xy, flavor = op
        return ("flip", xy, ab, "32" if flavor == "23" else "23")
    if op[0] == "delete":
        return ("attach", op[2], op[1])
    return ("delete", op[2], op[1])


def dz_of(op):
    return {"delete": -1, "attach": +1}.get(op[0], 0)


def phi(L, v):
    """Concave spoke potential sum sqrt(d(v,u)-3): rewards concentrating
    degree so that spokes reach the deg-3 deletion threshold (a 4->3
    drop gains 1.0, the compensating 5->6 rise costs only ~0.32). A pure
    state function -- telescopes along any path, so shaping with it
    leaves only boundary terms in the acceptance."""
    return sum(math.sqrt(d - 3) for e, d in L.edeg.items()
               if v in e and d > 3)


def price(L, s, v, op):
    """Exact dS (and potential change) of op by do-undo through the
    sampler; None if illegal."""
    cen, coc = op_moves(v, op)
    S0 = s.current_objective
    try:
        L.do(cen, coc)
    except Exception:
        return None
    dS = s.current_objective - S0
    dphi = phi(L, v)
    L.do(coc, cen)
    return dS, dphi


def weights(L, s, v):
    """Priced alphabet + STOP weight under the ONE shared rule
    w = exp(-beta*dS - eta*dZ). A single rule makes the per-step
    normalizers W(s) cancel exactly between a forward path and its
    retraced reverse (same states), so acceptance carries only boundary
    terms; the -eta*dZ shaping is a state POTENTIAL (telescopes along
    any path), funneling collapse walks downhill in Z without biasing
    the balance. STOP is always available; a removal walk that stops at
    Z != 4 is simply aborted (no move -- harmless for balance)."""
    cands = []
    phi0 = phi(L, v)
    for op in enum_ops(L, v):
        pr = price(L, s, v, op)
        if pr is None:
            continue
        dS, phi1 = pr
        w = math.exp(max(-60.0, min(60.0,
                                    -BETA * dS - ETA * dz_of(op)
                                    - GAMMA * (phi1 - phi0))))
        cands.append((op, dS, w))
    return cands, WSTOP


# ------------------------------------------------------------ constructors
def sample_reshape(L, s, v, rng):
    """Run one stochastic reshape; returns (path, logq) or None on
    step-cap abort. logq includes all w/W factors and the final STOP."""
    logq, path = 0.0, []
    for _ in range(MAXSTEPS + 1):
        cands, wstop = weights(L, s, v)
        Wtot = wstop + sum(w for _, _, w in cands)
        if Wtot <= 0.0:
            return None
        r = rng.random() * Wtot
        if r < wstop:
            return path, logq + math.log(wstop / Wtot)
        r -= wstop
        for op, dS, w in cands:
            if r < w:
                cen, coc = op_moves(v, op)
                L.do(cen, coc)
                path.append(op)
                logq += math.log(w / Wtot)
                break
            r -= w
        if len(path) > MAXSTEPS:
            return None
    return None


def retrace_reshape(L, s, v, ops):
    """Apply `ops` in order, accumulating the log-probability the shared
    constructor would have assigned (incl. final STOP). Returns logq or
    None if any op is unsupported."""
    logq = 0.0
    for op in ops:
        cands, wstop = weights(L, s, v)
        Wtot = wstop + sum(w for _, _, w in cands)
        hit = next((w for o, _, w in cands if o == op), None)
        if hit is None or Wtot <= 0.0:
            return None
        cen, coc = op_moves(v, op)
        L.do(cen, coc)
        logq += math.log(hit / Wtot)
    cands, wstop = weights(L, s, v)
    Wtot = wstop + sum(w for _, _, w in cands)
    if wstop <= 0.0:
        return None
    return logq + math.log(wstop / Wtot)


def label_pool_size(L):
    return LMAX - len({x for e in L.edeg for x in e})


# --------------------------------------------------------------- proposals
def propose(m, rng, force_dir=None):
    """One channel attempt on manifold m. Returns a dict with the
    proposed manifold (or None on abort) + diagnostics. force_dir is a
    testing hook only -- production draws the direction coin."""
    work = dup_of(m)               # never mutate the caller's state
    s, L = fresh(work)
    S0 = s.current_objective
    verts = sorted({x for e in L.edeg for x in e})
    f0, f3 = len(verts), int(m.f_vector[3])
    direction = force_dir or ("rm" if rng.random() < 0.5 else "ins")
    info = {"dir": direction}

    if direction == "rm":
        v = verts[rng.randrange(f0)]
        info["seed"] = v
        logq_f = -math.log(f0)
        res = sample_reshape(L, s, v, rng)
        if res is None:
            return {**info, "abort": "reshape"}
        path, lq = res
        logq_f += lq
        if zdeg(L, v) != 4:
            return {**info, "abort": "stopz", "steps": len(path)}
        nb = sorted({x for t in L.v2t[v] for x in t} - {v})
        try:
            L.do([v], nb)              # 4->1 through the sampler
        except Exception:
            return {**info, "abort": "4to1"}
        dS = s.current_objective - S0
        # reverse: insertion on y, seed tet nb, label v, retraced path
        logq_r = (-math.log(int(s.manifold.f_vector[3]))
                  - math.log(label_pool_size(L)))
        y = dup_of(s.manifold)         # proposal snapshot
        L.do(sorted(nb), [v])          # 1->4 back, in place
        lqr = retrace_reshape(L, s, v,
                              [rev_op(op) for op in reversed(path)])
        if lqr is None:
            return {**info, "abort": "retrace"}
        logq_r += lqr
        if VALIDATE:
            assert facet_set(s.manifold) == facet_set(m), \
                "reverse retrace does not reconstruct x"
        return {**info, "m_new": y, "dS": dS, "logq_f": logq_f,
                "logq_r": logq_r, "steps": len(path)}

    # insertion
    all_tets = sorted({t for ts in L.v2t.values() for t in ts})
    t_seed = all_tets[rng.randrange(len(all_tets))]
    pool = [x for x in range(LMAX) if x not in set(verts)]
    v_new = pool[rng.randrange(len(pool))]
    info["seed"] = list(t_seed)
    info["label"] = v_new
    logq_f = -math.log(len(all_tets)) - math.log(len(pool))
    try:
        L.do(sorted(t_seed), [v_new])  # 1->4 through the sampler
    except Exception:
        return {**info, "abort": "1to4"}
    res = sample_reshape(L, s, v_new, rng)
    if res is None:
        return {**info, "abort": "reshape"}
    path, lq = res
    logq_f += lq
    dS = s.current_objective - S0
    y = dup_of(s.manifold)           # proposed state (reshaped, v_new in)
    # reverse: removal on y seeded at v_new, retraced reversed path
    logq_r = -math.log(f0 + 1)
    lqr = retrace_reshape(L, s, v_new,
                          [rev_op(op) for op in reversed(path)])
    if lqr is None:
        return {**info, "abort": "retrace"}
    logq_r += lqr
    if VALIDATE:
        nb = sorted({x for t in L.v2t[v_new] for x in t} - {v_new})
        L.do([v_new], nb)            # 4->1 completes the reverse
        assert facet_set(s.manifold) == facet_set(m), \
            "reverse retrace does not reconstruct x"
    return {**info, "m_new": y, "dS": dS, "logq_f": logq_f,
            "logq_r": logq_r, "steps": len(path)}


def mh_step(m, rng):
    """One MH attempt; returns (m_or_new, accepted, info)."""
    info = propose(m, rng)
    if "m_new" not in info:
        return m, False, info
    la = -info["dS"] + info["logq_r"] - info["logq_f"]
    info["log_alpha"] = la
    if math.log(max(rng.random(), 1e-300)) < la:
        return info.pop("m_new"), True, info
    info.pop("m_new")
    return m, False, info


# --------------------------------------------------------------------- lab
def main():
    start, outbase = sys.argv[1], sys.argv[2]
    cycles = int(os.environ.get("CYCLES", "300"))
    relax = float(os.environ.get("RELAX", "2"))
    attempts = int(os.environ.get("ATTEMPTS", "5"))
    seed = int(os.environ.get("SEED", "42"))
    rng = random.Random(seed)
    ddg.set_random_seed(seed)
    m = ddg.Manifold.load(start, 3)
    log = open(outbase + ".chan.jsonl", "w")
    log.write(json.dumps({
        "start": start, "etarget": ETARGET, "cimp": CIMP, "f3t": F3T,
        "beta": BETA, "eta": ETA, "wstop": WSTOP, "maxsteps": MAXSTEPS,
        "lmax": LMAX, "cycles": cycles, "relax": relax,
        "attempts": attempts, "seed": seed}) + "\n")
    t0, n_acc = time.time(), {"rm": 0, "ins": 0}
    for cyc in range(cycles):
        s, _ = fresh(m)
        if relax > 0:
            s.run(sweeps=relax)
        m = s.manifold
        rows = []
        for _ in range(attempts):
            m, acc, info = mh_step(m, rng)
            if acc:
                n_acc[info["dir"]] += 1
            rows.append({k: v for k, v in info.items()
                         if k in ("dir", "seed", "label", "abort",
                                  "dS", "log_alpha", "steps")}
                        | {"acc": acc})
        s, _ = fresh(m)
        fv = [int(x) for x in m.f_vector]
        gap = fv[1] - 6.0 * fv[3] / ETARGET
        log.write(json.dumps({
            "cyc": cyc, "f": fv, "gap": round(gap, 3),
            "S": round(s.current_objective, 3), "att": rows}) + "\n")
        log.flush()
        if cyc % 10 == 0 or cyc == cycles - 1:
            print(f"cyc {cyc:4d}  f0={fv[0]} f3={fv[3]} gap={gap:+6.2f} "
                  f"S={s.current_objective:8.2f} acc rm/ins = "
                  f"{n_acc['rm']}/{n_acc['ins']} "
                  f"({time.time() - t0:.0f}s)", flush=True)
    m.save(outbase + ".final.mfd")
    log.close()
    print(f"done: {outbase}.final.mfd  +  {outbase}.chan.jsonl", flush=True)


if __name__ == "__main__":
    main()
