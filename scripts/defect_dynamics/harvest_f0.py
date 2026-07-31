"""Force-mode f0 harvest: planner-built vertex removals, no Hastings.

Repeated rounds of: cost-optimal (link, collar) plan (link_planner) ->
lift through the sampler -> 4->1 through the sampler (the capi now
supports targeted 1<->4 with label-pool bookkeeping) -> short thermal
relax of the fast sector. Benchmark (bilocal design doc 3): on the
down-quenched m^2-only C15 gas, close the +10 pin gap, f0 1536 -> 1522,
each removal net-favorable with the net price decaying to ~0 as the gap
closes.

Usage:
  python scripts/defect_dynamics/harvest_f0.py START.mfd OUT.mfd \
      [ETARGET [CIMP [F3T]]]
Env: ROUNDS RELAX NODECAP TRYMAX SEED (see defaults below).

RESULT 2026-07-30 (quench_down5q_wOFF): 14/14 removals, zero lift
failures, planned==executed dS exact, gap +10.12 -> +0.49 in 20 s.
"""
import os, sys, time

_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, os.path.join(_ROOT, "python"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import discrete_differential_geometry as ddg
import worm_deg4_slide as W
import link_planner as LP

START = sys.argv[1]
OUT = sys.argv[2]
ETARGET = float(sys.argv[3]) if len(sys.argv) > 3 else 5.1067907
CIMP = float(sys.argv[4]) if len(sys.argv) > 4 else 0.7
F3T = int(sys.argv[5]) if len(sys.argv) > 5 else 8704
ROUNDS = int(os.environ.get("ROUNDS", "18"))
RELAX = float(os.environ.get("RELAX", "10"))     # sweeps between removals
NODECAP = int(os.environ.get("NODECAP", "200000"))
TRYMAX = int(os.environ.get("TRYMAX", "8"))
SEED = int(os.environ.get("SEED", "1234"))


def fresh(m):
    p = ddg.SamplerParams(
        num_facets_target=F3T, num_facets_coef=0.1,
        hinge_degree_target=ETARGET, num_hinges_coef=1.0,
        hinge_degree_variance_coef=0.0, codim3_degree_variance_coef=0.0,
        hinge_degree_target_coef=0.0, codim3_degree_target_coef=0.0)
    s = ddg.ManifoldSampler(m, p)
    s.set_n6_potential(0.0, CIMP, tilt=[0.0] * 5)
    return s, W.Live(s.manifold, driver=s.do_bistellar_move)


def op_to_move(L, v, op):
    """Translate one 2D plan op into the 3D (center, cocenter)."""
    if op[0] == "delete":
        u = op[1]
        tets = [t for t in L.v2t[v] if u in t]
        lk = sorted({x for t in tets for x in t} - {v, u})
        return (v, u) if v < u else (u, v), tuple(lk)
    _, ab, xy, flavor = op
    if flavor == "23":
        return tuple(sorted((v,) + ab)), xy
    return ab, tuple(sorted((v,) + xy))


def main():
    m = ddg.Manifold.load(START, 3)
    s, L = fresh(m)
    S00 = s.current_objective
    fv = [int(x) for x in m.f_vector]
    print(f"start f={fv} S00={S00:.3f} "
          f"gap={fv[1] - 6.0 * fv[3] / ETARGET:+.2f}", flush=True)

    removed, t_start = 0, time.time()
    for rnd in range(ROUNDS):
        s, L = fresh(m)
        S0r = s.current_objective
        fv = [int(x) for x in m.f_vector]
        f1, f3 = fv[1], fv[3]
        gap = f1 - 6.0 * f3 / ETARGET
        if gap <= 0.5:
            print(f"\nGAP CLOSED at round {rnd}: f={fv} gap={gap:+.2f}",
                  flush=True)
            break
        Z = lambda vv: sum(1 for e in L.edeg if vv in e)
        targets = sorted({x for e in L.edeg for x in e},
                         key=lambda u: (Z(u), u))
        ok = False
        for tv in targets[:TRYMAX]:
            t0 = time.time()
            ops, pds = LP.plan_collapse(L, tv, CIMP, 1.0, ETARGET, f3, f1,
                                        F3T, nodecap=NODECAP,
                                        optimize="cost")
            if ops is None:
                print(f" r{rnd} v={tv} Z={Z(tv)}: no plan", flush=True)
                continue
            applied, barrier = [], 0.0
            try:
                for op in ops:
                    cen, coc = op_to_move(L, tv, op)
                    L.do(cen, coc)
                    applied.append((cen, coc))
                    barrier = max(barrier, s.current_objective - S0r)
            except Exception as ex:
                for cen, coc in reversed(applied):
                    L.do(coc, cen)
                print(f" r{rnd} v={tv}: lift failed ({ex}) -- rolled back "
                      f"(drift {s.current_objective - S0r:+.1e})",
                      flush=True)
                continue
            dS_exec = s.current_objective - S0r
            nb = sorted({x for t in L.v2t[tv] for x in t} - {tv})
            if len(nb) != 4:
                for cen, coc in reversed(applied):
                    L.do(coc, cen)
                print(f" r{rnd} v={tv}: Z_end={len(nb)} != 4 -- rolled "
                      f"back", flush=True)
                continue
            L.do([tv], nb)                     # 4->1 through the sampler
            net = s.current_objective - S0r
            fv2 = [int(x) for x in s.manifold.f_vector]
            gap2 = fv2[1] - 6.0 * fv2[3] / ETARGET
            match = abs(dS_exec - pds) < 1e-6
            print(f" r{rnd} v={tv} REMOVED: {len(ops)}+1 moves, planned "
                  f"{pds:+.3f} exec {dS_exec:+.3f} match={match} barrier "
                  f"{barrier:+.2f} refund {net - dS_exec:+.2f} NET "
                  f"{net:+.3f} | f={fv2} gap={gap2:+.2f} "
                  f"({time.time() - t0:.1f}s)", flush=True)
            ddg.set_random_seed(SEED + rnd)
            if RELAX > 0:
                s.run(sweeps=RELAX)
            dS_rel = s.current_objective - S0r - net
            m = s.manifold
            fv3 = [int(x) for x in m.f_vector]
            print(f"   relax {RELAX:g}sw dS={dS_rel:+.3f} -> f={fv3} "
                  f"gap={fv3[1] - 6.0 * fv3[3] / ETARGET:+.2f} | "
                  f"cum S-S00 = {s.current_objective - S00:+.3f}",
                  flush=True)
            removed += 1
            ok = True
            break
        if not ok:
            print(f"round {rnd}: no removable vertex in first {TRYMAX} "
                  f"targets -- stopping", flush=True)
            break

    s, _ = fresh(m)
    fv = [int(x) for x in m.f_vector]
    m.save(OUT)
    print(f"\nDONE: {removed} removals in {time.time() - t_start:.0f}s | "
          f"f={fv} gap={fv[1] - 6.0 * fv[3] / ETARGET:+.2f} "
          f"S-S00={s.current_objective - S00:+.3f}\nsnapshot: {OUT}",
          flush=True)


if __name__ == "__main__":
    main()
