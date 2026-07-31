"""f0 worm: extended-ensemble vertex removal/insertion (design doc 3.2).

Scheme C. Alongside closed states T (weight e^-S), open states (T, v)
with a marked head vertex and weight zeta * e^{-S + U(T,v)}. U and zeta
shape only the auxiliary open sector -- the closed-sector conditional is
e^-S exactly, so measurements (gated on closed) are unbiased for ANY U.
Choose U ~ the measured collapse-cost staircase and the open walk sees
S - U ~ flat: the +20 barrier is crossed diffusively by many O(1)
local accepts instead of one e^-20 proposal.

Moves (fixed mixture weights; invalid draws auto-reject = null move):
  open-flag    T -> (T,v):  alpha = zeta e^U f0 * BCF/AOF
  close-flag   (T,v) -> T:  alpha = e^-U/(zeta f0) * AOF/BCF
  reshape      (T,v) -> (T',v): any 2<->3 move whose support meets
               {v} u lk(v)  (the union-of-stars region; membership is
               EXACTLY move-symmetric -- a move not containing v
               changes no edge at v), uniform over the optimistic
               candidate set A: alpha = e^{-(dS-dU)} |A|/|A'|
  open-insert  T -> (T+,v): tet ~ U(f3), label ~ U(pool), 1->4:
               alpha = zeta e^{-dS14+U} f3*|pool| * BC4/AOI
  close-41     (T,v) at Z=4 -> T-: 4->1:
               alpha = e^{-dS41-U}/zeta * AOI/(BC4 f3(T-)|pool(T-)|)

Episodes are contiguous (ordinary sweeps run only while closed);
close-flag is the escape hatch (accepts freely where U ~ 0), so no
abort exists. Episode debris is physical and balanced. The 1<->4
sector crossings run through the sampler (capi label-pool support).

v1 umbrella: U = ETA*(Z0 - Z(v)) + GAMMA*(PHI0 - Phi(v)), Phi =
sum sqrt(d(v,u)-3); Z0/PHI0 calibrated from the start state so
U(typical star) ~ 0. Staircase table = v2.

Usage: python scripts/defect_dynamics/f0_worm.py START.mfd OUTBASE \
           [ETARGET [CIMP [F3T]]]
Env: ETA GAMMA ZETA AOF BR BCF BC4 CYCLES RELAX EPISODES MAXSTEP SEED
     LMAX  (AOI = 1-AOF)
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
import discrete_differential_geometry as ddg
import worm_deg4_slide as W

ETARGET = float(sys.argv[3]) if len(sys.argv) > 3 else 5.1067907
CIMP = float(sys.argv[4]) if len(sys.argv) > 4 else 0.7
F3T = int(sys.argv[5]) if len(sys.argv) > 5 else 8704
ETA = float(os.environ.get("ETA", "1.2"))
GAMMA = float(os.environ.get("GAMMA", "1.0"))
ZETA = float(os.environ.get("ZETA", "1e-3"))
AOF = float(os.environ.get("AOF", "0.5"))          # open-flag | closed
AOI = 1.0 - AOF                                    # open-insert | closed
BR = float(os.environ.get("BR", "0.90"))           # reshape | open
BCF = float(os.environ.get("BCF", "0.05"))         # close-flag | open
BC4 = 1.0 - BR - BCF                               # close-41 | open
MAXSTEP = int(os.environ.get("MAXSTEP", "3000"))
LMAX = int(os.environ.get("LMAX", "4096"))


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


def nb_of(L, v):
    return {x for t in L.v2t[v] for x in t} - {v}


def phi_of(L, v):
    out = 0.0
    for u in nb_of(L, v):
        d = L.edeg[_e(v, u)]
        if d > 3:
            out += math.sqrt(d - 3)
    return out


# Z0/PHI0 are fixed at startup (calibration constants, NOT state).
Z0 = 12.0
PHI0 = 18.0


def U_of(L, v):
    return ETA * (Z0 - len(nb_of(L, v))) + GAMMA * (PHI0 - phi_of(L, v))


def region_candidates(L, v):
    """All optimistically-valid 2<->3 moves whose support meets
    {v} u lk(v). Support-membership is move-symmetric (doc 3.2 lemma);
    optimistic validity (cheap static checks, D-core rejection of a
    drawn move = auto-reject) is applied identically on both sides."""
    head = {v} | nb_of(L, v)
    tets = set()
    for u in head:
        tets |= L.v2t[u]
    out = []
    seen_f, seen_e = set(), set()
    for t in tets:
        for f in combinations(t, 3):
            if f in seen_f:
                continue
            seen_f.add(f)
            ts = [tt for tt in L.v2t[f[0]] if
                  f[1] in tt and f[2] in tt]
            if len(ts) != 2:
                continue
            ap = tuple(sorted((set(ts[0]) | set(ts[1])) - set(f)))
            if len(ap) != 2 or ap in L.edeg:
                continue
            if not (set(f) | set(ap)) & head:
                continue
            out.append((f, ap))                       # 2->3
        for e in combinations(t, 2):
            if e in seen_e:
                continue
            seen_e.add(e)
            if L.edeg.get(e) != 3:
                continue
            t3 = [tt for tt in L.v2t[e[0]] if e[1] in tt]
            if len(t3) != 3:
                continue
            lk = sorted({x for tt in t3 for x in tt} - set(e))
            if len(lk) != 3:
                continue
            if not (set(e) | set(lk)) & head:
                continue
            out.append((e, tuple(lk)))                # 3->2
    return out


class Worm:
    """Single-head f0 worm on one live sampler."""

    def __init__(self, s, L, rng):
        self.s, self.L, self.rng = s, L, rng
        self.acc = {k: [0, 0] for k in
                    ("of", "oi", "re", "cf", "c4")}   # [tries, accepts]
        self.forced = 0

    def _hit(self, kind, ok):
        self.acc[kind][0] += 1
        self.acc[kind][1] += ok
        return ok

    def _mh(self, kind, log_alpha):
        return self._hit(kind, math.log(max(self.rng.random(), 1e-300))
                         < log_alpha)

    def verts(self):
        return {x for e in self.L.edeg for x in e}

    def try_open(self):
        """One opening attempt from closed. Returns head vertex or None."""
        s, L, rng = self.s, self.L, self.rng
        vs = sorted(self.verts())
        f0 = len(vs)
        if rng.random() < AOF:
            v = vs[rng.randrange(f0)]
            la = (math.log(ZETA) + U_of(L, v) + math.log(f0)
                  + math.log(BCF / AOF))
            return v if self._mh("of", la) else None
        tets = sorted({t for u in vs for t in L.v2t[u]})
        f3 = len(tets)
        tet = tets[rng.randrange(f3)]
        vset = set(vs)
        pool = [x for x in range(LMAX) if x not in vset]
        vn = pool[rng.randrange(len(pool))]
        S0 = s.current_objective
        try:
            L.do(sorted(tet), [vn])
        except Exception:
            self._hit("oi", 0)
            return None
        dS = s.current_objective - S0
        la = (math.log(ZETA) - dS + U_of(L, vn) + math.log(f3)
              + math.log(len(pool)) + math.log(BC4 / AOI))
        if self._mh("oi", la):
            return vn
        L.do([vn], sorted(tet))
        return None

    def reshape(self, v):
        s, L, rng = self.s, self.L, self.rng
        A = region_candidates(L, v)
        if not A:
            return
        cen, coc = A[rng.randrange(len(A))]
        S0, U0 = s.current_objective, U_of(L, v)
        try:
            L.do(cen, coc)
        except Exception:
            self._hit("re", 0)
            return
        dS = s.current_objective - S0
        dU = U_of(L, v) - U0
        n2 = len(region_candidates(L, v))
        if n2 == 0:                    # reverse unproposable: reject
            self._hit("re", 0)
            L.do(coc, cen)
            return
        la = -dS + dU + math.log(len(A)) - math.log(n2)
        if not self._mh("re", la):
            L.do(coc, cen)

    def try_close(self, v):
        """One closing attempt. Returns 'cf'/'c4' on success, None else."""
        s, L, rng = self.s, self.L, self.rng
        r = rng.random()
        if r < BR:
            self.reshape(v)
            return None
        f0 = len(self.verts())
        if r < BR + BCF:                               # close-flag
            la = (-U_of(L, v) - math.log(ZETA) - math.log(f0)
                  + math.log(AOF / BCF))
            return "cf" if self._mh("cf", la) else None
        nb = sorted(nb_of(L, v))                       # close-41
        if len(nb) != 4:
            self._hit("c4", 0)
            return None
        S0, U0 = s.current_objective, U_of(L, v)
        try:
            L.do([v], nb)
        except Exception:
            self._hit("c4", 0)
            return None
        dS = s.current_objective - S0
        f3m = int(s.manifold.f_vector[3])
        poolm = LMAX - (f0 - 1)
        la = (-dS - U0 - math.log(ZETA) + math.log(AOI / BC4)
              - math.log(f3m) - math.log(poolm))
        if self._mh("c4", la):
            return "c4"
        L.do(sorted(nb), [v])
        return None

    def episode(self):
        """One full attempt: open, walk, close. Returns a record."""
        t0 = time.time()
        S0 = self.s.current_objective
        v = self.try_open()
        if v is None:
            return {"opened": 0}
        steps, closed, umax = 0, None, U_of(self.L, v)
        while closed is None:
            steps += 1
            if steps > MAXSTEP:                        # emergency only:
                self.forced += 1                       # breaks exactness,
                closed = "forced"                      # must stay 0
                break
            closed = self.try_close(v)
            umax = max(umax, U_of(self.L, v))
        return {"opened": 1, "head": int(v), "steps": steps,
                "closed": closed, "dS": round(
                    self.s.current_objective - S0, 3),
                "umax": round(umax, 2), "t": round(time.time() - t0, 2)}


def main():
    global Z0, PHI0
    start, outbase = sys.argv[1], sys.argv[2]
    cycles = int(os.environ.get("CYCLES", "200"))
    relax = float(os.environ.get("RELAX", "2"))
    episodes = int(os.environ.get("EPISODES", "3"))
    seed = int(os.environ.get("SEED", "42"))
    rng = random.Random(seed)
    ddg.set_random_seed(seed)
    m = ddg.Manifold.load(start, 3)
    s, L = fresh(m)
    # calibrate U's zero point from the start state (fixed constants)
    vs = sorted({x for e in L.edeg for x in e})
    samp = vs[:: max(1, len(vs) // 200)]
    Z0 = sum(len(nb_of(L, u)) for u in samp) / len(samp)
    PHI0 = sum(phi_of(L, u) for u in samp) / len(samp)
    log = open(outbase + ".chan.jsonl", "w")
    log.write(json.dumps({
        "start": start, "etarget": ETARGET, "cimp": CIMP, "f3t": F3T,
        "eta": ETA, "gamma": GAMMA, "zeta": ZETA, "z0": round(Z0, 3),
        "phi0": round(PHI0, 3), "aof": AOF, "br": BR, "bcf": BCF,
        "bc4": BC4, "cycles": cycles, "relax": relax,
        "episodes": episodes, "maxstep": MAXSTEP, "seed": seed}) + "\n")
    print(f"U calib: Z0={Z0:.2f} PHI0={PHI0:.2f}  "
          f"U(Z=4 star)~{ETA * (Z0 - 4) + GAMMA * PHI0:.1f}", flush=True)
    worm = Worm(s, L, rng)
    t0 = time.time()
    for cyc in range(cycles):
        if relax > 0:
            s.run(sweeps=relax)
            worm.L = L = W.Live(s.manifold, driver=s.do_bistellar_move)
        eps = [worm.episode() for _ in range(episodes)]
        fv = [int(x) for x in s.manifold.f_vector]
        gap = fv[1] - 6.0 * fv[3] / ETARGET
        log.write(json.dumps({
            "cyc": cyc, "f": fv, "gap": round(gap, 3),
            "S": round(s.current_objective, 3), "ep": eps}) + "\n")
        log.flush()
        if cyc % 10 == 0 or cyc == cycles - 1:
            a = worm.acc
            print(f"cyc {cyc:4d} f0={fv[0]} gap={gap:+6.2f} "
                  f"S={s.current_objective:8.2f} | "
                  + " ".join(f"{k}:{x[1]}/{x[0]}" for k, x in a.items())
                  + f" forced={worm.forced} ({time.time() - t0:.0f}s)",
                  flush=True)
    s.manifold.save(outbase + ".final.mfd")
    log.close()
    print(f"done: {outbase}.final.mfd (+ .chan.jsonl) "
          f"forced={worm.forced}", flush=True)


if __name__ == "__main__":
    main()
