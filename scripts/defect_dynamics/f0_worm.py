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

v3 umbrella: FROZEN spoke-multiset table harvested at startup from
planner collapse corridors (exact on-corridor; per-spoke linear
least-squares fallback off-corridor). Any constant table is a valid
umbrella -- provenance never enters the balance.

Usage: python scripts/defect_dynamics/f0_worm.py START.mfd OUTBASE \
           [ETARGET [CIMP [F3T]]]
Reshape proposals are TWO-CLASS (v2 fix): with prob PH uniform over
candidates whose support contains v (the collapse-relevant class H --
class membership is move-symmetric since the support set is invariant),
else uniform over the rest; Hastings uses the drawn class's counts on
both sides (PH cancels). Episodes that hit MAXSTEP are EXACTLY UNDONE
(full accepted-move unwind, recorded dS must be 0) -- never committed.

Env: DSIDE PHD UNSEEDS UNC WLCYC WLD0 UTAB_LOAD ZETA AOF BR BCF BC4 PH
     CYCLES RELAX EPISODES MAXSTEP SEED LMAX  (AOI = 1-AOF)

WLCYC > 0 runs a Wang-Landau calibration phase first (visited
multisets penalized by a staged decaying delta, periodic re-anchor),
then FREEZES the table (saved to OUTBASE.utab.json, reloadable via
UTAB_LOAD), reloads the start state, and runs production with full
exactness (any frozen table is a valid umbrella).
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
UNSEEDS = int(os.environ.get("UNSEEDS", "16"))
UNC = int(os.environ.get("UNC", "150000"))
WLCYC = int(os.environ.get("WLCYC", "0"))          # WL calibration cycles
WLD0 = float(os.environ.get("WLD0", "0.5"))        # initial WL delta
UTAB_LOAD = os.environ.get("UTAB_LOAD", "")        # frozen table json
ZETA = float(os.environ.get("ZETA", "0.05"))
AOF = float(os.environ.get("AOF", "0.5"))          # open-flag | closed
AOI = 1.0 - AOF                                    # open-insert | closed
BR = float(os.environ.get("BR", "0.94"))           # reshape | open
BCF = float(os.environ.get("BCF", "0.01"))         # close-flag | open
BC4 = 1.0 - BR - BCF                               # close-41 | open
PH = float(os.environ.get("PH", "0.5"))           # head-class proposal prob
MAXSTEP = int(os.environ.get("MAXSTEP", "3000"))
DSIDE = bool(int(os.environ.get("DSIDE", "0")))   # D-core episode engine
PHD = float(os.environ.get("PHD", "0.45"))        # D: head-kernel share
MU = float(os.environ.get("MU", "3.0"))           # D: seed bias e^-mu*Z
ZETA_D = float(os.environ.get("ZETA_D", "50.0"))
BCF_D = float(os.environ.get("BCF_D", "1e-4"))
BC4_D = float(os.environ.get("BC4_D", "0.1"))
OFFPEN = float(os.environ.get("OFFPEN", "-3.0"))  # off-tube U
UCAP_HI = float(os.environ.get("UCAP_HI", "45.0"))
UCAP_LO = float(os.environ.get("UCAP_LO", "-15.0"))
LMAX = int(os.environ.get("LMAX", "4096"))
RETUBE_EVERY = int(os.environ.get("RETUBE_EVERY", "10"))  # commits/rebuild


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


def spoke_ms(L, v):
    """Sorted spoke-degree multiset of the head star."""
    return tuple(sorted(L.edeg[_e(v, u)] for u in nb_of(L, v)))


# The umbrella: a FROZEN table over spoke multisets, seeded at startup
# from planner-harvested collapse corridors (each corridor state has a
# distinct multiset, so the measured staircase is represented EXACTLY
# on-corridor -- a per-spoke SUM provably cannot do this, residuals
# +-7). Off-table states fall back to a least-squares per-spoke linear
# form fitted to the same corridor data. Any constant table is a valid
# umbrella; provenance does not enter the balance.
UTAB = {}
UFB = [0.0] * 6      # (n3, n4, n5, n6, n7plus, Zdeficit) coefficients
UFBC = 0.0           # fallback constant (shifts with re-anchoring)
Z0 = 12.0


def U_fallback(ms):
    n = [0] * 5
    for d in ms:
        n[min(max(d, 3), 7) - 3] += 1
    z = Z0 - len(ms)
    return (UFB[0] * n[0] + UFB[1] * n[1] + UFB[2] * n[2]
            + UFB[3] * n[3] + UFB[4] * n[4] + UFB[5] * z + UFBC)


def U_of(L, v):
    ms = spoke_ms(L, v)
    u = UTAB.get(ms)
    return u if u is not None else U_fallback(ms)


def build_umbrella(s, L, nseeds, nodecap):
    """Harvest collapse corridors from a stratified vertex sample and
    freeze the multiset->cost table + linear fallback. Deterministic;
    exact rollback after each replay."""
    import numpy as np
    import link_planner as LP
    fv = [int(x) for x in s.manifold.f_vector]
    vs = sorted({x for e in L.edeg for x in e},
                key=lambda u: (len(nb_of(L, u)), u))
    seeds = vs[:: max(1, len(vs) // nseeds)][:nseeds]
    S0 = s.current_objective
    data = {}                      # ms -> [cum values]
    for v in seeds:
        ops, c = LP.plan_collapse(L, v, CIMP, 1.0, ETARGET, fv[3],
                                  fv[1], F3T, nodecap=nodecap,
                                  optimize="cost")
        if ops is None:
            continue
        # only PLANNED corridors enter the tube (an unplanned start
        # state would be a dead-end entrance at U = 0)
        data.setdefault(spoke_ms(L, v), []).append(0.0)
        applied = []
        try:
            for op in ops:
                if op[0] == "delete":
                    u = op[1]
                    tets = [t for t in L.v2t[v] if u in t]
                    lk = sorted({x for t in tets for x in t} - {v, u})
                    cen, coc = _e(v, u), tuple(lk)
                else:
                    _, ab, xy, flavor = op
                    cen, coc = ((tuple(sorted((v,) + ab)), xy)
                                if flavor == "23"
                                else (ab, tuple(sorted((v,) + xy))))
                L.do(cen, coc)
                applied.append((cen, coc))
                data.setdefault(spoke_ms(L, v), []).append(
                    s.current_objective - S0)
        except Exception:
            pass
        for cen, coc in reversed(applied):
            L.do(coc, cen)
        assert abs(s.current_objective - S0) < 1e-6
    # MIN over seeds: where corridors of different orbits share a
    # multiset, undercrediting is safe (umbrella too weak = slow walk)
    # while overcrediting inflates sector-crossing rates.
    tab = {ms: min(cs) for ms, cs in data.items()}
    # linear fallback fit on the corridor data
    global Z0
    Z0 = sum(len(nb_of(L, u)) for u in seeds) / len(seeds)
    rows, targ = [], []
    for ms, u in tab.items():
        n = [0] * 5
        for d in ms:
            n[min(max(d, 3), 7) - 3] += 1
        rows.append(n + [Z0 - len(ms)])
        targ.append(u)
    coef, *_ = np.linalg.lstsq(np.array(rows, float),
                               np.array(targ), rcond=None)
    return tab, list(coef)


def region_split(L, v):
    """All optimistically-valid 2<->3 moves whose support meets
    {v} u lk(v), split into (H, R): H = support contains v (collapse-
    relevant), R = rest. Support-membership AND the H/R split are
    move-symmetric (doc 3.2 lemma; the support set is invariant);
    optimistic validity (cheap static checks, D-core rejection of a
    drawn move = auto-reject) is applied identically on both sides."""
    head = {v} | nb_of(L, v)
    tets = set()
    for u in head:
        tets |= L.v2t[u]
    outH, outR = [], []
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
            sup = set(f) | set(ap)
            if not sup & head:
                continue
            (outH if v in sup else outR).append((f, ap))   # 2->3
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
            sup = set(e) | set(lk)
            if not sup & head:
                continue
            (outH if v in sup else outR).append((e, tuple(lk)))  # 3->2
    return outH, outR


class Worm:
    """Single-head f0 worm on one live sampler."""

    def __init__(self, s, L, rng):
        self.s, self.L, self.rng = s, L, rng
        self.acc = {k: [0, 0] for k in
                    ("of", "oi", "re", "cf", "c4")}   # [tries, accepts]
        self.undone = 0
        self.wl_delta = 0.0        # >0 only during WL calibration

    def _hit(self, kind, ok):
        self.acc[kind][0] += 1
        self.acc[kind][1] += ok
        return ok

    def _mh(self, kind, log_alpha):
        return self._hit(kind, math.log(max(self.rng.random(), 1e-300))
                         < log_alpha)

    def verts(self):
        return {x for e in self.L.edeg for x in e}

    def try_open(self, applied):
        """One opening attempt from closed. Returns head vertex or
        None; appends the opening move (if any) to `applied`."""
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
            applied.append((sorted(tet), [vn]))
            return vn
        L.do([vn], sorted(tet))
        return None

    def reshape(self, v):
        """One two-class reshape attempt. Returns the accepted (cen,
        coc) or None. The drawn class's counts enter the Hastings on
        both sides; PH cancels (class assignment is move-invariant)."""
        s, L, rng = self.s, self.L, self.rng
        H, R = region_split(L, v)
        pick_h = rng.random() < PH
        A = H if pick_h else R
        if not A:
            self._hit("re", 0)
            return None
        cen, coc = A[rng.randrange(len(A))]
        n1 = len(A)
        S0, U0 = s.current_objective, U_of(L, v)
        try:
            L.do(cen, coc)
        except Exception:
            self._hit("re", 0)
            return None
        dS = s.current_objective - S0
        dU = U_of(L, v) - U0
        H2, R2 = region_split(L, v)
        n2 = len(H2) if pick_h else len(R2)
        if n2 == 0:                    # reverse unproposable: reject
            self._hit("re", 0)
            L.do(coc, cen)
            return None
        la = -dS + dU + math.log(n1) - math.log(n2)
        if self._mh("re", la):
            return (cen, coc)
        L.do(coc, cen)
        return None

    def close_flag(self, v):
        L = self.L
        f0 = len(self.verts())
        la = (-U_of(L, v) - math.log(ZETA) - math.log(f0)
              + math.log(AOF / BCF))
        return self._mh("cf", la)

    def close_41(self, v):
        s, L = self.s, self.L
        nb = sorted(nb_of(L, v))
        if len(nb) != 4:
            self._hit("c4", 0)
            return False
        f0 = len(self.verts())
        S0, U0 = s.current_objective, U_of(L, v)
        try:
            L.do([v], nb)
        except Exception:
            self._hit("c4", 0)
            return False
        dS = s.current_objective - S0
        f3m = int(s.manifold.f_vector[3])
        poolm = LMAX - (f0 - 1)
        la = (-dS - U0 - math.log(ZETA) + math.log(AOI / BC4)
              - math.log(f3m) - math.log(poolm))
        if self._mh("c4", la):
            return True
        L.do(sorted(nb), [v])
        return False

    def episode(self):
        """One full attempt: open, walk, close. A capped episode is
        EXACTLY UNDONE (recorded dS must be 0) -- never committed."""
        t0 = time.time()
        S0 = self.s.current_objective
        applied = []
        v = self.try_open(applied)
        if v is None:
            return {"opened": 0}
        rng = self.rng
        steps, closed, umax = 0, None, U_of(self.L, v)
        while closed is None:
            steps += 1
            if steps > MAXSTEP:
                for cen, coc in reversed(applied):
                    self.L.do(coc, cen)
                self.undone += 1
                return {"opened": 1, "head": int(v), "steps": steps,
                        "closed": "undone", "dS": round(
                            self.s.current_objective - S0, 6),
                        "umax": round(umax, 2),
                        "t": round(time.time() - t0, 2)}
            r = rng.random()
            if r < BR:
                mv = self.reshape(v)
                if mv is not None:
                    applied.append(mv)
            elif r < BR + BCF:
                if self.close_flag(v):
                    closed = "cf"
            else:
                if self.close_41(v):
                    closed = "c4"
            if closed is None and self.wl_delta:
                # Wang-Landau: penalize the visited multiset (flatten
                # the visit histogram; table is frozen for production)
                ms = spoke_ms(self.L, v)
                UTAB[ms] = U_of(self.L, v) - self.wl_delta
            umax = max(umax, U_of(self.L, v))
        return {"opened": 1, "head": int(v), "steps": steps,
                "closed": closed, "dS": round(
                    self.s.current_objective - S0, 3),
                "umax": round(umax, 2), "t": round(time.time() - t0, 2)}


def build_orbit_tube(s, L, ntry=8, nodecap=200000):
    """Single-orbit tube umbrella: the cost-optimal collapse corridor
    of the cheapest low-Z seed, replayed exactly (planner cost model ==
    executed dS to machine precision). Multiset values from ONE
    coherent context -- no cross-orbit min-agg traps. Returns
    ({multiset: (cum dS, df1, df3)}, (f1, f3)) -- each corridor state's
    cumulative dS, its exact f-vector offset from the corridor start,
    and the build-time f reference. The D engine reprices the global-pin
    part of each value at every episode open, so ONE build stays valid
    as f drifts (retube only to refresh the local m^2 content)."""
    import link_planner as LP
    fv = [int(x) for x in s.manifold.f_vector]
    S0 = s.current_objective
    vs = sorted({x for e in L.edeg for x in e},
                key=lambda u: (len(nb_of(L, u)), u))
    best = None
    for v in vs[:ntry]:
        ops, c = LP.plan_collapse(L, v, CIMP, 1.0, ETARGET, fv[3],
                                  fv[1], F3T, nodecap=nodecap,
                                  optimize="cost")
        if ops and (best is None or c < best[2]):
            best = (v, ops, c)
    if best is None:
        return None
    v, ops, c = best
    tab = {spoke_ms(L, v): (0.0, 0, 0)}
    applied = []
    try:
        for op in ops:
            if op[0] == "delete":
                u = op[1]
                tets = [t for t in L.v2t[v] if u in t]
                lk = sorted({x for t in tets for x in t} - {v, u})
                cen, coc = _e(v, u), tuple(lk)
            else:
                _, ab, xy, flavor = op
                cen, coc = ((tuple(sorted((v,) + ab)), xy)
                            if flavor == "23"
                            else (ab, tuple(sorted((v,) + xy))))
            L.do(cen, coc)
            applied.append((cen, coc))
            fvn = s.manifold.f_vector
            tab[spoke_ms(L, v)] = (
                round(s.current_objective - S0, 6),
                int(fvn[1]) - fv[1], int(fvn[3]) - fv[3])
    except Exception:
        pass
    for cen, coc in reversed(applied):
        L.do(coc, cen)
    assert abs(s.current_objective - S0) < 1e-6
    return tab, (fv[1], fv[3])


def anchor_utab(L, sample):
    """Uniform-shift the umbrella so U(typical closed star) ~ 0. A
    global shift is a gauge choice (trades against zeta only)."""
    global UFBC
    off = sum(U_of(L, u) for u in sample) / len(sample)
    for k in list(UTAB):
        UTAB[k] -= off
    UFBC -= off
    return off


def save_utab(path):
    with open(path, "w") as f:
        json.dump({"z0": Z0, "ufb": list(UFB), "ufbc": UFBC,
                   "tab": {",".join(map(str, k)): v
                           for k, v in UTAB.items()}}, f)


def load_utab(path):
    global Z0, UFBC
    d = json.load(open(path))
    Z0 = d["z0"]
    UFB[:] = d["ufb"]
    UFBC = d["ufbc"]
    UTAB.update({tuple(int(x) for x in k.split(",")): v
                 for k, v in d["tab"].items()})


def main():
    
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
    t_cal = time.time()
    if UTAB_LOAD:
        load_utab(UTAB_LOAD)
        print(f"umbrella: loaded {len(UTAB)} multisets from "
              f"{UTAB_LOAD}", flush=True)
    else:
        tab, coef = build_umbrella(s, L, UNSEEDS, UNC)
        UTAB.update(tab)
        UFB[:] = coef
        print(f"umbrella: {len(UTAB)} corridor multisets from "
              f"{UNSEEDS} seeds ({time.time() - t_cal:.0f}s); "
              f"Z0={Z0:.2f}; U[3,3,3,3]="
              f"{UTAB.get((3, 3, 3, 3), float('nan')):+.1f}",
              flush=True)
    log = open(outbase + ".chan.jsonl", "w")
    log.write(json.dumps({
        "start": start, "etarget": ETARGET, "cimp": CIMP, "f3t": F3T,
        "zeta": ZETA, "z0": round(Z0, 3), "unseeds": UNSEEDS,
        "utab": len(UTAB), "aof": AOF, "br": BR, "bcf": BCF,
        "bc4": BC4, "cycles": cycles, "relax": relax,
        "episodes": episodes, "maxstep": MAXSTEP, "seed": seed}) + "\n")

    if DSIDE:
        pg = 1.0 - PHD - BCF_D - BC4_D

        def retube():
            r = build_orbit_tube(s, L)
            if r is None:
                return 0
            tube, f0ref = r
            s.set_worm_f0(tube, [0.0] * 6, OFFPEN, Z0, lmax=LMAX,
                          zeta=ZETA_D, aof=AOF, ph=PHD, pg=pg,
                          bcf=BCF_D, bc4=BC4_D, maxstep=MAXSTEP,
                          ucap_hi=UCAP_HI, ucap_lo=UCAP_LO, mu=MU,
                          f0_ref=f0ref)
            return len(tube)
        nt = retube()
        print(f"D-side episodes: single-orbit tube ({nt} states), "
              f"zeta={ZETA_D} ph={PHD} pg={pg:.4f} bcf={BCF_D} "
              f"bc4={BC4_D} mu={MU} maxstep={MAXSTEP}", flush=True)
        acc = {"of": [0, 0], "oi": [0, 0], "re": [0, 0], "cf": [0, 0],
               "c4": [0, 0]}
        undone = 0
        commits = 0
        t0 = time.time()
        for cyc in range(cycles):
            if relax > 0:
                s.run(sweeps=relax)
            if relax > 0:
                worm_L = W.Live(s.manifold,
                                driver=s.do_bistellar_move)
                globals()  # (Live rebuilt for tube builds below)
                L = worm_L
            eps = []
            for _ in range(episodes):
                r = s.worm_f0_episode()
                if r["changed"]:
                    # the tube's pin part self-adjusts (f-adaptive
                    # skeleton, recompiled at each episode open), so a
                    # commit no longer forces a rebuild -- retube only
                    # every RETUBE_EVERY commits to refresh the local
                    # m^2 content (each episode is balanced under its
                    # own frozen compiled U; the closed-sector measure
                    # is U-independent, so any cadence is unbiased)
                    L = W.Live(s.manifold, driver=s.do_bistellar_move)
                    commits += 1
                    if commits % RETUBE_EVERY == 0:
                        retube()
                acc["re"][0] += r["nH"] + r["nG"]
                acc["re"][1] += r["accH"] + r["accG"]
                if r["opened"]:
                    kind = "of" if r["opened"] == 1 else "oi"
                    acc[kind][0] += 1
                    acc[kind][1] += 1
                    if r["closed"] in ("cf", "c4"):
                        acc[r["closed"]][1] += 1
                        acc[r["closed"]][0] += 1
                    elif r["closed"] == "undone":
                        undone += 1
                eps.append({k: r[k] for k in
                            ("opened", "head", "steps", "closed", "dS",
                             "umax", "accH", "nH", "accG", "nG",
                             "zmin", "nZ4", "changed")})
            fv = [int(x) for x in s.manifold.f_vector]
            gap = fv[1] - 6.0 * fv[3] / ETARGET
            log.write(json.dumps({
                "cyc": cyc, "f": fv, "gap": round(gap, 3),
                "S": round(s.current_objective, 3), "ep": eps}) + "\n")
            log.flush()
            if cyc % 10 == 0 or cyc == cycles - 1:
                print(f"cyc {cyc:4d} f0={fv[0]} gap={gap:+6.2f} "
                      f"S={s.current_objective:8.2f} | "
                      + " ".join(f"{k}:{x[1]}/{x[0]}"
                                 for k, x in acc.items())
                      + f" undone={undone} ({time.time() - t0:.0f}s)",
                      flush=True)
        s.manifold.save(outbase + ".final.mfd")
        log.close()
        print(f"done: {outbase}.final.mfd (+ .chan.jsonl) "
              f"undone={undone}", flush=True)
        return

    worm = Worm(s, L, rng)
    t0 = time.time()
    if WLCYC and not UTAB_LOAD:
        vs0 = sorted({x for e in L.edeg for x in e})
        anchor_sample = vs0[:: max(1, len(vs0) // 100)]
        worm.wl_delta = WLD0
        stage = max(1, WLCYC // 6)
        for cyc in range(WLCYC):
            if relax > 0:
                s.run(sweeps=relax)
                worm.L = L = W.Live(s.manifold,
                                    driver=s.do_bistellar_move)
            for _ in range(episodes):
                worm.episode()
            if (cyc + 1) % stage == 0:
                worm.wl_delta /= 2.0
                off = anchor_utab(L, anchor_sample)
                a = worm.acc
                print(f"WL {cyc + 1:3d}/{WLCYC} delta="
                      f"{worm.wl_delta:.3f} |UTAB|={len(UTAB)} "
                      f"anchor{off:+.2f} | "
                      + " ".join(f"{k}:{x[1]}/{x[0]}"
                                 for k, x in a.items())
                      + f" undone={worm.undone} "
                      f"({time.time() - t0:.0f}s)", flush=True)
        worm.wl_delta = 0.0
        save_utab(outbase + ".utab.json")
        print(f"WL frozen -> {outbase}.utab.json; reloading start "
              f"state for production", flush=True)
        m = ddg.Manifold.load(start, 3)
        s, L = fresh(m)
        worm = Worm(s, L, rng)
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
                  + f" undone={worm.undone} ({time.time() - t0:.0f}s)",
                  flush=True)
    s.manifold.save(outbase + ".final.mfd")
    log.close()
    print(f"done: {outbase}.final.mfd (+ .chan.jsonl) "
          f"undone={worm.undone}", flush=True)


if __name__ == "__main__":
    main()
