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
     CHORD CPH CPG CBCP CMAXSTEP CHAINK  (the chord carrier, off by
     default; CHORD = strict-closure chord episodes per cycle)

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
import discrete_differential_geometry as ddg
from . import worm_deg4_slide as W

# Positional CLI knobs (argv[3..5]) apply only when f0_worm itself is the
# script being run -- as a library import argv belongs to someone else
# (pytest, a driver CLI) and must not be parsed.
_AS_SCRIPT = os.path.basename(sys.argv[0]).startswith("f0_worm")
ETARGET = float(sys.argv[3]) if _AS_SCRIPT and len(sys.argv) > 3 else 5.1067907
CIMP = float(sys.argv[4]) if _AS_SCRIPT and len(sys.argv) > 4 else 0.7
F3T = int(sys.argv[5]) if _AS_SCRIPT and len(sys.argv) > 5 else 8704
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
AUTOZETA = bool(int(os.environ.get("AUTOZETA", "1")))  # calibrate zeta
ANCHOR_U = float(os.environ.get("ANCHOR_U", "20.6"))   # canonical U[3^4]
# --- CHORD carrier (the second bilocal carrier), interleaved per cycle.
# It conserves f0 EXACTLY (Delta f0 == 0 structurally), so it can never
# move the f0 read-out by itself; it is here as the fixed-f0 f3 mobility
# that unblocks the vertex channel's reshaper. Defaults reproduce the
# configuration whose two-sided balance was certified in
# scripts/defect_dynamics/twosided_chord.py -- do not drift from them
# without re-certifying.
CHORD = int(os.environ.get("CHORD", "0"))          # chord episodes/cycle
CPH = float(os.environ.get("CPH", "0.5"))          # chord head share
CPG = float(os.environ.get("CPG", "0.3"))          # chord global share
CBCP = float(os.environ.get("CBCP", "0.05"))       # chord close share
CMAXSTEP = int(os.environ.get("CMAXSTEP", "800"))  # chord step cap
CHAINK = int(os.environ.get("CHAINK", "20"))       # chord site window
BCP = float(os.environ.get("BCP", "0.002"))        # pair p_close
UFBSEEDS = int(os.environ.get("UFBSEEDS", "24"))   # seeds for the ufb fit
UFB_ON = bool(int(os.environ.get("UFB_ON", "1")))  # off-tube gradient


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
    from . import link_planner as LP
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
    from . import link_planner as LP
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
    # ANCHOR: pin U([3,3,3,3]) to a canonical value. A uniform shift of U
    # is a pure gauge choice -- it trades exactly against zeta -- but
    # leaving it free let U([3,3,3,3]) wander up to 10 units between
    # rebuilt corridors (each retube may find a different seed/orbit),
    # and log zeta* = dS14 - U4 - const inherited every bit of it:
    # measured log10 zeta spreads of 1.4 (quenched) and 2.2 (harvested),
    # i.e. acceptances swinging e^+-3 around 1 and alternating between
    # insert- and remove-favouring regimes. Anchored, zeta* tracks only
    # dS14 (sd 1.7). Returns the shift so the caller moves OFFPEN with
    # it and the tube's RELATIVE structure is untouched.
    shift = 0.0
    if (3, 3, 3, 3) in tab:
        shift = ANCHOR_U - tab[(3, 3, 3, 3)][0]
        tab = {ms: (cum + shift, d1, d3)
               for ms, (cum, d1, d3) in tab.items()}
    return tab, (fv[1], fv[3]), shift


def fit_ufb(s, L, nseeds=None, nodecap=UNC):
    """Fit the per-spoke linear fallback that prices OFF-TUBE stars.

    The orbit tube is a SINGLE collapse corridor -- 12 states from one
    vertex -- and the D config has always been called with ufb = 0, so
    wf0U returned the flat constant ufbc for everything off it. That is
    fine for the single-ball channel, which opens by inserting its own
    fresh z=4 head (already on the corridor), but it is fatal for the
    pair channel: the adopted ball is drawn from the whole crystal,
    lands on an ordinary Z12 star (the dominant TCP coordination), and
    sees a perfectly flat umbrella with no reason to collapse. Measured:
    the adopted ball's zmin stayed at its starting 12 in EVERY episode,
    so the transport closure -- which deletes that ball -- was
    unreachable no matter how zeta2 was tuned.

    So: keep the single-orbit tube as the EXACT on-corridor table (that
    is what avoids the cross-orbit min-agg traps build_umbrella warns
    about) and use a many-corridor least-squares fit only for the
    off-tube gradient. Both halves are free choices of U, so neither
    can bias the closed measure.

    Returns the 6 coefficients (n3, n4, n5, n6, n7plus, Zdeficit);
    build_umbrella also sets the module-level Z0 reference as a side
    effect. NOT f-adaptive (unlike the tube, whose pin part is
    recompiled per episode) -- again a free choice, so it costs walk
    efficiency as f drifts, never correctness."""
    _tab, coef = build_umbrella(s, L, nseeds or UFBSEEDS, nodecap)
    return list(coef)


def pin_part(f1, f3):
    """The action terms NONLINEAR in the f-vector (volume + edge pins) --
    the f-dependent half of every tube value (see sampler.wf0PinPart)."""
    return 0.1 * (f3 - F3T) ** 2 + (f1 - 6.0 * f3 / ETARGET) ** 2


def tube_u(tube, fref, f1, f3, ms, offpen=OFFPEN, ufb=None, z0=None):
    """Tube value at multiset ``ms``, compiled at the live f-vector
    exactly as the D engine does.

    Off-tube states take the SAME per-spoke linear fallback wf0U uses,
    ``ufbc + ufb[5](z0 - z) + sum ufb[i] n[i]`` with ``offpen`` as the
    constant -- pass ``ufb=None`` for the flat-plateau behaviour. This
    has to track wf0U exactly: once fit_ufb turned the plateau into a
    gradient on the D side, a Python model still returning the flat
    constant misprices every off-tube head, and calib_zeta2 is built
    entirely out of those values."""
    if ms not in tube:
        if ufb is None:
            return offpen
        n = [0] * 5
        for d in ms:
            n[min(max(int(d), 3), 7) - 3] += 1
        z0 = Z0 if z0 is None else z0
        u = offpen + ufb[5] * (z0 - len(ms))
        for i in range(5):
            u += ufb[i] * n[i]
        return max(min(u, UCAP_HI), UCAP_LO)
    ent = tube[ms]
    if not isinstance(ent, (tuple, list)):
        return ent          # plain-float table (build_umbrella): not
                            # f-adaptive, so no pin recompilation
    cum, d1, d3 = ent
    return (cum + (pin_part(f1 + d1, f3 + d3) - pin_part(f1, f3))
            - (pin_part(fref[0] + d1, fref[1] + d3) - pin_part(*fref)))


def calib_zeta(s, L, tube, fref, offpen=OFFPEN, ntet=24, ufb=None):
    """Auto-calibrate the open-sector fugacity zeta.

    From sampler.d (close-41 / open-insert, exact inverses):
      log a_c4 = -dS41 - U - log z + log(aoi/bc4) - log f3 - log lmax
      log a_oi = +log z - dS14 + U + log f3 + log lmax + log(bc4/aoi)
    so both legs are O(1) -- the condition for the insert/remove sector
    to actually TURN OVER -- iff
      log z* = dS14 - U([3,3,3,3]) - log f3 - log lmax - log(bc4/aoi).
    dS14 (the 1->4 insert cost at a uniformly random tet, which is what
    the engine proposes) is strongly state-dependent: it is dominated by
    the pins, so it falls steeply as the gap closes (~+40 at gap +10,
    ~+19 at gap 0). A fixed zeta therefore prices the sector correctly
    over a narrow band only -- outside it one leg saturates at alpha ~ 1
    and the other dies (measured: alpha_c4 ~ 1e-6 at gap 0 with
    zeta = 50, giving 207/207 abandoned insert episodes and ZERO
    committed f0 moves). zeta is a free ensemble parameter -- the
    closed-sector marginal is zeta-independent, the same theorem that
    licenses retubing -- so recalibrating it between episodes is exactly
    unbiased. Returns (zeta*, mean dS14, U[3,3,3,3])."""
    fv = [int(x) for x in s.manifold.f_vector]
    tets = [t for ts in L.v2t.values() for t in ts]
    top = max(x for t in tets for x in t)
    rng = random.Random(0xC0FFEE ^ fv[3])
    ds = []
    for i in range(ntet):
        t = rng.choice(tets)
        S0 = s.current_objective
        try:
            s.do_bistellar_move(list(t), [top + 1 + i])
        except Exception:
            continue
        ds.append(s.current_objective - S0)
        s.do_bistellar_move([top + 1 + i], sorted(t))
        assert abs(s.current_objective - S0) < 1e-6
    if not ds:
        return ZETA_D, float("nan"), float("nan")
    dS14 = sum(ds) / len(ds)
    u4 = tube_u(tube, fref, fv[1], fv[3], (3, 3, 3, 3), offpen, ufb)
    logz = (dS14 - u4 - math.log(fv[3]) - math.log(LMAX)
            - math.log(BC4_D / (1.0 - AOF)))
    return math.exp(max(min(logz, 30.0), -30.0)), dS14, u4


def calib_zeta2(s, L, tube, fref, offpen=OFFPEN, ntry=48,
                ph=None, pg=None, bcp=None, ufb=None):
    """Auto-calibrate the PAIR channel's fugacity zeta2.

    From sampler.d, the pair open is
      log a_open = log z2 - dS14 + u0 + u1
                   + log f3 + log lmax + log(W/wv) + log p_close
    with u0 the CREATED ball's umbrella at its fresh (3,3,3,3) star and
    u1 the ADOPTED ball's -- role-signed, so u0 = -U[3,3,3,3] and
    u1 = +U(star of the drawn vertex). The close carries the exact
    mirror, so both legs are O(1) iff

      log z2* = dS14 - u0 - u1
                - log f3 - log lmax - log(W/wv) - log p_close.

    W/wv is the drawn vertex's inverse seed probability under the
    engine's biased seed p(v) ~ exp(-mu(2 + deg/2)); it is sampled here
    the same way, so the average is over the SAME law the engine uses.
    Like zeta, zeta2 is a free ensemble parameter (the closed-sector
    marginal does not depend on it), so recalibrating between episodes
    is exactly unbiased -- and necessary, because the bracket is
    dominated by dS14, which is pin-driven and drifts as the gap closes.

    Returns (zeta2*, mean bracket, sd bracket). The sd is the number
    that matters: a scalar can only price the mean, so a draw sitting
    k units off has one leg at ~1 and its inverse at ~e^-k."""
    ph = PHD if ph is None else ph
    pg = (1.0 - PHD - BCF_D - BC4_D) if pg is None else pg
    pcl = BCP if bcp is None else bcp     # p_close IS bcp for the pair
    if not (0.0 < pcl < 1.0) or not (ph + pg > 0.0):
        return ZETA_D, float("nan"), float("nan")
    fv = [int(x) for x in s.manifold.f_vector]
    f1, f3 = fv[1], fv[3]
    tets = [t for ts in L.v2t.values() for t in ts]
    top = max(x for t in tets for x in t)
    verts = sorted({x for e in L.edeg for x in e})
    # the engine's seed weight is exp(-mu(2 + deg/2)) with deg =
    # Manifold._vertexDegrees[v], which counts FACETS in the star (it is
    # bumped once per tet on insert), NOT link vertices. len(L.v2t[v]) is
    # that same count. Using the neighbour count instead (mean 2f1/f0 =
    # 13.4 against the tet mean 4f3/f0 = 22.8) models a different
    # proposal than the one being priced and threw log(W/wv) off by ~17,
    # which is exactly the shortfall the close was showing.
    wt = {v: math.exp(-MU * (2.0 + 0.5 * len(L.v2t[v]))) for v in verts}
    W = sum(wt.values())
    u0 = -tube_u(tube, fref, f1, f3, (3, 3, 3, 3), offpen, ufb)
    rng = random.Random(0xBEEF ^ f3)
    brk = []
    for i in range(ntry):
        t = rng.choice(tets)
        S0 = s.current_objective
        try:
            s.do_bistellar_move(list(t), [top + 1 + i])
        except Exception:
            continue
        dS14 = s.current_objective - S0
        s.do_bistellar_move([top + 1 + i], sorted(t))
        assert abs(s.current_objective - S0) < 1e-6
        r, acc, vf = rng.random() * W, 0.0, None
        for v in verts:
            acc += wt[v]
            if acc >= r:
                vf = v
                break
        if vf is None or vf in t:
            continue
        u1 = tube_u(tube, fref, f1, f3, spoke_ms(L, vf), offpen, ufb)
        brk.append(dS14 - u0 - u1 - math.log(f3) - math.log(LMAX)
                   - math.log(W / wt[vf]) - math.log(pcl))
    if not brk:
        return ZETA_D, float("nan"), float("nan")
    mean = sum(brk) / len(brk)
    sd = (sum((x - mean) ** 2 for x in brk) / len(brk)) ** 0.5
    return math.exp(max(min(mean, 300.0), -300.0)), mean, sd


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

        state = {"tube": None, "fref": None, "zeta": ZETA_D,
                 "dS14": float("nan"), "u4": float("nan"),
                 "offpen": OFFPEN}

        def push(ph_, pg_, bcf_, bc4_, maxstep_):
            """Push the current tube under a given step mixture. The two
            carriers share one WormF0Params, so each channel's block
            swaps the mixture in and the other swaps it back. That is
            safe: every episode is detailed-balanced under whatever
            (U, zeta, mixture) is frozen at its open, and the closed
            measure does not depend on any of them -- the same master
            theorem that licenses between-episode recalibration."""
            s.set_worm_f0(state["tube"], [0.0] * 6, state["offpen"], Z0,
                          lmax=LMAX, zeta=state["zeta"], aof=AOF,
                          ph=ph_, pg=pg_, bcf=bcf_, bc4=bc4_,
                          maxstep=maxstep_, ucap_hi=UCAP_HI,
                          ucap_lo=UCAP_LO, mu=MU,
                          f0_ref=state["fref"])

        def push_vertex():
            push(PHD, pg, BCF_D, BC4_D, MAXSTEP)

        def push_chord():
            push(CPH, CPG, 1e-4, 1.0 - CPH - CPG - 1e-4, CMAXSTEP)

        def reconfig():
            """Push the current tube at a freshly calibrated zeta. The
            insert cost drifts with the state, so zeta is recalibrated
            here as well as at every rebuild (both are exactly
            unbiased: the closed measure is U- and zeta-independent)."""
            if AUTOZETA:
                z, d14, u4 = calib_zeta(s, L, state["tube"],
                                        state["fref"],
                                        state["offpen"])
                state.update(zeta=z, dS14=d14, u4=u4)
            push_vertex()

        def retube():
            r = build_orbit_tube(s, L)
            if r is None:
                return 0
            state["tube"], state["fref"], sh = r
            # move the off-tube flat value with the anchor shift so the
            # tube's relative structure is untouched (uniform U shifts
            # are gauge; only the tube-vs-off-tube gap is physical)
            state["offpen"] = OFFPEN + sh
            reconfig()
            return len(state["tube"])
        nt = retube()
        print(f"D-side episodes: single-orbit tube ({nt} states), "
              f"zeta={state['zeta']:.3e} (auto={int(AUTOZETA)}, "
              f"dS14={state['dS14']:+.2f} U4={state['u4']:+.2f}) "
              f"ph={PHD} pg={pg:.4f} bcf={BCF_D} "
              f"bc4={BC4_D} mu={MU} maxstep={MAXSTEP}", flush=True)
        if CHORD:
            s.set_worm_pair(zeta2=float("nan"), bcp=CBCP,
                            chain_k=CHAINK)
            print(f"chord carrier: {CHORD} episodes/cycle, ph={CPH} "
                  f"pg={CPG} bcp={CBCP} maxstep={CMAXSTEP} "
                  f"chain_k={CHAINK} (auto zeta2 + p_close)",
                  flush=True)
        acc = {"of": [0, 0], "oi": [0, 0], "re": [0, 0], "cf": [0, 0],
               "c4": [0, 0]}
        undone = 0
        commits = 0
        ccommit = 0
        ccensus = {}
        t0 = time.time()
        for cyc in range(cycles):
            if relax > 0:
                s.run(sweeps=relax)
            if relax > 0:
                worm_L = W.Live(s.manifold,
                                driver=s.do_bistellar_move)
                globals()  # (Live rebuilt for tube builds below)
                L = worm_L
            # CHORD block first: it is the fixed-f0 solvent, so let it
            # act before the vertex channel tries to remove anything.
            if CHORD:
                push_chord()
                cdirty = False
                for _ in range(CHORD):
                    f0_before = int(s.manifold.f_vector[0])
                    rc = s.worm_chord_strict_episode()
                    if rc["changed"]:
                        # CONTROL, asserted every commit: this channel
                        # conserves f0 exactly. If it ever does not, the
                        # f0 read-out below is contaminated and the run
                        # is void -- so fail loudly rather than quietly
                        # attribute a chord artefact to the vertex
                        # channel's free energy.
                        assert int(s.manifold.f_vector[0]) == f0_before, \
                            "chord channel changed f0 -- channel is bad"
                        ccommit += 1
                        ccensus[rc["df"]] = ccensus.get(rc["df"], 0) + 1
                        L = W.Live(s.manifold,
                                   driver=s.do_bistellar_move)
                        cdirty = True
                # a chord commit moves f3, hence the pin gap, hence the
                # vertex channel's insert cost -- recalibrate zeta
                if cdirty and AUTOZETA:
                    reconfig()
                else:
                    push_vertex()
            eps = []
            for _ in range(episodes):
                r = s.worm_f0_episode()
                if r["changed"]:
                    # the tube's pin part self-adjusts (f-adaptive
                    # skeleton, recompiled at each episode open), so a
                    # commit no longer forces a rebuild -- retube only
                    # every RETUBE_EVERY commits to refresh the local
                    # m^2 content. zeta IS recalibrated every commit
                    # (cheap: 24 insert/undo probes) -- the insert cost
                    # drifts continuously with the pin gap, and a stale
                    # zeta silently kills one leg of the insert/remove
                    # sector. Both are exactly unbiased: the closed
                    # measure is U- and zeta-independent.
                    L = W.Live(s.manifold, driver=s.do_bistellar_move)
                    commits += 1
                    if commits % RETUBE_EVERY == 0:
                        retube()
                    elif AUTOZETA:
                        reconfig()
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
                "S": round(s.current_objective, 3),
                "zeta": state["zeta"], "dS14": state["dS14"],
                "u4": state["u4"], "cc": ccommit, "ep": eps}) + "\n")
            log.flush()
            if cyc % 10 == 0 or cyc == cycles - 1:
                print(f"cyc {cyc:4d} f0={fv[0]} gap={gap:+6.2f} "
                      f"S={s.current_objective:8.2f} | "
                      + " ".join(f"{k}:{x[1]}/{x[0]}"
                                 for k, x in acc.items())
                      + f" undone={undone} z={state['zeta']:.1e}"
                      + (f" ch={ccommit}" if CHORD else "")
                      + f" ({time.time() - t0:.0f}s)", flush=True)
        s.manifold.save(outbase + ".final.mfd")
        if CHORD:
            log.write(json.dumps({"chord_census":
                                  {str(k): v
                                   for k, v in ccensus.items()},
                                  "chord_commits": ccommit}) + "\n")
        log.close()
        print(f"done: {outbase}.final.mfd (+ .chan.jsonl) "
              f"undone={undone}"
              + (f" chord_commits={ccommit}" if CHORD else ""),
              flush=True)
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
