"""FPKMC infrastructure (notes/FPKMC_DESIGN.md, M1).

Library-grade pieces of the first-passage / heat-bath samplers that are
pure bookkeeping over the D core:

  - the REQUIRED startup guard: face -> apex-pair injectivity on the
    reference crystal (a duplicate would make chord keys ambiguous; by
    no-halo, reference-level injectivity covers every clear region of
    every state, so one scan suffices permanently);
  - face <-> apex-pair maps (chord key <-> creation face);
  - single-site surveys (any pure-knot state is a creation state: knot
    state <-> chord <-> unique face, BY the injectivity guard);
  - the clean slide graph builder: BFS over pure-knot states with exact
    per-edge dS from ddg_sampler_site_survey. Nodes are chords; edges are
    clean slides (species-preserving, so the state stays a pure knot).
    Verified properties this relies on (2026-07-26, R m4): every legal
    slide has a legal inverse with dS' = -dS (24/24 spot checks), and
    detailed balance holds with pi ~ e^{-S} (M0 proof, given the exact
    per-edge dS).
  - the proposal-rate constant nu: P(propose a specific (chord, slot)) =
    slide_prob * (3/N3) * (1/6) * (1/12) per attempted move, read off
    sampler.d's slide branch (a deg-3 chord lies in exactly 3 facets).

The dirty-slide graph (species-changing edges; the physical default
channel) needs live-state menus and is deliberately NOT built here -- it
arrives with the D-side graph scan (design M2+).
"""
from __future__ import annotations

import ctypes
from collections import deque
from itertools import combinations

import numpy as np

from ._manifold import Manifold
from ._sampler import ManifoldSampler, SLIDE_SLOTS


def face_apex_maps(manifold: Manifold):
    """(face -> apex pair, apex pair -> face) over all interior triangles.

    Raises RuntimeError if the apex-pair map is NOT injective -- the
    required M1 guard. On a crystal that fails, chord keys are ambiguous
    and the samplers must key states by (chord, face) instead; that mode
    is not implemented (design R3, documented fallback).
    """
    F = np.asarray(manifold.facets())
    face_ap: dict[tuple, list] = {}
    for t in F:
        ts = sorted(int(x) for x in t)
        for f in combinations(ts, 3):
            face_ap.setdefault(f, []).append(
                next(v for v in ts if v not in f))
    pair_of = {}
    face_of = {}
    for f, aps in face_ap.items():
        if len(aps) != 2:
            raise RuntimeError(f"face {f} has {len(aps)} apexes "
                               f"(not a closed 3-manifold?)")
        p = tuple(sorted(aps))
        if p in face_of:
            raise RuntimeError(
                f"apex-pair injectivity FAILED: pair {p} is the apex set "
                f"of both {face_of[p]} and {f}. Chord keys are ambiguous "
                f"on this crystal; the samplers require the (chord, face) "
                f"re-keying fallback (notes/FPKMC_DESIGN.md R3), which is "
                f"not implemented.")
        pair_of[f] = p
        face_of[p] = f
    return pair_of, face_of


def nu_per_attempt(slide_prob: float, n3: int) -> float:
    """P(propose one specific (chord, slot)) per attempted move."""
    return slide_prob * (3.0 / n3) / 6.0 / SLIDE_SLOTS


def survey_chord(sampler: ManifoldSampler, face, pair) -> dict:
    """site_survey of the single pure-knot state with chord `pair` (the
    knot created by the 2->3 on `face`). Returns the one-window row."""
    chain5 = [pair[0], face[0], face[1], face[2], pair[1]]
    sv = sampler.site_survey(chain5)
    return {k: v[0] for k, v in sv.items()}


class CleanKnotGraph:
    """The pure-knot slide graph, built lazily by BFS from a seed face.

    Nodes: chords (sorted pairs), each the unique apex pair of a crystal
    face. Edges: CLEAN slides, with exact dS from the D survey. The
    stationary law on the graph is pi(chord) ~ e^{-S1(chord)} (M0 + the
    verified edge antisymmetry), where S1 is dS_create of the node.

    `blocked` (optional set of vertices): a node whose 5 knot vertices
    tet-touch the blocked set is a boundary/docking node -- expanded but
    not traversed. This is the segment machinery: pass the other defects'
    vertex reach to get the clear component around the seed.
    """

    def __init__(self, manifold: Manifold, sampler: ManifoldSampler,
                 blocked: set | None = None):
        self.pair_of, self.face_of = face_apex_maps(manifold)
        self.sampler = sampler
        self.blocked = blocked or set()
        self.S1 = {}          # chord -> creation action (site energy)
        self.edges = {}       # chord -> list of (dest chord, dS, slot)
        self.boundary = set()

    def _is_blocked(self, chord) -> bool:
        if not self.blocked:
            return False
        face = self.face_of[chord]
        verts = set(chord) | set(face)
        return bool(verts & self.blocked)

    def expand(self, seed_chord, max_nodes: int = 2000) -> int:
        """BFS the clean component containing `seed_chord`."""
        q = deque([tuple(sorted(seed_chord))])
        while q and len(self.S1) < max_nodes:
            ch = q.popleft()
            if ch in self.S1 or ch in self.boundary:
                continue
            if self._is_blocked(ch):
                self.boundary.add(ch)
                continue
            face = self.face_of.get(ch)
            if face is None:
                self.boundary.add(ch)     # not a creation state (shouldn't
                continue                  # happen for clean destinations)
            row = survey_chord(self.sampler, face, ch)
            if np.isnan(row["dS_create"]):
                self.boundary.add(ch)
                continue
            self.S1[ch] = float(row["dS_create"])
            out = []
            for slot in range(SLIDE_SLOTS):
                if row["slot_clean"][slot] != 1.0:
                    continue
                d = tuple(sorted(int(x) for x in row["slot_dest"][slot]))
                out.append((d, float(row["slot_dS"][slot]), slot))
                if d not in self.S1 and d not in self.boundary:
                    q.append(d)
            self.edges[ch] = out
        return len(self.S1)

    def check_antisymmetry(self, tol: float = 1e-9) -> int:
        """Every internal edge must appear in both directions with
        dS' = -dS. Returns the number of violations (0 required)."""
        bad = 0
        for a, outs in self.edges.items():
            for b, dS, _ in outs:
                if b not in self.edges:
                    continue
                back = [x for x in self.edges[b] if x[0] == a]
                if not back or all(abs(x[1] + dS) > tol for x in back):
                    bad += 1
        return bad

    def check_consistency(self, tol: float = 1e-9) -> int:
        """dS along an edge must equal S1(dest) - S1(src) (the slide is a
        pure site-energy difference on the pure-knot graph). Violations
        would mean the graph model is missing state. Returns count."""
        bad = 0
        for a, outs in self.edges.items():
            for b, dS, _ in outs:
                if b in self.S1 and abs((self.S1[b] - self.S1[a]) - dS) > tol:
                    bad += 1
        return bad


# ---------------------------------------------------------------------------
# HB: the Metropolized-ball heat-bath driver (design M2)
# ---------------------------------------------------------------------------

class HBDriver:
    """Metropolized-ball independence sampler on the slide graph.

    One step: scan the ball B(x) around the current state (dS_max window W,
    depth D), draw y from pi restricted to B(x), move there by replaying the
    slide path, scan B(y), and accept with

        alpha = min(1, e^{d} * Z(x) / Z(y)),   d = dS_x(y)

    (independence-sampler Hastings ratio; requires x identifiable in B(y),
    else the proposal is rejected). On rejection the path is walked back
    through the verified inverse slides. Exactness rests on the verified
    edge antisymmetry and the per-edge dS being exact.

    v1 scope: single-defect sector -- every path node must be single-chord
    (asserted). A running-action audit (`audit_every`) recomputes the true
    action change against the session start and errors on drift, so a
    mis-identified inverse cannot corrupt silently.
    """

    def __init__(self, sampler, chord, dS_max=8.0, max_depth=5,
                 rng=None, audit_every=50, oracle=None):
        self.s = sampler
        self.chord = tuple(sorted(int(c) for c in chord))
        self.W = float(dS_max)
        self.D = int(max_depth)
        self.rng = rng or np.random.default_rng()
        self.audit_every = audit_every
        self.oracle = oracle          # callable() -> true action rel. start
        self.S_rel = 0.0              # running action vs session start
        self.nstep = self.naccept = self.nself = 0

    # -- helpers ------------------------------------------------------------

    def _scan(self):
        return self.s.slide_graph_scan(self.chord, dS_max=self.W,
                                       max_depth=self.D)

    @staticmethod
    def _parents(g):
        par = {}
        for i in range(len(g["edge_dS"])):
            a, b = int(g["edge_src"][i]), int(g["edge_dst"][i])
            if b not in par and b != 0:
                par[b] = (a, i)
        return par

    @staticmethod
    def _path(g, par, t):
        """Edge indices root -> node t."""
        out = []
        while t != 0:
            a, i = par[t]
            out.append(i)
            t = a
        return out[::-1]

    def _chord_deg3(self, ch):
        import ctypes as _ct
        from ._dlang import _lib as _l
        arr = (_ct.c_int * 2)(int(ch[0]), int(ch[1]))
        try:
            return int(_l.ddg_sampler_degree(self.s._handle, arr, 2)) == 3
        except RuntimeError:
            return False

    def _replay(self, g, path):
        """Apply the slides along `path` (edge indices in g). Returns the
        total committed dS."""
        tot = 0.0
        for i in path:
            ch = (int(g["edge_chord"][i][0]), int(g["edge_chord"][i][1]))
            dS = self.s.slide_at(ch[0], ch[1], int(g["edge_slot"][i]),
                                 commit=True)
            if dS is None:
                raise RuntimeError("HB replay: recorded slide invalid")
            tot += dS
        return tot

    def _walk_back(self, g, path):
        """Invert the slides of `path` in reverse order. The exact inverse
        at each step is identified BEFORE committing: it is the slot at the
        arrival chord whose own arrival (from the frame decode, via
        slide_at2) is the predecessor's chord and whose dS is the exact
        negation. No guessing, no impostors."""
        for i in reversed(path):
            src = int(g["edge_src"][i])
            dst = int(g["edge_dst"][i])
            arr_ch = (int(g["chord"][dst][0]), int(g["chord"][dst][1]))
            src_ch = tuple(sorted((int(g["chord"][src][0]),
                                   int(g["chord"][src][1]))))
            want = -float(g["edge_dS"][i])
            done = False
            for slot in range(SLIDE_SLOTS):
                dS, arr = self.s.slide_at2(arr_ch[0], arr_ch[1], slot,
                                           commit=False)
                if dS is None or arr is None:
                    continue
                if tuple(sorted(arr)) != src_ch or abs(dS - want) > 1e-9:
                    continue
                self.s.slide_at(arr_ch[0], arr_ch[1], slot, commit=True)
                done = True
                break
            if not done:
                raise RuntimeError("HB walk-back: inverse not found")

    # -- one HB step --------------------------------------------------------

    def step(self):
        gx = self._scan()
        w = np.exp(-np.asarray(gx["dS"]))
        Zx = float(w.sum())
        i = int(self.rng.choice(len(w), p=w / Zx))
        self.nstep += 1
        if i == 0:
            self.nself += 1
            self.naccept += 1
            return True
        if int(gx["n_chords"][i]) != 1:
            return False            # outside v1 scope: treat as reject
        d = float(gx["dS"][i])
        par = self._parents(gx)
        path = self._path(gx, par, i)
        if any(int(gx["n_chords"][int(gx["edge_dst"][j])]) != 1
               for j in path):
            return False
        committed = self._replay(gx, path)
        assert abs(committed - d) < 1e-6, "path dS != node dS"
        y_chord = (int(gx["chord"][i][0]), int(gx["chord"][i][1]))
        old_chord = self.chord
        self.chord = tuple(sorted(y_chord))
        gy = self._scan()
        # identify x in B(y)
        cand = [j for j in range(len(gy["dS"]))
                if tuple(sorted((int(gy["chord"][j][0]),
                                 int(gy["chord"][j][1])))) == old_chord
                and abs(float(gy["dS"][j]) + d) < 1e-9]
        if len(cand) > 1:
            raise RuntimeError("HB: ambiguous x in B(y)")
        if not cand:
            alpha = 0.0
        else:
            Zy = float(np.exp(-np.asarray(gy["dS"])).sum())
            alpha = min(1.0, np.exp(d) * Zx / Zy)
        if self.rng.random() < alpha:
            self.naccept += 1
            self.S_rel += d
            accepted = True
        else:
            self._walk_back(gx, path)
            self.chord = old_chord
            accepted = False
        if self.audit_every and self.nstep % self.audit_every == 0 \
                and self.oracle is not None:
            true_rel = self.oracle()
            if abs(true_rel - self.S_rel) > 1e-6:
                raise RuntimeError(
                    f"HB audit FAILED: running action {self.S_rel:.9f} vs "
                    f"true {true_rel:.9f}")
        return accepted


# ---------------------------------------------------------------------------
# FP: first-passage flights on the slide graph (design M3, v1 FROZEN)
# ---------------------------------------------------------------------------

ABSORB_DOCK = "dock"
ABSORB_DS = "dS_frontier"
ABSORB_DEPTH = "depth_frontier"
ABSORB_MULTI = "multichord"


class FPFlight:
    """One scanned FP domain: node classification + exact flight sampling.

    Built from one blocked-aware graph scan. Interior nodes are exactly
    the fully expanded ones (depth < max_depth, dS <= dS_max, single
    chord, dock == 0) -- by the scan's complete-interior guarantee, every
    legal slide out of an interior node appears in the edge list, so the
    sampler's slide channel restricted to the domain is EXACTLY: per
    attempted move, proposal edge e (a specific (chord, slot)) commits
    with probability nu * min(1, e^{-dS_e}), mutually exclusively; all
    other outcomes hold. A flight simulates the embedded jump chain with
    geometric holding times (exact, attempted-move units) to absorption.
    Every non-interior node is absorbing, labeled dock / dS_frontier /
    depth_frontier / multichord.

    Scan-tree paths to any node have interior-only intermediates (only
    expanded nodes create edges), so absorbing states can always be
    materialized by replaying tree slides.
    """

    def __init__(self, g, dS_max, max_depth, nu):
        self.g = g
        self.nu = float(nu)
        n = len(g["dS"])
        interior = ((g["depth"] < max_depth) & (g["dS"] <= dS_max)
                    & (g["n_chords"] == 1) & (g["dock"] == 0))
        self.interior = interior
        reason = np.full(n, "", dtype=object)
        single = g["n_chords"] == 1
        reason[(g["dock"] == 1)] = ABSORB_DOCK
        reason[(g["dock"] == 0) & ~single] = ABSORB_MULTI
        reason[(g["dock"] == 0) & single & (g["dS"] > dS_max)] = ABSORB_DS
        reason[(g["dock"] == 0) & single & (g["dS"] <= dS_max)
               & (g["depth"] >= max_depth)] = ABSORB_DEPTH
        reason[interior] = ""
        self.reason = reason
        acc = self.nu * np.minimum(1.0, np.exp(-np.asarray(g["edge_dS"])))
        self.out = {int(i): [] for i in np.nonzero(interior)[0]}
        for k in range(len(g["edge_src"])):
            i = int(g["edge_src"][k])
            if interior[i]:
                self.out[i].append((int(g["edge_dst"][k]), float(acc[k])))
        for i, es in self.out.items():
            if not es:
                raise RuntimeError(f"interior node {i} has no out-edges")
        self.escape = {i: sum(r for _, r in es)
                       for i, es in self.out.items()}

    def sample(self, rng, start=0, max_jumps=1_000_000):
        """One exact flight from `start` to absorption.
        Returns (absorbing node, time in attempted moves, n jumps)."""
        j, t, jumps = int(start), 0, 0
        while self.interior[j]:
            es = self.out[j]
            t += int(rng.geometric(self.escape[j]))
            u = rng.random() * self.escape[j]
            c = 0.0
            for dst, r in es:
                c += r
                if u <= c:
                    j = dst
                    break
            jumps += 1
            if jumps > max_jumps:
                raise RuntimeError("FP flight: jump budget exceeded")
        return j, t, jumps

    def _jump_matrices(self):
        idx = np.nonzero(self.interior)[0]
        pos = {int(v): c for c, v in enumerate(idx)}
        m = len(idx)
        absorbing = [int(v) for v in np.nonzero(~self.interior)[0]]
        apos = {v: c for c, v in enumerate(absorbing)}
        P = np.zeros((m, m))
        R = np.zeros((m, len(absorbing)))
        for i in idx:
            i = int(i)
            esc = self.escape[i]
            for dst, r in self.out[i]:
                if dst in pos:
                    P[pos[i], pos[dst]] += r / esc
                else:
                    R[pos[i], apos[dst]] += r / esc
        return idx, pos, absorbing, P, R

    def splitting_exact(self, start=0):
        """Exact absorption law from `start` (dense solve on the jump
        chain; holding probabilities drop out -- M0 proof E).
        Returns {absorbing node: probability}."""
        idx, pos, absorbing, P, R = self._jump_matrices()
        B = np.linalg.solve(np.eye(len(idx)) - P, R)
        row = B[pos[int(start)]]
        return {a: float(p) for a, p in zip(absorbing, row) if p > 0}

    def mean_time_exact(self, start=0):
        """Exact mean absorption time from `start`, in attempted moves:
        sum over interior visits of the geometric holding mean."""
        idx, pos, absorbing, P, R = self._jump_matrices()
        # expected visits to each interior node starting from `start`
        v = np.linalg.solve((np.eye(len(idx)) - P).T,
                            np.eye(len(idx))[pos[int(start)]])
        return float(sum(v[c] / self.escape[int(i)]
                         for c, i in enumerate(idx)))


class FPDriver:
    """v1 FROZEN first-passage driver (design M3): repeated FP flights
    with the other defect's vertices blocked; between flights the exit
    state is materialized by replaying its scan-tree path. Stops when a
    flight absorbs anywhere that is not a clean single-chord state
    (dock / multichord) -- contact resolution is the caller's job (§8).

    nu must be the per-(chord,slot) proposal probability per attempted
    move of the reference dynamics (fpkmc.nu_per_attempt).
    """

    def __init__(self, sampler, chord, blocked_verts, nu, dS_max=5.0,
                 max_depth=3, rng=None, oracle=None):
        self.s = sampler
        self.chord = tuple(sorted(int(c) for c in chord))
        self.blocked = [int(v) for v in blocked_verts]
        self.nu = float(nu)
        self.W = float(dS_max)
        self.D = int(max_depth)
        self.rng = rng or np.random.default_rng()
        self.oracle = oracle
        self.S_rel = 0.0
        self.t = 0                  # attempted-move clock
        self.nflight = 0
        self.recs = []              # (src chord, dst chord, dS) of every
                                    # materialized slide, for walk-back

    def flight(self):
        """One flight. Returns (reason, exit node dict) where reason ''
        never occurs (exits are absorbing) and dict carries chord/dS."""
        g = self.s.slide_graph_scan(self.chord, dS_max=self.W,
                                    max_depth=self.D,
                                    blocked_verts=self.blocked)
        fl = FPFlight(g, self.W, self.D, self.nu)
        j, t, _ = fl.sample(self.rng)
        par = HBDriver._parents(g)
        path = HBDriver._path(g, par, j)
        for i in path:
            ch = (int(g["edge_chord"][i][0]), int(g["edge_chord"][i][1]))
            dS, arr = self.s.slide_at2(ch[0], ch[1],
                                       int(g["edge_slot"][i]), commit=True)
            if dS is None or arr is None:
                raise RuntimeError("FP materialize: recorded slide invalid")
            self.recs.append((tuple(sorted(ch)),
                              tuple(sorted(int(x) for x in arr)),
                              float(dS)))
        self.t += t
        self.S_rel += float(g["dS"][j])
        self.nflight += 1
        reason = fl.reason[j]
        if int(g["n_chords"][j]) == 1:
            self.chord = tuple(sorted((int(g["chord"][j][0]),
                                       int(g["chord"][j][1]))))
        if self.oracle is not None:
            true_rel = self.oracle()
            if abs(true_rel - self.S_rel) > 1e-6:
                raise RuntimeError(
                    f"FP audit FAILED: running action {self.S_rel:.9f} vs "
                    f"true {true_rel:.9f}")
        return reason, {"node": int(j), "chord": self.chord,
                        "dS": float(g["dS"][j]), "t": t}
