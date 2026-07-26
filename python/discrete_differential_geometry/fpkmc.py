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
