"""Build the CHORD carrier's umbrella by REPLAYING measured catalysed
paths (the analogue of build_orbit_tube for the vertex collapse).

Why replay rather than design. Catalysis is reachable -- 95.6% of wide
head moves leave the flicker closable, and paths with net dS -4 to -8
exist at depth 2-3 -- but the walk almost never gets there under plain
Metropolis (accH 2/2015), because a catalysed path OPENS by creating a
helper flicker, uphill by +6..+23 on the rung ladder. That is a barrier
to a destination we have actually measured, which is exactly what an
umbrella is for. Every hand-designed umbrella in the f0 campaign failed
(entropy plateaus, deg-3 farming, one-way traps); the one that worked
replayed a real corridor. So: find catalysed paths, replay them, freeze
the cumulative dS keyed on the chord's local signature.

The key is the sorted multiset of degrees of every edge at either chord
endpoint. That is the right key because U must change exactly when the
wide head class fires: an edge (x,u) changes degree only when a tet on
it is touched, which needs both x and u in the support, so U changes iff
an endpoint is in the support -- restoring dU == 0 for every move
outside the head class.
"""
import json
import os
import sys

_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, os.path.join(_ROOT, "python"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import discrete_differential_geometry as ddg  # noqa: E402


def chord_sig(L, c):
    """The chord's local signature: sorted degrees of every edge at
    either endpoint (the chord's own edge included)."""
    nb = set()
    for e in c:
        for t in L.v2t[e]:
            nb |= set(t)
    nb -= set(c)
    out = []
    for u in nb:
        for e in c:
            d = L.edeg.get((min(e, u), max(e, u)))
            if d is not None:
                out.append(d)
    return tuple(sorted(out))


def apexes(L, face):
    fs = set(face)
    ts = [t for t in L.v2t[face[0]] if fs <= set(t)]
    return None if len(ts) != 2 else tuple(sorted(
        {x for t in ts for x in t} - fs))


def wide_head_moves(L, c):
    """Moves whose support MEETS the chord (>= 1 endpoint), over the
    union of the two endpoint stars. Excludes the 3->2 on the chord
    itself, which is the episode's close."""
    out, seenf, seene = [], set(), set()
    tets = {tuple(sorted(t)) for e in c for t in L.v2t[e]}
    for t in tets:
        for i in range(4):
            face = tuple(sorted(x for j, x in enumerate(t) if j != i))
            if face in seenf:
                continue
            seenf.add(face)
            ap = apexes(L, face)
            if ap is None or len(ap) != 2 or ap in L.edeg:
                continue
            if not (set(c) & (set(face) | set(ap))):
                continue
            out.append((face, ap))
        for i in range(4):
            for j in range(i + 1, 4):
                e = (t[i], t[j])
                if e in seene or e == tuple(sorted(c)):
                    continue
                seene.add(e)
                if L.edeg.get(e, 0) != 3:
                    continue
                ts = [q for q in L.v2t[e[0]] if e[0] in q and e[1] in q]
                lk = tuple(sorted({x for q in ts for x in q} - set(e)))
                if len(lk) != 3 or not (set(c) & (set(lk) | set(e))):
                    continue
                out.append((e, lk))
    return out


def closable(L, c):
    e = tuple(sorted(c))
    if L.edeg.get(e, 0) != 3:
        return None
    ts = [q for q in L.v2t[e[0]] if e[0] in q and e[1] in q]
    lk = tuple(sorted({x for q in ts for x in q} - set(e)))
    return lk if len(lk) == 3 else None


def build_chord_tube(s, L, seeds_faces, depth=3, nodecap=6000,
                     nseed=30, clamp_negative=True):
    """Search for catalysed paths and replay the best ones into a
    {signature: cumulative dS} table. Returns (table, n_paths)."""
    tab = {}
    npath = 0
    for face, ap in sorted(seeds_faces)[:nseed]:
        S0 = s.current_objective
        try:
            L.do(face, ap)
        except Exception:
            continue
        chord = tuple(sorted(ap))
        best = [None]
        nodes = [0]

        def rec(path, depth_left):
            nodes[0] += 1
            if nodes[0] > nodecap:
                return
            lk = closable(L, chord)
            if lk is not None and path:
                Sa = s.current_objective
                L.do(chord, lk)
                net = s.current_objective - S0
                L.do(lk, chord)
                assert abs(s.current_objective - Sa) < 1e-9
                if net < -1e-9 and (best[0] is None or net < best[0][1]):
                    best[0] = (list(path), net)
            if depth_left <= 0:
                return
            for cen, coc in wide_head_moves(L, chord):
                try:
                    L.do(cen, coc)
                except Exception:
                    continue
                path.append((cen, coc))
                rec(path, depth_left - 1)
                path.pop()
                L.do(coc, cen)

        rec([], depth)
        if best[0] is not None:
            # REPLAY: walk the winning path recording (signature, cum dS).
            # The reference is S1, the BARE-FLICKER state -- not S0. The
            # open sector begins and ends at the bare flicker, so that is
            # where U must vanish. Measuring from S0 instead loads the
            # flicker's own creation cost (+6..+23) onto every entry and
            # puts the umbrella's zero at the catalysed END of the path,
            # i.e. exactly backwards: the walk is then drawn away from
            # the one state it can close in, and the close (carrying -U)
            # is suppressed. Path SELECTION still scores against S0,
            # since a catalysed path has to pay for its own flicker.
            path, net = best[0]
            S1 = s.current_objective
            tab.setdefault(chord_sig(L, chord), 0.0)
            for cen, coc in path:
                L.do(cen, coc)
                cum = s.current_objective - S1
                sig = chord_sig(L, chord)
                # keep the CHEAPEST route to each signature
                if sig not in tab or cum < tab[sig]:
                    tab[sig] = cum
            for cen, coc in reversed(path):
                L.do(coc, cen)
            npath += 1
        lk = closable(L, chord)
        if lk is not None:
            L.do(chord, lk)
        assert abs(s.current_objective - S0) < 1e-6, "seed not restored"
    # No min-anchoring: the bare-flicker signature is already the zero,
    # by construction of the S1 reference above. A uniform shift is a
    # gauge (it trades against zeta2), but pinning it to min() would move
    # the zero off the closable state again.
    #
    # CLAMP the below-zero entries at that same zero. They come from
    # catalysed END states (net dS < 0), which have no business pricing a
    # just-born flicker: the open ratio carries +U, so a fresh flicker
    # whose signature happened to match a -24 entry had its open
    # suppressed by e^-24. Measured: unopened episodes rose 135 -> 175
    # once the negative entries were in play. Clamping keeps every
    # signature at or above the bare-flicker reference, so the umbrella
    # can only ever HELP the walk climb, never tax the opening.
    if clamp_negative:
        tab = {k: (v if v > 0.0 else 0.0) for k, v in tab.items()}
    return tab, npath


def save(tab, path):
    json.dump({",".join(map(str, k)): v for k, v in tab.items()},
              open(path, "w"))


def load(path):
    return {tuple(int(x) for x in k.split(",")): v
            for k, v in json.load(open(path)).items()}
