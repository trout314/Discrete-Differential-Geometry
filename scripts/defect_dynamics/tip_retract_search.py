"""Endpoint-targeted Pachner search for the deg-4 tip RETRACT move, run on a
REAL gas snapshot (not pristine crystal).

Targets: (a) a lone deg-4 edge (valence 1-1, the dominant mobile species) and
(b) the tip edge of a longer deg-4 line.  Goal predicate: the illegal-edge set
in the core equals the initial set MINUS the tip edge -- i.e. the tip edge
healed and nothing new illegal.  Intermediates may be arbitrarily dirty
(dS depends only on endpoints); we cap only the transient illegal excess to
bound branching.  Every goal state is scored with the FULL shape action
(EDQ + zleg*U(n6) + cimp*m^2), in acceptance units (x LAM).

The inverse of a found RETRACT is the EXTEND at the arrival site; the hop
(transport) is retract o extend, whose dS ~ cancels by endpoint congruence.
"""
import os, sys, itertools, time
from collections import Counter, defaultdict
import numpy as np

_ROOT = os.path.normpath(os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "..", ".."))
sys.path.insert(0, os.path.join(_ROOT, "python"))
sys.path.insert(0, os.path.join(_ROOT, "scripts/defect_dynamics"))
import discrete_differential_geometry as ddg
import worm_moves as wm

SNAP = os.path.join(_ROOT, "data/mgas/lam35r_snap15000.mfd")
LAM, ESTAR, ZLEG, CIMP = 0.35, 5.105025, 0.6, 1.0
R_BALL, R_CORE = 3, 2
MAX_DEPTH = int(os.environ.get("MAXD", "6"))
ILL_EXCESS = int(os.environ.get("ILLX", "3"))    # transient illegal-count cap
NODE_CAP = int(os.environ.get("NODES", "4000000"))

def E(u, v): return (u, v) if u < v else (v, u)
def _legal(d): return d in (5, 6)

class Patch:
    def __init__(self, tets):
        self.tets = set(tets)
        self.edeg = Counter()
        for t in self.tets:
            for e in itertools.combinations(sorted(t), 2):
                self.edeg[e] += 1
        self.touched = set()

    def _rm(self, t):
        self.tets.discard(frozenset(t))
        for e in itertools.combinations(sorted(t), 2):
            self.edeg[e] -= 1; self.touched.add(e)
            if self.edeg[e] == 0: del self.edeg[e]

    def _add(self, t):
        self.tets.add(frozenset(t))
        for e in itertools.combinations(sorted(t), 2):
            self.edeg[e] += 1; self.touched.add(e)

    def apply_23(self, tri, x, y):
        a, b, c = tri
        self._rm((a, b, c, x)); self._rm((a, b, c, y))
        self._add((x, y, a, b)); self._add((x, y, b, c)); self._add((x, y, a, c))

    def undo_23(self, tri, x, y):
        a, b, c = tri
        self._rm((x, y, a, b)); self._rm((x, y, b, c)); self._rm((x, y, a, c))
        self._add((a, b, c, x)); self._add((a, b, c, y))

    def apply_32(self, u, w, a, b, c):
        self._rm((u, w, a, b)); self._rm((u, w, b, c)); self._rm((u, w, a, c))
        self._add((a, b, c, u)); self._add((a, b, c, w))

    def undo_32(self, u, w, a, b, c):
        self._rm((a, b, c, u)); self._rm((a, b, c, w))
        self._add((u, w, a, b)); self._add((u, w, b, c)); self._add((u, w, a, c))

    def moves(self, core, focus):
        """2->3 and 3->2 confined to core, touching `focus` vertices."""
        out = []
        face_ap = defaultdict(list)
        for t in self.tets:
            if not t <= core: continue
            st = tuple(sorted(t))
            for tri in itertools.combinations(st, 3):
                face_ap[tri].append(next(v for v in st if v not in tri))
        for tri, aps in face_ap.items():
            if len(aps) == 2:
                x, y = sorted(aps)
                if self.edeg.get(E(x, y), 0) == 0 and (set(tri) | {x, y}) & focus:
                    out.append(("23", tri, x, y))
        for e, d in list(self.edeg.items()):
            if d != 3 or e[0] not in core or e[1] not in core: continue
            u, w = e
            link = sorted({v for t in self.tets if u in t and w in t
                           for v in t if v not in (u, w)})
            if len(link) != 3: continue
            a, b, c = link
            if not (a in core and b in core and c in core): continue
            if any({a, b, c} <= t for t in self.tets): continue
            if {u, w, a, b, c} & focus:
                out.append(("32", (u, w), a, b, c))
        return out

def ball(F, seeds, R):
    adj = defaultdict(set)
    for t in F:
        for u, v in itertools.combinations(t, 2):
            adj[u].add(v); adj[v].add(u)
    seen = set(seeds); frontier = set(seeds)
    for _ in range(R):
        nxt = set()
        for v in frontier: nxt |= adj[v]
        nxt -= seen; seen |= nxt; frontier = nxt
    return seen

def dS_full(init_edeg, cur_edeg, touched):
    """Full shape action delta over the touched region, lam=1 units."""
    x = ESTAR - int(ESTAR)
    dS = 0.0; changed = set()
    for e in touched:
        a, b = init_edeg.get(e, 0), cur_edeg.get(e, 0)
        if a == b: continue
        changed |= set(e)
        if a: dS -= (ESTAR / 6.0) * ((a - ESTAR) ** 2 - x * (1 - x))
        if b: dS += (ESTAR / 6.0) * ((b - ESTAR) ** 2 - x * (1 - x))
    if not changed:
        return 0.0
    def counters(ed):
        n6 = defaultdict(int); mm = defaultdict(int)
        for e, d in ed.items():
            for v in e:
                if v in changed:
                    if d >= 6: n6[v] += 1
                    if not _legal(d): mm[v] += 1
        return n6, mm
    n60, mm0 = counters(init_edeg); n61, mm1 = counters(cur_edeg)
    for v in changed:
        dS += ZLEG * (wm.U_zleg(n61.get(v, 0)) - wm.U_zleg(n60.get(v, 0)))
        dS += CIMP * (mm1.get(v, 0) ** 2 - mm0.get(v, 0) ** 2)
    return dS

def search(P, core, tip_edge, label):
    init_edeg = dict(P.edeg)
    init_tets = len(P.tets)
    def ill_core():
        return frozenset(e for e, d in P.edeg.items()
                         if not _legal(d) and e[0] in core and e[1] in core)
    init_ill = ill_core()
    mode_goal = os.environ.get("GOAL", "tip")
    goal_ill = (frozenset() if mode_goal == "empty"
                else init_ill - {tip_edge})
    # translate mode: forbidden landing zone = the hinge-flip octahedron
    # (tip edge + its 4-cycle link) -- landing there is reachable by plain
    # 4->4 flips and is not transport.
    octa = set(tip_edge)
    for e, d in P.edeg.items():
        pass
    if mode_goal == "translate":
        for t in P.tets:
            if set(tip_edge) <= t:
                octa |= t
    print(f"\n=== {label}: tip edge {tip_edge} ===")
    print(f"  core |V|={len(core)}  patch tets={init_tets}  "
          f"illegal-in-core={len(init_ill)} -> goal {len(goal_ill)}", flush=True)

    sols = []
    visited = {}
    nodes = [0]
    t0 = time.time()

    def focus_verts():
        f = set(tip_edge)
        for e in P.touched:
            if P.edeg.get(e, 0) != init_edeg.get(e, 0):
                f |= set(e)
        return f

    def dfs(depth, seq):
        nodes[0] += 1
        if nodes[0] > NODE_CAP: return
        if nodes[0] % 500000 == 0:
            print(f"    ...{nodes[0]} nodes, {len(sols)} sols, "
                  f"{time.time()-t0:.0f}s", flush=True)
        cur_ill = ill_core()
        hit = False
        if seq and mode_goal == "translate":
            # every OTHER initial illegal edge must survive: deletions of
            # bundle partners are a different reaction (fusion), not
            # translation
            if tip_edge not in cur_ill and not (goal_ill - cur_ill):
                extra = cur_ill - goal_ill
                if len(extra) == 1:
                    f = next(iter(extra))
                    if P.edeg.get(f, 0) == 4 and not set(f) <= octa:
                        hit = True
        elif seq and cur_ill == goal_ill:
            hit = True
        if hit:
            ds = LAM * dS_full(init_edeg, P.edeg, frozenset(P.touched))
            land = (sorted(cur_ill - goal_ill)[0]
                    if mode_goal == "translate" else None)
            sols.append((list(seq), ds, len(P.tets) - init_tets, land))
            print(f"  [d{len(seq)}] GOAL  dS={ds:+.3f}  net_tets="
                  f"{len(P.tets)-init_tets:+d}  land={land}  seq={seq}",
                  flush=True)
            return                     # don't extend past a goal
        if depth == 0: return
        if len(cur_ill) > len(init_ill) + ILL_EXCESS: return
        sig = frozenset((e, P.edeg[e]) for e in P.touched
                        if P.edeg.get(e, 0) != init_edeg.get(e, 0))
        if visited.get(sig, -1) >= depth: return
        visited[sig] = depth
        for mm in P.moves(core, focus_verts()):
            if mm[0] == "23":
                _, tri, x, y = mm
                P.apply_23(tri, x, y)
                dfs(depth - 1, seq + [mm]); P.undo_23(tri, x, y)
            else:
                _, (u, w), a, b, c = mm
                P.apply_32(u, w, a, b, c)
                dfs(depth - 1, seq + [mm]); P.undo_32(u, w, a, b, c)
            if nodes[0] > NODE_CAP: return

    for d in range(2, MAX_DEPTH + 1):
        visited.clear()
        print(f"  depth {d} ... (nodes so far {nodes[0]})", flush=True)
        dfs(d, [])
        if sols:
            break                       # minimal length found; stop deepening
    print(f"  done: {nodes[0]} nodes, {time.time()-t0:.0f}s, "
          f"{len(sols)} solutions", flush=True)
    sols.sort(key=lambda s: s[1])
    for seq, ds, dn, land in sols[:12]:
        print(f"    dS={ds:+.3f} net_tets={dn:+d} len={len(seq)} "
              f"land={land}: {seq}", flush=True)
    return sols

def main():
    m = ddg.Manifold.load(SNAP, 3)
    pairs, degs = m.illegal_edges()
    ill = {tuple(sorted(int(x) for x in p)): int(d)
           for p, d in zip(pairs, degs)}
    F = [tuple(sorted(int(x) for x in t)) for t in np.asarray(m.facets())]
    d4 = [e for e, d in ill.items() if d == 4]
    adjv = defaultdict(list)
    for a, b in d4:
        adjv[a].append(b); adjv[b].append(a)
    val = {v: len(nb) for v, nb in adjv.items()}
    lone = [e for e in d4 if val[e[0]] == 1 and val[e[1]] == 1]
    line_tips = [e for e in d4
                 if (val[e[0]] == 1) != (val[e[1]] == 1)
                 and max(val[e[0]], val[e[1]]) == 2]
    print(f"{os.path.basename(SNAP)}: {len(d4)} deg-4 "
          f"({len(lone)} lone, {len(line_tips)} clean line tips)", flush=True)

    mode = os.environ.get("TARGETS", "first")
    targets = []
    if mode == "bundle":
        # deg-4 fragment-tip edges, biggest (longest-lived) complexes first
        import defect_state as dsm
        st = dsm.DefectState(m)
        comps = sorted(st.components(), key=lambda c: -len(c.sig))
        for cx in comps:
            if not cx.sig:
                continue
            cv = set(cx.verts)
            tips_c = [e for e in d4 if e[0] in cv and e[1] in cv
                      and (val[e[0]] == 1 or val[e[1]] == 1)]
            for e in tips_c:
                targets.append((f"B{len(cx.sig)}.{len(targets)}", e))
        targets = targets[:int(os.environ.get("NTARG", "25"))]
    elif mode == "sweep":
        targets += [(f"LONE[{i}]", e) for i, e in enumerate(lone)]
        targets += [(f"TIP[{i}]", e) for i, e in enumerate(line_tips)]
    elif mode.startswith("lone:"):
        targets = [(f"LONE[{mode[5:]}]", lone[int(mode[5:])])]
    else:
        if lone: targets.append(("LONE deg-4 edge", lone[0]))
        if line_tips: targets.append(("LINE tip", line_tips[0]))

    summary = []
    for label, e in targets:
        Vball = ball(F, e, R_BALL)
        Vcore = frozenset(ball(F, e, R_CORE))
        region = [t for t in F if set(t) <= Vball]
        P = Patch([frozenset(t) for t in region])
        sols = search(P, Vcore, e, label)
        best = min((s[1] for s in sols), default=None)
        summary.append((label, e, len(sols), best))
        print(f"SUMMARY {label} {e}: nsol={len(sols)}"
              + (f" best_dS={best:+.3f}" if best is not None else ""), flush=True)
    nfound = sum(1 for *_, n, _b in summary if n)
    print(f"\nSWEEP RESULT: {nfound}/{len(summary)} targets retractable "
          f"(depth<={MAX_DEPTH}, excess<={ILL_EXCESS})", flush=True)

if __name__ == "__main__":
    main()
