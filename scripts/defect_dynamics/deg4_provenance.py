"""Watch the ACTUAL dynamics create/destroy deg-4 edges: run the production
sampler on the m4 R gas with the accepted-move event log on, replay every move
against an incrementally-maintained tet/degree state (audited vs the manifold),
and record every degree-4 boundary crossing with its provenance window.

Answers: which moves birth a LONE deg-4 edge in the wild, which kill one, and
how long the causal episodes are -- the real mechanism, read off the record.
"""
import os, sys, json, itertools, time
from collections import Counter, defaultdict
import numpy as np

_ROOT = os.path.normpath(os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "..", ".."))
sys.path.insert(0, os.path.join(_ROOT, "python"))
import discrete_differential_geometry as ddg

SNAP = os.path.join(_ROOT, "data/mgas/lam35r_snap15000.mfd")
LAM, ET, ZLEG, CIMP = 0.35, 5.105025, 0.6, 1.0
SWEEPS = int(os.environ.get("SWEEPS", "1000"))
AUDIT_EVERY = 100
OUT = os.path.join(_ROOT, "data/rxn_lam035_m4/deg4_provenance")

def E(a, b): return (a, b) if a < b else (b, a)

def main():
    ddg.set_random_seed(20260729)
    m = ddg.Manifold.load(SNAP, 3)
    p = ddg.SamplerParams(
        num_facets_target=m.num_facets, hinge_degree_target=ET,
        num_facets_coef=0.1, num_hinges_coef=0.0,
        hinge_degree_variance_coef=0.0, codim3_degree_variance_coef=0.0,
        hinge_degree_target_coef=LAM * ET / 6.0)
    s = ddg.ManifoldSampler(m, p)
    s.set_n6_potential(ZLEG * LAM, CIMP * LAM, tilt=[0.0] * 5)
    v = s.manifold
    s.enable_event_log(64.0)
    s.drain_event_log()

    tets = {frozenset(int(x) for x in t) for t in np.asarray(v.facets())}
    edeg = Counter()
    for t in tets:
        for e in itertools.combinations(sorted(t), 2):
            edeg[e] += 1
    illv = Counter()                      # vertex -> #illegal edges at it
    for e, d in edeg.items():
        if d not in (5, 6):
            illv[e[0]] += 1; illv[e[1]] += 1

    events = []          # (clock, typ, labels) for every accepted move
    vlast = {}           # vertex -> last event index touching it
    crossings = []       # dicts for every event with a deg-4/legality crossing
    typc = Counter()

    def rmadd(typ, lab):
        """removed/added tets for an event record."""
        if typ == 1:      # 2->3: center face, coCenter apexes
            a, b, c, x, y = lab[:5]
            rm = [(a, b, c, x), (a, b, c, y)]
            ad = [(x, y, a, b), (x, y, b, c), (x, y, a, c)]
        elif typ == 2:    # 3->2: center edge, coCenter link
            u, w, a, b, c = lab[:5]
            rm = [(u, w, a, b), (u, w, b, c), (u, w, a, c)]
            ad = [(a, b, c, u), (a, b, c, w)]
        elif typ == 0:    # 1->4: center facet, coCenter new vertex
            a, b, c, d, x = lab[:5]
            rm = [(a, b, c, d)]
            ad = [(a, b, c, x), (a, b, d, x), (a, c, d, x), (b, c, d, x)]
        elif typ == 3:    # 4->1: center vertex, coCenter 4 neighbours
            x, a, b, c, d = lab[:5]
            rm = [(a, b, c, x), (a, b, d, x), (a, c, d, x), (b, c, d, x)]
            ad = [(a, b, c, d)]
        else:             # 4->4: removedEdge + linkCycle; diagonal (c0,c2)
            u, w, c0, c1, c2, c3 = lab
            rm = [(u, w, c0, c1), (u, w, c1, c2), (u, w, c2, c3),
                  (u, w, c3, c0)]
            ad = [(c0, c2, u, c1), (c0, c2, c1, w), (c0, c2, w, c3),
                  (c0, c2, c3, u)]
        return rm, ad

    def apply_event(idx, clock, typ, lab):
        rm, ad = rmadd(typ, lab)
        touched = set()
        old = {}
        for t in rm:
            tt = frozenset(t)
            if tt not in tets:
                raise RuntimeError(f"event {idx} type {typ}: removed tet "
                                   f"missing {t}")
            tets.discard(tt)
            for e in itertools.combinations(sorted(t), 2):
                if e not in old: old[e] = edeg[e]
                edeg[e] -= 1; touched.add(e)
        for t in ad:
            tets.add(frozenset(t))
            for e in itertools.combinations(sorted(t), 2):
                if e not in old: old[e] = edeg.get(e, 0)
                edeg[e] += 1; touched.add(e)
        cross = []
        for e in touched:
            d0, d1 = old[e], edeg.get(e, 0)
            if d1 == 0: del edeg[e]
            if d0 == d1: continue
            leg0, leg1 = d0 in (5, 6), d1 in (5, 6)
            if (not leg0) != (not leg1) or 4 in (d0, d1) or d0 == 0 or d1 == 0:
                cross.append((e, d0, d1))
            if not leg0:
                illv[e[0]] -= 1; illv[e[1]] -= 1
            if not leg1:
                illv[e[0]] += 1; illv[e[1]] += 1
        return cross

    t0 = time.time()
    audit_fail = 0
    for sw in range(1, SWEEPS + 1):
        s.run(sweeps=1)
        ev = s.drain_event_log()
        if s.event_log_overflowed():
            print(f"sw{sw}: EVENT LOG OVERFLOWED", flush=True)
        for rec in ev:
            idx = len(events)
            clock, typ = int(rec["clock"]), int(rec["type"])
            lab = [int(x) for x in rec["labels"]]
            typc[typ] += 1
            cross = apply_event(idx, clock, typ, lab)
            events.append((clock, typ, tuple(lab)))
            if cross:
                # loneliness of each crossing deg-4 edge, AFTER the event
                info = []
                for e, d0, d1 in cross:
                    lone = (illv[e[0]] + illv[e[1]] == 2
                            and edeg.get(e, 0) == d1 and d1 == 4) if d1 == 4 \
                        else None
                    prior = max((vlast.get(x, -1) for x in e), default=-1)
                    info.append({"e": list(e), "d0": d0, "d1": d1,
                                 "lone": lone, "prior_touch": prior})
                crossings.append({"i": idx, "sw": sw, "clock": clock,
                                  "typ": typ, "lab": lab, "cross": info})
            rm_, ad_ = rmadd(typ, lab)
            for x in set(v for tt in rm_ + ad_ for v in tt):
                vlast[x] = idx
        if sw % AUDIT_EVERY == 0:
            live = {frozenset(int(x) for x in t)
                    for t in np.asarray(v.facets())}
            if live != tets:
                audit_fail += 1
                print(f"sw{sw}: AUDIT FAIL diff {len(live ^ tets)} "
                      f"(type-4 diagonal convention?)", flush=True)
                if audit_fail >= 2:
                    sys.exit(1)
                tets.clear(); tets.update(live)
            n4 = sum(1 for d in edeg.values() if d == 4)
            print(f"sw{sw} ({time.time()-t0:.0f}s): {len(events)} events, "
                  f"{len(crossings)} crossings, n4={n4}, "
                  f"types {dict(typc)}", flush=True)

    with open(OUT + ".jsonl", "w") as fh:
        for c in crossings:
            fh.write(json.dumps(c) + "\n")
    np.save(OUT + "_events.npy",
            np.array([(c, t, *l) for c, t, l in events], dtype=np.int64))
    print(f"\nwrote {OUT}.jsonl ({len(crossings)} crossings) and "
          f"_events.npy ({len(events)} events)", flush=True)

    # ---- inline mechanism summary -------------------------------------
    TYP = {0: "1->4", 1: "2->3", 2: "3->2", 3: "4->1", 4: "4->4"}
    born = [c for c in crossings
            for x in c["cross"] if x["d1"] == 4 and x["lone"]]
    died = [c for c in crossings
            for x in c["cross"] if x["d0"] == 4 and x["d1"] != 4]
    print(f"\naccepted moves: {len(events)}  by type "
          f"{ {TYP[k]: n for k, n in sorted(typc.items())} }")
    print(f"LONE deg-4 births: {len(born)}  by move type "
          f"{Counter(TYP[c['typ']] for c in born)}")
    print(f"deg-4 exits (any): {len(died)}  by move type "
          f"{Counter(TYP[c['typ']] for c in died)}")
    ex = Counter()
    for c in crossings:
        for x in c["cross"]:
            if x["d0"] == 4 and x["d1"] != 4:
                ex[x["d1"]] += 1
    print(f"deg-4 exit destinations (new degree): {dict(sorted(ex.items()))}")

if __name__ == "__main__":
    main()
