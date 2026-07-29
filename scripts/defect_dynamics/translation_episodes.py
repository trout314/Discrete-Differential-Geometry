"""Natural deg-4 translation episodes, v2 (fixes: same-event relocations were
never checked; the restored-condition was global so any unrelated flicker in
the window aborted the episode).

Causal-chain scan: from each deg-4 death, grow a support set from the death
move's vertex labels; only later moves TOUCHING the support belong to the
episode (unrelated activity elsewhere is ignored).  Episode completes when the
net change over its own records is exactly {e lost, one new deg-4 f gained,
all other touched edges restored}.  Classified by whether f escaped the death
move's support (for a 4->4 death the labels are the octahedron -- escape means
genuine transport out of the flip cavity)."""
import json, os
from collections import Counter

_ROOT = os.path.normpath(os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "..", ".."))
SRC = os.path.join(_ROOT, "data/rxn_lam035_m4/deg4_provenance.jsonl")
W = 40                                # max related records per episode
TYP = {0: "1->4", 1: "2->3", 2: "3->2", 3: "4->1", 4: "4->4"}

recs = [json.loads(l) for l in open(SRC)]
print(f"{len(recs)} crossing records, window {W} related records")

episodes = []
for i, c in enumerate(recs):
    for dx in c["cross"]:
        if not (dx["d0"] == 4 and dx["d1"] != 4):
            continue
        e = tuple(dx["e"])
        death_support = frozenset(c["lab"]) - {-1}
        support = set(death_support)
        net = {}
        types = []
        nrel = 0

        def upd(cr):
            for x in cr["cross"]:
                k = tuple(x["e"])
                if k not in net:
                    net[k] = [x["d0"], x["d1"]]
                else:
                    net[k][1] = x["d1"]

        def status():
            """(complete, f): net == {e lost, exactly one new deg-4 f}."""
            if net.get(e, [4, 4])[1] == 4:
                return "dead", None
            f = None
            for k, (a, b) in net.items():
                if k == e or a == b:
                    continue
                if b == 4 and a != 4:
                    if f is not None:
                        return "open", None
                    f = k
                elif b != a:
                    return "open", None
            return ("done", f) if f else ("open", None)

        upd(c)
        types.append(TYP[c["typ"]])
        st, f = status()
        j = i
        while st == "open" and nrel < W and j + 1 < len(recs):
            j += 1
            lab = set(recs[j]["lab"]) - {-1}
            if not (lab & support):
                continue                        # unrelated activity: ignore
            support |= lab
            nrel += 1
            upd(recs[j])
            types.append(TYP[recs[j]["typ"]])
            st, f = status()
        if st == "done":
            episodes.append(dict(
                sw=c["sw"], e=e, f=f, nmv=len(types), types=tuple(types),
                escaped=not set(f) <= death_support,
                share=len(set(f) & set(e))))

print(f"\ntranslation episodes: {len(episodes)}")
for tag, sel in [("IN-CAVITY (hinge-orbit)", [p for p in episodes if not p["escaped"]]),
                 ("ESCAPED (genuine transport)", [p for p in episodes if p["escaped"]])]:
    print(f"\n{tag}: {len(sel)}")
    if not sel:
        continue
    print(f"  length histogram: {dict(sorted(Counter(p['nmv'] for p in sel).items()))}")
    print(f"  move strings (top 10):")
    for t, n in Counter(p["types"] for p in sel).most_common(10):
        print(f"    {n:6d}  {'|'.join(t)}")
    print(f"  e-f vertex sharing: {dict(sorted(Counter(p['share'] for p in sel).items()))}")

esc = [p for p in episodes if p["escaped"]]
print("\nshortest escaped episodes:")
for p in sorted(esc, key=lambda p: p["nmv"])[:10]:
    print(f"  sw{p['sw']}: {p['e']} -> {p['f']}  {p['nmv']} moves "
          f"({'|'.join(p['types'])})  share={p['share']}")
