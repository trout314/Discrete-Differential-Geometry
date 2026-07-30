"""Maximalist flight recorder for sampler runs.

Every research script that drives :class:`ManifoldSampler` should go through
:func:`record_run` (simple case) or :class:`Recorder` (custom loops) instead
of bare ``sampler.run()``. Recording is UNCONDITIONAL and cheap; equilibrium
gating / certification is post-processing on the recorded series -- the two
concerns are deliberately decoupled (see memory: maximalist-recording).

Outputs, for a given ``out`` prefix:

  <out>.rec.jsonl   one header record, then one row per chunk:
      O(1) tier (every chunk): sweep, wall time, objective, f-vector,
        mean edge degree 6*f3/f1, per-move-type acceptance deltas, channel
        stats (slide / nonlocal / worm), D GC heap used/free (leak
        sentinel).
      census tier (every ``census_every`` chunks): n_ill, illegal degree
        histogram, defect components (count + top sizes) -- via the D-side
        O(E) ``illegal_edges`` dump.
  <out>.mid.mfd(+.cocycle.npz)    snapshot at half the requested sweeps
  <out>.final.mfd(+.cocycle.npz)  final state, written by ``finish()``

The header captures SamplerParams, potential couplings and channel
probabilities (as recorded by their setters), the RNG seed (if set through
``discrete_differential_geometry.set_random_seed``), argv, git commit, and
the initial f-vector.
"""
from __future__ import annotations

import dataclasses
import json
import os
import subprocess
import sys
import time
from collections import Counter

import numpy as np

from . import gc_stats
from . import cocycle as _coc


def _git_commit() -> str | None:
    try:
        return subprocess.run(
            ["git", "rev-parse", "--short", "HEAD"],
            capture_output=True, text=True, timeout=5,
            cwd=os.path.dirname(os.path.abspath(__file__)),
        ).stdout.strip() or None
    except Exception:
        return None


def _components(pairs: np.ndarray) -> list[int]:
    """Sizes (in edges, descending) of vertex-connected illegal components."""
    parent: dict = {}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    v2e: dict = {}
    for i, (a, b) in enumerate(map(tuple, pairs)):
        parent[i] = i
        for v in (a, b):
            if v in v2e:
                ra, rb = find(i), find(v2e[v])
                if ra != rb:
                    parent[ra] = rb
            v2e[v] = i
    if not parent:
        return []
    return sorted(Counter(find(i) for i in parent).values(), reverse=True)


class Recorder:
    """Chunked run loop with unconditional instrumentation.

    Use :meth:`run` for the common case; scripts with their own loop call
    :meth:`step` per chunk and :meth:`finish` at the end.
    """

    def __init__(self, sampler, out: str, chunk: int = 25,
                 census_every: int = 1, snap_mid: bool = True,
                 cocycle_box=None, extra_meta: dict | None = None):
        self.s = sampler
        self.out = out
        self.chunk = int(chunk)
        self.census_every = int(census_every)
        self.snap_mid = snap_mid
        self.cocycle_box = cocycle_box
        self.sw = 0
        self.nchunk = 0
        self._t0 = time.time()
        self._mid_done = not snap_mid
        self._target = None
        os.makedirs(os.path.dirname(os.path.abspath(out)) or ".",
                    exist_ok=True)
        self._fh = open(out + ".rec.jsonl", "w")
        self._prev = self._counters()
        self._write_header(extra_meta or {})

    # -- internals ---------------------------------------------------------
    def _counters(self) -> dict:
        st = self.s.get_stats()
        d = {
            "tried": int(st.total_tried),
            "accepted": int(st.total_accepted),
            "hinge_tries": int(st.hinge_tries),
            "hinge_accepts": int(st.hinge_accepts),
            "bis_tries": [int(x) for x in st.bistellar_tries],
            "bis_accepts": [int(x) for x in st.bistellar_accepts],
        }
        for name in ("slide_stats", "nonlocal_slide_stats", "worm_stats"):
            try:
                d[name] = [int(x) for x in getattr(self.s, name)()]
            except Exception:
                d[name] = None
        return d

    def _write(self, rec: dict) -> None:
        self._fh.write(json.dumps(rec) + "\n")
        self._fh.flush()

    def _write_header(self, extra: dict) -> None:
        import __main__
        p = getattr(self.s, "_params", None)
        params = dataclasses.asdict(p) if dataclasses.is_dataclass(p) else None
        v = self.s.manifold
        fv = [int(x) for x in v.f_vector]
        self._write({
            "kind": "header",
            "time": time.strftime("%Y-%m-%d %H:%M:%S"),
            "argv": sys.argv,
            "main": getattr(__main__, "__file__", None),
            "git": _git_commit(),
            "params": params,
            "n6_potential": getattr(self.s, "_recorded_n6_potential", None),
            "slide_prob": getattr(self.s, "_recorded_slide_prob", None),
            "nonlocal_slide": getattr(self.s, "_recorded_nonlocal_slide",
                                      None),
            "worm_prob": getattr(self.s, "_recorded_worm_prob", None),
            "seed": self._last_seed(),
            "chunk": self.chunk,
            "census_every": self.census_every,
            "f_vector0": fv,
            "objective0": float(self.s.current_objective),
            **extra,
        })

    @staticmethod
    def _last_seed():
        from . import _last_random_seed
        return _last_random_seed

    def _save_state(self, tag: str) -> str:
        path = f"{self.out}.{tag}.mfd"
        self.s.manifold.save(path)
        try:
            e, om = self.s.read_cocycle()
            if len(e):
                meta = {"sweeps": self.sw}
                if self.cocycle_box is not None:
                    meta["box"] = self.cocycle_box
                _coc.save_cocycle(f"{self.out}.{tag}.cocycle.npz", e, om,
                                  **meta)
        except Exception:
            pass                      # no cocycle attached
        return path

    # -- public ------------------------------------------------------------
    def step(self, sweeps: float | None = None) -> dict:
        """Run one chunk and append its row. Returns the row."""
        n = self.chunk if sweeps is None else sweeps
        self.s.run(sweeps=n)
        self.sw += n
        self.nchunk += 1
        cur = self._counters()
        prev, self._prev = self._prev, cur

        def delta(key):
            a, b = cur.get(key), prev.get(key)
            if a is None or b is None:
                return None
            if isinstance(a, list):
                return [x - y for x, y in zip(a, b)]
            return a - b

        v = self.s.manifold
        fv = [int(x) for x in v.f_vector]
        used, free = gc_stats()
        row = {
            "kind": "row", "sw": self.sw,
            "t": round(time.time() - self._t0, 2),
            "obj": float(self.s.current_objective),
            "f": fv,
            "e_mean": 6.0 * fv[3] / fv[1] if len(fv) > 3 and fv[1] else None,
            "d_tried": delta("tried"), "d_accepted": delta("accepted"),
            "d_hinge": [delta("hinge_tries"), delta("hinge_accepts")],
            "d_bis_tries": delta("bis_tries"),
            "d_bis_accepts": delta("bis_accepts"),
            "slide": cur["slide_stats"],
            "nonlocal": cur["nonlocal_slide_stats"],
            "worm": cur["worm_stats"],
            "gc_used_mb": round(used / 1e6, 1),
            "gc_free_mb": round(free / 1e6, 1),
        }
        if self.census_every and self.nchunk % self.census_every == 0:
            pairs, degs = v.illegal_edges()
            comps = _components(pairs)
            row["n_ill"] = int(len(degs))
            row["ill_hist"] = {int(k): int(c)
                               for k, c in sorted(Counter(degs).items())}
            row["ncomp"] = len(comps)
            row["top_comps"] = comps[:5]
        self._write(row)
        if (not self._mid_done and self._target
                and self.sw >= self._target / 2):
            self._mid_done = True
            self._write({"kind": "snapshot", "sw": self.sw,
                         "path": self._save_state("mid")})
        return row

    def run(self, sweeps: float) -> None:
        """Chunked run to ``sweeps`` total (mid snapshot halfway)."""
        self._target = sweeps
        while self.sw < sweeps:
            self.step(min(self.chunk, sweeps - self.sw))

    def finish(self) -> str:
        """Final snapshot + closing record. Returns the final .mfd path."""
        path = self._save_state("final")
        self._write({"kind": "final", "sw": self.sw,
                     "t": round(time.time() - self._t0, 2),
                     "path": path, "counters": self._prev})
        self._fh.close()
        return path


def record_run(sampler, sweeps: float, out: str, **kw) -> str:
    """Run ``sweeps`` with full instrumentation; returns final .mfd path.

    Drop-in replacement for ``sampler.run(sweeps=...)`` in research
    scripts. Keyword args are forwarded to :class:`Recorder`.
    """
    rec = Recorder(sampler, out, **kw)
    rec.run(sweeps)
    return rec.finish()
