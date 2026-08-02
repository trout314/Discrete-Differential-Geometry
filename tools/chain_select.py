"""Shared BC-chain class selection and provenance.

A BC chain is the closed sliding-window walk of ``worm_helix.bc_orbit``. The
frames of a triangulation partition into cycles -- the complete chain set --
and ``Aut(K)`` acts on it, so chains fall into finitely many CLASSES that are
properties of the crystal, not of the supercell. R has 14 of them, with
lengths spanning 99 to 2439 at m=3.

WHY THIS EXISTS. Nearly every driver in ``scripts/defect_dynamics/`` starts
from ``bc_orbit(m, F[0])`` -- the chain through tet 0, in whatever vertex order
that tet happens to be stored in. That picks one of the 14 classes by accident,
does not record which, and cannot be reproduced or compared across runs. This
module makes the class an explicit, recorded input:

    cc = ChainClasses(ref_path)
    k  = cc.select("axis")            # or an index, or "longest"
    seq = cc.vertices(k)              # same convention as bc_orbit
    meta = cc.provenance(k)           # stamp this in run metadata

Positions (hence windings) are only available for crystals built by
``tcp_reference``; everything else degrades to combinatorial selection, and
winding-based selectors RAISE rather than quietly fall back.
"""
from __future__ import annotations

import os
import re

import numpy as np

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

__all__ = ["ChainClasses", "parse_ref_name", "reference_positions",
           "certify_positions", "chain_winding", "minimg", "add_chain_args",
           "resolve_chain"]


def parse_ref_name(path):
    """(structure, m) from a tcp_reference filename, else (None, None)."""
    mm = re.match(r"T3_([A-Z0-9]+)_m(\d+)_N\d+\.mfd", os.path.basename(path))
    return (mm.group(1).lower(), int(mm.group(2))) if mm else (None, None)


def reference_positions(struct, m, nvert):
    """Fractional torus coordinates (period m) per vertex id, or None.

    Reproduces ``tcp_reference``'s site perturbation and id scheme exactly
    (``v = cell*ns + site``); see ``cocycle_check.reference_frac_positions``.
    """
    import sys
    sys.path.insert(0, os.path.join(_ROOT, "scripts"))
    from tcp_reference import STRUCTURES
    if struct not in STRUCTURES:
        return None
    _, sites, _, _ = STRUCTURES[struct]
    ns = len(sites)
    if ns * m ** 3 != nvert:
        return None
    rng = np.random.default_rng(12345)
    sites = sites + 1e-6 * rng.standard_normal(sites.shape)
    v = np.arange(nvert)
    c = v // ns
    return (sites[v % ns] + np.stack([c // (m * m), (c // m) % m, c % m],
                                     axis=1)) % m


def minimg(d, m):
    return d - m * np.round(d / m)


def certify_positions(rp, sym, m, tol=1e-9):
    """Certify that minimum-image reconstructs the TRUE edge displacements.

    Stored positions give a vertex's spot in the fundamental domain, not which
    periodic image an edge runs to, so every displacement is RECONSTRUCTED as
    the minimum image. That is right only if the true displacement is already
    inside (-m/2, m/2) componentwise; otherwise the reconstruction is off by
    exactly m on that edge and the winding is off by exactly +-1 -- still a
    tidy integer vector, just the wrong one, and integrality cannot catch it.

    "The minimum image is under m/2" is NOT a test: minimg returns components
    in that range by construction. The certificate is CLOSURE -- the
    reconstruction is a 1-cochain and the truth is closed, so requiring the
    minimum images to sum to zero around every triangle pins it up to a
    cocycle, and one mis-wrapped edge breaks closure on every triangle
    containing it. The wrap MARGIN is also required (> 2x), which rules out
    the residual loophole of a coherent mis-wrapped cut surface.

    Returns (worst edge component, closure residual, margin, ok).
    """
    view = sym.view
    ext = view.labels.astype(np.int64)
    eu = np.array([[ext[u], ext[w]] for (u, w) in view.edges])
    worst = float(np.abs(minimg(rp[eu[:, 1]] - rp[eu[:, 0]], m)).max())
    fa = np.array([[ext[a], ext[b], ext[c]] for (a, b, c) in view.faces])
    cyc = (minimg(rp[fa[:, 1]] - rp[fa[:, 0]], m)
           + minimg(rp[fa[:, 2]] - rp[fa[:, 1]], m)
           + minimg(rp[fa[:, 0]] - rp[fa[:, 2]], m))
    closure = float(np.abs(cyc).max())
    margin = (m / 2.0) / worst if worst > 0 else float("inf")
    return worst, closure, margin, (closure < tol and margin > 2.0)


def chain_winding(sym, chain, rp, m):
    """Winding of ONE chain, in supercell periods (its class in H_1(T^3))."""
    view, lab = sym.view, sym.view.labels
    seq = [int(lab[view.frame_window(int(f))[0]]) for f in chain]
    p = rp[seq]
    return minimg(np.roll(p, -1, axis=0) - p, m).sum(axis=0) / m


class ChainClasses:
    """The BC chain classes of one reference crystal, selectable and citable.

    ``positions`` is optional: without it (a non-tcp_reference file, or a
    melt) the combinatorial selectors still work and the winding-based ones
    raise. Nothing silently falls back to an arbitrary chain.
    """

    def __init__(self, path, positions=True):
        import sys
        sys.path.insert(0, os.path.join(_ROOT, "python"))
        from discrete_differential_geometry import CrystalSymmetry
        self.path = path
        self.sym = CrystalSymmetry.for_manifold_path(path)
        self.view = self.sym.view
        _, self.chain_of, self.chains = self.view.chain_tables()
        self._oid, self.members, self.reps = self.sym.chain_orbits()
        self.struct, self.m = parse_ref_name(path)
        self.rp = None
        self.position_note = "not requested"
        if positions and self.struct:
            rp = reference_positions(self.struct, self.m, self.view.V)
            if rp is None:
                self.position_note = (
                    f"no tcp_reference positions for {self.struct} m={self.m}")
            else:
                worst, closure, margin, ok = certify_positions(rp, self.sym,
                                                               self.m)
                if ok:
                    self.rp = rp
                    self.position_note = (
                        f"certified (closure {closure:.1e}, margin "
                        f"{margin:.1f}x)")
                else:
                    self.position_note = (
                        f"REJECTED (closure {closure:.1e}, margin "
                        f"{margin:.1f}x) -- winding unavailable")
        elif positions:
            self.position_note = "filename is not a tcp_reference crystal"

    # -- basics --------------------------------------------------------------

    @property
    def n_classes(self):
        return len(self.members)

    def length(self, k):
        return len(self.chains[self.reps[k]])

    def frame(self, k):
        """Seed window (ordered vertex 4-tuple, original labels) of class k.

        ``bc_orbit(manifold, frame(k))`` reproduces ``vertices(k)`` exactly.
        """
        lab = self.view.labels
        w = self.view.frame_window(int(self.chains[self.reps[k]][0]))
        return tuple(int(lab[x]) for x in w)

    def vertices(self, k):
        """Cyclic vertex sequence of class k's representative, in the same
        convention as ``worm_helix.bc_orbit`` (window[0] per step)."""
        lab, view = self.view.labels, self.view
        return [int(lab[view.frame_window(int(f))[0]])
                for f in self.chains[self.reps[k]]]

    def class_of_frame(self, window):
        """Which class an ordered vertex 4-tuple lies on -- e.g. to record the
        class a legacy ``bc_orbit(m, F[0])`` seed actually landed in."""
        w = tuple(self.view.index[int(x)] for x in window)
        return int(self._oid[int(self.chain_of[self.view.frame_id(w)])])

    # -- winding -------------------------------------------------------------

    def _require_positions(self, what):
        if self.rp is None:
            raise ValueError(
                f"{what} needs vertex positions, which are unavailable for "
                f"{os.path.basename(self.path)}: {self.position_note}")

    def windings(self, k):
        """Distinct winding vectors realized across class k.

        A class does not have "a" winding: lattice translations preserve it
        but the point group ROTATES it, so this is the orbit, not one
        representative's vector.
        """
        self._require_positions("winding")
        out = set()
        for i in self.members[k]:
            w = chain_winding(self.sym, self.chains[i], self.rp, self.m)
            out.add(tuple(int(x) for x in np.round(w)))
        return sorted(out)

    def representative_winding(self, k):
        """Winding of the exact chain ``vertices(k)`` returns (not merely one
        of the class's windings -- the point group rotates them)."""
        self._require_positions("winding")
        w = chain_winding(self.sym, self.chains[self.reps[k]], self.rp, self.m)
        return np.round(w).astype(int)

    def pure_axis_classes(self):
        """Classes some of whose chains wind along a single lattice axis."""
        self._require_positions("axis selection")
        return [k for k in range(self.n_classes)
                if any(sum(1 for x in w if x) == 1 for w in self.windings(k))]

    # -- selection -----------------------------------------------------------

    def select(self, selector):
        """Resolve a selector to a class index.

        int / numeric string : that class index
        "axis"               : the SHORTEST pure-axis class (raises if none)
        "longest"/"shortest" : by chain length
        "w=a,b,c"            : the class realizing that winding vector
        """
        if isinstance(selector, (int, np.integer)):
            k = int(selector)
        elif isinstance(selector, str) and selector.strip().lstrip("-").isdigit():
            k = int(selector)
        elif selector == "axis":
            cands = self.pure_axis_classes()
            if not cands:
                raise ValueError(
                    f"{os.path.basename(self.path)} has no pure-axis BC chain "
                    f"class (checked all {self.n_classes} exactly)")
            return min(cands, key=self.length)
        elif selector == "longest":
            return max(range(self.n_classes), key=self.length)
        elif selector == "shortest":
            return min(range(self.n_classes), key=self.length)
        elif isinstance(selector, str) and selector.startswith("w="):
            want = tuple(int(x) for x in selector[2:].split(","))
            hits = [k for k in range(self.n_classes) if want in self.windings(k)]
            if len(hits) != 1:
                raise ValueError(f"winding {want} matches {len(hits)} classes "
                                 f"(need exactly 1): {hits}")
            return hits[0]
        else:
            raise ValueError(f"unrecognised chain selector {selector!r}")
        if not 0 <= k < self.n_classes:
            raise ValueError(f"chain class {k} out of range "
                             f"(0..{self.n_classes - 1})")
        return k

    # -- provenance ----------------------------------------------------------

    def provenance(self, k):
        """Metadata block to stamp into run output. A transport number without
        this cannot be compared to another run."""
        out = dict(crystal=os.path.basename(self.path),
                   aut_order=self.sym.order,
                   n_chain_classes=self.n_classes,
                   chain_class=int(k),
                   chain_length=self.length(k),
                   chains_in_class=len(self.members[k]),
                   seed_frame=list(self.frame(k)),
                   positions=self.position_note)
        if self.rp is not None:
            out["windings"] = [list(w) for w in self.windings(k)]
        return out

    def summary_line(self, k):
        w = ""
        if self.rp is not None:
            ws = self.windings(k)
            w = ("  winding " + " ".join("(" + ",".join(f"{x:+d}" for x in v)
                                         + ")" for v in ws[:2])
                 + (" ..." if len(ws) > 2 else ""))
        return (f"chain class {k}/{self.n_classes}: L={self.length(k)}, "
                f"{len(self.members[k])} chains, seed frame {self.frame(k)}{w}")


def add_chain_args(ap, default="axis"):
    """Add the standard ``--chain-class`` selector to an argparse parser."""
    ap.add_argument("--chain-class", default=default,
                    help="BC chain class: an index, 'axis', 'longest', "
                         "'shortest', or 'w=a,b,c'. Recorded in output "
                         "metadata (default: %(default)s)")
    return ap


def resolve_chain(path, selector, verbose=True):
    """(ChainClasses, class index, vertex sequence, provenance)."""
    cc = ChainClasses(path)
    k = cc.select(selector)
    if verbose:
        print(f"  {cc.summary_line(k)}")
    return cc, k, cc.vertices(k), cc.provenance(k)
