#!/usr/bin/env python3
"""M0: exact symbolic verification of the FPKMC chain formulas (SymPy).

Verifies, before any implementation (notes/FPKMC_DESIGN.md §3):

  A  the Metropolis min-identity used to telescope rate ratios;
  B  detailed balance / stationarity pi_j ∝ e^{−λ S_j} given slot
     symmetry I1 (edge multiplicities n_k shared by both directions);
  C  the Kolmogorov cycle condition on a closed orbit (lone-knot case:
     no absorbing sites, stationary distribution exists and is Boltzmann);
  D  the splitting-probability closed form
         P_R(s) = Σ_{k<s} e^{λ M_k}/n_k  /  Σ_{k<n} e^{λ M_k}/n_k ,
     M_k = max(S_k, S_{k+1}), against the exact harmonic solve;
  E  discrete-time invariance: per-attempt probabilities with arbitrary
     holding probabilities give the SAME splitting probabilities (the
     harmonic equation is scale-invariant row by row);
  F  conditional mean exit times: the w-equation construction
     (Q w = −h, w = E[T·1_R], E[T|R] = w/h) against the fundamental
     matrix of an exact rational instance;
  G  the phase-type tail P(T > t) = e_s^T P_int^t 1 against brute
     path enumeration on a small rational instance.

Everything is exact: symbolic where feasible (A–E), exact rationals where
symbolic case-splitting would explode (F, G). Any FAIL aborts with the
offending expression. Run time: seconds.
"""
import sys

import sympy as sp


def hdr(s):
    print(f"\n=== {s}")


ok_all = True


def check(name, expr_zero):
    global ok_all
    z = sp.simplify(expr_zero)
    good = z == 0 or z.equals(0)
    print(f"  {'PASS' if good else 'FAIL'}  {name}"
          + ("" if good else f"   residual: {z}"))
    if not good:
        ok_all = False


# ---------------------------------------------------------------- A
hdr("A: min(1, e^{-lam*d}) == e^{-lam*(max(d,0))}, both branches")
lam, d = sp.symbols("lam d", positive=True)
# branch d >= 0: min = e^{-lam d}, max(d,0) = d
check("d >= 0", sp.exp(-lam * d) - sp.exp(-lam * sp.Max(d, 0))
      .rewrite(sp.Piecewise).subs(sp.Max(d, 0), d))
# branch d <= 0 (substitute d -> -d, d >= 0): min = 1, max(-d,0) = 0
check("d <= 0", 1 - sp.exp(-lam * sp.Max(-d, 0))
      .subs(sp.Max(-d, 0), 0))
# Consequence used throughout: with M = max(a, b),
#   min(1, e^{-lam(a-b)}) = e^{-lam(M-b)}   (verified per branch)
a, b = sp.symbols("a b", real=True)
for cond, M in (("a>=b", a), ("a<=b", b)):
    dd = a - b
    lhs = sp.exp(-lam * dd) if cond == "a>=b" else sp.Integer(1)
    check(f"min-form, {cond}", lhs - sp.exp(-lam * (M - b)))

# ---------------------------------------------------------------- B
hdr("B: detailed balance with pi_j = e^{-lam S_j} (needs I1)")
# rates in M-form: q(j -> j+1) = nu * n_j * e^{-lam (M_j - S_j)}
#                  q(j+1 -> j) = nu * n_j * e^{-lam (M_j - S_{j+1})}
# (same edge multiplicity n_j both ways: THIS is invariant I1)
nu = sp.symbols("nu", positive=True)
Sj, Sj1, Mj, nj = sp.symbols("S_j S_j1 M_j n_j", real=True)
piL = sp.exp(-lam * Sj)
piR = sp.exp(-lam * Sj1)
qf = nu * nj * sp.exp(-lam * (Mj - Sj))
qb = nu * nj * sp.exp(-lam * (Mj - Sj1))
check("pi_j q(j->j+1) == pi_{j+1} q(j+1->j)", piL * qf - piR * qb)
print("  NOTE: without I1 (different multiplicities each way) the "
      "stationary law gains multiplicity factors -- V1 must verify I1.")

# ---------------------------------------------------------------- C
hdr("C: Kolmogorov cycle condition on a closed orbit (L=5)")
L = 5
S = sp.symbols(f"S0:{L}", real=True)
n = sp.symbols(f"n0:{L}", positive=True)
M = [sp.Symbol(f"M{k}", real=True) for k in range(L)]   # M_k on edge k,k+1
prod = sp.Integer(1)
for k in range(L):
    kp = (k + 1) % L
    qf = n[k] * sp.exp(-lam * (M[k] - S[k]))
    qb = n[k] * sp.exp(-lam * (M[k] - S[kp]))
    prod *= qf / qb
check("prod of q+/q- around the cycle == 1", prod - 1)

# ---------------------------------------------------------------- D
hdr("D: splitting probability closed form vs exact harmonic solve (n=6)")
N = 6                       # sites 0..5; 0 and 5 absorbing
S = sp.symbols(f"S0:{N}", real=True)
n = sp.symbols(f"n0:{N - 1}", positive=True)
M = [sp.Symbol(f"M{k}", real=True) for k in range(N - 1)]
qf = {j: n[j] * sp.exp(-lam * (M[j] - S[j])) for j in range(N - 1)}
qb = {j: n[j - 1] * sp.exp(-lam * (M[j - 1] - S[j])) for j in range(1, N)}
h = {0: sp.Integer(0), N - 1: sp.Integer(1)}
hs = sp.symbols(f"h1:{N - 1}")           # unknowns h_1 .. h_4
hh = {0: h[0], N - 1: h[N - 1]}
for j in range(1, N - 1):
    hh[j] = hs[j - 1]
eqs = [sp.Eq(qf[j] * (hh[j + 1] - hh[j]) + qb[j] * (hh[j - 1] - hh[j]), 0)
       for j in range(1, N - 1)]
sol = sp.solve(eqs, hs, dict=True)[0]
phi = [sp.exp(lam * M[k]) / n[k] for k in range(N - 1)]
Z = sum(phi)
for s_ in range(1, N - 1):
    closed = sum(phi[:s_]) / Z
    check(f"P_R({s_}) closed form", sp.simplify(sol[hs[s_ - 1]] - closed))

# ---------------------------------------------------------------- E
hdr("E: discrete-time invariance (holding probabilities drop out)")
# per-attempt chain: p_j^+ = c_j * qf_j, p_j^- = c_j * qb_j with ARBITRARY
# positive per-site scalings c_j (covers holding probability bookkeeping):
# the harmonic equation p^+ (h_{j+1}-h_j) + p^- (h_{j-1}-h_j) = 0 is
# row-scale invariant, so splitting probabilities are unchanged.
c = sp.symbols(f"c1:{N - 1}", positive=True)
eqs2 = [sp.Eq(c[j - 1] * qf[j] * (hh[j + 1] - hh[j])
              + c[j - 1] * qb[j] * (hh[j - 1] - hh[j]), 0)
        for j in range(1, N - 1)]
sol2 = sp.solve(eqs2, hs, dict=True)[0]
for s_ in range(1, N - 1):
    check(f"P_R({s_}) with per-site scalings",
          sp.simplify(sol2[hs[s_ - 1]] - sol[hs[s_ - 1]]))
print("  NOTE: exit TIMES do scale with c_j -- only probabilities are "
      "invariant; time bookkeeping must use the true per-attempt rates.")

# ---------------------------------------------------------------- F
hdr("F: conditional mean exit time, w-equation vs fundamental matrix")
# exact rational instance, sites 0..4, absorbing 0 and 4
pf = [None, sp.Rational(3, 10), sp.Rational(1, 5), sp.Rational(2, 5)]
pb = [None, sp.Rational(1, 10), sp.Rational(3, 10), sp.Rational(1, 10)]
Q = sp.Matrix([[-(pf[1] + pb[1]), pf[1], 0],
               [pb[2], -(pf[2] + pb[2]), pf[2]],
               [0, pb[3], -(pf[3] + pb[3])]])
vR = sp.Matrix([0, 0, pf[3]])           # absorption into right
vL = sp.Matrix([pb[1], 0, 0])
hvec = -Q.inv() * vR                    # splitting probabilities
w = -Q.inv() * hvec                     # w = E[T 1_R] (CTMC, unit rates)
# fundamental-matrix cross-check via the embedded DTMC with holding:
# P_int over interior with self-loops 1 - pf - pb
Pi = sp.Matrix([[1 - pf[1] - pb[1], pf[1], 0],
                [pb[2], 1 - pf[2] - pb[2], pf[2]],
                [0, pb[3], 1 - pf[3] - pb[3]]])
Nf = (sp.eye(3) - Pi).inv()
hv2 = Nf * sp.Matrix([0, 0, pf[3]])
check("splitting: resolvent vs fundamental matrix", (hvec - hv2).norm())
# E[steps 1_R](s) = sum_t t * (P^{t-1} v_R)_s = N^2 v_R (standard)
w2 = Nf * Nf * sp.Matrix([0, 0, pf[3]])
check("E[T 1_R]: w-equation vs N^2 v_R", (w - w2).norm())

# ---------------------------------------------------------------- G
hdr("G: phase-type tail vs brute path enumeration (T = first exit step)")
# P(T > t | start s) = (P_int^t 1)_s ; verify t = 1..6 by direct
# enumeration over paths of the 3-interior-site chain above
one = sp.Matrix([1, 1, 1])
import itertools
trans = {1: [(1, 1 - pf[1] - pb[1]), (2, pf[1]), ("L", pb[1])],
         2: [(2, 1 - pf[2] - pb[2]), (3, pf[2]), (1, pb[2])],
         3: [(3, 1 - pf[3] - pb[3]), ("R", pf[3]), (2, pb[3])]}
for s_ in (1, 2, 3):
    probs_alive = {s_: sp.Integer(1)}
    for t in range(1, 7):
        nxt = {}
        for st, pr in probs_alive.items():
            for st2, p2 in trans[st]:
                if st2 in ("L", "R"):
                    continue
                nxt[st2] = nxt.get(st2, sp.Integer(0)) + pr * p2
        probs_alive = nxt
        brute = sum(probs_alive.values())
        mat = ((Pi ** t) * one)[s_ - 1]
        check(f"P(T>{t} | start {s_})", sp.simplify(brute - mat))

print(f"\n{'ALL PASS' if ok_all else 'FAILURES PRESENT'}")
sys.exit(0 if ok_all else 1)
