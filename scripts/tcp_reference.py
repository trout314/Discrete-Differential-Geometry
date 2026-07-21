#!/usr/bin/env python3
"""Build T^3 triangulations of real tetrahedrally-close-packed (TCP) crystals.

A Frank-Kasper / TCP phase is BY DEFINITION a structure whose Delaunay
decomposition contains only tetrahedra, with every edge shared by 5 or 6 of
them. This script takes published Wyckoff positions of classical TCP phases,
builds an m x m x m supercell, computes the periodic Delaunay triangulation,
and quotients it to a simplicial triangulation of the 3-torus, saved as .mfd.

Structures (atoms/cell, mean coordination CN, mean edge degree q = 6-12/CN):
    a15   Cr3Si type,  8 atoms, CN 13.500, q 5.1111  (Z12_2 Z14_6)
    c15   MgCu2 Laves, 24 atoms, CN 13.333, q 5.1000  (Z12_16 Z16_8)
    c14   MgZn2 Laves, 12 atoms, CN 13.333, q 5.1000  (Z12_8 Z16_4)
    c36   MgNi2 Laves, 24 atoms, CN 13.333, q 5.1000  (Z12_16 Z16_8)
    sigma CrFe sigma,  30 atoms, CN 13.467, q 5.1090  (Z12_10 Z14_16 Z15_4)
          (dodecagonal-QC approximant; internal parameters from refinements)
    z     Zr4Al3,      7 atoms, CN 13.429, q 5.1064  (Z12_3 Z14_2 Z15_2)
    mu    W6Fe7 type, 39 atoms (hex cell), CN 13.385, q 5.1034
          (Z12_21 Z14_6 Z15_6 Z16_6; refined Co7Mo6 isotype parameters)
    r     Co-Cr-Mo R phase, 159 atoms (hex cell), CN 13.396, q 5.1042
          (Z12_81 Z14_36 Z15_18 Z16_24) — q matches the flat pin 5.1043
    p     Cr-Ni-Mo P phase, 56 atoms, CN 13.429, q 5.1064
          (Z12_24 Z14_20 Z15_8 Z16_4)
    delta MoNi delta, 56 atoms, CN 13.429, q 5.1064  (Z12_24 Z14_20 Z15_8 Z16_4)

Coordinates for c36/z/mu/r/p/delta are full expanded cells from published
refinements (COD/AFLOW; see per-block comments), independently validated by
periodic Voronoi CN census before inclusion. Two COD typos corrected: R-phase
A7 z = 0.3969 (COD 2310299 has 0.3696, transposed) and the P-phase cell
setting (COD 2310313 permutes the cell without permuting coordinates; we use
the original Pbnm cell a=9.070, b=16.983, c=4.752 of Shoemaker 1955/57).

Method notes:
  * Site positions get one deterministic ~1e-6 perturbation (fixed RNG seed)
    BEFORE tiling, so every periodic image is perturbed identically: this
    resolves co-spherical Delaunay degeneracies CONSISTENTLY across images
    (a per-point joggle like Qhull QJ could triangulate two images of the
    same degeneracy differently and break the quotient).
  * The supercell is tiled with a 2-unit-cell margin, Delaunay is computed on
    the extended cloud, and exactly one representative per translation class
    of tets is kept (centroid in the fundamental domain).
  * Validation: expected f-vector, no duplicate facets, Euler characteristic 0,
    orientable, the exact link sum rule (sum(6-deg)=12 at every vertex, i.e.
    a closed 3-manifold), edge degrees all in {5,6}, and the literature
    Z-class census. m >= 2 required (m >= 3 recommended) for a simplicial
    quotient.

Usage:
    python scripts/tcp_reference.py a15 -m 3                  # validate+save
    python scripts/tcp_reference.py c15 -m 6 --out data/tcp_reference/
    python scripts/tcp_reference.py all -m 3                  # all structures
"""
import argparse
import os
import sys

import numpy as np
from scipy.spatial import Delaunay

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_ROOT, "python"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from discrete_differential_geometry import Manifold
from fk_skeleton import edges_from_facets, vertex_class_census

SQ3 = np.sqrt(3.0)


def _fcc(base):
    """Apply +(0,0,0),(0,.5,.5),(.5,0,.5),(.5,.5,0) translations."""
    tr = np.array([[0, 0, 0], [0, .5, .5], [.5, 0, .5], [.5, .5, 0]])
    return np.vstack([(np.asarray(base) + t) % 1.0 for t in tr])


def _sigma_sites():
    """CrFe sigma phase, P4_2/mnm (136), 30 atoms; refined internal params."""
    def orbit_4f(x):
        return [(x, x, 0), (-x, -x, 0), (.5 + x, .5 - x, .5), (.5 - x, .5 + x, .5)]

    def orbit_8i(x, y):
        return [(x, y, 0), (-x, -y, 0), (y, x, 0), (-y, -x, 0),
                (.5 + x, .5 - y, .5), (.5 - x, .5 + y, .5),
                (.5 + y, .5 - x, .5), (.5 - y, .5 + x, .5)]

    def orbit_8j(x, z):
        return [(x, x, z), (x, x, -z), (-x, -x, z), (-x, -x, -z),
                (.5 + x, .5 - x, .5 + z), (.5 + x, .5 - x, .5 - z),
                (.5 - x, .5 + x, .5 + z), (.5 - x, .5 + x, .5 - z)]

    sites = [(0, 0, 0), (.5, .5, .5)]                     # 2a
    sites += orbit_4f(0.39864)
    sites += orbit_8i(0.46349, 0.13122)
    sites += orbit_8i(0.73933, 0.06609)
    sites += orbit_8j(0.18267, 0.25202)
    return np.array(sites) % 1.0


def _c14_sites(z_mg=0.06286, x_h=0.83048):
    """MgZn2, P6_3/mmc (194), 12 atoms (hexagonal axes)."""
    mg = [(1/3, 2/3, z_mg), (2/3, 1/3, z_mg + .5),
          (2/3, 1/3, -z_mg), (1/3, 2/3, .5 - z_mg)]
    zn2a = [(0, 0, 0), (0, 0, .5)]
    x = x_h
    zn6h = [(x, 2 * x, .25), (-2 * x, -x, .25), (x, -x, .25),
            (-x, -2 * x, .75), (2 * x, x, .75), (-x, x, .75)]
    return np.array(mg + zn2a + zn6h) % 1.0


def _hex(coa):
    return np.array([[1, 0, 0], [-.5, SQ3 / 2, 0], [0, 0, coa]])


def _parse(block):
    return np.array([[float(t) for t in ln.split()]
                     for ln in block.strip().splitlines()]) % 1.0


# C36 Laves MgNi2, P6_3/mmc, a=4.824 c=15.826; Komura & Tokunaga 1980, COD 2106100.
_C36 = """
0.00000 0.00000 0.09400
0.00000 0.00000 0.40600
0.00000 0.00000 0.59400
0.00000 0.00000 0.90600
0.00000 0.50000 0.00000
0.00000 0.50000 0.50000
0.16429 0.32858 0.25000
0.16429 0.83571 0.25000
0.32858 0.16429 0.75000
0.33333 0.66667 0.12514
0.33333 0.66667 0.37486
0.33333 0.66667 0.65580
0.33333 0.66667 0.84420
0.50000 0.00000 0.00000
0.50000 0.00000 0.50000
0.50000 0.50000 0.00000
0.50000 0.50000 0.50000
0.66667 0.33333 0.15580
0.66667 0.33333 0.34420
0.66667 0.33333 0.62514
0.66667 0.33333 0.87486
0.67142 0.83571 0.25000
0.83571 0.16429 0.75000
0.83571 0.67142 0.75000
"""

# Z phase Zr4Al3, P6/mmm, a=5.433 c=5.390; Wilson/Thomas/Spooner 1960 +
# Cenzual et al. 1991 symmetry revision (AFLOW A3B4_hP7_191_f_de, ICSD 150529).
_Z = """
0.00000 0.00000 0.25000
0.00000 0.00000 0.75000
0.00000 0.50000 0.00000
0.33333 0.66667 0.50000
0.50000 0.00000 0.00000
0.50000 0.50000 0.00000
0.66667 0.33333 0.50000
"""

# mu phase, R-3m hexagonal setting, 39 atoms; refined Co7Mo6 isotype
# parameters (Forsyth & d'Alte da Veiga 1962, COD 2310288), a=4.762 c=25.615.
_MU = """
0.00000 0.00000 0.00000
0.00000 0.00000 0.16550
0.00000 0.00000 0.34830
0.00000 0.00000 0.45180
0.00000 0.00000 0.54820
0.00000 0.00000 0.65170
0.00000 0.00000 0.83450
0.00007 0.50003 0.58953
0.16663 0.83337 0.92287
0.16663 0.33327 0.92287
0.16670 0.83330 0.74380
0.16670 0.33340 0.74380
0.33327 0.16663 0.07713
0.33333 0.66667 0.01497
0.33333 0.66667 0.11847
0.33333 0.66667 0.21487
0.33333 0.66667 0.31837
0.33333 0.66667 0.50117
0.33333 0.66667 0.66667
0.33333 0.66667 0.83217
0.33340 0.16670 0.25620
0.49997 0.50003 0.58953
0.49997 0.99993 0.58953
0.50003 0.00007 0.41047
0.50003 0.49997 0.41047
0.66660 0.83330 0.74380
0.66667 0.33333 0.16783
0.66667 0.33333 0.33333
0.66667 0.33333 0.49883
0.66667 0.33333 0.68163
0.66667 0.33333 0.78513
0.66667 0.33333 0.88153
0.66667 0.33333 0.98503
0.66673 0.83337 0.92287
0.83330 0.16670 0.25620
0.83330 0.66660 0.25620
0.83337 0.16663 0.07713
0.83337 0.66673 0.07713
0.99993 0.49997 0.41047
"""

# R phase Co-Cr-Mo, R-3 hexagonal setting, 159 atoms; Komura, Sly & Shoemaker
# 1960 (COD 2310299 with the A7 z typo corrected to 0.3969), a=10.903 c=19.342.
_R = """
0.00000 0.00000 0.07350
0.00000 0.00000 0.30440
0.00000 0.00000 0.50000
0.00000 0.00000 0.69560
0.00000 0.00000 0.92650
0.02120 0.13930 0.19620
0.02810 0.22500 0.73150
0.03300 0.25790 0.31830
0.04530 0.26710 0.87780
0.04940 0.17590 0.60310
0.05090 0.27900 0.10000
0.05433 0.43857 0.76667
0.06463 0.51117 0.13187
0.06623 0.44487 0.54447
0.07543 0.44177 0.98497
0.08633 0.57927 0.66867
0.08740 0.84040 0.00200
0.10523 0.71757 0.56667
0.10833 0.46977 0.39817
0.10843 0.69967 0.34837
0.11153 0.71197 0.78887
0.11320 0.26870 0.46520
0.11810 0.97880 0.19620
0.12650 0.95060 0.60310
0.13643 0.69477 0.93517
0.13930 0.11810 0.80380
0.15550 0.88680 0.46520
0.15743 0.54017 0.26977
0.15960 0.24700 0.00200
0.17373 0.41967 0.66467
0.17590 0.12650 0.39690
0.17783 0.77987 0.20147
0.19403 0.54857 0.86287
0.19690 0.97190 0.73150
0.20683 0.71607 0.06357
0.21523 0.68787 0.47047
0.22013 0.39797 0.20147
0.22180 0.95470 0.87780
0.22490 0.96700 0.31830
0.22500 0.19690 0.26850
0.22810 0.94910 0.10000
0.24593 0.82627 0.66467
0.24700 0.08740 0.99800
0.25790 0.22490 0.68170
0.26710 0.22180 0.12220
0.26870 0.15550 0.53480
0.27900 0.22810 0.90000
0.28243 0.38767 0.56667
0.28393 0.49077 0.06357
0.28803 0.39957 0.78887
0.30033 0.40877 0.34837
0.30523 0.44167 0.93517
0.31213 0.52737 0.47047
0.33333 0.66667 0.16667
0.33333 0.66667 0.36227
0.33333 0.66667 0.59317
0.33333 0.66667 0.74017
0.33333 0.66667 0.97107
0.35453 0.80597 0.86287
0.36143 0.89167 0.39817
0.36633 0.92457 0.98497
0.37863 0.93377 0.54447
0.38273 0.84257 0.26977
0.38423 0.94567 0.76667
0.38767 0.10523 0.43333
0.39797 0.17783 0.79853
0.39957 0.11153 0.21113
0.40877 0.10843 0.65163
0.41967 0.24593 0.33533
0.42073 0.50707 0.66867
0.43857 0.38423 0.23333
0.44167 0.13643 0.06483
0.44177 0.36633 0.01503
0.44487 0.37863 0.45553
0.44653 0.93537 0.13187
0.45143 0.64547 0.86287
0.45983 0.61727 0.26977
0.46977 0.36143 0.60183
0.47263 0.78477 0.47047
0.48883 0.55347 0.13187
0.49077 0.20683 0.93643
0.49293 0.91367 0.66867
0.50707 0.08633 0.33133
0.50923 0.79317 0.06357
0.51117 0.44653 0.86813
0.52737 0.21523 0.52953
0.53023 0.63857 0.39817
0.54017 0.38273 0.73023
0.54857 0.35453 0.13713
0.55347 0.06463 0.86813
0.55513 0.62137 0.54447
0.55823 0.63367 0.98497
0.55833 0.86357 0.93517
0.56143 0.61577 0.76667
0.57927 0.49293 0.33133
0.58033 0.75407 0.66467
0.59123 0.89157 0.34837
0.60043 0.88847 0.78887
0.60203 0.82217 0.20147
0.61233 0.89477 0.56667
0.61577 0.05433 0.23333
0.61727 0.15743 0.73023
0.62137 0.06623 0.45553
0.63367 0.07543 0.01503
0.63857 0.10833 0.60183
0.64547 0.19403 0.13713
0.66667 0.33333 0.02893
0.66667 0.33333 0.25983
0.66667 0.33333 0.40683
0.66667 0.33333 0.63773
0.66667 0.33333 0.83333
0.68787 0.47263 0.52953
0.69477 0.55833 0.06483
0.69967 0.59123 0.65163
0.71197 0.60043 0.21113
0.71607 0.50923 0.93643
0.71757 0.61233 0.43333
0.72100 0.77190 0.10000
0.73130 0.84450 0.46520
0.73290 0.77820 0.87780
0.74210 0.77510 0.31830
0.75300 0.91260 0.00200
0.75407 0.17373 0.33533
0.77190 0.05090 0.90000
0.77500 0.80310 0.73150
0.77510 0.03300 0.68170
0.77820 0.04530 0.12220
0.77987 0.60203 0.79853
0.78477 0.31213 0.52953
0.79317 0.28393 0.93643
0.80310 0.02810 0.26850
0.80597 0.45143 0.13713
0.82217 0.22013 0.79853
0.82410 0.87350 0.60310
0.82627 0.58033 0.33533
0.84040 0.75300 0.99800
0.84257 0.45983 0.73023
0.84450 0.11320 0.53480
0.86070 0.88190 0.19620
0.86357 0.30523 0.06483
0.87350 0.04940 0.39690
0.88190 0.02120 0.80380
0.88680 0.73130 0.53480
0.88847 0.28803 0.21113
0.89157 0.30033 0.65163
0.89167 0.53023 0.60183
0.89477 0.28243 0.43333
0.91260 0.15960 0.99800
0.91367 0.42073 0.33133
0.92457 0.55823 0.01503
0.93377 0.55513 0.45553
0.93537 0.48883 0.86813
0.94567 0.56143 0.23333
0.94910 0.72100 0.90000
0.95060 0.82410 0.39690
0.95470 0.73290 0.12220
0.96700 0.74210 0.68170
0.97190 0.77500 0.26850
0.97880 0.86070 0.80380
"""

# P phase Cr-Ni-Mo, Pbnm, a=9.070 b=16.983 c=4.752; Shoemaker, Shoemaker &
# Wilson 1957 (COD 2310313 coordinates in the ORIGINAL Pbnm cell — the COD
# entry's permuted cell is wrong).
_P = """
0.02020 0.46450 0.75000
0.02540 0.95360 0.25000
0.06170 0.63500 0.75000
0.07370 0.11340 0.25000
0.10580 0.31810 0.75000
0.11320 0.78830 0.00080
0.11320 0.78830 0.49920
0.13630 0.25470 0.25000
0.16500 0.14470 0.75000
0.17430 0.65780 0.25000
0.18480 0.92200 0.75000
0.19880 0.40470 0.25000
0.24960 0.03750 0.50140
0.24960 0.03750 0.99860
0.25040 0.53750 0.50140
0.25040 0.53750 0.99860
0.30120 0.90470 0.25000
0.31520 0.42200 0.75000
0.32570 0.15780 0.25000
0.33500 0.64470 0.75000
0.36370 0.75470 0.25000
0.38680 0.28830 0.00080
0.38680 0.28830 0.49920
0.39420 0.81810 0.75000
0.42630 0.61340 0.25000
0.43830 0.13500 0.75000
0.47460 0.45360 0.25000
0.47980 0.96450 0.75000
0.52020 0.03550 0.25000
0.52540 0.54640 0.75000
0.56170 0.86500 0.25000
0.57370 0.38660 0.75000
0.60580 0.18190 0.25000
0.61320 0.71170 0.50080
0.61320 0.71170 0.99920
0.63630 0.24530 0.75000
0.66500 0.35530 0.25000
0.67430 0.84220 0.75000
0.68480 0.57800 0.25000
0.69880 0.09530 0.75000
0.74960 0.46250 0.00140
0.74960 0.46250 0.49860
0.75040 0.96250 0.00140
0.75040 0.96250 0.49860
0.80120 0.59530 0.75000
0.81520 0.07800 0.25000
0.82570 0.34220 0.75000
0.83500 0.85530 0.25000
0.86370 0.74530 0.75000
0.88680 0.21170 0.50080
0.88680 0.21170 0.99920
0.89420 0.68190 0.25000
0.92630 0.88660 0.75000
0.93830 0.36500 0.25000
0.97460 0.04640 0.75000
0.97980 0.53550 0.25000
"""

# delta phase MoNi, P2_1 2_1 2_1, a=9.108 b=9.108 c=8.852; Shoemaker &
# Shoemaker 1963, COD 2310246.
_DELTA = """
0.00290 0.19690 0.67670
0.03380 0.33980 0.18070
0.04810 0.88470 0.03220
0.05760 0.63380 0.09710
0.06800 0.14420 0.95290
0.10310 0.41920 0.91330
0.11180 0.94770 0.77480
0.12320 0.56420 0.35670
0.13370 0.07070 0.21570
0.17630 0.48320 0.64250
0.18140 0.75360 0.57400
0.18850 0.01570 0.49600
0.22890 0.28650 0.40980
0.23520 0.80070 0.24860
0.26480 0.19930 0.74860
0.27110 0.71350 0.90980
0.31150 0.98430 0.99600
0.31860 0.24640 0.07400
0.32370 0.51680 0.14250
0.36630 0.92930 0.71570
0.37680 0.43580 0.85670
0.38820 0.05230 0.27480
0.39690 0.58080 0.41330
0.43200 0.85580 0.45290
0.44240 0.36620 0.59710
0.45190 0.11530 0.53220
0.46620 0.66020 0.68070
0.49710 0.80310 0.17670
0.50290 0.30310 0.32330
0.53380 0.16020 0.81930
0.54810 0.61530 0.96780
0.55760 0.86620 0.90290
0.56800 0.35580 0.04710
0.60310 0.08080 0.08670
0.61180 0.55230 0.22520
0.62320 0.93580 0.64330
0.63370 0.42930 0.78430
0.67630 0.01680 0.35750
0.68140 0.74640 0.42600
0.68850 0.48430 0.50400
0.72890 0.21350 0.59020
0.73520 0.69930 0.75140
0.76480 0.30070 0.25140
0.77110 0.78650 0.09020
0.81150 0.51570 0.00400
0.81860 0.25360 0.92600
0.82370 0.98320 0.85750
0.86630 0.57070 0.28430
0.87680 0.06420 0.14330
0.88820 0.44770 0.72520
0.89690 0.91920 0.58670
0.93200 0.64420 0.54710
0.94240 0.13380 0.40290
0.95190 0.38470 0.46780
0.96620 0.83980 0.31930
0.99710 0.69690 0.82330
"""

STRUCTURES = {
    # name: (lattice rows (unit a=1), fractional sites, CN_mean, census {Z: count})
    "a15": (np.eye(3),
            np.array([(0, 0, 0), (.5, .5, .5),
                      (.25, 0, .5), (.75, 0, .5),
                      (.5, .25, 0), (.5, .75, 0),
                      (0, .5, .25), (0, .5, .75)]),
            13.5, {"Z12": 2, "Z14": 6, "Z15": 0, "Z16": 0}),
    "c15": (np.eye(3),
            np.vstack([_fcc([(0, 0, 0), (.25, .25, .25)]),          # 8a (Mg)
                       _fcc([(5/8, 5/8, 5/8), (5/8, 7/8, 7/8),      # 16d (Cu)
                             (7/8, 5/8, 7/8), (7/8, 7/8, 5/8)])]),
            40 / 3, {"Z12": 16, "Z14": 0, "Z15": 0, "Z16": 8}),
    "c14": (np.array([[1, 0, 0], [-.5, SQ3 / 2, 0], [0, 0, np.sqrt(8 / 3)]]),
            _c14_sites(),
            40 / 3, {"Z12": 8, "Z14": 0, "Z15": 0, "Z16": 4}),
    "sigma": (np.diag([1, 1, 0.51759]),
              _sigma_sites(),
              404 / 30, {"Z12": 10, "Z14": 16, "Z15": 4, "Z16": 0}),
    "c36": (_hex(15.826 / 4.824), _parse(_C36),
            40 / 3, {"Z12": 16, "Z14": 0, "Z15": 0, "Z16": 8}),
    "z": (_hex(5.390 / 5.433), _parse(_Z),
          94 / 7, {"Z12": 3, "Z14": 2, "Z15": 2, "Z16": 0}),
    "mu": (_hex(25.615 / 4.762), _parse(_MU),
           174 / 13, {"Z12": 21, "Z14": 6, "Z15": 6, "Z16": 6}),
    "r": (_hex(19.342 / 10.903), _parse(_R),
          710 / 53, {"Z12": 81, "Z14": 36, "Z15": 18, "Z16": 24}),
    "p": (np.diag([1, 16.983 / 9.070, 4.752 / 9.070]), _parse(_P),
          94 / 7, {"Z12": 24, "Z14": 20, "Z15": 8, "Z16": 4}),
    "delta": (np.diag([1, 1, 8.852 / 9.108]), _parse(_DELTA),
              94 / 7, {"Z12": 24, "Z14": 20, "Z15": 8, "Z16": 4}),
}


def build_t3_triangulation(name, m, perturb=1e-6, pad=2):
    """Return (facets ndarray (f3,4), n_vertices) for an m^3 supercell of the
    named TCP structure, as a triangulation of T^3."""
    if m < 2:
        raise SystemExit("need m >= 2 for a simplicial quotient (m >= 3 safer)")
    L, sites, cn, _ = STRUCTURES[name]
    ns = len(sites)
    rng = np.random.default_rng(12345)
    sites = sites + perturb * rng.standard_normal(sites.shape)  # once, pre-tiling

    # Tile cells in [-pad, m+pad)^3; canonical id = site + m^3-cell (mod m).
    rng_cells = np.arange(-pad, m + pad)
    cells = np.array(np.meshgrid(rng_cells, rng_cells, rng_cells,
                                 indexing="ij")).reshape(3, -1).T   # (nc, 3)
    frac = (sites[None, :, :] + cells[:, None, :]).reshape(-1, 3)   # (nc*ns, 3)
    pts = frac @ L
    cell_mod = (cells[:, None, :] % m).repeat(ns, axis=1).reshape(-1, 3)
    site_idx = np.tile(np.arange(ns), len(cells))
    canon = ((cell_mod[:, 0] * m + cell_mod[:, 1]) * m + cell_mod[:, 2]) * ns + site_idx

    tri = Delaunay(pts)
    tets = tri.simplices                                            # (nt, 4)
    # keep one representative per translation class: centroid in [0, m)^3
    cent = pts[tets].mean(axis=1) @ np.linalg.inv(L)                # unit-cell frac
    keep = np.all((cent >= 0) & (cent < m), axis=1)
    fac = canon[tets[keep]]
    fac.sort(axis=1)
    if np.any(fac[:, :-1] == fac[:, 1:]):
        raise SystemExit(f"{name} m={m}: quotient tet with repeated vertex — increase m")
    fac = fac[np.lexsort(fac.T[::-1])]
    uniq = np.unique(fac, axis=0)
    if len(uniq) != len(fac):
        raise SystemExit(f"{name} m={m}: duplicate facets after quotient "
                         f"({len(fac)} vs {len(uniq)} unique) — degeneracy inconsistency")
    f3_expected = round(ns * m ** 3 * (cn / 2 - 1))
    if len(fac) != f3_expected:
        raise SystemExit(f"{name} m={m}: f3 = {len(fac)} != expected {f3_expected} "
                         f"(bad positions or selection)")
    return fac, ns * m ** 3


def validate_and_save(name, m, out_dir):
    fac, f0_expected = build_t3_triangulation(name, m)
    _, _, cn, census_per_cell = STRUCTURES[name]

    eu, edeg, V = edges_from_facets(fac)
    assert V == f0_expected, f"vertex count {V} != {f0_expected}"
    hist = np.bincount(edeg)
    ok_56 = hist.sum() == (hist[5] if len(hist) > 5 else 0) + (hist[6] if len(hist) > 6 else 0)
    fz, n_broken = vertex_class_census(eu, edeg, V)   # asserts link sum rule
    m3 = m ** 3

    mfd = Manifold(3, fac.tolist())
    euler = mfd.euler_characteristic
    orient = mfd.is_orientable

    print(f"[{name} m={m}] f0={V} f1={len(eu)} f3={len(fac)} "
          f"q̄={edeg.mean():.4f} (expect {6 - 12 / cn:.4f})")
    print(f"  edge degrees: " + ", ".join(f"deg{d}:{c}" for d, c in enumerate(hist) if c)
          + ("  [all {5,6} OK]" if ok_56 else "  [NON-{5,6} EDGES — positions off]"))
    print(f"  Z census (per cell, expect {census_per_cell}): "
          + ", ".join(f"{k}:{fz[k] * V / m3:.2f}" for k in ("Z12", "Z14", "Z15", "Z16"))
          + f"  fFK={fz['Z12'] + fz['Z14'] + fz['Z15'] + fz['Z16']:.4f}")
    print(f"  euler={euler} orientable={orient} broken56={n_broken}")
    good = (ok_56 and euler == 0 and orient
            and all(abs(fz[k] * V / m3 - census_per_cell[k]) < 1e-6
                    for k in census_per_cell))
    os.makedirs(out_dir, exist_ok=True)
    path = os.path.join(out_dir, f"T3_{name.upper()}_m{m}_N{len(fac)}.mfd")
    from seed_utils import get_git_info, make_leg, history_fields
    commit, dirty = get_git_info()
    qbar = float(edeg.mean())
    # Root leg: this crystal is a from-scratch lineage ORIGIN (built exactly from
    # Wyckoff positions, not sampled). Stamping it means downstream melt/dope/
    # quench runs append to it via read_history() and can never lose provenance.
    # See the history section in tools/seed_utils.py (is_root_leg_from).
    root_leg = make_leg(
        "build",
        {"struct": name, "m": m, "cn": cn, "qbar": round(qbar, 6),
         "source": "wyckoff-delaunay"},
        sweeps=0, from_=f"crystal:{name}@wyckoff",
        commit=commit, dirty=dirty, tried=0, accepted=0)
    mfd.save(path, comments=[
        f"topology = T3", f"structure = {name}", f"supercell = {m}",
        f"source = scripts/tcp_reference.py (periodic Delaunay of Wyckoff positions)",
        f"git_commit = {commit}", f"git_dirty = {'true' if dirty else 'false'}",
        f"mean_edge_degree = {qbar:.6f}",
        f"validation = {'PASS' if good else 'PARTIAL (see census)'}",
    ] + history_fields([root_leg]))
    print(f"  {'VALIDATED' if good else 'PARTIAL'} -> {path}")
    return good


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("structure", choices=list(STRUCTURES) + ["all"])
    ap.add_argument("-m", type=int, default=3, help="supercell multiplier (>=2)")
    ap.add_argument("--out", default=os.path.join(_ROOT, "data", "tcp_reference"))
    args = ap.parse_args()
    sys.path.insert(0, os.path.join(_ROOT, "tools"))
    names = list(STRUCTURES) if args.structure == "all" else [args.structure]
    results = {n: validate_and_save(n, args.m, args.out) for n in names}
    if not all(results.values()):
        sys.exit(1)


if __name__ == "__main__":
    main()
