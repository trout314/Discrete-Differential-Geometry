#!/usr/bin/env python3
"""Generate the JSON parameter file for the S3 seed triangulation library.

Outputs a JSON array of parameter dictionaries to stdout.
Each entry specifies one seed job (one replica at one parameter combo).

Usage:
    python3 tools/generate_seed_params.py > tools/seed_params.json
"""

import json
import math
import sys


def main():
    # --- Parameter grid ---
    # 21 sizes: 10^2, 10^2.25, 10^2.5, ..., 10^7 (geometric, 10^(1/4) spacing)
    sizes = [round(10 ** (2 + 0.25 * i)) for i in range(21)]

    hinge_degree_targets = [5.0043, 5.1043, 5.2043]
    codim3_degree_variance_coefs = [0, 0.1, 0.5]

    # Fixed parameters
    num_facets_coef = 0.1
    num_hinges_coef = 0.1
    hinge_degree_variance_coef = 0.0

    # --- Replicas per size ---
    def replicas_for_size(n):
        if n <= 562:
            return 128
        elif n <= 5623:
            return 64
        elif n <= 56234:
            return 32
        elif n <= 562341:
            return 16
        elif n <= 5623413:
            return 8
        else:
            return 4

    # --- Build parameter list ---
    params = []
    for n in sizes:
        growth_step_size = max(50, n // 20)
        num_replicas = replicas_for_size(n)
        for hdt in hinge_degree_targets:
            for cdvc in codim3_degree_variance_coefs:
                for replica in range(num_replicas):
                    params.append({
                        "topology": "S3",
                        "num_facets_target": n,
                        "num_facets_coef": num_facets_coef,
                        "hinge_degree_target": hdt,
                        "num_hinges_coef": num_hinges_coef,
                        "hinge_degree_variance_coef": hinge_degree_variance_coef,
                        "codim3_degree_variance_coef": cdvc,
                        "growth_step_size": growth_step_size,
                        "seed_index": replica,
                    })

    # --- Summary to stderr ---
    print(f"Generated {len(params)} seed jobs", file=sys.stderr)
    size_summary = []
    for n in sizes:
        r = replicas_for_size(n)
        count = r * len(hinge_degree_targets) * len(codim3_degree_variance_coefs)
        size_summary.append(f"  N={n:>10,}: {r:>3} replicas x {len(hinge_degree_targets) * len(codim3_degree_variance_coefs)} combos = {count}")
    print("\n".join(size_summary), file=sys.stderr)

    # --- Output JSON to stdout ---
    json.dump(params, sys.stdout, indent=2)
    print()  # trailing newline


if __name__ == "__main__":
    main()
