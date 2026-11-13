#!/usr/bin/env python3
"""How irregular can the true generators be before H, started from the
targets, fails? Realizable cases with equal volumes on the unit square:

  Y   true generators, from one of these families, in increasing regularity
        uniform   uniform random points
        lloydK    uniform points after K Lloyd iterations at equal volumes
                  (K = 1, 3, 10, 30: towards a centroidal diagram)
        jitterS   a jittered square grid, each point moved by up to S cell
                  sizes (S = 0.5, 0.3, 0.1)
  b   centroids of the equal-volume power diagram of Y (exact targets)

For each case the solver (H from sites at the targets, our start) runs to
1e-10, and the table reports the start distance |b - Y| (homothety removed,
in cell sizes h = sqrt(1/n)), the irregularity of Y (the coefficient of
variation of the nearest-neighbor distance), and convergence.

Usage:
  external/bourne/.venv/bin/python benchmarks/bourne/irregularity_sweep.py [--n 243 500 728] [--seeds 1 2 3]
"""

import argparse
import copy
import math
import os
import sys

import numpy as np
from scipy.spatial import cKDTree

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import diagnose_failures as df  # noqa: E402
import bourne_cases as rb  # noqa: E402

ref = rb.ref


def lloyd(Y, v, om, iterations):
    """Lloyd iterations at the prescribed volumes: sites <- centroids."""
    for _ in range(iterations):
        _, c, _ = ref.solve_OT_problem(Y, v, om, relTol=1e-8)
        Y = c
    return Y


def jittered_grid(n, s, rng):
    """n points on a jittered grid (k x k, k = ceil(sqrt(n)), first n cells),
    each moved uniformly by up to s cell sizes."""
    k = int(math.ceil(math.sqrt(n)))
    cell = 1.0 / k
    idx = rng.permutation(k * k)[:n]
    x = (idx % k + 0.5) * cell
    y = (idx // k + 0.5) * cell
    P = np.stack([x, y], 1) + s * cell * rng.uniform(-1, 1, (n, 2))
    return np.clip(P, 1e-3, 1 - 1e-3)


def start_distance(b, Y, h):
    n = len(Y)
    X = np.vstack([Y.ravel(), np.tile([1, 0], n), np.tile([0, 1], n)]).T
    coef = np.linalg.lstsq(X, b.ravel(), rcond=None)[0]
    return np.abs(X @ coef - b.ravel()).max() / abs(coef[0]) / h


def irregularity(Y):
    d, _ = cKDTree(Y).query(Y, k=2)
    nn = d[:, 1]
    return nn.std() / nn.mean()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", nargs="+", type=int, default=[243, 500, 728])
    ap.add_argument("--seeds", nargs="+", type=int, default=[1, 2, 3])
    args = ap.parse_args()
    base = [c for c in rb.make_cases() if c.name == "2026-voronoi"][0]
    families = [("uniform", None), ("lloyd1", 1), ("lloyd3", 3), ("lloyd10", 10), ("lloyd30", 30),
                ("jitter0.5", 0.5), ("jitter0.3", 0.3), ("jitter0.1", 0.1)]
    print("%-5s %-10s %-4s %-12s %-12s %s" % ("n", "family", "seed", "|b-Y|/h", "NN cv", "H from targets"))
    for n in args.n:
        v = np.ones(n) / n
        h = math.sqrt(1.0 / n)
        for fam, p in families:
            for seed in args.seeds:
                rng = np.random.default_rng(seed)
                if fam.startswith("jitter"):
                    Y = jittered_grid(n, p, rng)
                else:
                    Y = rng.random((n, 2))
                    if p:
                        Y = lloyd(Y, v, base.om, p)
                try:
                    _, b, _ = ref.solve_OT_problem(Y, v, base.om, relTol=1e-8)
                except Exception as e:  # noqa: BLE001
                    print("%-5d %-10s %-4d  OT failed: %s" % (n, fam, seed, e), flush=True)
                    continue
                c = copy.copy(base)
                c.v, c.b, c.sites, c.truth = v, b, b.copy(), Y
                r = df.solve(c, "H", math.inf, maxiter=300, tol=1e-10)
                ok = max(r["volumeResidual"], r["shiftResidual"]) <= 1e-8
                print("%-5d %-10s %-4d %-12.2f %-12.3f %s" % (
                    n, fam, seed, start_distance(b, Y, h), irregularity(Y),
                    ("converged %d it" % r["iterations"]) if ok else ("FAIL (%s)" % r["status"])), flush=True)


if __name__ == "__main__":
    main()
