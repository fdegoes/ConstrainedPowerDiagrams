#!/usr/bin/env python3
"""Is the failure of H from the targets caused by close generator pairs?

Realizable equal-volume cases on the unit square whose true generators Y are
uniform random points with a minimum spacing r_min (dart throwing, i.e.
Poisson-disk sampling), for r_min from 0 (uniform) up to a fraction of the
cell size h = sqrt(1/n). Everything else about Y stays random, so a change in
outcome isolates the close pairs.

Reported per case: the smallest generator distance, the nearest-neighbor
coefficient of variation, the start distance |b - Y| (homothety removed), and
whether full H from the targets converges.

Usage:
  external/bourne/.venv/bin/python benchmarks/bourne/spacing_sweep.py [--n 243 500 728] [--seeds 1 2 3]
      [--rmin 0 0.1 0.2 0.3 0.4 0.5]
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
from irregularity_sweep import irregularity, start_distance  # noqa: E402

ref = rb.ref


def dart_throwing(n, rmin, rng, max_tries=200000):
    """n uniform points in the unit square at pairwise distance >= rmin
    (sequential rejection sampling). Returns None if it cannot fit them."""
    pts = []
    tree_pts = np.empty((0, 2))
    tries = 0
    while len(pts) < n and tries < max_tries:
        tries += 1
        p = rng.random(2)
        if rmin > 0 and len(pts):
            if np.min(np.sum((tree_pts - p) ** 2, 1)) < rmin * rmin:
                continue
        pts.append(p)
        tree_pts = np.array(pts)
    return np.array(pts) if len(pts) == n else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", nargs="+", type=int, default=[243, 500, 728])
    ap.add_argument("--seeds", nargs="+", type=int, default=[1, 2, 3])
    ap.add_argument("--rmin", nargs="+", type=float, default=[0.0, 0.1, 0.2, 0.3, 0.4, 0.5],
                    help="minimum spacing in cell sizes h = sqrt(1/n)")
    args = ap.parse_args()
    base = [c for c in rb.make_cases() if c.name == "2026-voronoi"][0]
    print("%-5s %-6s %-4s %-12s %-8s %-10s %s" % ("n", "rmin/h", "seed", "min gap/h", "NN cv", "|b-Y|/h", "H from targets"))
    for n in args.n:
        v = np.ones(n) / n
        h = math.sqrt(1.0 / n)
        for rmin in args.rmin:
            for seed in args.seeds:
                Y = dart_throwing(n, rmin * h, np.random.default_rng(seed))
                if Y is None:
                    print("%-5d %-6g %-4d  could not place the points" % (n, rmin, seed), flush=True)
                    continue
                _, b, _ = ref.solve_OT_problem(Y, v, base.om, relTol=1e-8)
                c = copy.copy(base)
                c.v, c.b, c.sites, c.truth = v, b, b.copy(), Y
                r = df.solve(c, "H", math.inf, maxiter=1000, tol=1e-10)
                ok = max(r["volumeResidual"], r["shiftResidual"]) <= 1e-8
                gap = cKDTree(Y).query(Y, k=2)[0][:, 1].min() / h
                print("%-5d %-6g %-4d %-12.3f %-8.3f %-10.2f %s" % (
                    n, rmin, seed, gap, irregularity(Y), start_distance(b, Y, h),
                    ("converged %d it" % r["iterations"]) if ok else ("FAIL (%s)" % r["status"])), flush=True)


if __name__ == "__main__":
    main()
