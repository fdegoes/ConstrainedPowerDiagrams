#!/usr/bin/env python3
"""Why do my solvers fail on some Bourne cases? For each case and solver:

  1. the trajectory: residuals, smallest cell, and smallest generator distance
     (relative to the mean cell size h = sqrt(mean Vbar)) per iterate;
  2. an alpha sweep of H_alpha, as multiples of the distinctness bound
     eq:regularized-distinct, and at 0.99 times the per-cell estimate
     eq:regularized-practical, with the smallest generator distance at the end;
  3. for the aligned seeds, whether the centroids are invariant under
     per-axis warps of the generators (flat directions of H beyond the
     homotheties, so that the maximizer is not unique).

Usage:
  external/bourne/.venv/bin/python benchmarks/bourne/diagnose_failures.py [--case NAME ...]
"""

import argparse
import json
import math
import os
import subprocess
import sys
import tempfile

import numpy as np
from scipy.spatial.distance import pdist

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import bourne_cases as rb  # noqa: E402

ref = rb.ref


def solve(case, solver, alpha, start=0.0, maxiter=200, tol=1e-13, sites=None, eta=-1.0, hessmin=-1.0, pairreg=None,
          restarts=None, trajectory=True):
    """One benchBourne solve (H or H_alpha from MakeJointStart at s = start)."""
    sites = case.sites if sites is None else sites
    lines = ["box %r %r %r %r" % tuple(float(x) for x in case.box), "n %d" % case.n, "solver %s" % solver,
             "alpha %r" % (float(alpha) if math.isfinite(alpha) else 1e300), "start %r" % float(start),
             "volumetol %r" % tol, "shifttol %r" % tol, "maxiter %d" % maxiter, "eta %r" % float(eta),
            ] + (["noclamp"] if hessmin == 0 else [])
    if not trajectory:
        lines.append("notrajectory")
    if restarts is not None:  # (count, growth of the core radius)
        lines.append("restarts %d %r" % (restarts[0], float(restarts[1])))
    if pairreg is not None:  # R, scale
        lines.append("pairreg " + " ".join(repr(float(x)) for x in pairreg))
    lines += ["volumes " + " ".join(repr(float(x)) for x in case.v),
             "targets " + " ".join(repr(float(x)) for x in case.b.ravel()),
             "sites " + " ".join(repr(float(x)) for x in np.asarray(sites).ravel())]
    with tempfile.TemporaryDirectory() as d:
        cp, op = os.path.join(d, "c.txt"), os.path.join(d, "o.json")
        open(cp, "w").write("\n".join(lines) + "\n")
        subprocess.run([rb.BENCH, cp, op], check=True)
        return json.load(open(op))


def h_of(case):
    return math.sqrt(case.v.mean())


def gen_dist(p, h):
    return pdist(np.reshape(p, (-1, 2))).min() / h


def trajectory(case, solver, alpha, start=0.0):
    r = solve(case, solver, alpha, start)
    h = h_of(case)
    print("  %s alpha=%s start=%g: %s after %d it (%d builds)" % (
        solver, "inf" if not math.isfinite(alpha) else "%.3g" % alpha, start, r["status"], r["iterations"],
        r["numBuilds"]))
    T = r["trajectory"]
    step = max(1, len(T) // 12)
    for k in list(range(0, len(T), step)) + [len(T) - 1]:
        it = T[k]
        V = np.asarray(it["volumes"])
        print("    it %3d  vol %.1e  shift %.1e  min V/Vbar %.1e  min |p_i-p_j|/h %.1e" % (
            k, it["volumeResidual"], it["shiftResidual"], (V / case.v).min(), gen_dist(it["generators"], h)))
    return r


def alpha_sweep(case, factors):
    bound, practical = rb.distinct_alpha(case), rb.practical_alpha(case)
    h = h_of(case)
    print("  H_alpha sweep (bound %.3g, practical %.3g, mean cell size h %.3g):" % (bound, practical, h))
    for f in sorted(list(factors) + [0.99 * practical / bound]):
        a = f * bound
        r = solve(case, "Halpha", a)
        T = r["trajectory"]
        best = min(T, key=lambda it: max(it["volumeResidual"], it["shiftResidual"]))
        # Generators of H_alpha scale like alpha: compare their distance with
        # the separation of the sites (the bound is a sufficient condition).
        print("    alpha %9.3g (%7.3g x bound): %-20s %3d it  best res %.1e  final min|p_i-p_j|/h %.1e" % (
            a, f, r["status"], r["iterations"], max(best["volumeResidual"], best["shiftResidual"]),
            gen_dist(r["generators"], h) / max(a, 1.0)))


def fw_reference(case):
    """Their FW answer as an H_alpha point (alpha = t - 1): generators
    alpha (s b - B), and whether B is the centroid set of the hard power
    diagram of those generators (else B is a fractional partition, i.e. the
    maximizer of H_alpha has tied generators)."""
    alpha = case.t - 1
    y = case.t / (case.t - 1) * case.b
    B, k, st = ref.convex_order_projection_frank_wolfe(case.om, y, case.start, R0=case.R0, masses=case.v,
                                                       niter=case.niter, tol=case.gap)
    P = alpha * (y - B)
    h = h_of(case)
    _, c, _ = ref.solve_OT_problem(y - B, case.v, case.om, relTol=1e-6)
    L = math.sqrt(case.om.measure())
    d = np.sort(pdist(P) / h)
    print("  FW (t=%g, %s after %d it): generators alpha(sb - B), smallest distances / h: %s" % (
        case.t, st, k, " ".join("%.1e" % x for x in d[:5])))
    print("    max |B - c(hard diagram)| / L = %.1e  (round-off if B is a power diagram)" % (np.abs(B - c).max() / L))


def warp_test(case):
    """Centroids under a monotone per-axis warp x -> x + e sin(2 pi x) of the
    true generators (prescribed volumes kept by the OT solve)."""
    Y = case.truth
    _, c0, _ = ref.solve_OT_problem(Y, case.v, case.om, relTol=1e-8)
    L = math.sqrt(case.om.measure())
    for e in (1e-3, 1e-2, 3e-2):
        W = Y.copy()
        W[:, 0] = Y[:, 0] + e * np.sin(2 * np.pi * Y[:, 0])
        _, c1, _ = ref.solve_OT_problem(W, case.v, case.om, relTol=1e-8)
        print("    x-warp %.0e: generators move %.1e h, centroids move %.1e L" % (
            e, np.abs(W - Y).max() / h_of(case), np.abs(c1 - c0).max() / L))


def main():
    cases = {c.name: c for c in rb.make_cases()}
    default = ["2025-invert", "2025-perturbed-large", "ebsd", "2026-aligned"]
    ap = argparse.ArgumentParser()
    ap.add_argument("--case", nargs="+", choices=list(cases), default=default)
    args = ap.parse_args()
    factors = [0.99, 3, 10, 30, 100, 300, 1000, 3000]
    for name in args.case:
        c = cases[name]
        print("== %s (n = %d, %s)" % (name, c.n, c.group))
        if name == "2025-invert":
            for s in (0.0, 0.5, 1.0):
                trajectory(c, "H", math.inf, s)
            continue
        trajectory(c, "H", math.inf)
        if c.paper == "2026":
            trajectory(c, "Halpha", c.t - 1)
            fw_reference(c)
        alpha_sweep(c, factors)
        if name == "2026-aligned":
            print("  flat directions of the aligned seeds:")
            warp_test(c)


if __name__ == "__main__":
    main()
