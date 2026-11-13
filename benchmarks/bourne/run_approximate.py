#!/usr/bin/env python3
"""Benchmark on the approximate cases of Bourne et al. (no power diagram has
the prescribed volumes and centroids; only a fit is possible):
2025-perturbed-small, 2025-perturbed-large and ebsd (the EBSD map of both
papers: our H_alpha once, in the unit square; their 2025 least squares and
2026 Frank-Wolfe with each paper's scaling, as theirs-2025 and
theirs-2026). Runs
their method (SLSQP, 2025; Frank-Wolfe, 2026) and our max H_alpha over a sweep
of stiffnesses, under their termination rule and the common residual of the
problem each method solves (bourne_cases.py). Their 2025 methods have no
common rule here (H is infeasible, and their least squares has none).

Usage (all cases, methods, rules and tolerances by default):
  external/bourne/.venv/bin/python benchmarks/bourne/run_approximate.py
      [--case NAME ...] [--method theirs Halpha]
      [--tol 1e-8] [--alpha 10 100 1000] [--alpha-bound 0.99 10 100]
      [--start S ...] [--max-iter-mine N] [--max-iter-theirs N]
      [--max-seconds T] [--out DIR] [--list]
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bourne_cases as bc  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    bc.add_common_args(ap, bc.APPROXIMATE, ["theirs", "Halpha"], "approximate")
    ap.add_argument("--alpha", nargs="*", type=float, default=[10.0, 100.0, 1000.0], help="stiffnesses of H_alpha")
    ap.add_argument("--alpha-bound", nargs="*", type=float, default=[0.99, 10.0, 100.0],
                    help="further stiffnesses, as multiples of the distinctness bound min|q_i - q_j| / (2 diam)")
    bc.run(ap.parse_args(), bc.APPROXIMATE, bc.approximate_runs)


if __name__ == "__main__":
    main()
