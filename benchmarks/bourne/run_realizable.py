#!/usr/bin/env python3
"""Benchmark on the realizable cases of Bourne et al. (an exact power diagram
has the prescribed volumes and centroids): 2025-invert, 2026-voronoi and
2026-aligned. Runs their method (SLSQP on their H, 2025; Frank-Wolfe, 2026)
and our max H, under the SLSQP rule of 2025 and the common residual of H
(bourne_cases.py). H_alpha runs only on the approximate cases.

Usage (all cases, methods, rules and tolerances by default):
  external/bourne/.venv/bin/python benchmarks/bourne/run_realizable.py
      [--case NAME ...] [--method theirs H]
      [--tol 1e-8] [--start S ...] [--max-iter-mine N]
      [--max-iter-theirs N] [--max-seconds T] [--out DIR] [--list]
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bourne_cases as bc  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    bc.add_common_args(ap, bc.REALIZABLE, ["theirs", "H"], "realizable")
    bc.run(ap.parse_args(), bc.REALIZABLE, bc.realizable_runs)


if __name__ == "__main__":
    main()
