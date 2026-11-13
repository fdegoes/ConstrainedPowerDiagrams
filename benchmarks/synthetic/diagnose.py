#!/usr/bin/env python3
"""Why a synthetic case does not converge: reruns it with the trajectory and
reports, per attempt (the start and each JointSolver restart),
  - how far the residuals got and where they stalled,
  - the closest generator pair of the final iterate, its distance at the
    start and in the truth Y (cell sizes), and the volumes of the two cells,
  - the cells in the wrong place (the diagram mismatch, by cell),
and classifies the failure:
  collapse    a pair of generators ends at < 1e-2 of their distance in Y
  stall       no collapse; the residuals stop decreasing
Usage:
  external/bourne/.venv/bin/python benchmarks/synthetic/diagnose.py <case> [<case> ...]
      [--dim 2] [--volumetol 1e-6 --shifttol 1e-6] [--results benchmarks/synthetic/results(3d)]
"""

import argparse
import json
import math
import os
import subprocess
import sys
import tempfile

import numpy as np
from scipy.spatial import cKDTree

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import generate as G  # noqa: E402
from run_synthetic import BENCH, write_case_file  # noqa: E402


def attempts(T):
    """Split the trajectory at the restarts (the value drops back to the
    start's)."""
    starts = [0] + [k for k in range(1, len(T)) if T[k]["value"] < T[k - 1]["value"] - 1e-12 * abs(T[0]["value"])
                    and abs(T[k]["value"] - T[0]["value"]) <= 1e-12 * abs(T[0]["value"])]
    return [T[a:b] for a, b in zip(starts, starts[1:] + [len(T)])]


def closest_pair(P):
    d, j = cKDTree(P).query(P, k=2)
    i = int(np.argmin(d[:, 1]))
    return i, int(j[i, 1]), float(d[i, 1])


def diagnose(name, args):
    level, vol, n, seed = name.split("-")
    noise, n, seed = int(level[len("noise"):]), int(n), int(seed)
    c = G.make_case(args.results, n, noise, vol, seed, dim=args.dim)
    Y, V = c["truth"], c["volumes"]
    h = (1.0 / n) ** (1.0 / c["dim"])
    with tempfile.TemporaryDirectory() as d:
        cp = os.path.join(d, "c.txt")
        write_case_file(cp, c, args.maxiter, args.volumetol, args.shifttol)
        lines = [l for l in open(cp) if l.strip() != "notrajectory"]
        open(cp, "w").writelines(lines)
        subprocess.run([BENCH, cp, os.path.join(d, "o.json")], check=True)
        r = json.load(open(os.path.join(d, "o.json")))
    T = r["trajectory"]
    print("%s: %s after %d iterations in %d attempts (%.2fs)" % (name, r["status"], r["iterations"],
                                                                  r.get("attempts", 1), r["seconds"]))
    dY = cKDTree(Y).query(Y, k=2)[0][:, 1]
    print("  truth: min generator spacing %.2f h, cell volumes %.2f .. %.2f of the mean" % (
        dY.min() / h, V.min() / V.mean(), V.max() / V.mean()))
    kinds = []
    for a, A in enumerate(attempts(T)):
        res = np.array([max(t["volumeResidual"] / args.volumetol, t["shiftResidual"] / args.shifttol) for t in A])
        best = int(np.argmin(res))
        P0 = np.reshape(A[0]["generators"], (-1, c["dim"]))
        P = np.reshape(A[-1]["generators"], (-1, c["dim"]))
        i, j, dmin = closest_pair(P)
        dTrue = float(np.linalg.norm(Y[i] - Y[j]))
        d0 = float(np.linalg.norm(P0[i] - P0[j]))
        # First iterate where the pair is below 10% of its true distance.
        k_close = next((k for k, t in enumerate(A) if np.linalg.norm(
            np.reshape(t["generators"], (-1, c["dim"]))[i] - np.reshape(t["generators"], (-1, c["dim"]))[j]) < 0.1 * dTrue), None)
        kind = "collapse" if dmin < 1e-2 * dTrue else "stall"
        kinds.append(kind)
        print("  attempt %d: %4d iterations, best residual/tol %.1e at iteration %d, final volume %.1e shift %.1e"
              % (a + 1, len(A) - 1, res[best], best, A[-1]["volumeResidual"], A[-1]["shiftResidual"]))
        print("    closest pair (%d, %d): %.1e h at the end, %.2f h at the start, %.2f h in the truth;"
              " cells %.2f and %.2f of the mean volume; below 10%% of the truth from iteration %s  -> %s"
              % (i, j, dmin / h, d0 / h, dTrue / h, V[i] / V.mean(), V[j] / V.mean(), k_close, kind))
    return kinds


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cases", nargs="+")
    ap.add_argument("--volumetol", type=float, default=1e-6)
    ap.add_argument("--shifttol", type=float, default=1e-6)
    ap.add_argument("--maxiter", type=int, default=300)
    ap.add_argument("--dim", type=int, default=2, choices=[2, 3])
    ap.add_argument("--results", default=None, help="default: results (2D) or results3d (3D)")
    args = ap.parse_args()
    if args.results is None:
        args.results = os.path.join(HERE, "results" if args.dim == 2 else "results3d")
    for name in args.cases:
        diagnose(name, args)


if __name__ == "__main__":
    main()
