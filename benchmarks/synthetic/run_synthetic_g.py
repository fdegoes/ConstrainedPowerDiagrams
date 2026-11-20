#!/usr/bin/env python3
"""Runs the centroid dual G (max over the shifts, fixed weights; CentroidSolver,
LM Newton) on the synthetic benchmark.

The sites are the targets, q = cbar, and the fixed weights are those for
which the shifted weights of the solution are the weights wY of the
synthetic diagram. With q = cbar, eq:shifted-site-weight is
lambda_i = w_i + |u_i|^2, and the solution has u_i = Y_i - cbar_i, so
    w_i = wY_i - |Y_i - cbar_i|^2.
With fixed weights G has no homothety invariance: Y and wY are its unique
solution, and the errors are measured without removing one (the weights up
to their mean, which does not change the diagram).

The solve starts from zero shifts (the power diagram of (cbar, w)) and
stops on the shift residual of DualSolver (<= 1e-8 by default); the volumes
are not constrained by G, but match at the solution (volumeResidual,
max |V - Vbar| / Vbar measured on the final cells).

The start, the power diagram of (cbar, w), can have empty cells: they
have no gradient and the solver accepts no step that keeps them, so G
cannot run from there. benchBourne counts them exactly (startEmptyCells),
and each case gets an outcome:
  converged      G converged
  invalid start  the start has empty cells (a property of the data and
                 the fixed weights, not a failure of the solve)
  failed         valid start, but the solve did not converge
Each run writes
  <out>/<case>/oursG.json   benchBourne's result
and one row per case to <out>/oursG.csv.

Usage:
  external/bourne/.venv/bin/python benchmarks/synthetic/run_synthetic_g.py
      [--dim 2] [--n 10 25 ...] [--noise 0 1 2] [--volumes uniform voronoi]
      [--seeds 1 2 3 4 5] [--maxiter 100] [--out benchmarks/synthetic/results(3d)]
      [--force]
"""

import argparse
import csv
import json
import os
import subprocess
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import generate as G  # noqa: E402
import run_synthetic as RS  # noqa: E402

FIELDS = ["case", "n", "noise", "volumes", "seed", "outcome", "status", "converged", "iterations", "numBuilds", "seconds",
          "buildSeconds", "hessianSeconds", "factorSeconds", "volumeResidual", "shiftResidual", "generatorError",
          "weightError", "diagramMismatch", "minSpacing", "startEmptyCells", "load"]


def fixed_weights(c):
    """w = wY - |Y - cbar|^2: the shifted weights of the solution are wY."""
    return c["weights"] - ((c["truth"] - c["targets"]) ** 2).sum(1)


def write_case_file(path, c, maxiter, volume_tol, shift_tol):
    d = c["dim"]
    xy = " ".join(repr(float(x)) for x in c["targets"].ravel())
    lines = ["dim %d" % d, "box " + " ".join(["0"] * d + ["1"] * d), "n %d" % c["n"], "solver G",
             "volumetol %r" % volume_tol, "shifttol %r" % shift_tol, "maxiter %d" % maxiter, "notrajectory",
             "volumes " + " ".join(repr(float(x)) for x in c["volumes"]), "targets " + xy, "sites " + xy,
             "weights " + " ".join(repr(float(x)) for x in fixed_weights(c))]
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")


OUTCOMES = ("converged", "invalid start", "failed", "error")


def outcome(r):
    if r["converged"]:
        return "converged"
    return "invalid start" if r["startEmptyCells"] > 0 else "failed"


def run_case(c, args):
    out = os.path.join(c["dir"], "oursG%s.json" % args.tag)
    if args.force or not os.path.exists(out):
        cp = os.path.join(c["dir"], "caseG%s.txt" % args.tag)
        write_case_file(cp, c, args.maxiter, args.volumetol, args.shifttol)
        subprocess.run([RS.BENCH, cp, out], check=True)
    r = json.load(open(out))
    if "startEmptyCells" not in r:
        # A result of an older benchBourne: rerun to get the count.
        subprocess.run([RS.BENCH, os.path.join(c["dir"], "caseG%s.txt" % args.tag), out], check=True)
        r = json.load(open(out))
    # G has no volume gradient, so its volume residual is not computed by
    # the solver: measure it on the cells.
    r["volumeResidual"] = float(np.max(np.abs(np.asarray(r["volumes"]) - c["volumes"]) / c["volumes"]))
    h = RS.cell_size(c["n"], c["dim"])
    P = np.reshape(r["generators"], (-1, c["dim"]))
    lam = np.asarray(r["shiftedWeights"])
    gen = float(np.abs(P - c["truth"]).max() / h)
    dw = (lam - lam.mean()) - (c["weights"] - c["weights"].mean())
    return r, gen, float(np.abs(dw).max() / h ** 2), RS.diagram_mismatch(P, lam, c["truth"], c["weights"])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dim", type=int, default=2, choices=[2, 3])
    ap.add_argument("--n", nargs="+", type=int, default=list(RS.NS))
    ap.add_argument("--noise", nargs="+", type=int, default=list(G.NOISE))
    ap.add_argument("--volumes", nargs="+", default=list(G.VOLUMES))
    ap.add_argument("--seeds", nargs="+", type=int, default=[1, 2, 3, 4, 5])
    ap.add_argument("--maxiter", type=int, default=100, help="Newton iterations")
    ap.add_argument("--out", default=None, help="default: results (2D) or results3d (3D)")
    ap.add_argument("--volumetol", type=float, default=RS.VOLUME_TOL)
    ap.add_argument("--shifttol", type=float, default=RS.SHIFT_TOL)
    ap.add_argument("--force", action="store_true", help="rerun cases that already have a result")
    args = ap.parse_args()
    if args.out is None:
        args.out = os.path.join(HERE, "results" if args.dim == 2 else "results3d")
    default = args.volumetol == RS.VOLUME_TOL and args.shifttol == RS.SHIFT_TOL
    args.tag = "" if default else "-v%g-s%g" % (args.volumetol, args.shifttol)

    rows = []
    csv_path = os.path.join(args.out, "oursG%s.csv" % args.tag)
    os.makedirs(args.out, exist_ok=True)
    for n in args.n:
        for noise in args.noise:
            for vol in args.volumes:
                for seed in args.seeds:
                    c = G.make_case(args.out, n, noise, vol, seed, dim=args.dim)
                    load = os.getloadavg()[0]
                    t = time.time()
                    try:
                        r, gen, werr, mismatch = run_case(c, args)
                    except subprocess.CalledProcessError as e:
                        # A crash of benchBourne on one case: record it and go on.
                        rows.append({"case": c["name"], "n": n, "noise": noise, "volumes": vol, "seed": seed,
                                     "outcome": "error", "status": "exit %d" % e.returncode, "converged": False,
                                     "load": load})
                        print("%-26s benchBourne failed (exit %d)" % (c["name"], e.returncode), flush=True)
                        continue
                    row = {"case": c["name"], "n": n, "noise": noise, "volumes": vol, "seed": seed,
                           "outcome": outcome(r), "status": r["status"], "converged": r["converged"], "iterations": r["iterations"],
                           "numBuilds": r["numBuilds"], "seconds": r["seconds"],
                           "buildSeconds": r["buildSeconds"], "hessianSeconds": r["hessianSeconds"],
                           "factorSeconds": r["factorSeconds"], "volumeResidual": r["volumeResidual"],
                           "shiftResidual": r["shiftResidual"], "generatorError": gen, "weightError": werr,
                           "diagramMismatch": mismatch, "minSpacing": RS.min_spacing(c["truth"], n),
                           "startEmptyCells": r["startEmptyCells"], "load": load}
                    rows.append(row)
                    print("%-26s %-13s %-15s %5d it %8.3fs  res %.1e/%.1e  gen %.1e  w %.1e  empty0 %d  (%.0fs)" % (
                        c["name"], row["outcome"], r["status"], r["iterations"], r["seconds"], r["volumeResidual"],
                        r["shiftResidual"], gen, werr, row["startEmptyCells"], time.time() - t), flush=True)
                    with open(csv_path, "w", newline="") as f:
                        w = csv.DictWriter(f, FIELDS)
                        w.writeheader()
                        w.writerows(rows)
    count = {o: sum(r["outcome"] == o for r in rows) for o in OUTCOMES}
    valid = len(rows) - count["invalid start"] - count["error"]
    print("converged %d / %d with a valid start (%d invalid starts, %d failed, %d errors); summary in %s" % (
        count["converged"], valid, count["invalid start"], count["failed"], count["error"], csv_path))


if __name__ == "__main__":
    main()
