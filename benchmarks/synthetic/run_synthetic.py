#!/usr/bin/env python3
"""Runs our solver (max H, LM Newton) on the synthetic benchmark.

For every case (generate.py) the solve starts from the sites q = cbar
(MakeJointStart with s = 0) and stops on the common residual rule: volume
residual <= 1e-8 and shift residual <= 1e-8 (max-norm per cell, see
DualSolver). Each run writes
  <out>/<case>/ours.json   benchBourne's result (status, iterations, restarts,
                           timers, residuals, generators, weights, shifts,
                           volumes, centroids)
and one row per case to <out>/ours.csv.

Runs are sequential by default: the solver is multithreaded, and concurrent
runs would distort the timings.

The 3D cases (--dim 3) are those of results3d/, in the unit cube.

Usage:
  external/bourne/.venv/bin/python benchmarks/synthetic/run_synthetic.py
      [--dim 2] [--n 10 25 ...] [--noise 0 1 2] [--volumes uniform voronoi]
      [--seeds 1 2 3 4 5] [--maxiter 100] [--out benchmarks/synthetic/results(3d)]
      [--force] [--case-extra 'stretchlimit 0' --label nolimit]
"""

import argparse
import csv
import json
import math
import os
import subprocess
import sys
import tempfile
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import generate as G  # noqa: E402
import their_rules as TR  # noqa: E402

BENCH = os.path.join(G.ROOT, "build", "benchBourne")
NS = (10, 14, 21, 30, 43, 62, 89, 127, 183, 264, 379, 546, 785, 1129, 1624, 2336, 3360, 4833, 6952, 10000)  # 20 sizes, log-spaced: round(10 * 1000^(k / 19))
VOLUME_TOL = 1e-8
SHIFT_TOL = 1e-8

FIELDS = ["case", "n", "noise", "volumes", "seed", "status", "converged", "iterations", "attempts",
          "attemptIterations", "numBuilds", "seconds", "startSeconds", "buildSeconds", "hessianSeconds",
          "factorSeconds", "volumeResidual", "shiftResidual", "generatorError", "diagramMismatch", "minSpacing", "load",
          "rule1Ratio", "rule1Met"]


def write_case_file(path, c, maxiter, volume_tol=VOLUME_TOL, shift_tol=SHIFT_TOL, trajectory="notrajectory", extra=()):
    d = c["dim"]
    xy = " ".join(repr(float(x)) for x in c["targets"].ravel())
    lines = ["dim %d" % d, "box " + " ".join(["0"] * d + ["1"] * d), "n %d" % c["n"], "solver H", "alpha 1e300", "start 0",
             "volumetol %r" % volume_tol, "shifttol %r" % shift_tol, "maxiter %d" % maxiter, trajectory,
             "volumes " + " ".join(repr(float(x)) for x in c["volumes"]), "targets " + xy, "sites " + xy]
    lines += list(extra)
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")


def cell_size(n, d):
    """h = (1/n)^(1/d), the mean cell size in the unit square or cube."""
    return (1.0 / n) ** (1.0 / d)


def generator_error(P, Y, n):
    """max_i |p_i - Y_i| / h with the homothety about the centroid removed
    (H is invariant along it)."""
    Pc, Yc = P - P.mean(0), Y - Y.mean(0)
    a = float((Pc * Yc).sum() / (Yc * Yc).sum())
    return float(np.abs(Pc - a * Yc).max() / cell_size(n, Y.shape[1]))


def diagram_mismatch(P, w, Y, wY, m=1000, samples=10 ** 6):
    """Fraction of the unit square (an m x m grid of cell centres) or cube
    (samples uniform random points) assigned to a different cell by the power
    diagram (P, w) than by the true one (Y, wY): 0 when the cells are
    recovered, whichever generators represent them."""
    from scipy.spatial import cKDTree
    if Y.shape[1] == 2:
        g = (np.arange(m) + 0.5) / m
        X = np.stack(np.meshgrid(g, g), -1).reshape(-1, 2)
    else:
        X = np.random.default_rng(0).random((samples, Y.shape[1]))
    Xl = np.c_[X, np.zeros(len(X))]

    def labels(Q, wq):
        # Nearest lifted site (q, sqrt(W - w)): the power cell of x.
        return cKDTree(np.c_[Q, np.sqrt(wq.max() - wq)]).query(Xl)[1]

    return float(np.mean(labels(P, np.asarray(w)) != labels(Y, np.asarray(wY))))


def min_spacing(Y, n):
    from scipy.spatial import cKDTree
    d, _ = cKDTree(Y).query(Y, k=2)
    return float(d[:, 1].min() / cell_size(n, Y.shape[1]))



def _trim_to_solution(traj, generators):
    """The trajectory up to the iterate returned by the solver. With restarts,
    the trajectory chains all the attempts, each from the same start, and the
    returned solution (the converged attempt, or else the one with the
    smallest residual) need not be the last iterate."""
    if not traj:
        return traj
    g = np.asarray(generators, dtype=float)
    for k in range(len(traj) - 1, -1, -1):
        p = np.asarray(traj[k]["generators"], dtype=float)
        if p.shape == g.shape and np.array_equal(p, g):
            return traj[:k + 1]
    return traj

def run_case(c, args):
    if args.no_save:
        # Status and counts only: the case and the result live in a temporary
        # directory, deleted afterwards, and the costly measures are skipped.
        with tempfile.TemporaryDirectory() as d:
            cp, out = os.path.join(d, "case.txt"), os.path.join(d, "ours.json")
            write_case_file(cp, c, args.maxiter, args.volumetol, args.shifttol, "notrajectory", args.case_extra)
            subprocess.run([BENCH, cp, out], check=True)
            r = json.load(open(out))
        r["rule1"] = None
        return r, float("nan"), float("nan")
    out = os.path.join(c["dir"], "ours%s.json" % args.tag)
    if args.force or not os.path.exists(out):
        cp = os.path.join(c["dir"], "case%s.txt" % args.tag)
        write_case_file(cp, c, args.maxiter, args.volumetol, args.shifttol, "lighttrajectory", args.case_extra)
        subprocess.run([BENCH, cp, out], check=True)
    r = json.load(open(out))
    # Rule 1 of Bourne et al. (2025, their_rules.py) measured at our
    # solution; the trajectory is then dropped from memory.
    traj = r.pop("trajectory", None)
    if traj:
        traj = _trim_to_solution(traj, r["generators"])
    r["rule1"] = TR.slsqp_at_solution(c["targets"], c["volumes"], traj) if traj else None
    P = np.reshape(r["generators"], (-1, c["dim"]))
    return r, generator_error(P, c["truth"], c["n"]), diagram_mismatch(P, r["shiftedWeights"], c["truth"],
                                                                         c["weights"])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dim", type=int, default=2, choices=[2, 3])
    ap.add_argument("--n", nargs="+", type=int, default=list(NS))
    ap.add_argument("--noise", nargs="+", type=int, default=list(G.NOISE))
    ap.add_argument("--volumes", nargs="+", default=list(G.VOLUMES))
    ap.add_argument("--seeds", nargs="+", type=int, default=[1, 2, 3, 4, 5])
    ap.add_argument("--maxiter", type=int, default=100, help="Newton iterations per attempt")
    ap.add_argument("--out", default=None, help="default: results (2D) or results3d (3D)")
    ap.add_argument("--volumetol", type=float, default=VOLUME_TOL)
    ap.add_argument("--shifttol", type=float, default=SHIFT_TOL)
    ap.add_argument("--force", action="store_true", help="rerun cases that already have a result")
    ap.add_argument("--case-extra", nargs="*", default=[],
                    help="extra benchBourne case lines, e.g. 'stretchlimit 0' 'pairreg 0.3 2' (use with --label)")
    ap.add_argument("--label", default="", help="suffix of the result files (ours-<label>.json, .csv)")
    ap.add_argument("--no-save", action="store_true",
                    help="keep no per-case files (temporary directory): status, iterations and counts in the CSV only")
    args = ap.parse_args()
    if args.out is None:
        args.out = os.path.join(HERE, "results" if args.dim == 2 else "results3d")
    # Results of non-default tolerances get their own files.
    default = args.volumetol == VOLUME_TOL and args.shifttol == SHIFT_TOL
    args.tag = "" if default else "-v%g-s%g" % (args.volumetol, args.shifttol)
    if args.label:
        args.tag += "-" + args.label

    rows = []
    csv_path = os.path.join(args.out, "ours%s.csv" % args.tag)
    os.makedirs(args.out, exist_ok=True)
    for n in args.n:
        for noise in args.noise:
            for vol in args.volumes:
                for seed in args.seeds:
                    c = G.make_case(args.out, n, noise, vol, seed, dim=args.dim)
                    load = os.getloadavg()[0]
                    t = time.time()
                    try:
                        r, err, mismatch = run_case(c, args)
                    except subprocess.CalledProcessError as e:
                        # A crash of benchBourne on one case: record it and go on.
                        rows.append({"case": c["name"], "n": n, "noise": noise, "volumes": vol, "seed": seed,
                                     "status": "error (exit %d)" % e.returncode, "converged": False, "load": load})
                        print("%-26s benchBourne failed (exit %d)" % (c["name"], e.returncode), flush=True)
                        continue
                    row = {"case": c["name"], "n": n, "noise": noise, "volumes": vol, "seed": seed,
                           "status": r["status"], "converged": r["converged"], "iterations": r["iterations"],
                           "attempts": r.get("attempts", 1), "attemptIterations": r.get("attemptIterations"),
                           "numBuilds": r["numBuilds"], "seconds": r["seconds"],
                           "startSeconds": r["startSeconds"], "buildSeconds": r["buildSeconds"],
                           "hessianSeconds": r["hessianSeconds"], "factorSeconds": r["factorSeconds"],
                           "volumeResidual": r["volumeResidual"], "shiftResidual": r["shiftResidual"],
                           "generatorError": err, "diagramMismatch": mismatch, "minSpacing": min_spacing(c["truth"], n), "load": load}
                    if r["rule1"]:
                        row.update(rule1Ratio=r["rule1"]["ratio"], rule1Met=r["rule1"]["met"])
                    rows.append(row)
                    print("%-26s %-15s %5d it %d att %8.3fs  res %.1e/%.1e  mismatch %.1e  load %.1f  (%.0fs)" % (
                        c["name"], r["status"], r["iterations"], row["attempts"], r["seconds"],
                        r["volumeResidual"], r["shiftResidual"], mismatch, load, time.time() - t), flush=True)
                    with open(csv_path, "w", newline="") as f:
                        w = csv.DictWriter(f, FIELDS)
                        w.writeheader()
                        w.writerows(rows)
    ok = sum(r["converged"] for r in rows)
    print("converged %d / %d; summary in %s" % (ok, len(rows), csv_path))
    failed = [r for r in rows if not r["converged"]]
    if failed:
        # How far the failed solves are: the residuals of the returned
        # (best) attempt, the volume error relative to Vbar and the centroid
        # error in cell sizes.
        print("failed: %d; residuals of the best attempt (volume / centroid in cell sizes):" % len(failed))
        for r in failed:
            print("  %-26s %.1e / %.1e" % (r["case"], r["volumeResidual"], r["shiftResidual"]))


if __name__ == "__main__":
    main()
