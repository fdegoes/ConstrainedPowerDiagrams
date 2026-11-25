#!/usr/bin/env python3
"""Summarizes runPP2D runs: one row per run directory (energy.csv, run.txt, frame CSVs if any).

Columns:
  steps     steps completed (frames of energy.csv - 1)
  K/K0      kinetic energy at the last step over the first, and its min and max over the run
  dE/E0     max relative change of the total energy (energy.csv's total: K + the integrator's
            potential; for the constrained integrator kappa (Q + R))
  volErr    max relative volume error over the steps
  cGap      max centroid gap |c - q| over the steps, in mean spacings h
  its       median and max Newton iterations per step
  pair      closest particle pair at the last saved frame, in h
  |u|       max and rms generator offset |p - q| at the last saved frame, in h

Usage:
  external/bourne/.venv/bin/python experiments/pp2d/summarize_runs.py <run_dir> [<run_dir> ...]
"""

import csv
import glob
import os
import statistics
import sys

import numpy as np
from scipy.spatial import cKDTree


def summarize(d):
    rows = list(csv.DictReader(open(os.path.join(d, "energy.csv"))))
    col = lambda k: np.array([float(r[k]) for r in rows])
    K, E = col("kinetic"), col("total")
    n = None
    out = dict(run=os.path.relpath(d), steps=len(rows) - 1)
    out["K"] = "%.3f [%.3f, %.3f]" % (K[-1] / K[0], K.min() / K[0], K.max() / K[0])
    out["dE"] = "%.1e" % (np.abs(E - E[0]).max() / abs(E[0]))
    out["vol"] = "%.1e" % col("maxVolumeError").max()
    its = [int(r["iterations"]) for r in rows[1:]] or [0]
    out["its"] = "%d / %d" % (statistics.median(its), max(its))
    frames = sorted(glob.glob(os.path.join(d, "frame_*.csv")))
    out["pair"] = out["u"] = "-"
    if frames:
        R = list(csv.DictReader(open(frames[-1])))
        a = lambda k: np.array([float(r[k]) for r in R])
        q = np.c_[a("x"), a("y")]
        n = len(q)
        h = 1 / np.sqrt(n)
        dd, _ = cKDTree(q).query(q, k=2)
        U = np.linalg.norm(np.c_[a("ux"), a("uy")], axis=1)
        out["pair"] = "%.2f" % (dd[:, 1].min() / h)
        out["u"] = "%.2f / %.3f" % (U.max() / h, np.sqrt((U ** 2).mean()) / h)
    h = 1 / np.sqrt(n or 5000)
    out["gap"] = "%.1e" % (col("maxCentroidGap").max() / h)
    status = ""
    if os.path.exists(os.path.join(d, "status.txt")):
        status = open(os.path.join(d, "status.txt")).read().strip()
    out["status"] = status
    return out


def main():
    keys = ["run", "steps", "K", "dE", "vol", "gap", "its", "pair", "u", "status"]
    print("| " + " | ".join(["run", "steps", "K/K0 [min, max]", "dE/E0", "volErr", "cGap (h)", "its med/max",
                             "pair (h)", "|u| max/rms (h)", "status"]) + " |")
    print("|" + "---|" * len(keys))
    for d in sys.argv[1:]:
        if not os.path.exists(os.path.join(d, "energy.csv")):
            print("| %s | no energy.csv |" % d)
            continue
        s = summarize(d)
        print("| " + " | ".join(str(s[k]) for k in keys) + " |")


if __name__ == "__main__":
    main()
