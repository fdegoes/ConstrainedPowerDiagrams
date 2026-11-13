#!/usr/bin/env python3
"""Generates every case of the synthetic benchmark (generate.py) and reports
the statistics of the data, per noise level and volume distribution (pooled
over n and seeds) and per n (pooled over noise, volumes and seeds):
  nn min, nn cv    the distance of each generator to its nearest neighbour, in
                   cell sizes h = (1/n)^(1/dim): the minimum, and the coefficient
                   of variation (std / mean) per case
  V max/min        the ratio of the largest to the smallest cell volume
  |w| / h^2        the largest weight of the solution (weights shifted to zero
                   mean), in squared cell sizes
  |Y - cbar| / h   the distance from the start q = cbar to the solution Y: the
                   mean and the maximum over the cells
Ranges are min-max over the cases.

Usage:
  external/bourne/.venv/bin/python benchmarks/synthetic/stats.py
      [--dim 2] [--out benchmarks/synthetic/results(3d)] [--md stats.md]
"""

import argparse
import math
import os
import sys

import numpy as np
from scipy.spatial import cKDTree

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import generate as G  # noqa: E402
from run_synthetic import NS  # noqa: E402

SEEDS = (1, 2, 3, 4, 5)


def case_stats(c):
    n = c["n"]
    h = (1.0 / n) ** (1.0 / c["dim"])
    Y, cbar, V, w = c["truth"], c["targets"], c["volumes"], c["weights"]
    nn = cKDTree(Y).query(Y, k=2)[0][:, 1] / h
    d = np.linalg.norm(Y - cbar, axis=1) / h
    return {"nnMin": nn.min(), "nnCv": nn.std() / nn.mean(), "vRatio": V.max() / V.min(),
            "w": np.abs(w - w.mean()).max() / h ** 2, "dMean": d.mean(), "dMax": d.max()}


def rng(rows, key, fmt):
    v = [r[key] for r in rows]
    return (fmt + " - " + fmt) % (min(v), max(v))


def table(groups, label):
    lines = ["| %s | cases | nn min | nn cv | V max/min | \\|w\\| / h^2 | mean \\|Y - cbar\\| / h | max \\|Y - cbar\\| / h |"
             % label, "|---|---|---|---|---|---|---|---|"]
    for name, rows in groups:
        lines.append("| %s | %d | %.3f | %s | %s | %s | %s | %s |" % (
            name, len(rows), min(r["nnMin"] for r in rows), rng(rows, "nnCv", "%.2f"), rng(rows, "vRatio", "%.2f"),
            rng(rows, "w", "%.3f"), rng(rows, "dMean", "%.3f"), rng(rows, "dMax", "%.3f")))
    return "\n".join(lines)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dim", type=int, default=2, choices=[2, 3])
    ap.add_argument("--out", default=None, help="default: results (2D) or results3d (3D)")
    ap.add_argument("--md", default=None, help="also write the tables to this markdown file")
    args = ap.parse_args()
    if args.out is None:
        args.out = os.path.join(HERE, "results" if args.dim == 2 else "results3d")
    rows = []
    for n in NS:
        for noise in G.NOISE:
            for vol in G.VOLUMES:
                for seed in SEEDS:
                    r = case_stats(G.make_case(args.out, n, noise, vol, seed, dim=args.dim))
                    r.update(n=n, noise=noise, vol=vol)
                    rows.append(r)
        print("n = %d done" % n, file=sys.stderr)
    by_level = [("%d (rho %.2f), %s" % (k, G.NOISE_SPACING[k], v), [r for r in rows if r["noise"] == k and r["vol"] == v])
                for k in G.NOISE for v in G.VOLUMES]
    by_n = [("%d" % n, [r for r in rows if r["n"] == n]) for n in NS]
    text = "%d cases\n\n%s\n\n%s\n" % (len(rows), table(by_level, "noise, volumes"), table(by_n, "n"))
    print(text)
    if args.md:
        open(args.md, "w").write(text)


if __name__ == "__main__":
    main()
