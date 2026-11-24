#!/usr/bin/env python3
"""A power diagram with scattered prescribed centroids: target centroids cbar placed by
dart throwing in the unit square (kept half a dart radius from its boundary), volumes as
close to uniform as the conservation identities eq:partition-feasibility allow (the
max-entropy volumes of latex/appendix-maxent.tex, eq:maxent-dual), and H (or H_alpha with
--alpha) solved from sites at the targets, q = cbar (introFigure).
  latex/figures/intro-maxent.pdf
Needs numpy, matplotlib and the introFigure tool (ninja -C build introFigure):
  external/bourne/.venv/bin/python benchmarks/figures/plot_maxent.py [--seed 1] [--alpha -0.99]
"""

import argparse
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import plot_intro as PI  # noqa: E402  (sets MPLCONFIGDIR and the Agg backend)

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.patches import Polygon  # noqa: E402

GL = PI.GL
OUT = os.path.join(PI.ROOT, "latex", "figures", "intro-maxent.pdf")


def dart_throwing(rng, radius, rmax, tries=1500000):
    """Points of the unit square at least max(radius(y), radius(z)) from each other
    (Poisson-disk sampling with a varying radius), with a grid of cells of size rmax."""
    m = int(np.ceil(1.0 / rmax))
    grid = [[[] for _ in range(m)] for _ in range(m)]
    Y, R = [], []
    for _ in range(tries):
        y = rng.random(2)
        r = radius(y)
        gi, gj = min(int(y[0] * m), m - 1), min(int(y[1] * m), m - 1)
        ok = True
        for a in range(max(gi - 1, 0), min(gi + 2, m)):
            for b in range(max(gj - 1, 0), min(gj + 2, m)):
                for k in grid[a][b]:
                    if np.hypot(*(y - Y[k])) < max(r, R[k]):
                        ok = False
                        break
                if not ok:
                    break
            if not ok:
                break
        if ok:
            grid[gi][gj].append(len(Y))
            Y.append(y)
            R.append(r)
    return np.array(Y), np.array(R)


def maxent_volumes(cbar, area, centroid, tol=1e-14, maxiter=100):
    """The volumes V_i = |Omega| exp(-y . l_i) with the lifts l_i = (cbar_i, 1), for the
    minimizer y of |Omega| sum_i exp(-y . l_i) + y . l_Omega, l_Omega = |Omega| (c_Omega, 1)
    (eq:maxent-dual), by Newton's method: then sum_i V_i l_i = l_Omega."""
    L = np.column_stack([cbar, np.ones(len(cbar))])
    lO = area * np.append(centroid, 1.0)
    y = np.append(np.zeros(cbar.shape[1]), np.log(len(cbar)))  # uniform volumes area / n
    for _ in range(maxiter):
        V = area * np.exp(-L @ y)
        g = lO - V @ L
        if np.linalg.norm(g) <= tol * area:
            return V
        y -= np.linalg.solve((L * V[:, None]).T @ L, g)   # Newton: the gradient is g, the Hessian sum_i V_i l_i l_i^T
    raise RuntimeError("the max-entropy volumes did not converge")


def run_tool(cbar, V, tool, alpha, centroid_only=0, log=False):
    """introFigure with given targets and volumes, sites at the targets."""
    n = len(cbar)
    import subprocess
    import tempfile
    with tempfile.TemporaryDirectory() as tmp:
        src, dst = os.path.join(tmp, "in.txt"), os.path.join(tmp, "out.txt")
        with open(src, "w") as f:
            f.write("n %d\n" % n)
            f.write("targets " + " ".join("%.17g" % x for x in cbar.ravel()) + "\n")
            f.write("volumes " + " ".join("%.17g" % x for x in V) + "\n")
            f.write("sitesAtTargets 1\nclosest 1\n")
            if alpha:
                f.write("alpha %.17g\n" % alpha)
            if log:
                f.write("log 1\n")
            if centroid_only:
                f.write("centroidOnly %d\n" % centroid_only)
        run = subprocess.run([tool, src, dst], capture_output=True, text=True)
        if run.returncode != 0:
            return None, (run.stdout + run.stderr).strip()
        rec = {}
        with open(dst) as f:
            for line in f:
                key, *vals = line.split()
                rec[key] = np.array([float(v) for v in vals])
    for k in ("targets", "sites", "generators", "shifts"):
        rec[k] = rec[k].reshape(n, 2)
    return rec, run.stdout.strip()


def draw(rec, cols, out, dpi):
    cbar, p, lam = rec["targets"], rec["generators"], rec["weights"]
    cells = GL.power_cells(p, lam)
    fig = plt.figure(figsize=(2.0, 2.0))
    ax = fig.add_axes([0, 0, 1, 1])
    for P, c in zip(cells, cols):
        ax.add_patch(Polygon(P, closed=True, facecolor=GL.tint(c, GL.CELL_TINT), edgecolor="black",
                             linewidth=0.6 * PI.LW, joinstyle="round"))
    for i in range(len(p)):                             # the shifts, from the centroids to the generators
        ax.plot(*zip(cbar[i], p[i]), color="black", linewidth=0.4, solid_capstyle="round", zorder=3)
    ax.plot(cbar[:, 0], cbar[:, 1], "s", markerfacecolor="white", markeredgecolor="black", markersize=1.6,
            markeredgewidth=0.35, zorder=4)
    ax.plot(p[:, 0], p[:, 1], "o", color="black", markersize=1.3, markeredgewidth=0, zorder=5, clip_on=False)
    ax.add_patch(Polygon([(0, 0), (1, 0), (1, 1), (0, 1)], closed=True, fill=False, edgecolor="black",
                         linewidth=1.2 * PI.LW))
    pad = 0.006
    ax.set_xlim(-pad, 1 + pad)
    ax.set_ylim(-pad, 1 + pad)
    ax.set_aspect("equal")
    ax.axis("off")
    fig.savefig(out, dpi=dpi)
    plt.close(fig)
    return cells


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--r", type=float, default=0.07, help="dart radius")
    ap.add_argument("--margin", type=float, default=0.5, help="distance to the boundary, in dart radii")
    ap.add_argument("--shrink", type=float, default=1.0, help="scale of the targets about the centre")
    ap.add_argument("--alpha", type=float, default=0.0,
                    help="H_alpha with this stiffness (< 0: that fraction of the practical bound); 0: H")
    ap.add_argument("--g", type=int, default=0, choices=(0, 1, 2),
                    help="G instead of H: 1 with w = 0, 2 with w from F at the sites")
    ap.add_argument("--log", action="store_true", help="print the solver's iterates")
    ap.add_argument("--tool", default=os.path.join(PI.ROOT, "build", "introFigure"))
    ap.add_argument("--out", default=OUT)
    ap.add_argument("--dpi", type=int, default=300, help="for raster previews")
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)
    m = args.margin * args.r
    Y, _ = dart_throwing(rng, lambda y: args.r, args.r)    # in the unit square, mapped to the inset one
    cbar = m + (1 - 2 * m) * Y
    cbar = 0.5 + args.shrink * (cbar - 0.5)
    n = len(cbar)
    V = maxent_volumes(cbar, 1.0, np.array([0.5, 0.5]))
    print("n = %d; max-entropy volumes / (1/n): %.3f to %.3f; mean of the samples - c_Omega: (%.1e, %.1e)" % (
        n, V.min() * n, V.max() * n, *(cbar.mean(0) - 0.5)))
    rec, msg = run_tool(cbar, V, args.tool, args.alpha, args.g, args.log)
    print(msg)
    if rec is None:
        raise SystemExit(1)
    cells = draw(rec, GL.colors(n, args.seed), args.out, args.dpi)
    A = np.array([PI.area(P) for P in cells])
    C = np.array([PI.centroid(P) for P in cells])
    r = np.sqrt(V)
    print("drawn cells: max |V - Vbar| / Vbar %.1e, |c - cbar| / r_i: max %.1e median %.1e" % (
        np.max(np.abs(A - V) / V), np.max(np.linalg.norm(C - cbar, axis=1) / r),
        np.median(np.linalg.norm(C - cbar, axis=1) / r)))
    print("wrote " + args.out)


if __name__ == "__main__":
    main()
