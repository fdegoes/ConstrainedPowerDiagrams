#!/usr/bin/env python3
"""Draws the power diagram of every case of the synthetic benchmark (the
known solution Y, w), for browsing the data; the paper uses gallery.py. One
directory per noise level and volume distribution:
  <out>/noise<k>-<volumes>/n<n>-seed<seed>.<format>
The style is that of gallery.py: cells in a light tint of a random colour,
generators as dots inside a schematic disk of radius
h (R0 + R1 (w_i - mean w) / max_j |w_j - mean w|), h = sqrt(1/n) the cell
size, so the disks show the variation of the weights within a case (equal
disks: equal weights). Line widths and dots shrink with n.

Each cell is clipped to the unit square and against its neighbours in the
regular triangulation of (Y, w) (the lower hull of the lifted points
(Y_i, |Y_i|^2 - w_i)), which contain its neighbours in the clipped diagram.
The cell areas are checked against the stored volumes.

Usage:
  external/bourne/.venv/bin/python benchmarks/synthetic/render_all.py
      [--n N ...] [--noise 0 1 2] [--volumes uniform voronoi]
      [--seeds 1 2 3 4 5] [--format pdf] [--out benchmarks/synthetic/images]
      [--no-disks] [--force]
"""

import argparse
import os
import sys
import time

os.environ.setdefault("MPLCONFIGDIR", "/tmp/mplcfg")
import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.collections import PatchCollection, PolyCollection  # noqa: E402
from matplotlib.patches import Circle, Polygon  # noqa: E402
from scipy.spatial import ConvexHull  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import gallery as GL  # noqa: E402
import generate as G  # noqa: E402
from run_synthetic import NS  # noqa: E402

DISK_R0 = 0.13  # disk radius of the mean weight, in cell sizes h (gallery.py at n = 10)
DISK_R1 = 0.10  # radius change of the largest |weight| of the case, in cell sizes h
SIZE = 4.0  # inches


def neighbours(Y, w):
    """Neighbours of each site in the regular triangulation of (Y, w)."""
    L = np.c_[Y, (Y * Y).sum(1) - w]
    hull = ConvexHull(L)
    nb = [set() for _ in range(len(Y))]
    for simplex, eq in zip(hull.simplices, hull.equations):
        if eq[2] < 0:  # lower facet (outward normal pointing down)
            for i in simplex:
                nb[i].update(int(j) for j in simplex if j != i)
    return nb


def power_cells(Y, w):
    """The cells of the power diagram of (Y, w) clipped to the unit square."""
    nb = neighbours(Y, w)
    cells = []
    for i in range(len(Y)):
        poly = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]
        for j in nb[i]:
            if not poly:
                break
            poly = GL.clip(poly, 2.0 * (Y[j] - Y[i]), Y[j] @ Y[j] - w[j] - Y[i] @ Y[i] + w[i])
        cells.append(poly)
    return cells


def draw(path, cells, Y, w, cols, disks):
    n = len(Y)
    h = np.sqrt(1.0 / n)
    k = min(1.0, np.sqrt(30.0 / n))  # thinner lines and smaller dots for many cells
    fig = plt.figure(figsize=(SIZE, SIZE))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.add_collection(PolyCollection([p for p in cells if p], facecolors=[GL.tint(c, GL.CELL_TINT) for p, c in
                                                                          zip(cells, cols) if p],
                                     edgecolors="black", linewidths=0.6 * k, joinstyle="round"))
    if disks:
        dw = w - w.mean()
        s = np.abs(dw).max()
        radii = h * (DISK_R0 + DISK_R1 * dw / s) if s > 0 else np.full(n, h * DISK_R0)
        order = np.argsort(-radii)
        pc = PatchCollection([Circle(tuple(Y[i]), radii[i]) for i in order],
                             facecolors=[(*GL.tint(cols[i], GL.DISK_TINT), GL.DISK_ALPHA) for i in order],
                             edgecolors=[tuple(GL.DISK_SHADE * cols[i]) for i in order], linewidths=0.8 * k)
        pc.set_clip_path(Polygon([(0, 0), (1, 0), (1, 1), (0, 1)], closed=True, transform=ax.transData))
        ax.add_collection(pc)
    ax.add_patch(Polygon([(0, 0), (1, 0), (1, 1), (0, 1)], closed=True, fill=False, edgecolor="black",
                         linewidth=1.2))
    ax.plot(Y[:, 0], Y[:, 1], "o", color="black", markersize=2.6 * k, markeredgewidth=0, clip_on=False)
    pad = 0.004
    ax.set_xlim(-pad, 1 + pad)
    ax.set_ylim(-pad, 1 + pad)
    ax.set_aspect("equal")
    ax.axis("off")
    fig.savefig(path, dpi=200)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", nargs="+", type=int, default=list(NS))
    ap.add_argument("--noise", nargs="+", type=int, default=list(G.NOISE))
    ap.add_argument("--volumes", nargs="+", default=list(G.VOLUMES))
    ap.add_argument("--seeds", nargs="+", type=int, default=[1, 2, 3, 4, 5])
    ap.add_argument("--format", default="pdf", choices=["pdf", "png", "svg"])
    ap.add_argument("--results", default=os.path.join(HERE, "results"))
    ap.add_argument("--out", default=os.path.join(HERE, "images"))
    ap.add_argument("--no-disks", dest="disks", action="store_false")
    ap.add_argument("--force", action="store_true", help="redraw existing images")
    args = ap.parse_args()
    # 2D only: the drawing clips the cells to the unit square.
    t0 = time.time()
    for n in args.n:
        for noise in args.noise:
            for vol in args.volumes:
                d = os.path.join(args.out, "noise%d-%s" % (noise, vol))
                os.makedirs(d, exist_ok=True)
                for seed in args.seeds:
                    path = os.path.join(d, "n%d-seed%d.%s" % (n, seed, args.format))
                    if os.path.exists(path) and not args.force:
                        continue
                    c = G.make_case(args.results, n, noise, vol, seed)
                    Y, w, V = c["truth"], c["weights"], c["volumes"]
                    cells = power_cells(Y, w)
                    err = max(abs(GL.area(p) - v) if p else v for p, v in zip(cells, V)) / V.min()
                    assert err < 1e-6, "%s: cell volumes off by %g of the smallest" % (c["name"], err)
                    draw(path, cells, Y, w, GL.colors(n, seed), args.disks)
        print("n = %d done (%.0f s)" % (n, time.time() - t0), flush=True)


if __name__ == "__main__":
    main()
