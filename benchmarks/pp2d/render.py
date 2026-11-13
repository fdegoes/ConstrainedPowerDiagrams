#!/usr/bin/env python3
"""Renders the power diagrams of a runPP2D run (frame_NNNN.csv), with cells
colored by the initial position of their particle, so that the mixing of
the flow is visible (as the weight-optimization figure of Power Particles,
Before2015/pp/siggraph2015/images/wopt).

Colors (--color): band (band2 cases, default) the initial x, the middle
band red and the sides blue; radius (other cases, default) the initial
distance to the domain center, blue to red; each particle keeps its color
in every frame. area colors the cells by their current area, as the
figure of the paper (small blue, large red). The cells are the power cells of the state (generators q + u, power
weights w + |u|^2) clipped to the square domain [-0.5, 0.5]^2.

Usage:
  external/bourne/.venv/bin/python benchmarks/pp2d/render.py <run_dir>
      [--frames 0 100 300 | --every 50] [--color band|radius] [--points]
      [--format png|pdf] [--size 6] [--out <run_dir>/images]
"""

import argparse
import csv
import glob
import os
import sys

os.environ.setdefault("MPLCONFIGDIR", "/tmp/mplcfg")
import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.collections import PolyCollection  # noqa: E402
from scipy.spatial import ConvexHull  # noqa: E402

BOX = (-0.5, -0.5, 0.5, 0.5)


def read_frame(path):
    rows = list(csv.DictReader(open(path)))
    f = lambda k: np.array([float(r[k]) for r in rows])
    q = np.c_[f("x"), f("y")]
    u = np.c_[f("ux"), f("uy")]
    w = f("w")
    return q, u, w


def polygon_area(c):
    return 0.5 * abs(np.dot(c[:, 0], np.roll(c[:, 1], -1)) - np.dot(c[:, 1], np.roll(c[:, 0], -1)))


def power_cells(p, lam):
    """Power cells of generators p with power weights lam, clipped to BOX: the
    lower hull of the lifted points (p, |p|^2 - lam) gives the neighbors, each
    cell is the box clipped by its neighbors' bisectors."""
    L = np.c_[p, (p * p).sum(1) - lam]
    hull = ConvexHull(L)
    nb = [set() for _ in range(len(p))]
    for s, eq in zip(hull.simplices, hull.equations):
        if eq[2] < 0:
            for i in s:
                nb[i].update(int(j) for j in s if j != i)
    x0, y0, x1, y1 = BOX
    cells = []
    for i in range(len(p)):
        poly = np.array([[x0, y0], [x1, y0], [x1, y1], [x0, y1]])
        for j in nb[i]:
            # Half-plane 2 (p_j - p_i) . x <= |p_j|^2 - lam_j - |p_i|^2 + lam_i.
            a = 2.0 * (p[j] - p[i])
            b = p[j] @ p[j] - lam[j] - p[i] @ p[i] + lam[i]
            poly = clip(poly, a, b)
            if len(poly) == 0:
                break
        cells.append(poly)
    return cells


def clip(poly, a, b, eps=1e-14):
    out = []
    n = len(poly)
    for k in range(n):
        P, Q = poly[k], poly[(k + 1) % n]
        fp, fq = a @ P - b, a @ Q - b
        if fp <= eps:
            out.append(P)
        if (fp < -eps and fq > eps) or (fp > eps and fq < -eps):
            out.append(P + fp / (fp - fq) * (Q - P))
    return np.array(out)


def colors_for(q0, mode):
    if mode == "band":
        # The middle band [-0.15, 0.2] of band2 cases red, the rest blue.
        t = ((q0[:, 0] > -0.17) & (q0[:, 0] < 0.2)).astype(float)
    else:
        r = np.hypot(q0[:, 0], q0[:, 1])
        t = 1.0 - np.clip(r / 0.5, 0.0, 1.0)
    return plt.get_cmap("jet")(0.1 + 0.8 * t)


def draw(path, cells, q, colors, points, size):
    fig = plt.figure(figsize=(size, size))
    ax = fig.add_axes([0, 0, 1, 1])
    polys = [c for c in cells if len(c) >= 3]
    cols = [colors[i] for i, c in enumerate(cells) if len(c) >= 3]
    lw = 0.25 * min(1.0, np.sqrt(2000.0 / max(len(q), 1)))
    ax.add_collection(PolyCollection(polys, facecolors=cols, edgecolors="black", linewidths=lw))
    if points:
        ax.plot(q[:, 0], q[:, 1], ".", color="black", markersize=0.8)
    ax.set_xlim(BOX[0], BOX[2])
    ax.set_ylim(BOX[1], BOX[3])
    ax.set_aspect("equal")
    ax.axis("off")
    fig.savefig(path, dpi=200)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run")
    ap.add_argument("--frames", nargs="+", type=int, default=None)
    ap.add_argument("--every", type=int, default=None)
    ap.add_argument("--color", choices=["band", "radius", "area"], default=None,
                    help="default: band for band2 cases, radius otherwise")
    ap.add_argument("--points", action="store_true", help="also draw the particles")
    ap.add_argument("--format", choices=["png", "pdf"], default="png")
    ap.add_argument("--size", type=float, default=6.0)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    files = sorted(glob.glob(os.path.join(args.run, "frame_*.csv")))
    if not files:
        sys.exit("no frame_*.csv in " + args.run)
    index = {int(os.path.basename(f)[6:10]): f for f in files}
    if args.frames is not None:
        frames = [k for k in args.frames if k in index]
    elif args.every is not None:
        frames = [k for k in sorted(index) if k % args.every == 0]
    else:
        last = max(index)
        frames = sorted({0, last // 3, 2 * last // 3, last})
    mode = args.color or ("band" if "band" in os.path.abspath(args.run) else "radius")
    out = args.out or os.path.join(args.run, "images")
    os.makedirs(out, exist_ok=True)

    q0, _, _ = read_frame(index[0])
    colors = colors_for(q0, mode) if mode != "area" else None
    lo = hi = None
    for k in frames:
        q, u, w = read_frame(index[k])
        p = q + u
        cells = power_cells(p, w + (u * u).sum(1))
        if mode == "area":
            a = np.array([polygon_area(c) if len(c) >= 3 else 0.0 for c in cells])
            if lo is None:  # the range of the first frame drawn, for all frames
                lo, hi = a.min(), a.max()
            colors = plt.get_cmap("jet")(np.clip((a - lo) / max(hi - lo, 1e-300), 0.0, 1.0))
        path = os.path.join(out, "frame_%04d.%s" % (k, args.format))
        draw(path, cells, q, colors, args.points, args.size)
        area = sum(polygon_area(c) for c in cells if len(c) >= 3)
        print("frame %4d -> %s  (%d cells, area %.6f)" % (k, os.path.relpath(path), len(cells), area), flush=True)


if __name__ == "__main__":
    main()
