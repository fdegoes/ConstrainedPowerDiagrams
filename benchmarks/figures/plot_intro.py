#!/usr/bin/env python3
"""Draws the showcase figure of the introduction (latex/intro.tex, fig:intro): power
diagrams of the unit square with n cells of uniform volume and prescribed centroids
cbar (squares) in the rectangle [0, width] x [0, 1] (width 2 by default, so that the figure
fits the column at half its height), solved by H (introFigure) from sites q (circles) that are not the
centroids: the generators p = q + u (dots), with the shifts u as arrows from the sites.
The cells, colours and weight disks are those of the synthetic gallery
(benchmarks/synthetic/gallery.py, fig:synthetic-gallery), at the same printed line
widths: a schematic disk about each generator, not its power circle, of radius
(0.04 + 0.03 (w_i - mean w) / w_max) / sqrt(1/10) h for the cell size h = sqrt(1/n),
with w_max the largest |w_j - mean w|.

The generators Y of the targets are random points of the square, and the sites are
q = Y + a random offset. The solution of H is unique up to a homothety of the
generators, which leaves the cells unchanged; the figure shows the generators p
closest to the sites (--gauge closest), which differ from the sites and from the
centroids.
  latex/figures/intro.pdf
Needs numpy, matplotlib and the introFigure tool (ninja -C build introFigure):
  external/bourne/.venv/bin/python benchmarks/figures/plot_intro.py [--seed 2]
"""

import argparse
import os
import subprocess
import sys
import tempfile

os.environ.setdefault("MPLCONFIGDIR", "/tmp/mplcfg")
import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.patches import Circle, FancyArrowPatch, Polygon  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(ROOT, "benchmarks", "synthetic"))
import gallery as GL  # noqa: E402

OUT = os.path.join(ROOT, "latex", "figures", "intro.pdf")
# The gallery's panels are 2 in wide, printed at 0.29 \columnwidth; these are printed at
# 0.49 \columnwidth, so their lines are thinner by this factor to print alike.
LW = 0.29 / 0.49
H_GALLERY = np.sqrt(1.0 / 10)  # the cell size of the gallery's n = 10


def area(P):
    x, y = np.asarray(P).T
    return 0.5 * (x @ np.roll(y, -1) - y @ np.roll(x, -1))


def centroid(P):
    x, y = np.asarray(P).T
    c = x * np.roll(y, -1) - np.roll(x, -1) * y
    return np.array([((x + np.roll(x, -1)) * c).sum(), ((y + np.roll(y, -1)) * c).sum()]) / (3.0 * c.sum())


def random_square(rng, n, margin, sep, width=1.0):
    """n random points of [margin, width - margin] x [margin, 1 - margin], at least sep apart."""
    Y = []
    while len(Y) < n:
        y = margin + np.array([width - 2 * margin, 1 - 2 * margin]) * rng.random(2)
        if all(np.linalg.norm(y - z) >= sep for z in Y):
            Y.append(y)
    return np.array(Y)


def rectangle(width):
    return np.array([(0.0, 0.0), (width, 0.0), (width, 1.0), (0.0, 1.0)])


def power_cells(p, lam, width=1.0):
    """The power cells of (p, lam) clipped to the rectangle [0, width] x [0, 1]."""
    cells = []
    for i in range(len(p)):
        poly = [tuple(x) for x in rectangle(width)]
        for j in range(len(p)):
            if j != i and poly:
                poly = GL.clip(poly, 2.0 * (p[j] - p[i]), p[j] @ p[j] - p[i] @ p[i] - lam[j] + lam[i])
        cells.append(np.array(poly))
    return cells


def solve(Y, q, tool, closest, volumes=None, width=1.0):
    """Runs introFigure from the sites q, or from sites at the targets if q is None, with
    uniform volumes or the given ones."""
    n = len(Y)
    with tempfile.TemporaryDirectory() as tmp:
        src, dst = os.path.join(tmp, "in.txt"), os.path.join(tmp, "out.txt")
        with open(src, "w") as f:
            f.write("n %d\n" % n)
            if width != 1.0:
                f.write("domain " + " ".join("%.17g" % x for x in rectangle(width).ravel()) + "\n")
            f.write("generators " + " ".join("%.17g" % x for x in Y.ravel()) + "\n")
            if q is None:
                f.write("sitesAtTargets 1\n")
            else:
                f.write("sites " + " ".join("%.17g" % x for x in q.ravel()) + "\n")
            if volumes is not None:
                f.write("volumes " + " ".join("%.17g" % x for x in volumes) + "\n")
            f.write("closest %d\n" % closest)
        run = subprocess.run([tool, src, dst], capture_output=True, text=True)
        if run.returncode != 0:
            return None, run.stderr.strip()
        rec = {}
        with open(dst) as f:
            for line in f:
                key, *vals = line.split()
                rec[key] = np.array([float(v) for v in vals])
    for k in ("targets", "sites", "generators", "shifts"):
        rec[k] = rec[k].reshape(n, 2)
    return rec, run.stdout.strip()


def panel(path, cells, centers, radii, cols, dpi, draw_more, width=1.0):
    """The gallery's panel (GL.draw) at the printed widths of LW: cells, weight disks about
    centers, the rectangle [0, width] x [0, 1], then draw_more(ax) for the points."""
    fig = plt.figure(figsize=(2.0 * width, 2.0))
    ax = fig.add_axes([0, 0, 1, 1])
    for poly, c in zip(cells, cols):
        ax.add_patch(Polygon(poly, closed=True, facecolor=GL.tint(c, GL.CELL_TINT), edgecolor="black",
                             linewidth=0.6 * LW, joinstyle="round"))
    square = Polygon(rectangle(width), closed=True, transform=ax.transData)
    for i in np.argsort(-radii):                        # largest first, so that the small disks are on top
        d = Circle(tuple(centers[i]), radii[i], facecolor=(*GL.tint(cols[i], GL.DISK_TINT), GL.DISK_ALPHA),
                   edgecolor=tuple(GL.DISK_SHADE * cols[i]), linewidth=0.8 * LW)
        d.set_clip_path(square)
        ax.add_patch(d)
    ax.add_patch(Polygon(rectangle(width), closed=True, fill=False, edgecolor="black", linewidth=1.2 * LW))
    draw_more(ax)
    pad = 0.006
    ax.set_xlim(-pad, width + pad)
    ax.set_ylim(-pad, 1 + pad)
    ax.set_aspect("equal")
    ax.axis("off")
    fig.savefig(path, dpi=dpi)
    plt.close(fig)


def sites(ax, q):
    ax.plot(q[:, 0], q[:, 1], "o", markerfacecolor="white", markeredgecolor="black", markersize=3.0,
            markeredgewidth=0.5, zorder=7, clip_on=False)


def targets(ax, cbar):
    ax.plot(cbar[:, 0], cbar[:, 1], "s", markerfacecolor="white", markeredgecolor="black", markersize=3.0,
            markeredgewidth=0.5, zorder=8)


def draw(rec, cols, out, dpi, width=1.0):
    n = len(rec["volumes"])
    cbar, q, p, lam = rec["targets"], rec["sites"], rec["generators"], rec["weights"]
    h = np.sqrt(width / n)
    # The gallery's disks, scaled with the cell size.
    wmax = np.abs(lam - lam.mean()).max()
    radii = h / H_GALLERY * (GL.DISK_R0 + GL.DISK_R1 * (lam - lam.mean()) / wmax)
    cells = power_cells(p, lam, width)

    def points(ax):
        for i in range(n):
            ax.add_patch(FancyArrowPatch(q[i], p[i], arrowstyle="-|>", mutation_scale=6, shrinkA=1.5, shrinkB=1.5,
                                         color="black", linewidth=0.5, zorder=5))
        sites(ax, q)
        targets(ax, cbar)
        ax.plot(p[:, 0], p[:, 1], "o", color="black", markersize=2.4, markeredgewidth=0, zorder=9, clip_on=False)

    panel(out, cells, p, radii, cols, dpi, points, width)
    return cells, radii


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--n", type=int, default=8)
    ap.add_argument("--width", type=float, default=2.0, help="the domain is [0, width] x [0, 1]")
    ap.add_argument("--seed", type=int, default=2)
    ap.add_argument("--offset", type=float, default=0.35, help="|q - Y| relative to the cell size")
    ap.add_argument("--gauge", choices=("closest", "solver"), default="closest",
                    help="the generators closest to the sites (up to a homothety), or those the solver reaches")
    ap.add_argument("--tool", default=os.path.join(ROOT, "build", "introFigure"))
    ap.add_argument("--out", default=OUT)
    ap.add_argument("--dpi", type=int, default=300, help="for raster previews")
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)
    w = args.width
    h = np.sqrt(w / args.n)
    Y = random_square(rng, args.n, 0.1, 0.6 * h, w)
    q = np.empty_like(Y)
    for i in range(args.n):                             # sites inside the rectangle, off its boundary
        while True:
            b = 2 * np.pi * rng.random()
            q[i] = Y[i] + args.offset * h * np.array([np.cos(b), np.sin(b)])
            if 0.04 <= q[i][0] <= w - 0.04 and 0.04 <= q[i][1] <= 0.96:
                break
    rec, msg = solve(Y, q, args.tool, args.gauge == "closest", width=w)
    if rec is None:
        raise SystemExit("introFigure failed: " + msg)
    print(msg)

    cols = GL.colors(args.n, args.seed)
    cells, radii = draw(rec, cols, args.out, args.dpi, w)
    # Cross-check the drawn cells against the solver's targets.
    V = np.array([area(P) for P in cells])
    C = np.array([centroid(P) for P in cells])
    print("drawn cells: max |V - Vbar| / Vbar %.1e, max |c - cbar| / h %.1e" % (
        np.max(np.abs(V - rec["volumes"]) / rec["volumes"]), np.max(np.linalg.norm(C - rec["targets"], axis=1)) / h))
    print("|q - cbar| / h: min %.2f;  |p - cbar| / h: min %.2f max %.2f;  |u| / h: max %.2f;  radii / h: %.2f to %.2f" % (
        np.min(np.linalg.norm(rec["sites"] - rec["targets"], axis=1)) / h,
        np.min(np.linalg.norm(rec["generators"] - rec["targets"], axis=1)) / h,
        np.max(np.linalg.norm(rec["generators"] - rec["targets"], axis=1)) / h,
        np.max(np.linalg.norm(rec["shifts"], axis=1)) / h, radii.min() / h, radii.max() / h))
    print("wrote " + args.out)


if __name__ == "__main__":
    main()
