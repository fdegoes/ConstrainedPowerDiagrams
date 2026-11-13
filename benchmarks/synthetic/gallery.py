#!/usr/bin/env python3
"""Draws the power diagrams of the synthetic benchmark at n = 10 (seed 1),
one PDF per noise level and volume distribution, for the 2 x 3 gallery of the
paper:
  <out>/synthetic-noise<k>-<volumes>.pdf
The cells are clipped to the unit square and filled with a light tint of a
random colour; the colour belongs to the generator, so the two volume distributions of one
noise level (which share Y) show the same cell in the same colour. The
generators Y are drawn as dots inside a disk that illustrates its weight
(schematic, not the power circle): the weights of a case are shifted to zero
mean (they are defined up to a constant), and the radius is
R0 + (R1 - R0) w_i / max |w|, the max over all six panels, so that the
panels share one scale. The voronoi panels have equal weights and are drawn
without disks. As in the lifting inset
(latex/figures/lifting.tex), a disk is filled with a medium tint of its cell's
colour and outlined in a dark shade of it.

Usage:
  external/bourne/.venv/bin/python benchmarks/synthetic/gallery.py
      [--n 10] [--seed 1] [--out latex/figures/synthetic] [--no-generators] [--no-disks]
"""

import argparse
import colorsys
import os
import sys

os.environ.setdefault("MPLCONFIGDIR", "/tmp/mplcfg")
import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.patches import Circle, Polygon  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import generate as G  # noqa: E402

EPS = 1e-12


def power_cells(Y, w):
    """The cells of the power diagram of (Y, w) clipped to the unit square:
    cell i = {x : |x - Y_i|^2 - w_i <= |x - Y_j|^2 - w_j for all j}, i.e. the
    half-planes 2 (Y_j - Y_i) . x <= |Y_j|^2 - w_j - |Y_i|^2 + w_i."""
    cells = []
    for i in range(len(Y)):
        poly = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]
        for j in range(len(Y)):
            if j == i or not poly:
                continue
            a = 2.0 * (Y[j] - Y[i])
            b = Y[j] @ Y[j] - w[j] - Y[i] @ Y[i] + w[i]
            poly = clip(poly, a, b)
        cells.append(poly)
    return cells


def clip(poly, a, b):
    """The part of the polygon with a . x <= b (Sutherland-Hodgman)."""
    out = []
    for k in range(len(poly)):
        P, Q = np.asarray(poly[k]), np.asarray(poly[(k + 1) % len(poly)])
        fp, fq = a @ P - b, a @ Q - b
        if fp <= EPS:
            out.append(tuple(P))
        if (fp < -EPS and fq > EPS) or (fp > EPS and fq < -EPS):
            out.append(tuple(P + fp / (fp - fq) * (Q - P)))
    return out


def area(poly):
    x, y = np.array(poly).T
    return 0.5 * abs(x @ np.roll(y, -1) - y @ np.roll(x, -1))


CELL_TINT = 0.30  # cell fill: this much of the colour, the rest white
DISK_TINT = 0.55  # disk fill
DISK_ALPHA = 0.45  # of the disk fill, so that overlapping disks stay readable
DISK_SHADE = 0.60  # disk outline: the colour darkened to this fraction
DISK_R0 = 0.04  # disk radius of the mean weight
DISK_R1 = 0.03  # radius change of the largest |weight| over all panels


def colors(n, seed):
    """n random saturated colours: evenly spaced hues in a random order (so
    that no two cells share a hue), with random lightness and saturation."""
    rng = np.random.default_rng(seed)
    hues = (rng.permutation(n) + rng.random()) / n
    return [np.array(colorsys.hls_to_rgb(h, 0.45 + 0.1 * rng.random(), 0.6 + 0.3 * rng.random())) for h in hues]


def tint(c, t):
    return tuple(t * c + (1.0 - t))


def draw(path, cells, Y, radii, cols, generators, disks):
    fig = plt.figure(figsize=(2.0, 2.0))
    ax = fig.add_axes([0, 0, 1, 1])
    for poly, c in zip(cells, cols):
        if poly:
            ax.add_patch(Polygon(poly, closed=True, facecolor=tint(c, CELL_TINT), edgecolor="black",
                                 linewidth=0.6, joinstyle="round"))
    if disks:
        # Largest first, so that the small disks are on top.
        square = Polygon([(0, 0), (1, 0), (1, 1), (0, 1)], closed=True, transform=ax.transData)
        for i in np.argsort(-radii):
            d = Circle(tuple(Y[i]), radii[i], facecolor=(*tint(cols[i], DISK_TINT), DISK_ALPHA),
                       edgecolor=tuple(DISK_SHADE * cols[i]), linewidth=0.8)
            d.set_clip_path(square)
            ax.add_patch(d)
    ax.add_patch(Polygon([(0, 0), (1, 0), (1, 1), (0, 1)], closed=True, fill=False, edgecolor="black",
                         linewidth=1.2))
    if generators:
        ax.plot(Y[:, 0], Y[:, 1], "o", color="black", markersize=2.2, markeredgewidth=0,
                clip_on=False)
    pad = 0.006
    ax.set_xlim(-pad, 1 + pad)
    ax.set_ylim(-pad, 1 + pad)
    ax.set_aspect("equal")
    ax.axis("off")
    fig.savefig(path)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=10)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--results", default=os.path.join(HERE, "results"))
    ap.add_argument("--out", default=os.path.join(G.ROOT, "latex", "figures", "synthetic"))
    ap.add_argument("--no-generators", dest="generators", action="store_false")
    ap.add_argument("--no-disks", dest="disks", action="store_false", help="omit the weight disks")
    args = ap.parse_args()
    # 2D only: the drawing clips the cells to the unit square.
    os.makedirs(args.out, exist_ok=True)
    cols = colors(args.n, args.seed)
    cases = {(noise, vol): G.make_case(args.results, args.n, noise, vol, args.seed)
             for noise in G.NOISE for vol in G.VOLUMES}
    wmax = max(np.abs(c["weights"] - c["weights"].mean()).max() for c in cases.values())
    for (noise, vol), c in cases.items():
        Y, w, V = c["truth"], c["weights"], c["volumes"]
        radii = DISK_R0 + DISK_R1 * (w - w.mean()) / wmax if wmax > 0 else np.full(len(w), DISK_R0)
        cells = power_cells(Y, w)
        err = max(abs(area(p) - v) if p else v for p, v in zip(cells, V))
        assert err < 1e-9, "%s: cell volumes off by %g" % (c["name"], err)
        path = os.path.join(args.out, "synthetic-noise%d-%s.pdf" % (noise, vol))
        draw(path, cells, Y, radii, cols, args.generators, args.disks and vol != "voronoi")
        print("%s  V max/min %.2f  |w|/h^2 %.2f  ->  %s" % (
            c["name"], V.max() / V.min(), np.abs(w - w.mean()).max() * args.n, os.path.relpath(path)))


if __name__ == "__main__":
    main()
