#!/usr/bin/env python3
"""Draws the second panel of the showcase figure of the introduction (latex/intro.tex,
fig:intro): the SIGGRAPH logo hidden in the generators of a power diagram with uniform
volumes. Generators Y are placed by dart throwing inside the logo, filled from its
outlines (data/siggraph-logo.png), and the power diagram of Y with uniform volumes (F at
Y) gives the targets (Vbar, cbar), whose cells and centroids cover the whole square
evenly. H is then solved from these volumes and centroids alone, from sites at the
targets, q = cbar (introFigure), and its generators trace the logo. They are unique up
to a homothety, which leaves the cells unchanged; the panel shows the one whose bounding
box is centred and 0.7 wide. The generators are coloured by the part of the logo they
come from.
  latex/figures/intro-logo.pdf
Needs numpy, matplotlib, scipy and the introFigure tool (ninja -C build introFigure):
  external/bourne/.venv/bin/python benchmarks/figures/plot_logo.py [--seed 1]
"""

import argparse
import os
import sys
from collections import deque

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import plot_intro as PI  # noqa: E402  (sets MPLCONFIGDIR and the Agg backend)

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.patches import Polygon  # noqa: E402

GL = PI.GL
LOGO = os.path.join(HERE, "data", "siggraph-logo.png")
OUT = os.path.join(PI.ROOT, "latex", "figures", "intro-logo.pdf")
BLUE = np.array([0.0, 0.0, 1.0])     # the outline colours of the logo file
RED = np.array([0.9, 0.0, 0.0])


def logo_mask(path):
    """The parts of the logo, filled from their closed outlines: an image of labels, 0 outside,
    1 for the blue part and 2 for the red one, with the row index flipped so that y points up."""
    rgb = plt.imread(path)[..., :3]
    blue = (rgb[..., 2] > 0.6) & (rgb[..., 0] < 0.5)
    red = (rgb[..., 0] > 0.6) & (rgb[..., 2] < 0.5)
    line = blue | red
    H, W = line.shape
    outside = np.zeros_like(line)
    outside[0, 0] = True
    queue = deque([(0, 0)])
    while queue:                                        # flood fill of the outside from a corner
        y, x = queue.popleft()
        for a, b in ((y + 1, x), (y - 1, x), (y, x + 1), (y, x - 1)):
            if 0 <= a < H and 0 <= b < W and not outside[a, b] and not line[a, b]:
                outside[a, b] = True
                queue.append((a, b))
    # Each inside pixel takes the colour of the nearest outline pixel along its row and column.
    label = np.zeros(line.shape, dtype=int)
    inside = ~outside
    for y in range(H):
        xs = np.nonzero(inside[y])[0]
        if len(xs) == 0:
            continue
        lx = np.nonzero(line[y])[0]
        for x in xs:
            k = lx[np.argmin(np.abs(lx - x))]
            label[y, x] = 1 if blue[y, k] else 2
    return label[::-1]


class Logo:
    """The logo mask mapped to the box [lo, hi]^2 of the unit square, with the distance to
    the logo."""

    def __init__(self, path, lo, hi):
        self.label = logo_mask(path)
        rows, cols = np.nonzero(self.label)
        # the bounding square of the logo pixels, centred in [lo, hi]^2
        c = np.array([cols.mean(), rows.mean()])
        half = 0.5 * max(cols.max() - cols.min(), rows.max() - rows.min())
        self.origin, self.scale = c - half, (hi - lo) / (2 * half)
        self.lo = lo
        # distance of each pixel to the logo, in units of the unit square
        from scipy.ndimage import distance_transform_edt
        self.dist = distance_transform_edt(self.label == 0) * self.scale

    def _pixel(self, x):
        pix = self.origin + (np.asarray(x) - self.lo) / self.scale
        j, i = int(round(pix[0])), int(round(pix[1]))
        H, W = self.label.shape
        return (i, j) if 0 <= i < H and 0 <= j < W else None

    def distance(self, x):
        """The distance of x to the logo (far outside the image: a large value)."""
        ij = self._pixel(x)
        return self.dist[ij] if ij is not None else 1.0

    def part(self, x):
        """0 outside the logo, 1 or 2 in its blue or red part."""
        ij = self._pixel(x)
        return int(self.label[ij]) if ij is not None else 0


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


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--r", type=float, default=0.02, help="dart radius of the generators in the logo")
    ap.add_argument("--width", type=float, default=0.7, help="width of the generators' bounding box")
    ap.add_argument("--tool", default=os.path.join(PI.ROOT, "build", "introFigure"))
    ap.add_argument("--out", default=OUT)
    ap.add_argument("--dpi", type=int, default=300, help="for raster previews")
    args = ap.parse_args()

    logo = Logo(LOGO, 0.15, 0.85)
    rng = np.random.default_rng(args.seed)
    Y, _ = dart_throwing(rng, lambda y: args.r if logo.part(y) else 10.0, args.r)
    Y = np.array([y for y in Y if logo.part(y)])        # a dart outside the logo blocks nothing (radius 10 > 1)
    n = len(Y)
    rec, msg = PI.solve(Y, None, args.tool, True)
    if rec is None:
        raise SystemExit("introFigure failed: " + msg)
    print("n = %d; %s" % (n, msg))

    # The homothety p -> s p + t of the solution with the generators' bounding box centred and
    # args.width wide; the shifted weights become |s p + t|^2 - s (|p|^2 - lambda) (same cells).
    cbar, p, lam = rec["targets"], rec["generators"], rec["weights"]
    lo, hi = p.min(0), p.max(0)
    sc = args.width / (hi - lo).max()
    p2 = sc * p + (0.5 - sc * 0.5 * (lo + hi))
    lam = (p2 ** 2).sum(1) - sc * ((p ** 2).sum(1) - lam)
    p = p2
    cells = GL.power_cells(p, lam)
    parts = [logo.part(y) for y in Y]
    hues = {1: np.array([0.25, 0.35, 0.85]), 2: np.array([0.85, 0.2, 0.2])}

    fig = plt.figure(figsize=(2.0, 2.0))
    ax = fig.add_axes([0, 0, 1, 1])
    for P in cells:
        ax.add_patch(Polygon(P, closed=True, facecolor="white", edgecolor="0.45", linewidth=0.5 * PI.LW,
                             joinstyle="round"))
    ax.scatter(p[:, 0], p[:, 1], s=4.0, c=[hues[k] for k in parts], linewidths=0, zorder=5)
    ax.add_patch(Polygon([(0, 0), (1, 0), (1, 1), (0, 1)], closed=True, fill=False, edgecolor="black",
                         linewidth=1.2 * PI.LW))
    pad = 0.006
    ax.set_xlim(-pad, 1 + pad)
    ax.set_ylim(-pad, 1 + pad)
    ax.set_aspect("equal")
    ax.axis("off")
    fig.savefig(args.out, dpi=args.dpi)
    plt.close(fig)

    # Cross-check the drawn cells (after the homothety) against the targets.
    V = np.array([PI.area(P) for P in cells])
    C = np.array([PI.centroid(P) for P in cells])
    print("drawn cells: max |V - Vbar| / Vbar %.1e, max |c - cbar| / h %.1e" % (
        np.max(np.abs(V - rec["volumes"]) / rec["volumes"]), np.max(np.linalg.norm(C - cbar, axis=1)) * np.sqrt(n)))
    print("wrote " + args.out)


if __name__ == "__main__":
    main()
