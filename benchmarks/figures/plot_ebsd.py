#!/usr/bin/env python3
"""Draws the grain-map figure (fig:ebsd): the EBSD map of Bourne et al. (2025, 2026; external/bourne/Laguerre-fitting),
243 grains with measured areas and centroids, fitted by H_alpha (introFigure). No power
diagram has these areas and centroids exactly, so the fit keeps the areas exact and moves
the centroids as little as the stiffness allows. The data are scaled to the unit square
and made compatible with eq:partition-feasibility as in their code (make_data_compatible),
the sites are the targets, q = cbar, and alpha is a fraction of the practical bound
(eq:regularized-practical). Each cell takes the colour of the micrograph
(EBSD_image.png, whose rows are the data's y) at its target centroid, or just above
its scale bar.
One column-wide image: the measured map (cropped to the data box) and our fit side by side,
with the square of the largest centroid misses in the middle of the maps (above the inset)
outlined in both, and a close-up of the fit
on that square, with the targets (dots) and the fitted centroids (circles), as an inset
centred between them that hangs a little below (the labels of the maps go beside it).
  latex/figures/ebsd.pdf
Needs numpy, matplotlib and the introFigure tool (ninja -C build introFigure):
  external/bourne/.venv/bin/python benchmarks/figures/plot_ebsd.py [--alpha -0.99]
"""

import argparse
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import plot_intro as PI  # noqa: E402  (sets MPLCONFIGDIR and the Agg backend)
import plot_maxent as PM  # noqa: E402

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.patches import Polygon  # noqa: E402

GL = PI.GL
DATA = os.path.join(PI.ROOT, "external", "bourne", "Laguerre-fitting")
OUT = os.path.join(PI.ROOT, "latex", "figures", "ebsd.pdf")


def load(path):
    """Areas and centroids of the grains in the unit square (y up), made compatible with the
    conservation identities as in their make_data_compatible, and the centroids in pixels of
    the micrograph of width W."""
    E = np.loadtxt(os.path.join(path, "EBSD_data.txt"))
    L = np.sqrt(E[:, 3].sum())                          # the side of the box of their 2025 paper (microns)
    V = E[:, 3] / E[:, 3].sum()
    c = np.column_stack([E[:, 1] / L, 1.0 - E[:, 2] / L])
    c = c + 0.5 - V @ c                                 # sum_i V_i c_i = the centroid of the square
    return V, c, E[:, 1:3] / L


SCALE_BAR = (930, 990, 12, 384)  # rows and columns of the scale bar box of the micrograph (pixels)


def colours(path, pix):
    """The micrograph's colour at each centroid (rows are the data's y), the median of a
    small window; under the scale bar, just above it."""
    img = plt.imread(os.path.join(path, "EBSD_image.png"))[..., :3]
    H, W = img.shape[:2]
    out = []
    for x, y in pix:
        i, j = int(y * W), int(x * W)
        r0, r1, c0, c1 = SCALE_BAR
        if r0 - 3 <= i <= r1 + 3 and c0 - 3 <= j <= c1 + 3:
            i = r0 - 6
        win = img[max(i - 2, 0):min(i + 3, H), max(j - 2, 0):min(j + 3, W)].reshape(-1, 3)
        out.append(np.median(win, axis=0))
    return np.array(out)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--alpha", type=float, default=-0.99,
                    help="stiffness of H_alpha (< 0: that fraction of the practical bound)")
    ap.add_argument("--zoom", type=float, default=0.25, help="side of the zoomed square")
    ap.add_argument("--zoom-ymin", type=float, default=0.56, help="lowest bottom edge of the zoomed square")
    ap.add_argument("--side", type=float, default=1.62, help="printed side of each map (inches)")
    ap.add_argument("--gap", type=float, default=0.09, help="printed gap between the maps (inches)")
    ap.add_argument("--inset", type=float, default=1.0, help="printed side of the close-up (inches)")
    ap.add_argument("--hang", type=float, default=0.2,
                    help="how far the close-up hangs below the maps (inches); the labels go in that strip")
    ap.add_argument("--tool", default=os.path.join(PI.ROOT, "build", "introFigure"))
    ap.add_argument("--out", default=OUT)
    ap.add_argument("--dpi", type=int, default=300, help="for raster previews")
    args = ap.parse_args()

    V, cbar, pix = load(DATA)
    n = len(V)
    rec, msg = PM.run_tool(cbar, V, args.tool, args.alpha)
    print("n = %d; %s" % (n, msg))
    if rec is None:
        raise SystemExit(1)
    p, lam = rec["generators"], rec["weights"]
    cells = GL.power_cells(p, lam)
    A = np.array([PI.area(P) for P in cells])
    C = np.array([PI.centroid(P) for P in cells])
    r = np.sqrt(V)
    err = np.linalg.norm(C - cbar, axis=1) / r
    print("drawn cells: max |V - Vbar| / Vbar %.1e; |c - cbar| / r_i: median %.3f, 90%% %.3f, max %.3f" % (
        np.max(np.abs(A - V) / V), np.median(err), np.quantile(err, 0.9), err.max()))

    cols = colours(DATA, pix)
    # The zoom: the square of side args.zoom with the largest sum of centroid misses of the
    # cells whose target it contains.
    z = args.zoom
    # The square is searched in the middle of the maps and above the part that the inset covers
    # (centre x in [0.35, 0.65], bottom y >= args.zoom_ymin), so it stays visible in both maps.
    best = max(((sum(err[np.all((cbar >= (x, y)) & (cbar <= (x + z, y + z)), axis=1)]), x, y)
                for x in np.linspace(0.35 - 0.5 * z, 0.65 - 0.5 * z, 31)
                for y in np.linspace(args.zoom_ymin, 1 - z, 31)))
    box = (best[1], best[2], best[1] + z, best[2] + z)
    print("zoom box (%.3f, %.3f) to (%.3f, %.3f)" % box)
    base, ext = os.path.splitext(args.out)

    def frame(ax, lim):
        x0, y0, x1, y1 = lim
        pad = 0.006 * (x1 - x0)
        ax.set_xlim(x0 - pad, x1 + pad)
        ax.set_ylim(y0 - pad, y1 + pad)
        ax.set_aspect("equal")
        ax.axis("off")

    def outline(ax, lim, lw, color="black"):
        x0, y0, x1, y1 = lim
        ax.add_patch(Polygon([(x0, y0), (x1, y0), (x1, y1), (x0, y1)], closed=True, fill=False, edgecolor=color,
                             linewidth=lw, zorder=6))

    def diagram(ax, lw):
        for P, c in zip(cells, cols):
            ax.add_patch(Polygon(P, closed=True, facecolor=c, edgecolor="black", linewidth=lw, joinstyle="round"))

    # One column-wide image: the measured map and our fit side by side, and the close-up as an
    # inset centred on the gap that hangs a little (--hang) below them; the paper sets the
    # labels of the maps in that strip, beside the close-up. Sizes are the printed ones (inches),
    # so the line widths are in printed points.
    img = plt.imread(os.path.join(DATA, "EBSD_image.png"))[..., :3]
    W = img.shape[1]
    side, gap, inset, mat, hang = args.side, args.gap, args.inset, 0.03, args.hang
    width, height = 2 * side + gap, side + hang
    fig = plt.figure(figsize=(width, height))

    def axes(x, y, w, h):                               # in inches from the bottom left
        return fig.add_axes([x / width, y / height, w / width, h / height])

    ax = axes(0, hang, side, side)                      # the measured map, cropped to the data box
    ax.imshow(img[:W], extent=(0, 1, 0, 1), origin="upper", interpolation="nearest")
    outline(ax, (0, 0, 1, 1), 0.8)
    outline(ax, box, 0.8)
    frame(ax, (0, 0, 1, 1))

    ax = axes(side + gap, hang, side, side)             # our fit
    diagram(ax, 0.2)
    outline(ax, (0, 0, 1, 1), 0.8)
    outline(ax, box, 0.8)
    frame(ax, (0, 0, 1, 1))

    x0 = 0.5 * width - 0.5 * inset                      # the close-up, on a white mat
    mat_ax = axes(x0 - mat, 0, inset + 2 * mat, inset + 2 * mat)
    mat_ax.set_facecolor("white")
    mat_ax.set_xticks([])
    mat_ax.set_yticks([])
    for sp in mat_ax.spines.values():
        sp.set_visible(False)
    ax = axes(x0, mat, inset, inset)
    diagram(ax, 0.35)
    for i in range(n):                                  # targets (dots) and fitted centroids (circles)
        ax.plot(*zip(cbar[i], C[i]), color="black", linewidth=0.6, solid_capstyle="round", zorder=3)
    ax.plot(cbar[:, 0], cbar[:, 1], "o", color="black", markersize=1.6, markeredgewidth=0, zorder=4)
    ax.plot(C[:, 0], C[:, 1], "o", markerfacecolor="white", markeredgecolor="black", markersize=1.6,
            markeredgewidth=0.4, zorder=5)
    outline(ax, box, 0.8)
    frame(ax, box)
    for a_ in ax.patches + ax.lines:
        a_.set_clip_box(ax.bbox)
    fig.savefig(args.out, dpi=args.dpi)
    plt.close(fig)
    print("wrote " + args.out)


if __name__ == "__main__":
    main()
