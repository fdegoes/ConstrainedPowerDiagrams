#!/usr/bin/env python3
"""Draws the reducible diagram of the null-space appendix (appendix-nullspace.tex):
a centroidal Voronoi diagram of the unit square whose generators are symmetric
under the eight reflections of the square, and its coarsening into eight
sectors. The symmetry makes the diagonals and the medians unions of facets,
with vertices shared by four cells along them and one shared by eight cells at
the centre, as in the aligned example of Bourne et al. 2026, but with every
generator at the centroid of its cell. The eight sectors are the Voronoi
diagram of their centroids, whose cells are unions of the dense cells.
  nullspace/reducible-dense.pdf          the dense diagram, tinted by sector
  nullspace/reducible-dense-closeup.pdf  (--closeups) the disk around a vertex shared by four
                                         cells, with the dual graph and its 4-cycle
  nullspace/reducible-coarse.pdf         the eight sectors
  nullspace/reducible-coarse-closeup.pdf (--closeups) the same disk on the eight sectors
  nullspace/reducible-cut.pdf            a diagram in general position with a cut cell (orange)
  nullspace/reducible-cut-coarse.pdf     its coarsening, which merges the cut cell with the cells on its right
Needs numpy and matplotlib:
  external/bourne/.venv/bin/python latex/figures/nullspace_reducible.py [--per-sector 12] [--seed 3]
"""

import argparse
import math
import os

os.environ.setdefault("MPLCONFIGDIR", "/tmp/mplcfg")
import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.patches import Circle, Polygon, Rectangle  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "nullspace")
EPS = 1e-12
VTOL = 1e-7                 # vertices closer than this are one vertex
MARGIN = 0.02               # generators keep this distance from the mirror lines
BLUE = np.array([66, 113, 184]) / 255      # nsBlue, nsRed, nsGreen, nsOrange, nsGray of nullspace-styles.tex
RED = np.array([205, 50, 50]) / 255
GREEN = np.array([60, 160, 90]) / 255
ORANGE = np.array([230, 140, 40]) / 255
GRAY = np.array([150, 150, 150]) / 255
SECTOR_TINTS = [tuple(0.16 * c + 0.84) for c in (BLUE, GREEN)]   # alternating around the centre
FOCUS, FOCUS_R = (0.5, 0.8), 0.12           # the close-up: a disk around the vertex nearest FOCUS on the median


def images(p):
    """The eight images of p under the reflections of the square about its centre."""
    u, v = p[0] - 0.5, p[1] - 0.5
    return [(0.5 + a, 0.5 + b) for a, b in ((u, v), (v, u), (-v, u), (-u, v), (-u, -v), (-v, -u), (v, -u), (u, -v))]


def fold(p):
    """The image of p in the fundamental triangle 0 <= v <= u (sector 0), kept off the mirror lines."""
    u, v = abs(p[0] - 0.5), abs(p[1] - 0.5)
    u, v = max(u, v), min(u, v)
    v = min(max(v, MARGIN), u - MARGIN)
    return (0.5 + u, 0.5 + v)


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


def cell(Y, i, w=None, box=(0.0, 0.0, 1.0, 1.0)):
    """The power cell of site i (weights w, zero by default) clipped to the box."""
    x0, y0, x1, y1 = box
    w = np.zeros(len(Y)) if w is None else w
    poly = [(x0, y0), (x1, y0), (x1, y1), (x0, y1)]
    for j in np.argsort(np.linalg.norm(Y - Y[i], axis=1))[1:]:
        if not poly:
            break
        poly = clip(poly, 2.0 * (Y[j] - Y[i]), Y[j] @ Y[j] - Y[i] @ Y[i] - w[j] + w[i])
    return np.array(poly)


def voronoi(Y, w=None, box=(0.0, 0.0, 1.0, 1.0)):
    return [cell(Y, i, w, box) for i in range(len(Y))]


def area(P):
    x, y = np.asarray(P).T
    return 0.5 * abs(x @ np.roll(y, -1) - y @ np.roll(x, -1))


def centroid(P):
    x, y = np.asarray(P).T
    c = x * np.roll(y, -1) - np.roll(x, -1) * y
    return np.array([((x + np.roll(x, -1)) * c).sum(), ((y + np.roll(y, -1)) * c).sum()]) / (3.0 * c.sum())


def symmetric_cvt(per_sector, seed, iters=300):
    """Lloyd iterations on the orbits of per_sector points of sector 0."""
    rng = np.random.default_rng(seed)
    F = []
    while len(F) < per_sector:
        u, v = 0.5 * rng.random(2)
        if v < u:
            F.append(fold((0.5 + u, 0.5 + v)))
    for _ in range(iters):
        Y = np.array([q for p in F for q in images(p)])
        F = [fold(centroid(cell(Y, 8 * k))) for k in range(len(F))]
    return np.array([q for p in F for q in images(p)])


def vertices(cells):
    pts, owners = [], []
    for i, P in enumerate(cells):
        for p in P:
            for k, q in enumerate(pts):
                if abs(q[0] - p[0]) + abs(q[1] - p[1]) < VTOL:
                    owners[k].add(i)
                    break
            else:
                pts.append(p)
                owners.append({i})
    return np.array(pts), owners


def dual_edges(owners):
    shared = {}
    for k, o in enumerate(owners):
        for i in o:
            for j in o:
                if i < j:
                    shared.setdefault((i, j), set()).add(k)
    return [e for e, ks in shared.items() if len(ks) >= 2]


def degenerate(pts, owners, Y):
    """(vertex, its cells in cyclic order) for the interior vertices shared by k > 3 cells."""
    out = []
    for p, o in zip(pts, owners):
        if len(o) > 3 and min(p[0], p[1], 1 - p[0], 1 - p[1]) > VTOL:
            idx = sorted(o, key=lambda i: math.atan2(Y[i][1] - p[1], Y[i][0] - p[0]))
            out.append((p, idx))
    return out


def sector_of(P):
    c = centroid(P) - 0.5
    return int(math.floor(math.atan2(c[1], c[0]) / (math.pi / 4))) % 8


def new_axes(size=2.0):
    fig = plt.figure(figsize=(size, size))
    return fig, fig.add_axes([0, 0, 1, 1])


def finish(fig, ax, name, box):
    x0, y0, x1, y1 = box
    pad = 0.008 * (x1 - x0)
    ax.set_xlim(x0 - pad, x1 + pad)
    ax.set_ylim(y0 - pad, y1 + pad)
    ax.set_aspect("equal")
    ax.axis("off")
    fig.savefig(os.path.join(OUT, name))
    plt.close(fig)


def full(name, cells, Y, cycles, focus, lw):
    fig, ax = new_axes()
    for P in cells:
        ax.add_patch(Polygon(P, closed=True, facecolor=SECTOR_TINTS[sector_of(P) % 2], edgecolor="black",
                             linewidth=lw, joinstyle="round"))
    ax.plot(Y[:, 0], Y[:, 1], "o", color="black", markersize=1.8, markeredgewidth=0, zorder=4)
    for p, idx in cycles:
        big = len(idx) > 4
        ax.plot(*p, "o", color=RED if big else BLUE, markersize=4.6 if big else 3.0, markeredgewidth=0, zorder=5)
    if focus is not None:
        ax.add_patch(Circle(focus[0], focus[1], fill=False, edgecolor=ORANGE, linewidth=1.6, zorder=7))
    ax.add_patch(Rectangle((0, 0), 1, 1, fill=False, edgecolor="black", linewidth=1.2, zorder=6))
    finish(fig, ax, name, (0, 0, 1, 1))


def closeup(name, cells, Y, edges, cycle, focus):
    (cx, cy), r = focus
    fig, ax = new_axes(1.0)
    disk = Circle((cx, cy), r, transform=ax.transData)
    artists = [ax.add_patch(Polygon(P, closed=True, facecolor=SECTOR_TINTS[sector_of(P) % 2], edgecolor="black",
                                    linewidth=0.8, joinstyle="round")) for P in cells]
    for i, j in edges:
        artists += ax.plot(*zip(Y[i], Y[j]), color=GRAY, linewidth=1.0, solid_capstyle="round", zorder=2)
    if cycle is not None:
        p, idx = cycle
        artists.append(ax.add_patch(Polygon(Y[idx], closed=True, facecolor=(*BLUE, 0.18), edgecolor=BLUE,
                                            linewidth=1.8, joinstyle="round", zorder=3)))
        artists += ax.plot(*p, "o", color=BLUE, markersize=4.0, markeredgewidth=0, zorder=4)
    artists += ax.plot(Y[:, 0], Y[:, 1], "o", color="black", markersize=3.2, markeredgewidth=0, zorder=5)
    for a in artists:
        a.set_clip_path(disk)
    ax.add_patch(Circle((cx, cy), r, fill=False, edgecolor=ORANGE, linewidth=2.0, zorder=6))
    finish(fig, ax, name, (cx - r, cy - r, cx + r, cy + r))


CUT_TINT = tuple(0.4 * ORANGE + 0.6)       # nsCut of nullspace-styles.tex


CUT_BOX = (0.0, 0.0, 2.5, 1.0)             # the rectangle of the cut-cell example


def lloyd_box(n, box, rng, iters=300):
    """A centroidal Voronoi diagram of n sites in the box."""
    x0, y0, x1, y1 = box
    Y = np.column_stack([x0 + (x1 - x0) * rng.random(n), y0 + (y1 - y0) * rng.random(n)])
    for _ in range(iters):
        Y = np.array([centroid(P) for P in voronoi(Y, None, box)])
    return Y


def inside(P, y):
    """Whether y lies in the counter-clockwise convex polygon P."""
    Q = np.roll(P, -1, axis=0)
    return bool(np.all((Q[:, 0] - P[:, 0]) * (y[1] - P[:, 1]) - (Q[:, 1] - P[:, 1]) * (y[0] - P[:, 0]) >= -EPS))


def cut_cell_example(n_side, seed, w_cut=0.12, gap=0.35, iters=300):
    """A power diagram in general position with a cut cell, in the rectangle CUT_BOX. Two
    centroidal Voronoi diagrams of n_side sites fill its left and right parts, up to a gap around
    the middle, where a cut site of weight w_cut has a cell that spans the rectangle from its bottom
    to its top edge. Following the arm mode of Appendix E until the generators of the right part
    reach the cut site merges them with the cut cell; the coarser diagram is the power diagram of
    the left sites and the cut site alone."""
    rng = np.random.default_rng(seed)
    x0, y0, x1, y1 = CUT_BOX
    xm = 0.5 * (x0 + x1)
    L = lloyd_box(n_side, (x0, y0, xm - gap, y1), rng, iters)
    R_ = lloyd_box(n_side, (xm + gap, y0, x1, y1), rng, iters)
    c = np.array([[xm, 0.5 * (y0 + y1)]])
    Y = np.vstack([L, c, R_])
    side = np.array([-1] * n_side + [0] + [1] * n_side)
    w = np.where(side == 0, w_cut, 0.0)
    cells = voronoi(Y, w, CUT_BOX)
    assert all(len(P) >= 3 for P in cells), "empty cell"
    pts, owners = vertices(cells)
    edges = dual_edges(owners)
    k = n_side
    assert min(cells[k][:, 1]) < y0 + VTOL and max(cells[k][:, 1]) > y1 - VTOL, "the middle cell is not a cut cell"
    assert not any(side[i] * side[j] < 0 for i, j in edges), "a left and a right site are neighbors"
    inner = [o for p, o in zip(pts, owners) if min(p[0] - x0, p[1] - y0, x1 - p[0], y1 - p[1]) > VTOL]
    assert all(len(o) == 3 for o in inner), "not in general position"
    assert all(inside(P, y) for P, y in zip(cells, Y)), "a generator lies outside its cell"
    Yc, wc = Y[:k + 1], w[:k + 1]
    coarse = voronoi(Yc, wc, CUT_BOX)
    merged = sum(area(cells[i]) for i in range(k, len(Y)))
    assert abs(area(coarse[k]) - merged) < 1e-9, "the coarser cell is not the union"
    assert all(abs(area(coarse[i]) - area(cells[i])) < 1e-9 for i in range(k)), "the left cells changed"
    print("cut cell: n = %d, %d dual edges, %d interior vertices, all shared by 3 cells" % (
        len(Y), len(edges), len(inner)))
    return Y, cells, side, Yc, coarse


def cut_panel(name, cells, Y, tints, lw):
    x0, y0, x1, y1 = CUT_BOX
    # same width as the square panels, and the aspect of the padded box of finish(), so that
    # both rows are scaled alike in the paper and the frame fills the figure
    pad = 0.008 * (x1 - x0)
    fig = plt.figure(figsize=(2.0, 2.0 * (y1 - y0 + 2 * pad) / (x1 - x0 + 2 * pad)))
    ax = fig.add_axes([0, 0, 1, 1])
    for P, t in zip(cells, tints):
        ax.add_patch(Polygon(P, closed=True, facecolor=t, edgecolor="black", linewidth=lw, joinstyle="round"))
    ax.plot(Y[:, 0], Y[:, 1], "o", color="black", markersize=1.8, markeredgewidth=0, zorder=4)
    ax.add_patch(Rectangle((x0, y0), x1 - x0, y1 - y0, fill=False, edgecolor="black", linewidth=1.2, zorder=6))
    finish(fig, ax, name, CUT_BOX)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--per-sector", type=int, default=12, help="generators per sector")
    ap.add_argument("--seed", type=int, default=3)
    ap.add_argument("--closeups", action="store_true", help="also draw the close-ups (not used by the paper)")
    ap.add_argument("--n-side", type=int, default=9, help="sites on each side of the cut cell")
    ap.add_argument("--cut-seed", type=int, default=1)
    args = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)

    # dense: the symmetric centroidal Voronoi diagram
    Y = symmetric_cvt(args.per_sector, args.seed)
    cells = voronoi(Y)
    assert all(np.linalg.norm(centroid(P) - y) < 1e-3 for P, y in zip(cells, Y)), "Lloyd has not converged"
    sec = [sector_of(P) for P in cells]
    assert all(abs(sum(area(P) for P, s in zip(cells, sec) if s == k) - 0.125) < 1e-9 for k in range(8))
    pts, owners = vertices(cells)
    edges = dual_edges(owners)
    cycles = degenerate(pts, owners, Y)
    print("dense: n = %d, %d dual edges, %d vertices shared by 4 cells, %d by 8" % (
        len(Y), len(edges), sum(len(i) == 4 for _, i in cycles), sum(len(i) == 8 for _, i in cycles)))

    # coarse: the Voronoi diagram of the sector centroids, which is the eight sectors
    C = np.array([centroid(np.array([(0.5, 0.5), images((1.0, 0.5))[k], images((1.0, 1.0))[k]])) for k in range(8)])
    coarse = voronoi(C)
    assert all(abs(area(P) - 0.125) < 1e-12 for P in coarse)
    cpts, cowners = vertices(coarse)
    cedges = dual_edges(cowners)
    ccycles = degenerate(cpts, cowners, C)

    # the close-up: the vertex shared by four cells on the median nearest FOCUS, centred
    # between the vertex and its cycle
    p, idx = min(((p, i) for p, i in cycles if len(i) == 4), key=lambda c: np.hypot(*(c[0] - np.asarray(FOCUS))))
    focus = (0.5 * (p + Y[idx].mean(0)), FOCUS_R)
    full("reducible-dense.pdf", cells, Y, cycles, None, 0.5)
    full("reducible-coarse.pdf", coarse, C, ccycles, None, 0.8)
    if args.closeups:
        full("reducible-dense-circled.pdf", cells, Y, cycles, focus, 0.35)
        closeup("reducible-dense-closeup.pdf", cells, Y, edges, (p, idx), focus)
        closeup("reducible-coarse-closeup.pdf", coarse, C, cedges, None, focus)

    # cut cell: a diagram in general position, reducible through its cut cell
    Y, cells, side, Yc, cut_coarse = cut_cell_example(args.n_side, args.cut_seed)
    tint = {-1: SECTOR_TINTS[0], 0: CUT_TINT, 1: SECTOR_TINTS[1]}
    cut_panel("reducible-cut.pdf", cells, Y, [tint[s] for s in side], 0.5)
    cut_panel("reducible-cut-coarse.pdf", cut_coarse, Yc, [tint[-1]] * (len(Yc) - 1) + [SECTOR_TINTS[1]], 0.8)


if __name__ == "__main__":
    main()
