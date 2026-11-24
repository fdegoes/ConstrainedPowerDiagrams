#!/usr/bin/env python3
"""Draws the 3D example: a power diagram of the unit cube with prescribed
volumes and centroids, solved by H.

The case is a 3D synthetic case (benchmarks/synthetic/results3d/<case>/case.txt,
start from the sites q = cbar). The script runs benchBourne (solver H) on it,
takes the generators p and the shifted weights lambda of the result, and
computes the cells itself: cell i is the intersection of the cube with the
half-spaces |x - p_i|^2 - lambda_i <= |x - p_j|^2 - lambda_j (scipy
HalfspaceIntersection). It checks their volumes and centroids against the
targets and prints the statistics of the solve. The cells are unchanged by a
homothety of the generators about the cube centre (with the weights changed
to match), along which H is invariant; by default the figure shows the
generators closest to the sites (as in plot_intro.py, --gauge closest).

The drawing is an orthographic cutaway: the cells whose target centroid lies
in the octant nearest to the viewer are removed; their targets cbar (squares)
and generators p (dots) are drawn in the hollow, joined by the shift segment.
The other cells are filled with a light tint of a random colour
(benchmarks/synthetic/gallery.py) and shaded by a directional light. They are
drawn back to front in the order of p . e (e the direction to the viewer): for
two neighbours i, j the facet normal is p_j - p_i, so this order is a
visibility order of a power diagram; within a cell only the facets that face
the viewer are drawn.

  external/bourne/.venv/bin/python benchmarks/figures/plot_3d.py
      [--case noise1-uniform-62-1] [--cut 0.4] [--gauge closest|solver] [--out latex/figures/example-3d.pdf]
"""

import argparse
import json
import os
import subprocess
import sys
import tempfile

os.environ.setdefault("MPLCONFIGDIR", "/tmp/mplcfg")
import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.collections import PolyCollection  # noqa: E402
from scipy.optimize import linprog  # noqa: E402
from scipy.spatial import ConvexHull, HalfspaceIntersection  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(ROOT, "benchmarks", "synthetic"))
import gallery as GL  # noqa: E402

RESULTS = os.path.join(ROOT, "benchmarks", "synthetic", "results3d")
BENCH = os.path.join(ROOT, "build", "benchBourne")
OUT = os.path.join(ROOT, "latex", "figures", "example-3d.pdf")

AZIM, ELEV = 32.0, 24.0  # view direction (degrees)
LIGHT = np.array([0.35, 0.55, 1.0])  # direction to the light
AMBIENT = 0.8  # shade = AMBIENT + (1 - AMBIENT) max(0, n . light)


def read_case(path):
    rec = {}
    with open(path) as f:
        for line in f:
            t = line.split()
            if t:
                rec[t[0]] = t[1:]
    n = int(rec["n"][0])
    vol = np.array(rec["volumes"], float)
    cbar = np.array(rec["targets"], float).reshape(n, 3)
    return n, vol, cbar


def solve(case):
    with tempfile.TemporaryDirectory() as d:
        out = os.path.join(d, "out.json")
        subprocess.run([BENCH, case, out], check=True, stdout=subprocess.DEVNULL)
        with open(out) as f:
            return json.load(f)


def power_cell(i, p, lam):
    """Vertices and facets (outward normal, vertex loop) of cell i in [0,1]^3."""
    # Halfspaces A x + b <= 0.
    others = np.arange(len(p)) != i
    A = 2.0 * (p[others] - p[i])
    b = (p[i] @ p[i] - lam[i]) - (np.sum(p[others] ** 2, axis=1) - lam[others])
    A = np.vstack([A, -np.eye(3), np.eye(3)])
    b = np.concatenate([b, np.zeros(3), -np.ones(3)])
    # Chebyshev centre: an interior point.
    norm = np.linalg.norm(A, axis=1)
    lp = linprog(np.r_[0, 0, 0, -1], A_ub=np.c_[A, norm], b_ub=-b, bounds=[(None, None)] * 3 + [(0, None)])
    if not lp.success or lp.x[3] <= 1e-12:
        return None
    hs = HalfspaceIntersection(np.c_[A, b], lp.x[:3])
    hull = ConvexHull(hs.intersections)
    V = hull.points
    # Volume and centroid from the hull's triangles.
    o = V.mean(axis=0)
    vol, mom = 0.0, np.zeros(3)
    for s in hull.simplices:
        t = V[s]
        v = abs(np.linalg.det(t - o)) / 6.0
        vol += v
        mom += v * (t.sum(axis=0) + o) / 4.0
    # Merge coplanar triangles into facets.
    facets = []
    eq = hull.equations
    used = np.zeros(len(eq), bool)
    for k in range(len(eq)):
        if used[k]:
            continue
        same = np.where(np.linalg.norm(eq - eq[k], axis=1) < 1e-9)[0]
        used[same] = True
        idx = np.unique(hull.simplices[same])
        P = V[idx]
        nrm = eq[k, :3]
        c = P.mean(axis=0)
        u = P[0] - c
        u /= np.linalg.norm(u)
        w = np.cross(nrm, u)
        ang = np.arctan2((P - c) @ w, (P - c) @ u)
        facets.append((nrm, P[np.argsort(ang)]))
    return vol, mom / vol, facets


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--case", default="noise1-uniform-62-1")
    ap.add_argument("--cut", type=float, default=0.4, help="remove cells with cbar > cut in all coordinates")
    ap.add_argument("--seed", type=int, default=1, help="seed of the colours")
    ap.add_argument("--size", type=float, default=3.2, help="width of the figure (inches)")
    ap.add_argument("--gauge", choices=("closest", "solver"), default="closest")
    ap.add_argument("--out", default=OUT)
    a = ap.parse_args()

    case = os.path.join(RESULTS, a.case, "case.txt")
    n, vbar, cbar = read_case(case)
    res = solve(case)
    p0 = np.array(res["generators"]).reshape(n, 3)
    lam0 = np.array(res["shiftedWeights"])
    h = (1.0 / n) ** (1.0 / 3.0)
    # The cells are invariant under a homothety of the generators about the
    # centroid c0 of the cube, with the weights changed accordingly; take the
    # generators closest to the sites q = cbar (--gauge closest) or those the
    # solver reaches (--gauge solver).
    c0 = np.full(3, 0.5)
    d0 = p0 - c0
    s = np.sum(d0 * (cbar - c0)) / np.sum(d0 * d0) if a.gauge == "closest" else 1.0
    p = c0 + s * d0
    lam = s * lam0 + (s * s - s) * np.sum(d0 * d0, axis=1)

    cells = [power_cell(i, p, lam) for i in range(n)]
    if s != 1.0:
        V0 = np.array([power_cell(i, p0, lam0)[0] for i in range(n)])
        print(f"gauge: homothety {s:.6f} about c0; max cell-volume change {np.max(np.abs(V0 - [c[0] for c in cells])):.1e}")
    vol = np.array([c[0] for c in cells])
    cen = np.array([c[1] for c in cells])
    with open(os.path.join(RESULTS, a.case, "truth.json")) as f:
        truth = json.load(f)
    Y = np.array(truth["generators"]).reshape(n, 3)
    print(f"case {a.case}: n = {n}, h = {h:.4f}")
    print(f"H: status {res['status']}, {res['iterations']} iterations, {res['numBuilds']} builds, "
          f"{res['seconds']:.3f} s (start: {res['startIterations']} iterations of F)")
    print(f"   volume residual {res['volumeResidual']:.2e}, shift residual {res['shiftResidual']:.2e}")
    print(f"check (halfspace cells): max |V - Vbar|/Vbar = {np.max(np.abs(vol - vbar) / vbar):.2e}, "
          f"max |c - cbar|/h = {np.max(np.linalg.norm(cen - cbar, axis=1)) / h:.2e}, "
          f"sum V = {vol.sum():.15f}")
    # Y and p agree up to a homothety about the centroid (H is invariant along it).
    pc, Yc = p - p.mean(axis=0), Y - Y.mean(axis=0)
    sY = np.sum(pc * Yc) / np.sum(pc * pc)
    print(f"   max |p - Y|/h after the best homothety (scale {sY:.6f}) = "
          f"{np.max(np.linalg.norm(sY * pc - Yc, axis=1)) / h:.2e}")
    shift = np.linalg.norm(p - cbar, axis=1) / h
    print(f"   |p - cbar|/h: mean {shift.mean():.3f}, max {shift.max():.3f}; "
          f"generators outside their cells: {sum(not inside(f, p[i]) for i, (_, _, f) in enumerate(cells))}")

    # View.
    az, el = np.radians(AZIM), np.radians(ELEV)
    e = np.array([np.cos(el) * np.cos(az), np.cos(el) * np.sin(az), np.sin(el)])
    right = np.array([-np.sin(az), np.cos(az), 0.0])
    up = np.cross(e, right)
    proj = lambda X: np.c_[X @ right, X @ up]  # noqa: E731
    light = LIGHT / np.linalg.norm(LIGHT)

    removed = np.all(cbar > a.cut, axis=1)
    print(f"cutaway: {removed.sum()} cells removed (cbar > {a.cut} in all coordinates)")
    cols = GL.colors(n, a.seed)

    polys, faces, edges = [], [], []
    for i in np.argsort(p @ e):
        if removed[i]:
            continue
        for nrm, P in cells[i][2]:
            if nrm @ e <= 1e-12:
                continue
            shade = AMBIENT + (1.0 - AMBIENT) * max(0.0, nrm @ light)
            polys.append(proj(P))
            faces.append(np.clip(np.array(GL.tint(cols[i], GL.CELL_TINT)) * shade, 0, 1))
    fig = plt.figure(figsize=(a.size, a.size))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.add_collection(PolyCollection(polys, facecolors=faces, edgecolors="black", linewidths=0.4,
                                     joinstyle="round"))
    # The cube edges at the removed corner.
    c = np.ones(3) * 1.0
    for k in range(3):
        q = c.copy()
        q[k] = 0.0
        S = proj(np.array([q, c]))
        ax.plot(S[:, 0], S[:, 1], color="0.45", linewidth=0.45, dashes=(2, 2), zorder=2)
    # Targets and generators of the removed cells, back to front.
    for i in np.where(removed)[0][np.argsort(p[removed] @ e)]:
        S = proj(np.array([cbar[i], p[i]]))
        dark = tuple(GL.DISK_SHADE * cols[i])
        ax.plot(S[:, 0], S[:, 1], color=dark, linewidth=0.8, zorder=3, solid_capstyle="round")
        ax.plot(*S[0], "s", color=dark, markersize=2.6, markeredgewidth=0, zorder=4)
        ax.plot(*S[1], "o", color="black", markersize=2.4, markeredgewidth=0, zorder=4)
    allp = np.vstack(polys)
    pad = 0.01
    ax.set_xlim(allp[:, 0].min() - pad, allp[:, 0].max() + pad)
    ax.set_ylim(allp[:, 1].min() - pad, allp[:, 1].max() + pad)
    ax.set_aspect("equal")
    ax.axis("off")
    # Tight figure height for the aspect of the drawing.
    wx = np.ptp(allp[:, 0]) + 2 * pad
    wy = np.ptp(allp[:, 1]) + 2 * pad
    fig.set_size_inches(a.size, a.size * wy / wx)
    fig.savefig(a.out)
    print(f"wrote {a.out}")


def inside(facets, x):
    return all(nrm @ (x - P[0]) <= 1e-12 for nrm, P in facets)


if __name__ == "__main__":
    main()
