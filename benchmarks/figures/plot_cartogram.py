#!/usr/bin/env python3
"""Draws the cartogram panel of the showcase figure of the introduction (latex/intro.tex,
fig:intro): the lower 48 US states (DC with Maryland) as a power diagram of the outline of the lower
48, with areas proportional to their 2020 population and centroids as close as possible
to their 2020 centres of population, by H_alpha (introFigure). No power diagram has these
areas and centroids exactly, so the fit keeps the areas exact and moves the centroids.

Data (public domain), in data/:
  CenPop2020_Mean_ST.txt  US Census Bureau, centers of population by state (2020):
                          population, latitude and longitude of the centre
  us-states.json          state boundaries (Leaflet choropleth example), only for the
                          outline: the union of the states, i.e. their unshared edges
Coordinates are longitude cos(39 deg) and latitude (degrees). The centres are kept as
they are, although they do not verify eq:partition-feasibility-centroid with these areas;
the sites are the targets, q = cbar, and alpha is a fraction of the practical bound
(eq:regularized-practical).
  latex/figures/intro-cartogram.pdf
Needs numpy, matplotlib and the introFigure tool (ninja -C build introFigure):
  external/bourne/.venv/bin/python -I benchmarks/figures/plot_cartogram.py [--alpha -0.99]
"""

import argparse
import json
import math
import os
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import plot_intro as PI  # noqa: E402  (sets MPLCONFIGDIR and the Agg backend)

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.patches import Polygon  # noqa: E402

GL = PI.GL
DATA = os.path.join(HERE, "data")
OUT = os.path.join(PI.ROOT, "latex", "figures", "intro-cartogram.pdf")
SKIP = {"Alaska", "Hawaii", "Puerto Rico"}
COS = math.cos(math.radians(39.0))


def project(lon, lat):
    return np.array([lon * COS, lat])


def outline(path):
    """The largest loop of the union of the state polygons (the lower 48): the edges that
    belong to a single state, chained, in counter-clockwise order."""
    d = json.load(open(path))
    edges = Counter()
    for f in d["features"]:
        if f["properties"]["name"] in SKIP:
            continue
        g = f["geometry"]
        rings = [g["coordinates"][0]] if g["type"] == "Polygon" else [p[0] for p in g["coordinates"]]
        for ring in rings:
            pts = [(round(x, 4), round(y, 4)) for x, y in ring]
            if pts[0] == pts[-1]:
                pts = pts[:-1]
            for a, b in zip(pts, pts[1:] + pts[:1]):
                if a != b:
                    edges[tuple(sorted((a, b)))] += 1
    adj = defaultdict(list)
    for (a, b), k in edges.items():
        if k == 1:
            adj[a].append(b)
            adj[b].append(a)
    assert all(len(v) == 2 for v in adj.values()), "the unshared edges do not form simple loops"
    seen, loops = set(), []
    for s in adj:
        if s in seen:
            continue
        loop, prev, cur = [s], None, s
        seen.add(s)
        while True:
            nxt = [v for v in adj[cur] if v != prev and v not in seen]
            if not nxt:
                break
            prev, cur = cur, nxt[0]
            loop.append(cur)
            seen.add(cur)
        loops.append(np.array([project(*p) for p in loop]))
    L = max(loops, key=lambda P: abs(PI.area(P)))
    return L if PI.area(L) > 0 else L[::-1]


def centres(path):
    """Name, population and projected centre of population of the states kept."""
    rows = []
    with open(path, encoding="utf-8-sig") as f:
        next(f)
        for line in f:
            fp, name, pop, lat, lon = line.strip().split(",")
            if name not in SKIP:
                rows.append((name, float(pop), project(float(lon), float(lat))))
    # DC joins Maryland (population-weighted centre): its centre is a few km from
    # Maryland's, which would bound alpha by eq:regularized-practical to almost zero.
    dc = [r for r in rows if r[0] == "District of Columbia"][0]
    md = [r for r in rows if r[0] == "Maryland"][0]
    rows = [r for r in rows if r[0] not in ("District of Columbia", "Maryland")]
    rows.append(("Maryland + DC", md[1] + dc[1], (md[1] * md[2] + dc[1] * dc[2]) / (md[1] + dc[1])))
    return rows


def clip(poly, a, b):
    """The part of the polygon with a . x <= b (Sutherland-Hodgman: exact areas and centroids
    for a non-convex polygon clipped by a half-plane)."""
    return GL.clip(poly, a, b)


def power_cells(p, lam, domain):
    """The power cells of (p, lam) clipped to the domain polygon."""
    cells = []
    for i in range(len(p)):
        poly = [tuple(x) for x in domain]
        for j in range(len(p)):
            if j != i and poly:
                poly = clip(poly, 2.0 * (p[j] - p[i]), p[j] @ p[j] - p[i] @ p[i] - lam[j] + lam[i])
        cells.append(np.array(poly))
    return cells


def run_tool(domain, cbar, V, alpha, tool):
    n = len(cbar)
    with tempfile.TemporaryDirectory() as tmp:
        src, dst = os.path.join(tmp, "in.txt"), os.path.join(tmp, "out.txt")
        with open(src, "w") as f:
            f.write("n %d\n" % n)
            f.write("domain " + " ".join("%.17g" % x for x in domain.ravel()) + "\n")
            f.write("targets " + " ".join("%.17g" % x for x in cbar.ravel()) + "\n")
            f.write("volumes " + " ".join("%.17g" % x for x in V) + "\n")
            f.write("sitesAtTargets 1\nalpha %.17g\n" % alpha)
            if os.environ.get("CARTO_LOG"):
                f.write("log 1\n")
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


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--alpha", type=float, default=-0.99,
                    help="stiffness of H_alpha (< 0: that fraction of the practical bound)")
    ap.add_argument("--tool", default=os.path.join(PI.ROOT, "build", "introFigure"))
    ap.add_argument("--out", default=OUT)
    ap.add_argument("--dpi", type=int, default=300, help="for raster previews")
    args = ap.parse_args()

    domain = outline(os.path.join(DATA, "us-states.json"))
    shape = os.environ.get("CARTO_DOMAIN", "outline")
    rows = centres(os.path.join(DATA, "CenPop2020_Mean_ST.txt"))
    names = [r[0] for r in rows]
    pop = np.array([r[1] for r in rows])
    cbar = np.array([r[2] for r in rows])
    area = PI.area(domain)
    V = area * pop / pop.sum()
    if shape == "rectangle":
        # a rectangle of the outline's aspect and area, centred at the population-weighted
        # mean of the centres, so that eq:partition-feasibility-centroid holds
        w0 = np.ptp(domain[:, 0]); h0 = np.ptp(domain[:, 1])
        k = math.sqrt(area / (w0 * h0)); m = V @ cbar / area
        hw, hh = 0.5 * k * w0, 0.5 * k * h0
        domain = np.array([[m[0] - hw, m[1] - hh], [m[0] + hw, m[1] - hh], [m[0] + hw, m[1] + hh], [m[0] - hw, m[1] + hh]])
    n = len(V)
    # The population-weighted mean of the centres (near Missouri) is not the centroid of the
    # outline (in Kansas), so the targets violate eq:partition-feasibility-centroid; H_alpha
    # needs no such test, its slacks absorb the mismatch.
    print("n = %d, domain %d vertices, area %.1f; mean of the centres - centroid of the outline (%.2f, %.2f) deg" % (
        n, len(domain), area, *(V @ cbar / area - PI.centroid(domain))))
    rec, msg = run_tool(domain, cbar, V, args.alpha, args.tool)
    print(msg)
    if rec is None:
        raise SystemExit(1)
    p, lam = rec["generators"], rec["weights"]
    cells = power_cells(p, lam, domain)
    A = np.array([abs(PI.area(P)) if len(P) > 2 else 0.0 for P in cells])
    C = np.array([PI.centroid(P) if len(P) > 2 else cbar[i] for i, P in enumerate(cells)])
    err = np.linalg.norm(C - cbar, axis=1)
    print("drawn cells: max |V - Vbar| / Vbar %.1e; |c - cbar| (degrees): median %.2f, max %.2f (%s)" % (
        np.max(np.abs(A - V) / V), np.median(err), err.max(), names[int(np.argmax(err))]))

    cols = GL.colors(n, 7)
    w = domain[:, 0].max() - domain[:, 0].min()
    hgt = domain[:, 1].max() - domain[:, 1].min()
    fig = plt.figure(figsize=(2.0, 2.0 * hgt / w))
    ax = fig.add_axes([0, 0, 1, 1])
    for P, c in zip(cells, cols):
        if len(P) > 2:
            ax.add_patch(Polygon(P, closed=True, facecolor=GL.tint(c, GL.CELL_TINT), edgecolor="black",
                                 linewidth=0.6 * PI.LW, joinstyle="round"))
    for i in range(n):                                  # the centroid residual: the target to the fitted centroid
        ax.plot(*zip(cbar[i], C[i]), color="black", linewidth=0.5, solid_capstyle="round", zorder=3)
    ax.plot(cbar[:, 0], cbar[:, 1], "o", color="black", markersize=1.4, markeredgewidth=0, zorder=4)
    ax.add_patch(Polygon(domain, closed=True, fill=False, edgecolor="black", linewidth=1.2 * PI.LW,
                         joinstyle="round"))
    pad = 0.006 * w
    ax.set_xlim(domain[:, 0].min() - pad, domain[:, 0].max() + pad)
    ax.set_ylim(domain[:, 1].min() - pad, domain[:, 1].max() + pad)
    ax.set_aspect("equal")
    ax.axis("off")
    fig.savefig(args.out, dpi=args.dpi)
    plt.close(fig)
    print("wrote " + args.out)


if __name__ == "__main__":
    main()
