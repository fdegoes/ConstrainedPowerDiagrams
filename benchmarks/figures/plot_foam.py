#!/usr/bin/env python3
"""The foam of Numerow et al. (2024, Sec. 4.4, Fig. 6; external/VoronoiCellSim, Apache-2.0,
Projects/VoronoiFoam/resource/ronan.jpg and ronan.txt): a photo of a 2D soap foam and its
hand annotation (152 vertices, 75 cells). We fit a power diagram to the areas and centroids
of the annotated cells by H_alpha (introFigure), from sites at the targets, and score it
with their objective: the sum over the annotated vertices of the squared distance to the
vertex the generators define (the triple point of the three cells, or the point where the
boundary of two cells meets the side of the image; the four corners are skipped), in their
units (half the image height; the image centred at the origin, y up).

  latex/figures/foam-{photo,fit}.pdf
Needs numpy, matplotlib and the introFigure tool (ninja -C build introFigure):
  external/bourne/.venv/bin/python -I benchmarks/figures/plot_foam.py [--alpha -0.99 ...]
"""

import argparse
import os
import re
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import plot_intro as PI  # noqa: E402  (sets MPLCONFIGDIR and the Agg backend)

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

GL = PI.GL
RES = os.path.join(PI.ROOT, "external", "VoronoiCellSim", "Projects", "VoronoiFoam", "resource")
OUT = os.path.join(PI.ROOT, "latex", "figures", "foam.pdf")


def load():
    """The image, the annotated vertices in their units, the cells (vertex indices), and the
    half extents (dx, dy) of the image."""
    img = plt.imread(os.path.join(RES, "ronan.jpg"))
    H, W = img.shape[:2]
    s = 0.5 * H                                         # their loadImageNetwork
    txt = open(os.path.join(RES, "ronan.txt")).read()
    head, tail = txt.split("Cells:")
    px = np.array([tuple(map(float, m)) for m in re.findall(r"\((\d+),\s*(\d+)\)", head)])
    verts = np.column_stack([(px[:, 0] - 0.5 * W) / s, -(px[:, 1] - 0.5 * H) / s])
    cells = [list(map(int, l.split(":")[1].split(","))) for l in tail.strip().splitlines() if ":" in l]
    return img, verts, cells, (0.5 * W / s, 0.5 * H / s)


def run(domain, cbar, V, alpha, tool, max_iter=100):
    n = len(cbar)
    with tempfile.TemporaryDirectory() as tmp:
        src, dst = os.path.join(tmp, "in.txt"), os.path.join(tmp, "out.txt")
        with open(src, "w") as f:
            f.write("n %d\n" % n)
            f.write("domain " + " ".join("%.17g" % x for x in domain.ravel()) + "\n")
            f.write("targets " + " ".join("%.17g" % x for x in cbar.ravel()) + "\n")
            f.write("volumes " + " ".join("%.17g" % x for x in V) + "\n")
            f.write("sitesAtTargets 1\nalpha %.17g\nmaxIterations %d\n" % (alpha, max_iter))
        r = subprocess.run([tool, src, dst], capture_output=True, text=True)
        if r.returncode != 0:
            return None, (r.stdout + r.stderr).strip()
        rec = {}
        with open(dst) as f:
            for line in f:
                key, *vals = line.split()
                rec[key] = np.array([float(v) for v in vals])
    for k in ("targets", "sites", "generators", "shifts"):
        rec[k] = rec[k].reshape(n, 2)
    return rec, r.stdout.strip()


def their_objective(verts, cells, p, lam, ext):
    """Their image-match objective for the power diagram (p, lam): |x - p_i|^2 - lam_i."""
    dx, dy = ext
    owners = [[] for _ in verts]
    for i, c in enumerate(cells):
        for v in c:
            owners[v].append(i)
    total, count = 0.0, 0
    for v, o in enumerate(owners):
        o = sorted(o)
        if len(o) == 3:                                 # the triple point of the three cells
            a, b, c = o
            M = 2.0 * np.array([p[b] - p[a], p[c] - p[a]])
            r = np.array([p[b] @ p[b] - p[a] @ p[a] - lam[b] + lam[a], p[c] @ p[c] - p[a] @ p[a] - lam[c] + lam[a]])
            x = np.linalg.solve(M, r)
        elif len(o) == 2:                               # the boundary of the two cells on the nearest side
            a, b = o
            X, Y = verts[v]
            d = {"bottom": abs(Y + dy), "right": abs(X - dx), "top": abs(Y - dy), "left": abs(X + dx)}
            side = min(d, key=d.get)
            n_ = 2.0 * (p[b] - p[a])
            r = p[b] @ p[b] - p[a] @ p[a] - lam[b] + lam[a]
            if side in ("bottom", "top"):
                y = -dy if side == "bottom" else dy
                x = np.array([(r - n_[1] * y) / n_[0], y])
            else:
                xx = dx if side == "right" else -dx
                x = np.array([xx, (r - n_[0] * xx) / n_[1]])
        else:
            continue
        total += float(np.sum((x - verts[v]) ** 2))
        count += 1
    return total, count


def clip_cells(p, lam, domain):
    cells = []
    for i in range(len(p)):
        poly = [tuple(x) for x in domain]
        for j in range(len(p)):
            if j != i and poly:
                poly = GL.clip(poly, 2.0 * (p[j] - p[i]), p[j] @ p[j] - p[i] @ p[i] - lam[j] + lam[i])
        cells.append(np.array(poly))
    return cells


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--alpha", type=float, nargs="+", default=[1e-6, -0.99, -3, -10],
                    help="stiffnesses of H_alpha (< 0: that fraction of the practical bound); the first is the baseline")
    ap.add_argument("--max-iter", type=int, default=100, help="Newton iterations per attempt")
    ap.add_argument("--draw", type=float, default=None, help="the alpha of the drawn fit (default: the last that converges)")
    ap.add_argument("--tool", default=os.path.join(PI.ROOT, "build", "introFigure"))
    ap.add_argument("--out", default=OUT)
    ap.add_argument("--dpi", type=int, default=300)
    args = ap.parse_args()

    img, verts, cells, (dx, dy) = load()
    domain = np.array([[-dx, -dy], [dx, -dy], [dx, dy], [-dx, dy]])
    polys = [verts[c] for c in cells]
    V = np.array([abs(PI.area(P)) for P in polys])
    cbar = np.array([PI.centroid(P) for P in polys])
    print("n = %d cells, annotated area %.4f of %.4f" % (len(V), V.sum(), 4 * dx * dy))

    results = []
    for a in args.alpha:
        rec, msg = run(domain, cbar, V, a, args.tool, args.max_iter)
        if rec is None:
            print("alpha %g: failed: %s" % (a, msg))
            continue
        L, m = their_objective(verts, cells, rec["generators"], rec["weights"], (dx, dy))
        it = int(msg.split(" iterations")[0].split()[-1])
        sec = float(msg.split(" s (start")[0].split()[-1])
        al = float(msg.split("alpha ")[1].split(";")[0]) if "alpha " in msg else a
        C = np.array([PI.centroid(P) for P in clip_cells(rec["generators"], rec["weights"], domain)])
        miss = np.linalg.norm(C - cbar, axis=1) / np.sqrt(V)
        print("alpha %-9.4g %2d iterations %.4f s  their objective %.4f (%d vertices)  centroid miss / cell size: median %.3f max %.3f" % (
            al, it, sec, L, m, np.median(miss), miss.max()))
        results.append((al, rec, L))
    al, rec, L = results[-1] if args.draw is None else min(results, key=lambda r: abs(r[0] - args.draw))

    base, ext = os.path.splitext(args.out)
    # Drawn: the photo inside its black border (pixels 50 to 877 x 50 to 1157); the solve
    # uses their whole domain.
    H_, W_ = img.shape[:2]
    s = 0.5 * H_
    crop = (((50 - 0.5 * W_) / s, (878 - 0.5 * W_) / s), (-(1158 - 0.5 * H_) / s, -(50 - 0.5 * H_) / s))
    figsize = (2.2, 2.2 * (crop[1][1] - crop[1][0]) / (crop[0][1] - crop[0][0]))
    for name, draw_fit in (("photo", False), ("fit", True)):
        fig = plt.figure(figsize=figsize)
        ax = fig.add_axes([0, 0, 1, 1])
        ax.imshow(img, extent=(-dx, dx, -dy, dy), origin="upper", interpolation="bilinear")
        if draw_fit:
            for P in clip_cells(rec["generators"], rec["weights"], domain):
                Q = np.vstack([P, P[:1]])
                ax.plot(Q[:, 0], Q[:, 1], color="#38d430", linewidth=0.8, solid_joinstyle="round")
        else:
            for c in cells:
                Q = verts[c + c[:1]]
                ax.plot(Q[:, 0], Q[:, 1], color="#f0c020", linewidth=0.8, solid_joinstyle="round")
        ax.set_xlim(*crop[0])
        ax.set_ylim(*crop[1])
        ax.set_aspect("equal")
        ax.axis("off")
        fig.savefig(base + "-" + name + ext, dpi=args.dpi)
        plt.close(fig)
    print("drew alpha %.4g (their objective %.4f): %s-{photo,fit}%s" % (al, L, base, ext))


if __name__ == "__main__":
    main()
