#!/usr/bin/env python3
"""Draws the softening figure (fig:alpha-example): one realizable synthetic case solved by H
and by H_alpha at increasing stiffness alpha (introFigure), drawn as fig:intro (cells and
schematic weight disks as in the synthetic gallery, targets as squares, sites as circles,
generators as dots, shifts as arrows), with the centroid misses of H_alpha as red ticks.

The generators Y of the targets are dart-thrown in the unit square (mid noise, 0.5 h apart),
the volumes are uniform, the targets are the centroids of the diagram of Y with these
volumes, and the sites are q = Y + a random offset of 0.35 h (as fig:intro), so that the
solution of H is the diagram of Y, while that of H_alpha also depends on q.
  latex/figures/alpha-example-{H,guaranteed,practical,large}.pdf
Needs numpy, matplotlib and the introFigure tool (ninja -C build introFigure):
  external/bourne/.venv/bin/python benchmarks/figures/plot_alpha_example.py [--n 20] [--seed 1]
"""

import argparse
import os
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import plot_intro as PI  # noqa: E402  (sets MPLCONFIGDIR and the Agg backend)

import numpy as np  # noqa: E402
from matplotlib.patches import FancyArrowPatch  # noqa: E402

GL = PI.GL
OUT = os.path.join(PI.ROOT, "latex", "figures", "alpha-example.pdf")
RED = "#c83232"


def run(Y, q, tool, alpha=None, guaranteed=False):
    """introFigure from the sites q with the targets of the diagram of Y; H if alpha is None,
    else H_alpha (alpha < 0: that fraction of a bound)."""
    n = len(Y)
    with tempfile.TemporaryDirectory() as tmp:
        src, dst = os.path.join(tmp, "in.txt"), os.path.join(tmp, "out.txt")
        with open(src, "w") as f:
            f.write("n %d\n" % n)
            f.write("generators " + " ".join("%.17g" % x for x in Y.ravel()) + "\n")
            f.write("sites " + " ".join("%.17g" % x for x in q.ravel()) + "\n")
            if alpha is None:
                f.write("closest 1\n")
            else:
                f.write("alpha %.17g\n" % alpha)
                if guaranteed:
                    f.write("guaranteed 1\n")
        r = subprocess.run([tool, src, dst], capture_output=True, text=True)
        if r.returncode != 0:
            raise SystemExit("introFigure failed: " + (r.stdout + r.stderr).strip())
        rec = {}
        with open(dst) as f:
            for line in f:
                key, *vals = line.split()
                rec[key] = np.array([float(v) for v in vals])
    for k in ("targets", "sites", "generators", "shifts"):
        rec[k] = rec[k].reshape(n, 2)
    rec["message"] = r.stdout.strip()
    return rec


def stats(rec):
    """Iterations, seconds, alpha and the centroid misses (RMS and max, over the cell size)
    from the tool's output and the drawn cells."""
    msg = rec["message"]
    it = int(msg.split(" iterations")[0].split()[-1])
    sec = float(msg.split(" s (start")[0].split()[-1])
    alpha = float(msg.split("alpha ")[1].split(";")[0]) if "alpha " in msg else float("inf")
    cells = GL.power_cells(rec["generators"], rec["weights"])
    C = np.array([PI.centroid(P) for P in cells])
    V = np.array([PI.area(P) for P in cells])
    n = len(V)
    miss = np.linalg.norm(C - rec["targets"], axis=1) * np.sqrt(n)
    return dict(iterations=it, seconds=sec, alpha=alpha, rms=float(np.sqrt(np.mean(miss ** 2))),
                max=float(miss.max()), volume=float(np.max(np.abs(V - rec["volumes"]) / rec["volumes"])),
                cells=cells, centroids=C)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--n", type=int, default=20)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--large", type=float, default=30.0, help="the last alpha, in practical bounds")
    ap.add_argument("--tool", default=os.path.join(PI.ROOT, "build", "introFigure"))
    ap.add_argument("--out", default=OUT)
    ap.add_argument("--dpi", type=int, default=300, help="for raster previews")
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)
    h = np.sqrt(1.0 / args.n)
    Y = PI.random_square(rng, args.n, 0.06, 0.5 * h)
    q = np.empty_like(Y)
    for i in range(args.n):                             # sites inside the square, off its boundary
        while True:
            b = 2 * np.pi * rng.random()
            q[i] = Y[i] + 0.35 * h * np.array([np.cos(b), np.sin(b)])
            if np.all(q[i] >= 0.03) and np.all(q[i] <= 0.97):
                break

    runs = [("H", run(Y, q, args.tool)),
            ("guaranteed", run(Y, q, args.tool, -0.99, guaranteed=True)),
            ("practical", run(Y, q, args.tool, -0.99)),
            ("large", run(Y, q, args.tool, -args.large))]
    print(runs[1][1]["message"].split(";")[0])
    S = {k: stats(r) for k, r in runs}
    for k, r in runs:
        s = S[k]
        print("%-10s alpha %-9.4g %2d iterations, %.4f s, centroid miss / h: RMS %.2e, max %.2e, volume error %.1e" % (
            k, s["alpha"], s["iterations"], s["seconds"], s["rms"], s["max"], s["volume"]))

    # one weight scale for the four panels, as the gallery
    lams = [r["weights"] for _, r in runs]
    wmax = max(np.abs(l - l.mean()).max() for l in lams)
    cols = GL.colors(args.n, args.seed)
    base, ext = os.path.splitext(args.out)
    for k, r in runs:
        p, lam, cbar = r["generators"], r["weights"], r["targets"]
        radii = h / PI.H_GALLERY * (GL.DISK_R0 + GL.DISK_R1 * (lam - lam.mean()) / wmax)
        C = S[k]["centroids"]

        def points(ax, r=r, p=p, cbar=cbar, C=C, k=k):
            for i in range(len(p)):
                ax.add_patch(FancyArrowPatch(r["sites"][i], p[i], arrowstyle="-|>", mutation_scale=5, shrinkA=1.2,
                                             shrinkB=1.2, color="black", linewidth=0.45, zorder=5))
                if k != "H":                            # the centroid miss: target to fitted centroid
                    ax.plot(*zip(cbar[i], C[i]), color=RED, linewidth=0.9, solid_capstyle="round", zorder=6)
            PI.sites(ax, r["sites"])
            PI.targets(ax, cbar)
            ax.plot(p[:, 0], p[:, 1], "o", color="black", markersize=2.2, markeredgewidth=0, zorder=9, clip_on=False)

        PI.panel(base + "-" + k + ext, S[k]["cells"], p, radii, cols, args.dpi, points)
    print("wrote %s-{H,guaranteed,practical,large}%s" % (base, ext))


if __name__ == "__main__":
    main()
