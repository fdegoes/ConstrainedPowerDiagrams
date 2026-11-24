#!/usr/bin/env python3
"""Draws the ablation of the stiffness alpha (fig:alpha-ablation): one synthetic case with
n = 20 cells of equal area, solved by G (the joint dual, H in the code) and by G_alpha at
the guaranteed bound (eq:regularized-distinct) and at 1, 4 and 16 times the practical bound
(eq:regularized-practical), all with alpha at 0.99 times the value, from the same sites;
each column uses the same alpha in both rows, the smaller of the two rows' values.
Two rows: the realizable targets (the centroids of a power diagram, as plot_alpha_example.py),
and the same targets with each centroid moved by up to 0.2 cell sizes in a random direction
and translated back to verify eq:partition-feasibility, which no power diagram realizes.
Each panel shows the cells (colours as fig:intro), the targets as squares and the fitted
centroids as crosses, joined by red segments (the centroid residual); a solve that fails
leaves an empty panel with the targets. Prints the statistics of every
panel (alpha, iterations, time, RMS and max centroid miss over the cell size).
  latex/figures/alpha-ablation-{real,pert}-{G,guaranteed,1,4,16}.pdf
Needs numpy, matplotlib and the introFigure tool (ninja -C build introFigure):
  external/bourne/.venv/bin/python benchmarks/figures/plot_alpha_ablation.py
"""

import argparse
import json
import os
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import plot_intro as PI  # noqa: E402  (sets MPLCONFIGDIR and the Agg backend)
import plot_alpha_example as PA  # noqa: E402

import numpy as np  # noqa: E402
from matplotlib.patches import FancyArrowPatch  # noqa: E402

GL = PI.GL
OUT = os.path.join(PI.ROOT, "latex", "figures", "alpha-ablation.pdf")
COLUMNS = [("G", None), ("guaranteed", "g"), ("1", 1.0), ("4", 4.0), ("16", 16.0)]


def run(cbar, V, q, tool, column):
    """introFigure with given targets, volumes and sites; G or G_alpha for the column."""
    n = len(cbar)
    with tempfile.TemporaryDirectory() as tmp:
        src, dst = os.path.join(tmp, "in.txt"), os.path.join(tmp, "out.txt")
        with open(src, "w") as f:
            f.write("n %d\n" % n)
            f.write("targets " + " ".join("%.17g" % x for x in cbar.ravel()) + "\n")
            f.write("volumes " + " ".join("%.17g" % x for x in V) + "\n")
            f.write("sites " + " ".join("%.17g" % x for x in q.ravel()) + "\n")
            key, k = column
            if k is None:
                f.write("closest 1\n")
            elif k == "g":
                f.write("alpha -0.99\nguaranteed 1\n")
            elif isinstance(k, tuple):                  # ("alpha", value): this stiffness
                f.write("alpha %.17g\n" % k[1])
            else:
                f.write("alpha %.17g\n" % (-0.99 * k))
        r = subprocess.run([tool, src, dst], capture_output=True, text=True)
        msg = (r.stdout + r.stderr).strip()
        if r.returncode != 0:
            return None, msg
        rec = {}
        with open(dst) as f:
            for line in f:
                name, *vals = line.split()
                rec[name] = np.array([float(v) for v in vals])
    for name in ("targets", "sites", "generators", "shifts"):
        rec[name] = rec[name].reshape(n, 2)
    rec["message"] = msg
    return rec, msg


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--n", type=int, default=20)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--perturb", type=float, default=0.2, help="largest centroid move, in cell sizes")
    ap.add_argument("--tool", default=os.path.join(PI.ROOT, "build", "introFigure"))
    ap.add_argument("--out", default=OUT)
    ap.add_argument("--dpi", type=int, default=300)
    args = ap.parse_args()

    # the realizable case of plot_alpha_example.py
    rng = np.random.default_rng(args.seed)
    n, h = args.n, np.sqrt(1.0 / args.n)
    Y = PI.random_square(rng, n, 0.06, 0.5 * h)
    q = np.empty_like(Y)
    for i in range(n):
        while True:
            b = 2 * np.pi * rng.random()
            q[i] = Y[i] + 0.35 * h * np.array([np.cos(b), np.sin(b)])
            if np.all(q[i] >= 0.03) and np.all(q[i] <= 0.97):
                break
    base_rec = PA.run(Y, q, args.tool)
    cbar, V = base_rec["targets"], base_rec["volumes"]
    # the perturbed targets: moved, then translated back so that sum V c = |Omega| c_Omega
    prng = np.random.default_rng(7)
    a, m = 2 * np.pi * prng.random(n), args.perturb * h * prng.random(n)
    pert = cbar + np.column_stack([m * np.cos(a), m * np.sin(a)])
    pert = pert + 0.5 - V @ pert / V.sum()

    rows = {"real": cbar, "pert": pert}
    cols = GL.colors(n, args.seed)
    base, ext = os.path.splitext(args.out)
    table = {}
    # The same alpha in both rows: for each column, the smaller of the two rows' values (their
    # bounds depend on the targets), so that every alpha is within the bounds of both rows.
    def bounds(targets):
        _, msg = run(targets, V, q, args.tool, ("1", 1.0))
        g = float(msg.split("guaranteed bound ")[1].split(",")[0])
        p = float(msg.split("practical bound ")[1].split(",")[0])
        return g, p
    B = {row: bounds(t) for row, t in rows.items()}
    alphas = {key: min(0.99 * (B[r][0] if k == "g" else k * B[r][1]) for r in rows)
              for key, k in COLUMNS if k is not None}
    print("bounds (guaranteed, practical):", B, " alphas:", alphas)
    for row, targets in rows.items():
        recs = {key: run(targets, V, q, args.tool, (key, None if k is None else ("alpha", alphas[key])))
                for key, k in COLUMNS}
        good = [r for r, _ in recs.values() if r is not None]
        wmax = max(np.abs(r["weights"] - r["weights"].mean()).max() for r in good)
        table[row] = {"bounds": "guaranteed bound %g, practical bound %g" % B[row]}
        for key, _ in COLUMNS:
            rec, msg = recs[key]
            path = "%s-%s-%s%s" % (base, row, key, ext)
            if rec is None:
                table[row][key] = {"failed": msg.split(":")[-1].strip()[:80]}

                def points(ax, targets=targets):
                    PI.targets(ax, targets)

                PI.panel(path, [], np.zeros((0, 2)), np.zeros(0), [], args.dpi, points)
                continue
            s = PA.stats(rec)
            if key in alphas:                           # passed to the solver, not printed by it
                s["alpha"] = alphas[key]
            table[row][key] = {k: s[k] for k in ("alpha", "iterations", "seconds", "rms", "max")}
            C = s["centroids"]

            def points(ax, C=C, targets=targets):
                # the centroid residual: from each target (square) to the fitted centroid (cross)
                for i in range(n):
                    ax.plot(*zip(targets[i], C[i]), color=PA.RED, linewidth=0.9, solid_capstyle="round", zorder=6)
                PI.targets(ax, targets)
                ax.plot(C[:, 0], C[:, 1], "x", color="black", markersize=3.0, markeredgewidth=0.7, zorder=9)

            # no weight disks: they are centred on the generators, which are not drawn
            PI.panel(path, s["cells"], np.zeros((0, 2)), np.zeros(0), cols, args.dpi, points)
    print(json.dumps(table, indent=1))
    with open(base + ".json", "w") as f:
        json.dump(table, f, indent=1)
    print("wrote %s-{real,pert}-{G,guaranteed,1,4,16}%s" % (base, ext))


if __name__ == "__main__":
    main()
