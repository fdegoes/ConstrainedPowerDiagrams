#!/usr/bin/env python3
"""Timing of the solvers on the 2D synthetic benchmark against n.

One color per (volumes, noise) pair. Every converged case is a dot (one per
seed), and a line joins the medians over the seeds at each n. Line styles
tell the methods apart:
  solid   our H (results/ours.csv, run_synthetic.py), time to converge
  dashed  their SLSQP on H, time to stop on rule 1 (its own)
  dotted  their FW (t = 100), time to stop on their FW gap (their tolerance);
          it stalls at about 1/t in our residual, so it is usually not drawn
  dash-dot their subgradient method, time to stop on their rule (change of the
          running average below 1e-5)
Their methods come from <results>/theirs.csv (run_theirs.py) and are drawn
with --theirs. A (method, pair) curve is drawn only where at least
--min-converged of its cases met the rule; the cases that did not are
counted in the printed summary, so that a comparison of the timings that
cannot be drawn is reported in the text or a table instead.

Usage:
  external/bourne/.venv/bin/python benchmarks/synthetic/plot_timing.py
      [--results benchmarks/synthetic/results] [--theirs]
      [--out latex/figures/synthetic/timing.pdf]
"""

import argparse
import collections
import csv
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.ticker import NullFormatter  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))

def shade(color, amount):
    """The colour mixed with white by amount > 0 (lighter), or with black by -amount (darker)."""
    target = 255 if amount > 0 else 0
    rgb = [int(color[i:i + 2], 16) for i in (1, 3, 5)]
    return "#" + "".join("%02x" % round(c + (target - c) * abs(amount)) for c in rgb)

# Two hues for the volumes (blue uniform, orange Voronoi: the two most
# separable hues for colour-vision deficiencies), and within each hue three
# lightness steps for the noise level (dark low, base mid, light high). Markers
# give a second encoding of the volumes.
PAIRS = [("uniform", 0), ("uniform", 1), ("uniform", 2), ("voronoi", 0), ("voronoi", 1), ("voronoi", 2)]
VOLUME_HUE = {"uniform": "#2a78d6", "voronoi": "#eb6834"}
NOISE_SHADE = {0: -0.3, 1: 0.0, 2: 0.4}
COLORS = [shade(VOLUME_HUE[vol], NOISE_SHADE[noise]) for vol, noise in PAIRS]
MARKERS = ["o", "s", "^", "D", "v", "P"]
NOISE_NAME = {0: "low", 1: "mid", 2: "high"}
METHODS = [("ours", "H (ours)", "-"), ("slsqp", "SLSQP (Bourne 2025)", "--"),
           ("fw", "Frank–Wolfe (Bourne 2026)", ":"), ("subgrad", "Subgradient (Bourne 2026)", "-.")]


def read(path):
    return list(csv.DictReader(open(path))) if os.path.exists(path) else []


def seconds(r, method, tol):
    """Seconds to stop on the method's own rule, or None if it did not:
    ours on rule 3; SLSQP on rule 1; FW and the subgradient on their own rules (their gap, their change)."""
    if method == "ours":
        return float(r["seconds"]) if r["converged"] == "True" else None
    rule = "rule 1" if method == "slsqp" else "their rule"
    return float(r["seconds"]) if r.get("stopped") == rule else None


def collect(rows, method, tol):
    """{pair: {n: [seconds of converged cases]}}, and the counts (met, all)."""
    data = collections.defaultdict(lambda: collections.defaultdict(list))
    count = collections.defaultdict(lambda: [0, 0])
    for r in rows:
        if method != "ours" and r["method"] != method:
            continue
        pair = (r["volumes"], int(r["noise"]))
        s = seconds(r, method, tol)
        count[pair][1] += 1
        if s is not None:
            count[pair][0] += 1
            data[pair][int(r["n"])].append(s)
    return data, count


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results", default=os.path.join(HERE, "results"))
    ap.add_argument("--theirs", action="store_true", help="add their SLSQP and FW (results/theirs.csv)")
    ap.add_argument("--min-converged", type=float, default=0.5,
                    help="draw a curve of theirs only if this fraction of its cases met the rule")
    ap.add_argument("--out", default=os.path.join(ROOT, "latex", "figures", "synthetic", "timing.pdf"))
    args = ap.parse_args()

    ours = read(os.path.join(args.results, "ours.csv"))
    theirs = read(os.path.join(args.results, "theirs.csv")) if args.theirs else []
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    drawn = set()
    print("method, pair: cases meeting the rule / cases")
    for method, label, style in METHODS:
        rows = ours if method == "ours" else theirs
        if not rows:
            continue
        data, count = collect(rows, method, None)
        for k, pair in enumerate(PAIRS):
            met, total = count[pair]
            if total == 0:
                continue
            print("  %-6s %-8s noise %d: %d / %d" % (method, pair[0], pair[1], met, total))
            if method != "ours" and met < args.min_converged * total:
                continue
            ns = sorted(data[pair])
            if not ns:
                continue
            for n in ns:
                ax.scatter([n] * len(data[pair][n]), data[pair][n], s=10, color=COLORS[k], alpha=0.35,
                           linewidths=0, zorder=2)
            med = [np.median(data[pair][n]) for n in ns]
            ax.plot(ns, med, style, color=COLORS[k], lw=1.6, marker=MARKERS[k], ms=4.5,
                    mfc=COLORS[k] if method == "ours" else "white", mec=COLORS[k], zorder=3)
            drawn.add((method, k))

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.set_xlabel("number of cells $n$")
    ax.set_ylabel("seconds to converge")
    ax.grid(True, which="major", color="#e6e5e1", lw=0.6)
    ax.grid(True, which="minor", color="#f1f0ec", lw=0.4)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)

    pair_handles = [Line2D([], [], color=COLORS[k], marker=MARKERS[k], lw=1.6, ms=4.5,
                           label="%s volumes, %s noise" % (v.capitalize() if v == "uniform" else "Voronoi",
                                                           NOISE_NAME[z]))
                    for k, (v, z) in enumerate(PAIRS) if any((m, k) in drawn for m, _, _ in METHODS)]
    # Legends outside the axes, so that they cover no curve.
    leg = ax.legend(handles=pair_handles, loc="upper left", bbox_to_anchor=(1.02, 1.0), fontsize=7.5, frameon=False)
    methods_drawn = [(m, label, style) for m, label, style in METHODS if any((m, k) in drawn for k in range(6))]
    if len(methods_drawn) > 1:
        ax.add_artist(leg)
        ax.legend(handles=[Line2D([], [], color="#52514e", ls=style, lw=1.6, label=label)
                           for _, label, style in methods_drawn], loc="lower left", bbox_to_anchor=(1.02, 0.0),
                  fontsize=7.5, frameon=False)
    fig.tight_layout()
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    fig.savefig(args.out)
    if not args.out.endswith(".png"):
        fig.savefig(os.path.splitext(args.out)[0] + ".png", dpi=150)
    print("wrote", args.out)


if __name__ == "__main__":
    main()
