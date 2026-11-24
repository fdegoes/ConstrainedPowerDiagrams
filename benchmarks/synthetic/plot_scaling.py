#!/usr/bin/env python3
"""Time and iterations of our H against the number of sites n, on the 600
2D synthetic cases (log-synthetic.csv of run_overnight.sh), and the summary
statistics for the text.

One panel per noise level (low, mid, high), the same y-scale in all panels.
In each panel, the Voronoi volumes (circles, solid) and the uniform volumes
(squares, dashed) are drawn as the median over the seeds at each n, with a
shaded band between the 25th and 75th percentiles. The lines are labeled at
their right end. The time figure has a grey reference line proportional to n.
The figures are drawn at the width of the text (the paper's figure*), so the
type is at its printed size. Time is in milliseconds.

Usage:
  external/bourne/.venv/bin/python benchmarks/synthetic/plot_scaling.py
      [--log benchmarks/overnight/2026-11-20/log-synthetic.csv]
      [--out latex/figures/synthetic] [--md stats.md]
"""

import argparse
import collections
import csv
import os
import statistics

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))

# Colour is the volume: blue for uniform, orange for Voronoi (two hues, separable
# for colour-vision deficiencies). Each panel is one noise level.
NOISE = {"0": ("low", None), "1": ("mid", None), "2": ("high", None)}  # panel titles
VOLUMES = {"voronoi": ("o", "-", "Voronoi", "#eb6834"), "uniform": ("s", "--", "uniform", "#2a78d6")}  # marker, line style, label, colour
INK = "#2b2a27"  # text ink, not a series color

GRID = "#e4e2dc"
TOL = 1e-8

plt.rcParams.update({"font.size": 8, "axes.labelsize": 8, "xtick.labelsize": 7, "ytick.labelsize": 7})


def load(path):
    return list(csv.DictReader(open(path)))


def ours(rows):
    return [r for r in rows if r["method"] == "our H" and r["converged"] == "True"]


def groups(rows):
    """{(noise, volumes): {n: [(ms, iterations), ...]}} for our H."""
    g = collections.defaultdict(lambda: collections.defaultdict(list))
    for r in rows:
        g[(r["noise"], r["volumes"])][int(r["n"])].append((1000 * float(r["seconds"]), int(r["iterations"])))
    return g


def band(values):
    """Median, 25th and 75th percentiles of the seeds."""
    q1, med, q3 = np.percentile(values, [25, 50, 75])
    return med, q1, q3


def plot(g, index, ylabel, path, reference):
    fig, axes = plt.subplots(1, 3, figsize=(7.0, 2.5), sharey=True)
    for ax, (noise, (name, _)) in zip(axes, NOISE.items()):
        for vol, (marker, ls, label, color) in VOLUMES.items():
            by_n = g[(noise, vol)]
            ns = np.array(sorted(by_n))
            stats = np.array([band([v[index] for v in by_n[n]]) for n in ns])
            med, lo, hi = stats[:, 0], stats[:, 1], stats[:, 2]
            c = color
            ax.fill_between(ns, lo, hi, color=c, alpha=0.2, lw=0, zorder=1)
            ax.plot(ns, med, ls, color=c, lw=1.4, zorder=3)
            ax.plot(ns, med, marker, color=c, ms=3.5, mec="white", mew=0.4, lw=0, zorder=4)
            ax.annotate(label, (ns[-1], med[-1]), xytext=(4, 0), textcoords="offset points",
                        va="center", ha="left", fontsize=7, color=INK)
        if reference and index == 0:
            # Linear scaling, through the pooled median at the size nearest to 1000.
            sizes = sorted({n for key in g if key[0] == noise for n in g[key]})
            n0 = min(sizes, key=lambda n: abs(n - 1000))
            c = statistics.median(v[0] for key in g if key[0] == noise for v in g[key][n0]) / n0
            xs = np.array([10, 10000])
            ax.plot(xs, c * xs, color="#8a8880", lw=0.9, ls=(0, (1, 1.5)), zorder=0)
            ax.annotate("$\\propto n$", (25, c * 25), xytext=(3, -8), textcoords="offset points",
                        fontsize=7, color="#8a8880")
        ax.set_xscale("log")
        ax.set_xlim(8, 40000)
        if index == 0:
            ax.set_yscale("log")
        ax.set_xlabel("number of sites $n$")
        ax.text(0.04, 0.96, f"{name} noise", transform=ax.transAxes, va="top", ha="left",
                fontsize=8, color=INK)
        ax.grid(True, which="major", color=GRID, lw=0.6)
        ax.set_axisbelow(True)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
    axes[0].set_ylabel(ylabel)
    fig.tight_layout(w_pad=1.2)
    fig.savefig(path)
    plt.close(fig)


NOISE_STYLE = {"0": ("-", "o", "low"), "1": ("-", "s", "mid"), "2": ("-", "^", "high")}  # line, marker, label


def plot_combined(g, path):
    """One column-wide figure: time (top) and iterations (bottom) against n, all six groups in each
    panel. Hue is the volume type and the marker the noise level, with one legend above both panels;
    the medians over the seeds only (the bands of six overlapping groups would blur)."""
    from matplotlib.lines import Line2D
    fig, axes = plt.subplots(2, 1, figsize=(3.3, 4.0), sharex=True)
    for ax, index, ylabel in ((axes[0], 0, "time (ms)"), (axes[1], 1, "iterations")):
        for vol, (_, _, _, color) in VOLUMES.items():
            for noise, (ls, marker, _) in NOISE_STYLE.items():
                by_n = g[(noise, vol)]
                ns = np.array(sorted(by_n))
                med = np.array([band([v[index] for v in by_n[n]])[0] for n in ns])
                ax.plot(ns, med, ls, color=color, lw=1.2, alpha=0.85, zorder=3)
                ax.plot(ns, med, marker, color=color, ms=3.6, mec="white", mew=0.4, lw=0, zorder=4)
        if index == 0:
            # Linear scaling, through the pooled median at the size nearest to 1000.
            sizes = sorted({n for key in g for n in g[key]})
            n0 = min(sizes, key=lambda n: abs(n - 1000))
            c = statistics.median(v[0] for key in g for v in g[key][n0]) / n0
            xs = np.array([10, 10000])
            ax.plot(xs, c * xs, color="#8a8880", lw=0.9, ls=(0, (1, 1.5)), zorder=0)
            ax.annotate("$\\propto n$", (25, c * 25), xytext=(3, -8), textcoords="offset points",
                        fontsize=7, color="#8a8880")
            ax.set_yscale("log")
        ax.set_xscale("log")
        ax.set_xlim(8, 14000)
        if index == 1:                                  # the shared n axis, labelled under the bottom panel only
            ax.set_xlabel("number of sites $n$")
        ax.set_ylabel(ylabel)
        ax.grid(True, which="major", color=GRID, lw=0.6)
        ax.set_axisbelow(True)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
    # one legend for both panels, above them: the volume types (hue) and the noise levels (marker)
    handles = [Line2D([], [], color=c, lw=1.6, label=lab) for _, (_, _, lab, c) in reversed(list(VOLUMES.items()))]
    handles += [Line2D([], [], color="#6b6a66", ls="none", marker=m, ms=3.8, label=lab)
                for _, m, lab in NOISE_STYLE.values()]
    fig.tight_layout(h_pad=0.6, rect=(0, 0, 1, 0.95))
    # the legend centred over the plot boxes (not the figure, whose left part holds the y labels),
    # and the saved figure cropped to its content, so that the paper centres the plots
    box = axes[0].get_position()
    fig.legend(handles=handles, frameon=False, fontsize=7, loc="upper center", ncol=5, handlelength=1.4,
               columnspacing=1.0, handletextpad=0.4, bbox_to_anchor=(0.5 * (box.x0 + box.x1), 1.0))
    fig.savefig(path, bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)


def slope(g, key, index):
    """Log-log slope of the median against n over the sizes of the group."""
    ns = sorted(g[key])
    x = np.log([n for n in ns])
    y = np.log([statistics.median(v[index] for v in g[key][n]) for n in ns])
    return np.polyfit(x, y, 1)[0]


def summary(rows, g):
    out = []
    ok = ours(rows)
    total = [r for r in rows if r["method"] == "our H"]
    out.append(f"Our H converged on {len(ok)} of {len(total)} synthetic cases "
               f"(rule 3 at {TOL:g}); the iterations range from {min(int(r['iterations']) for r in ok)} "
               f"to {max(int(r['iterations']) for r in ok)}, with a median of "
               f"{statistics.median(int(r['iterations']) for r in ok):g}.")
    out.append("")
    out.append("| n | median time (ms) | median iterations |")
    out.append("|---|---|---|")
    for n in sorted({int(r["n"]) for r in ok}):
        t = [1000 * float(r["seconds"]) for r in ok if int(r["n"]) == n]
        it = [int(r["iterations"]) for r in ok if int(r["n"]) == n]
        out.append(f"| {n} | {statistics.median(t):.3g} | {statistics.median(it):g} |")
    out.append("")
    out.append("Log-log slope of the median time against n (time ~ n^p), per group:")
    for noise, (name, _) in NOISE.items():
        for vol in VOLUMES:
            out.append(f"- {name} noise, {vol}: p = {slope(g, (noise, vol), 0):.2f}")
    pooled = collections.defaultdict(list)
    for r in ok:
        pooled[int(r["n"])].append(1000 * float(r["seconds"]))
    ns = sorted(pooled)
    p = np.polyfit(np.log(ns), np.log([statistics.median(pooled[n]) for n in ns]), 1)[0]
    out.append(f"- all cases pooled: p = {p:.2f}")
    out.append("")
    return out


def comparison(rows):
    """Our H and their methods on the synthetic cases they ran: FW and the subgradient on all cases, SLSQP on
    seed 1 and n <= 500. Their times are to their own stop; our residual is checked at their solution."""
    out = ["| method | stopped by | converged | median time (ms) | our residual met at its solution |",
           "|---|---|---|---|---|"]
    for method in ["our H", "their SLSQP", "their FW (t = 100)", "their subgradient (t = 100)"]:
        s = [r for r in rows if r["method"] == method
             and (method != "their SLSQP" or (r["seed"] == "1" and int(r["n"]) <= 500))]
        if not s:
            continue
        stops = collections.Counter(r["stoppedBy"] for r in s)
        conv = sum(r["converged"] == "True" for r in s)
        t = [1000 * float(r["seconds"]) for r in s]
        if method == "our H":
            met = "by construction"
        else:
            m = sum(float(r["volumeResidual"]) <= TOL and float(r["shiftResidual"]) <= TOL for r in s)
            met = f"{m}/{len(s)}"
        stop_txt = ", ".join(f"{k}: {v}" for k, v in sorted(stops.items()))
        out.append(f"| {method} | {stop_txt} | {conv}/{len(s)} | {statistics.median(t):.3g} | {met} |")
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--log", default=os.path.join(ROOT, "benchmarks/overnight/2026-11-20/log-synthetic.csv"))
    ap.add_argument("--out", default=os.path.join(ROOT, "latex/figures/synthetic"))
    ap.add_argument("--md", default=None, help="also write the statistics to this markdown file")
    ap.add_argument("--panels", action="store_true", help="also write the per-noise-level figures")
    args = ap.parse_args()

    rows = load(args.log)
    g = groups(ours(rows))
    os.makedirs(args.out, exist_ok=True)
    plot_combined(g, os.path.join(args.out, "scaling.pdf"))
    if args.panels:   # the former two figures, one panel per noise level
        plot(g, 0, "time (ms)", os.path.join(args.out, "scaling-time.pdf"), reference=True)
        plot(g, 1, "iterations", os.path.join(args.out, "scaling-iterations.pdf"), reference=False)

    lines = ["## Our H on the synthetic benchmark", ""] + summary(rows, g)
    lines += ["## Matched cases (seed 1, n <= 500)", ""] + comparison(rows) + [""]
    text = "\n".join(lines)
    print(text)
    if args.md:
        open(args.md, "w").write(text)


if __name__ == "__main__":
    main()
