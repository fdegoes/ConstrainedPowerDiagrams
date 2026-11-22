#!/usr/bin/env python3
"""Time and iterations of our H against the number of sites n, on the 600
2D synthetic cases (log-synthetic.csv of run_overnight.sh), and the summary
statistics for the text.

One color per noise level (low, mid, high). One marker per volume type:
a dot for the Voronoi volumes, a square for the uniform volumes. Every case
is a marker (one per seed); a line joins the medians over the seeds at each n,
solid for Voronoi and dashed for uniform. Time is in milliseconds.

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
from matplotlib.lines import Line2D  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))

# Categorical slots 1-3 of the dataviz reference palette (as plot_timing.py).
NOISE = {"0": ("low", "#2a78d6"), "1": ("mid", "#eb6834"), "2": ("high", "#1baf7a")}
VOLUMES = {"voronoi": ("o", "-"), "uniform": ("s", "--")}  # marker, line style
TOL = 1e-8


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


def plot(g, index, ylabel, title, path):
    fig, ax = plt.subplots(figsize=(4.6, 3.4))
    for noise, (name, color) in NOISE.items():
        for vol, (marker, ls) in VOLUMES.items():
            by_n = g[(noise, vol)]
            ns = sorted(by_n)
            for n in ns:
                ax.scatter([n] * len(by_n[n]), [v[index] for v in by_n[n]], marker=marker, s=12,
                           color=color, alpha=0.55, linewidths=0, zorder=2)
            meds = [statistics.median(v[index] for v in by_n[n]) for n in ns]
            ax.plot(ns, meds, ls, color=color, lw=1.4, zorder=3)
    ax.set_xscale("log")
    if index == 0:
        ax.set_yscale("log")
    ax.set_xlabel("number of sites $n$")
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontsize=9)
    ax.grid(True, which="major", color="#e4e2dc", lw=0.6)
    ax.set_axisbelow(True)
    handles = [Line2D([], [], color=c, lw=0, marker="o", ms=5, label=f"{name} noise")
               for name, c in NOISE.values()]
    handles += [Line2D([], [], color="#52514e", lw=0, marker=m, ms=5, label=f"{v} volumes")
                for v, (m, _) in VOLUMES.items()]
    handles += [Line2D([], [], color="#52514e", lw=1.4, ls=ls, label=f"{v}, median")
                for v, (_, ls) in VOLUMES.items()]
    ax.legend(handles=handles, fontsize=7, frameon=False, loc="upper left", bbox_to_anchor=(1.01, 1.0))
    fig.tight_layout()
    fig.savefig(path)
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
    """Matched cases of Bourne et al. (seed 1, n <= 500), with our H on the same
    cases. Their times are to their own stop; the rule-3 check is at their solution."""
    out = ["| method | stopped by | converged | median time (ms) | rule 3 met at its solution |",
           "|---|---|---|---|---|"]
    matched = [r for r in rows if r["seed"] == "1" and int(r["n"]) <= 500]
    for method in ["our H", "their SLSQP", "their FW (t = 100)"]:
        s = [r for r in matched if r["method"] == method]
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
    args = ap.parse_args()

    rows = load(args.log)
    g = groups(ours(rows))
    os.makedirs(args.out, exist_ok=True)
    plot(g, 0, "time (ms)", "Our H, time to converge", os.path.join(args.out, "scaling-time.pdf"))
    plot(g, 1, "iterations", "Our H, iterations to converge", os.path.join(args.out, "scaling-iterations.pdf"))

    lines = ["## Our H on the synthetic benchmark", ""] + summary(rows, g)
    lines += ["## Matched cases (seed 1, n <= 500)", ""] + comparison(rows) + [""]
    text = "\n".join(lines)
    print(text)
    if args.md:
        open(args.md, "w").write(text)


if __name__ == "__main__":
    main()
