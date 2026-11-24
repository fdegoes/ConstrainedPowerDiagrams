#!/usr/bin/env python3
"""Ablation of the safeguards of our H solver (latex/implementation.tex) on the
mid- and high-noise synthetic cases with n <= 1129.

Variants (benchBourne case lines, defaults otherwise):
  full        the solver as in the paper (eta = 1e-3, pairing term R = 1.2,
              truncation 1/2, up to 3 restarts with R doubled)
  no pairing  pairreg 0 2       (R = 0; the solver then makes no restart)
  no trunc.   stretchlimit 0
  no restarts restarts 0 2
  no pairing, no truncation    pairreg 0 2, stretchlimit 0
  no truncation, no restarts   stretchlimit 0, restarts 0 2
  eta 1e-6    eta 1e-6          (eta also sets the distance clamp)
  eta 0       eta 0

Success: status "converged" (volume and shift residuals <= 1e-8 within 100
Newton iterations per attempt). Iterations are the total over the attempts.
Wall-clock time is not recorded.

The runs are cached in benchmarks/figures/cache/ablation.csv (one row per
variant and case); cached rows are not rerun unless --force.

Usage:
  MPLCONFIGDIR=/tmp/mplcfg external/bourne/.venv/bin/python benchmarks/figures/plot_ablation.py
      [--force] [--jobs 4] [--plot-only] [--out latex/figures/ablation.pdf]
"""

import argparse
import collections
import concurrent.futures as cf
import csv
import json
import os
import subprocess
import sys
import tempfile

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
SYN = os.path.join(ROOT, "benchmarks", "synthetic")
sys.path.insert(0, SYN)
import generate as G  # noqa: E402
import run_synthetic as RS  # noqa: E402

CACHE = os.path.join(HERE, "cache", "ablation.csv")
SIZES = [n for n in RS.NS if n <= 1129]
NOISES = (1, 2)
VOLS = ("uniform", "voronoi")
SEEDS = (1, 2, 3, 4, 5)
MAXITER = 100

VARIANTS = collections.OrderedDict([
    ("full", ("full solver", [])),
    ("nopair", ("no\npairing", ["pairreg 0 2"])),
    ("notrunc", ("no\ntruncation", ["stretchlimit 0"])),
    ("norestart", ("no\nrestarts", ["restarts 0 2"])),
    ("nopairtrunc", ("no pairing,\nno trunc.", ["pairreg 0 2", "stretchlimit 0"])),
    ("notruncrestart", ("no trunc.,\nno restarts", ["stretchlimit 0", "restarts 0 2"])),
    ("eta1e-6", (r"$\eta = 10^{-6}$", ["eta 1e-6"])),
    ("eta0", (r"$\eta = 0$", ["eta 0"])),
])

FIELDS = ["variant", "case", "n", "noise", "volumes", "seed", "status", "converged", "iterations", "attempts",
          "attemptIterations", "volumeResidual", "shiftResidual"]

NOISE_NAME = {1: "mid", 2: "high"}
VOL = {"uniform": ("uniform", "#2a78d6"), "voronoi": ("Voronoi", "#eb6834")}
INK = "#2b2a27"
GRID = "#e4e2dc"

plt.rcParams.update({"font.size": 8, "axes.labelsize": 8, "xtick.labelsize": 7, "ytick.labelsize": 7})


def lighten(hex_color, t):
    """Mix the colour with white (t = 0: unchanged, 1: white)."""
    c = np.array([int(hex_color[i:i + 2], 16) for i in (1, 3, 5)]) / 255.0
    return tuple(c + (1 - c) * t)


def run_one(variant, c):
    with tempfile.TemporaryDirectory() as d:
        cp, out = os.path.join(d, "case.txt"), os.path.join(d, "ours.json")
        RS.write_case_file(cp, c, MAXITER, extra=VARIANTS[variant][1])
        p = subprocess.run([RS.BENCH, cp, out], capture_output=True, text=True)
        n, noise, vol, seed = c["_meta"]
        row = {"variant": variant, "case": c["name"], "n": n, "noise": noise, "volumes": vol,
               "seed": seed}
        if p.returncode != 0 or not os.path.exists(out):
            row.update(status="error (exit %d)" % p.returncode, converged=False)
            return row
        r = json.load(open(out))
    row.update(status=r["status"], converged=r["converged"], iterations=r["iterations"],
               attempts=r.get("attempts", 1), attemptIterations=r.get("attemptIterations"),
               volumeResidual=r["volumeResidual"], shiftResidual=r["shiftResidual"])
    return row


def run_all(force, jobs):
    rows = {}
    if os.path.exists(CACHE) and not force:
        for r in csv.DictReader(open(CACHE)):
            rows[(r["variant"], r["case"])] = r
    cases = []
    for n in SIZES:
        for noise in NOISES:
            for vol in VOLS:
                for seed in SEEDS:
                    c = G.make_case(os.path.join(SYN, "results"), n, noise, vol, seed)
                    c.update(_meta=(n, noise, vol, seed))
                    cases.append(c)
    todo = [(v, c) for v in VARIANTS for c in cases if (v, c["name"]) not in rows]
    print("%d runs to do, %d cached" % (len(todo), len(rows)), flush=True)
    os.makedirs(os.path.dirname(CACHE), exist_ok=True)

    def save():
        with open(CACHE, "w", newline="") as f:
            w = csv.DictWriter(f, FIELDS)
            w.writeheader()
            for v in VARIANTS:
                for c in cases:
                    if (v, c["name"]) in rows:
                        w.writerow({k: rows[(v, c["name"])].get(k) for k in FIELDS})

    with cf.ThreadPoolExecutor(jobs) as ex:
        futs = {ex.submit(run_one, v, c): (v, c) for v, c in todo}
        for k, f in enumerate(cf.as_completed(futs)):
            r = f.result()
            rows[(r["variant"], r["case"])] = {k2: str(v2) for k2, v2 in r.items()}
            print("%4d/%d %-10s %-26s %-18s %s it" % (k + 1, len(todo), r["variant"], r["case"], r["status"],
                                                      r.get("iterations")), flush=True)
            if k % 50 == 0:
                save()
    save()
    return [r for r in rows.values() if r["variant"] in VARIANTS]


def summary(rows):
    print("\n| variant | noise | volumes | converged | median it. (converged) | statuses of failures |")
    print("|---|---|---|---|---|---|")
    for v in VARIANTS:
        for noise in NOISES:
            for vol in VOLS + (None,):
                s = [r for r in rows if r["variant"] == v and int(r["noise"]) == noise
                     and (vol is None or r["volumes"] == vol)]
                ok = [int(r["iterations"]) for r in s if r["converged"] == "True"]
                fail = collections.Counter(r["status"] for r in s if r["converged"] != "True")
                print("| %s | %s | %s | %d/%d | %s | %s |" % (
                    v, NOISE_NAME[noise], vol or "all", len(ok), len(s), np.median(ok) if ok else "-",
                    ", ".join("%s: %d" % kv for kv in sorted(fail.items()))))
    print("\nper variant, all cases:")
    for v in VARIANTS:
        s = [r for r in rows if r["variant"] == v]
        ok = [int(r["iterations"]) for r in s if r["converged"] == "True"]
        rest = sum(int(r["attempts"]) > 1 for r in s if r["converged"] == "True")
        print("  %-10s %3d/%d converged (%.1f%%), median it %s, IQR %s, max %s, converged after a restart %d" % (
            v, len(ok), len(s), 100 * len(ok) / len(s), np.median(ok) if ok else "-",
            np.percentile(ok, [25, 75]) if ok else "-", max(ok) if ok else "-", rest))
        # Failures by size
        f = collections.Counter(int(r["n"]) for r in s if r["converged"] != "True")
        if f:
            print("             failures by n:", dict(sorted(f.items())))


def plot(rows, path):
    keys = list(VARIANTS)
    groups = [(noise, vol) for noise in NOISES for vol in VOLS]
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(7.0, 3.6), sharex=True,
                                 gridspec_kw={"height_ratios": [1, 1.2]})
    width = 0.8 / len(groups)
    x = np.arange(len(keys))
    for j, (noise, vol) in enumerate(groups):
        label, color = VOL[vol]
        col = color if noise == 2 else lighten(color, 0.5)
        off = (j - (len(groups) - 1) / 2) * width
        rates, its = [], []
        for v in keys:
            s = [r for r in rows if r["variant"] == v and int(r["noise"]) == noise and r["volumes"] == vol]
            rates.append(100.0 * sum(r["converged"] == "True" for r in s) / max(1, len(s)))
            its.append([int(r["iterations"]) for r in s if r["converged"] == "True"])
        a1.bar(x + off, rates, width * 0.92, color=col, label="%s noise, %s" % (NOISE_NAME[noise], label), zorder=2)
        for xi, r in zip(x + off, rates):
            if r < 99.95:
                a1.text(xi, r + 1.5, "%.0f" % r, ha="center", va="bottom", fontsize=6, color=INK)
        bp = a2.boxplot([i if i else [np.nan] for i in its], positions=x + off, widths=width * 0.8,
                        patch_artist=True, showfliers=True, whis=1.5, manage_ticks=False)
        for b in bp["boxes"]:
            b.set(facecolor=col, edgecolor=INK, lw=0.5)
        for k in ("whiskers", "caps"):
            for l in bp[k]:
                l.set(color=INK, lw=0.5)
        for l in bp["medians"]:
            l.set(color=INK, lw=0.9)
        for l in bp["fliers"]:
            l.set(marker="o", ms=1.6, mfc=col, mec=INK, mew=0.3)
    a1.set_ylabel("converged (%)")
    a2.set_ylabel("iterations (converged)")
    a2.set_yscale("log")
    a2.set_xticks(x)
    a2.set_xticklabels([VARIANTS[v][0] for v in keys], fontsize=7)
    for ax in (a1, a2):
        ax.tick_params(axis="x", length=0)
        ax.set_xlim(-0.5, len(keys) - 0.5)
    a1.set_ylim(0, 112)
    a1.set_yticks([0, 25, 50, 75, 100])
    fig.legend(*a1.get_legend_handles_labels(), fontsize=7, ncol=4, frameon=False, loc="upper center",
               bbox_to_anchor=(0.5, 1.0), handlelength=1.2, columnspacing=1.5)
    for ax in (a1, a2):
        ax.grid(True, axis="y", which="major", color=GRID, lw=0.6)
        ax.set_axisbelow(True)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
    fig.tight_layout(h_pad=0.4, rect=(0, 0, 1, 0.95))
    fig.savefig(path)
    plt.close(fig)
    print("wrote", path)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--force", action="store_true", help="rerun the cached runs")
    ap.add_argument("--jobs", type=int, default=4, help="concurrent benchBourne runs")
    ap.add_argument("--plot-only", action="store_true", help="only read the cache")
    ap.add_argument("--out", default=os.path.join(ROOT, "latex", "figures", "ablation.pdf"))
    args = ap.parse_args()
    if args.plot_only:
        rows = list(csv.DictReader(open(CACHE)))
    else:
        rows = run_all(args.force, args.jobs)
    summary(rows)
    plot(rows, args.out)


if __name__ == "__main__":
    main()
