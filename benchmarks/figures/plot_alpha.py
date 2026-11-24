#!/usr/bin/env python3
"""What the stiffness bounds on alpha mean in practice: H_alpha (eq:regularized-dual)
maximized over a log sweep of alpha on the approximate (unrealizable) cases of Bourne et
al. (benchmarks/bourne/bourne_cases.py: 2025-perturbed-small, 2025-perturbed-large and
the EBSD map, all in the unit square), from sites at the targets (q = cbar), with the
introFigure tool (it starts from the volume solve at the sites).

Top: the centroid residual of the result, the volume-weighted RMS of |c_i - cbar_i| over
the mean cell size h = sqrt(|Omega| / n); failed solves (damping failed, max iterations)
are drawn as crosses in rows above the data. Bottom: Newton iterations.
Vertical lines: the guaranteed bound eq:regularized-distinct (dashed) and the practical
estimate eq:regularized-practical (dotted) of each case (bourne_cases.distinct_alpha and
practical_alpha).

The solves are cached in benchmarks/figures/data/alpha-sweep.json (--rerun to redo them).
Needs pysdot (external/bourne/.venv) and the introFigure tool (build/introFigure):
  MPLCONFIGDIR=/tmp/mplcfg external/bourne/.venv/bin/python benchmarks/figures/plot_alpha.py
      [--rerun] [--amin 1e-3] [--amax 1e3] [--num 31] [--out latex/figures/alpha-sweep.pdf]
"""

import argparse
import json
import math
import os
import re
import subprocess
import sys
import tempfile

os.environ.setdefault("MPLCONFIGDIR", "/tmp/mplcfg")
import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(ROOT, "benchmarks", "bourne"))
import bourne_cases as BC  # noqa: E402
import bourne_ref as ref  # noqa: E402
from pysdot import PowerDiagram  # noqa: E402

CASES = {  # name: (label, colour, marker)
    "2025-perturbed-small": ("perturbed, small", "#2a78d6", "o"),
    "2025-perturbed-large": ("perturbed, large", "#eb6834", "s"),
    "ebsd": ("EBSD", "#3a9a5b", "D"),
}
INK = "#2b2a27"
GRID = "#e4e2dc"
CACHE = os.path.join(HERE, "data", "alpha-sweep.json")
OUT = os.path.join(ROOT, "latex", "figures", "alpha-sweep.pdf")
COLLIDE = 1e-6  # min |p_i - p_j| / h below which converged generators count as merged

plt.rcParams.update({"font.size": 8, "axes.labelsize": 8, "xtick.labelsize": 7, "ytick.labelsize": 7})


def run_tool(cbar, V, alpha, tool):
    """introFigure on H_alpha from sites at the targets, with the per-iterate log.
    Returns (status, iterations, generators, weights, last min |p_i - p_j| / h)."""
    n = len(cbar)
    with tempfile.TemporaryDirectory() as tmp:
        src, dst = os.path.join(tmp, "in.txt"), os.path.join(tmp, "out.txt")
        with open(src, "w") as f:
            f.write("n %d\n" % n)
            f.write("targets " + " ".join("%.17g" % x for x in cbar.ravel()) + "\n")
            f.write("volumes " + " ".join("%.17g" % x for x in V) + "\n")
            f.write("sitesAtTargets 1\nalpha %.17g\nlog 1\n" % alpha)
        run = subprocess.run([tool, src, dst], capture_output=True, text=True)
        its = re.findall(r"^\s+it (\d+) .*min \|p_i - p_j\| / h (\S+)", run.stderr, re.M)
        dmin = float(its[-1][1]) if its else math.nan
        if run.returncode != 0:
            m = re.search(r"H_alpha solve: (.+?) \((\d+) iterations\)", run.stderr)
            status, it = (m.group(1), int(m.group(2))) if m else (run.stderr.strip().splitlines()[-1], 0)
            return status, it, None, None, dmin
        rec = {}
        with open(dst) as f:
            for line in f:
                key, *vals = line.split()
                rec[key] = np.array([float(v) for v in vals])
        it = int(re.search(r"H_alpha: (\d+) iterations", run.stdout).group(1))
    return "converged", it, rec["generators"].reshape(n, 2), rec["weights"], dmin


def residuals(c, cbar, V, h):
    e = np.linalg.norm(c - cbar, axis=1)
    return math.sqrt(np.sum(V * e ** 2) / np.sum(V)) / h, float(e.max()) / h


def sweep(args):
    cases = {c.name: c for c in BC.make_cases()}
    out = {}
    for name in CASES:
        case = cases[name]
        cbar, V = np.asarray(case.b), np.asarray(case.v)
        n, h = case.n, math.sqrt(V.sum() / case.n)
        # the volume-only diagram at the sites (alpha -> 0), as in the helpers
        with ref.Evaluation():
            _, c0, _ = ref.solve_OT_problem(cbar, V, case.om, relTol=BC.OT_MIN_REL_TOL)
        rms0, max0 = residuals(c0, cbar, V, h)
        entry = dict(n=n, distinct=BC.distinct_alpha(case), practical=BC.practical_alpha(case),
                     rms0=rms0, max0=max0, runs=[])
        for a in np.logspace(math.log10(args.amin), math.log10(args.amax), args.num):
            status, it, p, lam, dmin = run_tool(cbar, V, float(a), args.tool)
            r = dict(alpha=float(a), status=status, iterations=it, minGeneratorDistance=dmin)
            if p is not None:
                pd = PowerDiagram(positions=p, weights=lam, domain=case.om)
                c, A = pd.centroids(), pd.integrals()
                r["rms"], r["max"] = residuals(c, cbar, V, h)
                r["volumeError"] = float(np.max(np.abs(A - V) / V))
                if not dmin > COLLIDE:
                    r["status"] = "generators merged"
            out[name] = entry
            entry["runs"].append(r)
            print("%-22s alpha %9.3g  %-18s %3d it  rms %s  min|p_i-p_j|/h %.2g" % (
                name, a, r["status"], it, "%.3g" % r["rms"] if "rms" in r else "-", dmin), flush=True)
    os.makedirs(os.path.dirname(args.cache), exist_ok=True)
    json.dump(dict(amin=args.amin, amax=args.amax, num=args.num, cases=out), open(args.cache, "w"), indent=1)
    return out


def plot(data, path):
    fig, (top, bot) = plt.subplots(2, 1, figsize=(3.3, 3.4), sharex=True,
                                   gridspec_kw=dict(height_ratios=[2.0, 1.0]))
    fail_y = 0.96  # failed solves: one row per case above the data (axes fraction)
    for name, (label, color, marker) in CASES.items():
        e = data[name]
        ok = [r for r in e["runs"] if r["status"] == "converged"]
        bad = [r for r in e["runs"] if r["status"] != "converged"]
        a = np.array([r["alpha"] for r in ok])
        top.plot(a, [r["rms"] for r in ok], "-", color=color, lw=1.4, zorder=3)
        top.plot(a, [r["rms"] for r in ok], marker, color=color, ms=3.2, mec="white", mew=0.4, lw=0,
                 zorder=4, label=label)
        bot.plot(a, [r["iterations"] for r in ok], "-", color=color, lw=1.2, zorder=3)
        bot.plot(a, [r["iterations"] for r in ok], marker, color=color, ms=2.8, mec="white", mew=0.4, lw=0,
                 zorder=4)
        if bad:
            top.plot([r["alpha"] for r in bad], [fail_y] * len(bad), "x", color=color, ms=3.5, mew=1.0,
                     transform=top.get_xaxis_transform(), clip_on=False, zorder=5)
            fail_y -= 0.055
        for ax in (top, bot):
            ax.axvline(e["distinct"], color=color, lw=0.9, ls=(0, (4, 2)), zorder=0)
            ax.axvline(e["practical"], color=color, lw=0.9, ls=(0, (1, 1.5)), zorder=0)
    top.set_xscale("log")
    top.set_yscale("log")
    top.set_ylim(7e-4, 0.4)
    top.set_ylabel(r"RMS $|c_i - \bar c_i|\,/\,h$")
    bot.set_ylabel("iterations")
    bot.set_xlabel(r"stiffness $\alpha$")
    bot.set_ylim(bottom=0)
    top.set_xlim(data[next(iter(CASES))]["runs"][0]["alpha"] / 1.4, data[next(iter(CASES))]["runs"][-1]["alpha"] * 1.4)
    # legend: the cases, and the line styles of the bounds and of the failed solves
    from matplotlib.lines import Line2D
    kw = dict(fontsize=6, frameon=True, facecolor="white", edgecolor="none", framealpha=1.0, handlelength=1.8)
    top.legend(loc="lower left", **kw)
    h = [Line2D([], [], color="#8a8880", lw=0.9, ls=(0, (4, 2))), Line2D([], [], color="#8a8880", lw=0.9, ls=(0, (1, 1.5))),
         Line2D([], [], color="#8a8880", marker="x", lw=0, ms=3.5, mew=1.0)]
    bot.legend(h, ["guaranteed bound", "practical bound", "failed (top)"], loc="upper left", **kw)
    for ax in (top, bot):
        ax.grid(True, which="major", color=GRID, lw=0.6)
        ax.set_axisbelow(True)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
    fig.tight_layout(h_pad=0.6)
    fig.savefig(path)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--rerun", action="store_true", help="redo the solves even if the cache exists")
    ap.add_argument("--amin", type=float, default=1e-3)
    ap.add_argument("--amax", type=float, default=1e3)
    ap.add_argument("--num", type=int, default=31)
    ap.add_argument("--tool", default=os.path.join(ROOT, "build", "introFigure"))
    ap.add_argument("--cache", default=CACHE)
    ap.add_argument("--out", default=OUT)
    args = ap.parse_args()
    if args.rerun or not os.path.exists(args.cache):
        data = sweep(args)
    else:
        data = json.load(open(args.cache))["cases"]
    for name, e in data.items():
        bad = [r["alpha"] for r in e["runs"] if r["status"] != "converged"]
        ok = [r for r in e["runs"] if r["status"] == "converged"]
        print("%s (n = %d): guaranteed bound %.4g, practical bound %.4g; alpha -> 0 RMS %.3g; "
              "converged RMS %.3g to %.3g; largest converged alpha %.3g; smallest failed alpha %s" % (
                  name, e["n"], e["distinct"], e["practical"], e["rms0"], min(r["rms"] for r in ok),
                  max(r["rms"] for r in ok), max(r["alpha"] for r in ok), "%.3g" % min(bad) if bad else "-"))
    plot(data, args.out)
    print("wrote " + args.out)


if __name__ == "__main__":
    main()
