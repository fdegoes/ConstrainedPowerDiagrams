#!/usr/bin/env python3
"""Convergence curves (fig:convergence): the residual of rule 3 (latex/implementation.tex,
eq. results-termination) against the charged time, for our Newton solver and the three
methods of Bourne et al., on the two benchmarks of Table 1 (benchmarks/bourne/run_table1.py):

  (top)    voronoi, n = 20, the data of 2025-invert and 2026-voronoi;
  (bottom) cross, n = 728, the data of 2026-aligned.

Methods and settings are those of run_table1.py: ours (max H, rule 3 at 1e-8), SLSQP
(Bourne 2025), Frank-Wolfe at t = 100 and the subgradient method (Bourne 2026, t = 1250 on
voronoi, 100 on cross). run_table1.py runs their methods on rule 3 or to its time cap,
with their own rules disabled and read off the history; here each of their curves ends
at the iterate where its own rule (SLSQP rule 1, FW gap, change of the running average)
stops it, which is the run with their termination (the iterates do not depend on the rule).
A curve that reaches the time cap first ends there; the stops are not marked.

The residual drawn is the centroid residual eps_u of the equation, max_i |c_i - cbar_i| /
Vbar_i^{1/d} (bourne_cases.common_residuals); --y max draws max(eps_w, eps_u), with the
volume residual eps_w = max_i |V_i - Vbar_i| / Vbar_i. report() prints the largest ratio
eps_w / eps_u over our iterates. Their iterates are diagrams solved for the
prescribed volumes (eps_w ~ 1e-10 at every iterate), ours move shifts and weights together
(eps_w falls below the tolerance only at convergence). Every curve starts at t = 0 with
iteration 0, the start generators on the diagram with the prescribed volumes (uncharged for
their methods; for ours, the initial volume solve, whose time is charged to iteration 1). The
time axis is linear on [0, 1 ms] and logarithmic above.

  MPLCONFIGDIR=/tmp/mplcfg external/bourne/.venv/bin/python benchmarks/figures/plot_convergence.py \
      --from benchmarks/bourne/results/table1-<date> [--y max|centroid] [--out latex/figures/convergence.pdf]

--ebsd prints (no figure) the former EBSD comparison from the cache in
benchmarks/figures/cache/convergence/ (--rerun to recompute): our H_alpha at 0.99 x the
practical estimate, and their FW at t = 10 and at t = 1 + alpha.
"""

import argparse
import json
import math
import os
import subprocess
import sys
import tempfile
import time

os.environ.setdefault("MPLCONFIGDIR", "/tmp/mplcfg")
import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(ROOT, "benchmarks", "bourne"))
import bourne_cases as bc  # noqa: E402
import bourne_ref as ref  # noqa: E402

CACHE = os.path.join(HERE, "cache", "convergence")
OUT = os.path.join(ROOT, "latex", "figures", "convergence.pdf")
TOL = 1e-8
LINTHRESH = 1e-3  # s: the time axis is linear below
EVAL_TOL = bc.OT_MIN_REL_TOL  # OT tolerance of the evaluations (relative to min V)

INK = "#2b2a27"
GRID = "#e4e2dc"
STYLE = {  # label, colour, line style, marker
    "ours": ("ours (Newton)", "#2a78d6", "-", "o"),
    "slsqp": ("SLSQP [2025]", "#eb6834", "-", None),
    "fw": ("FW [2026]", "#1a9e77", "-", None),
    "fw-same": ("FW [2026], our $\\alpha$", "#1a9e77", "--", None),
    "subgrad": ("subgradient [2026]", "#7570b3", "-", None),
}

plt.rcParams.update({"font.size": 8, "axes.labelsize": 8, "xtick.labelsize": 7, "ytick.labelsize": 7})


def _cases():
    return {c.name: c for c in bc.make_cases()}


def _evaluate(case, X, alpha):
    """(volume, shift) residuals of the power diagram of generators X with the
    prescribed volumes; X are FW's scaled generators when alpha is finite."""
    with ref.Evaluation():
        _, c, V = ref.solve_OT_problem(np.asarray(X), case.v, case.om, relTol=EVAL_TOL, strict=False)
    p = alpha * np.asarray(X) if math.isfinite(alpha) else np.asarray(X)
    return bc.common_residuals(case, V, c, p, alpha)


# --------------------------------------------------------------------------
# Runs
# --------------------------------------------------------------------------

def run_ours(case, solver, alpha):
    """benchBourne, sites at the targets, rule 3 at TOL; the residuals of every
    iterate, recomputed from the trajectory as in bourne_cases.run_mine."""
    lines = ["box %r %r %r %r" % tuple(float(x) for x in case.box), "n %d" % case.n, "solver %s" % solver,
             "alpha %r" % (float(alpha) if math.isfinite(alpha) else 1e300), "start 0",
             "volumetol %r" % TOL, "shifttol %r" % TOL, "maxiter 100", "lighttrajectory",
             "volumes " + " ".join(repr(float(x)) for x in case.v),
             "targets " + " ".join(repr(float(x)) for x in case.b.ravel()),
             "sites " + " ".join(repr(float(x)) for x in np.asarray(case.sites).ravel())]
    with tempfile.TemporaryDirectory() as d:
        cpath, opath = os.path.join(d, "case.txt"), os.path.join(d, "out.json")
        with open(cpath, "w") as f:
            f.write("\n".join(lines) + "\n")
        subprocess.run([bc.BENCH, cpath, opath], check=True)
        with open(opath) as f:
            r = json.load(f)
    traj = bc.trim_to_solution(r["trajectory"], r["generators"])
    hist = []
    for k, it in enumerate(traj):
        V, c = np.asarray(it["volumes"]), np.reshape(it["centroids"], (-1, 2))
        p = np.reshape(it["generators"], (-1, 2))
        vol, shift = bc.common_residuals(case, V, c, p, alpha)
        hist.append(dict(iteration=k, solves=it["numBuilds"], volume=vol, shift=shift))
    return dict(method="ours", solver=solver, alpha=alpha if math.isfinite(alpha) else None, status=r["status"],
                iterations=r["iterations"], attempts=r["attempts"], builds=r["numBuilds"], history=hist)


def run_slsqp(case):
    """Their SLSQP on H (Bourne 2025) with their settings, stopped on their rule."""
    ref.COUNTERS.reset()
    ref.OT_REL_TOL["value"] = 1e-3
    hist = []

    def record(X):
        vol, shift = _evaluate(case, X, math.inf)
        hist.append(dict(iteration=len(hist), solves=ref.COUNTERS.otSolves, volume=vol, shift=shift))
        return False

    record(case.start)
    t0 = time.perf_counter()
    res, _ = ref.max_H_solver(case.start, case.b, case.v, case.om, R0=case.R0, d=case.d, max_iter=case.maxiter,
                              ptol=case.ptol, callback=record)
    return dict(method="slsqp", status=str(res.message), success=bool(res.success), iterations=int(res.nit),
                otSolves=ref.COUNTERS.otSolves, seconds=time.perf_counter() - t0 - ref.COUNTERS.evalSeconds,
                history=hist)


def run_fw(case, t, alpha_ref, niter):
    """Their FW (Bourne 2026) at dilation t from case.start, run to niter
    iterations with their gap rule disabled; the gap of every iterate is
    recorded, so that the stop of their rule (gap < case.gap) is read off the
    history. Residual of H (alpha_ref = inf) or H_alpha."""
    ref.COUNTERS.reset()
    ref.OT_REL_TOL["value"] = 1e-4
    y = t / (t - 1) * case.b
    hist = []

    def cb(it, B, gap):
        vol, shift = _evaluate(case, y - B, alpha_ref)
        # The LMO of this iterate has been solved when the callback is called.
        hist.append(dict(iteration=int(it), solves=ref.COUNTERS.otSolves, volume=vol, shift=shift, gap=float(gap)))
        return False

    _, k, status = ref.convex_order_projection_frank_wolfe(case.om, y, case.start, R0=case.R0, masses=case.v,
                                                           niter=niter, tol=-1.0, callback=cb)
    return dict(method="fw", t=t, referenceAlpha=alpha_ref if math.isfinite(alpha_ref) else None, gapTol=case.gap,
                status=status, iterations=int(k), otSolves=ref.COUNTERS.otSolves, history=hist)


def run_subgrad(case, niter=600, tau=2.0, tol=1e-5):
    """Their subgradient method (Bourne 2026, convex_order_projection) at dilation t = 100 from case.start, run to
    niter iterations with their rule disabled; the change of the running average (their rule: < tol) is recorded
    as the gap of every iterate, so that their stop is read off the history. Residual of H."""
    ref.COUNTERS.reset()
    ref.OT_REL_TOL["value"] = 1e-4
    t = 100.0
    y = t / (t - 1) * case.b
    hist = []

    def cb(i, Z, error):
        # Generators of the iterate's diagram: the rescaled iterate, as in their code (run_theirs.subgrad_sites).
        X = (Z - np.min(Z) + .02) / (np.max(Z) - np.min(Z) + .02)
        vol, shift = _evaluate(case, X, math.inf)
        hist.append(dict(iteration=i + 1, solves=ref.COUNTERS.otSolves, volume=vol, shift=shift, gap=float(error)))
        return False

    _, k, status, _ = ref.convex_order_subgradient(case.om, case.start.copy(), y, case.v, niter=niter, tau=tau,
                                                   tol=-1.0, callback=cb)
    return dict(method="subgrad", t=t, referenceAlpha=None, gapTol=tol, status=status, iterations=int(k),
                otSolves=ref.COUNTERS.otSolves, history=hist)


def _cached(name, rerun, fn):
    path = os.path.join(CACHE, name + ".json")
    if not rerun and os.path.exists(path):
        with open(path) as f:
            return json.load(f)
    t0 = time.time()
    r = fn()
    r["wallSeconds"] = time.time() - t0
    os.makedirs(CACHE, exist_ok=True)
    with open(path, "w") as f:
        json.dump(r, f, indent=1)
    print("%-22s %-8s %-28s %4d iterations, final residual %.2e (%.0f s)" % (
        name, r["method"], r["status"][:28], len(r["history"]) - 1,
        max(r["history"][-1]["volume"], r["history"][-1]["shift"]), r["wallSeconds"]), flush=True)
    return r


def load_table1(d, benchmark):
    """The runs of benchmarks/bourne/run_table1.py in directory d: per method the history (charged seconds,
    volume and centroid residuals per iterate), the stop, and ownRuleIndex (the iterate where the method's own
    rule would stop it; None for ours, or when it never does)."""
    out = {}
    for m in ("ours", "slsqp", "fw", "subgrad"):
        path = os.path.join(d, "%s__%s.json" % (benchmark, m))
        if not os.path.exists(path):
            print("missing %s: %s not drawn" % (path, m))
            continue
        with open(path) as f:
            r = json.load(f)
        h = r["history"]
        own = r.get("ownRule") or {}
        out[m] = dict(method=m, stop=r["stop"], status="%s (%s)" % (r["stop"], r["status"]),
                      iterations=r.get("iterations"), ownRuleIndex=own.get("index"),
                      history=[dict(iteration=k, solves=n, seconds=t, volume=v, shift=c, gap=g)
                               for k, n, t, v, c, g in zip(h["iteration"], h["solves"], h["seconds"], h["volume"],
                                                           h["centroid"], h["own"])])
    return out


def ebsd_histories(args):
    """The former EBSD comparison (cached): our H_alpha and their FW at t = 10 and at our alpha."""
    cases = _cases()
    ebsd = cases["ebsd"]
    ebsd26 = [v for v in ebsd.variants if v.paper == "2026"][0]
    a = 0.99 * bc.practical_alpha(ebsd)
    return {"ebsd": {
        "ours": _cached("ebsd-ours", args.rerun, lambda: run_ours(ebsd, "Halpha", a)),
        "fw": _cached("ebsd-fw", args.rerun, lambda: run_fw(ebsd26, ebsd26.t, ebsd26.t - 1, args.fw_iter)),
        "fw-same": _cached("ebsd-fw-same", args.rerun, lambda: run_fw(ebsd26, 1 + a, a, args.fw_iter)),
    }}, a


# --------------------------------------------------------------------------
# Figure
# --------------------------------------------------------------------------

def residual(e, y):
    return e["shift"] if y == "centroid" else max(e["volume"], e["shift"])


def own_stop(r):
    """Index of the last drawn iterate: where the method's own rule stops it (ours: rule 3, its last iterate)."""
    if r["method"] == "ours" or r.get("ownRuleIndex") is None:
        return len(r["history"]) - 1
    return r["ownRuleIndex"]


def plot(runs, args):
    # One column, the two benchmarks stacked, sharing the time axis; one legend above both.
    fig, axes = plt.subplots(2, 1, figsize=(3.3, 4.0), sharex=True)
    panels = [("voronoi", "Voronoi, $n = 20$"), ("cross", "cross, $n = 728$")]
    tmax = 0.0
    for ax, (key, title) in zip(axes, panels):
        for m, r in runs[key].items():
            label, color, ls, marker = STYLE[m]
            k = own_stop(r)
            h = r["history"][:k + 1]
            if not h:
                continue
            xs = np.array([0.0 if e["iteration"] == 0 else e["seconds"] for e in h])   # every start at t = 0
            ys = np.maximum([residual(e, args.y) for e in h], 1e-16)
            tmax = max(tmax, xs[-1])
            ax.plot(xs, ys, "-", color=color, lw=1.3 if m == "ours" else 1.0, zorder=3 if m == "ours" else 2)
            if m == "ours":
                ax.plot(xs, ys, "o", color=color, ms=3.0, mec="white", mew=0.4, lw=0, zorder=4)
        ax.axhline(TOL, color="#8a8880", lw=0.9, ls=(0, (1, 1.5)), zorder=1)
        ax.annotate("tolerance $10^{-8}$", (0.98, TOL), xycoords=("axes fraction", "data"), xytext=(0, 3),
                    textcoords="offset points", ha="right", va="bottom", fontsize=7, color="#8a8880")
        ax.text(0.03, 0.04, title, transform=ax.transAxes, ha="left", va="bottom", fontsize=8, color=INK)
        ax.set_xscale("symlog", linthresh=LINTHRESH, linscale=0.4)   # linear on [0, 1 ms], log above
        ax.set_yscale("log")
        ax.set_ylim(1e-12, 1e2)
        ax.set_ylabel("$\\varepsilon_{\\mathbf{u}}$" if args.y == "centroid" else
                      "$\\max(\\varepsilon_w, \\varepsilon_{\\mathbf{u}})$")
        ax.grid(True, which="major", color=GRID, lw=0.6)
        ax.set_axisbelow(True)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
    axes[1].set_xlim(0, 1.6 * tmax)
    axes[1].set_xlabel("time (s)")
    from matplotlib.lines import Line2D
    handles = [Line2D([], [], color=STYLE[m][1], lw=1.3, marker="o" if m == "ours" else None, ms=3.0, mec="white",
                      label=STYLE[m][0]) for m in ("ours", "slsqp", "fw", "subgrad")]
    fig.tight_layout(h_pad=0.6, rect=(0, 0, 1, 0.9))
    # the legend centred over the plot boxes (not the figure, whose left part holds the y labels),
    # and the saved figure cropped to its content, so that the paper centres the plots
    box = axes[0].get_position()
    fig.legend(handles=handles, frameon=False, fontsize=7, loc="upper center", ncol=2, handlelength=1.4,
               columnspacing=1.0, handletextpad=0.4, bbox_to_anchor=(0.5 * (box.x0 + box.x1), 1.0))
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    fig.savefig(args.out, bbox_inches="tight", pad_inches=0.02)
    if args.png:
        fig.savefig(os.path.splitext(args.out)[0] + ".png", dpi=200, bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)


def report(runs, y="centroid"):
    for key, rs in runs.items():
        if "ours" in rs:
            h = rs["ours"]["history"]
            print("%-8s ours: largest eps_w / eps_u over its iterates %.2f" % (
                key, max(e["volume"] / max(e["shift"], 1e-300) for e in h)))
        for m, r in rs.items():
            h = r["history"]
            res = [residual(e, y) for e in h]
            k = own_stop(r) if "seconds" in h[0] else len(h) - 1
            first = next((e for e, x in zip(h, res) if x <= TOL), None)
            print("%-8s %-8s %-30s drawn to iteration %d (%s s): residual %.2e (vol %.1e, centroid %.1e); "
                  "rule 3 first met: %s" % (
                      key, m, r["status"][:30], h[k]["iteration"], "%.3g" % h[k]["seconds"] if "seconds" in h[k]
                      else "-", res[k], h[k]["volume"], h[k]["shift"],
                      "-" if first is None else "iteration %d" % first["iteration"]))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--from", dest="from_dir", default=None, help="output directory of benchmarks/bourne/run_table1.py")
    ap.add_argument("--y", choices=["max", "centroid"], default="centroid",
                    help="max(eps_w, eps_u) of rule 3, or the centroid residual eps_u alone")
    ap.add_argument("--out", default=OUT)
    ap.add_argument("--png", action="store_true", help="also write a PNG next to the PDF")
    ap.add_argument("--ebsd", action="store_true", help="print the former EBSD comparison (cache), no figure")
    ap.add_argument("--rerun", action="store_true", help="with --ebsd: recompute the cached histories")
    ap.add_argument("--fw-iter", type=int, default=500, help="with --ebsd: iteration cap of FW")
    args = ap.parse_args()
    if args.ebsd:
        runs, alpha = ebsd_histories(args)
        print("alpha (EBSD, 0.99 x practical) = %.4g" % alpha)
        report(runs)
        return
    if not args.from_dir:
        ap.error("--from DIR (the output of benchmarks/bourne/run_table1.py) is required")
    runs = {b: load_table1(args.from_dir, b) for b in ("voronoi", "cross")}
    report(runs, args.y)
    plot(runs, args)
    print("wrote", args.out)


if __name__ == "__main__":
    main()
