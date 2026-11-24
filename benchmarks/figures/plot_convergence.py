#!/usr/bin/env python3
"""Convergence curves: the residual of rule 3 (latex/implementation.tex,
eq. results-termination) against the iteration, for our Newton solver and the
methods of Bourne et al., on two cases of benchmarks/bourne/bourne_cases.py.

  (a) Voronoi, n = 20 (realizable): the data of 2025-invert and 2026-voronoi.
      our H            benchBourne, sites at the targets, tolerance 1e-8;
      SLSQP (2025)     their max_H_solver from their random start (2025-invert:
                       inner OT 1e-3 min V, ptol 1e-8, maxiter 1000), stopped on
                       their own rule;
      FW (2026)        t = 100 from the true seeds (2026-voronoi), capped at
                       --fw-iter iterations; a cross marks where their gap rule
                       (1e-14) stops;
      subgradient [2026] t = 100 from the same start, their step tau/(i+1)
                       with tau = 2, capped at 600 iterations; a cross marks
                       where their rule (change of the running average
                       below 1e-5) stops.
      Every curve is the residual of H.
  (b) EBSD map, n = 243 (approximate), unit square.
      our H_alpha      alpha = 0.99 x the practical estimate
                       (bourne_cases.practical_alpha), tolerance 1e-8;
      FW (2026)        t = 10 from 0.5 (their settings), capped at --fw-iter
                       iterations: residual of H_alpha, alpha = t - 1 = 9;
      FW at our alpha  t = 1 + alpha, same start: residual of the same
                       H_alpha as ours.
      FW runs to the cap with their gap rule disabled; a cross marks the
      iterate where their rule (gap < 1e-9) stops.
      Their 2025 least-squares SLSQP solves neither H nor H_alpha, so it has
      no rule-3 residual and is not drawn.

The residual of an iterate is max(volume, shift) of bourne_cases.common_residuals:
max_i |V_i - Vbar_i| / Vbar_i and max_i |g_u,i| / (2 Vbar_i^{3/2}); rule 3 holds
when it is <= 1e-8. For their methods it is evaluated at every iterate on the
power diagram of the iterate's generators with the prescribed volumes (a tight
OT solve, not charged); for ours it is that of the solver's own iterate (the
trajectory of benchBourne). Iteration 0 is the start.

The histories are cached as JSON in benchmarks/figures/cache/convergence/
(--rerun to recompute them). Needs build/benchBourne and the pysdot Python:
  MPLCONFIGDIR=/tmp/mplcfg external/bourne/.venv/bin/python benchmarks/figures/plot_convergence.py
      [--rerun] [--fw-iter 500] [--x iterations|solves] [--out latex/figures/convergence.pdf]
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
EVAL_TOL = bc.OT_MIN_REL_TOL  # OT tolerance of the evaluations (relative to min V)

INK = "#2b2a27"
GRID = "#e4e2dc"
STYLE = {  # label, colour, line style, marker
    "ours": ("ours (Newton)", "#2a78d6", "-", "o"),
    "slsqp": ("SLSQP [2025]", "#eb6834", "-", None),
    "fw": ("FW [2026]", "#1a9e77", "-", None),
    "fw-same": ("FW [2026], our $\\alpha$", "#1a9e77", "--", None),
    "subgrad": ("subgradient [2026]", "#7570b3", "-.", None),
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


def histories(args):
    cases = _cases()
    inv, vor, ebsd = cases["2025-invert"], cases["2026-voronoi"], cases["ebsd"]
    ebsd26 = [v for v in ebsd.variants if v.paper == "2026"][0]
    a = 0.99 * bc.practical_alpha(ebsd)
    out = {
        "voronoi": {
            "ours": _cached("voronoi-ours", args.rerun, lambda: run_ours(vor, "H", math.inf)),
            "slsqp": _cached("voronoi-slsqp", args.rerun, lambda: run_slsqp(inv)),
            "fw": _cached("voronoi-fw", args.rerun, lambda: run_fw(vor, vor.t, math.inf, args.fw_iter)),
            "subgrad": _cached("voronoi-subgrad", args.rerun, lambda: run_subgrad(vor)),
        },
        "ebsd": {
            "ours": _cached("ebsd-ours", args.rerun, lambda: run_ours(ebsd, "Halpha", a)),
            "fw": _cached("ebsd-fw", args.rerun, lambda: run_fw(ebsd26, ebsd26.t, ebsd26.t - 1, args.fw_iter)),
            "fw-same": _cached("ebsd-fw-same", args.rerun, lambda: run_fw(ebsd26, 1 + a, a, args.fw_iter)),
        },
    }
    return out, a


# --------------------------------------------------------------------------
# Figure
# --------------------------------------------------------------------------

def _series(r, x):
    h = r["history"]
    xs = np.array([e["iteration" if x == "iterations" else "solves"] for e in h], dtype=float)
    ys = np.array([max(e["volume"], e["shift"]) for e in h])
    return xs, np.maximum(ys, 1e-16)


def gap_stop(r):
    """Index of the iterate where their gap rule stops FW (None: not reached)."""
    if r["method"] == "slsqp":
        return len(r["history"]) - 1  # it stops on its own rule (rule 1)
    if r["method"] not in ("fw", "subgrad"):
        return None
    return next((k for k, e in enumerate(r["history"]) if e["gap"] < r["gapTol"]), None)


def plot(runs, alpha, args):
    # One column, the two cases stacked. On EBSD only FW at our alpha is drawn: with their
    # t = 10 (alpha = 9, far above the bound) it does not converge, which the text states.
    fig, axes = plt.subplots(2, 1, figsize=(3.3, 4.0))
    panels = [("voronoi", "(a) Voronoi, $n = 20$, realizable: $\\mathcal{G}$"),
              ("ebsd", "(b) EBSD, $n = 243$, approximate: $\\mathcal{G}_\\alpha$")]
    for ax, (key, title) in zip(axes, panels):
        for m, r in runs[key].items():
            if key == "ebsd" and m == "fw":
                continue
            label, color, ls, marker = STYLE[m]
            if m == "fw-same":
                label, ls = "FW [2026], same $\\alpha$", "-"
            xs, ys = _series(r, args.x)
            ax.plot(xs, ys, ls, color=color, lw=1.3 if m == "ours" else 1.0, label=label,
                    zorder=3 if m == "ours" else 2)
            if marker:
                ax.plot(xs, ys, marker, color=color, ms=3.2, mec="white", mew=0.4, lw=0, zorder=4)
            k = gap_stop(r)
            if k is not None:
                ax.plot([xs[k]], [ys[k]], "x", color=INK, ms=4.5, mew=1.0, zorder=5)
        ax.axhline(TOL, color="#8a8880", lw=0.9, ls=(0, (1, 1.5)), zorder=1)
        ax.annotate("tolerance $10^{-8}$", (0.62, TOL), xycoords=("axes fraction", "data"), xytext=(-2, 3),
                    textcoords="offset points", ha="right", va="bottom", fontsize=7, color="#8a8880")
        ax.set_xscale("symlog", linthresh=10, linscale=0.6)
        ax.set_yscale("log")
        ax.set_ylim(1e-15, 1e2)
        ax.set_xlim(0, 1.15 * max(len(r["history"]) for r in runs[key].values()))
        ax.set_xlabel("iteration" if args.x == "iterations" else "OT solves / power diagrams")
        ax.set_title(title, fontsize=8, color=INK, loc="left")
        ax.grid(True, which="major", color=GRID, lw=0.6)
        ax.set_axisbelow(True)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        ax.plot([], [], "x", color=INK, ms=4.5, mew=1.0, label="their stopping rule")
        ax.legend(frameon=False, fontsize=7, loc="lower right", handlelength=1.8, borderaxespad=0.2)
        ax.set_ylabel("residual (rule 3)")
    fig.tight_layout(h_pad=0.8)
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    fig.savefig(args.out)
    fig.savefig(os.path.splitext(args.out)[0] + ".png", dpi=200) if args.png else None
    plt.close(fig)


def report(runs, alpha):
    print("alpha (EBSD, 0.99 x practical) = %.4g" % alpha)
    for key, rs in runs.items():
        for m, r in rs.items():
            h = r["history"]
            res = [max(e["volume"], e["shift"]) for e in h]
            first = next((e for e, x in zip(h, res) if x <= TOL), None)
            k = gap_stop(r)
            if k is not None:
                print("%-8s %-8s their gap rule stops at iteration %d, residual %.2e" % (key, m, k, res[k]))
            print("%-8s %-8s %-26s %4d it, %5d solves, final %.2e (vol %.1e, shift %.1e), min %.2e, to 1e-8: %s" % (
                key, m, r["status"][:26], h[-1]["iteration"], h[-1]["solves"], res[-1], h[-1]["volume"],
                h[-1]["shift"], min(res), "-" if first is None else "%d it / %d solves" % (
                    first["iteration"], first["solves"])))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--rerun", action="store_true", help="recompute the cached histories")
    ap.add_argument("--fw-iter", type=int, default=500, help="iteration cap of FW")
    ap.add_argument("--x", choices=["iterations", "solves"], default="iterations")
    ap.add_argument("--out", default=OUT)
    ap.add_argument("--png", action="store_true", help="also write a PNG next to the PDF")
    args = ap.parse_args()
    runs, alpha = histories(args)
    report(runs, alpha)
    plot(runs, alpha, args)
    print("wrote", args.out)


if __name__ == "__main__":
    main()
