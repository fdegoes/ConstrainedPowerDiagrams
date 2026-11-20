#!/usr/bin/env python3
"""Head-to-head comparison on the approximate cases of Bourne et al.
(2025-perturbed-small, 2025-perturbed-large, ebsd): their methods and our
max H_alpha, each run until it terminates, with up to three stiffnesses per
case:

  bound      our bound (latex/contributions.tex, eq:regularized-distinct),
             bourne_cases.distinct_alpha():
               alpha < min_{i != j} |q_i - q_j| / sqrt(resid (1/Vbar_i + 1/Vbar_j)),
               resid = sum_i Vbar_i |c0_i - cbar_i|^2, c0 the centroids of max F
               at the sites q = cbar
  practical  our per-cell estimate (eq:regularized-practical),
             bourne_cases.practical_alpha():
               alpha < min_{i != j} |q_i - q_j| / (|c0_i - cbar_i| + |c0_j - cbar_j|)
             (both --bound-factor (0.99) times these values)
  theirs     their hand-picked value (Bourne 2026 sec. 6.3: "we set t = 10" for
             the EBSD map, alpha = t - 1 = 9); none for the perturbed cases,
             whose method (SLSQP on H, Bourne 2025) has no stiffness

Runs per case (FW with dilation t = 1 + alpha solves the same projection as
our H_alpha, so FW and ours run at the same alphas):
  perturbed-*  their SLSQP on H (no alpha), and ours at bound and practical
  ebsd         their FW and ours at bound, practical and theirs (t = 1 + alpha
               for FW), and their 2025 least-squares SLSQP (no alpha, in their
               251.25 um box)

Every run stops when its method converges by its own rule (SLSQP: their
relative objective change ptol; FW: their gap; ours: volume and shift
residuals <= --tol), or at --max-iter iterations, or after --max-seconds
seconds of wall time, whichever comes first. For each run it reports why it
stopped, the iterations, the time, the OT solves (theirs) or diagram builds
(ours), and the residuals of the returned solution, all evaluated in the
same way (benchmarks/bourne/bourne_cases.py):
  volume      max_i |V_i - Vbar_i| / Vbar_i
  centroid    max_i |c_i - cbar_i| / sqrt(|domain|)
  stationary  the per-cell shift gradient of the problem the run solves
              (H for SLSQP on H, H_alpha for FW and ours; none for the
              least-squares SLSQP, which solves neither), as in the common rule

Outputs in --out: one JSON per run, results.csv and summary.txt.

Long: the EBSD least-squares SLSQP (finite-difference gradients, 243 cells)
typically runs to the time cap. Run overnight, e.g.
  caffeinate -i nohup external/bourne/.venv/bin/python benchmarks/bourne/run_comparison.py \\
      > benchmarks/bourne/results/comparison.log 2>&1 &

Usage:
  external/bourne/.venv/bin/python benchmarks/bourne/run_comparison.py
      [--case NAME ...] [--max-iter 1000] [--max-iter-mine 100] [--max-seconds 1800] [--tol 1e-8]
      [--bound-factor 0.99] [--out benchmarks/bourne/results/comparison]
      [--list]
"""

import argparse
import csv
import json
import math
import os
import subprocess
import sys
import tempfile
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import bourne_cases as bc  # noqa: E402
import bourne_ref as ref  # noqa: E402

# Their hand-picked stiffness: Bourne 2026 sec. 6.3.2 ("we set t = 10").
THEIR_T = {"ebsd": 10.0}


def _solve_volumes(case, X):
    """Centroids and volumes of the power diagram of generators X with the
    prescribed volumes: the tightest OT solve that pysdot reaches (about
    1e-9 min(v); on some EBSD diagrams only about 1e-8), looser if needed."""
    for tol in (bc.OT_MIN_REL_TOL, 1e-8, 1e-7, 1e-6):
        try:
            with ref.Evaluation():
                _, c, V = ref.solve_OT_problem(np.asarray(X), case.v, case.om, relTol=tol)
            return c, V
        except RuntimeError:
            continue
    raise RuntimeError("the volumes of the final diagram cannot be solved to 1e-6 min(v)")


def residuals(case, X, alpha):
    """Volume, centroid and stationarity residuals of the power diagram with
    generators X and the prescribed volumes (stationarity of H_alpha, or H for
    alpha = inf, or None)."""
    try:
        c, V = _solve_volumes(case, X)
    except RuntimeError as e:
        return dict(residualError=str(e))
    out = dict(volume=float(np.max(np.abs(V - case.v) / case.v)),
               centroid=float(np.max(np.linalg.norm(c - case.b, axis=1)) / math.sqrt(case.om.measure())))
    if alpha is not None:
        p = alpha * np.asarray(X) if math.isfinite(alpha) else np.asarray(X)
        out["stationary"] = bc.common_residuals(case, V, c, p, alpha)[1]
    return out


def run_slsqp(case, args):
    """Their 2025 SLSQP (on H, or least squares), with their rule, capped."""
    ref.COUNTERS.reset()
    ref.OT_REL_TOL["value"] = 1e-3
    t0 = time.perf_counter()
    k = [0]
    timed_out = [False]

    def cb(X):
        k[0] += 1
        if time.perf_counter() - t0 - ref.COUNTERS.evalSeconds > args.max_seconds:
            timed_out[0] = True
            return True
        return False

    try:
        if case.objective == "H":
            res, X = ref.max_H_solver(case.start, case.b, case.v, case.om, R0=case.R0, d=case.d,
                                      max_iter=args.max_iter, ptol=case.ptol, callback=cb)
        else:
            res, X = ref.min_F_solver(case.start, case.b, case.v, case.om, d=case.d, max_iter=args.max_iter,
                                      ptol=case.ptol, callback=cb)
        stop = "time cap" if timed_out[0] else ("converged" if res.success else (
            "iteration cap" if "Iteration limit" in str(res.message) else str(res.message)))
        iters = int(res.nit)
    except Exception as e:  # noqa: BLE001
        return dict(stop="error: %s" % e, iterations=k[0], seconds=time.perf_counter() - t0 - ref.COUNTERS.evalSeconds,
                    otSolves=ref.COUNTERS.otSolves)
    seconds = time.perf_counter() - t0 - ref.COUNTERS.evalSeconds
    out = dict(stop=stop, iterations=iters, seconds=seconds, otSolves=ref.COUNTERS.otSolves)
    out.update(residuals(case, X, math.inf if case.objective == "H" else None))
    return out


def slsqp_rule(case):
    """Their 2025 rule: SLSQP's relative objective change, ptol times the
    objective at the start."""
    return "SLSQP objective change (Bourne 2025, %s)" % ("H" if case.objective == "H" else "least squares"), case.ptol


def run_fw(case, t, args):
    """Their 2026 FW at dilation t (= H_alpha with alpha = t - 1), with their
    gap rule, capped."""
    ref.COUNTERS.reset()
    ref.OT_REL_TOL["value"] = 1e-4
    y = t / (t - 1) * case.b
    t0 = time.perf_counter()
    timed_out = [False]

    def cb(it, B, gap):
        if time.perf_counter() - t0 - ref.COUNTERS.evalSeconds > args.max_seconds:
            timed_out[0] = True
            return True
        return False

    try:
        B, k, status = ref.convex_order_projection_frank_wolfe(case.om, y, case.start, R0=case.R0, masses=case.v,
                                                               niter=args.max_iter, tol=case.gap, callback=cb)
    except Exception as e:  # noqa: BLE001
        return dict(stop="error: %s" % e, seconds=time.perf_counter() - t0 - ref.COUNTERS.evalSeconds,
                    otSolves=ref.COUNTERS.otSolves)
    seconds = time.perf_counter() - t0 - ref.COUNTERS.evalSeconds
    stop = "time cap" if timed_out[0] else ("converged" if status == "converged" else "iteration cap")
    out = dict(stop=stop, iterations=int(k), seconds=seconds, otSolves=ref.COUNTERS.otSolves)
    # FW's generators Y = y - B correspond to our p / alpha (eq. 1 + 1/alpha dilation).
    out.update(residuals(case, y - B, t - 1))
    return out


OUR_RULE = "volume and shift residuals of H_alpha (eq. results-termination)"
FW_RULE = "FW gap (Bourne 2026)"


def run_ours(case, alpha, args):
    """Our max H_alpha (benchBourne) with sites and targets cbar, stopped on
    volume and shift residuals <= --tol, capped."""
    lines = ["box %r %r %r %r" % tuple(float(x) for x in case.box), "n %d" % case.n, "solver Halpha",
             "alpha %r" % float(alpha), "start 0", "volumetol %r" % args.tol, "shifttol %r" % args.tol,
             "maxiter %d" % args.max_iter_mine, "notrajectory",
             "volumes " + " ".join(repr(float(x)) for x in case.v),
             "targets " + " ".join(repr(float(x)) for x in case.b.ravel()),
             "sites " + " ".join(repr(float(x)) for x in np.asarray(case.sites).ravel())]
    with tempfile.TemporaryDirectory() as d:
        cpath, opath = os.path.join(d, "case.txt"), os.path.join(d, "out.json")
        with open(cpath, "w") as f:
            f.write("\n".join(lines) + "\n")
        try:
            subprocess.run([bc.BENCH, cpath, opath], check=True, timeout=args.max_seconds + 60)
        except subprocess.TimeoutExpired:
            return dict(stop="time cap", seconds=float(args.max_seconds))
        except subprocess.CalledProcessError as e:
            # A crash of benchBourne on one run: record it and go on.
            return dict(stop="error (exit %d)" % e.returncode)
        r = json.load(open(opath))
    stop = "converged" if r["status"] == "converged" else (
        "iteration cap" if r["status"] == "max iterations" else r["status"])
    out = dict(stop=stop, iterations=r["iterations"], seconds=r["seconds"], builds=r["numBuilds"])
    V, c = np.asarray(r["volumes"]), np.reshape(r["centroids"], (-1, 2))
    if V.size:
        p = np.reshape(r["generators"], (-1, 2))
        out.update(volume=float(np.max(np.abs(V - case.v) / case.v)),
                   centroid=float(np.max(np.linalg.norm(c - case.b, axis=1)) / math.sqrt(case.om.measure())),
                   stationary=bc.common_residuals(case, V, c, p, alpha)[1])
    return out


FIELDS = ["case", "method", "alphaChoice", "alpha", "rule", "tolerance", "stop", "iterations", "seconds", "otSolves",
          "builds", "volume", "centroid", "stationary"]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    names = [c.name for c in bc.make_cases() if c.group == bc.APPROXIMATE]
    ap.add_argument("--case", nargs="+", choices=names, default=names)
    ap.add_argument("--max-iter", type=int, default=1000, help="iteration cap of their methods")
    ap.add_argument("--max-iter-mine", type=int, default=100, help="Newton iterations per attempt of ours")
    ap.add_argument("--max-seconds", type=float, default=1800.0)
    ap.add_argument("--tol", type=float, default=1e-8, help="our volume and shift tolerance")
    ap.add_argument("--bound-factor", type=float, default=0.99,
                    help="alpha = factor * our bound and factor * our practical estimate")
    ap.add_argument("--out", default=os.path.join(HERE, "results", "comparison"))
    ap.add_argument("--list", action="store_true", help="print the runs and alphas, without running")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    plan = []
    for case in (c for c in bc.make_cases() if c.group == bc.APPROXIMATE and c.name in args.case):
        variants = getattr(case, "variants", [case])
        alphas = [("bound", args.bound_factor * bc.distinct_alpha(case)),
                  ("practical", args.bound_factor * bc.practical_alpha(case))]
        if case.name in THEIR_T:
            alphas.append(("theirs", THEIR_T[case.name] - 1.0))
        for v in variants:
            if v.paper == "2025":
                plan.append((case, v, "theirs-SLSQP" + ("-H" if v.objective == "H" else "-LSQ"), "none", None))
            else:
                for choice, a in alphas:
                    plan.append((case, v, "theirs-FW", choice, a))
        for choice, a in alphas:
            plan.append((case, case, "ours-Halpha", choice, a))
    for case, v, method, choice, a in plan:
        print("%-22s %-18s alpha %-9s %s" % (case.name, method, choice, "-" if a is None else "%.4g" % a), flush=True)
    if args.list:
        return

    rows = []
    csv_path = os.path.join(args.out, "results.csv")
    for case, v, method, choice, a in plan:
        t0 = time.time()
        if method.startswith("theirs-SLSQP"):
            r = run_slsqp(v, args)
        elif method == "theirs-FW":
            r = run_fw(v, 1.0 + a, args)
        else:
            r = run_ours(case, a, args)
        if method.startswith("theirs-SLSQP"):
            rule, tol = slsqp_rule(v)
        elif method == "theirs-FW":
            rule, tol = FW_RULE, v.gap
        else:
            rule, tol = OUR_RULE, args.tol
        row = dict(case=case.name, method=method, alphaChoice=choice, alpha=a, rule=rule, tolerance=tol, **r)
        rows.append(row)
        tag = "%s__%s__%s" % (case.name, method, choice)
        with open(os.path.join(args.out, tag + ".json"), "w") as f:
            json.dump(row, f, indent=1)
        with open(csv_path, "w", newline="") as f:
            w = csv.DictWriter(f, FIELDS, extrasaction="ignore")
            w.writeheader()
            w.writerows(rows)
        fmt = lambda k: "%.2e" % row[k] if isinstance(row.get(k), float) else "-"
        print("%-22s %-18s alpha %-9s %-14s it %5s  %8.2fs  vol %s  centroid %s  stationary %s  (wall %.0fs)" % (
            case.name, method, choice, row["stop"], row.get("iterations", "-"), row.get("seconds", float("nan")),
            fmt("volume"), fmt("centroid"), fmt("stationary"), time.time() - t0), flush=True)

    with open(os.path.join(args.out, "summary.txt"), "w") as f:
        for row in rows:
            f.write("%-22s %-18s alpha %-9s %-10s %-14s it %5s  %9.2fs  vol %s  centroid %s  stationary %s\n" % (
                row["case"], row["method"], row["alphaChoice"], "-" if row["alpha"] is None else "%.4g" % row["alpha"],
                row["stop"], row.get("iterations", "-"), row.get("seconds", float("nan")),
                *["%.2e" % row[k] if isinstance(row.get(k), float) else "-" for k in ("volume", "centroid",
                                                                                      "stationary")]))
    print("wrote", csv_path)


if __name__ == "__main__":
    main()
