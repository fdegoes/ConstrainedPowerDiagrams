#!/usr/bin/env python3
"""Runs the methods of Bourne et al. on the synthetic benchmark.

Methods (benchmarks/bourne/bourne_ref.py, ported from their code):
  fw     Frank-Wolfe on the convex-order projection (2026), t = 100, R0 = 1,
         from the targets b (as in their aligned-seeds example)
  slsqp  SLSQP on H with the pairwise-distance and ball constraints (2025),
         d = 1e-3, R0 = sqrt(n), from b. Its n(n-1)/2 constraints make it
         impractical beyond a few hundred cells (--slsqp-max-n).

Each method stops on its own rule at 1e-8, and the other rule is measured
at its solution (the last iterate):
  slsqp  stops on rule 1, their relative objective change ptol = 1e-8
         (their 2025 notebook: 1e-10); rule 3 is measured at its solution.
  fw     has no rule 1 (the FW gap is not used as a rule): it stops on
         rule 3 at 1e-8, the volume and shift residuals of H (max-norm per
         cell), evaluated at every iterate with a tight OT solve that is not
         charged. The data are realizable and FW solves H_alpha with
         alpha = t - 1 = 99, so its residual of H stalls at about 1/t and
         the run usually ends at the time cap.
Every run stops at the wall-time cap --max-seconds (not counting the
evaluations), and reports its last iterate.

Writes <out>/<case>/<method>.json and <out>/theirs.csv.
"""

import argparse
import csv
import json
import math
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "bourne"))
import bourne_ref as ref  # noqa: E402
import generate as G  # noqa: E402
import their_rules as TR  # noqa: E402
from run_synthetic import NS, generator_error  # noqa: E402

TOL = 1e-8  # rule 3 tolerance
OT_EVAL_TOL = 1e-9  # relative to min(v), the OT of the common-rule evaluations (about what pysdot reaches;
                    # a missed tolerance is kept, as volume residual, not raised)

FIELDS = ["case", "n", "noise", "volumes", "seed", "method", "stoppingRule", "status", "stopped", "iterations",
          "seconds", "otSolves", "finalVolumeResidual", "finalShiftResidual", "rule3Met", "finalGeneratorError",
          "load"]


def residuals(c, X, alpha=math.inf):
    """Volume and shift residuals at the diagram with seeds X and the
    prescribed volumes (uncharged OT), of H (alpha = inf) or of H_alpha with
    sites cbar and generators p = alpha X (the FW iterates: X = y - B,
    bourne_cases.common_residuals)."""
    with ref.Evaluation():
        _, cen, V = ref.solve_OT_problem(X, c["volumes"], c["om"], relTol=OT_EVAL_TOL, strict=False)
    vb, b = c["volumes"], c["targets"]
    vol = float(np.max(np.abs(V - vb) / vb))
    g = 2 * V[:, None] * (b - cen)
    if math.isfinite(alpha):
        g -= 2 / alpha * vb[:, None] * (alpha * np.asarray(X) - b)
    shift = float(np.max(np.linalg.norm(g, axis=1) / (2 * vb ** 1.5)))
    return vol, shift


class Recorder:
    """Iterates of one run: counts them, keeps the last, and (FW) stops at
    the first iterate meeting rule 3 at TOL. Evaluations are not charged."""

    def __init__(self, c, args, stop_on_rule3):
        self.c, self.args, self.stop_on_rule3 = c, args, stop_on_rule3
        self.t0 = time.perf_counter()
        self.k = 0
        self.last = None
        self.met = False

    def elapsed(self):
        return time.perf_counter() - self.t0 - ref.COUNTERS.evalSeconds

    def record(self, X):
        self.k += 1
        self.last = np.array(X, copy=True)
        if self.stop_on_rule3:
            vol, shift = residuals(self.c, X)
            self.met = vol <= TOL and shift <= TOL

    def done(self):
        return self.met or self.elapsed() > self.args.max_seconds


def run_fw(c, args):
    """FW with t = 100, stopped on rule 3 at TOL with the residual of H (it
    has no rule 1; the FW gap is not used as a rule). The data are
    realizable and FW solves H_alpha with alpha = t - 1, so its residual of
    H stalls at about 1/t and the run usually ends at the time cap."""
    t = 100.0
    rec = Recorder(c, args, stop_on_rule3=True)
    y = t / (t - 1) * c["targets"]

    def cb(it, B, gap):
        rec.record(y - B)
        return rec.done()

    try:
        _, _, status = ref.convex_order_projection_frank_wolfe(c["om"], y, c["targets"].copy(), R0=1.0,
                                                               masses=c["volumes"], niter=args.max_iter,
                                                               tol=-1.0, callback=cb)
    except Exception as e:  # noqa: BLE001
        status = "error: %s" % e
    stopped = "rule 3" if rec.met else ("time cap" if rec.elapsed() > args.max_seconds else status)
    return rec, status, stopped


def run_slsqp(c, args):
    """SLSQP on H stopped on rule 1 (ptol = SLSQP_PTOL), or the time cap."""
    rec = Recorder(c, args, stop_on_rule3=False)
    n = c["n"]
    try:
        res, X = ref.max_H_solver(c["targets"].copy(), c["targets"], c["volumes"], c["om"], R0=math.sqrt(n),
                                  d=1e-3, max_iter=args.max_iter, ptol=TR.SLSQP_PTOL,
                                  callback=lambda X: (rec.record(X), rec.done())[1])
        status = str(res.message)
        stopped = "rule 1" if res.success else ("time cap" if rec.elapsed() > args.max_seconds else status)
        rec.last = X
    except Exception as e:  # noqa: BLE001
        status = stopped = "error: %s" % e
    return rec, status, stopped


def summarize(c, method, rec, status, stopped, load):
    """One row; rule 3 measured at the solution (the last iterate)."""
    row = {"case": c["name"], "n": c["n"], "noise": c["noise"], "volumes": c["vol"], "seed": c["seed"],
           "method": method, "stoppingRule": "rule 3" if method == "fw" else "rule 1", "status": status,
           "stopped": stopped, "iterations": rec.k, "seconds": rec.elapsed(), "otSolves": ref.COUNTERS.otSolves,
           "load": load}
    if rec.last is not None:
        try:
            vol, shift = residuals(c, rec.last)
            row.update(finalVolumeResidual=vol, finalShiftResidual=shift, rule3Met=vol <= TOL and shift <= TOL)
        except Exception:  # noqa: BLE001
            pass
        row["finalGeneratorError"] = generator_error(np.asarray(rec.last), c["truth"], c["n"])
    return row


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--n", nargs="+", type=int, default=list(NS))
    ap.add_argument("--noise", nargs="+", type=int, default=list(G.NOISE))
    ap.add_argument("--volumes", nargs="+", default=list(G.VOLUMES))
    ap.add_argument("--seeds", nargs="+", type=int, default=[1])
    ap.add_argument("--method", nargs="+", choices=["fw", "slsqp"], default=["slsqp", "fw"])
    ap.add_argument("--max-n", type=int, default=None, help="only the sizes of --n up to this")
    ap.add_argument("--max-iter", type=int, default=100000)
    ap.add_argument("--max-seconds", type=float, default=1800.0)
    ap.add_argument("--slsqp-max-n", type=int, default=500)
    ap.add_argument("--out", default=os.path.join(HERE, "results"))
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()
    if args.max_n is not None:
        args.n = [n for n in args.n if n <= args.max_n]
    ref.OT_REL_TOL["value"] = 1e-4  # their inner OT tolerance (2026), relative to min(v)

    rows, csv_path = [], os.path.join(args.out, "theirs.csv")
    if os.path.exists(csv_path):
        rows = list(csv.DictReader(open(csv_path)))
    done = {(r["case"], r["method"]) for r in rows}
    for n in args.n:
        for noise in args.noise:
            for vol in args.volumes:
                for seed in args.seeds:
                    c = G.make_case(args.out, n, noise, vol, seed)
                    c.update(noise=noise, vol=vol, seed=seed, om=ref.make_box(0, 0, 1, 1))
                    for method in args.method:
                        if method.startswith("slsqp") and n > args.slsqp_max_n:
                            continue
                        if (c["name"], method) in done and not args.force:
                            continue
                        load = os.getloadavg()[0]
                        ref.COUNTERS.reset()
                        if method == "fw":
                            rec, status, stopped = run_fw(c, args)
                        else:
                            rec, status, stopped = run_slsqp(c, args)
                        row = summarize(c, method, rec, status, stopped, load)
                        json.dump(row, open(os.path.join(c["dir"], method + ".json"), "w"))
                        rows = [r for r in rows if (r["case"], r["method"]) != (c["name"], method)] + [row]
                        print("%-26s %-6s stopped by %-10s %6d it %8.2fs  rule 3 at solution %.1e/%.1e (%s)  "
                              "load %.1f" % (c["name"], method, stopped[:10], rec.k, row["seconds"],
                                             row.get("finalVolumeResidual", float("nan")),
                                             row.get("finalShiftResidual", float("nan")),
                                             "met" if row.get("rule3Met") else "not met", load), flush=True)
                        with open(csv_path, "w", newline="") as f:
                            w = csv.DictWriter(f, FIELDS, extrasaction="ignore")
                            w.writeheader()
                            w.writerows(rows)


if __name__ == "__main__":
    main()
