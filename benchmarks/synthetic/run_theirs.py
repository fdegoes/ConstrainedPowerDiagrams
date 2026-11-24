#!/usr/bin/env python3
"""Runs the methods of Bourne et al. with their own termination rules.

Sets:
  synthetic   the 2D synthetic benchmark (generate.py): sizes NS, noise 0-2, uniform and
              Voronoi volumes, seeds 1-5 (600 cases). Default.
  realizable  the realizable cases of 2026 (bourne_cases.py): 2026-voronoi and 2026-aligned.

Methods (benchmarks/bourne/bourne_ref.py, ported from their code; the subgradient port reproduces
the error sequence of their convex_order_projection exactly over 30 iterations on a random case):
  fw       Frank-Wolfe on the convex-order projection (2026) at dilation t = 100, from the start.
           Stops on their FW gap at their tolerance and iteration cap: realizable, the case's
           gap and niter (their notebook); synthetic, 1e-14 and 4000 (their voronoi example).
  subgrad  the subgradient method on the dual of the same projection (2026,
           ConvexOrderProjectionUtils.convex_order_projection): tau = 2, steps tau/(i+1),
           niter 600. Stops when the change of its running average of the iterates is below
           1e-5 (their rule). Same dilation and start as fw.
  slsqp    SLSQP on H (2025), from the targets, stopped on their rule 1 (relative change of the
           objective, ptol = 1e-8) or on rule 3 with --slsqp-stop rule3. Synthetic only, n up to
           --slsqp-max-n (its n(n-1)/2 constraints make it impractical beyond a few hundred cells).

Each run stops on its own rule, at the iteration cap, or at the time cap (case_cap: max(60 s, our time) on synthetic cases,
1800 s on realizable; --max-seconds overrides; not counting
the uncharged evaluations of our residuals). At its last iterate we report its own rule's value
(theirRuleValue: the FW gap, or the subgradient change) and our residuals of H (volume and shift,
max-norm per cell), and whether our rule (1e-8 on both) holds there.

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
import bourne_cases as bc  # noqa: E402
import bourne_ref as ref  # noqa: E402
import generate as G  # noqa: E402
import their_rules as TR  # noqa: E402
from run_synthetic import NS, generator_error  # noqa: E402

TOL = 1e-8  # our rule 3 tolerance
OT_EVAL_TOL = 1e-9  # relative to min(v), the OT of the evaluations (about what pysdot reaches; a missed
                    # tolerance is kept, as volume residual, not raised)
DILATION = 100.0  # t of the 2026 methods (the realizable cases and the synthetic benchmark)
SYNTH_FW_TOL, SYNTH_FW_NITER = 1e-14, 4000  # their 2026-voronoi settings, used for the synthetic cases
SUBGRAD_TAU, SUBGRAD_TOL, SUBGRAD_NITER = 2.0, 1e-5, 600  # their convex_order_projection defaults
CAPS = {"realizable": 1800.0}  # time caps of a run, in seconds
SYNTH_MIN_CAP = 60.0  # synthetic cases: cap = max(SYNTH_MIN_CAP, our time on the case) (results/ours.csv)
OURS = os.path.join(HERE, "results", "ours.csv")
OUT = {"synthetic": os.path.join(HERE, "results"), "realizable": os.path.join(os.path.dirname(HERE), "bourne",
                                                                             "results", "theirs-realizable")}
RULE = {"fw": "FW gap (Bourne 2026)", "subgrad": "change of the running average (Bourne 2026)",
        "slsqp": "rule 1 (SLSQP ftol, Bourne 2025)"}

FIELDS = ["case", "set", "n", "noise", "volumes", "seed", "method", "stoppingRule", "status", "stopped", "iterations",
          "seconds", "capSeconds", "otSolves", "theirRuleValue", "theirRuleTol", "finalVolumeResidual", "finalShiftResidual",
          "rule3Met", "finalGeneratorError", "load"]


def residuals(c, X, alpha=math.inf):
    """Volume and shift residuals at the diagram with seeds X and the prescribed volumes (uncharged OT),
    of H (alpha = inf) or of H_alpha with sites cbar and generators p = alpha X (bourne_cases.common_residuals)."""
    with ref.Evaluation():
        _, cen, V = ref.solve_OT_problem(X, c["volumes"], c["om"], relTol=OT_EVAL_TOL, strict=False)
    vb, b = c["volumes"], c["targets"]
    vol = float(np.max(np.abs(V - vb) / vb))
    g = 2 * V[:, None] * (b - cen)
    if math.isfinite(alpha):
        g -= 2 / alpha * vb[:, None] * (alpha * np.asarray(X) - b)
    shift = float(np.max(np.linalg.norm(g, axis=1) / (2 * vb ** 1.5)))
    return vol, shift


def subgrad_sites(Z):
    """The generators of the diagram of a subgradient iterate: its rescaling, as in their code."""
    return (Z - np.min(Z) + .02) / (np.max(Z) - np.min(Z) + .02)


def our_seconds():
    """Our time per synthetic case (results/ours.csv)."""
    return {r["case"]: float(r["seconds"]) for r in csv.DictReader(open(OURS))} if os.path.exists(OURS) else {}


def case_cap(c, override, ours):
    """Time cap of a run on case c: the override of --max-seconds, else for the realizable cases CAPS, and for
the synthetic cases max(SYNTH_MIN_CAP, our time on the case), so that theirs get as long as ours."""
    if override is not None:
        return override
    if c["set"] == "realizable":
        return CAPS["realizable"]
    return max(SYNTH_MIN_CAP, ours.get(c["name"], 0.0))


class Clock:
    """Wall time of one run, not counting the uncharged evaluations."""

    def __init__(self, max_seconds):
        self.t0, self.max_seconds = time.perf_counter(), max_seconds

    def elapsed(self):
        return time.perf_counter() - self.t0 - ref.COUNTERS.evalSeconds

    def out(self):
        return self.elapsed() > self.max_seconds


def their_stop(status):
    """The stop we report for their status (converged, stopped, max iterations, or an OT error of the subgradient)."""
    return {"converged": "their rule", "stopped": "time cap", "max iterations": "iteration cap"}.get(status, status)


def run_fw(c, args):
    """Their FW at dilation t, from the start, stopped by their gap rule, the cap or the time cap."""
    t = c["t"]
    y = t / (t - 1) * c["targets"]
    clock = Clock(args.max_seconds)
    gap = {"value": None}

    def cb(it, B, value):
        gap["value"] = value
        return clock.out()

    try:
        X, k, status = ref.convex_order_projection_frank_wolfe(c["om"], y, c["start"].copy(), R0=c["R0"],
                                                               masses=c["volumes"], niter=c["niter"], tol=c["gap"],
                                                               callback=cb)
        X = y - X
        stopped = their_stop(status)
    except Exception as e:  # noqa: BLE001
        X, k, stopped = None, None, "error: %s" % e
    return dict(X=X, k=k, stopped=stopped, value=gap["value"], tol=c["gap"], seconds=clock.elapsed())


def run_subgrad(c, args):
    """Their subgradient method on the dual of the projection of the dilated targets, from the start."""
    t = c["t"]
    y = t / (t - 1) * c["targets"]
    clock = Clock(args.max_seconds)
    try:
        Z, k, status, error = ref.convex_order_subgradient(c["om"], c["start"].copy(), y, c["volumes"],
                                                           niter=SUBGRAD_NITER, tau=SUBGRAD_TAU, tol=SUBGRAD_TOL,
                                                           callback=lambda i, Z, e: clock.out())
        X = subgrad_sites(Z)
        stopped = their_stop(status)
    except Exception as e:  # noqa: BLE001
        X, k, error, stopped = None, None, None, "error: %s" % e
    return dict(X=X, k=k, stopped=stopped, value=error, tol=SUBGRAD_TOL, seconds=clock.elapsed())


class Recorder:
    """Iterates of an SLSQP run: counts them, keeps the last, and (rule 3) stops at the first iterate meeting
    our rule at TOL. Evaluations are not charged."""

    def __init__(self, c, args, stop_on_rule3):
        self.c, self.args, self.stop_on_rule3 = c, args, stop_on_rule3
        self.clock = Clock(args.max_seconds)
        self.k = 0
        self.last = None
        self.met = False

    def elapsed(self):
        return self.clock.elapsed()

    def record(self, X):
        self.k += 1
        self.last = np.array(X, copy=True)
        if self.stop_on_rule3:
            vol, shift = residuals(self.c, X)
            self.met = vol <= TOL and shift <= TOL

    def done(self):
        return self.met or self.clock.out()


def run_slsqp(c, args):
    """SLSQP on H stopped on rule 1 (ptol = TR.SLSQP_PTOL), or on rule 3 at TOL with --slsqp-stop rule3."""
    rule3 = args.slsqp_stop == "rule3"
    rec = Recorder(c, args, stop_on_rule3=rule3)
    n = c["n"]
    # With rule 3 as the stop, rule 1 is switched off (ptol = 0), so that SLSQP runs until our rule holds.
    ptol = 0.0 if rule3 else TR.SLSQP_PTOL
    try:
        res, X = ref.max_H_solver(c["targets"].copy(), c["targets"], c["volumes"], c["om"], R0=math.sqrt(n),
                                  d=1e-3, max_iter=args.max_iter, ptol=ptol,
                                  callback=lambda X: (rec.record(X), rec.done())[1])
        stopped = "rule 3" if rec.met else ("rule 1" if res.success else
                                             ("time cap" if rec.elapsed() > args.max_seconds else str(res.message)))
        return dict(X=X, k=rec.k, stopped=stopped, value=None, tol=TR.SLSQP_PTOL, seconds=rec.elapsed(),
                    status=str(res.message))
    except Exception as e:  # noqa: BLE001
        return dict(X=None, k=rec.k, stopped="error: %s" % e, value=None, tol=None, seconds=rec.elapsed(),
                    status="error: %s" % e)


def summarize(c, method, r, load):
    """One row. The residuals are of H at the last iterate of the run, evaluated uncharged."""
    row = {"case": c["name"], "set": c["set"], "n": c["n"], "noise": c["noise"], "volumes": c["vol"],
           "seed": c["seed"], "method": method, "stoppingRule": RULE[method], "status": r.get("status", ""),
           "stopped": r["stopped"], "iterations": r["k"], "seconds": r["seconds"], "otSolves": ref.COUNTERS.otSolves,
           "theirRuleValue": r["value"], "theirRuleTol": r["tol"], "load": load}
    if r["X"] is not None:
        try:
            vol, shift = residuals(c, r["X"])
            row.update(finalVolumeResidual=vol, finalShiftResidual=shift, rule3Met=vol <= TOL and shift <= TOL)
        except Exception:  # noqa: BLE001
            pass
        if c.get("truth") is not None:
            row["finalGeneratorError"] = generator_error(np.asarray(r["X"]), c["truth"], c["n"])
    return row


def synthetic_cases(args):
    for n in args.n:
        for noise in args.noise:
            for vol in args.volumes:
                for seed in args.seeds:
                    c = G.make_case(args.out, n, noise, vol, seed)
                    c.update(set="synthetic", noise=noise, vol=vol, seed=seed, om=ref.make_box(0, 0, 1, 1),
                             t=DILATION, start=c["targets"].copy(), gap=SYNTH_FW_TOL, niter=SYNTH_FW_NITER, R0=1.0)
                    yield c


def realizable_cases(args):
    for case in bc.make_cases():
        if not (case.realizable and case.paper == "2026"):
            continue
        d = os.path.join(args.out, case.name)
        os.makedirs(d, exist_ok=True)
        yield dict(set="realizable", name=case.name, dir=d, n=case.n, noise=None, vol=None, seed=None, om=case.om,
                   volumes=case.v, targets=case.b, truth=case.truth, start=case.start, t=case.t, gap=case.gap,
                   niter=case.niter, R0=case.R0)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--set", choices=["synthetic", "realizable"], default="synthetic")
    ap.add_argument("--method", nargs="+", choices=["fw", "subgrad", "slsqp"], default=["fw", "subgrad"])
    ap.add_argument("--n", nargs="+", type=int, default=list(NS))
    ap.add_argument("--noise", nargs="+", type=int, default=list(G.NOISE))
    ap.add_argument("--volumes", nargs="+", default=list(G.VOLUMES))
    ap.add_argument("--seeds", nargs="+", type=int, default=[1, 2, 3, 4, 5])
    ap.add_argument("--max-n", type=int, default=None, help="only the sizes of --n up to this")
    ap.add_argument("--max-iter", type=int, default=100000, help="iteration cap of SLSQP")
    ap.add_argument("--max-seconds", type=float, default=None,
                    help="time cap of every run (default: max(60 s, our time) on synthetic cases, 1800 s realizable)")
    ap.add_argument("--slsqp-max-n", type=int, default=500)
    ap.add_argument("--slsqp-stop", choices=["rule1", "rule3"], default="rule1",
                    help="SLSQP stops on its own rule 1 (default) or on rule 3 at 1e-8 (our termination)")
    ap.add_argument("--out", default=None, help="default: benchmarks/synthetic/results or benchmarks/bourne/"
                                                "results/theirs-realizable")
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()
    args.out = args.out or OUT[args.set]
    override, ours = args.max_seconds, our_seconds() if args.set == "synthetic" else {}
    if args.set == "realizable":
        args.method = [m for m in args.method if m != "slsqp"]
    if args.max_n is not None:
        args.n = [n for n in args.n if n <= args.max_n]
    ref.OT_REL_TOL["value"] = 1e-4  # their inner OT tolerance (2026), relative to min(v)
    os.makedirs(args.out, exist_ok=True)

    rows, csv_path = [], os.path.join(args.out, "theirs.csv")
    if os.path.exists(csv_path):
        rows = list(csv.DictReader(open(csv_path)))
    done = {(r["case"], r["method"]) for r in rows}
    cases = synthetic_cases(args) if args.set == "synthetic" else realizable_cases(args)
    for c in cases:
        for method in args.method:
            if method == "slsqp" and c["n"] > args.slsqp_max_n:
                continue
            if (c["name"], method) in done and not args.force:
                continue
            args.max_seconds = case_cap(c, override, ours)
            load = os.getloadavg()[0]
            ref.COUNTERS.reset()
            if method == "fw":
                r = run_fw(c, args)
            elif method == "subgrad":
                r = run_subgrad(c, args)
            else:
                r = run_slsqp(c, args)
            row = summarize(c, method, r, load)
            row["capSeconds"] = args.max_seconds
            json.dump(row, open(os.path.join(c["dir"], method + ".json"), "w"))
            rows = [x for x in rows if (x["case"], x["method"]) != (c["name"], method)] + [row]
            print("%-26s %-8s %-14s %6s it %8.2fs (cap %.0f)  own rule %s  ours %.1e/%.1e (%s)  load %.1f" % (
                c["name"], method, r["stopped"], r["k"], row["seconds"], args.max_seconds,
                "-" if r["value"] is None else "%.1e" % r["value"],
                row.get("finalVolumeResidual", float("nan")), row.get("finalShiftResidual", float("nan")),
                "met" if row.get("rule3Met") else "not met", load), flush=True)
            with open(csv_path, "w", newline="") as f:
                w = csv.DictWriter(f, FIELDS, extrasaction="ignore")
                w.writeheader()
                w.writerows(rows)


if __name__ == "__main__":
    main()
