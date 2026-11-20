#!/usr/bin/env python3
"""Tables of an overnight run (benchmarks/run_overnight.sh). Each method
stops on its own rule at 1e-8, and the other rule is measured at its
solution:
  our H         stops on rule 3 (common per-cell residual of H); rule 1
                measured at its solution
  their SLSQP   stops on rule 1 (relative objective change, ptol = 1e-8);
                rule 3 measured at its solution
  their FW      has no rule of its own (the FW gap is not a rule): stops on
                rule 3, or at the time cap

Reads <dir>/bourne-realizable*/*.json, <dir>/synthetic/{ours,theirs}.csv
and <dir>/bourne-approximate/results.csv, writes <dir>/REPORT_TABLES.md,
draws <dir>/timing.pdf (our H only) and <dir>/timing-with-theirs.pdf
(synthetic/plot_timing.py), and writes one log per set of results, one row
per run, with the same columns (LOG_FIELDS):
  <dir>/log-synthetic.csv           the 600 synthetic 2D cases (ours, and theirs on seed 1, n <= 500)
  <dir>/log-bourne-realizable.csv   the 3 realizable examples of Bourne et al.
  <dir>/log-bourne-approximate.csv  the approximate examples of Bourne et al.

Usage:
  external/bourne/.venv/bin/python benchmarks/report_overnight.py --dir benchmarks/overnight/<date>
"""

import argparse
import collections
import csv
import glob
import json
import os
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
TOL = "1e-08"


def read_csv(path):
    return list(csv.DictReader(open(path))) if os.path.exists(path) else []


def med(xs):
    return "%.3g" % np.median(xs) if xs else "–"


def mx(xs):
    return "%.3g" % max(xs) if xs else "–"


# One row per run in every log: whether the run met its rule, how it stopped
# (its rule, the iteration or time cap, or an error), its iterations and
# seconds, the rule and tolerance it stopped on, and the volume and shift
# residuals of eq. results-termination at its solution when measured.
LOG_FIELDS = ["case", "method", "converged", "stoppedBy", "iterations", "seconds", "rule", "tolerance",
              "volumeResidual", "shiftResidual"]
SYNTHETIC_FIELDS = ["n", "noise", "volumes", "seed"]
RULE_3 = "rule 3: volume and shift residuals of H (eq. results-termination)"
RULE_1 = "rule 1: SLSQP objective change <= tol |H(start)|, volumes within 0.1 % (Bourne 2025)"


def write_log(path, rows, extra=()):
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, ["case", *extra, *LOG_FIELDS[1:]], extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    print("wrote", path, "(%d runs)" % len(rows))


def synthetic_log(d):
    sd = os.path.join(d, "synthetic")
    rows = []
    for r in read_csv(os.path.join(sd, "ours.csv")):
        ok = r.get("converged") == "True"
        rows.append(dict(r, method="our H", converged=ok, stoppedBy="rule 3" if ok else r.get("status"),
                         rule=RULE_3, tolerance=float(TOL)))
    for r in read_csv(os.path.join(sd, "theirs.csv")):
        rule = r.get("stoppingRule")
        rows.append(dict(r, method="their SLSQP" if r.get("method") == "slsqp" else "their FW (t = 100)",
                         converged=r.get("stopped") == rule, stoppedBy=r.get("stopped"),
                         rule=RULE_1 if rule == "rule 1" else RULE_3, tolerance=float(TOL),
                         volumeResidual=r.get("finalVolumeResidual"), shiftResidual=r.get("finalShiftResidual")))
    write_log(os.path.join(d, "log-synthetic.csv"), rows, SYNTHETIC_FIELDS)


def realizable_log(d):
    rows = []
    for f in sorted(glob.glob(os.path.join(d, "bourne-realizable*", "*__*.json"))):
        r = json.load(open(f))
        method = r["method"]
        if method.startswith("theirs") and r.get("termination") == "theirs":
            r3 = r.get("rule3AtSolution") or {}
            rows.append(dict(case=r["case"], method="their SLSQP", converged=bool(r.get("converged")),
                             stoppedBy="rule 1" if r.get("converged") else r.get("status"),
                             iterations=r.get("iterations"), seconds=r.get("seconds"), rule=RULE_1,
                             tolerance=float(TOL), volumeResidual=r3.get("volumeResidual"),
                             shiftResidual=r3.get("shiftResidual")))
        elif method.startswith("theirs"):
            s, fin = r.get("stops", {}).get(TOL), r.get("final") or {}
            rows.append(dict(case=r["case"], method="their FW (t = 100)", converged=bool(s),
                             stoppedBy="rule 3" if s else r.get("status"),
                             iterations=s["iteration"] if s else fin.get("iteration"),
                             seconds=s["seconds"] if s else r.get("seconds"), rule=RULE_3, tolerance=float(TOL),
                             volumeResidual=fin.get("volumeResidual"), shiftResidual=fin.get("shiftResidual")))
        else:
            s = r.get("common", {}).get("stops", {}).get(TOL)
            rows.append(dict(case=r["case"], method="our " + method, converged=bool(s),
                             stoppedBy="rule 3" if s else r.get("solverStatus"),
                             iterations=s["iteration"] if s else r.get("solverIterations"),
                             seconds=s["seconds"] if s else r.get("solverSeconds"), rule=RULE_3,
                             tolerance=float(TOL), volumeResidual=r.get("volumeResidual", r.get("finalVolumeError")),
                             shiftResidual=r.get("shiftResidual")))
    write_log(os.path.join(d, "log-bourne-realizable.csv"), rows)


def approximate_log(d):
    rows = []
    for r in read_csv(os.path.join(d, "bourne-approximate", "results.csv")):
        method = r["method"] + ("" if r.get("alphaChoice") in (None, "", "none") else
                                " (alpha %s = %.4g)" % (r["alphaChoice"], float(r["alpha"])))
        rows.append(dict(case=r["case"], method=method, converged=r.get("stop") == "converged",
                         stoppedBy=r.get("stop"), iterations=r.get("iterations"), seconds=r.get("seconds"),
                         rule=r.get("rule"), tolerance=r.get("tolerance"), volumeResidual=r.get("volume"),
                         shiftResidual=r.get("stationary")))
    write_log(os.path.join(d, "log-bourne-approximate.csv"), rows)


def bourne_table(d, out):
    files = sorted(glob.glob(os.path.join(d, "bourne-realizable*", "*__*.json")))
    if not files:
        out.append("_No Bourne results._\n")
        return
    out.append("| case | method | stopped by | iterations | seconds | other rule at the solution |")
    out.append("|---|---|---|---|---|---|")
    for f in files:
        r = json.load(open(f))
        name, method = r["case"], r["method"]
        if method.startswith("theirs") and r.get("termination") == "theirs":
            r3 = r.get("rule3AtSolution") or {}
            other = ("rule 3: %.1e / %.1e (%s)" % (r3["volumeResidual"], r3["shiftResidual"],
                                                    "met" if r3["met"] else "not met")) if "met" in r3 else "–"
            out.append("| %s | SLSQP | rule 1 (%s) | %s | %.3g | %s |" % (
                name, "met" if r.get("converged") else r.get("status", "?")[:30], r.get("iterations"),
                r.get("seconds", float("nan")), other))
        elif method.startswith("theirs"):
            s = r.get("stops", {}).get(TOL)
            f = r.get("final") or {}
            stopped = "rule 3" if s else "cap (%s)" % r.get("status", "?")[:30]
            where = ("rule 3 residuals at the cap: %.1e / %.1e" % (f["volumeResidual"], f["shiftResidual"])
                     if not s and f else "– (no rule 1)")
            out.append("| %s | FW (t = 100) | %s | %s | %.3g | %s |" % (
                name, stopped, s["iteration"] if s else f.get("iteration", "–"),
                s["seconds"] if s else r.get("seconds", 0), where))
        else:
            s = r.get("common", {}).get("stops", {}).get(TOL)
            th = r.get("theirs")
            other = "–" if th is None else "rule 1: |Δobj|/|obj₀| %.1e (%s)" % (
                th["ratio"], "met" if th["met"] else "not met")
            out.append("| %s | our %s | %s | %s | %.3g | %s |" % (
                name, method, "rule 3" if s else "not met (%s)" % r.get("solverStatus", "?"),
                s["iteration"] if s else r.get("solverIterations"),
                s["seconds"] if s else r.get("solverSeconds", 0), other))
    out.append("")


def synthetic_tables(d, out):
    sd = os.path.join(d, "synthetic")
    ours, theirs = read_csv(os.path.join(sd, "ours.csv")), read_csv(os.path.join(sd, "theirs.csv"))

    out.append("### Our H on all cases (stops on rule 3 at 1e-8)\n")
    if ours:
        ok = [r for r in ours if r["converged"] == "True"]
        r1 = sum(r.get("rule1Met") == "True" for r in ours)
        s = lambda m: [float(r["seconds"]) for r in ok if int(r["n"]) <= m]
        out.append("| converged | median iterations | median seconds (n ≤ 100 / ≤ 1000 / all) | rule 1 met at the solution |")
        out.append("|---|---|---|---|")
        out.append("| %d/%d | %s | %s / %s / %s | %d/%d |" % (
            len(ok), len(ours), med([int(r["iterations"]) for r in ok]), med(s(100)), med(s(1000)), med(s(10 ** 9)),
            r1, len(ours)))
        out.append("")

    if not theirs:
        out.append("_No results of their methods._\n")
        return
    keys = {r["case"] for r in theirs}
    sub = [r for r in ours if r["case"] in keys]
    out.append("### Matched cases (seed 1, n ≤ 500)\n")
    out.append("| method | stops on | stopped by its rule | median seconds | max seconds | other rule met at the solution |")
    out.append("|---|---|---|---|---|---|")
    rows_at_cap = []
    if sub:
        ok = [float(r["seconds"]) for r in sub if r["converged"] == "True"]
        r1 = sum(r.get("rule1Met") == "True" for r in sub)
        out.append("| our H | rule 3 | %d/%d | %s | %s | rule 1: %d/%d |" % (len(ok), len(sub), med(ok), mx(ok), r1, len(sub)))
    for method, label, rule in (("slsqp", "SLSQP", "rule 1"), ("fw", "FW (t = 100)", "rule 3")):
        rows = [r for r in theirs if r["method"] == method]
        if not rows:
            continue
        mine = [float(r["seconds"]) for r in rows if r["stopped"] == rule]
        other = ("rule 3: %d/%d" % (sum(r.get("rule3Met") == "True" for r in rows), len(rows))
                 if method == "slsqp" else "– (no rule 1)")
        out.append("| %s | %s | %d/%d | %s | %s | %s |" % (label, rule, len(mine), len(rows), med(mine), mx(mine), other))
        rows_at_cap += [(label, r) for r in rows if r["stopped"] != rule]
    out.append("")
    if rows_at_cap:
        # How far from the solution the runs are that did not stop on their
        # rule (time cap or failure): rule 3 residuals at the last iterate.
        out.append("#### Runs of theirs that hit the cap: rule 3 residuals at their last iterate\n")
        out.append("| method | cases | volume residual (median / max) | centroid residual (median / max) |")
        out.append("|---|---|---|---|")
        for label in sorted({l for l, _ in rows_at_cap}):
            rr = [r for l, r in rows_at_cap if l == label and r.get("finalShiftResidual")]
            v = [float(r["finalVolumeResidual"]) for r in rr]
            c = [float(r["finalShiftResidual"]) for r in rr]
            out.append("| %s | %d | %s / %s | %s / %s |" % (
                label, len(rr), "%.1e" % np.median(v) if v else "–", "%.1e" % max(v) if v else "–",
                "%.1e" % np.median(c) if c else "–", "%.1e" % max(c) if c else "–"))
        out.append("")
        out.append("Per case: `synthetic/theirs.csv`, columns `finalVolumeResidual` and `finalShiftResidual`.\n")
    stops = collections.Counter((r["method"], r["stopped"]) for r in theirs)
    out.append("How their runs ended: %s\n" % ", ".join("%s %s: %d" % (m, s, k) for (m, s), k in sorted(stops.items())))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", required=True)
    args = ap.parse_args()
    out = ["# Overnight run: %s\n" % os.path.basename(os.path.abspath(args.dir))]
    env = os.path.join(args.dir, "environment.txt")
    if os.path.exists(env):
        out += ["```", open(env).read().rstrip(), "```", ""]
    out.append("Each method stops on its own rule at 1e-8; the other rule is measured at its solution. "
               "Rule 1 = SLSQP objective change (Bourne 2025, ptol = 1e-8); rule 3 = common per-cell residual of H "
               "(the data are realizable). FW has no rule of its own (the gap is not a rule), so it stops on rule 3; "
               "it solves H_alpha with alpha = t - 1 = 99, so its residual of H stalls at about 1/t.\n")
    out.append("## Bourne et al., realizable cases\n")
    bourne_table(args.dir, out)
    out.append("## Synthetic 2D\n")
    synthetic_tables(args.dir, out)
    path = os.path.join(args.dir, "REPORT_TABLES.md")
    open(path, "w").write("\n".join(out) + "\n")
    print("wrote", path)
    synthetic_log(args.dir)
    realizable_log(args.dir)
    approximate_log(args.dir)
    # Two plots: our H only (the scaling of our method), and with their
    # methods added where at least half of their cases stopped on their rule.
    plot = [sys.executable, os.path.join(HERE, "synthetic", "plot_timing.py"), "--results",
            os.path.join(args.dir, "synthetic")]
    subprocess.run(plot + ["--out", os.path.join(args.dir, "timing.pdf")], check=False)
    if os.path.exists(os.path.join(args.dir, "synthetic", "theirs.csv")):
        subprocess.run(plot + ["--theirs", "--out", os.path.join(args.dir, "timing-with-theirs.pdf")], check=False)


if __name__ == "__main__":
    main()
