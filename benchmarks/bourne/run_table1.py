#!/usr/bin/env python3
"""Table 1 runs: four methods on the two realizable benchmarks, every run stopped
on rule 3 at 1e-8 (or at a time cap), with the residual history of every iterate.

Benchmarks (bourne_cases.make_cases()):
  voronoi  n = 20, the data of 2025-invert and 2026-voronoi (the same Voronoi
           diagram); each method with the case settings of its own paper.
  cross    n = 728, the data of 2026-aligned.

Methods (the paper calls our joint dual G; the code calls it H):
  ours     max H with build/benchBourne, as bourne_cases.run_mine: sites = targets,
           start 0, rule 3 at 1e-8, at most 100 Newton iterations per attempt,
           restarts as the solver default.
  slsqp    Bourne et al. 2025, SLSQP on their H (bourne_ref.max_H_solver), with
           the settings of 2025-invert: pairwise distance d = 1e-3, ball radius
           R0 = sqrt(n). Voronoi from the random start of 2025-invert; cross
           (never run by them) from the targets, as their perturbed examples and
           our synthetic runs. Their rule 1 is disabled (ptol = 0) and read off
           the history (first iterate where the change of H since the previous
           iterate is <= 1e-8 |H(start)|). The inner OT is solved to 1e-9 min V,
           as in bourne_cases.run_theirs_common (their 1e-3 caps the accuracy of
           the gradient at about 1e-3).
  fw       Bourne et al. 2026, Frank-Wolfe on the convex-order projection
           (bourne_ref.convex_order_projection_frank_wolfe), t = 100, R0 = 1, from
           the case start (voronoi: B0 = the true seeds, their notebook; cross:
           B0 = the targets, their initial condition 2). Gap rule disabled; the gap
           of every iterate is recorded and their stop (gap < 1e-14 on voronoi,
           1e-16 on cross, their notebook) read off the history.
  subgrad  Bourne et al. 2026, subgradient method on the dual of the same
           projection (bourne_ref.convex_order_subgradient, a port of their
           ConvexOrderProjectionUtils.convex_order_projection without the lower
           bound): steps tau/(k+1) with tau = 2 (their code's default, and the
           step 2/(k+1) of their Lemma 6.2), lengthscale 1, cluster tolerance 1e-6.
           Voronoi: t = 1250 from the random generators of 2025-invert (their
           Sec. 6.3.1, Fig. 8: t = 1250, a random initial distribution of
           generators; the notebook draws the same random start). Cross: t = 100
           from Y1 = the targets (their Sec. 6.3.1, Figs. 11-14: t = 100, initial
           condition 2, the same as FW). Their rule (change of the running average
           of the iterates < 1e-5) is disabled and read off the history; their
           stop on an OT error (|V - Vbar|_2 > 1e-6) is kept and reported as an
           error.

Termination of their methods: rule 3 at 1e-8 (bourne_cases.common_residuals, the
residual of H, evaluated at every iterate on the power diagram of the iterate's
generators with the prescribed volumes, an OT solve to 1e-9 min V that is not
charged: bourne_ref.Evaluation), a time cap (--max-seconds, default 1800 s, of
charged time, checked between iterations) or an iteration cap (--max-iter). The
time of a run excludes these evaluations. Each run is a child process; one that
overruns its cap inside an iteration (e.g. an SLSQP iteration on the cross) is
killed at cap + --grace seconds of charged time and reported as "time cap" with
its history up to then.

History: rule 3 is evaluated at every iterate (iteration 0 = the start); the history
stores every iterate up to 10^4, then every 10th up to 10^5, every 100th up to 10^6,
and so on (keep_iterate), plus the last iterate and those where rule 3 or the method's
own rule first holds. Columns: iteration, solves (OT solves of their methods, diagram
builds of ours, cumulative), seconds (charged, cumulative), volume and centroid
residuals of rule 3, and own (the value of the method's own rule: SLSQP |dH|/|H0|, FW
gap, subgradient change).

Output (--out DIR): <benchmark>__<method>.json per run (summary and history),
<benchmark>__<method>.log (the child's output), results.csv, table.tex (rows of
tab:realizable in latex/results.tex). Re-running a subset updates its runs only;
--report rewrites results.csv and table.tex from the JSON files.

  MPLCONFIGDIR=/tmp/mplcfg external/bourne/.venv/bin/python benchmarks/bourne/run_table1.py
      [--benchmark voronoi cross] [--method ours slsqp fw subgrad] [--max-seconds 1800]
      [--out DIR] [--quick] [--report]
"""

import argparse
import csv
import json
import math
import os
import platform
import signal
import subprocess
import sys
import tempfile
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import bourne_cases as bc  # noqa: E402
import bourne_ref as ref  # noqa: E402

TOL = 1e-8  # rule 3
EVAL_TOL = bc.OT_MIN_REL_TOL  # OT of the evaluations and of SLSQP's inner solve, relative to min V
BENCHMARKS = ["voronoi", "cross"]
METHODS = ["ours", "slsqp", "fw", "subgrad"]
LABEL = {"voronoi": "Voronoi", "cross": "cross"}
SUBGRAD = {"voronoi": dict(t=1250.0, start="2025-invert"), "cross": dict(t=100.0, start="2026-aligned")}
SUBGRAD_TAU, SUBGRAD_RULE_TOL = 2.0, 1e-5
FIELDS = ["benchmark", "method", "n", "stop", "status", "iterations", "solves", "seconds", "volumeResidual",
          "centroidResidual", "rule3Met", "ownRule", "ownRuleIteration", "ownRuleSeconds", "ownRuleVolumeResidual",
          "ownRuleCentroidResidual", "wallSeconds", "evalSeconds", "load"]
OUT = os.path.join(HERE, "results", "table1")


def cases():
    return {c.name: c for c in bc.make_cases()}


# --------------------------------------------------------------------------
# Child: one run
# --------------------------------------------------------------------------

def keep_iterate(k):
    """Iterates stored in the history (all are evaluated): every one up to KEEP_ALL, then
    about 9000 per decade (every 10^(m - 3)-th from 10^m on), so that the runs of their
    methods that last to the time cap (10^5 to 10^6 iterations on the Voronoi case) stay
    a few MB. The last iterate, and those where rule 3 or the method's own rule first holds,
    are always stored."""
    if k <= KEEP_ALL:
        return True
    return k % 10 ** (int(math.log10(k)) - 3) == 0


KEEP_ALL = 10000


class Run:
    """History of one run: appended to <tag>.history.jsonl as it goes (so that a
    killed run keeps it), with the charged clock and the uncharged evaluations.
    own_test(k, own): whether the method's own rule holds at iterate k."""

    def __init__(self, path, meta, max_seconds, grace, own_test=None):
        self.f = open(path, "w", buffering=1)
        self.f.write(json.dumps(dict(meta=meta)) + "\n")
        self.t0 = time.perf_counter()
        self.max_seconds, self.grace = max_seconds, grace
        self.cols = {k: [] for k in ("iteration", "solves", "seconds", "volume", "centroid", "own")}
        self.met = False
        self.own_test, self.own_met = own_test, False
        self.pending = None  # the last iterate, when it is not stored
        self.arm()

    def arm(self):
        """Hard stop at cap + grace of charged time: SIGALRM with its default action ends the
        process even inside a compiled call (one SLSQP iteration on the cross), which the cap,
        checked between iterations, cannot. Re-armed at every iterate, so that the uncharged
        evaluations do not count."""
        signal.setitimer(signal.ITIMER_REAL, max(1.0, self.max_seconds + self.grace - self.elapsed()))

    def elapsed(self):
        return time.perf_counter() - self.t0 - ref.COUNTERS.evalSeconds

    def out_of_time(self):
        return self.elapsed() > self.max_seconds

    def add(self, k, solves, seconds, vol, cen, own=None):
        row = dict(iteration=int(k), solves=int(solves), seconds=float(seconds), volume=float(vol),
                   centroid=float(cen), own=None if own is None else float(own))
        self.met = vol <= TOL and cen <= TOL
        first_own = (not self.own_met and own is not None and self.own_test is not None
                     and self.own_test(k, own))
        self.own_met = self.own_met or first_own
        if keep_iterate(k) or self.met or first_own:
            self._store(row)
            self.pending = None
        else:
            self.pending = row
        self.arm()
        return self.met

    def _store(self, row):
        for key in self.cols:
            self.cols[key].append(row[key])
        self.f.write(json.dumps(dict(row, evalSeconds=ref.COUNTERS.evalSeconds)) + "\n")

    def close(self):
        """Stores the last iterate."""
        if self.pending is not None:
            self._store(self.pending)
            self.pending = None

    def evaluate(self, case, X):
        """Rule-3 residuals (of H) of the diagram of generators X with the prescribed volumes, not charged."""
        with ref.Evaluation():
            _, c, V = ref.solve_OT_problem(np.asarray(X, dtype=float), case.v, case.om, relTol=EVAL_TOL, strict=False)
        vol, cen = bc.common_residuals(case, V, c, np.asarray(X), math.inf)
        return vol, cen, V, c


def run_ours(case, args, run):
    lines = ["box %r %r %r %r" % tuple(float(x) for x in case.box), "n %d" % case.n, "solver H", "alpha 1e300",
             "start 0", "volumetol %r" % TOL, "shifttol %r" % TOL, "maxiter %d" % args.max_iter_ours,
             "lighttrajectory",
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
    for k, it in enumerate(traj):
        V, c = np.asarray(it["volumes"]), np.reshape(it["centroids"], (-1, 2))
        vol, cen = bc.common_residuals(case, V, c, np.reshape(it["generators"], (-1, 2)), math.inf)
        run.add(k, it["numBuilds"], it["seconds"], vol, cen)
    stop = "rule 3" if r["converged"] else ("iteration cap" if r["status"] == "max iterations" else
                                             "error: %s" % r["status"])
    return dict(stop=stop, status=r["status"], iterations=int(r["iterations"]), solves=int(r["numBuilds"]),
                seconds=float(r["seconds"]), attempts=r.get("attempts"),
                solverResiduals=dict(volume=r["volumeResidual"], shift=r["shiftResidual"]))


def run_slsqp(case, start, args, run):
    ref.OT_REL_TOL["value"] = EVAL_TOL
    H = []

    def record(X):
        vol, cen, V, c = run.evaluate(case, X)
        H.append(float(np.sum(case.v * np.sum((case.b - c) * np.asarray(X), axis=1))))
        own = None if len(H) < 2 else abs(H[-1] - H[-2]) / abs(H[0])
        return run.add(len(H) - 1, ref.COUNTERS.otSolves, run.elapsed(), vol, cen, own)

    record(start)
    if run.met:
        return dict(stop="rule 3", status="met at the start")
    R0 = math.sqrt(case.n)
    res, _ = ref.max_H_solver(start, case.b, case.v, case.om, R0=R0, d=1e-3, max_iter=args.max_iter_slsqp, ptol=0.0,
                              callback=lambda X: record(X) or run.out_of_time())
    msg = str(res.message)
    if run.met:
        stop = "rule 3"
    elif run.out_of_time():
        stop = "time cap"
    elif "Iteration limit" in msg:
        stop = "iteration cap"
    else:
        stop = "error: %s" % msg
    return dict(stop=stop, status=msg)


def run_fw(case, args, run):
    y = case.t / (case.t - 1) * case.b

    def cb(it, B, gap):
        vol, cen, _, _ = run.evaluate(case, y - B)
        # The LMO of this iterate is solved when the callback is called.
        return run.add(it, ref.COUNTERS.otSolves, run.elapsed(), vol, cen, gap) or run.out_of_time()

    _, _, status = ref.convex_order_projection_frank_wolfe(case.om, y, case.start.copy(), R0=case.R0, masses=case.v,
                                                           niter=args.max_iter, tol=-1.0, callback=cb)
    stop = "rule 3" if run.met else ("time cap" if status == "stopped" else
                                     "iteration cap" if status == "max iterations" else status)
    return dict(stop=stop, status=status)


def subgrad_sites(Z):
    """The generators of the diagram of a subgradient iterate: its rescaling, as in their code."""
    return (Z - np.min(Z) + .02) / (np.max(Z) - np.min(Z) + .02)


def run_subgrad(case, start, t, args, run):
    y = t / (t - 1) * case.b
    vol, cen, _, _ = run.evaluate(case, subgrad_sites(start))
    if run.add(0, 0, run.elapsed(), vol, cen):
        return dict(stop="rule 3", status="met at the start")

    def cb(i, Z, error):
        vol, cen, _, _ = run.evaluate(case, subgrad_sites(Z))
        return run.add(i + 1, ref.COUNTERS.otSolves, run.elapsed(), vol, cen, error) or run.out_of_time()

    _, _, status, _ = ref.convex_order_subgradient(case.om, np.array(start, dtype=float), y, case.v,
                                                   niter=args.max_iter, tau=SUBGRAD_TAU, tol=-1.0, callback=cb)
    stop = "rule 3" if run.met else ("time cap" if status == "stopped" else
                                     "iteration cap" if status == "max iterations" else "error: %s" % status)
    return dict(stop=stop, status=status)


def settings(benchmark, method, C):
    """(case of the data, start case, settings) of a run."""
    if method == "ours":
        case = C["2026-voronoi"] if benchmark == "voronoi" else C["2026-aligned"]
        return case, dict(objective="H", sites="targets", start=0, tol=TOL, maxiterPerAttempt=None, restarts="default")
    if method == "slsqp":
        case = C["2025-invert"] if benchmark == "voronoi" else C["2026-aligned"]
        start = "2025-invert random start" if benchmark == "voronoi" else "targets"
        return case, dict(d=1e-3, R0="sqrt(n)", ptol=0.0, rule1Ptol=bc.SLSQP_PTOL, innerOTRelTol=EVAL_TOL, start=start)
    if method == "fw":
        case = C["2026-voronoi"] if benchmark == "voronoi" else C["2026-aligned"]
        start = "true seeds (B0)" if benchmark == "voronoi" else "targets (B0)"
        return case, dict(t=case.t, R0=case.R0, gapRuleTol=case.gap, start=start)
    case = C["2026-voronoi"] if benchmark == "voronoi" else C["2026-aligned"]
    s = SUBGRAD[benchmark]
    start = "2025-invert random generators" if benchmark == "voronoi" else "targets (Y1)"
    return case, dict(t=s["t"], tau=SUBGRAD_TAU, lengthscale=1.0, tolcluster=1e-6, ruleTol=SUBGRAD_RULE_TOL,
                      start=start)


OWN_RULE = {"slsqp": "rule 1 (SLSQP ftol, ptol 1e-8)", "fw": "FW gap", "subgrad": "change of the running average"}


def own_rule_index(method, cols, s):
    own = cols["own"]
    if method == "slsqp":
        tol = s["rule1Ptol"]
        return next((k for k, x in enumerate(own) if k > 0 and x is not None and x <= tol), None)
    if method == "fw":
        tol = s["gapRuleTol"]
        return next((k for k, x in enumerate(own) if x is not None and x < tol), None)
    if method == "subgrad":
        tol = s["ruleTol"]
        return next((k for k, x in enumerate(own) if x is not None and x < tol), None)
    return None


def finish(benchmark, method, n, s, cols, summary, wall, evals):
    """The JSON of a run from its history and the summary of the child (or of the parent, if killed)."""
    last = len(cols["iteration"]) - 1
    out = dict(benchmark=benchmark, method=method, n=n, settings=s, tol=TOL, wallSeconds=wall, evalSeconds=evals,
               load=os.getloadavg()[0], host=platform.node())
    out.update(summary)
    if last >= 0:
        out.setdefault("iterations", cols["iteration"][last])
        out.setdefault("solves", cols["solves"][last])
        out.setdefault("seconds", cols["seconds"][last])
        out.update(volumeResidual=cols["volume"][last], centroidResidual=cols["centroid"][last],
                   rule3Met=cols["volume"][last] <= TOL and cols["centroid"][last] <= TOL)
    k = own_rule_index(method, cols, s)
    out["ownRule"] = None if method == "ours" else dict(
        rule=OWN_RULE[method], index=k, met=k is not None,
        **({} if k is None else dict(iteration=cols["iteration"][k], solves=cols["solves"][k],
                                     seconds=cols["seconds"][k], volumeResidual=cols["volume"][k],
                                     centroidResidual=cols["centroid"][k])))
    out["history"] = cols
    return out


def worker(args):
    C = cases()
    b, m = args.benchmark[0], args.method[0]
    case, s = settings(b, m, C)
    if m == "ours":
        s["maxiterPerAttempt"] = args.max_iter_ours
    tag = os.path.join(args.out, "%s__%s" % (b, m))
    ref.COUNTERS.reset()
    own_test = {"slsqp": lambda k, x: k > 0 and x <= s["rule1Ptol"], "fw": lambda k, x: x < s["gapRuleTol"],
                "subgrad": lambda k, x: x < s["ruleTol"]}.get(m)
    run = Run(tag + ".history.jsonl", dict(benchmark=b, method=m, n=case.n, settings=s), args.max_seconds, args.grace,
              own_test)
    t0 = time.time()
    try:
        if m == "ours":
            summary = run_ours(case, args, run)
        elif m == "slsqp":
            summary = run_slsqp(case, case.start if b == "voronoi" else case.b.copy(), args, run)
        elif m == "fw":
            summary = run_fw(case, args, run)
        else:
            start = C[SUBGRAD[b]["start"]].start if b == "voronoi" else case.b.copy()
            summary = run_subgrad(case, start, SUBGRAD[b]["t"], args, run)
        if m != "ours":
            summary.update(seconds=run.elapsed())
    except Exception as e:  # noqa: BLE001
        summary = dict(stop="error: %s" % e, status="error: %s" % e)
        if m != "ours":
            summary.update(seconds=run.elapsed())
    signal.setitimer(signal.ITIMER_REAL, 0)  # the run is over: writing its JSON is not capped
    run.close()
    if m != "ours" and run.cols["iteration"]:
        summary.setdefault("iterations", run.cols["iteration"][-1])
        summary["solves"] = ref.COUNTERS.otSolves
    out = finish(b, m, case.n, s, run.cols, summary, time.time() - t0, ref.COUNTERS.evalSeconds)
    with open(tag + ".json", "w") as f:
        json.dump(out, f)
    run.f.close()
    os.remove(tag + ".history.jsonl")
    print("%-8s %-8s %-14s %7s it %8.3f s  residuals %.2e / %.2e" % (
        b, m, out["stop"][:14], out.get("iterations"), out.get("seconds", float("nan")),
        out.get("volumeResidual", float("nan")), out.get("centroidResidual", float("nan"))), flush=True)


# --------------------------------------------------------------------------
# Parent: runs the children, enforces the caps, writes the tables
# --------------------------------------------------------------------------

def _last_line(path):
    try:
        with open(path, "rb") as f:
            f.seek(0, 2)
            size = f.tell()
            f.seek(max(0, size - 8192))
            lines = f.read().splitlines()
        return json.loads(lines[-1]) if lines else None
    except (OSError, ValueError):
        return None


def _from_jsonl(path):
    meta, cols, evals = None, {k: [] for k in ("iteration", "solves", "seconds", "volume", "centroid", "own")}, 0.0
    with open(path) as f:
        for line in f:
            try:
                d = json.loads(line)
            except ValueError:
                continue  # a line cut by the kill
            if "meta" in d:
                meta = d["meta"]
                continue
            for k in cols:
                cols[k].append(d[k])
            evals = d["evalSeconds"]
    return meta, cols, evals


def spawn(b, m, args):
    tag = os.path.join(args.out, "%s__%s" % (b, m))
    cmd = [sys.executable, os.path.abspath(__file__), "--worker", "--benchmark", b, "--method", m,
           "--out", args.out, "--max-seconds", repr(args.max_seconds), "--max-iter", str(args.max_iter),
           "--max-iter-slsqp", str(args.max_iter_slsqp), "--max-iter-ours", str(args.max_iter_ours), "--grace", repr(args.grace)]
    for p in (tag + ".json", tag + ".history.jsonl"):
        if os.path.exists(p):
            os.remove(p)
    load = os.getloadavg()[0]
    print("== %s %s (cap %.0f s, load %.2f)" % (b, m, args.max_seconds, load), flush=True)
    t0 = time.time()
    killed = False
    with open(tag + ".log", "w") as log:
        p = subprocess.Popen(cmd, stdout=log, stderr=subprocess.STDOUT)
        while p.poll() is None:
            time.sleep(0.5)
            last = _last_line(tag + ".history.jsonl")
            evals = last.get("evalSeconds", 0.0) if last else 0.0
            # Fallback of the child's own alarm (Run.arm).
            if not killed and time.time() - t0 - evals > args.max_seconds + 2 * args.grace + 60:
                try:
                    p.kill()
                except PermissionError:  # a sandbox without signals: wait for the child's alarm
                    pass
                killed = True
    killed = killed or p.returncode == -signal.SIGALRM
    if os.path.exists(tag + ".json"):
        with open(tag + ".json") as f:
            out = json.load(f)
    else:
        # Killed at the cap inside an iteration, or crashed: the history written so far.
        meta, cols, evals = (_from_jsonl(tag + ".history.jsonl") if os.path.exists(tag + ".history.jsonl")
                             else (None, {k: [] for k in ("iteration", "solves", "seconds", "volume", "centroid",
                                                          "own")}, 0.0))
        C = cases()
        case, s = settings(b, m, C)
        wall = time.time() - t0
        if killed:
            summary = dict(stop="time cap", status="killed inside an iteration at %.0f s (cap %.0f s + grace %.0f s)" % (
                wall - evals, args.max_seconds, args.grace), seconds=wall - evals)
        else:
            summary = dict(stop="error: exit %d" % p.returncode, status="exit %d (see %s.log)" % (p.returncode, tag),
                           seconds=wall - evals)
        out = finish(b, m, case.n, s, cols, summary, wall, evals)
        with open(tag + ".json", "w") as f:
            json.dump(out, f)
        if os.path.exists(tag + ".history.jsonl"):
            os.remove(tag + ".history.jsonl")
        print("%-8s %-8s %-14s %7s it %8.3f s  residuals %s / %s  (%s)" % (
            b, m, out["stop"][:14], out.get("iterations"), out["seconds"], out.get("volumeResidual"),
            out.get("centroidResidual"), out["status"]), flush=True)
    out["load"] = load
    with open(tag + ".json", "w") as f:
        json.dump(out, f)
    with open(tag + ".log") as f:
        sys.stdout.write(f.read())
    return out


def load_results(d):
    rs = []
    for b in BENCHMARKS:
        for m in METHODS:
            p = os.path.join(d, "%s__%s.json" % (b, m))
            if os.path.exists(p):
                with open(p) as f:
                    rs.append(json.load(f))
    return rs


def _tex_num(x):
    if x is None or not math.isfinite(x):
        return "---"
    if x == 0:
        return "$0$"
    e = int(math.floor(math.log10(abs(x))))
    m = x / 10 ** e
    if round(m, 1) >= 10:
        m, e = m / 10, e + 1
    if e == 0:
        return "$%.1f$" % m
    return "$%.1f{\\times}10^{%d}$" % (m, e)


def _tex_time(s):
    if s is None:
        return "---"
    return "%#.3g" % s if s < 1000 else "%.0f" % s


TEX_METHOD = {"ours": "our solver", "slsqp": "SLSQP (2025)", "fw": "Frank--Wolfe (2026, $t = %g$)",
              "subgrad": "subgradient (2026, $t = %g$)"}
TEX_STOP = {"rule 3": "\\eqref{eq:results-termination}", "time cap": "time cap", "iteration cap": "iteration cap"}


def tex_method(r):
    m = TEX_METHOD[r["method"]]
    return m % r["settings"]["t"] if "%" in m else m


def write_reports(d):
    rs = load_results(d)
    rows = []
    for r in rs:
        own = r.get("ownRule") or {}
        rows.append(dict(benchmark=r["benchmark"], method=r["method"], n=r["n"], stop=r["stop"], status=r["status"],
                         iterations=r.get("iterations"), solves=r.get("solves"), seconds=r.get("seconds"),
                         volumeResidual=r.get("volumeResidual"), centroidResidual=r.get("centroidResidual"),
                         rule3Met=r.get("rule3Met"), ownRule=own.get("rule"), ownRuleIteration=own.get("iteration"),
                         ownRuleSeconds=own.get("seconds"), ownRuleVolumeResidual=own.get("volumeResidual"),
                         ownRuleCentroidResidual=own.get("centroidResidual"), wallSeconds=r.get("wallSeconds"),
                         evalSeconds=r.get("evalSeconds"), load=r.get("load")))
    with open(os.path.join(d, "results.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, FIELDS)
        w.writeheader()
        w.writerows(rows)

    lines = ["%% benchmarks/bourne/run_table1.py: rows of tab:realizable, generated %s on %s" % (
                 time.strftime("%Y-%m-%d %H:%M"), platform.node()),
             "% columns: case & method & stopped by & iterations & time (s) & volume & centroid",
             "% time: wall seconds, excluding the evaluations of rule 3; residuals of rule 3 at the last iterate"]
    for r in rs:
        stop = r["stop"]
        stop = TEX_STOP.get(stop, "failure" if stop.startswith("error") else stop)
        lines.append("    %-8s & %-30s & %-30s & %s & %s & %s & %s \\\\" % (
            LABEL[r["benchmark"]], tex_method(r), stop, r.get("iterations", "---"), _tex_time(r.get("seconds")),
            _tex_num(r.get("volumeResidual")), _tex_num(r.get("centroidResidual"))))
        own = r.get("ownRule")
        if r["method"] == "slsqp" and own:
            if own["met"]:
                lines.append("    %-8s & %-30s & %-30s & %s & %s & %s & %s \\\\" % (
                    LABEL[r["benchmark"]], tex_method(r), "rule 1", own["iteration"], _tex_time(own["seconds"]),
                    _tex_num(own["volumeResidual"]), _tex_num(own["centroidResidual"])))
            else:
                lines.append("    %% %s SLSQP: rule 1 not met within the run" % LABEL[r["benchmark"]])
        elif own:
            lines.append("    %% %s %s: their rule (%s) %s" % (
                LABEL[r["benchmark"]], r["method"], own["rule"],
                "at iteration %d, %.3g s, residuals %.1e / %.1e" % (own["iteration"], own["seconds"],
                                                                    own["volumeResidual"], own["centroidResidual"])
                if own["met"] else "not met within the run"))
    with open(os.path.join(d, "table.tex"), "w") as f:
        f.write("\n".join(lines) + "\n")
    for r in rows:
        print("%-8s %-8s %-16s %8s it %8s solves %10.3f s  vol %.2e  cen %.2e  own rule %s" % (
            r["benchmark"], r["method"], r["stop"][:16], r["iterations"], r["solves"], r["seconds"] or float("nan"),
            r["volumeResidual"] if r["volumeResidual"] is not None else float("nan"),
            r["centroidResidual"] if r["centroidResidual"] is not None else float("nan"),
            "-" if r["ownRule"] is None else ("not met" if r["ownRuleIteration"] is None else
                                              "it %s" % r["ownRuleIteration"])))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--benchmark", nargs="+", choices=BENCHMARKS, default=BENCHMARKS)
    ap.add_argument("--method", nargs="+", choices=METHODS, default=METHODS)
    ap.add_argument("--max-seconds", type=float, default=1800.0, help="time cap of their methods (charged time)")
    ap.add_argument("--grace", type=float, default=None,
                    help="a child still running at cap + grace (charged) is killed (default max(60, 5%% of the cap))")
    ap.add_argument("--max-iter", type=int, default=10000000, help="iteration cap of FW and the subgradient")
    ap.add_argument("--max-iter-slsqp", type=int, default=100000)
    ap.add_argument("--max-iter-ours", type=int, default=100, help="Newton iterations per attempt")
    ap.add_argument("--out", default=None, help="output directory (default %s, or /tmp/table1-quick)" % OUT)
    ap.add_argument("--quick", action="store_true", help="smoke test: 3 s caps, 2000 iterations, /tmp/table1-quick")
    ap.add_argument("--report", action="store_true", help="only rewrite results.csv and table.tex from --out")
    ap.add_argument("--worker", action="store_true", help=argparse.SUPPRESS)
    args = ap.parse_args()
    if args.quick:
        args.max_seconds = min(args.max_seconds, 3.0)
        args.max_iter = min(args.max_iter, 2000)
        args.grace = 15.0 if args.grace is None else args.grace
    if args.grace is None:
        args.grace = max(60.0, 0.05 * args.max_seconds)
    args.out = os.path.abspath(args.out or ("/tmp/table1-quick" if args.quick else OUT))
    os.makedirs(args.out, exist_ok=True)
    if args.worker:
        worker(args)
        return
    if not args.report:
        t0 = time.time()
        for b in args.benchmark:
            for m in args.method:
                spawn(b, m, args)
        print("total %.0f s" % (time.time() - t0))
    write_reports(args.out)
    print("wrote", os.path.join(args.out, "results.csv"), "and table.tex")


if __name__ == "__main__":
    main()
