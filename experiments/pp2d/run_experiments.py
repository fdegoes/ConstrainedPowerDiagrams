#!/usr/bin/env python3
"""Runs the three four-vortex examples (yesweights, Power Particles inputs made centroidal with
bnot_centroidal.py) to the final time T = 1, at two time steps, with:
  pp           Power Particles (cap_tol and cap_max_iters of param.txt);
  gm           Gallouet & Merigot (2018), stiffness kappa = 1/eps^2 with eps = 3 dt;
  softened     our relaxed integrator, at two stiffnesses alpha (ALPHAS);
  constrained  our constrained integrator (SHAKE).
All integrators start from the same centroidal particles and run the same number of steps T / dt.

Time steps: DTS = (0.01, 0.005), the 0.01 of the Power Particles inputs (param.txt) and half of it.

GM's eps: their experiments (Gallouet & Merigot 2018, Sec. 5) use tau / eps between 0.2 and 0.5
(Beltrami N = 900, tau = 0.02, eps = 0.1; Kelvin-Helmholtz N = 2e5, tau = 0.002, eps = 0.005;
Rayleigh-Taylor N = 5e4, tau = 0.001, eps = 0.002), with eps between 1 and 3 particle spacings h.
eps = 3 dt (tau / eps = 1/3) gives eps = 0.03 (2.1 h) and 0.015 (1.06 h) for h = 0.014 (n = 5000),
within both ranges. Their convergence theorem asks for tau = o(eps^2), which they report to be
pessimistic in practice.

Relaxed alpha: the cells of the relaxed scheme are those of the shifted generators q + u with
u = alpha (q - c), so alpha is bounded by the distinctness of these generators, not by the time
step (eq:softened-time-step allows alpha up to 4 / dt^2 - 1). ALPHAS = (0.028, 0.003):
  0.028  0.99 times the practical bound (eq:regularized-practical) at T = 1 along the run on the
         uniform case, the largest fixed alpha that runs (0.1 fails at step 100, 0.5 at 70, 1 at 55
         for dt = 0.01; 9999 = 1/dt^2 - 1, GM's stiffness at eps = dt, fails within 2 steps);
  0.003  a decade lower, towards the limit alpha -> 0 (GM with kappa = 1).
Both run to T = 1 at both time steps on the uniform case.

Each run writes <out>/dt<dt>/<case>/<run>/ (frame_NNNN.csv, energy.csv, run.txt), as runPP2D does,
with <run> one of pp, gm, softened-a<alpha>, constrained. The parameter file of each time step is
written to <out>/dt<dt>/param.txt.

Usage (from the repository root; needs build/runPP2D):
  external/bourne/.venv/bin/python experiments/pp2d/run_experiments.py
      [--dt 0.01 0.005] [--case uniform band2 radial] [--run pp gm softened constrained]
      [--alpha 0.028 0.003] [--final-time 1] [--out experiments/pp2d/runs] [--no-csv]
"""

import argparse
import os
import subprocess
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
DATA = os.path.join(HERE, "yesweights")
BIN = os.path.join(ROOT, "build", "runPP2D")
CASES = {
    "uniform": "uniform_5k_4vort_centroidal.pp",
    "band2": "band2_5k_4vort_centroidal.pp",
    "radial": "radial_5k_4vort_centroidal.pp",
}
RUNS = ("pp", "gm", "softened", "constrained")
DTS = (0.01, 0.005)
ALPHAS = (0.028, 0.003)
GM_EPS_PER_DT = 3.0
# Volume and shift tolerance of the Newton solve of every step (GM, relaxed, constrained; PP uses
# cap_tol and cap_max_iters). The integrators' default, 1e-4, leaves the relaxed volume errors near
# 1e-4; 1e-8 makes them comparable with GM's for about 45 % more Newton iterations.
STEP_TOL = 1e-8


def write_param(dt, path):
    """param.txt of the inputs with the time step dt."""
    lines = [("dt: %r\n" % dt) if line.startswith("dt:") else line for line in open(os.path.join(DATA, "param.txt"))]
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        f.writelines(lines)


def commands(dt, alphas):
    """(run name, integrator options) of one time step."""
    eps = GM_EPS_PER_DT * dt
    out = [("pp", ["--integrator", "powerparticles"]),
           ("gm", ["--integrator", "gm", "--stiffness", repr(1.0 / eps ** 2)])]
    for a in alphas:
        if not dt < 2.0 / np.sqrt(1.0 + a):  # eq:softened-time-step
            sys.exit("dt %g outside the time-step bound for alpha %g" % (dt, a))
        out.append(("softened-a%g" % a, ["--integrator", "softened", "--alpha", repr(a)]))
    out.append(("constrained", ["--integrator", "constrained"]))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dt", nargs="+", type=float, default=list(DTS))
    ap.add_argument("--case", nargs="+", choices=list(CASES), default=list(CASES))
    ap.add_argument("--run", nargs="+", choices=list(RUNS), default=list(RUNS))
    ap.add_argument("--alpha", nargs="+", type=float, default=list(ALPHAS), help="stiffnesses of the relaxed runs")
    ap.add_argument("--final-time", type=float, default=1.0)
    ap.add_argument("--out", default=os.path.join(HERE, "runs"))
    ap.add_argument("--no-csv", action="store_true")
    args = ap.parse_args()
    if not os.path.exists(BIN):
        sys.exit("missing " + BIN + " (cmake --build build --target runPP2D)")

    failed = []
    for dt in args.dt:
        steps = int(round(args.final_time / dt))
        root = os.path.join(args.out, "dt%g" % dt)
        param = os.path.join(root, "param.txt")
        write_param(dt, param)
        print("dt %g: %d steps, GM eps %g (kappa %.4g)" % (dt, steps, GM_EPS_PER_DT * dt,
                                                           1 / (GM_EPS_PER_DT * dt) ** 2), flush=True)
        for case in args.case:
            for name, opts in commands(dt, args.alpha):
                if name.split("-")[0] not in args.run:
                    continue
                cmd = [BIN, os.path.join(DATA, "square.line"), os.path.join(DATA, CASES[case]), param,
                       os.path.join(root, case, name), "--frames", str(steps), "--lloydtol", "1e-6",
                       "--volumetol", repr(STEP_TOL)] + opts
                if args.no_csv:
                    cmd.append("--no-csv")
                t0 = time.perf_counter()
                r = subprocess.run(cmd, capture_output=True, text=True)
                ok = r.returncode == 0
                last = ((r.stdout if ok else r.stderr).strip().splitlines() or ["no message"])[-1]
                if not ok:
                    failed.append("dt%g/%s/%s" % (dt, case, name))
                print("  %-8s %-17s %-7s %7.1f s  %s" % (case, name, "ok" if ok else "FAILED",
                                                         time.perf_counter() - t0, last), flush=True)
    if failed:
        print("failed:", ", ".join(failed), flush=True)


if __name__ == "__main__":
    main()
