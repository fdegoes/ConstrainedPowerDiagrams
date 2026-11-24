#!/usr/bin/env python3
"""Runs the three four-vortex examples (yesweights, Power Particles inputs made centroidal with
bnot_centroidal.py) with Power Particles (PP), the Gallouet-Merigot scheme (GM), our constrained
integrator and our relaxed (softened) integrator.

Inputs: experiments/pp2d/yesweights/<case>_centroidal.pp with the parameters of param.txt (dt 0.01,
cap_tol 1e-6, cap_max_iters 5); all four integrators start from the same centroidal particles.

Stiffness of the relaxed integrator: alpha = 0.028 (--alpha). It is the rule of the paper, 0.99 times
the practical distinctness bound (eq:regularized-practical),
    alpha < min_{i != j} |q_i - q_j| / (|c_i - q_i| + |c_j - q_j|),
evaluated along the run: the bound falls from about 2 at frame 20 to 0.029 at frame 100 on the uniform
case (the centroid gap grows with the flow, for GM as well). Fixed values tried on the uniform case:
alpha = 0.01 runs (bound stays above 0.029), 0.1 fails at frame 100, 0.5 at frame 70, 1 at frame 55,
and the Lloyd-step value 1/dt^2 - 1 = 9999 fails within 2 frames. The time-step bound
(eq:softened-time-step), dt < 2 / sqrt(1 + alpha), holds with a margin of about 200.
GM is the limit alpha -> 0 of the relaxed scheme, with kappa = 1 (--stiffness 1).

Each run writes <out>/<case>/<integrator>/ (frame_NNNN.csv, energy.csv, run.txt), as runPP2D does.

Usage (from the repository root; needs build/runPP2D):
  external/bourne/.venv/bin/python experiments/pp2d/run_experiments.py
      [--case uniform band2 radial] [--integrator powerparticles gm constrained softened]
      [--frames 100] [--out experiments/pp2d/runs] [--no-csv]
"""

import argparse
import itertools
import os
import re
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
INTEGRATORS = ("powerparticles", "gm", "constrained", "softened")
POS = re.compile(r"<pos\s+([-\d.eE+]+)\s+([-\d.eE+]+)\s*/>")
VEL = re.compile(r"<vel\s+([-\d.eE+]+)\s+([-\d.eE+]+)\s*/>")


def dt_of_param(path):
    for line in open(path):
        if line.startswith("dt:"):
            return float(line.split(":")[1])
    raise ValueError("no dt in " + path)


def time_step_ok(dt, alpha):
    """The time-step bound of the relaxed scheme (eq:softened-time-step)."""
    return dt < 2.0 / np.sqrt(1.0 + alpha)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--case", nargs="+", choices=list(CASES), default=list(CASES))
    ap.add_argument("--integrator", nargs="+", choices=list(INTEGRATORS), default=list(INTEGRATORS))
    ap.add_argument("--frames", type=int, default=100)
    ap.add_argument("--out", default=os.path.join(HERE, "runs"))
    ap.add_argument("--no-csv", action="store_true")
    ap.add_argument("--alpha", type=float, default=0.028, help="stiffness of the relaxed integrator")
    args = ap.parse_args()

    param = os.path.join(DATA, "param.txt")
    dt = dt_of_param(param)
    alpha = args.alpha
    print("dt %g, relaxed alpha %g, time-step bound 2/sqrt(1+alpha) = %.4g" % (dt, alpha, 2 / np.sqrt(1 + alpha)),
          flush=True)
    if not os.path.exists(BIN):
        sys.exit("missing " + BIN + " (cmake --build build --target runPP2D)")

    failed = []
    for case, integrator in itertools.product(args.case, args.integrator):
        particles = os.path.join(DATA, CASES[case])
        if integrator == "softened" and not time_step_ok(dt, args.alpha):
            sys.exit("time step %g outside the bound 2/sqrt(1+alpha) for alpha %g" % (dt, args.alpha))
        cmd = [BIN, os.path.join(DATA, "square.line"), particles, param,
               os.path.join(args.out, case, integrator), "--integrator", integrator,
               "--frames", str(args.frames)]
        cmd += ["--lloydtol", "1e-6"]  # the start of the centroidal integrators (default 1e-9 is too strict)
        if integrator == "softened":
            cmd += ["--alpha", repr(args.alpha)]
        if integrator == "gm":
            cmd += ["--stiffness", "1"]
        if args.no_csv:
            cmd += ["--no-csv"]
        t0 = time.perf_counter()
        r = subprocess.run(cmd, capture_output=True, text=True)
        status = "ok" if r.returncode == 0 else "FAILED"
        last = (r.stdout.strip().splitlines() or [""])[-1]
        if r.returncode != 0:
            failed.append((case, integrator))
            last = (r.stderr.strip().splitlines() or ["no message"])[-1]
        print("%-8s %-15s %-7s %7.1f s  %s" % (case, integrator, status, time.perf_counter() - t0, last), flush=True)
    if failed:
        print("failed:", ", ".join("%s/%s" % f for f in failed), flush=True)


if __name__ == "__main__":
    main()
