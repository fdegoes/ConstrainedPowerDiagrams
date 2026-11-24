#!/usr/bin/env python3
"""Runs the three four-vortex examples (yesweights, Power Particles inputs made centroidal with
bnot_centroidal.py) with Power Particles (PP), the Gallouet-Merigot scheme (GM), our constrained
integrator and our relaxed (softened) integrator.

Inputs: experiments/pp2d/yesweights/<case>_centroidal.pp with the parameters of param.txt (dt 0.01,
cap_tol 1e-6, cap_max_iters 5); all four integrators start from the same centroidal particles.

Stiffness of the relaxed integrator: alpha = 1 (--alpha). It is checked against the two bounds of
the paper at the start:
  - time step (eq:softened-time-step): dt < 2 / sqrt(1 + alpha);
  - distinct generators: dt (|v_i| + |v_j|) must stay below |q_i - q_j| for every pair of nearest
    neighbours (the residual of a step is dt |v| / alpha, so the bound does not depend on alpha).
alpha = 1/dt^2 - 1 (the full Lloyd step of the paper, (1 + alpha) dt^2 = 1) was tried first: the
Newton solve of the relaxed step fails at the cap of 100 iterations within 2 frames. alpha = 1 runs
20 frames without failure; alpha = 10 fails at frame 18, alpha = 30 at frame 10.
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


def check_bounds(path, dt, alpha):
    """Time-step bound and the distinctness check of the first step."""
    text = open(path).read()
    q = np.array(POS.findall(text), float)
    v = np.array(VEL.findall(text), float)
    stable = dt < 2.0 / np.sqrt(1.0 + alpha)
    speed = np.linalg.norm(v, axis=1)
    # Pairwise distance over the speeds of the pair, only the nearest neighbours matter for the check:
    # for each particle, its 8 nearest others.
    worst = np.inf
    from scipy.spatial import cKDTree
    tree = cKDTree(q)
    dist, idx = tree.query(q, k=9)
    for i in range(len(q)):
        for d, j in zip(dist[i, 1:], idx[i, 1:]):
            worst = min(worst, d / (dt * (speed[i] + speed[j]) + 1e-300))
    return stable, worst


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--case", nargs="+", choices=list(CASES), default=list(CASES))
    ap.add_argument("--integrator", nargs="+", choices=list(INTEGRATORS), default=list(INTEGRATORS))
    ap.add_argument("--frames", type=int, default=100)
    ap.add_argument("--out", default=os.path.join(HERE, "runs"))
    ap.add_argument("--no-csv", action="store_true")
    ap.add_argument("--alpha", type=float, default=1.0, help="stiffness of the relaxed integrator")
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
        if integrator == "softened":
            stable, worst = check_bounds(particles, dt, args.alpha)
            print("%-8s softened: time step %s the bound; distinctness ratio (>1 is safe) %.3g"
                  % (case, "within" if stable else "OUTSIDE", worst), flush=True)
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
