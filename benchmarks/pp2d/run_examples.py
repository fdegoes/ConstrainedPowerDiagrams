#!/usr/bin/env python3
"""Runs the fluid-only 2D examples of the original Power Particles code
(Before2015/pp/data/2D) with runPP2D and one of our integrators.

The cases are those whose fluid fills the domain (no air, no solids):
  yesweights/{uniform,band2,radial}_5k, band2_5k_4vort
  param: dt 0.01, cap_tol 1e-6, cap_max_iters 5
(5000 particles in [-0.5, 0.5]^2, no gravity, no viscosity). The noweights
cases are skipped: they are the same particles with cap_max_iters 0 (Power
Particles without the weight solve), which none of our schemes reproduces.

Each run writes <out>/<case>/<integrator>/ with frame_NNNN.csv (state per
frame: x, y, vx, vy, ux, uy, w, p, volume, cx, cy, cellVolume; not with
--no-csv), energy.csv (energies, momenta, substeps, iterations, errors per
frame; not with --no-energy) and run.txt. render.py needs the frame CSVs.

Usage:
  python3 benchmarks/pp2d/run_examples.py [--case CASE ...]
      [--integrator powerparticles softened gm constrained] [--frames 100]
      [--every 1] [--alpha 0.01] [--stiffness 1] [--lloydtol 1e-6] [--lloyditers N]
      [--data DIR] [--out DIR]
      [--no-csv] [--no-energy] [--list]
"""

import argparse
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
RUN = os.path.join(ROOT, "build", "runPP2D")
DATA = os.path.expanduser("~/Documents/Projects/Before2015/pp/data/2D")

# name: (directory, particles file)
CASES = {
    "yesweights-uniform": ("yesweights", "uniform_5k.pp"),
    "yesweights-band2": ("yesweights", "band2_5k.pp"),
    "yesweights-band2-4vort": ("yesweights", "band2_5k_4vort.pp"),
    "yesweights-radial": ("yesweights", "radial_5k.pp"),
}
# softened and constrained start from a centroidal diagram (Lloyd), which does
# not reach the default 1e-9 on 5000 particles, hence the looser --lloydtol.
# The constrained integrator does not yet run these cases (with that start its
# first step reaches the iteration cap): opt-in only.
INTEGRATORS = ("powerparticles", "softened", "gm", "constrained")
DEFAULT_INTEGRATORS = ("powerparticles", "softened", "gm")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--case", nargs="+", choices=list(CASES), default=list(CASES))
    ap.add_argument("--integrator", nargs="+", choices=INTEGRATORS, default=list(DEFAULT_INTEGRATORS))
    ap.add_argument("--frames", type=int, default=100)
    ap.add_argument("--every", type=int, default=1, help="advances of dt per saved frame")
    ap.add_argument("--alpha", type=float, default=0.01, help="softened integrator stiffness")
    ap.add_argument("--stiffness", type=float, default=1.0, help="kappa of GM (1 + alpha matches the softened integrator)")
    ap.add_argument("--lloydtol", type=float, default=1e-6, help="tolerance of the centroidal start (softened, constrained)")
    ap.add_argument("--lloyditers", type=int, default=None, help="iteration cap of the centroidal start")
    ap.add_argument("--data", default=DATA)
    ap.add_argument("--out", default=os.path.join(HERE, "results"))
    ap.add_argument("--no-csv", action="store_true", help="do not write the per-frame state CSVs")
    ap.add_argument("--no-energy", action="store_true", help="do not write energy.csv")
    ap.add_argument("--list", action="store_true")
    args = ap.parse_args()
    if args.list:
        for name, (d, f) in CASES.items():
            print("%-24s %s/%s" % (name, d, f))
        return
    if not os.path.exists(RUN):
        sys.exit("build runPP2D first: cmake --build build --target runPP2D")

    failed = []
    for name in args.case:
        d, f = CASES[name]
        src = os.path.join(args.data, d)
        for integ in args.integrator:
            out = os.path.join(args.out, name, integ)
            cmd = [RUN, os.path.join(src, "square.line"), os.path.join(src, f), os.path.join(src, "param.txt"), out,
                   "--integrator", integ, "--frames", str(args.frames), "--every", str(args.every),
                   "--alpha", repr(args.alpha), "--stiffness", repr(args.stiffness)]
            cmd += ["--lloydtol", repr(args.lloydtol)]
            if args.lloyditers is not None:
                cmd += ["--lloyditers", str(args.lloyditers)]
            cmd += ["--no-csv"] * args.no_csv + ["--no-energy"] * args.no_energy
            t = time.time()
            r = subprocess.run(cmd, capture_output=True, text=True)
            last = (r.stdout.strip().splitlines() or [""])[-1]
            status = "ok" if r.returncode == 0 else "FAILED (%s)" % (r.stderr.strip().splitlines() or ["?"])[-1]
            print("%-24s %-15s %6.1fs  %s  %s" % (name, integ, time.time() - t, status, last), flush=True)
            if r.returncode != 0:
                failed.append((name, integ))
    if failed:
        sys.exit("%d runs failed" % len(failed))


if __name__ == "__main__":
    main()
