#!/usr/bin/env python3
"""Makes the particles of a Power Particles input centroidal with the BNOT relaxation
(de Goes et al. 2012, "Blue noise through optimal transport").

The particles x_i keep their target volumes V_i (the input's targetVol). The energy is
the quantization energy of the semi-discrete optimal transport,
    E(x) = sum_i int_{cell_i} |y - x_i|^2 dy,
where the cells are the power cells of the weights that give the volumes (pysdot
OptimalTransport). Its gradient is dE/dx_i = 2 V_i (x_i - c_i), with c_i the centroid
of cell i, so E is minimal exactly when the particles are the centroids of their cells.
BNOT minimizes E with L-BFGS-B, which needs far fewer evaluations than the Lloyd
iterations of the centroidal start (linear convergence).

Output: a copy of the input .pp with the <pos> records replaced and the <weight> records set to the
weights of the new positions (solved to 1e-13 in the volumes), so that integrators that start from
the loaded weights (Power Particles, pp2) start from the cells with the prescribed volumes.
Velocities, masses and targetVol are unchanged. --maxiter 0 only rewrites the weights.

Usage (from the repository root):
  external/bourne/.venv/bin/python experiments/pp2d/bnot_centroidal.py
      <folder>/<name>.pp [--out <folder>/<name>_centroidal.pp] [--maxiter 2000]
"""

import argparse
import os
import re
import sys
import time

import numpy as np
from scipy.optimize import minimize
from pysdot import OptimalTransport

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(ROOT, "benchmarks", "bourne"))
import bourne_ref as ref  # noqa: E402

POS = re.compile(r"<pos\s+([-\d.eE+]+)\s+([-\d.eE+]+)\s*/>")
TARGET = re.compile(r"<targetVol\s+([-\d.eE+]+)\s*/>")
WEIGHT = re.compile(r"<weight\s+([-\d.eE+]+)\s*/>")
BOX = (-0.5, -0.5, 0.5, 0.5)  # square.line of the yesweights inputs
DIAM = np.hypot(BOX[2] - BOX[0], BOX[3] - BOX[1])


def energy(x, v, om):
    """E(x) and its gradient 2 V (x - c) at the power cells of the volumes v (OT weights solved here)."""
    ot = OptimalTransport(domain=om, positions=x, masses=v, weights=ref.init_w(x, om), verbosity=0,
                          linear_solver="Scipy")
    ot.adjust_weights()
    E = float(np.sum(np.asarray(ot.pd.second_order_moments())))
    c = np.asarray(ot.get_centroids())
    V = np.asarray(ot.pd.integrals())
    return E, (2 * V[:, None] * (x - c)).ravel(), c, V


def weights(x, v, om):
    """Weights of the power cells of x with the volumes v, solved tightly."""
    ot = OptimalTransport(domain=om, positions=x, masses=v, weights=ref.init_w(x, om), verbosity=0,
                          linear_solver="Scipy")
    ot.set_stopping_criterion(1e-13)
    ot.adjust_weights()
    return np.asarray(ot.get_weights()), float(np.max(np.abs(np.asarray(ot.pd.integrals()) - v) / v))


def residual(x, v, om):
    """Largest centroid miss max_i |c_i - x_i| / diameter of the domain."""
    _, _, c, _ = energy(x, v, om)
    return float(np.max(np.linalg.norm(c - x, axis=1))) / DIAM


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input")
    ap.add_argument("--out", default=None)
    ap.add_argument("--maxiter", type=int, default=2000)
    args = ap.parse_args()
    out = args.out or os.path.splitext(args.input)[0] + "_centroidal.pp"

    text = open(args.input).read()
    x0 = np.array(POS.findall(text), float)
    v = np.array(TARGET.findall(text), float)
    if len(v) != len(x0) or abs(v.sum() - 1.0) > 1e-4:
        sys.exit("unexpected input: %d positions, %d targetVol, sum %g" % (len(x0), len(v), v.sum()))
    v = v / v.sum()  # the targetVol values are written to 6 digits
    om = ref.make_box(*BOX)

    res0 = residual(x0, v, om)
    print("input: %d particles, centroid residual %.3e of the diagonal" % (len(v), res0), flush=True)

    t0 = time.perf_counter()

    def f(z):
        E, g, _, _ = energy(z.reshape(-1, 2), v, om)
        return E, g

    x = x0
    if args.maxiter > 0:
        sol = minimize(f, x0.ravel(), jac=True, method="L-BFGS-B",
                       options={"maxiter": args.maxiter, "maxfun": 10 * args.maxiter, "gtol": 1e-14, "ftol": 1e-16})
        x = sol.x.reshape(-1, 2)
        res = residual(x, v, om)
        disp = float(np.max(np.linalg.norm(x - x0, axis=1))) / DIAM
        print("BNOT: %s, %d iterations, %d evaluations, %.1f s" % (sol.message, sol.nit, sol.nfev,
                                                                 time.perf_counter() - t0), flush=True)
        print("output: centroid residual %.3e, max displacement %.3e of the diagonal" % (res, disp), flush=True)
    w, verr = weights(x, v, om)
    print("weights: max relative volume error %.1e" % verr, flush=True)

    it = iter(x.tolist())
    text = POS.sub(lambda m: "<pos %.10g %.10g />" % tuple(next(it)), text)
    iw = iter(w.tolist())
    text = WEIGHT.sub(lambda m: "<weight %.17g />" % next(iw), text)
    with open(out, "w") as f_out:
        f_out.write(text)
    print("wrote", out, flush=True)


if __name__ == "__main__":
    main()
