#!/usr/bin/env python3
"""Cases of the synthetic benchmark: realizable power diagrams in the unit square
(dim = 2) or the unit cube (dim = 3).

Axes:
  n         number of cells
  noise     0 (low noise), 1 (mid noise) or 2 (high noise): the generators Y
            are uniform random points in the square (cube) at pairwise
            distance >= rho h (dart throwing), rho = 0.8, 0.5, 0.25,
            h = (1/n)^(1/dim) the cell size; the smaller rho, the less
            uniform the spacing
  volumes   uniform: all 1 / n
            voronoi: the Voronoi cell volumes of Y (zero weights)
            lognormal50, lognormal100: lognormal volumes with a coefficient of
              variation of 0.5 or 1, at random over the generators (stress set,
              not in the default cases; written to volumes.txt)
  seed      of the generators

The targets are the volumes and centroids of the power diagram of Y with the
prescribed volumes (makeSynthetic), so Y and its weights are the solution.
The solver starts from the sites q = cbar.

The 3D cases go into their own root directory (results3d), with the same case
names. Each case is a directory with
  generators.txt   n (and dim in 3D), then Y (input of makeSynthetic)
  data.txt         n, volumes, targets (cbar) and weights of the solution
  truth.json       Y, weights, volumes, centroids (the known solution)
"""

import itertools
import json
import math
import os
import subprocess

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
MAKE = os.path.join(ROOT, "build", "makeSynthetic")

NOISE = (0, 1, 2)
NOISE_SPACING = (0.8, 0.5, 0.25)  # minimum distance of the generators per noise level, in cell sizes h
VOLUMES = ("uniform", "voronoi")
# Stress set (not in the default 600 cases): lognormal volumes, the usual model of grain sizes,
# with these coefficients of variation, assigned to the generators at random.
LOGNORMAL = {"lognormal50": 0.5, "lognormal100": 1.0}


def dart_throwing(n, rng, spacing=None, max_tries_per_point=1000, dim=2):
    """n uniform random points in the unit square or cube (sequential
    rejection sampling).
    spacing(p): with it, two points p, q are at distance >= (spacing(p) +
    spacing(q)) / 2 (Poisson-disk sampling with a varying radius; a background
    grid of cell max spacing finds the neighbors)."""
    pts = []
    smax = None
    if spacing is not None:
        m = 200 if dim == 2 else 40
        g = (np.arange(m) + 0.5) / m
        smax = float(spacing(np.stack(np.meshgrid(*[g] * dim), -1)).max()) * 1.01
        cell = 1.0 / max(1, int(1.0 / smax))
        offsets = list(itertools.product((-1, 0, 1), repeat=dim))
    grid = {}
    tries = 0
    while len(pts) < n:
        tries += 1
        if tries > max_tries_per_point * n:
            raise RuntimeError("dart throwing: cannot fit %d points" % n)
        p = rng.random(dim)
        if spacing is not None:
            sp = float(spacing(p))
            key = tuple(int(x / cell) for x in p)
            ok = True
            for off in offsets:
                for q, sq in grid.get(tuple(k + o for k, o in zip(key, off)), ()):
                    r = 0.5 * (sp + sq)
                    if ((p - q) ** 2).sum() < r * r:
                        ok = False
                        break
                if not ok:
                    break
            if not ok:
                continue
            grid.setdefault(key, []).append((p, sp))
        pts.append(p)
    return np.array(pts)


def generators(n, noise, seed, dim=2):
    # Seeds per (n, noise): uniform and voronoi volumes share Y. The 3D
    # streams get their own seed.
    rng = np.random.default_rng([seed, n, noise] if dim == 2 else [seed, n, noise, dim])
    r = NOISE_SPACING[noise] * (1.0 / n) ** (1.0 / dim)
    return dart_throwing(n, rng, spacing=lambda p: np.full(np.shape(p)[:-1], r), dim=dim)


def lognormal_volumes(n, volumes, seed, dim=2):
    """n volumes with a lognormal distribution of coefficient of variation LOGNORMAL[volumes]
    (sigma^2 = log(1 + cv^2)), normalized to sum 1 (makeSynthetic normalizes too); their own
    random stream, so that the generators are those of the uniform and Voronoi cases."""
    cv = LOGNORMAL[volumes]
    rng = np.random.default_rng([seed, n, dim, int(round(100 * cv)), 7])
    v = rng.lognormal(0.0, math.sqrt(math.log(1.0 + cv * cv)), n)
    return v / v.sum()


def case_name(n, noise, volumes, seed):
    return "noise%d-%s-%d-%d" % (noise, volumes, n, seed)


def read_data(path):
    rec = {}
    for line in open(path):
        key, *vals = line.split()
        rec[key] = np.array([float(v) for v in vals])
    n = int(rec["n"][0])
    dim = int(rec["dim"][0]) if "dim" in rec else 2
    return {"n": n, "dim": dim, "volumes": rec["volumes"], "targets": rec["targets"].reshape(n, dim),
            "weights": rec["weights"]}


def make_case(root, n, noise, volumes, seed, force=False, dim=2):
    """Writes the case directory (if needed) and returns its data."""
    d = os.path.join(root, case_name(n, noise, volumes, seed))
    data = os.path.join(d, "data.txt")
    if force or not os.path.exists(data):
        os.makedirs(d, exist_ok=True)
        Y = generators(n, noise, seed, dim)
        with open(os.path.join(d, "generators.txt"), "w") as f:
            f.write("%d\n" % n if dim == 2 else "%d %d\n" % (n, dim))
            f.write(" ".join(repr(float(x)) for x in Y.ravel()) + "\n")
        mode = volumes
        if volumes in LOGNORMAL:
            mode = os.path.join(d, "volumes.txt")
            with open(mode, "w") as f:
                f.write("%d\n" % n + " ".join(repr(float(x)) for x in lognormal_volumes(n, volumes, seed, dim)) + "\n")
        subprocess.run([MAKE, os.path.join(d, "generators.txt"), mode, data], check=True,
                       stdout=subprocess.DEVNULL)
        c = read_data(data)
        json.dump({"n": n, "dim": dim, "noise": noise, "volumes": volumes, "seed": seed,
                   "generators": Y.ravel().tolist(), "weights": c["weights"].tolist(),
                   "cellVolumes": c["volumes"].tolist(), "centroids": c["targets"].ravel().tolist()},
                  open(os.path.join(d, "truth.json"), "w"))
    c = read_data(data)
    c["dir"] = d
    c["name"] = os.path.basename(d)
    c["truth"] = np.array(json.load(open(os.path.join(d, "truth.json")))["generators"]).reshape(-1, c["dim"])
    return c
