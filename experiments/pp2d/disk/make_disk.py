#!/usr/bin/env python3
"""Disk test case without corners: a regular polygon of radius R (disk.line), n particles on a
jittered hexagonal grid inside it with equal volumes, and the differential rotation
    v_theta(r) = 4 U (r / R) (1 - r^2 / R^2),
divergence-free and tangent to the boundary, with its peak speed 8 U / (3 sqrt 3) at r = R / sqrt 3.
The particles are not centroidal: the integrators' Initialize() relaxes them (Lloyd).

Usage (from the repository root):
  external/bourne/.venv/bin/python experiments/pp2d/disk/make_disk.py [--n 2000] [--sides 128]
      [--speed 0.1] [--out experiments/pp2d/disk]
"""

import argparse
import os

import numpy as np


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--n", type=int, default=2000)
    ap.add_argument("--sides", type=int, default=128)
    ap.add_argument("--radius", type=float, default=0.5)
    ap.add_argument("--speed", type=float, default=0.1, help="peak speed")
    ap.add_argument("--out", default=os.path.dirname(os.path.abspath(__file__)))
    args = ap.parse_args()
    R, m = args.radius, args.sides
    t = 2 * np.pi * np.arange(m) / m
    poly = np.c_[R * np.cos(t), R * np.sin(t)]
    area = 0.5 * m * R * R * np.sin(2 * np.pi / m)
    with open(os.path.join(args.out, "disk.line"), "w") as f:
        for x, y in poly:
            f.write("v %.17g %.17g\n" % (x, y))
        for i in range(m):
            f.write("e %d %d\n" % (i + 1, (i + 1) % m + 1))
    # Jittered hexagonal grid with about n points inside the polygon's inscribed circle.
    rin = R * np.cos(np.pi / m) * 0.985
    a = np.sqrt(2 * np.pi * rin * rin / (np.sqrt(3) * args.n))
    pts = []
    for j in range(-int(R / a) - 2, int(R / a) + 3):
        for i in range(-int(R / a) - 2, int(R / a) + 3):
            x, y = (i + 0.5 * (j % 2)) * a, j * a * np.sqrt(3) / 2
            if x * x + y * y < rin * rin:
                pts.append((x, y))
    rng = np.random.default_rng(1)
    q = np.array(pts) + rng.uniform(-0.05 * a, 0.05 * a, (len(pts), 2))
    q = q[np.hypot(q[:, 0], q[:, 1]) < rin]
    n = len(q)
    U = args.speed * 3 * np.sqrt(3) / 8
    r = np.hypot(q[:, 0], q[:, 1])
    vt = 4 * U * (r / R) * (1 - (r / R) ** 2)
    v = np.c_[-vt * q[:, 1] / np.maximum(r, 1e-300), vt * q[:, 0] / np.maximum(r, 1e-300)]
    vol = area / n
    with open(os.path.join(args.out, "disk_%d_rotation.pp" % n), "w") as f:
        f.write("<Fluid>\n<size %d />\n<offset 0 />\n" % n)
        for (x, y), (vx, vy) in zip(q, v):
            f.write("<Particle>\n<pos %.17g %.17g />\n<vel %.17g %.17g />\n<mass %.17g />\n<pres 0 />\n"
                    "<weight 0 />\n<targetVol %.17g />\n</Particle>\n" % (x, y, vx, vy, vol, vol))
        f.write("</End>\n")
    print("disk: %d particles, area %.6f, peak speed %.3f, spacing %.4f" % (n, area, np.abs(vt).max(), a))


if __name__ == "__main__":
    main()
