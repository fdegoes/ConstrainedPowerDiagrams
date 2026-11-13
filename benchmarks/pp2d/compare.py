#!/usr/bin/env python3
"""Compares a runPP2D run with a run of the original Power Particles code
(pp2, `pp2 -r -i <in.pp> -d <domain.line> -p <param> -n N -f 1 -o <dir>`),
frame by frame, particle by particle (same order in both).

Reported per frame:
  dq      max and RMS |q_ours - q_pp2|, in mean particle spacings
          h = sqrt(|domain| / n)
  dv      max and RMS |v_ours - v_pp2| relative to the RMS speed of pp2
  K       kinetic energy sum_i m_i |v_i|^2 / 2 of both (pp2's energy.csv)
  dV      max relative cell-volume error of both

Writes <out>/compare.csv and, with --plot, compare.png (curves over time),
and with --frames, side-by-side power diagrams of both (render.py style).

Usage:
  external/bourne/.venv/bin/python benchmarks/pp2d/compare.py <ours_dir> <pp2_dir>
      [--out DIR] [--plot] [--frames 0 100 300] [--color band|radius]
"""

import argparse
import csv
import glob
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)


def read_pp(path):
    """pos, vel, mass, targetVol of the fluid particles of a .pp frame, and the
    volume of each cell (from its DualCell facets, about the site)."""
    pos, vel, mass, tvol, cellvol = [], [], [], [], []
    fluid = False
    area = None
    q = None
    for line in open(path):
        line = line.strip()
        if line == "<Fluid>":
            fluid = True
        elif line in ("<Air>", "<Solid>"):
            fluid = False
        elif line == "</End>":
            fluid = False
        if not fluid:
            continue
        if line.startswith("<pos "):
            q = [float(x) for x in line[5:-2].split()]
            pos.append(q)
        elif line.startswith("<vel "):
            vel.append([float(x) for x in line[5:-2].split()])
        elif line.startswith("<mass "):
            mass.append(float(line[6:-2]))
        elif line.startswith("<targetVol "):
            tvol.append(float(line[11:-2]))
        elif line == "<DualCell>":
            area = 0.0
        elif line.startswith("<segment ") and area is not None:
            x0, y0, x1, y1 = (float(x) for x in line[9:-2].split())
            # Signed area of the triangle (q, s0, s1); facets are CCW.
            area += 0.5 * ((x0 - q[0]) * (y1 - q[1]) - (x1 - q[0]) * (y0 - q[1]))
        elif line == "</DualCell>":
            cellvol.append(abs(area))
            area = None
    return np.array(pos), np.array(vel), np.array(mass), np.array(tvol), np.array(cellvol)


def read_ours(path):
    rows = list(csv.DictReader(open(path)))
    f = lambda k: np.array([float(r[k]) for r in rows])
    return np.c_[f("x"), f("y")], np.c_[f("vx"), f("vy")], f("volume"), f("cellVolume")


def frame_index(paths, pattern):
    out = {}
    for p in paths:
        m = re.search(pattern, os.path.basename(p))
        if m:
            out[int(m.group(1))] = p
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("ours")
    ap.add_argument("pp2")
    ap.add_argument("--out", default=None)
    ap.add_argument("--plot", action="store_true")
    ap.add_argument("--frames", nargs="*", type=int, default=None, help="render these frames side by side")
    ap.add_argument("--color", choices=["band", "radius"], default=None)
    args = ap.parse_args()
    out = args.out or os.path.join(args.ours, "compare")
    os.makedirs(out, exist_ok=True)

    ours = frame_index(glob.glob(os.path.join(args.ours, "frame_*.csv")), r"frame_(\d+)\.csv$")
    theirs = frame_index(glob.glob(os.path.join(args.pp2, "frame*.pp")), r"frame(\d+)\.pp$")
    common = sorted(set(ours) & set(theirs))
    if not common:
        sys.exit("no common frames")

    rows = []
    for k in common:
        q, v, Vbar, V = read_ours(ours[k])
        Q, W, M, T, C = read_pp(theirs[k])
        if len(Q) != len(q):
            sys.exit("frame %d: %d particles in ours, %d in pp2" % (k, len(q), len(Q)))
        h = np.sqrt(Vbar.sum() / len(q))
        dq = np.linalg.norm(q - Q, axis=1) / h
        vrms = np.sqrt((W ** 2).sum(1).mean()) or 1.0
        dv = np.linalg.norm(v - W, axis=1) / vrms
        r = dict(frame=k, dq_max=dq.max(), dq_rms=np.sqrt((dq ** 2).mean()), dv_max=dv.max(),
                 dv_rms=np.sqrt((dv ** 2).mean()),
                 K_ours=0.5 * (Vbar * (v ** 2).sum(1)).sum(), K_pp2=0.5 * (M * (W ** 2).sum(1)).sum(),
                 dV_ours=(np.abs(V - Vbar) / Vbar).max() if np.isfinite(V).all() else np.nan,
                 dV_pp2=(np.abs(C - T) / T).max() if len(C) == len(T) else np.nan)
        rows.append(r)
    with open(os.path.join(out, "compare.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    for r in rows:
        if r["frame"] % max(1, len(rows) // 10) == 0 or r["frame"] == rows[-1]["frame"]:
            print("frame %4d  dq/h max %.2e rms %.2e  dv/vrms max %.2e rms %.2e  K ours %.6g pp2 %.6g (rel %.1e)  "
                  "dV ours %.1e pp2 %.1e" % (r["frame"], r["dq_max"], r["dq_rms"], r["dv_max"], r["dv_rms"], r["K_ours"],
                                            r["K_pp2"], abs(r["K_ours"] - r["K_pp2"]) / r["K_pp2"], r["dV_ours"], r["dV_pp2"]))

    if args.plot:
        os.environ.setdefault("MPLCONFIGDIR", "/tmp/mplcfg")
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        t = np.array([r["frame"] for r in rows])
        fig, ax = plt.subplots(1, 3, figsize=(15, 4))
        ax[0].semilogy(t, [r["dq_rms"] for r in rows], label="rms")
        ax[0].semilogy(t, [r["dq_max"] for r in rows], label="max")
        ax[0].set_title("|q_ours - q_pp2| / h"); ax[0].legend()
        ax[1].plot(t, [r["K_ours"] for r in rows], label="ours")
        ax[1].plot(t, [r["K_pp2"] for r in rows], "--", label="pp2")
        ax[1].set_title("kinetic energy"); ax[1].legend()
        ax[2].semilogy(t, [r["dV_ours"] for r in rows], label="ours")
        ax[2].semilogy(t, [r["dV_pp2"] for r in rows], "--", label="pp2")
        ax[2].set_title("max relative volume error"); ax[2].legend()
        for a in ax:
            a.set_xlabel("frame")
        fig.tight_layout()
        fig.savefig(os.path.join(out, "compare.png"), dpi=100)
        print("wrote", os.path.join(out, "compare.png"))

    if args.frames:
        import render as R
        os.environ.setdefault("MPLCONFIGDIR", "/tmp/mplcfg")
        import matplotlib.pyplot as plt
        q0, _, _ = R.read_frame(ours[common[0]])
        mode = args.color or ("band" if "band" in os.path.abspath(args.ours) else "radius")
        colors = R.colors_for(q0, mode)
        for k in args.frames:
            if k not in ours or k not in theirs:
                continue
            q, u, w = R.read_frame(ours[k])
            R.draw(os.path.join(out, "ours_%04d.png" % k), R.power_cells(q + u, w + (u * u).sum(1)), q, colors, False, 6)
            # pp2's cells, from its DualCell segments of each particle.
            draw_pp(theirs[k], os.path.join(out, "pp2_%04d.png" % k), colors, R)
            print("frame %d: wrote ours_%04d.png and pp2_%04d.png" % (k, k, k))


def draw_pp(path, png, colors, R):
    """Draws the cells of a .pp frame (polygons chained from the segments)."""
    cells, segs, inside = [], [], False
    for line in open(path):
        line = line.strip()
        if line == "<Fluid>":
            inside = True
        elif line in ("<Air>", "<Solid>", "</End>"):
            inside = False
        if not inside:
            continue
        if line == "<DualCell>":
            segs = []
        elif line.startswith("<segment "):
            x0, y0, x1, y1 = (float(x) for x in line[9:-2].split())
            segs.append(((x0, y0), (x1, y1)))
        elif line == "</DualCell>":
            cells.append(_chain(segs))
    R.draw(png, cells, np.zeros((len(cells), 2)), colors, False, 6)


def _chain(segs):
    if not segs:
        return np.zeros((0, 2))
    nxt = {}
    for a, b in segs:
        nxt[(round(a[0], 12), round(a[1], 12))] = b
    start = segs[0][0]
    poly = [start]
    cur = segs[0][1]
    for _ in range(len(segs)):
        if abs(cur[0] - start[0]) < 1e-12 and abs(cur[1] - start[1]) < 1e-12:
            break
        poly.append(cur)
        cur = nxt.get((round(cur[0], 12), round(cur[1], 12)), start)
    return np.array(poly)


if __name__ == "__main__":
    main()
