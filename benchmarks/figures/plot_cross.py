#!/usr/bin/env python3
"""Draws the aligned (cross) example of Bourne et al. 2026, n = 728 equal cells in the unit
square (case 2026-aligned of benchmarks/bourne/bourne_cases.py), solved by our H and by their
Frank-Wolfe method (t = 100), side by side:
  latex/figures/cross.pdf
(a) our H (benchBourne, solver H, rule 3 at 1e-8, sites = targets, as in run_realizable.py);
(b) their FW (bourne_ref.convex_order_projection_frank_wolfe from the start of the case,
    gap rule disabled as in run_theirs_common, stopped at --fw-iter iterations or
    --fw-seconds seconds, whichever comes first).
The cells are tinted by sector: the true diagram has straight facet rays along the
diagonals and the medians from an 8-valent centre vertex, which cut it into eight sectors of
91 cells (the reducible diagram of appendix-reducible.tex); each cell keeps the sector of its
cell in panel (a) in both panels. Generators are black dots, targets small open circles, and
in (b) a red segment joins each centroid to its target where the error exceeds --arrow-min
cell sizes (rule 3 at 1e-8 is met by (a) only). FW's generators are drawn as (t - 1) Z,
the generators of the equivalent H_alpha (display_generators). To the right of each diagram,
two zooms: the centre (red) and the inner corner of the cross at the top left (orange), where
the centroid errors of FW are largest. The script also prints the residuals of both runs and
checks the cells against the solvers' volumes and centroids and against the true diagram.

The solutions are cached in benchmarks/figures/data/cross.json (--resolve to recompute).
Needs the built build/benchBourne and the pysdot environment:
  MPLCONFIGDIR=/tmp/mplcfg external/bourne/.venv/bin/python benchmarks/figures/plot_cross.py
      [--resolve] [--fw-iter 5000] [--fw-seconds 900] [--width 6.9]
"""

import argparse
import json
import math
import os
import subprocess
import sys
import tempfile
import time

os.environ.setdefault("MPLCONFIGDIR", "/tmp/mplcfg")
import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.collections import LineCollection, PolyCollection  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(ROOT, "benchmarks", "bourne"))
import bourne_cases as bc  # noqa: E402
import bourne_ref as ref  # noqa: E402

CACHE = os.path.join(HERE, "data", "cross.json")
OUT = os.path.join(ROOT, "latex", "figures", "cross.pdf")
EPS = 1e-12
TOL = 1e-8                                  # rule 3 tolerance of both runs
BLUE = np.array([66, 113, 184]) / 255       # nsBlue, nsGreen, nsRed of nullspace-styles.tex
GREEN = np.array([60, 160, 90]) / 255
RED = np.array([205, 50, 50]) / 255
ORANGE = np.array([230, 140, 40]) / 255
SECTOR_TINTS = [tuple(0.16 * c + 0.84) for c in (BLUE, GREEN)]   # alternating around the centre


# --------------------------------------------------------------------------
# Solves
# --------------------------------------------------------------------------

def the_case():
    return next(c for c in bc.make_cases() if c.name == "2026-aligned")


def solve_ours(case, max_iter=100):
    """benchBourne with solver H, stopped on rule 3 at TOL (the settings of run_mine)."""
    lines = ["box %r %r %r %r" % tuple(float(x) for x in case.box), "n %d" % case.n, "solver H",
             "alpha 1e300", "start 0.0", "volumetol %r" % TOL, "shifttol %r" % TOL, "maxiter %d" % max_iter,
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
    p = np.reshape(r["generators"], (-1, 2))
    V, c = np.asarray(r["volumes"]), np.reshape(r["centroids"], (-1, 2))
    vol, shift = bc.common_residuals(case, V, c, p, math.inf)
    return dict(method="our H", status=r["status"], converged=r["converged"], iterations=r["iterations"],
                attempts=r["attempts"], numBuilds=r["numBuilds"], seconds=r["seconds"],
                solverVolumeResidual=r["volumeResidual"], solverShiftResidual=r["shiftResidual"],
                volumeResidual=vol, shiftResidual=shift, generators=p.tolist(),
                weights=list(r["shiftedWeights"]), volumes=V.tolist(), centroids=c.tolist())


def solve_theirs(case, max_iter, max_seconds, every=250):
    """Their FW (t = 100) from the start of the case, with the gap rule disabled and the
    inner OT tolerance of run_theirs_common; stops at max_iter iterations, at max_seconds
    (evaluations not charged), or on rule 3 at TOL. Rule 3 is evaluated every `every`
    iterations and at the stop."""
    ref.COUNTERS.reset()
    ref.OT_REL_TOL["value"] = max(bc.OT_MIN_REL_TOL, min(1e-3, 0.1 * TOL))
    y = case.t / (case.t - 1) * case.b
    t0 = time.perf_counter()
    last, trace, stop = {}, [], ["iteration cap"]

    def elapsed():
        return time.perf_counter() - t0 - ref.COUNTERS.evalSeconds

    def evaluate(it, Z):
        with ref.Evaluation():
            w, c, V = ref.solve_OT_problem(Z, case.v, case.om, strict=False)
            vol, shift = bc.common_residuals(case, V, c, Z, math.inf)
        trace.append(dict(iteration=it, seconds=elapsed(), volumeResidual=vol, shiftResidual=shift))
        return vol, shift

    def cb(it, B, gap):
        last.update(iteration=it, B=B.copy(), gap=float(gap), seconds=elapsed())
        if it % every == 0:
            vol, shift = evaluate(it, y - B)
            print("   FW %5d  %.1f s  vol %.2e  shift %.2e  gap %.2e" % (it, elapsed(), vol, shift, gap), flush=True)
            if vol <= TOL and shift <= TOL:
                stop[0] = "rule 3"
                return True
        if elapsed() > max_seconds:
            stop[0] = "time cap"
            return True
        return False

    ref.convex_order_projection_frank_wolfe(case.om, y, case.start, R0=case.R0, masses=case.v,
                                            niter=max_iter, tol=-1.0, callback=cb)
    Z = y - last["B"]
    with ref.Evaluation():
        w, c, V = ref.solve_OT_problem(Z, case.v, case.om, strict=False)
    vol, shift = bc.common_residuals(case, V, c, Z, math.inf)
    if not trace or trace[-1]["iteration"] != last["iteration"]:
        trace.append(dict(iteration=last["iteration"], seconds=last["seconds"], volumeResidual=vol,
                          shiftResidual=shift))
    return dict(method="their FW (t = 100)", status=stop[0], iterations=last["iteration"],
                seconds=last["seconds"], otSolves=ref.COUNTERS.otSolves, gap=last["gap"],
                maxIter=max_iter, maxSeconds=max_seconds, innerOTRelTol=ref.OT_REL_TOL["value"],
                volumeResidual=vol, shiftResidual=shift, trace=trace, generators=Z.tolist(),
                weights=np.asarray(w).tolist(), volumes=np.asarray(V).tolist(), centroids=np.asarray(c).tolist())


# --------------------------------------------------------------------------
# Cells
# --------------------------------------------------------------------------

def clip(poly, a, b):
    """The part of the polygon with a . x <= b (Sutherland-Hodgman)."""
    out = []
    for k in range(len(poly)):
        P, Q = poly[k], poly[(k + 1) % len(poly)]
        fp, fq = a @ P - b, a @ Q - b
        if fp <= EPS:
            out.append(P)
        if (fp < -EPS and fq > EPS) or (fp > EPS and fq < -EPS):
            out.append(P + fp / (fp - fq) * (Q - P))
    return out


def power_cells(Y, w):
    """The cells {x : |x - Y_i|^2 - w_i <= |x - Y_j|^2 - w_j} clipped to the unit square, by
    half-plane clipping: first against the nearest generators, then against any generator
    whose half-plane still cuts the polygon, until none does."""
    sq = np.sum(Y * Y, axis=1) - w
    cells = []
    for i in range(len(Y)):
        poly = [np.array(p, dtype=float) for p in ((0, 0), (1, 0), (1, 1), (0, 1))]
        A = 2.0 * (Y - Y[i])
        B = sq - sq[i]
        order = list(np.argsort(np.linalg.norm(Y - Y[i], axis=1))[1:40])
        while order and poly:
            for j in order:
                poly = clip(poly, A[j], B[j])
                if not poly:
                    break
            if not poly:
                break
            P = np.array(poly)
            viol = (P @ A.T - B).max(axis=0)
            viol[i] = -1
            order = list(np.nonzero(viol > EPS)[0])
        cells.append(np.array(poly))
    return cells


def area(P):
    x, y = P.T
    return 0.5 * abs(x @ np.roll(y, -1) - y @ np.roll(x, -1))


def centroid(P):
    x, y = P.T
    c = x * np.roll(y, -1) - np.roll(x, -1) * y
    return np.array([((x + np.roll(x, -1)) * c).sum(), ((y + np.roll(y, -1)) * c).sum()]) / (3.0 * c.sum())


def sector_of(c):
    d = c - 0.5
    return int(math.floor(math.atan2(d[1], d[0]) / (math.pi / 4))) % 8


def centre_valence(cells, tol=1e-7):
    """The number of cells with a vertex within tol of the centre."""
    return sum(bool(len(P)) and np.min(np.linalg.norm(P - 0.5, axis=1)) < tol for P in cells)


# --------------------------------------------------------------------------
# Drawing
# --------------------------------------------------------------------------

def display_generators(sol, t):
    """Generators and weights to draw. Ours: the generators of H. FW: its iterate Z
    (y - B) is the generator set of the dilated problem, of size about 1/(t - 1); we draw
    p = (t - 1) Z, the generators of the equivalent H_alpha with alpha = t - 1 (as in
    run_theirs_common), with the weights that keep the same cells:
    |x - sZ|^2 - w' = s (|x - Z|^2 - w) + (1 - s)|x|^2 for w' = s^2 |Z|^2 - s |Z|^2 + s w."""
    Y, w = np.asarray(sol["generators"]), np.asarray(sol["weights"])
    if "FW" not in sol["method"]:
        return Y, w
    s = t - 1.0
    q = np.sum(Y * Y, axis=1)
    return s * Y, s * s * q - s * q + s * w


def draw_cells(ax, cells, sectors, Y, b, c, V, lw, ms, arrow_min):
    ax.add_collection(PolyCollection([P for P in cells if len(P)],
                                     facecolors=[SECTOR_TINTS[s % 2] for P, s in zip(cells, sectors) if len(P)],
                                     edgecolors="black", linewidths=lw, joinstyle="round", zorder=1))
    if arrow_min is not None:
        err = np.linalg.norm(c - b, axis=1) / np.sqrt(V)
        k = err > arrow_min
        if k.any():
            ax.add_collection(LineCollection(np.stack([c[k], b[k]], 1), colors=[RED], linewidths=3.0 * lw,
                                             capstyle="round", zorder=3))
    ax.plot(b[:, 0], b[:, 1], "o", markerfacecolor="white", markeredgecolor="black", markersize=1.7 * ms,
            markeredgewidth=0.25 * ms, zorder=4)
    ax.plot(Y[:, 0], Y[:, 1], "o", color="black", markersize=ms, markeredgewidth=0, zorder=5)


def frame(ax, color, lw):
    ax.set_xticks([])
    ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_edgecolor(color)
        sp.set_linewidth(lw)


def method_block(fig, x, y, h, gap, sol, cells, sectors, b, Vbar, label, zooms, args):
    """The full diagram (h x h, in figure fractions of the height) and, to its right, the
    zooms stacked in one column."""
    W, H = fig.get_figwidth(), fig.get_figheight()
    Y, wts = display_generators(sol, args.t)
    c = np.asarray(sol["centroids"])
    ax = fig.add_axes([x, y, h * H / W, h])
    draw_cells(ax, cells, sectors, Y, b, c, Vbar, 0.2, 0.9, args.arrow_min)
    ax.add_patch(Rectangle((0, 0), 1, 1, fill=False, edgecolor="black", linewidth=0.8, zorder=6))
    for (cx, cy), r, col in zooms:
        ax.add_patch(Rectangle((cx - r, cy - r), 2 * r, 2 * r, fill=False, edgecolor=col, linewidth=0.8, zorder=7))
    pad = 0.012
    ax.set_xlim(-pad, 1 + pad)
    ax.set_ylim(-pad, 1 + pad)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.text(-pad, 1 + 2 * pad, label, ha="left", va="bottom", fontsize=7)
    hz = (h - gap * (len(zooms) - 1)) / len(zooms)
    for k, ((cx, cy), r, col) in enumerate(zooms):
        za = fig.add_axes([x + (h + gap) * H / W, y + h - (k + 1) * hz - k * gap, hz * H / W, hz])
        draw_cells(za, cells, sectors, Y, b, c, Vbar, 0.45, 2.2, args.arrow_min)
        za.set_xlim(cx - r, cx + r)
        za.set_ylim(cy - r, cy + r)
        za.set_aspect("equal")
        frame(za, col, 1.0)
    return (h + gap + hz) * H / W


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--resolve", action="store_true", help="recompute the solutions (else read the cache)")
    ap.add_argument("--fw-iter", type=int, default=5000, help="iteration cap of FW (that of the overnight run)")
    ap.add_argument("--fw-seconds", type=float, default=900.0, help="time cap of FW (evaluations not charged)")
    ap.add_argument("--width", type=float, default=6.9, help="figure width (in)")
    ap.add_argument("--zoom", type=float, default=0.075, help="half-side of the zooms")
    ap.add_argument("--arrow-min", type=float, default=0.02,
                    help="draw centroid-target segments where |c - b| / sqrt(V) exceeds this")
    ap.add_argument("--out", default=OUT)
    args = ap.parse_args()

    case = the_case()
    args.t = case.t
    if args.resolve or not os.path.exists(CACHE):
        print("our H ...", flush=True)
        ours = solve_ours(case)
        print("   %s, %d it, %.3g s, residuals %.1e / %.1e" % (
            ours["status"], ours["iterations"], ours["seconds"], ours["volumeResidual"], ours["shiftResidual"]))
        print("their FW ...", flush=True)
        theirs = solve_theirs(case, args.fw_iter, args.fw_seconds)
        print("   %s, %d it, %.3g s, residuals %.1e / %.1e" % (
            theirs["status"], theirs["iterations"], theirs["seconds"], theirs["volumeResidual"],
            theirs["shiftResidual"]))
        os.makedirs(os.path.dirname(CACHE), exist_ok=True)
        with open(CACHE, "w") as f:
            json.dump(dict(case=case.name, n=case.n, tol=TOL, targets=case.b.tolist(), volumes=case.v.tolist(),
                           ours=ours, theirs=theirs), f)
    with open(CACHE) as f:
        data = json.load(f)
    b, Vbar = np.asarray(data["targets"]), np.asarray(data["volumes"])

    # The true diagram (the cross-shaped grid of the case with the weights of its volumes),
    # to tell whether a solution recovers its cells.
    wt, _, _ = ref.solve_OT_problem(case.truth, case.v, case.om, relTol=bc.OT_MIN_REL_TOL)
    true_cells = power_cells(case.truth, wt)
    sols = [data["ours"], data["theirs"]]
    allcells = []
    for sol in sols:
        Y, w = display_generators(sol, args.t)
        cells = power_cells(Y, w)
        A = np.array([area(P) if len(P) else 0.0 for P in cells])
        C = np.array([centroid(P) if len(P) else [np.nan, np.nan] for P in cells])
        dv = np.max(np.abs(A - np.asarray(sol["volumes"])) / Vbar)
        dc = np.max(np.linalg.norm(C - np.asarray(sol["centroids"]), axis=1) / np.sqrt(Vbar))
        assert dv < 1e-6 and dc < 1e-6, "%s: cells off the solver's (%g, %g)" % (sol["method"], dv, dc)
        e = np.linalg.norm(C - b, axis=1) / np.sqrt(Vbar)
        same = max(np.max(np.abs(np.sort(P, 0) - np.sort(T, 0))) if len(P) == len(T) else np.inf
                   for P, T in zip(cells, true_cells))
        print("%-20s %-14s %5d it  %7.3g s  vol %.1e  shift %.1e  (|c - b| / sqrt(V): max %.2e, median %.2e, "
              "%d cells > %g)  centre shared by %d cells  max vertex distance to the true cells %.1e" % (
                  sol["method"], sol["status"], sol["iterations"], sol["seconds"], sol["volumeResidual"],
                  sol["shiftResidual"], e.max(), np.median(e), (e > args.arrow_min).sum(), args.arrow_min,
                  centre_valence(cells), same))
        allcells.append(cells)
    sectors = [sector_of(centroid(P)) for P in allcells[0]]
    print("sector sizes (a):", np.bincount(sectors, minlength=8).tolist())

    # zooms: the centre, and the inner corner of the cross at the top left, where the
    # centroid errors of FW are largest
    zooms = [((0.5, 0.5), args.zoom, RED), ((0.31, 0.69), args.zoom, ORANGE)]
    gap, top = 0.03, 0.07           # fractions of the height
    h = 1 - top
    # width of one block = (h + gap + hz) H / W with hz = (h - gap) / 2; two blocks and a gap
    per = h + gap + (h - gap) / 2
    H = args.width / (2 * per + 2 * gap)
    fig = plt.figure(figsize=(args.width, H))
    x = 0.0
    for sol, cells, lab in zip(sols, allcells, ("(a) our H", "(b) Frank–Wolfe (Bourne et al. 2026), 5000 it.")):
        x += method_block(fig, x, 0.0, h, gap, sol, cells, sectors, b, Vbar, lab, zooms, args) + 2 * gap * H / args.width
    fig.savefig(args.out)
    plt.close(fig)
    print("->", os.path.relpath(args.out, ROOT), "(%.2f x %.2f in)" % (args.width, H))


if __name__ == "__main__":
    main()
