"""Shared parts of the benchmarks against Bourne, Pearce & Roper (2025) and
Bourne, Gallouet, Merigot & Natale (2026): the cases (generated exactly as in
their notebooks), their methods, ours (max H, max H_alpha through
benchBourne), the termination rules, and the driver loop. The drivers are
  run_realizable.py   realizable data (an exact power diagram exists): their
                      method and our H
  run_approximate.py  approximate data (fit only): their method and our
                      H_alpha

Each method stops on its own rule at 1e-8 (--tol), and the other rule is
measured at its solution:
  * rule 1, the SLSQP rule of Bourne 2025 ("theirs": relative change of
    their objective, ptol = 1e-8): stops their SLSQP; measured at our
    solution. The FW gap of Bourne 2026 is not used as a rule
    (latex/results.tex, "Termination"), so FW, which has no rule of its
    own, stops on rule 3;
  * a common per-cell L-infinity residual ("common"): of H on the
    realizable cases, for every method (their FW with dilation t solves
    H_alpha with alpha = t - 1, so it stalls at about 1/t there); of
    H_alpha on the approximate cases, with alpha = t - 1 for their FW and
    our H_alpha (the projection of the dilated targets):
        max_i |V_i - Vbar_i| / Vbar_i                 <= tol
        max_i |g_u,i| / Vbar_i^((d+1)/d)              <= tol
    with g_u = 2 V (cbar - c) - (2/alpha) Vbar (p - q) the shift gradient of
    the reference objective of the case (alpha = infinity: H).

Every method runs once per case (evaluations of the rules are not charged).
"""

import json
import math
import os
import subprocess
import sys
import tempfile
import time

import numpy as np
from scipy.spatial.distance import pdist, squareform

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, HERE)
import bourne_ref as ref  # noqa: E402

BENCH = os.path.join(ROOT, "build", "benchBourne")
DATA = {
    "2025": os.path.join(ROOT, "external", "bourne", "Laguerre-Polycrystalline-Microstructures", "Examples",
                         "Inverting_Laguerre_Tessellations_2024", "EBSD_data.txt"),
    "2026": os.path.join(ROOT, "external", "bourne", "Laguerre-fitting", "EBSD_data.txt"),
}


# --------------------------------------------------------------------------
# Cases (data generated exactly as in their notebooks)
# --------------------------------------------------------------------------

# Realizability of the data (v, b), independent of the solver:
#   REALIZABLE  b are the centroids of a power diagram with volumes v (the
#               data are generated from one), so the targets can be met
#               exactly.
#   APPROXIMATE no power diagram has volumes v and centroids b (perturbed or
#               measured data); only a fit is possible. For 2025-perturbed-*
#               this is stated by Bourne et al. (2025, Sec. 5.2.1: the
#               maximum of their H is positive); for the EBSD data the
#               centroids of real grains are not those of a power diagram.
REALIZABLE = "realizable"
# Relative objective change of the SLSQP rule of Bourne 2025 (their
# notebook: 1e-10).
SLSQP_PTOL = 1e-8
APPROXIMATE = "approximate"


class Case:
    """box, volumes v, targets b, start of their method, and the settings of
    their rule. group is REALIZABLE or APPROXIMATE. start is the start of
    their method; sites (= b) the sites of mine, whose initialization is
    our own (MakeJointStart). ref_alpha is the stiffness of the reference
    objective of the common rule for their method (inf: H, the volume and
    centroid match of the realizable cases; None: no common rule).
    variants: when the same data were used by both papers, one Case per
    paper with its own scaling and settings for their methods (named
    theirs-<paper>); our methods run once, on this Case."""

    def __init__(self, name, paper, box, v, b, **kw):
        self.name, self.paper, self.box, self.v, self.b = name, paper, box, v, b
        self.om = ref.make_box(*box)
        self.__dict__.update(kw)

    @property
    def n(self):
        return self.v.size

    @property
    def realizable(self):
        return self.group == REALIZABLE


def _voronoi(amp, p_seed):
    v, b, om, X0 = ref.generate_perturbed_voronoi_data(1, 1, 20, x_seed=999, p_seed=p_seed, amp=amp)
    return v, b, X0


def make_cases():
    cases = []

    # SLSQP ftol of their 2025 rule: ptol |H(X_init)|. Their notebook uses
    # ptol = 1e-10; we use 1e-8 (SLSQP_PTOL).
    # 2025 Sec. 5.1: recover a Voronoi diagram from a random start.
    v, b, X0 = _voronoi(0, None)
    xinit = np.random.rand(20, 2)  # continues the RNG stream, as in the notebook
    cases.append(Case("2025-invert", "2025", (0, 0, 1, 1), v, b, truth=X0, start=xinit, sites=b.copy(),
                      objective="H", ptol=SLSQP_PTOL, maxiter=1000, R0=math.sqrt(20), d=1e-3,
                      group=REALIZABLE, ref_alpha=math.inf))

    # 2025 Sec. 5.2.1 and 5.2.2: perturbed Voronoi centroids, start at b.
    for tag, amp in (("small", 1e-3), ("large", 5e-2)):
        v, b, X0 = _voronoi(amp, 1002)
        cases.append(Case("2025-perturbed-" + tag, "2025", (0, 0, 1, 1), v, b, truth=X0, start=b.copy(),
                          sites=b.copy(), objective="H", ptol=SLSQP_PTOL, maxiter=1000, R0=math.sqrt(20), d=1e-3,
                          group=APPROXIMATE, ref_alpha=None))


    # 2026 Sec. 5.2: Voronoi n = 20, t = 100, FW from the true seeds. On
    # realizable data every method is measured with the residual of H; FW
    # with dilation t solves H_alpha with alpha = t - 1 instead, so its
    # residual of H stalls at about 1/t.
    v, b, X0 = _voronoi(0, None)
    cases.append(Case("2026-voronoi", "2026", (0, 0, 1, 1), v, b, truth=X0, start=X0.copy(), sites=b.copy(),
                      t=100.0, gap=1e-14, niter=4000, R0=1.0, group=REALIZABLE, ref_alpha=math.inf))

    # 2026 Sec. 5.3: aligned seeds (cross), n = 728, t = 100.
    m = 40
    c0 = np.linspace(0, 1, m)
    gx, gy = np.meshgrid(c0, c0)
    P = np.stack([gx.ravel(), gy.ravel()], 1)
    cross = np.logical_or(np.logical_and(P[:, 0] > .33, P[:, 0] < .67), np.logical_and(P[:, 1] < .67, P[:, 1] > .33))
    centre = np.logical_or(np.logical_or(P[:, 0] < .33, P[:, 0] > .67), np.logical_or(P[:, 1] > .67, P[:, 1] < .33))
    Y = P[np.logical_and(cross, centre)].copy()
    vv = np.ones(len(Y)) / len(Y)
    om = ref.make_box(0, 0, 1, 1)
    # Their notebook solves this OT to 1e-4 min(v); solved to 1e-7 min(v)
    # here (about 1e-10 absolute, what pysdot reaches on this lattice). The
    # targets are the volumes and centroids of that diagram, so that they are
    # realizable up to round-off (as in makeSynthetic), not 1 / n, which the
    # diagram only meets to the OT tolerance.
    _, cc, vv = ref.solve_OT_problem(Y, vv, om, relTol=1e-7)
    cases.append(Case("2026-aligned", "2026", (0, 0, 1, 1), vv, cc, truth=Y, start=cc.copy(), sites=cc.copy(),
                      t=100.0, gap=1e-16, niter=1000, R0=1.0, group=REALIZABLE, ref_alpha=math.inf))

    # The EBSD map of both papers (the same file in both repositories), one
    # case: our methods run once, in the unit square; their methods run with
    # each paper's scaling and settings (variants).
    # 2025 Sec. 5.3: in a box of side sqrt(sum of areas) (microns), least
    # squares from b.
    E = np.loadtxt(DATA["2025"])
    areas, cents = E[:, 3], E[:, 1:3]
    box = math.sqrt(areas.sum())
    om = ref.make_box(0, 0, box, box)
    v, b = ref.make_data_compatible(areas * box ** 2 / areas.sum(), cents, om)
    ebsd2025 = Case("ebsd", "2025", (0, 0, box, box), v, b, truth=None, start=b.copy(), sites=b.copy(),
                    objective="F", ptol=1e-8, maxiter=10000, d=1e-3, group=APPROXIMATE, ref_alpha=None)
    # 2026 Sec. 5.4: in the unit box, t = 10, FW from 0.5.
    E = np.loadtxt(DATA["2026"])
    areas, cents = E[:, 3], E[:, 1:3]
    cents = cents / math.sqrt(areas.sum())
    areas = areas / areas.sum()
    om = ref.make_box(0, 0, 1, 1)
    v, b = ref.make_data_compatible(areas, cents, om)
    ebsd2026 = Case("ebsd", "2026", (0, 0, 1, 1), v, b, truth=None, start=.5 * np.ones(b.shape),
                    sites=b.copy(), t=10.0, gap=1e-9, niter=500, R0=1.0, group=APPROXIMATE, ref_alpha=9.0)
    ebsd2026.variants = [ebsd2025, ebsd2026]
    cases.append(ebsd2026)
    return cases


# --------------------------------------------------------------------------
# Rules
# --------------------------------------------------------------------------

def common_residuals(case, V, c, p, alpha):
    """(volume, shift) residuals of the reference objective at generators p
    (alpha = inf: H; otherwise H_alpha with sites q = b)."""
    vb = case.v
    vol = np.max(np.abs(V - vb) / vb)
    g = 2 * V[:, None] * (case.b - c)
    if math.isfinite(alpha):
        g -= 2 / alpha * vb[:, None] * (p - case.b)
    shift = np.max(np.linalg.norm(g, axis=1) / (2 * vb ** 1.5))
    return float(vol), float(shift)


def their_objective(case, V, c, p):
    """H of the 2025 paper (sum Vbar (b - c) . X) or their least-squares f."""
    if case.objective == "H":
        return float(np.sum(case.v * np.sum((case.b - c) * p, axis=1)))
    return float(np.linalg.norm(case.v[:, None] * (case.b - c)) ** 2 / ref.normalisation(case.om.measure(), case.n))


class Stops:
    """First iterate passing each rule."""

    def __init__(self, tols):
        self.tols = sorted(tols, reverse=True)
        self.common = {t: None for t in self.tols}
        self.theirs = None

    def record_common(self, k, seconds, solves, vol, shift):
        for t in self.tols:
            if self.common[t] is None and vol <= t and shift <= t:
                self.common[t] = dict(iteration=k, seconds=seconds, solves=solves, volumeResidual=vol,
                                      shiftResidual=shift)

    def all_common(self):
        return all(v is not None for v in self.common.values())


# --------------------------------------------------------------------------
# Their methods
# --------------------------------------------------------------------------

def run_theirs_default(case, args):
    """Their method with their exact settings and stopping rule (2026: the
    FW gap at their notebook tolerance, kept as a reference run only, not
    as a rule of the comparison)."""
    ref.COUNTERS.reset()
    ref.OT_REL_TOL["value"] = 1e-3 if case.paper == "2025" else 1e-4
    t0 = time.perf_counter()
    out = dict(case=case.name, method="theirs", termination="theirs")
    try:
        if case.paper == "2025":
            # Time cap, checked between iterations as in run_theirs_common.
            def capped(X):
                return time.perf_counter() - t0 > args.max_seconds
            if case.objective == "H":
                res, X = ref.max_H_solver(case.start, case.b, case.v, case.om, R0=case.R0, d=case.d,
                                          max_iter=case.maxiter, ptol=case.ptol, callback=capped)
            else:
                res, X = ref.min_F_solver(case.start, case.b, case.v, case.om, d=case.d, max_iter=case.maxiter,
                                          ptol=case.ptol, callback=capped)
            timedOut = time.perf_counter() - t0 > args.max_seconds
            out.update(converged=bool(res.success) and not timedOut,
                       status="time cap" if timedOut else str(res.message), iterations=int(res.nit))
        else:
            y = case.t / (case.t - 1) * case.b
            B, k, status = ref.convex_order_projection_frank_wolfe(case.om, y, case.start, R0=case.R0, masses=case.v,
                                                                   niter=case.niter, tol=case.gap)
            X = y - B
            out.update(converged=status == "converged", status=status, iterations=int(k))
    except Exception as e:  # noqa: BLE001
        out.update(converged=False, status="error: %s" % e, iterations=None)
        X = None
    out.update(seconds=time.perf_counter() - t0, otSolves=ref.COUNTERS.otSolves, otNewton=ref.COUNTERS.otNewton)
    if X is not None:
        out.update(_final_errors(case, X=X))
        if case.ref_alpha is not None:
            # Rule 3 at their solution, with the residual of the case.
            try:
                with ref.Evaluation():
                    _, c, V = ref.solve_OT_problem(np.asarray(X), case.v, case.om, relTol=OT_MIN_REL_TOL, strict=False)
                vol, shift = common_residuals(case, V, c, np.asarray(X), case.ref_alpha)
                tol = min(args.tol)
                out["rule3AtSolution"] = dict(volumeResidual=vol, shiftResidual=shift, tol=tol,
                                              met=vol <= tol and shift <= tol)
            except Exception as e:  # noqa: BLE001
                out["rule3AtSolution"] = dict(error="%s" % e)
    return out


def _final_errors(case, X=None, V=None, c=None):
    """Volume and centroid errors, and their least-squares error f, of the
    cells (V, c), or of the power diagram with generators X and the
    prescribed volumes."""
    if V is None:
        try:
            with ref.Evaluation():
                _, c, V = ref.solve_OT_problem(np.asarray(X), case.v, case.om, relTol=1e-6)
        except Exception as e:  # noqa: BLE001
            return dict(finalError="%s" % e)
    f = np.linalg.norm(case.v[:, None] * (case.b - c)) ** 2 / ref.normalisation(case.om.measure(), case.n)
    return dict(finalVolumeError=float(np.max(np.abs(V - case.v) / case.v)),
                finalCentroidError=float(np.max(np.linalg.norm(c - case.b, axis=1)) / math.sqrt(case.om.measure())),
                finalF=float(f))


OT_MIN_REL_TOL = 1e-9


def run_theirs_common(case, args):
    """Their method with the common rule: stop disabled, the residual of the
    reference objective evaluated at every iterate (not charged)."""
    if case.ref_alpha is None:
        return None
    stops = Stops(args.tol)
    ref.COUNTERS.reset()
    # Inner OT tight enough for the volume part of the common rule, within
    # what pysdot reaches (about 1e-10 min(v)); smaller tolerances are then
    # reported as not reached. The evaluations accept the diagram pysdot
    # reaches when it misses this tolerance (strict=False), so that a missed
    # evaluation shows up in the volume residual instead of ending the run.
    ref.OT_REL_TOL["value"] = max(OT_MIN_REL_TOL, min(1e-3, 0.1 * min(args.tol)))
    alpha = case.ref_alpha
    t0 = time.perf_counter()
    k = [0]

    def elapsed():
        return time.perf_counter() - t0 - ref.COUNTERS.evalSeconds

    def done():
        return stops.all_common() or elapsed() > args.max_seconds

    status = "max iterations"
    # Rule 3 residuals of the last iterate: how far the run is from the
    # solution when it stops at the cap.
    last = {}
    try:
        if case.paper == "2025":
            def cb(X):
                k[0] += 1
                with ref.Evaluation():
                    _, c, V = ref.solve_OT_problem(X, case.v, case.om, strict=False)
                res3 = common_residuals(case, V, c, X, alpha)
                last.update(iteration=k[0], volumeResidual=res3[0], shiftResidual=res3[1])
                stops.record_common(k[0], elapsed(), ref.COUNTERS.otSolves, *res3)
                return done()
            res, _ = ref.max_H_solver(case.start, case.b, case.v, case.om, R0=case.R0, d=case.d,
                                      max_iter=args.max_iter_theirs, ptol=0.0, callback=cb)
            status = str(res.message)
        else:
            s = case.t / (case.t - 1)
            y = s * case.b

            def cb(it, B, gap):
                # Generators of the current FW point and their diagram (the
                # next LMO, recomputed here as an evaluation).
                Z = y - B
                with ref.Evaluation():
                    _, c, V = ref.solve_OT_problem(Z, case.v, case.om, strict=False)
                p = alpha * Z if math.isfinite(alpha) else Z
                res3 = common_residuals(case, V, c, p, alpha)
                last.update(iteration=it, volumeResidual=res3[0], shiftResidual=res3[1])
                stops.record_common(it, elapsed(), ref.COUNTERS.otSolves, *res3)
                return done()
            _, _, status = ref.convex_order_projection_frank_wolfe(case.om, y, case.start, R0=case.R0, masses=case.v,
                                                                   niter=args.max_iter_theirs, tol=-1.0, callback=cb)
    except Exception as e:  # noqa: BLE001
        status = "error: %s" % e
    return dict(case=case.name, method="theirs", termination="common", referenceAlpha=alpha, status=status,
                stops={"%g" % t: v for t, v in stops.common.items()}, seconds=elapsed(),
                otSolves=ref.COUNTERS.otSolves, final=last or None)


# --------------------------------------------------------------------------
# My methods
# --------------------------------------------------------------------------


def trim_to_solution(traj, generators):
    """The trajectory up to the iterate returned by the solver. With restarts,
    the trajectory chains all the attempts, each from the same start, and the
    returned solution (the converged attempt, or else the one with the
    smallest residual) need not be the last iterate."""
    if not traj:
        return traj
    g = np.asarray(generators, dtype=float)
    for k in range(len(traj) - 1, -1, -1):
        p = np.asarray(traj[k]["generators"], dtype=float)
        if p.shape == g.shape and np.array_equal(p, g):
            return traj[:k + 1]
    return traj

def run_mine(case, solver, alpha, start, args):
    """One benchBourne solve, stopped on rule 3 at the smallest --tol (the
    larger ones are read off the trajectory); rule 1 is measured at its
    solution."""
    tight = min(args.tol)
    lines = ["box %r %r %r %r" % tuple(float(x) for x in case.box), "n %d" % case.n, "solver %s" % solver,
             "alpha %r" % (float(alpha) if math.isfinite(alpha) else 1e300), "start %r" % float(start),
             "volumetol %r" % tight, "shifttol %r" % tight, "maxiter %d" % args.max_iter_mine,
             "volumes " + " ".join(repr(float(x)) for x in case.v),
             "targets " + " ".join(repr(float(x)) for x in case.b.ravel()),
             "sites " + " ".join(repr(float(x)) for x in np.asarray(case.sites).ravel())]
    with tempfile.TemporaryDirectory() as d:
        cpath, opath = os.path.join(d, "case.txt"), os.path.join(d, "out.json")
        with open(cpath, "w") as f:
            f.write("\n".join(lines) + "\n")
        subprocess.run([BENCH, cpath, opath], check=True)
        with open(opath) as f:
            r = json.load(f)

    traj = trim_to_solution(r.pop("trajectory"), r["generators"])
    method = solver if solver == "H" else "Halpha(%.3g)" % alpha
    if len(args.start) > 1:
        method += " s=%g" % start
    out = dict(case=case.name, method=method, solverStatus=r["status"], solverIterations=r["iterations"],
               solverSeconds=r["seconds"], numBuilds=r["numBuilds"], startSeconds=r["startSeconds"],
               # The residuals of the solver's termination test at the returned solution.
               volumeResidual=r["volumeResidual"], shiftResidual=r["shiftResidual"],
               phaseSeconds=dict(build=r["buildSeconds"], hessian=r["hessianSeconds"], factor=r["factorSeconds"]))
    if r["volumes"]:
        out.update(_final_errors(case, V=np.asarray(r["volumes"]), c=np.reshape(r["centroids"], (-1, 2))))

    # Common rule, on the objective of the run (H or H_alpha): each method
    # is tested for stationarity of the problem it solves.
    ref_alpha = alpha
    stops = Stops(args.tol)
    for k, it in enumerate(traj):
        V, c = np.asarray(it["volumes"]), np.reshape(it["centroids"], (-1, 2))
        p = np.reshape(it["generators"], (-1, 2))
        stops.record_common(k, it["seconds"], it["numBuilds"], *common_residuals(case, V, c, p, ref_alpha))
    out["common"] = dict(referenceAlpha=ref_alpha, stops={"%g" % t: v for t, v in stops.common.items()})

    # Their rule 1, measured at our solution.
    out["theirs"] = _rule1_at_solution(case, traj)
    return out


def _rule1_at_solution(case, traj):
    """Rule 1 (SLSQP ftol of the 2025 paper) at the final iterate: the last
    change of their objective relative to its value at the start, and the
    0.1 % volume condition. None for the 2026 cases (no rule of theirs)."""
    if case.paper != "2025" or len(traj) < 2:
        return None
    obj = [their_objective(case, np.asarray(it["volumes"]), np.reshape(it["centroids"], (-1, 2)),
                           np.reshape(it["generators"], (-1, 2))) for it in (traj[0], traj[-2], traj[-1])]
    ratio = abs(obj[2] - obj[1]) / abs(obj[0])
    volOk = bool(np.max(np.abs(np.asarray(traj[-1]["volumes"]) - case.v)) < 1e-3 * np.min(case.v))
    return dict(rule="SLSQP ftol (ptol %g)" % case.ptol, ratio=float(ratio), met=bool(volOk and ratio <= case.ptol),
                iteration=len(traj) - 1, seconds=traj[-1]["seconds"], solves=traj[-1]["numBuilds"])


# --------------------------------------------------------------------------
# Driver
# --------------------------------------------------------------------------

def _start_residuals(case):
    """|c0_i - b_i| per cell, with c0 the centroids of the power diagram of
    the sites q = b with the prescribed volumes (max F at q): the tightest OT
    solve that pysdot reaches, looser if needed."""
    for tol in (OT_MIN_REL_TOL, 1e-8, 1e-7, 1e-6):
        try:
            with ref.Evaluation():
                _, c0, _ = ref.solve_OT_problem(np.asarray(case.sites), case.v, case.om, relTol=tol)
            return np.linalg.norm(c0 - case.b, axis=1)
        except RuntimeError:
            continue
    raise RuntimeError("the volume solve at the sites cannot be solved to 1e-6 min(v)")


def _site_distances(case):
    D = squareform(pdist(np.asarray(case.sites)))
    np.fill_diagonal(D, np.inf)
    return D


def distinct_alpha(case):
    """Bound eq:regularized-distinct at the sites q = b: below it, the
    maximizer of H_alpha is a power diagram (distinct generators),
        min_{i != j} |q_i - q_j| / sqrt(resid (1/Vbar_i + 1/Vbar_j)),
    resid = sum_i Vbar_i |c0_i - b_i|^2."""
    e, V = _start_residuals(case), case.v
    resid = float(np.sum(V * e ** 2))
    return float(np.min(_site_distances(case) / np.sqrt(resid * (1.0 / V[:, None] + 1.0 / V[None, :]))))


def practical_alpha(case):
    """Per-cell estimate eq:regularized-practical of the same bound, assuming
    |u_i| <= alpha |c0_i - b_i|:
        min_{i != j} |q_i - q_j| / (|c0_i - b_i| + |c0_j - b_j|).
    Not a guarantee: it held, with a margin of at least 1.5, on the
    synthetic and Bourne cases (2D and 3D), where the bound is up to 60 times
    smaller."""
    e = _start_residuals(case)
    return float(np.min(_site_distances(case) / (e[:, None] + e[None, :])))


def realizable_runs(case, methods, args):
    """(solver, alpha, start) of our runs on a realizable case: H only
    (H_alpha runs on the approximate cases)."""
    return [("H", math.inf, s) for s in args.start] if "H" in methods else []


def approximate_runs(case, methods, args):
    """(solver, alpha, start) of our runs on an approximate case: H_alpha over
    the alpha sweep (with alpha = t - 1 too for the 2026 cases, the same
    problem as their dilation t/(t-1))."""
    if "Halpha" not in methods:
        return []
    alphas = (list(args.alpha) + [r * distinct_alpha(case) for r in args.alpha_bound]
              + [r * practical_alpha(case) for r in args.alpha_practical])
    if case.paper == "2026":
        alphas.append(case.t - 1)
    return [("Halpha", a, s) for a in sorted(set(alphas)) for s in args.start]


def _fmt(s):
    return "-" if s is None else "%d it, %.3g s, %d solves" % (s["iteration"], s["seconds"], s["solves"])


def summarize(results, tols):
    lines = []
    for r in results:
        if r is None:
            continue
        head = "%-22s %-14s" % (r["case"], r["method"])
        if r["method"].startswith("theirs") and r["termination"] == "theirs":
            lines.append("%s rule 1   %s, %s it, %.3g s, %d OT solves" % (
                head, "met" if r["converged"] else r["status"], r["iterations"], r["seconds"], r["otSolves"]))
            r3 = r.get("rule3AtSolution")
            if r3 and "met" in r3:
                lines.append("%s rule 3 at solution: %.1e / %.1e (%s)" % (
                    head, r3["volumeResidual"], r3["shiftResidual"], "met" if r3["met"] else "not met"))
        elif r["method"].startswith("theirs"):
            for t in sorted(tols, reverse=True):
                lines.append("%s tol %-5g %s" % (head, t, _fmt(r["stops"]["%g" % t])))
            f = r.get("final")
            if f and not r["stops"]["%g" % min(tols)]:
                lines.append("%s stopped at iteration %d (%s): rule 3 residuals %.1e / %.1e" % (
                    head, f["iteration"], r["status"][:30], f["volumeResidual"], f["shiftResidual"]))
        else:
            # theirs (rule 1 at our solution) is absent on the 2026 cases.
            th = r.get("theirs")
            if th is not None:
                lines.append("%s rule 1 at solution: |dobj|/|obj0| %.1e (%s)" % (
                    head, th["ratio"], "met" if th["met"] else "not met"))
            if "common" in r:
                for t in sorted(tols, reverse=True):
                    lines.append("%s tol %-5g %s" % (head, t, _fmt(r["common"]["stops"]["%g" % t])))
    return "\n".join(lines)


def add_common_args(ap, group, methods, out):
    """Options shared by the drivers; the cases are those of group."""
    names = [c.name for c in make_cases() if c.group == group]
    ap.add_argument("--case", nargs="+", choices=names, default=names)
    ap.add_argument("--method", nargs="+", choices=methods, default=list(methods))
    ap.add_argument("--tol", nargs="+", type=float, default=[1e-8],
                    help="rule 3 tolerance (the smallest stops the runs; larger ones are read off the trajectory)")
    ap.add_argument("--start", nargs="+", type=float, default=[0.0],
                    help="s of MakeJointStart (irrelevant when the sites are the targets, the default)")
    ap.add_argument("--max-iter-mine", type=int, default=100)
    ap.add_argument("--max-iter-theirs", type=int, default=5000, help="iteration cap of the common-rule runs")
    ap.add_argument("--max-seconds", type=float, default=1800.0, help="time cap of the common-rule runs")
    ap.add_argument("--out", default=os.path.join(HERE, "results", out))
    ap.add_argument("--list", action="store_true", help="list the cases and exit")


def run(args, group, my_runs):
    """Runs their method and ours (my_runs) on the cases of args.case, writes
    one JSON per run and summary.txt to args.out, and prints the summary."""
    cases = [c for c in make_cases() if c.group == group and c.name in args.case]
    if args.list:
        for c in cases:
            print("%-22s n=%-4d paper %s" % (c.name, c.n, c.paper))
        return
    os.makedirs(args.out, exist_ok=True)
    results = []
    for case in cases:
        print("== %s (n = %d, %s)" % (case.name, case.n, case.group), flush=True)
        if "theirs" in args.method:
            variants = getattr(case, "variants", [case])
            for variant in variants:
                # Each method stops on its own rule: SLSQP (2025) on rule 1,
                # rule 3 measured at its solution; FW (2026), which has no
                # rule of its own (the gap is not a rule), on rule 3.
                if variant.paper == "2025":
                    runs = [run_theirs_default(variant, args)]
                else:
                    runs = [run_theirs_common(variant, args)]
                for r in runs:
                    if r is not None and len(variants) > 1:
                        r["method"] = "theirs-%s" % variant.paper
                    results.append(r)
        for solver, alpha, start in my_runs(case, args.method, args):
            # Ours always stops on rule 3 (common); rule 1 is measured at
            # its solution (None on the 2026 cases).
            try:
                r = run_mine(case, solver, alpha, start, args)
            except subprocess.CalledProcessError as e:
                # A crash of benchBourne on one run: record it and go on.
                print("   %s %s alpha %g start %g: benchBourne failed (exit %d)" % (
                    case.name, solver, alpha, start, e.returncode), flush=True)
                with open(os.path.join(args.out, "%s__%s__error.json" % (case.name, solver)), "w") as f:
                    json.dump(dict(case=case.name, method=solver, alpha=alpha, start=start,
                                   error="exit %d" % e.returncode), f, indent=1)
                continue
            if r["theirs"] is None:
                r.pop("theirs")
            results.append(r)
        for r in results:
            if r is not None and r["case"] == case.name:
                tag = "%s__%s__%s" % (r["case"], r["method"], r.get("termination", "all"))
                with open(os.path.join(args.out, tag + ".json"), "w") as f:
                    json.dump(r, f, indent=1)
        print(summarize([r for r in results if r is not None and r["case"] == case.name], args.tol), flush=True)
    text = summarize(results, args.tol)
    with open(os.path.join(args.out, "summary.txt"), "w") as f:
        f.write(text + "\n")
