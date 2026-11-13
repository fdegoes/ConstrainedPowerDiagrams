"""Reference methods of Bourne et al., ported from their code with
instrumentation (OT-solve counters, per-iterate callbacks, timing).

Sources (cloned in external/bourne):
  Laguerre-Polycrystalline-Microstructures/Examples/Inverting_Laguerre_Tessellations_2024/
      Inverting_Laguerre_Tessellations.ipynb  (BSD-3, D. Bourne)       -> 2025 paper
  Laguerre-fitting/ConvexOrderProjectionUtils.py, semidiscrete_convex_order.ipynb
      (MIT, A. Natale)                                                    -> 2026 paper

The numerical steps are unchanged; only the pandas tables and printing of the
notebooks are removed, and the iterates are reported to a callback.
"""

import time

import numpy as np
from pysdot import OptimalTransport, PowerDiagram
from pysdot.domain_types import ConvexPolyhedraAssembly
from scipy.optimize import minimize
from scipy.spatial.distance import cdist


class Counters:
    """OT solves and their damped Newton iterations, outside evaluation."""

    def __init__(self):
        self.reset()

    def reset(self):
        self.otSolves = 0
        self.otNewton = 0
        self.evalSeconds = 0.0
        self._evaluating = 0

    def add(self, ot):
        if self._evaluating == 0:
            self.otSolves += 1
            self.otNewton += max(len(ot.delta_m), len(ot.delta_w))


COUNTERS = Counters()


class Evaluation:
    """Context in which OT solves and time are not charged to the method."""

    def __enter__(self):
        COUNTERS._evaluating += 1
        self._t = time.perf_counter()

    def __exit__(self, *exc):
        COUNTERS._evaluating -= 1
        if COUNTERS._evaluating == 0:
            COUNTERS.evalSeconds += time.perf_counter() - self._t


def make_box(x0, y0, x1, y1):
    om = ConvexPolyhedraAssembly()
    om.add_box([x0, y0], [x1, y1])
    return om


# --------------------------------------------------------------------------
# 2025: Inverting Laguerre tessellations
# --------------------------------------------------------------------------

def init_w(x, om):
    """Meyron's rescaling initialization (their init_w)."""
    com = 0.5 * (om.max_position() + om.min_position())
    Lx = om.max_position()[0] - om.min_position()[0]
    Ly = om.max_position()[1] - om.min_position()[1]
    alpha = 0.49
    l = alpha * np.min([Lx, Ly]) / np.max(np.abs(x - x[0, :]))
    t = com - l * x[0, :]
    return -2. * np.sum(t * x, axis=1) - (l - 1.) * np.sum(x * x, axis=1)


# Relative volume tolerance of the 2025 OT solve: max|V - Vbar| < OT_REL_TOL min(Vbar).
OT_REL_TOL = {"value": 1e-3}


def solve_OT_problem(x, v, om, relTol=None, strict=True):
    """Their solve_OT_problem: weights and centroids with prescribed areas.
    strict=False returns the diagram pysdot reached when it misses the
    tolerance, instead of raising (for evaluations, whose residuals then
    include the volume error it left)."""
    err_tol = (OT_REL_TOL["value"] if relTol is None else relTol) * np.min(v)
    w0 = init_w(x, om)
    ot = OptimalTransport(positions=x, masses=v, weights=w0, domain=om,
                          obj_max_dm=err_tol, verbosity=0, linear_solver='Scipy')
    ot.adjust_weights()
    COUNTERS.add(ot)
    c = ot.get_centroids()
    av = ot.pd.integrals()
    w = ot.get_weights()
    dv = np.max(np.abs(v - av))
    if strict and dv > err_tol:
        raise RuntimeError("Volume tolerance not met (max error %g)" % dv)
    return w, c, av


def compute_H_dH(x, b, v, om):
    w, c, _ = solve_OT_problem(x, v, om)
    H = np.sum(v * np.sum((b - c) * x, axis=1))
    dH = v.reshape(-1, 1) * (b - c)
    return H, dH


def normalisation(vol, n):
    return (vol ** 3) / (n ** 2)


def compute_F(x, b, v, om):
    w, c, _ = solve_OT_problem(x, v, om)
    return (np.linalg.norm(v.reshape(-1, 1) * (b - c)) ** 2) / normalisation(om.measure(), v.size)


def pairwise(x, delta):
    tmp = cdist(x, x)
    indices = np.tril_indices_from(tmp, k=-1)
    pd = np.asarray(tmp[indices])
    return 0.5 * (pd ** 2 - delta ** 2)


def Dpairwise(x):
    tmp = cdist(x, x)
    n_seeds = np.shape(x)[0]
    indices = np.tril_indices_from(tmp, k=-1)
    n_constraints = indices[0].size
    J = np.zeros((n_constraints, n_seeds, 2))
    i0, i1 = indices
    J[np.arange(n_constraints), i0, :] = x[i0, :] - x[i1, :]
    J[np.arange(n_constraints), i1, :] = x[i1, :] - x[i0, :]
    return J


def constraint_ball(x, x0, R):
    return 0.5 * (R ** 2 - np.linalg.norm(x - x0) ** 2)


def Dconstraint_ball(x, x0):
    return x0 - x


def _slsqp(objective, x0, cons, ftol, maxiter, jac, callback):
    def cb(xk):
        if callback is not None and callback(xk):
            raise StopIteration
    return minimize(objective, x0, method='SLSQP', constraints=cons, callback=cb, jac=jac,
                    options={'ftol': ftol, 'disp': False, 'maxiter': maxiter})


def max_H_solver(xinit, b, v, om, R0=1, d=1e-3, max_iter=100, ptol=1e-6, callback=None):
    """Their max_H_solver (SLSQP on H with pairwise and ball constraints).
    callback(X) is called at every iterate; returning True stops."""
    n = v.size
    x0 = 0.5 * (om.max_position() + om.min_position())
    objective_fun = lambda x: -compute_H_dH(np.reshape(x, (n, 2)), b, v, om)[0]
    Dobjective_fun = lambda x: -np.reshape(compute_H_dH(np.reshape(x, (n, 2)), b, v, om)[1], -1)
    pair_con = {'type': 'ineq',
                'fun': lambda x: pairwise(np.reshape(x, (-1, 2)), d),
                'jac': lambda x: np.reshape(Dpairwise(np.reshape(x, (-1, 2))), (n * (n - 1) // 2, 2 * n))}
    ball_con = {'type': 'ineq',
                'fun': lambda x: constraint_ball(np.reshape(x, (n, 2)), x0, R0),
                'jac': lambda x: np.reshape(Dconstraint_ball(np.reshape(x, (n, 2)), x0), (1, -1))}
    Hinit = -objective_fun(xinit.flatten())
    ftol = np.abs(Hinit) * ptol
    cbw = None if callback is None else (lambda xk: callback(np.reshape(xk, (n, 2))))
    res = _slsqp(objective_fun, xinit.flatten(), [pair_con, ball_con], ftol, max_iter, Dobjective_fun, cbw)
    return res, np.reshape(res.x, (-1, 2))


def min_F_solver(xinit, b, v, om, d=1e-3, max_iter=100, ptol=1e-2, callback=None):
    """Their min_F_solver (SLSQP on the least-squares error, finite-difference
    gradient, pairwise constraints)."""
    n = v.size
    objective_fun = lambda x: compute_F(np.reshape(x, (n, 2)), b, v, om)
    pair_con = {'type': 'ineq',
                'fun': lambda x: pairwise(np.reshape(x, (-1, 2)), d),
                'jac': lambda x: np.reshape(Dpairwise(np.reshape(x, (-1, 2))), (n * (n - 1) // 2, 2 * n))}
    Finit = objective_fun(xinit.flatten())
    ftol = Finit * ptol
    cbw = None if callback is None else (lambda xk: callback(np.reshape(xk, (n, 2))))
    res = _slsqp(objective_fun, xinit.flatten(), pair_con, ftol, max_iter, None, cbw)
    return res, np.reshape(res.x, (-1, 2))


def make_data_compatible(vi, bi, om):
    total_vol = om.measure()
    vr = vi * total_vol / vi.sum()
    com = 0.5 * (om.max_position() + om.min_position())
    br = bi + com - np.sum(vr.reshape(-1, 1) * bi, axis=0) / total_vol
    return vr, br


def generate_perturbed_voronoi_data(Lx, Ly, n, x_seed=None, p_seed=None, amp=0):
    om = make_box(0, 0, Lx, Ly)
    Lmat = np.diag([Lx, Ly])
    if x_seed is not None:
        np.random.seed(x_seed)
    X0 = np.random.rand(n, 2) @ Lmat
    pd = PowerDiagram(positions=X0, weights=np.zeros(n), domain=om)
    b = pd.centroids()
    v = pd.integrals()
    if p_seed is not None:
        np.random.seed(p_seed)
    angles = 2 * np.pi * np.random.rand(n, 1)
    unit_vec = np.hstack([np.cos(angles), np.sin(angles)])
    amps = amp * np.random.rand(n, 1)
    br = b + unit_vec * amps
    va, ba = make_data_compatible(v, br, om)
    return va, ba, om, X0


# --------------------------------------------------------------------------
# 2026: Semi-discrete convex order (Frank-Wolfe)
# --------------------------------------------------------------------------

def get_centroids(domain, X, masses, rescale=True, R0=1.):
    """Their get_centroids: the linear minimization oracle (centroids of the
    power diagram of generators X with the prescribed masses)."""
    if rescale:
        X = X - np.reshape(np.mean(X, axis=0), (1, 2))
        X = X / np.max(np.linalg.norm(X, axis=0)) * .25 * R0 + .5 * R0
    ot = OptimalTransport(domain=domain, positions=X, masses=masses,
                          weights=np.zeros(len(masses)), verbosity=0, linear_solver='Scipy')
    ot.adjust_weights()
    COUNTERS.add(ot)
    return ot.get_centroids(), ot.pd.integrals()


def frank_wolfe(y, x0, lmo, dot_product, tolerance, t_max=100, callback=None):
    """Their frank_wolfe (no away steps). callback(t, xt, gap) after each
    gap evaluation; returning True stops."""
    xt = x0.copy()
    t = 0
    status = "max iterations"
    while t <= t_max:
        t += 1
        st = lmo(xt - y)
        dtFW = st - xt
        FWgap = dot_product(-(xt - y), dtFW)
        if callback is not None and callback(t - 1, xt, FWgap):
            status = "stopped"
            break
        if FWgap < tolerance:
            status = "converged"
            break
        gammat = np.max((np.min((dot_product(y - xt, dtFW) / dot_product(dtFW, dtFW), 1)), 0))
        xt = xt + gammat * dtFW
    return xt, t - 1, status


def convex_order_projection_frank_wolfe(domain, Y, X0, R0=1., masses=None, niter=600, tol=1e-5,
                                        callback=None):
    N = len(Y[:, 0])
    if masses is None:
        masses = 1 / N * np.ones(N)
    dot_product = lambda A, B: np.sum(A * B * masses.reshape((N, 1)))
    lmo = lambda X: get_centroids(domain, -X, masses, R0=R0)[0]
    return frank_wolfe(Y, X0, lmo, dot_product, tol, t_max=niter, callback=callback)
