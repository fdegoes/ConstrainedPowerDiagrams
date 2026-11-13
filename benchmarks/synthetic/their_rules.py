"""The stopping rules of Bourne et al., applied to the iterates of any method.

Neither paper stops on the centroid residual of the solution:

  SLSQP on H (2025, max_H_solver): scipy's SLSQP stops when the change of the
      objective between iterates is below ftol = ptol |H(X_init)|, with their
      realizable setting ptol = 1e-10 (1e-8 here, SLSQP_PTOL), and its inner OT solves meet the volumes
      to max |V - Vbar| < 1e-3 min Vbar. Here: the first iterate k >= 1 with
      max |V - Vbar| < 1e-3 min Vbar and |H_k - H_(k-1)| <= ptol |H_0|, with
      their objective H(X) = sum Vbar (b - c) . X at generators X and
      centroids c, and H_0 at our start.

  Frank-Wolfe (2026, convex_order_projection_frank_wolfe): stops when the FW
      gap sum Vbar (y - c) . (c' - c) is below tol (their realizable settings
      1e-14 and 1e-16; 1e-14 here), with y = t / (t - 1) b, t = 100, and c'
      the centroids of the linear minimization oracle at c (one OT solve, not
      charged to the method). The gap is only evaluated at iterates with the
      prescribed volumes (relative error < 1e-8), the points of the set FW
      works on.
"""

import numpy as np

SLSQP_PTOL = 1e-8  # their 2025 notebook: 1e-10
SLSQP_VOLUME_TOL = 1e-3
FW_GAP = 1e-14
FW_T = 100.0
FW_VOLUME_TOL = 1e-8


def their_H(b, v, V, c, X):
    return float(np.sum(v * np.sum((b - c) * X, axis=1)))


def slsqp_stop(b, v, traj, ptol=SLSQP_PTOL):
    """First iterate meeting the SLSQP rule (dict) or None."""
    H0 = their_H(b, v, np.asarray(traj[0]["volumes"]), np.reshape(traj[0]["centroids"], (-1, b.shape[1])),
                 np.reshape(traj[0]["generators"], (-1, b.shape[1])))
    ftol = ptol * abs(H0)
    prev = None
    for k, it in enumerate(traj):
        V, c = np.asarray(it["volumes"]), np.reshape(it["centroids"], (-1, b.shape[1]))
        H = their_H(b, v, V, c, np.reshape(it["generators"], (-1, b.shape[1])))
        if prev is not None and np.max(np.abs(V - v)) < SLSQP_VOLUME_TOL * np.min(v) and abs(H - prev) <= ftol:
            return dict(iteration=k, seconds=it["seconds"], builds=it["numBuilds"],
                        volumeResidual=it["volumeResidual"], shiftResidual=it["shiftResidual"])
        prev = H
    return None


def slsqp_at_solution(b, v, traj, ptol=SLSQP_PTOL):
    """Rule 1 measured at the solution (the final iterate) of a run: the
    last objective change |H_K - H_(K-1)| relative to |H_0|, and whether it
    and the volume condition pass (dict), or None without two iterates."""
    if len(traj) < 2:
        return None
    d = b.shape[1]
    H = [their_H(b, v, np.asarray(it["volumes"]), np.reshape(it["centroids"], (-1, d)),
                 np.reshape(it["generators"], (-1, d))) for it in (traj[0], traj[-2], traj[-1])]
    ratio = abs(H[2] - H[1]) / abs(H[0])
    volOk = bool(np.max(np.abs(np.asarray(traj[-1]["volumes"]) - v)) < SLSQP_VOLUME_TOL * np.min(v))
    return dict(ratio=float(ratio), met=bool(volOk and ratio <= ptol))


def fw_gap(om, b, v, c, R0=1.0, t=FW_T):
    """Their FW gap at centroids c (one LMO, uncharged)."""
    import bourne_ref as ref
    y = t / (t - 1) * b
    with ref.Evaluation():
        cp, _ = ref.get_centroids(om, -(c - y), v, R0=R0)
    return float(np.sum(v[:, None] * (y - c) * (cp - c)))


def fw_stop(om, b, v, traj, gap=FW_GAP, stride=1):
    """First iterate meeting the FW-gap rule (dict) or None. stride: evaluate
    the gap every stride iterates (it costs an OT solve)."""
    for k, it in enumerate(traj):
        if it["volumeResidual"] > FW_VOLUME_TOL or (k % stride and k != len(traj) - 1):
            continue
        g = fw_gap(om, b, v, np.reshape(it["centroids"], (-1, b.shape[1])))
        if abs(g) < gap:
            return dict(iteration=k, seconds=it["seconds"], builds=it["numBuilds"], gap=g,
                        volumeResidual=it["volumeResidual"], shiftResidual=it["shiftResidual"])
    return None


def centroid_error(v, V, c, b):
    """max_i |b_i - c_i| / Vbar_i^(1/2): the centroid error in cell sizes (2D)."""
    return float(np.max(np.linalg.norm(np.asarray(c) - b, axis=1) / np.sqrt(v)))

