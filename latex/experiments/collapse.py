"""Existence test for the centroid-constrained problem (Sec. 4.1) on a pixel grid.

The discretized dual max_u G(u) is an LP whose primal is the fractional-partition problem
  min sum_k a_k theta_ik c_ik  s.t.  sum_i theta_ik = 1,  sum_k a_k theta_ik (x_k - cbar_i) = 0.
Its centroid multipliers y_i equal 2 u_i. We check whether the hard power cells of the optimal
shifts realize the targets (no collapse) or not (collapse: coinciding affine functions)."""
import sys
import numpy as np
from scipy.optimize import linprog
from scipy.sparse import coo_matrix, vstack

def grid(dim, m):
    if dim == 1:
        x = (np.arange(m) + 0.5) / m
        return x[:, None], np.full(m, 1.0 / m)
    g = (np.arange(m) + 0.5) / m
    X, Y = np.meshgrid(g, g, indexing="ij")
    return np.c_[X.ravel(), Y.ravel()], np.full(m * m, 1.0 / m**2)

def solve(x, a, q, w, cb):
    n, d = cb.shape
    K = len(a)
    cost = ((x[None, :, :] - q[:, None, :]) ** 2).sum(-1) - w[:, None]  # n x K
    idx = np.arange(n * K).reshape(n, K)                                 # var (i,k)
    # partition-of-unity rows
    rows = np.tile(np.arange(K), n)
    B = coo_matrix((np.ones(n * K), (rows, idx.ravel())), shape=(K, n * K))
    # centroid rows: n*d
    r, c, v = [], [], []
    for i in range(n):
        for j in range(d):
            r.append(np.full(K, i * d + j)); c.append(idx[i]); v.append(a * (x[:, j] - cb[i, j]))
    A = coo_matrix((np.concatenate(v), (np.concatenate(r), np.concatenate(c))), shape=(n * d, n * K))
    Aeq = vstack([B, A]).tocsr()
    beq = np.r_[np.ones(K), np.zeros(n * d)]
    res = linprog((a[None, :] * cost).ravel(), A_eq=Aeq, b_eq=beq, bounds=(0, None), method="highs")
    if res.status != 0:
        return None
    theta = res.x.reshape(n, K)
    y = res.eqlin.marginals[K:].reshape(n, d)
    return cost, theta, y

def diagnose(x, a, q, w, cb, cost, theta, y):
    n, d = cb.shape
    # shifted power distance per pixel (per unit area), for both sign conventions of y; keep the
    # one consistent with complementary slackness (support of theta at the minimum)
    best = None
    for s in (+1, -1):
        f = cost - s * np.einsum("id,kd->ik", y, x) + s * (y * cb).sum(1)[:, None]
        gap = (f * (theta > 1e-9)).sum() - (f.min(0) * 1).sum() * 0  # placeholder
        viol = ((f - f.min(0)) * theta).sum()
        if best is None or viol < best[0]:
            best = (viol, s, f)
    viol, s, f = best
    u = s * y / 2.0
    lab = f.argmin(0)
    srt = np.sort(f, 0)
    tie = (srt[1] - srt[0]) < 1e-6 * (np.abs(srt[0]).max() + 1)
    hard_res, vols = [], []
    for i in range(n):
        msk = lab == i
        V = a[msk].sum()
        vols.append(V)
        c = (a[msk, None] * x[msk]).sum(0) / V if V > 0 else np.full(d, np.nan)
        hard_res.append(np.linalg.norm(c - cb[i]))
    # lifted generators (p_i, h_i): a_i(x) = 2 p.x - h ; from f = |x|^2 - a_i(x) up to a common term
    p = q + u
    lam = w + (u**2).sum(1) + 2 * (u * (q - cb)).sum(1)
    h = (p**2).sum(1) - lam
    G = np.c_[p, h]
    dmin = min(np.linalg.norm(G[i] - G[j]) for i in range(n) for j in range(i + 1, n))
    frac = ((theta > 1e-6) & (theta < 1 - 1e-6)).any(0).sum()
    return dict(cs_viol=viol, hard_res=np.nanmax(hard_res), tie_area=a[tie].sum(),
                min_gen_dist=dmin, frac_pixels=int(frac), min_vol=min(vols), unorm=np.abs(u).max())

def run(name, dim, m, q, w, cb):
    x, a = grid(dim, m)
    out = solve(x, a, q, w, cb)
    if out is None:
        print(f"{name:38s} INFEASIBLE"); return
    dg = diagnose(x, a, q, w, cb, *out)
    tag = "REALIZED" if dg["hard_res"] < 3.0 / m and dg["tie_area"] < 5.0 / m else "COLLAPSE"
    print(f"{name:38s} {tag:9s} " + " ".join(f"{k}={v:.3g}" for k, v in dg.items()))
    return tag

if __name__ == "__main__":
    rng = np.random.default_rng(0)
    # 1D sanity check: c = (0.3, 0.7) should collapse, (0.25, 0.75) realized
    for cb in ([0.3, 0.7], [0.25, 0.75], [0.2, 0.6]):
        cb = np.array(cb)[:, None]
        run(f"1D cbar={cb.ravel()}", 1, 2000, cb.copy(), np.zeros(2), cb)
    m = int(sys.argv[1]) if len(sys.argv) > 1 else 80
    print("--- 2D, centroidal setting q = cbar, w = 0, grid", m)
    # (a) centroids of a Voronoi diagram of random points (realizable by construction w/ q=sites, not q=cbar)
    for n in (4, 8, 16):
        for trial in range(3):
            s = rng.random((n, 2))
            x, a = grid(2, m)
            lab = ((x[None] - s[:, None]) ** 2).sum(-1).argmin(0)
            cb = np.array([x[lab == i].mean(0) for i in range(n)])
            run(f"(a) Voronoi centroids n={n} t={trial}", 2, m, cb.copy(), np.zeros(n), cb)
            for sig in (0.02, 0.05):
                cbp = np.clip(cb + sig * rng.standard_normal(cb.shape), 0.02, 0.98)
                run(f"(a') perturbed sig={sig} n={n} t={trial}", 2, m, cbp.copy(), np.zeros(n), cbp)
    # (b) two nearby targets among a jittered grid
    base = np.array([[i, j] for i in range(3) for j in range(3)], float) / 3 + 1 / 6
    for delta in (0.0, 0.01, 0.03, 0.06):
        cb = base.copy(); cb[4] = [0.5 - delta / 2, 0.5]; cb = np.vstack([cb, [0.5 + delta / 2, 0.5]])
        run(f"(b) close pair delta={delta}", 2, m, cb.copy(), np.zeros(len(cb)), cb)
    # (c) uniform random targets
    for n in (4, 8, 16):
        for trial in range(4):
            cb = 0.1 + 0.8 * rng.random((n, 2))
            run(f"(c) random n={n} t={trial}", 2, m, cb.copy(), np.zeros(n), cb)
