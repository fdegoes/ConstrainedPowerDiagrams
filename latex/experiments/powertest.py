"""Rigorous (up to discretization) realizability test.
If some power diagram with distinct generators realizes the targets, its shifts u* maximize G and,
by complementary slackness, EVERY primal optimum theta equals its hard cells. So: take the LP primal
theta and test whether its (rounded) cells are separable by a power diagram, via a margin LP over
affine functions a_i(x) = 2 p_i.x - h_i. Margin > 0  => realized; margin <= 0 => no realizing diagram."""
import sys
import numpy as np
from scipy.optimize import linprog
from scipy.sparse import coo_matrix
sys.path.insert(0, "/tmp/exp2d")
from collapse import grid, solve

def power_margin(x, lab, n):
    K, d = x.shape
    nv = n * (d + 1) + 1  # (p_i, h_i) ..., t
    r, c, v = [], [], []
    row = 0
    for k in range(K):
        i = lab[k]
        for j in range(n):
            if j == i: continue
            # -(2(p_i-p_j).x - (h_i-h_j)) + t <= 0
            for e in range(d):
                r += [row, row]; c += [i*(d+1)+e, j*(d+1)+e]; v += [-2*x[k, e], 2*x[k, e]]
            r += [row, row, row]; c += [i*(d+1)+d, j*(d+1)+d, nv-1]; v += [1.0, -1.0, 1.0]
            row += 1
    A = coo_matrix((v, (r, c)), shape=(row, nv)).tocsr()
    cost = np.zeros(nv); cost[-1] = -1
    bounds = [(-10, 10)] * (nv - 1) + [(None, 1.0)]
    res = linprog(cost, A_ub=A, b_ub=np.zeros(row), bounds=bounds, method="highs")
    return -res.fun if res.status == 0 else np.nan

def test(name, m, cb):
    n = len(cb)
    x, a = grid(2, m)
    out = solve(x, a, cb.copy(), np.zeros(n), cb)
    if out is None:
        print(f"{name:34s} INFEASIBLE (unbounded G)"); return "INF"
    cost, theta, y = out
    lab = theta.argmax(0)
    hard = theta.max(0) > 1 - 1e-6
    vols = np.array([a[lab == i].sum() for i in range(n)])
    keep = hard & (vols[lab] > 0)
    nonempty = (vols > 0).sum()
    # relabel nonempty cells
    marg = power_margin(x[keep], lab[keep], n)
    tag = "REALIZED" if marg > 1e-7 else "NOT-POWER"
    print(f"{name:34s} {tag:9s} margin={marg:.2e} empty_cells={n-nonempty} frac_px={int((~hard).sum())}")
    return tag

if __name__ == "__main__":
    m = int(sys.argv[1]) if len(sys.argv) > 1 else 40
    rng = np.random.default_rng(1)
    summary = {}
    def rec(group, tag): summary.setdefault(group, []).append(tag)
    for n in (5, 8, 12):
        for t in range(4):
            s = rng.random((n, 2)); x, a = grid(2, m)
            lab = ((x[None] - s[:, None])**2).sum(-1).argmin(0)
            cb = np.array([x[lab == i].mean(0) for i in range(n)])
            rec(f"Voronoi-centroids n={n}", test(f"(a) Voronoi centroids n={n} t={t}", m, cb))
            for sig in (0.005, 0.02):
                cbp = np.clip(cb + sig*rng.standard_normal(cb.shape), 0.02, 0.98)
                rec(f"perturbed sig={sig} n={n}", test(f"(a') sig={sig} n={n} t={t}", m, cbp))
    for n in (5, 8, 12):
        for t in range(4):
            cb = 0.1 + 0.8*rng.random((n, 2))
            rec(f"random n={n}", test(f"(c) random n={n} t={t}", m, cb))
    print("\nSUMMARY")
    for g, tags in summary.items():
        print(f"{g:28s}", {k: tags.count(k) for k in set(tags)})
