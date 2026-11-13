"""Classify the maximizer of the discretized centroid dual (q = cbar, w = 0, unit square):
  REALIZED : distinct lifted generators, negligible tie area, hard power cells hit the targets
  COLLAPSE : some lifted generators coincide (tie region of positive area)
  UNBOUNDED: LP infeasible (0 outside the residual set R)"""
import sys
import numpy as np
sys.path.insert(0, "/tmp/exp2d")
from collapse import grid, solve, diagnose

def classify(m, cb):
    n = len(cb); x, a = grid(2, m)
    out = solve(x, a, cb.copy(), np.zeros(n), cb)
    if out is None: return "UNBOUNDED", None
    dg = diagnose(x, a, cb.copy(), np.zeros(n), cb, *out)
    if dg["min_gen_dist"] < 1e-8 and dg["tie_area"] > 0.02: return "COLLAPSE", dg
    if dg["hard_res"] < 3.0 / m: return "REALIZED", dg
    return "OTHER", dg

if __name__ == "__main__":
    rng = np.random.default_rng(7)
    groups = {}
    for n in (5, 8, 12):
        for t in range(6):
            s = rng.random((n, 2)); x, a = grid(2, 200)
            lab = ((x[None] - s[:, None])**2).sum(-1).argmin(0)
            cb = np.array([x[lab == i].mean(0) for i in range(n)])
            groups.setdefault(f"Voronoi centroids n={n}", []).append(cb)
            groups.setdefault(f"  + noise 0.03  n={n}", []).append(np.clip(cb + 0.03*rng.standard_normal(cb.shape), .03, .97))
            groups.setdefault(f"uniform random  n={n}", []).append(0.1 + 0.8*rng.random((n, 2)))
    for m in (40, 80):
        print(f"grid {m}x{m}")
        for g, cbs in groups.items():
            tags = [classify(m, cb)[0] for cb in cbs]
            print(f"  {g:26s}", {k: tags.count(k) for k in ("REALIZED", "COLLAPSE", "UNBOUNDED", "OTHER")})
