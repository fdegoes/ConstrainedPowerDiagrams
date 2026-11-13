import sys
import numpy as np
sys.path.insert(0, "/tmp/exp2d")
from collapse import grid, solve, diagnose
rng = np.random.default_rng(11)
cases = []
for t in range(3):
    s = rng.random((8, 2)); x, a = grid(2, 200)
    lab = ((x[None] - s[:, None])**2).sum(-1).argmin(0)
    cb = np.array([x[lab == i].mean(0) for i in range(8)])
    cases.append((f"Voronoi t={t}", cb))
    for sig in (0.003, 0.01):
        cases.append((f"  noise {sig} t={t}", cb + sig*rng.standard_normal(cb.shape)))
cases.append(("random", 0.1 + 0.8*rng.random((8, 2))))
for name, cb in cases:
    row = []
    for m in (40, 80, 120):
        x, a = grid(2, m)
        out = solve(x, a, cb.copy(), np.zeros(8), cb)
        if out is None: row.append("unbounded"); continue
        dg = diagnose(x, a, cb.copy(), np.zeros(8), cb, *out)
        row.append(f"tie={dg['tie_area']:.3f} gd={dg['min_gen_dist']:.1e}")
    print(f"{name:20s}", " | ".join(row))
