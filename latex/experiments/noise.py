import sys
import numpy as np
sys.path.insert(0, "/tmp/exp2d")
from final import classify
from collapse import grid
rng = np.random.default_rng(11)
base = []
for t in range(6):
    s = rng.random((8, 2)); x, a = grid(2, 200)
    lab = ((x[None] - s[:, None])**2).sum(-1).argmin(0)
    base.append(np.array([x[lab == i].mean(0) for i in range(8)]))
for sig in (0.0, 0.002, 0.005, 0.01, 0.02):
    tags = []
    for cb in base:
        cbp = cb + sig * rng.standard_normal(cb.shape)
        tags.append(classify(80, cbp)[0])
    print(f"n=8 noise={sig:<6}", {k: tags.count(k) for k in ("REALIZED", "COLLAPSE", "UNBOUNDED", "OTHER")})
