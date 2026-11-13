# Realizable Bourne cases with our H (2026-11-12)

Our max H (`JointSolver`, from `MakeJointStart` with s = 0 and sites q = b̄)
on the three realizable cases of `bourne_cases.py`, under the common rule
at tolerances 1e-2, 1e-4, 1e-6 and 1e-8. Their methods are not run here.
Command:

```sh
$PY benchmarks/bourne/run_realizable.py --method H --tol 1e-2 1e-4 1e-6 1e-8 \
    --out benchmarks/bourne/results/realizable-H
```

Results in `results/realizable-H/`. Timings are ballpark: other jobs were
running (load average about 2.5).

## Convergence

| case | n | status | iterations to tol 1e-2 / 1e-4 / 1e-6 / 1e-8 | final volume / centroid error |
|---|---|---|---|---|
| 2025-invert | 20 | converged (18 it, 30 builds, 6 ms) | 6 / 9 / 11 / 13 | 8.7e-15 / 1.3e-15 |
| 2026-voronoi | 20 | converged (18 it, 30 builds, 5 ms) | 6 / 9 / 11 / 13 | 8.7e-15 / 1.3e-15 |
| 2026-aligned | 728 | **not converged** (4 attempts × 200 it, 989 builds, 3.3 s) | – / – / – / – | 1.7e-1 / 2.1e-2 |

- **2025-invert and 2026-voronoi have the same data** (Voronoi diagram of 20
  random seeds, seed 999). They differ only in their methods' starts. Ours
  starts from the targets in both, so the runs are identical. Under their
  rules: the SLSQP rule of 2025 is met at iteration 18. The FW gap rule of
  2026 (< 1e-14) is never met; the last gap is 1.5e-6, measured on our
  iterates.

## 2026-aligned: H fails on the cross

The data: equal volumes, generators on a regular 40 × 40 grid restricted to
a cross (56 of the 728 generators lie on the boundary of the square). The
tessellation is not irreducible (Bourne 2026, §3.4).

Our solve:

- The start (volume solve at q = b̄) converges, with no empty cell.
- Within about 10 iterations, pairs of generators come together: the
  smallest distance drops from 0.43 h to 0.01 h, and to 1e-8 h later.
- The residuals then oscillate around 0.1–1 for the rest of the
  iterations, with cells repeatedly shrinking to a few percent of V̄.
- It ends with 20 pairs closer than 0.01 h (the closest 7.6e-8 h).
- The pair regularization and the 3 restarts do not prevent it. Without
  them it fails earlier ("damping failed" after 115 iterations, residuals
  5.6e-2 / 6.8e-2).
- Starting with s = 0.5 or s = 1 gives the same result. Starting from the
  true generators fails at the start (the volume solve stops on a line
  search failure).

The other solvers on this case:

- Their FW, under their rule, also stops at its iteration cap (1000), with
  a gap of 3.4e-4 measured on our iterates (`results/realizable/summary.txt`).
- In their paper, the reconstructed generators are very different from the
  true ones, while the centroids fit (Figs. 11–12).

Our failure has the same root: the diagram is reducible, so the targets lie
on a face of positive dimension of the achievable set. Many generator
configurations, including ones with coincident generators, have nearly the
same centroids. The maximum of H is attained (the data are realizable),
but the iterates are drawn towards merging generators, which our
regularization was designed to prevent and here does not.

## For the paper

- The n = 20 Voronoi case is easy: 13 iterations to 1e-8, regardless of
  the start of their methods.
- The cross is a genuine failure of our H, shared in a different form by
  their FW. It should be reported, with the explanation (reducible
  diagram, generators not identifiable).
- This contrasts with the synthetic benchmark (598/600), whose generators
  are random, so their diagrams are generically irreducible (not checked case by case).
