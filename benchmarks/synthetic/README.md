# Synthetic benchmark

These are inverse problems in the unit square with a known solution. Each case
gives target volumes V̄ and centroids c̄. The solver has to find generators and
weights whose power diagram has those volumes and centroids. Every case is
built from a power diagram we already have, so a solution always exists and we
know it.

## Files

- `generate.py`: generates the cases (`make_case`).
- `makeSynthetic.cpp`: computes the power diagram of Y with the given volumes
  and writes its volumes, centroids and weights. Target `makeSynthetic` in
  `CMakeLists.txt`.
- `run_synthetic.py`: runs our solver (`benchBourne`, max H) on every case and
  writes `<out>/ours.csv`.
- `run_synthetic_g.py`: runs the centroid dual G (`benchBourne` with
  `solver G`, max over the shifts with fixed weights) on every case and
  writes `<out>/oursG.csv`. See "Max G with fixed weights" below.
- `run_theirs.py`: runs the methods of Bourne et al. (Frank–Wolfe, SLSQP) on
  every case and writes `<out>/theirs.csv`.
- `diagnose.py`: explains why a case fails to converge.
- `stats.py`: generates all the cases and prints the statistics below.
- `report.py`: the tables of `REPORT.md` (H vs. G) from `results/ours.csv`
  and `results/oursG.csv`.
- `REPORT.md`: results of H and G on the 600 2D cases.
- `gallery.py`: draws the n = 10 diagrams (seed 1), with weight disks, for each noise level and
  volume distribution, as PDFs in `latex/figures/synthetic/`.
- `render_all.py`: draws the diagram of every case, for browsing (not used in the paper),
  in `images/noise<k>-<volumes>/n<n>-seed<seed>.pdf`.

## How a case is generated

1. **Generators Y.** Draw n points by dart throwing in [0,1]². Candidates are
   uniform in the square. A candidate is rejected if it falls within ρ·h of an
   accepted point, where h = √(1/n) is the cell size. The radius ρ sets the
   noise level: the smaller ρ, the less uniform the spacing.

   | noise | meaning | ρ |
   |---|---|---|
   | 0 | low noise | 0.8 |
   | 1 | mid noise | 0.5 |
   | 2 | high noise | 0.25 |

   The random generator is seeded with `[seed, n, noise]`, so the uniform and
   voronoi cases of one (n, noise, seed) share the same Y.
2. **Target volumes V̄.**
   - `uniform`: V̄ᵢ = 1/n.
   - `voronoi`: V̄ᵢ is the area of Yᵢ's Voronoi cell. The weights of the
     solution are then zero.
3. **The solution diagram.** `makeSynthetic` runs `VolumeSolver` to find the
   weights w for which the power diagram of Y has volumes V̄ (tolerance 1e-13).
   It measures that diagram's volumes and centroids c̄. These realized volumes
   and c̄ become the targets, so both come from the same diagram exactly.
4. **Start.** The solvers start from the sites q = c̄ and have to recover Y. Y
   and w are determined only up to a homothety about the centroid, which H is
   invariant along. `run_synthetic.py` measures `generatorError` with that
   homothety removed, and `diagramMismatch` (the area assigned to the wrong
   cell).

Each case is stored in the directory `<out>/noise<k>-<volumes>-<n>-<seed>/`,
which contains:

- `generators.txt`: n, then Y (the input of `makeSynthetic`).
- `data.txt`: n, the volumes, the targets c̄ and the weights.
- `truth.json`: Y, the weights, the volumes and the centroids.

## Grid

- n: 20 sizes, log-spaced in [10, 10000] as round(10 · 1000^(k/19)), k = 0…19.
  The sizes are 10, 14, 21, 30, 43, 62, 89, 127, 183, 264, 379, 546, 785,
  1129, 1624, 2336, 3360, 4833, 6952 and 10000 (`NS` in `run_synthetic.py`).
- noise: 0, 1, 2.
- volumes: uniform, voronoi.
- seeds: 1–5.

That is 20 × 3 × 2 × 5 = **600 cases**. Generating all of them takes about
2.5 minutes. Dart throwing at ρ = 0.8 is the slowest step, about 13 s per case
at n = 10000.

## Usage

```sh
cmake --build build --target makeSynthetic benchBourne
external/bourne/.venv/bin/python benchmarks/synthetic/stats.py       # generate + statistics
external/bourne/.venv/bin/python benchmarks/synthetic/run_synthetic.py [--n ...] [--noise 0 1 2] [--volumes uniform voronoi] [--seeds ...]
external/bourne/.venv/bin/python benchmarks/synthetic/run_synthetic_g.py [same options] [--dim 3]
external/bourne/.venv/bin/python benchmarks/synthetic/run_theirs.py    [same options] [--method fw slsqp]
```

The output directory defaults to `benchmarks/synthetic/results`. Cases that
already exist are reused.

## Max G with fixed weights

`run_synthetic_g.py` solves the same cases with the centroid dual G
(`CentroidSolver`), whose unknowns are the shifts only. The sites are the
targets, q = c̄, and the weights are fixed so that the shifted weights of
the solution are the weights w_Y of the synthetic diagram. With q = c̄,
eq:shifted-site-weight is λ_i = w_i + |u_i|², and the solution has
u_i = Y_i − c̄_i, so

    w_i = w_Y,i − |Y_i − c̄_i|².

With fixed weights G has no homothety invariance, so Y and w_Y are its
unique solution. The errors are measured directly:

- `generatorError`: max |p − Y| / h.
- `weightError`: max |λ − w_Y| / h², both with their mean removed (a
  constant does not change the diagram).
- `volumeResidual`: max |V − V̄| / V̄ on the final cells (G does not
  constrain the volumes; they match at the solution).
- `startEmptyCells`: empty cells of the start, the diagram of (c̄, w),
  counted exactly by `benchBourne`.
- `outcome`: `converged`; `invalid start` (the start has empty cells, so
  G cannot take a step: an empty cell has no gradient, and no step that
  keeps it empty is accepted); or `failed` (valid start, no convergence).

The solve starts from zero shifts and stops on the shift residual (1e-8).
Results go to `<case>/oursG.json` and `<out>/oursG.csv`.

The results of H and G on all 600 2D cases, and the key findings, are in
`REPORT.md` (tables regenerated by `report.py`). In short: H converges on
598/600; G on 523/600, with 62 invalid starts (uniform volumes at noise
1–2) and 15 stalls where a cell collapses during the solve.

## Statistics of the data

Columns:

- **nn min**: the smallest distance from a generator to its nearest
  neighbour, over all cases in the row, in units of h.
- **nn cv**: the coefficient of variation (std / mean) of the
  nearest-neighbour distances within one case.
- **V max/min**: the ratio of the largest to the smallest cell volume.
- **|w| / h²**: the largest weight of the solution, after shifting the
  weights to zero mean.
- **|Y − c̄| / h**: the distance from the start c̄ to the solution Y, as the
  mean or the maximum over the cells of one case.

Ranges are min–max over the cases in the row. To regenerate the tables, run
`stats.py`.

### Per noise level and volume distribution (pooled over n and seeds)

| noise, volumes | cases | nn min | nn cv | V max/min | \|w\| / h^2 | mean \|Y - cbar\| / h | max \|Y - cbar\| / h |
|---|---|---|---|---|---|---|---|
| 0 (rho 0.80), uniform | 100 | 0.800 | 0.06 - 0.19 | 1.00 - 1.00 | 0.552 - 37.724 | 0.192 - 0.441 | 0.412 - 1.048 |
| 0 (rho 0.80), voronoi | 100 | 0.800 | 0.06 - 0.19 | 2.34 - 10.13 | 0.000 - 0.000 | 0.069 - 0.207 | 0.266 - 0.494 |
| 1 (rho 0.50), uniform | 100 | 0.500 | 0.10 - 0.43 | 1.00 - 1.00 | 0.328 - 66.413 | 0.189 - 0.698 | 0.316 - 2.244 |
| 1 (rho 0.50), voronoi | 100 | 0.500 | 0.10 - 0.43 | 1.77 - 24.28 | 0.000 - 0.000 | 0.149 - 0.279 | 0.288 - 0.861 |
| 2 (rho 0.25), uniform | 100 | 0.250 | 0.24 - 0.63 | 1.00 - 1.00 | 0.776 - 142.068 | 0.320 - 1.157 | 0.631 - 3.338 |
| 2 (rho 0.25), voronoi | 100 | 0.250 | 0.24 - 0.63 | 4.23 - 41.38 | 0.000 - 0.000 | 0.198 - 0.331 | 0.363 - 1.077 |

### Per n (pooled over noise, volumes and seeds)

| n | cases | nn min | nn cv | V max/min | \|w\| / h^2 | mean \|Y - cbar\| / h | max \|Y - cbar\| / h |
|---|---|---|---|---|---|---|---|
| 10 | 30 | 0.296 | 0.06 - 0.63 | 1.00 - 11.43 | 0.000 - 2.141 | 0.151 - 0.570 | 0.288 - 1.286 |
| 14 | 30 | 0.253 | 0.06 - 0.52 | 1.00 - 9.75 | 0.000 - 2.132 | 0.121 - 0.618 | 0.266 - 1.494 |
| 21 | 30 | 0.253 | 0.09 - 0.53 | 1.00 - 23.91 | 0.000 - 3.405 | 0.115 - 0.776 | 0.292 - 1.346 |
| 30 | 30 | 0.255 | 0.11 - 0.53 | 1.00 - 13.02 | 0.000 - 3.196 | 0.128 - 0.602 | 0.270 - 1.592 |
| 43 | 30 | 0.257 | 0.08 - 0.56 | 1.00 - 18.96 | 0.000 - 5.713 | 0.115 - 0.689 | 0.275 - 1.721 |
| 62 | 30 | 0.252 | 0.08 - 0.48 | 1.00 - 13.10 | 0.000 - 9.145 | 0.103 - 0.925 | 0.276 - 1.751 |
| 89 | 30 | 0.250 | 0.08 - 0.48 | 1.00 - 25.99 | 0.000 - 11.460 | 0.096 - 0.862 | 0.266 - 1.781 |
| 127 | 30 | 0.252 | 0.08 - 0.45 | 1.00 - 24.52 | 0.000 - 10.331 | 0.093 - 0.727 | 0.296 - 1.886 |
| 183 | 30 | 0.250 | 0.08 - 0.44 | 1.00 - 23.16 | 0.000 - 11.676 | 0.088 - 0.721 | 0.279 - 1.955 |
| 264 | 30 | 0.250 | 0.09 - 0.44 | 1.00 - 24.36 | 0.000 - 13.708 | 0.081 - 0.668 | 0.295 - 1.930 |
| 379 | 30 | 0.250 | 0.08 - 0.44 | 1.00 - 29.00 | 0.000 - 26.867 | 0.078 - 0.925 | 0.305 - 2.333 |
| 546 | 30 | 0.251 | 0.07 - 0.39 | 1.00 - 22.21 | 0.000 - 13.104 | 0.081 - 0.630 | 0.268 - 1.975 |
| 785 | 30 | 0.250 | 0.07 - 0.43 | 1.00 - 24.42 | 0.000 - 25.208 | 0.076 - 0.820 | 0.307 - 2.303 |
| 1129 | 30 | 0.250 | 0.07 - 0.40 | 1.00 - 32.80 | 0.000 - 48.165 | 0.074 - 1.019 | 0.273 - 2.647 |
| 1624 | 30 | 0.250 | 0.08 - 0.40 | 1.00 - 25.99 | 0.000 - 35.367 | 0.074 - 0.883 | 0.285 - 2.451 |
| 2336 | 30 | 0.250 | 0.07 - 0.41 | 1.00 - 41.38 | 0.000 - 42.868 | 0.071 - 0.903 | 0.309 - 2.468 |
| 3360 | 30 | 0.250 | 0.07 - 0.40 | 1.00 - 30.92 | 0.000 - 75.843 | 0.072 - 1.025 | 0.305 - 2.934 |
| 4833 | 30 | 0.250 | 0.07 - 0.40 | 1.00 - 32.42 | 0.000 - 70.935 | 0.070 - 0.834 | 0.322 - 2.822 |
| 6952 | 30 | 0.250 | 0.07 - 0.40 | 1.00 - 39.19 | 0.000 - 142.068 | 0.070 - 1.157 | 0.327 - 2.805 |
| 10000 | 30 | 0.250 | 0.07 - 0.39 | 1.00 - 38.12 | 0.000 - 102.758 | 0.069 - 0.891 | 0.308 - 3.338 |

### Observations

- **The noise levels work as intended.** The minimum spacing equals ρ·h at
  every level. The spread of the spacing (nn cv) grows with the level:
  0.06–0.19, 0.10–0.43 and 0.24–0.63. The ranges overlap a little, mostly at
  small n, where a case has few distances to measure.
- **Uniform volumes** make all cells the same size. The weights of these cases
  grow with the noise and with n (up to 142 h² at noise 2): the less evenly
  spaced Y is, the further its equal-volume diagram is from its Voronoi
  diagram. The start c̄ is also further from Y (max 1.0, 2.2 and 3.3 h).
- **Voronoi volumes** have zero weights. The cell volumes vary by a factor of
  2–10 at noise 0 and up to 41 at noise 2. c̄ stays within about 1 h of Y,
  because these diagrams are ordinary Voronoi diagrams.
- The hardest cases should be uniform volumes at noise 2 and large n. They have the largest weights and the farthest start.

## 3D cases

The same grid in the unit cube [0,1]³, in `results3d/` with the same case
names. Generate it with `stats.py --dim 3`, which takes about 4 minutes.

- **Cell size:** h = (1/n)^(1/3), and dart throwing uses the minimum distance
  ρ·h with the same three values of ρ.
- **Seeds:** the random generator is seeded with `[seed, n, noise, 3]`, so the
  3D generators are independent of the 2D ones.
- **Files:** `generators.txt` starts with `n 3`, and `makeSynthetic` writes
  a `dim 3` record into `data.txt`. 2D files are unchanged.
- **Domain:** the unit cube, tetrahedralized as in `tests/testUtils.h`.
- **Check:** a Monte Carlo check at n = 20 (4 × 10⁶ samples) reproduced the
  stored volumes and centroids to within the sampling error.

To run our solver on them, use `run_synthetic.py --dim 3`, which writes
`results3d/ours.csv`, or `diagnose.py --dim 3 <case>`. `benchBourne` takes a
`dim 3` key and a 6-value box. `run_theirs.py`, `render_all.py` and
`gallery.py` stay 2D, because Bourne's methods and the drawings are 2D.

### Per noise level and volume distribution (pooled over n and seeds)

| noise, volumes | cases | nn min | nn cv | V max/min | \|w\| / h^2 | mean \|Y - cbar\| / h | max \|Y - cbar\| / h |
|---|---|---|---|---|---|---|---|
| 0 (rho 0.80), uniform | 100 | 0.800 | 0.04 - 0.22 | 1.00 - 1.00 | 0.382 - 6.621 | 0.211 - 0.426 | 0.429 - 0.811 |
| 0 (rho 0.80), voronoi | 100 | 0.800 | 0.04 - 0.22 | 2.27 - 13.91 | 0.000 - 0.000 | 0.094 - 0.290 | 0.318 - 0.543 |
| 1 (rho 0.50), uniform | 100 | 0.500 | 0.12 - 0.30 | 1.00 - 1.00 | 0.323 - 6.472 | 0.277 - 0.634 | 0.472 - 1.283 |
| 1 (rho 0.50), voronoi | 100 | 0.500 | 0.12 - 0.30 | 2.40 - 29.59 | 0.000 - 0.000 | 0.163 - 0.403 | 0.349 - 0.873 |
| 2 (rho 0.25), uniform | 100 | 0.250 | 0.14 - 0.46 | 1.00 - 1.00 | 0.400 - 10.067 | 0.312 - 0.679 | 0.501 - 1.435 |
| 2 (rho 0.25), voronoi | 100 | 0.250 | 0.14 - 0.46 | 2.53 - 73.04 | 0.000 - 0.000 | 0.199 - 0.360 | 0.456 - 0.800 |

### Per n (pooled over noise, volumes and seeds)

| n | cases | nn min | nn cv | V max/min | \|w\| / h^2 | mean \|Y - cbar\| / h | max \|Y - cbar\| / h |
|---|---|---|---|---|---|---|---|
| 10 | 30 | 0.283 | 0.06 - 0.46 | 1.00 - 10.63 | 0.000 - 1.125 | 0.204 - 0.634 | 0.344 - 0.886 |
| 14 | 30 | 0.318 | 0.07 - 0.40 | 1.00 - 22.98 | 0.000 - 1.295 | 0.177 - 0.679 | 0.340 - 1.198 |
| 21 | 30 | 0.349 | 0.06 - 0.43 | 1.00 - 16.25 | 0.000 - 1.499 | 0.204 - 0.585 | 0.339 - 0.886 |
| 30 | 30 | 0.254 | 0.04 - 0.44 | 1.00 - 17.85 | 0.000 - 1.908 | 0.160 - 0.494 | 0.323 - 0.901 |
| 43 | 30 | 0.255 | 0.08 - 0.43 | 1.00 - 40.57 | 0.000 - 3.747 | 0.164 - 0.643 | 0.318 - 1.414 |
| 62 | 30 | 0.252 | 0.07 - 0.42 | 1.00 - 14.13 | 0.000 - 2.042 | 0.166 - 0.413 | 0.386 - 1.142 |
| 89 | 30 | 0.252 | 0.08 - 0.36 | 1.00 - 16.70 | 0.000 - 2.967 | 0.156 - 0.491 | 0.369 - 1.153 |
| 127 | 30 | 0.253 | 0.09 - 0.40 | 1.00 - 16.07 | 0.000 - 3.564 | 0.145 - 0.532 | 0.375 - 1.076 |
| 183 | 30 | 0.253 | 0.08 - 0.38 | 1.00 - 20.02 | 0.000 - 3.092 | 0.131 - 0.463 | 0.354 - 1.053 |
| 264 | 30 | 0.250 | 0.08 - 0.37 | 1.00 - 25.82 | 0.000 - 4.758 | 0.129 - 0.530 | 0.356 - 1.385 |
| 379 | 30 | 0.251 | 0.07 - 0.37 | 1.00 - 24.19 | 0.000 - 4.095 | 0.127 - 0.474 | 0.358 - 1.234 |
| 546 | 30 | 0.252 | 0.07 - 0.35 | 1.00 - 25.13 | 0.000 - 5.932 | 0.121 - 0.495 | 0.352 - 1.296 |
| 785 | 30 | 0.251 | 0.07 - 0.34 | 1.00 - 36.74 | 0.000 - 5.521 | 0.118 - 0.439 | 0.346 - 1.183 |
| 1129 | 30 | 0.250 | 0.07 - 0.34 | 1.00 - 40.01 | 0.000 - 4.678 | 0.113 - 0.435 | 0.397 - 1.435 |
| 1624 | 30 | 0.250 | 0.07 - 0.35 | 1.00 - 42.60 | 0.000 - 6.386 | 0.110 - 0.440 | 0.392 - 1.415 |
| 2336 | 30 | 0.250 | 0.07 - 0.34 | 1.00 - 34.18 | 0.000 - 6.658 | 0.107 - 0.439 | 0.378 - 1.275 |
| 3360 | 30 | 0.250 | 0.06 - 0.33 | 1.00 - 45.93 | 0.000 - 6.419 | 0.102 - 0.414 | 0.367 - 1.305 |
| 4833 | 30 | 0.250 | 0.06 - 0.33 | 1.00 - 45.07 | 0.000 - 8.830 | 0.100 - 0.415 | 0.370 - 1.273 |
| 6952 | 30 | 0.250 | 0.06 - 0.33 | 1.00 - 73.04 | 0.000 - 10.067 | 0.098 - 0.425 | 0.383 - 1.321 |
| 10000 | 30 | 0.250 | 0.06 - 0.33 | 1.00 - 64.29 | 0.000 - 9.345 | 0.094 - 0.407 | 0.359 - 1.314 |

Compared with 2D, the weights are smaller (up to 10 h² against 142 h²), the
start is closer to the solution (at most 1.4 h against 3.3 h), and the Voronoi
volume ratios are larger (up to 73 against 41).

