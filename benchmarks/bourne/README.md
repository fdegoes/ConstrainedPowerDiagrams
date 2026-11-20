# Benchmarks against Bourne et al. (2025, 2026)

Reproduces the experiments of

- Bourne, Pearce & Roper, *Inverting Laguerre tessellations* (M2AN 2025) — `Bourne2025Inverting`
- Bourne, Gallouët, Mérigot & Natale, *Semi-discrete convex order and Laguerre tessellation fitting* (2026) — `Bourne2026ConvexOrder`

with their methods (SLSQP on their H or least squares; Frank–Wolfe) and with
the dual solvers of this project (max H, max H_alpha). Each method stops on
its own rule (their SLSQP rule, or a common per-cell L∞ residual), and the
other rule is measured at its solution.

## Files

- `run_realizable.py` — driver for the realizable cases: their method and
  our H.
- `run_approximate.py` — driver for the approximate cases: their method and
  our H_alpha over a sweep of stiffnesses.
- `run_comparison.py` — head-to-head on the approximate cases (see
  `README_COMPARISON.md`): their
  methods and our H_alpha, at our alpha bound (eq:regularized-distinct), at
  our practical estimate of it, and at their hand-picked value (EBSD t = 10), each run to convergence or to the
  caps (1000 iterations, 1800 s); reports the stop reason, iterations, time
  and the volume, centroid and stationarity residuals (`results.csv`).
- `bourne_cases.py` — shared by both drivers: the cases, the runs of every
  method, the termination rules, and the summary.
- `bourne_ref.py` — their methods, ported from their code with
  instrumentation (OT-solve counters, per-iterate callbacks). Sources in
  `external/bourne/`.
- `../common/boxDomain.h` — the box domain (triangulated in 2D, tetrahedralized
  in 3D) shared by `benchBourne`, `checkHessian` and `makeSynthetic`.
- `benchBourne.cpp` — one dual solve, 2D or 3D (`dim` key, default 2)
  (`JointSolver` / `SoftenedSolver` from `MakeJointStart`, or
  `CentroidSolver` for `solver G`, from zero shifts with the fixed
  `weights` of the case file), writes the result and the trajectory as JSON. Target
  `benchBourne` in `CMakeLists.txt`.

## Setup

```sh
cmake --build build --target benchBourne
# Python with pysdot (built against RV64's Boost, see external/bourne/src/pysdot-0.2.39)
external/bourne/.venv/bin/python benchmarks/bourne/run_realizable.py --list
external/bourne/.venv/bin/python benchmarks/bourne/run_approximate.py --list
```

`external/bourne/` is not in git. To rebuild it:

- Their code, unmodified (`bourne_ref.py` ports it and adds the
  instrumentation):
  - `Laguerre-Polycrystalline-Microstructures`:
    https://github.com/DPBourne/Laguerre-Polycrystalline-Microstructures at
    `691aeedd94e6461ea876c86d48980579bfff12e9`
  - `Laguerre-fitting`: https://github.com/andnatale/Laguerre-fitting at
    `675360d109d98084d16c98bb3c9185ba1f161658`
- `src/pysdot-0.2.39`: the PyPI sdist of pysdot 0.2.39 with
  `pysdot-0.2.39.patch` (one Boost define in `setup.py`, for recent libc++).
  Its `setup.py` clones `ext/sdot` (`1b718ef`), `ext/pybind11` (`3beb376`) and
  `ext/eigen3` (`36b9596`) and downloads Boost 1.87.0.
- `.venv`: Python 3.13.3 with numpy 2.5.4, scipy 1.18.1, matplotlib 3.11.2,
  pandas 3.1.1, vtk 9.7.2, h5py 3.16.0, and pysdot installed from
  `src/pysdot-0.2.39`.

## Usage

Each driver runs everything in its group by default (all cases, methods,
rules, tolerances):

```sh
MPLCONFIGDIR=/tmp/mplcfg external/bourne/.venv/bin/python benchmarks/bourne/run_realizable.py
MPLCONFIGDIR=/tmp/mplcfg external/bourne/.venv/bin/python benchmarks/bourne/run_approximate.py
```

Options (both drivers, except where noted):

| option | default | meaning |
|---|---|---|
| `--case NAME...` | all in the driver's group | realizable: `2025-invert`, `2026-voronoi`, `2026-aligned`; approximate: `2025-perturbed-small`, `2025-perturbed-large`, `ebsd` |
| `--method M...` | `theirs H` / `theirs Halpha` | their method of the case's paper; our H (realizable); our H_alpha (approximate) |
| `--tol X...` | `1e-8` | tolerances of the common rule (the smallest stops the runs; the larger ones are read off the trajectory) |
| `--alpha A...` | `10 100 1000` | approximate only: H_alpha stiffnesses |
| `--alpha-bound R...` | `0.99 10 100` | approximate only: further stiffnesses, as multiples of the distinctness bound eq:regularized-distinct (`bourne_cases.distinct_alpha`) |
| `--alpha-practical R...` | `0.99` | approximate only: further stiffnesses, as multiples of the per-cell estimate eq:regularized-practical (`bourne_cases.practical_alpha`) |
| `--start S...` | `0` | `s` of `MakeJointStart` for our H when the sites differ from the targets (2025-invert) |
| `--max-iter-mine N` | 200 | Newton iteration cap of my solves |
| `--max-iter-theirs N` | 5000 | iteration cap of their common-rule runs |
| `--max-seconds T` | 1800 | time cap of their common-rule runs |
| `--out DIR` | `benchmarks/bourne/results/realizable` or `.../approximate` | one JSON per run and `summary.txt` |

## Datasets

Every numerical example of the two papers is a case of `make_cases()` in
`bourne_cases.py`, generated as in their notebooks
(`external/bourne/Laguerre-Polycrystalline-Microstructures/Examples/Inverting_Laguerre_Tessellations_2024/`
for 2025, `external/bourne/Laguerre-fitting/semidiscrete_convex_order.ipynb`
for 2026). Nothing is downloaded or cached: synthetic cases are rebuilt with
their fixed random seeds on every run, and the EBSD case reads the same
file (`EBSD_data.txt`, identical in both repositories; 243 grains; column 1
is the grain id, columns 2–3 the centroid x, y in µm and column 4 the area
in µm²; the rest is unused). The commands below assume

```sh
export MPLCONFIGDIR=/tmp/mplcfg
PY=external/bourne/.venv/bin/python
```

and run every method of the case; add `--method` to narrow them.

### 2025-invert — realizable, N = 20

Bourne 2025 §5.1, Figs. 1–2. Voronoi diagram of 20 random seeds in the unit
square (seed 999); data = its cell areas and centroids. Their SLSQP on H
starts from 20 further random seeds (the next draws of the same stream).

```sh
$PY benchmarks/bourne/run_realizable.py --case 2025-invert
```

### 2025-perturbed-small, 2025-perturbed-large — approximate, N = 20

Bourne 2025 §5.2.1 (Figs. 3–5) and §5.2.2 (Figs. 5–7). The centroids of
2025-invert, each moved by a random vector of length up to 1e-3 (small) or
5e-2 (large) (seed 1002), then made compatible (volumes summing to the box,
mean centroid at its centre). Their SLSQP starts from the targets.

```sh
$PY benchmarks/bourne/run_approximate.py --case 2025-perturbed-small 2025-perturbed-large
```

### 2026-voronoi — realizable, N = 20

Bourne 2026 §6.3.1, Figs. 8–10. Same Voronoi data as 2025-invert. FW with
t = 100 from the true seeds, as in the notebook (the paper's Fig. 8 uses
t = 1250 and a random start).

```sh
$PY benchmarks/bourne/run_realizable.py --case 2026-voronoi
```

### 2026-aligned — realizable, N = 728

Bourne 2026 §6.3.1, Figs. 11–14. Seeds on a 40 × 40 grid of the unit
square, kept in a cross (four 14 × 13 arms, the centre and corners removed),
equal volumes 1/728; data = centroids of the equal-volume power diagram of
these seeds (OT solved to 1e-7·min V̄). FW with t = 100 from the targets
(their initial condition 2). This tessellation is not irreducible (paper
§3.4): the centroids can be matched without recovering the seeds (Fig. 12).

```sh
$PY benchmarks/bourne/run_realizable.py --case 2026-aligned
```

### ebsd — approximate, N = 243

The Tata Steel EBSD map of both papers: grain areas and centroids, with the
same file in both repositories. It is one case: our H_alpha runs once, in
the unit square. Their methods run with each paper's own scaling and
settings, and are reported as `theirs-2025` and `theirs-2026`:
- Bourne 2025 §5.3, Figs. 8–9: least squares (finite-difference gradient),
  from the targets, in a square box of side sqrt(sum of areas) ≈ 251.25 µm;
  it has no common rule.
- Bourne 2026 §6.3.2, Figs. 15–16: in the unit square, FW with t = 10,
  started at 0.5 as in the notebook.

Our problem does not depend on the scale (α is dimensionless), so the two
earlier cases, `2025-ebsd` and `2026-ebsd`, were the same problem for our
solver.

```sh
$PY benchmarks/bourne/run_approximate.py --case ebsd
```

### All cases

```sh
$PY benchmarks/bourne/run_realizable.py     # 2025-invert, 2026-voronoi, 2026-aligned
$PY benchmarks/bourne/run_approximate.py    # 2025-perturbed-*, ebsd
$PY benchmarks/bourne/run_comparison.py     # head-to-head on the approximate cases (README_COMPARISON.md)
```

## Problem groups (independent of the solver)

| group | cases | why |
|---|---|---|
| **realizable**: a power diagram with volumes V̄ and centroids b̄ exists | `2025-invert`, `2026-voronoi`, `2026-aligned` | the data are the volumes and centroids of a known power diagram (Voronoi of random seeds; grid of aligned seeds) |
| **approximate**: no power diagram matches (V̄, b̄); only a fit is possible | `2025-perturbed-small`, `2025-perturbed-large`, `ebsd` | perturbed centroids (Bourne 2025, §5.2.1: the maximum of their H is positive, so no exact diagram) or measured grains |

Each group has its own driver: `run_realizable.py` or `run_approximate.py`.

### Which of our solvers

- **Realizable: H** (`JointSolver`). A power diagram with exactly V̄ and b̄
  exists, so the maximum of H is attained and H converges to it.
- **Approximate: H_alpha** (`SoftenedSolver`), the softened centroid match.
  No power diagram has V̄ and b̄ (perturbed centroids; real grains are not
  polygonal), so H has no maximizer to converge to; H_alpha always has one
  and gives the best fit for its stiffness. This is also what their 2026
  method solves: FW with dilation t is H_alpha with α = t − 1 and sites b,
  so `ebsd` (their t = 10) is compared with ours at α = 9 as well.

α sets the trade-off of the fit. Larger α matches the centroids more
closely but lets the generators move far from the sites and come close to
each other. Below the distinctness bound (eq:regularized-distinct,
`bourne_cases.distinct_alpha`)

    α < min_{i≠j} |q_i − q_j| / ( r (1/V̄_i + 1/V̄_j) )^{1/2},   r = Σ_i V̄_i |c⁰_i − b_i|²

with c⁰ the centroids of the volume solve (max F) at the sites q = b
(`--alpha-bound` multiples < 1), the generators are guaranteed distinct, at
the cost of a looser fit. The per-cell estimate (eq:regularized-practical,
`bourne_cases.practical_alpha`)

    α < min_{i≠j} |q_i − q_j| / ( |c⁰_i − b_i| + |c⁰_j − b_j| )

is larger but not guaranteed (`--alpha-practical` multiples). Both need one
OT solve per case. The default sweep (`--alpha 10 100 1000`,
`--alpha-bound 0.99 10 100`, `--alpha-practical 0.99`, plus t − 1 for 2026)
shows the whole range instead of a single value.

## Termination rules

Each method stops on its own rule, and the other rule is measured at its
solution (latex/results.tex, "Termination of the comparison"):

| case | method | stops on |
|---|---|---|
| 2025-invert, -perturbed-* | SLSQP on their H, pairwise (δ = 1e-3) and ball constraints | rule 1: SLSQP `ftol = 1e-8·|H(X_init)|` (their notebook: 1e-10), maxiter 1000; inner OT `max|V-V̄| < 1e-3·min V̄` |
| ebsd (theirs-2025) | SLSQP least squares, pairwise constraints, finite-difference gradient | rule 1: `ftol = 1e-8·f(X_init)`, maxiter 1e4 |
| 2026-voronoi, 2026-aligned (t = 100), ebsd (theirs-2026: t = 10, start 0.5) | FW | the common rule, or `--max-iter-theirs` / `--max-seconds` |
| all | our H or H_alpha | the common rule, or `--max-iter-mine` |

The FW gap of Bourne 2026 is not used as a rule: their paper sets no
tolerance, and their code's (1e-14 to 1e-16) is at round-off for a gap that
is quadratic in the error. FW therefore has no rule of its own.

Rule 1 on my trajectory (2025 cases only) is `|obj_k − obj_{k−1}| ≤ ftol` at
the final iterate, with the 0.1 % volume condition. The common rule is
measured at the solution of their SLSQP when the case has one
(`rule3AtSolution`).

**Common:** `max|V−V̄|/V̄ ≤ tol` and `max_i |g_u,i| / V̄_i^{3/2} ≤ tol`, the
shift gradient of the problem each method solves: H for SLSQP-on-H and my H;
H_alpha (α = t − 1, sites = b) for FW and my H_alpha. Their methods get
tightened inner OT solves (`≤ 0.1·tol·min V̄`) so the volume part can pass.
Their least-squares method (theirs-2025 on ebsd) has no common rule. The perturbed 2025
cases are infeasible for H, so SLSQP-on-H has no common rule there either.

Each method runs once; every rule stops at the first iterate that passes,
timed up to that iterate. Evaluation of the rules is not charged.

## Known issues (2026-11-12)

- **ebsd at large α:** our H_alpha converges at α ≤ 0.41 (5–6 iterations),
  but fails (damping or iteration cap) at α = 9, 10, 100 and 1000, all above
  our bound (≈ 0.015, eq:regularized-distinct). perturbed-large behaves the
  same above ≈ 0.30. This is the regime where the generators of the
  maximizer can collide, not a solver bug. The earlier "no H_alpha iterate"
  issue of the µm-box 2025-ebsd no longer applies: our H_alpha now runs in
  the unit square.
- **2026-aligned: their notebook's FW never meets its gap** (< 1e-16): it
  ends at its iteration cap (1000, about 36 s). This is one reason the gap is
  not used as a rule; the drivers stop FW on the common rule instead.
- **Long runs.** Their methods under the common rule run up to
  `--max-iter-theirs` (5000) or `--max-seconds` (1800 s) per case. Lower
  both for a quick check, e.g.
  `--tol 1e-1 --max-iter-theirs 20 --max-seconds 120 --max-iter-mine 50`.

## Further parameters worth exposing

- The minimum seed distance δ and ball radius of the 2025 SLSQP.
- The dilation t of the 2026 FW (and hence α = t − 1).
- The regularization η and the Levenberg–Marquardt damping of my solver.
- Random seeds of the cases (theirs are fixed: 999, 1002), to run several
  instances per case.
