# Head-to-head comparison on the approximate cases

`run_comparison.py` compares the methods of Bourne et al. with our max H_α on
the cases for which no power diagram matches the targets, so only a fit is
possible. Every method runs until it terminates, and every result is
measured with the same residuals.

## Cases

| case | n | domain | data |
|---|---|---|---|
| `2025-perturbed-small` | 20 | unit square | Voronoi cells of 20 random generators (seed 999), centroids moved by up to 0.001 (seed 1002) |
| `2025-perturbed-large` | 20 | unit square | same data, centroids moved by up to 0.05 |
| `ebsd` | 243 | unit square (ours, their FW); 251.25 µm square (their 2025 least squares) | the Tata Steel EBSD map used in both papers |

The cases are built in `bourne_cases.py`, exactly as in their notebooks.

## Stiffness α

Each case runs with up to three stiffnesses. FW and ours run at the same
values:

- **bound**: 0.99 × our bound (latex/contributions.tex,
  eq:regularized-distinct; `bourne_cases.distinct_alpha`)

      α < min_{i≠j} |q_i − q_j| / ( r (1/V̄_i + 1/V̄_j) )^{1/2},   r = Σ_i V̄_i |c⁰_i − c̄_i|²

  Here q = c̄ are the sites, and c⁰ are the centroids of the volume solve
  (max F) at q. Below this bound, the maximizer of H_α is a power diagram,
  i.e. its generators are distinct. It is computed for each case: 0.777
  (perturbed-small), 0.299 (perturbed-large) and 0.0152 (ebsd).
- **practical**: 0.99 × our per-cell estimate (`bourne_cases.practical_alpha`)

      α < min_{i≠j} |q_i − q_j| / ( |c⁰_i − c̄_i| + |c⁰_j − c̄_j| )

  It assumes that each shift is at most α times the residual of its cell,
  |u_i| ≤ α |c⁰_i − c̄_i|, which is not proven. In sweeps over the
  synthetic cases (2D and 3D, up to 1129 cells) and these three cases, the
  generators stayed distinct up to at least 1.5 times this value. It is 4.21
  (perturbed-small), 1.21 (perturbed-large) and 0.938 (ebsd), i.e. 4 to 60
  times the bound.
- **theirs**: their hand-picked value. Bourne 2026 (§6.3.2) sets t = 10 for
  the EBSD map, which is α = t − 1 = 9. Their papers give no rule for t: the
  theory only asks for t "large enough", and they observe that "numerically,
  the approach seems robust". The 2025 SLSQP on H has no stiffness, so the
  perturbed cases have no "theirs" value.

Change the factor of both of ours with `--bound-factor`.

Their Frank–Wolfe with dilation t solves the same projection as our H_α with
α = t − 1, so the two methods are compared at equal α.

## Runs

| case | method | α |
|---|---|---|
| perturbed-small, perturbed-large | their SLSQP on H (2025) | none |
| | our H_α | bound, practical |
| ebsd | their least-squares SLSQP (2025, finite-difference gradient, µm box) | none |
| | their FW (2026), started at 0.5 | bound, practical, theirs (t = 1 + α) |
| | our H_α, started at the sites q = c̄ | bound, practical, theirs |

These are 13 runs in total.

## Stopping

Each run stops at the first of:

- **its own convergence rule:**
  - their SLSQP: relative change of their objective below ptol (1e-8, or
    1e-8 of f for least squares);
  - their FW: their gap rule (1e-9 for EBSD);
  - ours: volume and shift residuals ≤ `--tol` (1e-8);
- **the iteration cap**, `--max-iter` (1000 iterations);
- **the time cap**, `--max-seconds` (1800 s of wall time).

For their methods the time cap is checked at every iterate. Our solver is
not given it: benchBourne runs to the iteration cap, and the process is
killed after `--max-seconds` + 60 s, which is reported as `time cap` at
`--max-seconds` seconds, with no residuals.

## What is reported

For every run:

| field | meaning |
|---|---|
| `stop` | `converged`, `iteration cap`, `time cap`, or a failure reason (e.g. `damping failed`) |
| `iterations`, `seconds` | iterations and wall time to the stop (evaluations of the residuals are not timed) |
| `otSolves` / `builds` | OT solves (theirs) or power diagrams built (ours) |
| `volume` | max_i \|V_i − V̄_i\| / V̄_i of the returned diagram |
| `centroid` | max_i \|c_i − c̄_i\| / √\|domain\|, the fitting error |
| `stationary` | the per-cell shift gradient of the problem the run solves, as in the common rule of `bourne_cases.py`: H for SLSQP on H, H_α for FW and ours. Empty for the least-squares SLSQP, which solves neither. Small means the run actually solved its problem. |

All three residuals are evaluated in the same way for all methods, on the
power diagram of the returned generators with the prescribed volumes (a
tight OT solve, not charged to the method).

## How to run

Requirements: build `benchBourne`, and use the Python with pysdot.

```sh
cmake --build build --target benchBourne
```

To see the plan (the runs and the α values) without running anything:

```sh
external/bourne/.venv/bin/python benchmarks/bourne/run_comparison.py --list
```

For a quick check (a few minutes, small caps):

```sh
external/bourne/.venv/bin/python benchmarks/bourne/run_comparison.py \
    --max-iter 30 --max-seconds 30 --out /tmp/comparison_check
```

For the full run (about 1–1.5 h; the EBSD least squares usually runs to its
30-minute cap), run it overnight. Use a normal terminal: the Claude Code
sandbox is not needed. `caffeinate` keeps the Mac awake and `nohup` keeps the
run alive if the terminal closes:

```sh
cd ~/Documents/Projects/PowerElasticity
mkdir -p benchmarks/bourne/results
caffeinate -i nohup external/bourne/.venv/bin/python benchmarks/bourne/run_comparison.py \
    > benchmarks/bourne/results/comparison.log 2>&1 &
```

To follow its progress:

```sh
tail -f benchmarks/bourne/results/comparison.log
```

### Options

| option | default | meaning |
|---|---|---|
| `--case NAME...` | all three | `2025-perturbed-small`, `2025-perturbed-large`, `ebsd` |
| `--max-iter N` | 1000 | iteration cap of every method |
| `--max-seconds T` | 1800 | wall-time cap of every run |
| `--tol T` | 1e-8 | our volume and shift tolerance |
| `--bound-factor F` | 0.99 | α = F × our bound, and F × our practical estimate |
| `--out DIR` | `benchmarks/bourne/results/comparison` | outputs |
| `--list` | | print the plan and exit |

## Outputs

In `--out`:

- `results.csv`: one row per run, with the fields above. It is rewritten
  after each run, so a partial run still leaves its finished rows.
- `summary.txt`: the same table, aligned for reading.
- `<case>__<method>__<alpha>.json`: the full record of each run.

## What to expect

From a run with small caps (30 iterations, 30 s):

- **At our bound:** our H_α converges in 3 Newton iterations, with
  stationarity residuals of 1e-10. Their FW on EBSD stops on its gap after 2
  iterations, but at a stationarity residual of 3e-4, so their gap rule fires
  before the projection is solved.
- **At our practical α:** our H_α converges in 4–6 iterations (stationarity
  5e-9 to 5e-12), with a smaller centroid error (EBSD: 1.6e-2 against
  1.9e-2). Their FW stops on its gap after 13 iterations, at a stationarity
  residual of 3e-2.
- **At their α = 9 (EBSD):** both FW and ours reach the cap. α = 9 is about
  600 times our bound and 10 times our practical estimate, where the
  generators of the maximizer can collide.
- **Their SLSQP on perturbed-large:** still far from stationary after 30
  iterations (1.4). The full caps show whether it gets there.
