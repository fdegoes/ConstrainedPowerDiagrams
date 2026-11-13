# Overnight comparison: our H vs. Bourne et al.

`run_overnight.sh` runs our solver and the methods of Bourne et al. on all
the realizable data, and writes the tables and timing plots for the paper.
The timings are the point of the run, so run it on an idle machine: nothing
else running, laptop plugged in, sleep disabled.

## How to run

From the project root, in a normal terminal (not the Claude Code sandbox):

```sh
# 1. Smoke test (about 4 minutes, writes only to /tmp/overnight-quick)
OVERNIGHT_QUICK=1 bash benchmarks/run_overnight.sh

# 2. The real run (worst case about 5 hours)
caffeinate -i bash benchmarks/run_overnight.sh
```

- `caffeinate -i` keeps macOS from sleeping during the run.
- The script prints a warning if the load average is already ≥ 2 when it
  starts.
- Re-running on the same day overwrites that day's folder.

Requirements:

- `cmake` and the build tree in `build/` (`sh config.sh` once;
  `PP_BUILD_HOUDINI=OFF sh config.sh` without Houdini).
- The Python environment with pysdot, `external/bourne/.venv`, made by
  `bash benchmarks/setup_bourne.sh`.
- Their code in `external/bourne/` (`Laguerre-fitting`,
  `Laguerre-Polycrystalline-Microstructures`), git submodules (clone with
  `--recursive`).

On a Pixar Linux machine, follow `OVERNIGHT_LINUX_PLAN.md`, which covers the
clone, the dependencies and the setup.

## What it runs

The steps run one after the other, never in parallel, so that the timings
are clean. **All our runs come first**, so our results are complete even
if the night is cut short.

| step | log | what | data |
|---|---|---|---|
| 1 | `build.log` | build `benchBourne` and `makeSynthetic` | — |
| 2 | `ours-bourne.log` | our H | the 3 realizable Bourne cases |
| 3 | `ours-synthetic.log` | our H | all 600 synthetic 2D cases |
| 4 | `theirs-bourne.log` | their SLSQP (2025 cases) and FW (2026 cases) | the 3 realizable Bourne cases |
| 5 | `theirs-synthetic.log` | their SLSQP and FW | synthetic, seed 1, n ≤ 500 (66 cases), 120 s cap per run |
| 6 | `report.log` | tables and plots | — |

A failed step is logged as `FAILED` in `overnight.log`, and the next step
still runs.

The realizable Bourne cases are `2025-invert` and `2026-voronoi` (the same
Voronoi diagram, n = 20) and `2026-aligned` (the cross, n = 728). The
synthetic cases are described in `synthetic/README.md`. All of these data
are realizable, so our runs use H (H_α is for approximate data only).

## How each method is stopped

Each method stops on **its own rule at 1e-8**, and the **other rule is
measured at its solution**:

| method | stops on | measured at its solution |
|---|---|---|
| our H | rule 3 | rule 1 |
| their SLSQP (Bourne 2025) | rule 1 | rule 3 |
| their FW (Bourne 2026) | rule 3, or the time cap | — (it has no rule 1) |

- **Rule 1**, the rule of Bourne 2025: the change of their objective
  between iterations is ≤ 1e-8 times its value at the start, with every
  cell volume within 0.1 % of its target. Their code uses 1e-10; we use
  1e-8.
- **Rule 3**, ours (`latex/results.tex`, eq. `results-termination`): for
  every cell, the relative volume error and the centroid error relative to
  the cell size are ≤ 1e-8. On realizable data it is always the residual of
  H.
- **The FW gap is not used as a rule.** Their paper sets no tolerance, the
  tolerances of their code (1e-14 to 1e-16) are at round-off, and FW
  reaches its iteration cap on the cross. FW with t = 100 solves H_α with
  α = t − 1 = 99, so its residual of H stalls at about 1/t and its runs
  usually end at the cap.
- **Runs that hit the cap** are reported with their rule-3 residuals
  (volume and centroid) at the last iterate, so the report shows how far
  from the solution they stopped.
- Evaluating a rule (an extra OT solve per iterate for their methods) is
  not charged to the method's time.

## Output

Everything goes to `benchmarks/overnight/<date>/`:

| file | contents |
|---|---|
| `overnight.log` | start, end, duration and exit code of every step |
| `environment.txt` | date, git commit, number of modified files, CPU, load at the start and the end |
| `<step>.log` | the output of each step |
| `bourne-realizable/` | our H on the Bourne cases: one JSON per run, `summary.txt` |
| `bourne-realizable-theirs/` | their methods on the Bourne cases: one JSON per run, `summary.txt` |
| `synthetic/ours.csv` | our H, one row per synthetic case (status, iterations, seconds, residuals, rule 1 at the solution) |
| `synthetic/theirs.csv` | their methods, one row per run (rule it stopped on, seconds, rule-3 residuals at the solution) |
| `synthetic/<case>/` | the case data and the per-run JSON files |
| `REPORT_TABLES.md` | the tables: Bourne cases one row each; synthetic, all cases (ours) and the matched cases (ours and theirs); residuals of the runs that hit the cap |
| `timing.pdf`, `.png` | seconds against n for our H, one curve per volume type and noise level (the scaling plot for the paper) |
| `timing-with-theirs.pdf`, `.png` | the same, with their methods where at least half their cases stopped on their own rule |

To regenerate the tables and plots from an existing run:

```sh
external/bourne/.venv/bin/python benchmarks/report_overnight.py --dir benchmarks/overnight/<date>
```

## Expected time

| step | idle machine |
|---|---|
| ours, Bourne | under 1 min |
| ours, synthetic (600 cases) | about 20 min (case generation included) |
| theirs, Bourne | up to 30 min per case (`--max-seconds 1800`) |
| theirs, synthetic | up to 66 cases × 2 methods × 120 s ≤ 4.4 h |

FW usually runs to the cap, so expect close to the worst case for it.
SLSQP mostly stops on rule 1 well before the cap.

## Changing the settings

The settings are at the top of `run_overnight.sh`:

- `BOURNE_ARGS`: `--max-seconds` for their methods on the Bourne cases.
- `THEIRS_ARGS`: `--max-n` (largest synthetic size for their methods) and
  `--max-seconds` (cap per run).
- `OURS_ARGS`: empty, i.e. all sizes and seeds. Use e.g. `--n 10 100` for a
  subset.

The quick mode (`OVERNIGHT_QUICK=1`) uses n = 10 and 14, two Bourne cases,
5 s caps, and writes to `/tmp/overnight-quick`.

## Scripts it calls

- `bourne/run_realizable.py`: the Bourne cases (`--method H` or
  `--method theirs`). The cases and rules are in `bourne/bourne_cases.py`,
  their methods in `bourne/bourne_ref.py`.
- `synthetic/run_synthetic.py`: our H on the synthetic cases.
- `synthetic/run_theirs.py`: their SLSQP and FW on the synthetic cases.
- `synthetic/their_rules.py`: rule 1, measured at our solution.
- `report_overnight.py`: the tables, and `synthetic/plot_timing.py` for the
  plots.
