# Overnight run: 2026-11-20

```
start: Fri Nov 20 11:50:40 AM PST 2026
git: ff47a20 (3 modified files)
host: inhuman (Linux x86_64)
cpu: AMD EPYC 7763 64-Core Processor (15 logical cores)
sparse solver: Cholmod
load: 0.37
```

Each method stops on its own rule at 1e-8; the other rule is measured at its solution. Rule 1 = SLSQP objective change (Bourne 2025, ptol = 1e-8); rule 3 = common per-cell residual of H (the data are realizable). FW has no rule of its own (the gap is not a rule), so it stops on rule 3; it solves H_alpha with alpha = t - 1 = 99, so its residual of H stalls at about 1/t.

## Bourne et al., realizable cases

| case | method | stopped by | iterations | seconds | other rule at the solution |
|---|---|---|---|---|---|
| 2025-invert | SLSQP | rule 1 (met) | 190 | 4.51 | rule 3: 1.7e-14 / 3.5e-03 (not met) |
| 2026-aligned | FW (t = 100) | cap (max iterations) | 5000 | 290 | rule 3 residuals at the cap: 8.4e-09 / 4.7e-01 |
| 2026-voronoi | FW (t = 100) | cap (max iterations) | 5000 | 51.7 | rule 3 residuals at the cap: 5.2e-13 / 3.5e-03 |
| 2025-invert | our H | rule 3 | 7 | 0.0125 | rule 1: |Δobj|/|obj₀| 1.4e-05 (not met) |
| 2026-aligned | our H | rule 3 | 27 | 0.243 | – |
| 2026-voronoi | our H | rule 3 | 7 | 0.0108 | – |

## Synthetic 2D

### Our H on all cases (stops on rule 3 at 1e-8)

| converged | median iterations | median seconds (n ≤ 100 / ≤ 1000 / all) | rule 1 met at the solution |
|---|---|---|---|
| 600/600 | 10 | 0.0126 / 0.0181 / 0.0449 | 93/600 |

### Matched cases (seed 1, n ≤ 500)

| method | stops on | stopped by its rule | median seconds | max seconds | other rule met at the solution |
|---|---|---|---|---|---|
| our H | rule 3 | 66/66 | 0.015 | 0.079 | rule 1: 2/66 |
| SLSQP | rule 1 | 57/66 | 9.76 | 102 | rule 3: 0/66 |
| FW (t = 100) | rule 3 | 0/66 | – | – | – (no rule 1) |

#### Runs of theirs that hit the cap: rule 3 residuals at their last iterate

| method | cases | volume residual (median / max) | centroid residual (median / max) |
|---|---|---|---|
| FW (t = 100) | 66 | 1.1e-12 / 5.2e-07 | 3.8e-03 / 3.9e-01 |
| SLSQP | 9 | 3.3e-13 / 1.9e-10 | 2.7e-03 / 2.1e-02 |

Per case: `synthetic/theirs.csv`, columns `finalVolumeResidual` and `finalShiftResidual`.

How their runs ended: fw time cap: 66, slsqp rule 1: 57, slsqp time cap: 9

