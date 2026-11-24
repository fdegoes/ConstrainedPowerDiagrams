# Plan: 3D synthetic runs of method G on the Linux machine

Hand this file to a Claude Code session on the Linux machine ("execute
benchmarks/OVERNIGHT_3D_PLAN.md"). It runs method G (the centroid dual,
run_synthetic_g.py) on the 3D synthetic benchmark. The figures are made
later, once their script is in the repository.

Follow the environment steps of benchmarks/OVERNIGHT_LINUX_PLAN.md first
(dependencies, build/benchBourne, external/bourne/.venv).

## 0. Rules for the session

- Timings are the output. Run nothing else heavy while the steps run.
- Work from the project root (the folder with `CMakeLists.txt`).
- Do not edit `src/`, `tests/` or any benchmark script. If a script fails,
  report it with the log instead.
- Long steps go in the background (`nohup`), and their logs are checked
  periodically.
- Record the machine load (`uptime`) at the start and end of each step.

## 1. Check the 3D data

The 600 cases are already in `benchmarks/synthetic/results3d/`, one folder
per case (data.txt, truth.json). Do not regenerate them: the script
creates any missing case by itself.

```sh
ls -d benchmarks/synthetic/results3d/*/ | wc -l    # expect 600
```

Report the count.

## 2. Small test of G in 3D

```sh
mkdir -p logs
external/bourne/.venv/bin/python benchmarks/synthetic/run_synthetic_g.py \
  --dim 3 --n 10 30 --noise 0 --volumes uniform --seeds 1 > logs/g3d-test.log 2>&1
tail -20 logs/g3d-test.log
```

Each case should print `converged` or `invalid start`. If a case prints
`failed`, stop here and send the log tail.

## 3. Full 3D run of G, all 600 cases

```sh
nohup external/bourne/.venv/bin/python benchmarks/synthetic/run_synthetic_g.py \
  --dim 3 > logs/g3d-all.log 2>&1 &
```

Runs are sequential. The script skips finished cases, so it can be
restarted after an interruption.

Outputs:
- `benchmarks/synthetic/results3d/oursG.csv` (one row per case)
- `benchmarks/synthetic/results3d/<case>/oursG.json` (one per case)

When it finishes, report:
- counts of `converged`, `invalid start` and `failed`;
- wall time of the step;
- the machine load at the start and end (`uptime`).

## 4. Report back

- the commit hash of this repository on the machine (`git log -1 --oneline`);
- the result of step 1 (count) and step 2 (test outcome);
- the counts and wall time of step 3;
- the median and maximum time of G per n (from `oursG.csv`);
- the log tail for any case whose status is neither `converged` nor
  `invalid start`.
