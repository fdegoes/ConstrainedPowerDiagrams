# Plan: four-vortex experiments (PP, GM, constrained, relaxed) on the Linux machine

Hand this file to a Claude Code session on the Linux machine ("execute
experiments/pp2d/RUN_4VORT.md"). It runs the three four-vortex examples
(uniform, band2, radial; Power Particles inputs made centroidal with BNOT)
with four integrators. The inputs are already in the repository
(`experiments/pp2d/yesweights/`), so BNOT does not need to run again.

Follow the environment steps of `benchmarks/OVERNIGHT_LINUX_PLAN.md` first
(dependencies, `build/runPP2D`, `external/bourne/.venv`).

## 0. Rules for the session

- Timings are not the output here; the results are. Still, run nothing else
  heavy at the same time.
- Work from the project root. Do not edit any script or input.
- Run long steps with `nohup` in the background and check the log.
- If a run fails, do not retry it with different settings: record the message.

## 1. Check the inputs and the binary

```sh
mkdir -p logs
ls experiments/pp2d/yesweights/
ls -la build/runPP2D
```

Expect `uniform_5k_4vort_centroidal.pp`, `band2_5k_4vort_centroidal.pp`,
`radial_5k_4vort_centroidal.pp`, `param.txt`, `square.line`, and the
original `uniform_5k.pp`, `radial_5k.pp`, `band2_5k_4vort.pp`.

If `build/runPP2D` is missing: `cmake --build build --target runPP2D`.

## 2. Run the four integrators on the three cases (100 frames)

```sh
nohup external/bourne/.venv/bin/python experiments/pp2d/run_experiments.py \
  --frames 100 > logs/4vort-all.log 2>&1 &
```

The script runs, for each case: `powerparticles`, `gm` (stiffness 1),
`constrained`, and `softened` with alpha = 1. It starts all of them from the
same centroidal particles, with the start tolerance `--lloydtol 1e-6`. It
checks the time-step bound and the distinctness check before the relaxed
runs and prints both. Outputs: `experiments/pp2d/runs/<case>/<integrator>/`
(frame CSVs, `energy.csv`, `run.txt`). `runs/` is not committed.

Expected from the Mac run (for comparison, not a requirement):

| integrator | uniform | band2 | radial |
|---|---|---|---|
| powerparticles | 100 frames, max volume error 3.9e-3 | 100 frames, 2.5e-2 | 100 frames, 2.2e-2 |
| gm | 100 frames, 6.4e-7 | 100 frames, 4.2e-7 | 100 frames, 2.6e-6 |
| constrained | fails at frame 1 (max iterations) | same | same |
| softened, alpha = 1 | fails at frame 58 (damping failed) | frame 62 | frame 55 |

Send me the log, and report any difference from this table.

## 3. Summaries for the figures

For each run, read `energy.csv` and report:

- the number of frames completed and the status;
- the maximum relative change of the total energy over the run,
  `max |E(t) - E(0)| / E(0)`;
- the maximum change of the momentum (`momentum_x`, `momentum_y`);
- the maximum volume error (`maxVolumeError`) and the centroid gap
  (`maxCentroidGap`) over the frames;
- the total number of substeps and Newton iterations.

Write the summary to `experiments/pp2d/runs/summary.csv` (one row per run)
with a short Python script (`external/bourne/.venv/bin/python`), and send it
to me.

## 4. Report back

- the commit hash (`git log -1 --oneline`);
- the log of step 2 (the status line of every run);
- `summary.csv` from step 3;
- for the failed runs, the last frame written and the error message.
