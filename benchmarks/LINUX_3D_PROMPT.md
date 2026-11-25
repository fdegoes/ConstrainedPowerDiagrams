# Prompt: 3D synthetic runs of our G on the Linux machine

For the user, before pasting: commit and push this file first, so that
`git pull` on the Linux machine gets it. The results come back by `git pull`. Run it after the Table 1 run
(`benchmarks/LINUX_TABLE1_PROMPT.md`) has finished, never at the same time.

Paste everything below the line into a Claude Code session on the Linux
machine (the one that ran `benchmarks/OVERNIGHT_LINUX_PLAN.md`), from the
project root.

---

You are running a benchmark for a research paper on power diagrams with
prescribed volumes and centroids. The paper calls our objective G (the joint
dual over shifts and weights); the code and the scripts call it H
(`solver H` in the case files). Only our method runs here: the codes of
Bourne et al. are 2D only. The timings are the point of the run, so the
machine must be idle while it runs. Work from the project root (the folder
with `CMakeLists.txt`). Follow the steps in order.

## Rules

- Do not edit `src/`, `tests/` or any script under `benchmarks/`. If a step
  fails, stop and report the command, its output and the log. Do not work
  around it.
- No `sudo`, and no system packages. If a download or network access is
  needed and fails, stop and report it.
- Never run two benchmark steps at the same time, and run nothing else heavy
  during the full run (step 5). Check first that the Table 1 run is not
  still going (`pgrep -af run_table1`).
- Long steps go in the background (`nohup`). Check their logs every 20 to 30
  minutes, not more often.
- Do not use `benchmarks/OVERNIGHT_3D_PLAN.md` if an old copy exists: it ran
  a different method (the fixed-weights dual) and is retired.

## 1. Update the repository

```sh
git status --short | head      # report any local changes; do not discard them
git pull
git submodule update --init --recursive
git log -1 --oneline           # report the commit
```

## 2. Build

The build tree comes from `benchmarks/OVERNIGHT_LINUX_PLAN.md` (step 5).
If `build/` does not exist, follow that step first, then come back here.

```sh
if [ -f build/build.ninja ]; then ninja -C build benchBourne makeSynthetic
else make -C build -j"$(nproc)" benchBourne makeSynthetic; fi
grep '^PP_SPARSE_SOLVER:' build/CMakeCache.txt      # expect Cholmod
ldd build/benchBourne | grep -E 'cholmod|not found' # a path, no "not found"
```

If the build fails, stop and report the end of its output.

## 3. Generate the 3D cases (not timed)

The 600 cases (20 sizes n from 10 to 10000, three noise levels, uniform and
Voronoi volumes, five seeds) are not in git. Generate them into the output
directory first, so that the timed run does not include their generation.
The generation is deterministic (seeded) and runs `build/makeSynthetic`.

```sh
export MPLCONFIGDIR=/tmp/mplcfg
PY=external/bourne/.venv/bin/python
OUT=benchmarks/synthetic/results3d-linux-$(date +%Y-%m-%d)
mkdir -p "$OUT"
df -h .                        # the run needs about 3 GB free; report it if less
nohup "$PY" -c "
import sys; sys.path.insert(0, 'benchmarks/synthetic')
import generate as G, run_synthetic as RS
for n in RS.NS:
    for noise in G.NOISE:
        for vol in G.VOLUMES:
            for seed in (1, 2, 3, 4, 5):
                G.make_case('$OUT', n, noise, vol, seed, dim=3)
    print('n =', n, 'done', flush=True)
" > "$OUT/generate.log" 2>&1 &
```

Wait for `n = 10000 done` in `$OUT/generate.log` (minutes to an hour), then:

```sh
ls -d "$OUT"/*/ | wc -l        # expect 600
```

If the count is not 600, or the log has a traceback, stop and report it.

## 4. Smoke test (about 1 minute)

```sh
"$PY" benchmarks/synthetic/run_synthetic.py --dim 3 --n 10 30 127 --noise 0 2 \
    --volumes uniform voronoi --seeds 1 --out /tmp/g3d-quick
```

Each of the 12 lines should show `converged` (in about 5 to 15 iterations),
with residuals below 1e-8 and a mismatch near 0. Report the output. If a case
shows `benchBourne failed` or a status other than `converged`, stop and report
it.

## 5. The full run (about 1 to 3 hours)

Check that the machine is idle first: `uptime` should show a load average
below 1. If it does not, report the load and the top processes to the user
before you start.

```sh
uptime > "$OUT/load-start.txt"
if command -v systemd-inhibit >/dev/null; then
    nohup systemd-inhibit --what=sleep:idle --why="3D synthetic benchmark" \
        "$PY" benchmarks/synthetic/run_synthetic.py --dim 3 --force --out "$OUT" > "$OUT/run.log" 2>&1 &
else
    nohup "$PY" benchmarks/synthetic/run_synthetic.py --dim 3 --force --out "$OUT" > "$OUT/run.log" 2>&1 &
fi
echo $! > "$OUT/run.pid"
```

Runs are sequential, by increasing n. Each case prints one line (status,
iterations, solver seconds, residuals, mismatch, load). Check progress with
`tail -n 3 "$OUT/run.log"`. The CSV `$OUT/ours.csv` is rewritten after every
case. If the run is interrupted, restart the same command without `--force`
(finished cases are then read back, not rerun); tell the user that you did.

For reference, a run on a laptop (not idle) took 0.7 s per case at n = 1129
and 13 to 21 s at n = 10000 (median 12.8 s). The solver time per case is the
`seconds` column. The wall time per case is longer, because the script also
measures the diagram mismatch on 10^6 random points. Report it if the run
goes past 5 hours, or if the load column rises above 2.

The run ends with the line `converged ... / 600; summary in .../ours.csv`.

```sh
uptime > "$OUT/load-end.txt"
tail -n 20 "$OUT/run.log"
```

## 6. Report and copy back

Report to the user:

- the commit (step 1), and `uptime` at the start and end of the run;
- the count of converged cases, and for every other case its line from
  `run.log`;
- the median and maximum of `seconds` and `iterations` per n, from
  `$OUT/ours.csv` (a short table);
- the largest `load` value in `ours.csv`;
- the hardware: CPU model (`lscpu | grep 'Model name'`), `nproc`, and RAM
  (`free -g`).

Commit the CSV and the logs (a few hundred kB; the case files and the per-case
`ours.json` files stay on the machine). The folder is ignored, so add the files
with `-f`, and commit only them:

```sh
git add -f "$OUT/ours.csv" "$OUT/run.log" "$OUT/generate.log" "$OUT/load-start.txt" "$OUT/load-end.txt"
git commit -m "3D synthetic benchmark of our G on the Linux machine" -- "$OUT"
git push
```

If the push fails (no credentials or no network), stop and tell the user;
leave the commit in place. Keep `$OUT` on the machine.

---

## On the Mac, after `git pull`

The figure and the statistics for the paper are made from
`benchmarks/synthetic/results3d-linux-<date>/ours.csv`.
