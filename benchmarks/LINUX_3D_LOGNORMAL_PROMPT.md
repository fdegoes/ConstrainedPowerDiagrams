# Prompt: 3D stress set with lognormal volumes, on the Linux machine

For the user, before pasting: commit and push this file and
`benchmarks/synthetic/generate.py` (its lognormal volumes) first. The results
come back by `git pull`. Do not run it at the same time as another benchmark.

Paste everything below the line into a Claude Code session on the Linux
machine (the one that ran `benchmarks/LINUX_3D_PROMPT.md`), from the project
root.

---

You are running a benchmark for a research paper on power diagrams with
prescribed volumes and centroids. It runs our solver (the paper's G,
`solver H` in the code) on 48 synthetic 3D cases whose cell volumes follow a
lognormal distribution (coefficient of variation 0.5 and 1), a harder
variant of the 3D benchmark already run on this machine. The timings are
part of the result, so the machine must be idle while it runs. Work from the
project root (the folder with `CMakeLists.txt`). Follow the steps in order.

## Rules

- Do not edit `src/`, `tests/` or any script under `benchmarks/`. If a step
  fails, stop and report the command, its output and the log. Do not work
  around it.
- No `sudo`, and no system packages.
- Never run two benchmark steps at the same time. Check first that no other
  benchmark is running (`pgrep -af 'run_table1|run_synthetic|run_theirs'`).
- Long steps go in the background (`nohup`). Check their logs every 10 to 15
  minutes, not more often.

## 1. Update the repository and build

```sh
git status --short | head      # report any local changes; do not discard them
git pull
git log -1 --oneline           # report the commit
grep -c lognormal benchmarks/synthetic/generate.py   # expect 5 or more
if [ -f build/build.ninja ]; then ninja -C build benchBourne makeSynthetic
else make -C build -j"$(nproc)" benchBourne makeSynthetic; fi
```

If the count is 0, the pull did not get the new script. Stop and tell the
user. If the build fails, stop and report the end of its output.

## 2. Generate the cases (not timed)

```sh
export MPLCONFIGDIR=/tmp/mplcfg
PY=external/bourne/.venv/bin/python
OUT=benchmarks/synthetic/results3d-lognormal-linux-$(date +%Y-%m-%d)
mkdir -p "$OUT"
"$PY" -c "
import sys; sys.path.insert(0, 'benchmarks/synthetic')
import generate as G
for n in (10, 30, 89, 264, 785, 2336, 6952, 10000):
    for vol in ('lognormal50', 'lognormal100'):
        for seed in (1, 2, 3):
            G.make_case('$OUT', n, 1, vol, seed, dim=3)
    print('n =', n, 'done', flush=True)
" > "$OUT/generate.log" 2>&1
tail -2 "$OUT/generate.log"
ls -d "$OUT"/*/ | wc -l        # expect 48
```

This takes a few minutes. If the count is not 48, or the log has a traceback,
stop and report it.

## 3. The run (about 10 to 20 minutes)

Check that the machine is idle first: `uptime` should show a load average
below 1. If it does not, report the load and the top processes to the user
before you start.

```sh
uptime > "$OUT/load-start.txt"
nohup "$PY" benchmarks/synthetic/run_synthetic.py --dim 3 --force \
    --n 10 30 89 264 785 2336 6952 10000 --noise 1 --volumes lognormal50 lognormal100 --seeds 1 2 3 \
    --out "$OUT" > "$OUT/run.log" 2>&1 &
echo $! > "$OUT/run.pid"
```

Each case prints one line (status, iterations, seconds, residuals, mismatch,
load). The run ends with `converged ... / 48; summary in .../ours.csv`. For
reference, the uniform-volume 3D cases took 13 to 44 s each at n = 10000.

```sh
uptime > "$OUT/load-end.txt"
tail -n 20 "$OUT/run.log"
```

## 4. Commit the results and report

Commit the CSV and the logs only (the case files stay on the machine). The
folder is ignored, so add the files with `-f`:

```sh
git add -f "$OUT/ours.csv" "$OUT/run.log" "$OUT/generate.log" "$OUT/load-start.txt" "$OUT/load-end.txt"
git commit -m "3D lognormal-volume stress set of our G on the Linux machine" -- "$OUT"
git push
```

If the push fails, stop and tell the user; leave the commit in place.

Report to the user: the commits, `uptime` at the start and end, the count of
converged cases and the line from `run.log` of any other case, and the median
and maximum of `seconds` and `iterations` per n and volume distribution (a
short table).
