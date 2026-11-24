# Prompt: Table 1 and Fig. 6 with their own stopping rules, on the Linux machine

For the user, before pasting: commit and push this file and
`benchmarks/bourne/run_table1.py` (its `--own-rules` option) first, so that
`git pull` on the Linux machine gets them. This run replaces the Table 1 run of
`benchmarks/LINUX_TABLE1_PROMPT.md`; do not run that one. Do not run it at the
same time as another benchmark (e.g. `benchmarks/LINUX_3D_PROMPT.md`).

Paste everything below the line into a Claude Code session on the Linux
machine (the one that ran `benchmarks/OVERNIGHT_LINUX_PLAN.md`), from the
project root.

---

You are running a benchmark for a research paper on power diagrams with
prescribed volumes and centroids. It runs four methods on two cases: our
Newton solver (the paper's G, `solver H` in the code) and three methods of
Bourne et al. (SLSQP of 2025; Frank–Wolfe and a subgradient method of 2026).
Their methods run with the settings and stopping rules of their papers and
code (`--own-rules`). The timings are the point of the run, so the machine
must be idle while it runs. Work from the project root (the folder with
`CMakeLists.txt`). Follow the steps in order.

## Rules

- Do not edit `src/`, `tests/` or any script under `benchmarks/`. If a step
  fails, stop and report the command, its output and the log. Do not work
  around it.
- No `sudo`, and no system packages. If a download or network access is
  needed and fails, stop and report it.
- Never run two benchmark steps at the same time, and run nothing else heavy
  during the full run (step 5). Check first that no other benchmark is
  running (`pgrep -af 'run_table1|run_synthetic|run_theirs'`).
- Long steps go in the background (`nohup`). Check their logs every 10 to 15
  minutes, not more often.

## 1. Update the repository

```sh
git status --short | head      # report any local changes; do not discard them
git pull
git submodule update --init --recursive
git log -1 --oneline           # report the commit
grep -c -- '--own-rules' benchmarks/bourne/run_table1.py   # expect 2 or more
```

If the count is 0, the pull did not get the new script. Stop and tell the user.

## 2. Build benchBourne

The build tree comes from `benchmarks/OVERNIGHT_LINUX_PLAN.md` (step 5).

```sh
if [ -f build/build.ninja ]; then ninja -C build benchBourne
else make -C build -j"$(nproc)" benchBourne; fi
grep '^PP_SPARSE_SOLVER:' build/CMakeCache.txt      # expect Cholmod
ldd build/benchBourne | grep -E 'cholmod|not found' # a path, no "not found"
```

If the build fails, stop and report the end of its output.

## 3. Check the Python environment

```sh
export MPLCONFIGDIR=/tmp/mplcfg
PY=external/bourne/.venv/bin/python
$PY -c "import numpy, scipy, matplotlib, pysdot; print(numpy.__version__, scipy.__version__, matplotlib.__version__)"
```

If it fails, run `benchmarks/setup_bourne.sh` as in
`OVERNIGHT_LINUX_PLAN.md` step 4. Report the versions.

## 4. Smoke test (about 2 minutes)

```sh
$PY benchmarks/bourne/run_table1.py --own-rules --max-seconds 20 --grace 15 --out /tmp/table1-own-quick
cat /tmp/table1-own-quick/table.tex
```

Expect (the times mean nothing here):

- `ours` on both cases: `rule 3`, residuals below 1e-8.
- On voronoi: `slsqp` and `fw` stop by `own rule` within a few seconds;
  `subgrad` by `iteration cap` (600 iterations).
- On cross: `slsqp`, `fw` and `subgrad` stop by `time cap` at 20 s. One SLSQP
  iteration on the cross takes 20 to 60 s and about 5 GB of RAM, so it may be
  stopped inside its first iteration ("killed inside an iteration"). This is
  expected.
- No `error` in the output or in `table.tex`.

If a run reports `error`, read `/tmp/table1-own-quick/<case>__<method>.log`,
then stop and report it.

## 5. The full run (about 35 to 45 minutes)

Check that the machine is idle first: `uptime` should show a load average
below 1. If it does not, report the load and the top processes to the user
before you start.

```sh
OUT=benchmarks/bourne/results/table1-own-$(date +%Y-%m-%d)
mkdir -p "$OUT"
uptime > "$OUT/load-start.txt"
if command -v systemd-inhibit >/dev/null; then
    nohup systemd-inhibit --what=sleep:idle --why="table 1 benchmark" \
        "$PY" benchmarks/bourne/run_table1.py --own-rules --out "$OUT" > "$OUT/run.log" 2>&1 &
else
    nohup "$PY" benchmarks/bourne/run_table1.py --own-rules --out "$OUT" > "$OUT/run.log" 2>&1 &
fi
echo $! > "$OUT/run.pid"
```

The time cap is 1800 s per run. Runs are sequential: voronoi then cross, and
ours, slsqp, fw, subgrad within each. A line `== <case> <method>` starts a
run and a result line ends it. Check progress with `tail -n 3 "$OUT/run.log"`.

Expected on an idle machine:

| run | stops by | wall time |
|---|---|---|
| ours, both cases | rule 3 | under a second |
| voronoi slsqp, fw | own rule | seconds |
| voronoi subgrad | iteration cap (600) | seconds |
| cross fw | iteration cap (1000) | about 1 minute |
| cross subgrad | iteration cap (600) | about 30 s |
| cross slsqp | time cap (1800 s), or its own rule | up to about 32 minutes |

Report it if the run goes past 1 hour. When it finishes (the log ends with
`wrote .../results.csv and table.tex`):

```sh
uptime > "$OUT/load-end.txt"
cat "$OUT/table.tex"; cat "$OUT/results.csv"
du -sh "$OUT"
```

## 6. Commit the results and report

The output directory is a few MB (JSON histories, logs, the table). It is in
an ignored folder, so add it with `-f`. Commit only this directory:

```sh
git add -f "$OUT"
git commit -m "Table 1 with their own stopping rules (run_table1.py --own-rules)" -- "$OUT"
git push
```

If the push fails (no credentials or no network), stop and tell the user;
leave the commit in place.

Report to the user:

- the commit of step 1 and the commit of the results;
- the versions (step 3), and `uptime` at the start and end of the run;
- per run: stop reason, iterations, seconds and the two residuals (from
  `results.csv`);
- any run whose stop is `error`, with the end of its `.log`;
- the hardware: CPU model (`lscpu | grep 'Model name'`), `nproc`, and RAM
  (`free -g`).

---

## On the Mac, after `git pull`

```sh
export MPLCONFIGDIR=/tmp/mplcfg
PY=external/bourne/.venv/bin/python
DIR=benchmarks/bourne/results/table1-own-<date>
$PY benchmarks/bourne/run_table1.py --report --out $DIR      # rewrites $DIR/results.csv and $DIR/table.tex
$PY benchmarks/figures/plot_convergence.py --from $DIR --out latex/figures/convergence.pdf
```

`table.tex` has the rows of `tab:realizable` (`latex/results.tex`), and the
figure is `fig:convergence`.
