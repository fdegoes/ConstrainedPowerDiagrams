> **Superseded** by `benchmarks/LINUX_TABLE1_OWN_PROMPT.md` (their methods on their own stopping rules). Do not run this one.

# Prompt: Table 1 runs (realizable comparison) on the Linux machine

For the user, before pasting: `benchmarks/bourne/run_table1.py` and the
`--from` option of `benchmarks/figures/plot_convergence.py` must be committed
and pushed first, so that `git pull` on the Linux machine gets them.

Paste everything below the line into a Claude Code session on the Linux
machine (the one that ran `benchmarks/OVERNIGHT_LINUX_PLAN.md`), from the
project root.

---

You are running a benchmark for a research paper on power diagrams with
prescribed volumes and centroids. The paper calls our objective G; the code
and the scripts call it H. The timings are the point of the run, so the
machine must be idle while it runs. Work from the project root (the folder
with `CMakeLists.txt`). Follow the steps in order.

## Rules

- Do not edit `src/`, `tests/` or any script under `benchmarks/`. If a step
  fails, stop and report the command, its output and the log. Do not work
  around it.
- No `sudo`, and no system packages. If a download or network access is
  needed and fails, stop and report it.
- Never run two benchmark steps at the same time, and run nothing else heavy
  during the full run (step 5).
- Long steps go in the background (`nohup`). Check their logs every 20 to 30
  minutes, not more often.

## 1. Update the repository

```sh
git status --short | head      # report any local changes; do not discard them
git pull
git submodule update --init --recursive
git log -1 --oneline           # report the commit
ls benchmarks/bourne/run_table1.py benchmarks/figures/plot_convergence.py
```

If `run_table1.py` is missing, the pull did not get it. Stop and tell the
user.

## 2. Build benchBourne

The build tree comes from `benchmarks/OVERNIGHT_LINUX_PLAN.md` (step 5). If
`build/` does not exist, follow that step first (it uses `config.sh` and the
dependencies of its step 3), then come back here.

```sh
if [ -f build/build.ninja ]; then ninja -C build benchBourne
else make -C build -j"$(nproc)" benchBourne; fi
grep '^PP_SPARSE_SOLVER:' build/CMakeCache.txt      # expect Cholmod
ldd build/benchBourne | grep -E 'cholmod|not found' # a path, no "not found"
```

If the build fails, stop and report the end of its output.

## 3. Check the Python environment

The environment is `external/bourne/.venv`, made by
`benchmarks/setup_bourne.sh` (`OVERNIGHT_LINUX_PLAN.md`, step 4).

```sh
export MPLCONFIGDIR=/tmp/mplcfg
PY=external/bourne/.venv/bin/python
$PY -c "import numpy, scipy, matplotlib, pysdot; print(numpy.__version__, scipy.__version__, matplotlib.__version__)"
$PY benchmarks/bourne/run_realizable.py --list
# expect: 2025-invert (n=20), 2026-voronoi (n=20), 2026-aligned (n=728)
```

If the environment is missing or broken, run
`benchmarks/setup_bourne.sh` as in `OVERNIGHT_LINUX_PLAN.md` step 4. Report
the versions.

## 4. Smoke test (about 1 minute)

```sh
$PY benchmarks/bourne/run_table1.py --quick
cat /tmp/table1-quick/table.tex
```

This runs all 8 runs (2 benchmarks x 4 methods) with 3 s caps and writes
to `/tmp/table1-quick`. The timings mean nothing here. Expect:

- `ours` on both benchmarks: `rule 3`, residuals below 1e-8.
- `slsqp` on voronoi: `rule 3` or `time cap`.
- `slsqp` on cross: `time cap`, 0 iterations, with "killed inside an
  iteration". One SLSQP iteration on the cross takes 30 to 60 s, so the run
  is stopped at the cap plus a grace of 15 s. This is expected.
- `fw` and `subgrad` on both: `time cap`.
- No `error` or `failure` in the output or in `table.tex`.

If a run reports `error`, read `/tmp/table1-quick/<benchmark>__<method>.log`,
then stop and report it.

Also check that the figure renders:

```sh
$PY benchmarks/figures/plot_convergence.py --from /tmp/table1-quick --out /tmp/table1-quick/convergence.pdf
```

The figure is made on the user's Mac; here it only checks that the JSON files load.

## 5. The full run (about 4 hours)

Check that the machine is idle first: `uptime` should show a load average
below 1. If it does not, report the load and the top processes to the user
before you start.

```sh
OUT=benchmarks/bourne/results/table1-$(date +%Y-%m-%d)
mkdir -p "$OUT"
uptime > "$OUT/load-start.txt"
if command -v systemd-inhibit >/dev/null; then
    nohup systemd-inhibit --what=sleep:idle --why="table 1 benchmark" \
        "$PY" benchmarks/bourne/run_table1.py --out "$OUT" > "$OUT/run.log" 2>&1 &
else
    nohup "$PY" benchmarks/bourne/run_table1.py --out "$OUT" > "$OUT/run.log" 2>&1 &
fi
echo $! > "$OUT/run.pid"
```

Defaults: a time cap of 1800 s per run, which does not count the
evaluations of rule 3. Runs are sequential, in the order voronoi then cross,
and ours, slsqp, fw, subgrad within each. Check progress with
`tail -n 3 "$OUT/run.log"`. A line `== <benchmark> <method>` starts a run,
and a result line ends it.

Expected durations on an idle machine. They were measured on a laptop and
include the uncharged evaluations of rule 3 at every iterate.

| run | stops by | wall time |
|---|---|---|
| voronoi ours, cross ours | rule 3 | seconds |
| voronoi slsqp | rule 3 (about 290 iterations) | seconds |
| voronoi fw | time cap | about 50 min (about 0.6 M iterations) |
| voronoi subgrad | time cap | about 40 min |
| cross slsqp | time cap (30 to 60 s per iteration, about 5 GB of RAM) | about 31 min |
| cross fw | time cap | about 50 min |
| cross subgrad | time cap | about 50 min |

The total is about 3.5 to 4 hours. Report it if it runs past 5 hours.

When it finishes (the log ends with `wrote .../results.csv and table.tex`):

```sh
uptime > "$OUT/load-end.txt"
cat "$OUT/table.tex"; cat "$OUT/results.csv"
```

## 6. Report and copy back

Report to the user:

- the commit (step 1), the versions (step 3), and `uptime` at the start and
  end of the run;
- per run: stop reason, iterations, seconds and the two residuals (from
  `results.csv`), and the `ownRule` columns (where SLSQP's rule 1, the FW
  gap and the subgradient's change would have stopped);
- any run whose stop is `error` or `failure`, with the end of its `.log`;
- the hardware: CPU model (`lscpu | grep 'Model name'`), `nproc`, and RAM
  (`free -g`).

Package the output directory for the user's Mac. It is a few MB to a few
tens of MB. Report the path and the size of the tar file.

```sh
tar czf "table1-$(date +%Y-%m-%d)-linux.tar.gz" -C benchmarks/bourne/results "$(basename "$OUT")"
ls -la table1-*-linux.tar.gz
```

---

## On the Mac, after copying the tar file back

Unpack it into `benchmarks/bourne/results/`, then regenerate the table and
the figure from it:

```sh
export MPLCONFIGDIR=/tmp/mplcfg
PY=external/bourne/.venv/bin/python
DIR=benchmarks/bourne/results/table1-<date>
$PY benchmarks/bourne/run_table1.py --report --out $DIR      # rewrites $DIR/results.csv and $DIR/table.tex
$PY benchmarks/figures/plot_convergence.py --from $DIR --out latex/figures/convergence.pdf
```

`table.tex` has the rows of `tab:realizable` (`latex/results.tex`). The figure
(`fig:convergence`) shows both benchmarks and all four methods against the
charged time, each of their methods cut at the iterate where its own rule
stops it (SLSQP rule 1, FW gap, change of the running average), or at the
time cap if that comes first.
