#!/bin/bash
# Overnight comparison of our solvers with the methods of Bourne et al., on
# an idle machine (timings): the realizable Bourne cases and the 2D
# synthetic benchmark, under rule 1 (SLSQP ftol, Bourne 2025) and rule 3
# (common per-cell residual). The FW gap is not used as a rule.
#
#   bash benchmarks/run_overnight.sh          # full run, worst case ~8 h
#   OVERNIGHT_QUICK=1 bash benchmarks/run_overnight.sh   # smoke test, minutes
#
# Each method stops on its own rule at 1e-8, and the other rule is measured
# at its solution: our H stops on rule 3; their SLSQP on rule 1; their FW,
# which has no rule of its own (the gap is not a rule), on rule 3. All the
# data are realizable: our runs use H, and rule 3 is the residual of H.
#
# Steps, sequential (never in parallel, so that the timings are clean), ours
# first so that our results are complete even if the run is cut short:
#   1. build benchBourne and makeSynthetic
#   2. ours: H on the realizable Bourne cases
#   3. ours: H on all 600 synthetic 2D cases
#   4. theirs: SLSQP and FW on the realizable Bourne cases
#   5. theirs: SLSQP and FW on the synthetic cases, seed 1, n <= 500,
#      120 s per run
#   6. the approximate examples of Bourne et al. (perturbed small and large,
#      EBSD map), ours and theirs, a separate set of results
#      (bourne/run_comparison.py, its own results.csv and summary.txt)
#   7. report tables and timing plot (report_overnight.py)
# Output: benchmarks/overnight/<date>/ (quick: /tmp/overnight-quick/), one
# log per step. A failed step is logged and the next one runs.
set -u
cd "$(dirname "$0")/.."
PY=external/bourne/.venv/bin/python
export MPLCONFIGDIR=/tmp/mplcfg
if [ ! -x "$PY" ]; then
    echo "No $PY: run bash benchmarks/setup_bourne.sh first" >&2
    exit 1
fi
if [ ! -f build/CMakeCache.txt ]; then
    echo "No build tree: run PP_BUILD_HOUDINI=OFF sh config.sh first" >&2
    exit 1
fi

if [ "${OVERNIGHT_QUICK:-0}" = 1 ]; then
    OUT=/tmp/overnight-quick
    BOURNE_ARGS=(--case 2025-invert 2026-voronoi --max-seconds 5 --max-iter-theirs 50)
    OURS_ARGS=(--n 10 14)
    THEIRS_ARGS=(--max-n 14 --max-seconds 5)
    APPROX_ARGS=(--case 2025-perturbed-small --max-seconds 5 --max-iter 50)
else
    OUT=benchmarks/overnight/$(date +%Y-%m-%d)
    BOURNE_ARGS=(--max-seconds 1800)
    OURS_ARGS=()
    THEIRS_ARGS=(--max-n 500 --max-seconds 120)
    APPROX_ARGS=(--max-seconds 1800)
fi
mkdir -p "$OUT"
LOG="$OUT/overnight.log"

note() { echo "[$(date '+%F %T')] $*" | tee -a "$LOG"; }

# macOS and Linux.
cpu_name() {
    sysctl -n machdep.cpu.brand_string 2> /dev/null ||
        sed -n 's/^model name[[:space:]]*: //p' /proc/cpuinfo 2> /dev/null | head -1
}
cpu_count() { getconf _NPROCESSORS_ONLN; }
load_average() { "$PY" -c "import os; print('%.2f' % os.getloadavg()[0])"; }

environment() {
    {
        echo "$1: $(date)"
        echo "git: $(git rev-parse --short HEAD 2>/dev/null) ($(git status --porcelain 2>/dev/null | wc -l | tr -d ' ') modified files)"
        echo "host: $(hostname) ($(uname -sm))"
        echo "cpu: $(cpu_name) ($(cpu_count) logical cores)"
        echo "sparse solver: $(grep '^PP_SPARSE_SOLVER:' build/CMakeCache.txt 2>/dev/null | cut -d= -f2)"
        echo "load: $(load_average)"
    } >> "$OUT/environment.txt"
}

# step <name> <command...>: runs the command, logs to <name>.log, and
# records the time and the exit code, which it returns (the caller decides
# whether to continue).
step() {
    local name=$1
    shift
    note "start $name: $*"
    local t0=$SECONDS
    "$@" > "$OUT/$name.log" 2>&1
    local code=$?
    note "end   $name: exit $code, $((SECONDS - t0)) s"
    if [ $code -ne 0 ]; then
        note "FAILED $name (see $OUT/$name.log)"
    fi
    return $code
}

# Idle check: the timings need a quiet machine. Refuses to start when the 1-
# or 5-minute load average is >= MAX_LOAD, or when another user is logged in.
# Skipped by the smoke test, whose timings are not kept.
if [ "${OVERNIGHT_QUICK:-0}" != 1 ]; then
    MAX_LOAD=${MAX_LOAD:-1.0}
    read -r LOAD1 LOAD5 < <("$PY" -c "import os; print('%.2f %.2f' % os.getloadavg()[:2])")
    OTHERS=$(who | awk '{print $1}' | sort -u | grep -v -x "${USER:-$(id -un)}" || true)
    if awk -v a="$LOAD1" -v b="$LOAD5" -v m="$MAX_LOAD" 'BEGIN { exit !(a >= m || b >= m) }'; then
        echo "refusing to start: load average $LOAD1 (1 min), $LOAD5 (5 min); need both < $MAX_LOAD (set MAX_LOAD to override)" >&2
        exit 1
    fi
    if [ -n "$OTHERS" ]; then
        echo "refusing to start: other users logged in: $(echo $OTHERS)" >&2
        exit 1
    fi
fi

# Keep the machine awake for the whole run, through systemd-inhibit (not on
# macOS, where caffeinate is used instead). The re-run skips this step.
if [ "${OVERNIGHT_QUICK:-0}" != 1 ] && [ "${OVERNIGHT_INHIBITED:-0}" != 1 ] &&
    command -v systemd-inhibit > /dev/null && systemd-inhibit --what=sleep:idle true > /dev/null 2>&1; then
    export OVERNIGHT_INHIBITED=1
    exec systemd-inhibit --what=sleep:idle --why="overnight benchmark" bash "$PWD/benchmarks/$(basename "$0")" "$@"
fi

environment start
note "output in $OUT; worst case ~8 h for the full run (theirs on synthetic <= 66 cases x 2 x 120 s = 4.4 h, theirs on the approximate Bourne examples <= 6 runs x 1800 s = 3 h)"
note "load at the start: ${LOAD1:-n/a} (1 min), ${LOAD5:-n/a} (5 min); sleep and idle inhibited: ${OVERNIGHT_INHIBITED:-0}"

# Only these targets: the build tree may also configure the Houdini plugin,
# which is not needed here (configure with PP_BUILD_HOUDINI=OFF sh config.sh
# on a machine without Houdini). Stale binaries would give wrong timings, so
# a failed build stops the run.
# The cmake that configured the tree (the cmake on PATH may be too old for --build).
CMAKE_BIN=$(sed -n 's/^CMAKE_COMMAND:INTERNAL=//p' build/CMakeCache.txt)
step build "$CMAKE_BIN" --build build --target benchBourne makeSynthetic || { note "stopping: the build failed"; exit 1; }
step ours-bourne "$PY" benchmarks/bourne/run_realizable.py --method H --tol 1e-8 "${BOURNE_ARGS[@]}" \
    --out "$OUT/bourne-realizable"
step ours-synthetic "$PY" benchmarks/synthetic/run_synthetic.py ${OURS_ARGS[@]+"${OURS_ARGS[@]}"} --force --out "$OUT/synthetic"
step theirs-bourne "$PY" benchmarks/bourne/run_realizable.py --method theirs --tol 1e-8 "${BOURNE_ARGS[@]}" \
    --out "$OUT/bourne-realizable-theirs"
step theirs-synthetic "$PY" benchmarks/synthetic/run_theirs.py --seeds 1 --method slsqp fw \
    "${THEIRS_ARGS[@]}" --force --out "$OUT/synthetic"
step bourne-approximate "$PY" benchmarks/bourne/run_comparison.py "${APPROX_ARGS[@]}" \
    --out "$OUT/bourne-approximate"
step report "$PY" benchmarks/report_overnight.py --dir "$OUT"

environment end
note "done; tables in $OUT/REPORT_TABLES.md, plots $OUT/timing.pdf (ours) and $OUT/timing-with-theirs.pdf"
