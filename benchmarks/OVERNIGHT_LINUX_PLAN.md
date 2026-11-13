# Plan: overnight realizable comparison on a Pixar Linux machine

This file is meant to be handed to a Claude Code session on the Linux
machine ("execute benchmarks/OVERNIGHT_LINUX_PLAN.md"). It clones the
repository, finds the dependencies in Pixar's builds, sets up the Python
environment, builds, runs a smoke test, then the full overnight comparison
of our H with the methods of Bourne et al. on all the realizable data, and
writes the report. Follow the steps in order. Stop and report to the user
if a step fails in a way the notes below do not cover.

What it compares is described in `benchmarks/OVERNIGHT.md`: the 3 realizable
Bourne cases and the 600 synthetic 2D cases, ours (H) and theirs (SLSQP and
Frank–Wolfe), each method on its own rule at 1e-8 with the other rule
measured at its solution, plus tables and timing plots. Read that file
first.

## 0. Rules for the session

- The timings are the point of the run. Do not run anything else heavy
  while it runs, and never run two benchmark steps in parallel.
- Work from the project root (the folder with `CMakeLists.txt`).
- Do not edit `src/`, `tests/` or the benchmark scripts to make things run.
  If a script fails, report it with the log instead. The only changes
  allowed are build configuration (cmake flags, environment variables),
  the Python environment, and dependencies built under `$HOME`.
- Do not use `sudo`, and do not install system packages. Pixar workstations
  are managed: use what Pixar's builds provide (step 3), and ask the user
  when something is missing.
- Network access (GitHub, PyPI, archives.boost.io) may go through a proxy
  or be blocked. If a download fails, stop and report the command and the
  error. Do not try other mirrors or routes on your own.
- Long steps go in the background (`nohup`), and their logs are checked
  periodically.
- Write down every path you choose (compiler, Python, Eigen, Boost,
  SuiteSparse). They go in the final report.

## 1. Check the machine

```sh
uname -sm; nproc; free -g; df -h "$HOME"
gcc --version; cmake --version; ninja --version || make --version
python3 --version; git --version; curl --version | head -1
uptime
```

The machine needs:
- a C++17 compiler (gcc ≥ 9);
- CMake ≥ 3.18, and ninja or make;
- Python ≥ 3.12 with venv and pip (3.10 or 3.11 work with `UNPINNED=1`,
  step 4);
- git and curl;
- 10 GB of free disk;
- and it should be idle, with a load average below 1.

If the system `gcc`, `cmake` or `python3` is too old, look for newer ones
in Pixar's builds (step 3) before asking the user. If the machine is
shared, report the load to the user before the full run.

## 2. Clone the repository

The repository is private:
`git@github.com:fdegoes/ConstrainedPowerDiagrams.git`. Cloning it needs the
user's SSH key for GitHub on this machine. If `git clone` fails with a
permission error, stop and ask the user to set up the key (or to clone it
themselves).

```sh
cd ~/Documents/Projects        # or wherever the user wants it
git clone --recursive git@github.com:fdegoes/ConstrainedPowerDiagrams.git
cd ConstrainedPowerDiagrams
git log -1 --oneline           # report the commit
```

`--recursive` fetches the submodules: geogram (and its own submodules) and
the two Bourne et al. repositories, pinned at the commits of the paper.
Check that they are there:

```sh
ls src/core/integrator.h benchmarks/run_overnight.sh benchmarks/setup_bourne.sh
ls external/geogram/CMakeLists.txt \
   external/bourne/Laguerre-Polycrystalline-Microstructures/Examples \
   external/bourne/Laguerre-fitting/EBSD_data.txt
```

If the repository was cloned without `--recursive`, run
`git submodule update --init --recursive`.

## 3. Find the dependencies in Pixar's builds

The build needs Eigen and SuiteSparse (CHOLMOD), and pysdot needs Boost
headers. The machine may not have them as system packages. Pixar's
software builds usually do. Search them in this order:

1. the VFX reference sets, `/pixar/d2/sets/` (e.g. `vfx2018.006`, which
   `config.sh` already uses for Eigen; prefer the newest set);
2. the software tree `/dist/sw/`;
3. the system (`/usr`, `/usr/local`), and `ldconfig -p`.

```sh
ls -d /pixar/d2/sets/*/ /dist/sw/*/ 2>/dev/null

# SuiteSparse: a prefix with include/cholmod.h (or include/suitesparse/cholmod.h)
# and lib/libcholmod.so (or lib64/)
find /pixar/d2/sets /dist/sw /usr /usr/local -maxdepth 5 -name cholmod.h 2>/dev/null
find /pixar/d2/sets /dist/sw /usr /usr/local -maxdepth 5 -name 'libcholmod.so*' 2>/dev/null
ldconfig -p | grep -i cholmod

# Boost headers: a dir <dir> with <dir>/boost/version.hpp
find /pixar/d2/sets /dist/sw /usr /usr/local -maxdepth 5 -path '*/boost/version.hpp' 2>/dev/null

# Eigen: a dir with Eigen/Core
find /pixar/d2/sets /dist/sw /usr /usr/local -maxdepth 5 -path '*/Eigen/Core' 2>/dev/null

# Newer compilers and Python, if the system ones are too old
ls /pixar/d2/sets/*/bin/{gcc,g++,cmake,python3*} /dist/sw/*/bin/python3* 2>/dev/null
```

`find` over these trees can be slow. Raise `-maxdepth` only if nothing
shows up.

Known good (surveyed on 2026-11-13, RHEL 9 workstation): everything is in
`/pixar/d2/sets/vfx2025.002`: SuiteSparse 5.10 (prefix = the set), Boost 1.85
(`include/`), `bin/python3.11` (use `UNPINNED=1`), cmake 3.30 and ninja. `config.sh`
runs that set's cmake itself (the `cmake` on `PATH` may be 3.2.1, which has no
`--build -j`: build with `make -C build -j"$(nproc)" <targets>`). The VFX sets only have Eigen 3.2.8, too old: use
`$HOME/mylibs/eigen337` (3.3.7). `config.sh` picks these by default; check
that its first line shows them, and still record them in the report.

How to choose:

- **SuiteSparse.** The prefix must have `cholmod.h` and all of
  `libcholmod`, `libamd`, `libcamd`, `libcolamd`, `libccolamd` and
  `libsuitesparseconfig`. CMake checks this. Check the version in
  `SuiteSparse_config.h`; any 5.x or later works. That prefix is
  `SUITESPARSE_PATH` in step 5.
- **Boost.** Any version from about 1.70 works with the patch. Write down
  its version (`grep BOOST_LIB_VERSION <dir>/boost/version.hpp`). That
  `<dir>` is `BOOST_INCLUDE` in step 4. Houdini's bundled Boost is renamed
  to `hboost` and cannot be used.
- **Eigen.** Needs 3.3 or later (`grep EIGEN_WORLD_VERSION -A2
  <dir>/Eigen/src/Core/util/Macros.h`). `config.sh` takes the VFX set's
  Eigen if it exists; otherwise set `EIGEN3_PATH`.

If nothing is found:

- **SuiteSparse:** report it to the user and ask before choosing.
  - The first option is to build SuiteSparse 5.x under `$HOME`. This needs
    a BLAS and LAPACK; look for OpenBLAS or MKL in the same trees.
  - The fallback is the portable Eigen backend (`PP_SPARSE_SOLVER=Eigen`).
    It is slower, so record the choice: it goes in the report, and the
    timings are not comparable with a CHOLMOD run.
- **Boost:** without `BOOST_INCLUDE`, `setup_bourne.sh` downloads Boost
  1.87.0 from archives.boost.io (140 MB). If that download is blocked, stop
  and ask the user.

## 4. The Python environment with pysdot

```sh
PYTHON=python3.12 BOOST_INCLUDE=<boost dir> bash benchmarks/setup_bourne.sh
```

- **`PYTHON`:** a Python ≥ 3.12 (from the system or step 3). With 3.10 or
  3.11, add `UNPINNED=1`. That installs the latest numpy, scipy and
  matplotlib instead of the pinned versions, so record it.
- **`BOOST_INCLUDE`:** the Boost dir from step 3. If you leave it out, the
  script downloads Boost.

The script does these steps:
1. fetches the submodules;
2. downloads the pysdot 0.2.39 sdist from PyPI and applies
   `benchmarks/bourne/pysdot-0.2.39.patch`;
3. clones pysdot's dependencies (sdot, pybind11, Eigen) at pinned commits;
4. creates `external/bourne/.venv` with `benchmarks/bourne/requirements.txt`;
5. builds pysdot, which takes a few minutes;
6. checks the result with one OT solve.

It ends with `[setup_bourne] done`, after a line with the package versions
and an OT volume error below 1e-6. You can re-run it, and it skips the
steps already done. `FORCE=1` rebuilds pysdot and the venv.

If pysdot fails to compile, check the compiler (C++17), then the Boost
version. If the patch does not apply, stop and report it.

Then check that the drivers see their cases:

```sh
export MPLCONFIGDIR=/tmp/mplcfg
external/bourne/.venv/bin/python benchmarks/bourne/run_realizable.py --list
# expect: 2025-invert (n=20), 2026-voronoi (n=20), 2026-aligned (n=728)
```

## 5. Configure and build (no Houdini)

`config.sh` detects Linux. Pass it the paths from step 3:

```sh
PP_BUILD_HOUDINI=OFF \
SUITESPARSE_PATH=<suitesparse prefix> \
EIGEN3_PATH=<eigen dir> \
CC=<gcc> CXX=<g++> \
sh config.sh
make -C build -j"$(nproc)" pp_core benchBourne makeSynthetic \
    testPowerDiagram testVolumeSolver testJointDual testDualSolvers
```

- **Optional variables:** `EIGEN3_PATH` can be left out if the VFX set's
  Eigen is the one to use. `CC`/`CXX` can be left out if the system gcc is
  new enough.
- **Missing SuiteSparse:** if `SUITESPARSE_PATH` is wrong, CMake stops with
  "SuiteSparse: lib… not found" or "cholmod.h not found". Fix the path. Do
  not switch backends without asking the user (step 3).
- **Generator:** `config.sh` uses make on Linux. If ninja is available, it
  is faster: add `-G Ninja` by running the `cmake -S . -B build` line from
  `config.sh` by hand.

Check that the binaries use CHOLMOD and find it at run time:

```sh
grep '^PP_SPARSE_SOLVER:' build/CMakeCache.txt     # Cholmod
ldd build/benchBourne | grep -E 'cholmod|not found' # a path, no "not found"
```

Then run the unit tests. All four must print `PASSED`:

```sh
for t in testPowerDiagram testVolumeSolver testJointDual testDualSolvers; do build/$t | tail -1; done
```

If one fails, stop and report its output. Do not continue to the
benchmarks.

## 6. Smoke test (about 5 minutes)

```sh
OVERNIGHT_QUICK=1 bash benchmarks/run_overnight.sh
cat /tmp/overnight-quick/overnight.log
cat /tmp/overnight-quick/environment.txt
```

Expect every step to end with `exit 0`, `/tmp/overnight-quick/REPORT_TABLES.md`
to exist, and `environment.txt` to show the CPU and `sparse solver:
Cholmod`. If a step says `FAILED`, read its log in
`/tmp/overnight-quick/<step>.log`:

- **build:** fix the cmake configuration (step 5).
- **ours-* or theirs-*:** a Python error. Check the environment (step 4).
  If it is in our scripts, stop and report it.
- **report:** matplotlib. Make sure `MPLCONFIGDIR` is writable.

## 7. The full run (worst case about 5 h)

The machine should be idle: nothing else running, no sleep. On Linux there
is no `caffeinate`, so use `systemd-inhibit` if it is available, otherwise
plain `nohup`:

```sh
mkdir -p benchmarks/overnight
if command -v systemd-inhibit >/dev/null; then
    nohup systemd-inhibit --what=sleep:idle --why="overnight benchmark" \
        bash benchmarks/run_overnight.sh > benchmarks/overnight/run.log 2>&1 &
else
    nohup bash benchmarks/run_overnight.sh > benchmarks/overnight/run.log 2>&1 &
fi
echo $! > benchmarks/overnight/run.pid
```

Monitor it every 20–30 minutes (not more often, to keep the machine idle):

```sh
tail -n 5 benchmarks/overnight/$(date +%Y-%m-%d)/overnight.log
```

The steps and their usual durations are:
1. build;
2. ours on the Bourne cases (under 1 min);
3. ours on the synthetic cases (about 20 min);
4. theirs on the Bourne cases (up to 30 min per case);
5. theirs on the synthetic cases (up to 4.4 h);
6. the report.

If the run started just before midnight, the output folder has the date of
the start.

## 8. When it finishes

Output: `benchmarks/overnight/<date>/`, described in `benchmarks/OVERNIGHT.md`.

1. Check `overnight.log`: every step should have `exit 0`. List any
   `FAILED` steps and the end of their logs.
2. Check `environment.txt`: the host, the CPU, the sparse solver, and the
   load at the start and at the end.
3. Summarize `REPORT_TABLES.md` for the user:
   - the Bourne realizable cases: iterations, time and residuals, ours and
     theirs;
   - the synthetic cases: convergence of our H over the 600 cases, and the
     matched comparison with theirs;
   - the runs that hit their caps, with their residuals.
4. Point to `timing.pdf` and `timing-with-theirs.pdf`.
5. List the paths and versions chosen in steps 3–5:
   - compiler, Python, Eigen, Boost, SuiteSparse;
   - the sparse solver;
   - whether `UNPINNED=1` was used.
6. Package the results for the user's Mac. Everything except the per-case
   data, which can be regenerated, is small:

   ```sh
   cd benchmarks/overnight
   tar czf overnight-$(date +%Y-%m-%d)-linux.tar.gz \
       --exclude='synthetic/noise*' "$(ls -d 20* | tail -1)"
   ```

   Report the tar file's path and size.

## Notes

- **Timings differ from the Mac runs** (another CPU, BLAS and compiler, even
  with the same CHOLMOD backend). Compare methods within this run, not
  against numbers from the Mac.
- **Known behaviour, not failures:**
  - `2026-aligned`: our H does not converge (it ends with "damping failed"
    or at its iteration cap). Its targets are realizable to round-off, so
    this is a solver issue under investigation. Report its residuals.
  - Their FW: it usually runs to its cap, because FW with t = 100 solves
    H_α with α = 99, whose H residual stalls at about 1/t.
  - `testIntegrators` is not run here. It has one known failing check (a
    2D constrained step), unrelated to this comparison.
- **Approximate cases:** this plan does not include them. Their
  head-to-head is `benchmarks/bourne/run_comparison.py`
  (`benchmarks/bourne/README_COMPARISON.md`), which needs the same build
  and Python environment and can be run afterwards the same way.
