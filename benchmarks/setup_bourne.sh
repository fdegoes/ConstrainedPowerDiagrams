#!/bin/bash
# Sets up what the comparisons with Bourne et al. need and git does not
# hold: the Python environment external/bourne/.venv with pysdot 0.2.39
# (patched, built from source) and the Bourne et al. repositories (git
# submodules). Linux or macOS; needs git, curl, python3 >= 3.12 (3.10 with
# UNPINNED=1) and a C++17 compiler.
#
#   bash benchmarks/setup_bourne.sh
#
# Environment:
#   PYTHON         Python used to create the venv (default: python3)
#   BOOST_INCLUDE  a dir with Boost headers (<dir>/boost/version.hpp), used
#                  instead of downloading Boost 1.87.0 (about 140 MB)
#   UNPINNED=1     install the latest numpy, scipy and matplotlib instead of
#                  benchmarks/bourne/requirements.txt (for an older Python)
#   FORCE=1        rebuild pysdot and the venv even if they exist
#
# Idempotent: the steps already done are skipped.
set -eu
cd "$(dirname "$0")/.."
ROOT=$PWD
PY=${PYTHON:-python3}
DIR=external/bourne
PYSDOT=$DIR/src/pysdot-0.2.39
VENV=$DIR/.venv

note() { echo "[setup_bourne] $*"; }

# The pinned packages (requirements.txt) need Python >= 3.12.
if ! "$PY" -c "import sys; sys.exit(sys.version_info < (3, 10))" 2> /dev/null; then
    echo "[setup_bourne] $PY is $("$PY" --version 2>&1); need Python >= 3.10 (3.12 for the pinned versions)." >&2
    echo "[setup_bourne] Set PYTHON=/path/to/python3.12 (or UNPINNED=1 with 3.10/3.11)." >&2
    exit 1
fi
if [ "${UNPINNED:-0}" != 1 ] && ! "$PY" -c "import sys; sys.exit(sys.version_info < (3, 12))"; then
    echo "[setup_bourne] $PY is $("$PY" --version 2>&1); the pinned versions need >= 3.12." >&2
    echo "[setup_bourne] Set PYTHON=/path/to/python3.12, or UNPINNED=1 for the latest versions." >&2
    exit 1
fi

# Pinned commits of the dependencies pysdot's setup.py clones (it would take
# their HEAD otherwise); those of the macOS environment of the paper.
SDOT_COMMIT=1b718eff36f910c1d9ce30c99d0ea9187cfb2add
PYBIND11_COMMIT=3beb3764ae967fa39dba1fe43b5a867f73da19ef
EIGEN_COMMIT=36b95962756c1fce8e29b1f8bc45967f30773c00

clone_at() {  # clone_at <url> <dir> <commit>
    if [ ! -d "$2/.git" ]; then
        git clone --quiet "$1" "$2"
    fi
    git -C "$2" checkout --quiet "$3"
}

# 1. Submodules: geogram (and its own) and the two Bourne et al. repos.
note "submodules"
git submodule update --init --recursive

if [ "${FORCE:-0}" = 1 ]; then
    rm -rf "$PYSDOT" "$VENV"
fi

# 2. pysdot 0.2.39 from its PyPI sdist, with our patch.
if [ ! -f "$PYSDOT/setup.py" ]; then
    note "pysdot 0.2.39 source"
    mkdir -p "$DIR/src"
    tmp=$(mktemp -d)
    "$PY" -m pip download --quiet --no-deps --no-binary :all: pysdot==0.2.39 -d "$tmp"
    tar xzf "$tmp/pysdot-0.2.39.tar.gz" -C "$DIR/src"
    rm -rf "$tmp"
    (cd "$PYSDOT" && patch -p1 --quiet < "$ROOT/benchmarks/bourne/pysdot-0.2.39.patch")
fi
mkdir -p "$PYSDOT/ext"
clone_at https://github.com/sd-ot/sdot.git "$PYSDOT/ext/sdot" "$SDOT_COMMIT"
clone_at https://github.com/pybind/pybind11.git "$PYSDOT/ext/pybind11" "$PYBIND11_COMMIT"
clone_at https://github.com/eigenteam/eigen-git-mirror.git "$PYSDOT/ext/eigen3" "$EIGEN_COMMIT"
if [ ! -e "$PYSDOT/boost_1_87_0/boost" ]; then
    if [ -n "${BOOST_INCLUDE:-}" ]; then
        note "Boost headers from $BOOST_INCLUDE"
        mkdir -p "$PYSDOT/boost_1_87_0"
        ln -s "$BOOST_INCLUDE/boost" "$PYSDOT/boost_1_87_0/boost"
    else
        note "Boost 1.87.0 headers (download)"
        curl -fsSL -o "$PYSDOT/boost_1_87_0.tar.gz" \
            https://archives.boost.io/release/1.87.0/source/boost_1_87_0.tar.gz
        tar xzf "$PYSDOT/boost_1_87_0.tar.gz" -C "$PYSDOT" boost_1_87_0/boost
        rm "$PYSDOT/boost_1_87_0.tar.gz"
    fi
fi

# 3. The venv, its packages, and pysdot built in it.
if [ ! -x "$VENV/bin/python" ]; then
    note "venv ($("$PY" --version))"
    "$PY" -m venv "$VENV"
fi
"$VENV/bin/python" -m pip install --quiet --upgrade pip setuptools wheel
if [ "${UNPINNED:-0}" = 1 ]; then
    "$VENV/bin/python" -m pip install --quiet numpy scipy matplotlib
else
    "$VENV/bin/python" -m pip install --quiet -r benchmarks/bourne/requirements.txt
fi
if ! "$VENV/bin/python" -c "import pysdot" 2> /dev/null; then
    note "build pysdot (a few minutes)"
    "$VENV/bin/python" -m pip install --quiet --no-build-isolation "./$PYSDOT"
fi

# 4. Check: the packages import and one OT solve runs.
export MPLCONFIGDIR=${MPLCONFIGDIR:-/tmp/mplcfg}
"$VENV/bin/python" - << 'EOF'
import numpy as np
import scipy, matplotlib
from pysdot import OptimalTransport
from pysdot.domain_types import ConvexPolyhedraAssembly
om = ConvexPolyhedraAssembly()
om.add_box([0, 0], [1, 1])
x = np.random.default_rng(0).random((20, 2))
ot = OptimalTransport(positions=x, masses=np.full(20, 1 / 20), domain=om, verbosity=0)
ot.adjust_weights()
err = np.max(np.abs(ot.pd.integrals() - 1 / 20))
print("numpy %s, scipy %s, matplotlib %s; pysdot OT volume error %.1e"
      % (np.__version__, scipy.__version__, matplotlib.__version__, err))
assert err < 1e-6
EOF
note "done: $VENV"
