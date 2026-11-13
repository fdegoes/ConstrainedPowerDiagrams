#!/bin/sh
# Configure the build tree in ./build (Release), on macOS or Linux.
#
#   sh config.sh             # installs into the Houdini user dir
#                            # (~/Library/Preferences/houdini/22.0)
#   sh config.sh <user_dir>  # installs into <user_dir>/dso and <user_dir>/help
#   cmake --build build -j8
#   ./build/testPowerDiagram
#
# The sparse Cholesky backend is CHOLMOD on macOS and Linux; set
# PP_SPARSE_SOLVER=Accelerate (or Eigen) in the environment to override it,
# and BUILD_DIR to configure another tree (e.g. to compare backends).
# PP_BUILD_HOUDINI=OFF builds without the Houdini plugin (tests and
# benchmarks only; no Houdini install needed). EIGEN3_PATH and
# SUITESPARSE_PATH override the Eigen and SuiteSparse locations.
set -e

EXTRA=""
if [ -n "${1:-}" ]; then
    EXTRA="-DPP_HOUDINI_HOME=$1"
fi
: "${BUILD_DIR:=build}"
: "${PP_BUILD_HOUDINI:=ON}"
: "${EIGEN3_PATH:=}"
EXTRA="$EXTRA -DPP_BUILD_HOUDINI=${PP_BUILD_HOUDINI}"

if [ "$(uname -s)" = "Linux" ]
then
    # Linux (Pixar workstations). The compilers are CC/CXX if set, else the
    # system's; they need C++17 (geogram). Eigen defaults to the VFX set's.
    : "${HFS:=/dist/sw/houdini/22.0}"
    VFX="/pixar/d2/sets/vfx2018.006"
    if [ -z "${EIGEN3_PATH}" ] && [ -f "${VFX}/include/eigen3/Eigen/Core" ]; then
        EIGEN3_PATH="${VFX}/include/eigen3"
    fi
    # TODO: SuiteSparse install prefix (include/cholmod.h, lib/libcholmod.so).
    : "${SUITESPARSE_PATH:=}"
    : "${PP_SPARSE_SOLVER:=Cholmod}"
    export HFS
    cmake -S . -B "${BUILD_DIR}" \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_EXPORT_COMPILE_COMMANDS=ON \
        -DEIGEN3_INCLUDE_DIR="${EIGEN3_PATH:-$HOME/MyLibs/eigen-3.4.0}" \
        -DPP_SPARSE_SOLVER="${PP_SPARSE_SOLVER}" \
        -DSUITESPARSE_DIR="${SUITESPARSE_PATH}" \
        $EXTRA \
        -Wno-dev
else
    # macOS config (arm64).
    : "${HFS:=/Applications/Houdini/Current/Frameworks/Houdini.framework/Versions/Current/Resources}"
    : "${SUITESPARSE_PATH:=$HOME/MyLibs/SuiteSparse-5.4-arm64}"
    : "${PP_SPARSE_SOLVER:=Cholmod}"
    : "${EIGEN3_PATH:=$HOME/MyLibs/eigen-3.4.0}"
    export HFS
    cmake -S . -B "${BUILD_DIR}" -G Ninja \
        -DEIGEN3_INCLUDE_DIR="${EIGEN3_PATH}" \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_OSX_ARCHITECTURES=arm64 \
        -DCMAKE_EXPORT_COMPILE_COMMANDS=ON \
        -DPP_SPARSE_SOLVER="${PP_SPARSE_SOLVER}" \
        -DSUITESPARSE_DIR="${SUITESPARSE_PATH}" \
        $EXTRA \
        -Wno-dev
fi
