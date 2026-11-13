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
    # system's; they need C++17 (geogram). SuiteSparse and Eigen default to the
    # first match below; SUITESPARSE_PATH and EIGEN3_PATH override them.
    : "${HFS:=/dist/sw/houdini/22.0}"
    : "${PP_SPARSE_SOLVER:=Cholmod}"
    VFX="/pixar/d2/sets/vfx2025.002"

    # SuiteSparse prefix: needs include/cholmod.h and lib/libcholmod.so.
    if [ "${PP_SPARSE_SOLVER}" = "Cholmod" ] && [ -z "${SUITESPARSE_PATH:-}" ]; then
        for d in "${VFX}" "$HOME/mylibs/SuiteSparse"; do
            if [ -f "$d/include/cholmod.h" ] && { [ -f "$d/lib/libcholmod.so" ] || [ -f "$d/lib64/libcholmod.so" ]; }; then
                SUITESPARSE_PATH="$d"
                break
            fi
        done
        if [ -z "${SUITESPARSE_PATH:-}" ]; then
            echo "config.sh: no SuiteSparse found (tried ${VFX}, ~/mylibs/SuiteSparse)." >&2
            echo "config.sh: set SUITESPARSE_PATH=<prefix>, or PP_SPARSE_SOLVER=Eigen." >&2
            exit 1
        fi
    fi
    : "${SUITESPARSE_PATH:=}"

    # Eigen >= 3.3. The VFX sets only have 3.2.8, so they are used only if newer.
    if [ -z "${EIGEN3_PATH}" ]; then
        for d in "${VFX}/include/eigen3" "$HOME/mylibs/eigen337"; do
            m="$d/Eigen/src/Core/util/Macros.h"
            if [ -f "$m" ] && [ "$(sed -n 's/^#define EIGEN_WORLD_VERSION //p' "$m")" -ge 3 ] \
                && [ "$(sed -n 's/^#define EIGEN_MAJOR_VERSION //p' "$m")" -ge 3 ]; then
                EIGEN3_PATH="$d"
                break
            fi
        done
    fi
    echo "config.sh: SuiteSparse = ${SUITESPARSE_PATH:-none}, Eigen = ${EIGEN3_PATH:-default}"

    # cmake >= 3.18: the VFX set's, unless CMAKE is set (the cmake on PATH may be 3.2).
    if [ -z "${CMAKE:-}" ] && [ -x "${VFX}/bin/cmake" ]; then
        CMAKE="${VFX}/bin/cmake"
    fi
    : "${CMAKE:=cmake}"
    echo "config.sh: cmake = $CMAKE"
    export HFS
    "${CMAKE}" -S . -B "${BUILD_DIR}" \
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
