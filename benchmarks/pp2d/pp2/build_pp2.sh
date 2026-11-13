#!/bin/sh
# Builds the original 2D Power Particles code (Before2015/pp/code/pp2)
# headless, with its bundled voro++, into bin/pp2. The sources are copied to
# a scratch directory and not modified. See ../README.md.
#
#   sh benchmarks/pp2d/pp2/build_pp2.sh [pp2_src] [cgal_dir] [boost_include]
set -e
HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=$(cd "$HERE/../../.." && pwd)
SRC=${1:-$HOME/Documents/Projects/Before2015/pp/code/pp2}
CGAL=${2:-$ROOT/external/cgal-6.0.1}
BOOST=${3:-/Applications/RV64.app/Contents/include}
WORK=$(mktemp -d)
cp -R "$SRC/." "$WORK/src"
INC=""
for d in "$CGAL"/*/include; do INC="$INC -I$d"; done
# -ffp-contract=off: without it voro++'s radical-cell bounds test stops the
# neighbor search too early (overlapping cells, diverging solves).
FLAGS="-std=c++17 -O3 -ffp-contract=off -DNDEBUG -DCGAL_NDEBUG -D_LIBCPP_ENABLE_CXX17_REMOVED_UNARY_BINARY_FUNCTION -DCGAL_DISABLE_GMP -DCGAL_NO_CORE -I$WORK/src -isystem $BOOST $INC -w"
mkdir -p "$WORK/obj"
for f in voropp/voro++_2d.cc Dual.cpp FastBN.cpp Grid.cpp KdTree.cpp Marker.cpp Parameters.cpp Particles.cpp \
         ParticlesIO.cpp Polygon.cpp Sampler.cpp Sim.cpp SimAir.cpp SimBN.cpp SimDEC.cpp SimDyn.cpp SimInit.cpp \
         SimIO.cpp SimMesh.cpp main.cpp; do
    c++ $FLAGS -c "$WORK/src/$f" -o "$WORK/obj/$(basename "$f").o"
done
mkdir -p "$HERE/bin"
c++ "$WORK"/obj/*.o -framework Accelerate -o "$HERE/bin/pp2"
echo "$FLAGS" > "$HERE/bin/FLAGS"
rm -rf "$WORK"
echo "built $HERE/bin/pp2"
