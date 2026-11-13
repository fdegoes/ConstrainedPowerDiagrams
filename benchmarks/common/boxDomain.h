#ifndef PP_BENCHMARKS_BOXDOMAIN_H
#define PP_BENCHMARKS_BOXDOMAIN_H

// The clipping domain of an axis-aligned box, for the benchmark tools:
// a triangulated rectangle in 2D, a tetrahedralized box in 3D.

#include "core/domain.h"
#include "core/tetrahedralize.h"
#include "core/triangulate2D.h"

#include <string>
#include <vector>

namespace PPBench
{

// box  2 dim values: the lower corner, then the upper corner.
inline bool
MakeBoxDomain(int dim, const std::vector<double>& box, PP::Domain& domain, std::string& err)
{
    if (box.size() != size_t(2 * dim)) {
        err = "the box needs 2 dim values";
        return false;
    }
    if (dim == 2) {
        const double x0 = box[0], y0 = box[1], x1 = box[2], y1 = box[3];
        std::vector<double> pts;
        std::vector<int> tris;
        return PP::Triangulate2D({x0, y0, x1, y0, x1, y1, x0, y1}, {4}, pts, tris, err) and
               PP::Domain::Make(2, pts, tris, domain, err);
    }
    if (dim == 3) {
        const double x0 = box[0], y0 = box[1], z0 = box[2], x1 = box[3], y1 = box[4], z1 = box[5];
        // The 12 boundary triangles of the box (as tests/testUtils.h's unit cube).
        const std::vector<double> pts = {x0, y0, z0, x1, y0, z0, x1, y1, z0, x0, y1, z0,
                                         x0, y0, z1, x1, y0, z1, x1, y1, z1, x0, y1, z1};
        const std::vector<int> tris = {0, 2, 1, 0, 3, 2, 4, 5, 6, 4, 6, 7, 0, 1, 5, 0, 5, 4,
                                       2, 3, 7, 2, 7, 6, 1, 2, 6, 1, 6, 5, 0, 4, 7, 0, 7, 3};
        std::vector<double> tetPts;
        std::vector<int> tets;
        return PP::Tetrahedralize(pts, tris, PP::TetrahedralizeOptions(), tetPts, tets, err) and
               PP::Domain::Make(3, tetPts, tets, domain, err);
    }
    err = "the dimension must be 2 or 3";
    return false;
}

} // namespace PPBench

#endif // PP_BENCHMARKS_BOXDOMAIN_H
