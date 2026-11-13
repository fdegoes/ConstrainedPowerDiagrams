#ifndef PP_CORE_TRIANGULATE2D_H
#define PP_CORE_TRIANGULATE2D_H

#include <string>
#include <vector>

namespace PP
{

// Constrained Delaunay triangulation of the region bounded by closed
// polygonal loops (geogram's CDT2d).
//
// The region is defined by the even-odd rule: a point is inside if a ray
// from it crosses the loops an odd number of times. An outer loop plus
// loops inside it gives a domain with holes; several disjoint loops give
// several components. Loop orientation does not matter.
//
// The loops must not cross each other or themselves. Vertices shared by
// several loops are merged, but a loop must not repeat its first vertex
// at the end.
//
// loopPts   2 doubles per loop vertex, all loops concatenated
// loopSizes number of vertices of each loop
// triPts    2 doubles per output vertex (only the loop vertices
//           are used; no points are added)
// tris      3 indices per output triangle, counter-clockwise
bool Triangulate2D(
    const std::vector<double>& loopPts,
    const std::vector<int>& loopSizes,
    std::vector<double>& triPts,
    std::vector<int>& tris,
    std::string& err);

} // namespace PP

#endif // PP_CORE_TRIANGULATE2D_H
