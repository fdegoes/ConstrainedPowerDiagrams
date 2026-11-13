#ifndef PP_CORE_TETRAHEDRALIZE_H
#define PP_CORE_TETRAHEDRALIZE_H

#include <string>
#include <vector>

namespace PP
{

struct TetrahedralizeOptions
{
    // Merge duplicated vertices, fill small holes, remove intersections and
    // keep the outer boundary before meshing.
    bool _preprocess = false;

    // Insert Steiner points to improve tet quality.
    bool _refine = false;
    
    // Target radius-edge ratio when refining, typically in [1, 2]
    // (smaller is better quality and more tets).
    double _quality = 2.0;
};

// Fills the volume bounded by a triangulated surface with tetrahedra
// (TetGen through geogram's mesh_tetrahedralize).
//
// Without preprocessing, the surface must be closed, 2-manifold, welded
// (shared edges reference the same vertex indices) and free of
// self-intersections and degenerate triangles. The surface may be
// non-convex, have several components or nested shells (cavities).
//
// surfPts  3 doubles per surface vertex
// surfTris 3 indices per triangle
// tetPts   3 doubles per output vertex. The input vertices come
//          first, in order, unless preprocessing modified them.
// tets     4 indices per output tet, positively oriented
bool Tetrahedralize(
    const std::vector<double>& surfPts,
    const std::vector<int>& surfTris,
    const TetrahedralizeOptions& options,
    std::vector<double>& tetPts,
    std::vector<int>& tets,
    std::string& err);

} // namespace PP

#endif // PP_CORE_TETRAHEDRALIZE_H
