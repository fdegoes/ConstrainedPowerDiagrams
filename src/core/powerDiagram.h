#ifndef PP_CORE_POWERDIAGRAM_H
#define PP_CORE_POWERDIAGRAM_H

#include <string>
#include <vector>

namespace PP
{

class Domain;

// Polygonal output of a clipped power diagram.
//
// Points are stored with `dim` coordinates each (2 for the 2D diagram,
// 3 for the 3D one). Points are shared among the polygons of a same cell,
// never across cells. Polygons are oriented counter-clockwise (2D) or with
// outward normals by the right-hand rule (3D).
struct PolyMesh
{
    int                 _dim = 0;
    std::vector<double> _points;
    std::vector<int>    _polySizes;  // number of vertices per polygon
    std::vector<int>    _polyVerts;  // concatenated point indices
    std::vector<int>    _polyCell;   // site index for each polygon

    size_t GetNumPoints() const { return _dim ? _points.size() / _dim : 0; }
    size_t GetNumPolygons() const { return _polySizes.size(); }
};

// Power diagram of weighted sites clipped to a prebuilt domain, of either
// dimension: the 2D or 3D output below, with the domain's requirements.
// sites          d doubles per site (d = domain.GetDimension())
// weights        1 double per site
// boundaryAngle  3D only, as in ComputePowerDiagram3D()
bool ComputePowerDiagram(
    const Domain& domain,
    const std::vector<double>& sites,
    const std::vector<double>& weights,
    double boundaryAngle,
    PolyMesh& out,
    std::string& err);

// Power diagram of 2D weighted sites clipped to a triangulated domain. The
// power distance of x to site i is |x - sites_i|^2 - weights_i. Each power
// cell is clipped against the triangles it overlaps (powerCells2D.h), and
// the pieces of a cell are merged across shared edges. Use Triangulate2D() to build the domain from polygons.
//
// Requirements on the triangle mesh (not checked; the caller must
// guarantee them, otherwise the output is undefined):
//  - conforming: two triangles are either disjoint or share a whole vertex
//    or edge (no T-junctions, no overlaps);
//  - welded: triangles sharing an edge reference the same vertex indices;
//  - no degenerate (zero-area) triangles.
// Triangle orientation does not matter. The domain may be non-convex,
// have holes or several components.
//
// Sites may lie outside the domain and may coincide; sites whose cell
// does not meet the domain (outside, or hidden by heavier neighbors)
// produce no output.
//
// sites   2 doubles per site
// weights 1 double per site
// triPts  2 doubles per domain vertex
// tris    3 vertex indices per domain triangle
// out     one polygon per boundary loop of each non-empty cell:
//         counter-clockwise for outer boundaries, clockwise for the
//         boundaries of holes inside a cell
bool ComputePowerDiagram2D(
    const std::vector<double>& sites,
    const std::vector<double>& weights,
    const std::vector<double>& triPts,
    const std::vector<int>& tris,
    PolyMesh& out,
    std::string& err);

// Power diagram of 3D weighted sites clipped to the volume of a tetrahedral
// mesh. Each power cell is clipped against the tetrahedra it overlaps
// (powerCells3D.h), and the pieces of a cell are merged across shared tet
// faces.
//
// Requirements on the tetrahedral mesh (not checked; the caller must
// guarantee them, otherwise the output is undefined):
//  - conforming: two tets are either disjoint or share a whole vertex, edge
//    or face (no T-junctions, no overlaps);
//  - welded: tets sharing a face reference the same vertex indices;
//  - no degenerate (zero-volume) tets.
// Tet orientation does not matter (tets are reoriented internally). The
// domain may be non-convex, have several components, or cavities.
//
// Sites may lie outside the volume; sites whose cell does not meet the
// volume produce no output.
//
// sites     3 doubles per site
// weights   1 double per site
// tetPts    3 doubles per tet mesh vertex
// tets      4 vertex indices per tetrahedron
// boundaryAngle  degrees: the boundary pieces of a cell (its clipped
//           boundary tet faces) that share an edge and whose normals
//           differ by at most this angle are merged into one polygon.
//           0 merges coplanar pieces only; a positive value also merges
//           across edges of curved domains, into polygons that need not
//           be planar. Pieces that do not merge into a single loop are
//           kept apart.
// out       closed boundary polygons of each non-empty cell, with
//           outward normals (right-hand rule)
bool ComputePowerDiagram3D(
    const std::vector<double>& sites,
    const std::vector<double>& weights,
    const std::vector<double>& tetPts,
    const std::vector<int>& tets,
    double boundaryAngle,
    PolyMesh& out,
    std::string& err);

} // namespace PP

#endif // PP_CORE_POWERDIAGRAM_H
