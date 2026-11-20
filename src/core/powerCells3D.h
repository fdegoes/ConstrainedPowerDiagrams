#ifndef PP_CORE_POWERCELLS3D_H
#define PP_CORE_POWERCELLS3D_H

// 3D power cells (internal to the core library): each cell is the domain's
// bounding box cut by the power half-spaces of its neighbors in the regular
// triangulation, a convex polyhedron, then clipped against each domain
// tetrahedron it overlaps, independently of the other cells. Errors are
// bounded by round-off: no fragment of a cell can be missed or duplicated.

#include <vector>

namespace PP
{

class Accumulator;
class Domain;
struct PowerTriangulation;
struct PolyMesh;

// Integrates the cells and facets of the regular triangulation pd
// (MakePowerTriangulation()) of the sites into acc, in parallel over the
// sites with parallel. Requires a 3D domain.
void ComputePowerCells3D(
    const Domain& domain,
    const std::vector<double>& sites,
    const PowerTriangulation& pd,
    bool parallel,
    Accumulator& acc);

// Boundary polygons of the cells of pd (numSites sites), clipped to the 3D
// domain, appended to out (see PolyMesh): one polygon per facet between two
// cells, and per patch of the domain boundary in a cell, where a patch joins
// boundary pieces across edges whose dihedral angle is at most
// boundaryAngle degrees (0: coplanar pieces only). A facet or patch whose
// pieces do not merge into a single loop is output as its pieces. Polygons
// are oriented outwards.
void ComputePowerPolygons3D(
    const Domain& domain,
    int numSites,
    const PowerTriangulation& pd,
    double boundaryAngle,
    PolyMesh& out);

} // namespace PP

#endif // PP_CORE_POWERCELLS3D_H
