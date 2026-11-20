#ifndef PP_CORE_POWERCELLS2D_H
#define PP_CORE_POWERCELLS2D_H

// 2D power cells (internal to the core library): each cell is read from the
// regular triangulation as a convex polygon and clipped against the domain
// triangles it overlaps, independently of the other cells, so that no
// fragment can be missed or duplicated, even when diagram vertices lie on
// the edges of domain triangles (e.g. lattice sites), and without walking
// through the domain mesh as geogram's restricted Voronoi diagram does.

#include "core/powerCells.h"

#include <vector>

namespace PP
{

class Accumulator;
class Domain;
struct PowerTriangulation;
struct PolyMesh;

// Integrates the cells and facets of the regular triangulation pd
// (MakePowerTriangulation()) of the sites into acc, in parallel over the
// sites with parallel. Requires a 2D domain. False if a cell is unbounded
// or a power vertex is undefined (not expected: the sentinel sites of
// LiftSites() bound every cell).
bool ComputePowerCells2D(
    const Domain& domain,
    const std::vector<double>& sites,
    const PowerTriangulation& pd,
    bool parallel,
    Accumulator& acc);

// Boundary loops of the cells of the regular triangulation pd of numSites
// sites, clipped to the 2D domain, appended to out (see PolyMesh): one loop
// per boundary component of each non-empty cell, counter-clockwise for outer
// boundaries and clockwise for holes. False as ComputePowerCells2D().
bool ComputePowerPolygons2D(const Domain& domain, int numSites, const PowerTriangulation& pd, PolyMesh& out);

} // namespace PP

#endif // PP_CORE_POWERCELLS2D_H
