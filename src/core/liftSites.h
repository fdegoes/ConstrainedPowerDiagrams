#ifndef PP_CORE_LIFTSITES_H
#define PP_CORE_LIFTSITES_H

// Weighted sites as geogram Delaunay / restricted Voronoi diagram inputs
// (internal to the core library). Callers must hold a GeogramLock.

#include "core/domain.h"

#include <geogram/delaunay/delaunay.h>
#include <geogram/voronoi/RVD.h>

#include <string>
#include <vector>

namespace PP
{

// Builds the lifted sites (s_i, h_i) with h_i = sqrt(W - w_i), W = max_i w_i:
// the squared distance between (x, 0) and (s_i, h_i) is the power distance
// plus the constant W.
//
// dim+1 sentinel sites are appended at the corners of a simplex enclosing
// the sites and the domain, 10x larger than their bounding box. They keep
// the regular triangulation non-degenerate (geogram's RVD requires at least
// one simplex), and never own any part of the domain: with weight min_i w_i,
// their power distance is at least 81 L^2 - wmin whereas any real site's is
// at most L^2 - wmin (L = bounding box diagonal).
std::vector<double> LiftSites(
    const std::vector<double>& sites,
    const std::vector<double>& weights,
    const std::vector<double>& domainPts,
    int dim);

// Checks that there is one weight per site and that all values are finite.
bool CheckSites(
    const std::vector<double>& sites,
    const std::vector<double>& weights,
    int dim,
    std::string& err);

// Regular (weighted Delaunay) triangulation of the lifted sites and the
// restricted Voronoi diagram of the domain, i.e. the power diagram of the
// sites clipped to the domain. Seeds >= the number of sites are sentinels.
struct ClippedPowerDiagram
{
    // Lifted sites. geogram's Delaunay keeps a pointer to them (it does not
    // copy), so they must live as long as the triangulation.
    std::vector<double>               _lifted;
    GEO::Delaunay_var                 _delaunay;
    GEO::RestrictedVoronoiDiagram_var _rvd;
    // With parts: one diagram per part of the walk mesh
    // (Domain::GetWalkMesh()), restricted to its element range, instead of
    // _rvd.
    std::vector<GEO::RestrictedVoronoiDiagram_var> _parts;
    // Sites that are not vertices of the regular triangulation (their power
    // cell is empty), and whether one of them would seed geogram's RVD walk
    // (nearest to a domain vertex), which crashes the walk.
    int                               _numHidden = 0;
    bool                              _hiddenSeed = false;

    // Exact or inexact predicates for the walks of _rvd and _parts (exact by
    // default, see ComputePowerCells()).
    void SetExactPredicates(bool exact) const;
};

// Requires: CheckSites() succeeded and domain.IsValid(). Callers must not
// walk the RVD when _hiddenSeed is set. With parts, the diagrams are those
// of the parts of the walk mesh (_parts), which can be walked concurrently.
ClippedPowerDiagram MakeClippedPowerDiagram(
    const Domain& domain,
    const std::vector<double>& sites,
    const std::vector<double>& weights,
    bool parts = false);

} // namespace PP

#endif // PP_CORE_LIFTSITES_H
