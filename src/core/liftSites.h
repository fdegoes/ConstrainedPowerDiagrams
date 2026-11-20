#ifndef PP_CORE_LIFTSITES_H
#define PP_CORE_LIFTSITES_H

// Weighted sites as inputs of geogram's regular triangulation (internal to
// the core library). Callers must hold a GeogramLock.

#include "core/domain.h"

#include <geogram/delaunay/delaunay.h>

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
// the regular triangulation non-degenerate (at least one simplex, with every
// real site inside its hull), and never own any part of the domain: with weight min_i w_i,
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

// Regular (weighted Delaunay) triangulation of the lifted sites, whose power
// cells powerCells2D.h and powerCells3D.h clip to the domain. Vertices >= the
// number of sites are sentinels.
struct PowerTriangulation
{
    // Lifted sites. geogram's Delaunay keeps a pointer to them (it does not
    // copy), so they must live as long as the triangulation.
    std::vector<double> _lifted;
    GEO::Delaunay_var   _delaunay;
};

// Requires: CheckSites() succeeded and domain.IsValid().
PowerTriangulation MakePowerTriangulation(
    const Domain& domain,
    const std::vector<double>& sites,
    const std::vector<double>& weights);

} // namespace PP

#endif // PP_CORE_LIFTSITES_H
