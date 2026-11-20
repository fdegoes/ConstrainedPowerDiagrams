#include "core/liftSites.h"
#include "core/cellBackend.h"
#include "core/geometry.h"

#include <geogram/basic/numeric.h>

#include <algorithm>
#include <cmath>

namespace PP
{

namespace
{

// Unit directions to the corners of a regular simplex centered at 0.
std::vector<std::vector<double>>
_GetSimplexDirections(int dim)
{
    if (dim == 2) {
        std::vector<std::vector<double>> dirs;
        for (int k = 0; k < 3; ++k) {
            const double a = 2.0 * M_PI * k / 3.0;
            dirs.push_back({std::cos(a), std::sin(a)});
        }
        return dirs;
    }
    const double r = 1.0 / std::sqrt(3.0);
    return {{r, r, r}, {r, -r, -r}, {-r, r, -r}, {-r, -r, r}};
}

} // namespace

std::vector<double>
LiftSites(
    const std::vector<double>& sites,
    const std::vector<double>& weights,
    const std::vector<double>& domainPts,
    int dim)
{
    const double wmin = *std::min_element(weights.begin(), weights.end());
    const double wmax = *std::max_element(weights.begin(), weights.end());

    BoundingBox box;
    box.Add(sites, dim);
    box.Add(domainPts, dim);
    const double L = std::max(box.GetDiagonal(dim), 1e-6);

    std::vector<double> lifted;
    lifted.reserve((dim + 1) * (weights.size() + dim + 1));
    for (size_t i = 0; i < weights.size(); ++i) {
        lifted.insert(lifted.end(), &sites[dim * i], &sites[dim * i] + dim);
        lifted.push_back(std::sqrt(std::max(0.0, wmax - weights[i])));
    }
    for (const std::vector<double>& d : _GetSimplexDirections(dim)) {
        for (int c = 0; c < dim; ++c) {
            lifted.push_back(box.GetCenter(c) + 10.0 * L * d[c]);
        }
        lifted.push_back(std::sqrt(std::max(0.0, wmax - wmin)));
    }
    return lifted;
}

bool
CheckSites(
    const std::vector<double>& sites,
    const std::vector<double>& weights,
    int dim,
    std::string& err)
{
    if (weights.empty()) {
        err = "No sites.";
        return false;
    }
    if (sites.size() != size_t(dim) * weights.size()) {
        err = "Number of sites and weights do not match.";
        return false;
    }
    if (not AreAllFinite(sites) or not AreAllFinite(weights)) {
        err = "Sites or weights are not finite.";
        return false;
    }
    return true;
}

PowerTriangulation
MakePowerTriangulation(
    const Domain& domain,
    const std::vector<double>& sites,
    const std::vector<double>& weights)
{
    const int dim = domain.GetDimension();
    PowerTriangulation pd;
    pd._lifted = LiftSites(sites, weights, domain.GetPoints(), dim);
    GEO::Numeric::random_reset(-1);  // same insertion order at every build
    pd._delaunay = GEO::Delaunay::create(GEO::coord_index_t(dim + 1), domain.GetBackend().GetDelaunayName());
    // The cells are read through the neighbor and incidence lists, which
    // must be requested before the vertices are set.
    pd._delaunay->set_stores_neighbors(true);
    pd._delaunay->set_vertices(GEO::index_t(pd._lifted.size() / (dim + 1)), pd._lifted.data());
    return pd;
}

} // namespace PP
