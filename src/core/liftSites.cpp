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

void
ClippedPowerDiagram::SetExactPredicates(bool exact) const
{
    if (_rvd) {
        _rvd->set_exact_predicates(exact);
    }
    for (const GEO::RestrictedVoronoiDiagram_var& rvd : _parts) {
        rvd->set_exact_predicates(exact);
    }
}

ClippedPowerDiagram
MakeClippedPowerDiagram(
    const Domain& domain,
    const std::vector<double>& sites,
    const std::vector<double>& weights,
    bool parts)
{
    const int dim = domain.GetDimension();
    const CellBackend& backend = domain.GetBackend();
    ClippedPowerDiagram pd;
    pd._lifted = LiftSites(sites, weights, domain.GetPoints(), dim);
    GEO::Numeric::random_reset(-1);  // same insertion order at every build
    pd._delaunay = GEO::Delaunay::create(GEO::coord_index_t(dim + 1), backend.GetDelaunayName());
    // The RVD walks the triangulation through neighbor and incidence lists,
    // which must be requested before the vertices are set.
    pd._delaunay->set_stores_neighbors(true);
    pd._delaunay->set_vertices(GEO::index_t(pd._lifted.size() / (dim + 1)), pd._lifted.data());

    // Sites that are not vertices of the triangulation (empty cells).
    std::vector<char> present(pd._lifted.size() / (dim + 1), 0);
    for (GEO::index_t c = 0; c < pd._delaunay->nb_cells(); ++c) {
        for (GEO::index_t lv = 0; lv < pd._delaunay->cell_size(); ++lv) {
            const GEO::signed_index_t v = pd._delaunay->cell_vertex(c, lv);
            if (v >= 0 and size_t(v) < present.size()) {
                present[size_t(v)] = 1;
            }
        }
    }
    pd._numHidden = int(std::count(present.begin(), present.begin() + weights.size(), 0));
    // geogram's RVD walk seeds each connected piece of the domain mesh at its
    // first element (facet in 2D, tet in 3D), from the triangulation vertex
    // nearest to that element's first vertex, then propagates through
    // Delaunay neighbors, which are never hidden. A hidden seed makes the
    // walk read a missing incident cell and crash: replicate the seeds and
    // flag that case so that callers do not walk. With parts, each part's
    // walk seeds the connected pieces of its element range.
    GEO::Mesh& M = parts ? domain.GetWalkMesh() : domain.GetMesh();
    const bool volumetric = backend.IsVolumetric();
    const GEO::index_t nbElements = backend.GetNumElements(M);
    const std::vector<unsigned> ranges =
        parts ? domain.GetWalkRanges() : std::vector<unsigned>{0u, unsigned(nbElements)};
    if (pd._numHidden > 0) {
        std::vector<unsigned> part(nbElements, 0);
        for (size_t r = 0; r + 1 < ranges.size(); ++r) {
            std::fill(part.begin() + ranges[r], part.begin() + ranges[r + 1], unsigned(r));
        }
        std::vector<char> visited(nbElements, 0);
        std::vector<GEO::index_t> stack;
        std::vector<double> p(size_t(dim + 1), 0.0);
        for (GEO::index_t e = 0; e < nbElements and not pd._hiddenSeed; ++e) {
            if (visited[e]) {
                continue;
            }
            const GEO::index_t v0 = backend.GetElementVertex(M, e, 0);
            std::copy(M.vertices.point_ptr(v0), M.vertices.point_ptr(v0) + dim, p.begin());
            const GEO::index_t v = pd._delaunay->nearest_vertex(p.data());
            pd._hiddenSeed = size_t(v) < present.size() and not present[v];
            // Mark the connected piece of e.
            visited[e] = 1;
            stack.assign(1, e);
            while (not stack.empty()) {
                const GEO::index_t x = stack.back();
                stack.pop_back();
                const GEO::index_t nb = backend.GetNumAdjacent(M, x);
                for (GEO::index_t le = 0; le < nb; ++le) {
                    const GEO::index_t y = backend.GetAdjacent(M, x, le);
                    if (y != GEO::NO_INDEX and part[y] == part[x] and not visited[y]) {
                        visited[y] = 1;
                        stack.push_back(y);
                    }
                }
            }
        }
    }
    if (not parts) {
        pd._rvd = GEO::RestrictedVoronoiDiagram::create(pd._delaunay, &M);
        pd._rvd->set_volumetric(volumetric);
        pd.SetExactPredicates(true);
        return pd;
    }
    for (size_t r = 0; r + 1 < ranges.size(); ++r) {
        GEO::RestrictedVoronoiDiagram_var rvd = GEO::RestrictedVoronoiDiagram::create(pd._delaunay, &M);
        rvd->set_volumetric(volumetric);
        backend.SetElementRange(*rvd, ranges[r], ranges[r + 1]);
        pd._parts.push_back(rvd);
    }
    pd.SetExactPredicates(true);
    return pd;
}

} // namespace PP
