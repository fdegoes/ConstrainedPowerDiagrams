#include "core/powerCells.h"
#include "core/cellMoments.h"
#include "core/domain.h"
#include "core/geogramContext.h"
#include "core/liftSites.h"
#include "core/powerCells2D.h"
#include "core/powerCells3D.h"


namespace PP
{

namespace
{

// Number of sites from which the cells are integrated in parallel: below,
// the threads cost more than they save.
constexpr size_t _minParallelSites = 2000;

} // namespace

int
PowerCells::GetNumEmptyCells() const
{
    int count = 0;
    for (double v : _volumes) {
        count += (v <= 0.0);
    }
    return count;
}

bool
ComputePowerCells(
    const Domain& domain,
    const std::vector<double>& sites,
    const std::vector<double>& weights,
    PowerCells& out,
    std::string& err)
{
    out = PowerCells();
    if (not domain.IsValid()) {
        err = "Invalid clipping domain.";
        return false;
    }
    const int dim = domain.GetDimension();
    if (not CheckSites(sites, weights, dim, err)) {
        return false;
    }

    // Cells read from the regular triangulation and clipped one by one
    // against the domain elements (powerCells2D.h, powerCells3D.h).
    GeogramLock lock;
    const PowerTriangulation pd = MakePowerTriangulation(domain, sites, weights);
    Accumulator acc(dim, sites);
    const bool parallel = weights.size() >= _minParallelSites;
    if (dim == 2) {
        if (not ComputePowerCells2D(domain, sites, pd, parallel, acc)) {
            err = "The power cells could not be built (unbounded cell, or undefined power vertex of nearly collinear sites).";
            return false;
        }
    } else {
        ComputePowerCells3D(domain, sites, pd, parallel, acc);
    }
    acc.Fill(out);
    return true;
}

} // namespace PP
