#include "core/powerCells.h"
#include "core/cellBackend.h"
#include "core/cellMoments.h"
#include "core/domain.h"
#include "core/geogramContext.h"
#include "core/liftSites.h"

#include <cmath>
#include <thread>

namespace PP
{

namespace
{

// Number of sites from which the cells are integrated in parallel, over the
// parts of the domain's walk mesh: below, the threads cost more than they
// save.
constexpr size_t _minParallelSites = 2000;

// geogram's RVD walk can miss or duplicate a (cell, element) fragment at
// degenerate diagram vertices, e.g. where several cells meet on an edge of
// a domain element: with inexact predicates on lattices whose vertices are
// shared by four cells, and with exact predicates on some other diagrams.
// Such an error shows as a total cell volume different from the domain's.
constexpr double _volumeTolerance = 1e-12;

// Walks the RVD of pd (in parallel over its parts with parallel) into acc;
// false if the cells do not tile the domain.
bool
_WalkCells(
    const Domain& domain,
    const std::vector<double>& sites,
    const ClippedPowerDiagram& pd,
    bool parallel,
    Accumulator& acc)
{
    const CellBackend& backend = domain.GetBackend();
    if (not parallel) {
        backend.Walk(*pd._rvd, acc);
    } else {
        // One walk per part of the domain (they only read the
        // triangulation), merged in part order so that the result does not
        // depend on the thread timing.
        std::vector<Accumulator> parts(pd._parts.size() - 1, Accumulator(domain.GetDimension(), sites));
        std::vector<std::thread> threads;
        for (size_t p = 1; p < pd._parts.size(); ++p) {
            threads.emplace_back([&, p]() { backend.Walk(*pd._parts[p], parts[p - 1]); });
        }
        backend.Walk(*pd._parts[0], acc);
        for (std::thread& t : threads) {
            t.join();
        }
        for (const Accumulator& part : parts) {
            acc.Merge(part);
        }
    }
    return std::fabs(acc.GetTotalMass() - domain.GetVolume()) <= _volumeTolerance * domain.GetVolume();
}

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

    GeogramLock lock;
    const bool parallel = weights.size() >= _minParallelSites and domain.GetWalkRanges().size() > 2;
    ClippedPowerDiagram pd = MakeClippedPowerDiagram(domain, sites, weights, parallel);
    if (pd._hiddenSeed) {
        // geogram cannot walk this diagram; some cell is empty, so report
        // all cells as empty (the dual solvers reject such trial points) and
        // flag it.
        Accumulator(dim, sites).Fill(out);
        out._complete = false;
        return true;
    }
    // Each predicate mode fails on some degenerate diagrams: try exact
    // predicates, then inexact ones, and keep the walk closest to tiling the
    // domain.
    const double V = domain.GetVolume();
    double bestError = HUGE_VAL;
    for (bool exact : {true, false}) {
        pd.SetExactPredicates(exact);
        Accumulator trial(dim, sites);
        const bool tiles = _WalkCells(domain, sites, pd, parallel, trial);
        const double error = std::fabs(trial.GetTotalMass() - V);
        if (error < bestError) {
            bestError = error;
            trial.Fill(out);
        }
        if (tiles) {
            break;
        }
    }
    return true;
}

} // namespace PP
