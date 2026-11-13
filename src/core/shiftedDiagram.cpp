#include "core/shiftedDiagram.h"
#include "core/domain.h"
#include "core/geometry.h"

#include <algorithm>
#include <cmath>

namespace PP
{

namespace
{

bool
_CheckShiftedInputs(
    const Domain& domain,
    const std::vector<double>& sites,
    const std::vector<double>& targets,
    const std::vector<double>& shifts,
    const std::vector<double>& weights,
    std::string& err)
{
    if (not domain.IsValid()) {
        err = "Invalid clipping domain.";
        return false;
    }
    const size_t n = weights.size();
    const size_t dim = size_t(domain.GetDimension());
    if (n == 0) {
        err = "No sites.";
        return false;
    }
    if (sites.size() != dim * n or targets.size() != dim * n or shifts.size() != dim * n) {
        err = "Number of sites, targets, shifts and weights do not match.";
        return false;
    }
    if (not AreAllFinite(targets) or not AreAllFinite(shifts)) {
        err = "Targets or shifts are not finite.";
        return false;
    }
    return true;
}

} // namespace

bool
ShiftedDiagram::Build(
    const Domain& domain,
    const std::vector<double>& sites,
    const std::vector<double>& targets,
    const std::vector<double>& shifts,
    const std::vector<double>& weights,
    std::string& err)
{
    *this = ShiftedDiagram();
    if (not _CheckShiftedInputs(domain, sites, targets, shifts, weights, err)) {
        return false;
    }
    const int dim = domain.GetDimension();
    const int n = int(weights.size());
    _targets = targets;
    _generators.resize(dim * n);
    _shiftedWeights.resize(n);
    for (int i = 0; i < n; ++i) {
        const double* q = &sites[dim * i];
        const double* u = &shifts[dim * i];
        const double* cbar = &targets[dim * i];
        double lambda = weights[i];
        for (int c = 0; c < dim; ++c) {
            _generators[dim * i + c] = q[c] + u[c];
            lambda += u[c] * u[c] + 2.0 * u[c] * (q[c] - cbar[c]);
        }
        _shiftedWeights[i] = lambda;
    }
    if (not ComputePowerCells(domain, _generators, _shiftedWeights, _cells, err)) {
        return false;
    }
    _distances.resize(_cells._facets.size());
    for (size_t f = 0; f < _cells._facets.size(); ++f) {
        const Facet& F = _cells._facets[f];
        _distances[f] = std::sqrt(ComputeDistance2(GetGenerator(F._i), GetGenerator(F._j), dim));
    }
    return true;
}

double
ShiftedDiagram::GetMinVolume() const
{
    return *std::min_element(_cells._volumes.begin(), _cells._volumes.end());
}

double
ShiftedDiagram::GetEnergy() const
{
    // _secondMoments are about the sites given to ComputePowerCells, i.e.
    // the generators.
    double e = 0.0;
    for (int i = 0; i < GetNumSites(); ++i) {
        e += _cells._secondMoments[i] - _shiftedWeights[i] * _cells._volumes[i];
    }
    return e;
}

void
ShiftedDiagram::GetShiftShiftBlock(size_t f, int a, int b, double* B, double minDistance) const
{
    // ∫ (x - ca)(x - cb)^T = S + A (m - ca)(m - cb)^T, with the facet
    // centroid m and central second moment S.
    const int dim = GetDimension();
    const Facet& F = _cells._facets[f];
    const double* ca = _GetTarget(a);
    const double* cb = _GetTarget(b);
    const double k = 2.0 * _GetHessianInverseDistance(f, minDistance);
    for (int r = 0; r < dim; ++r) {
        for (int c = 0; c < dim; ++c) {
            const double outer = (F._centroid[r] - ca[r]) * (F._centroid[c] - cb[c]);
            B[dim * r + c] = k * (F._secondMoment[dim * r + c] + F._area * outer);
        }
    }
}

void
ShiftedDiagram::GetShiftWeightBlock(size_t f, int a, double* b, double minDistance) const
{
    // ∫ (x - ca) = A (m - ca).
    const int dim = GetDimension();
    const Facet& F = _cells._facets[f];
    const double* ca = _GetTarget(a);
    for (int c = 0; c < dim; ++c) {
        b[c] = F._area * (F._centroid[c] - ca[c]) * _GetHessianInverseDistance(f, minDistance);
    }
}

double
ShiftedDiagram::GetWeightWeightCoef(size_t f, double minDistance) const
{
    return 0.5 * _cells._facets[f]._area * _GetHessianInverseDistance(f, minDistance);
}

} // namespace PP
