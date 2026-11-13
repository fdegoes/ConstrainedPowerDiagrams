#ifndef PP_CORE_SHIFTEDDIAGRAM_H
#define PP_CORE_SHIFTEDDIAGRAM_H

#include "core/powerCells.h"

#include <algorithm>
#include <string>
#include <vector>

namespace PP
{

class Domain;

// Power diagram of the shifted weighted sites (latex/contributions.tex,
// eq:shifted-site), clipped to a domain:
//   p_i      = q_i + u_i                              (generators)
//   lambda_i = w_i + |u_i|^2 + 2 u_i^T (q_i - cbar_i)  (shifted weights)
// for sites q, prescribed centroids cbar, shifts u and weights w. Its cells
// are the power cells of (p, lambda).
//
// Besides the cells, it gives the facet coefficients of the Hessian of the
// joint dual H (latex/appendix-joint.tex, eq:appendix-blocks-off). For a
// facet A_ij with generator distance d_ij = |p_i - p_j| and a, b in {i, j}:
//   GetShiftShiftBlock(f, a, b)  = 2 / d_ij ∫_{A_ij} (x - cbar_a)(x - cbar_b)^T dA
//   GetShiftWeightBlock(f, a)    = 1 / d_ij ∫_{A_ij} (x - cbar_a) dA
//   GetWeightWeightCoef(f)       = A_ij / (2 d_ij)
// so that, for i ~ j,
//   d2H/du_i du_j^T = GetShiftShiftBlock(f, i, j)
//   d2H/du_i dw_j   = GetShiftWeightBlock(f, i)
//   d2H/dw_i dw_j   = GetWeightWeightCoef(f)
// and the diagonal blocks are minus the sums of the (a, a) terms over the
// facets of cell i (eq:appendix-blocks-diag).
class ShiftedDiagram
{
public:
    using Facet = PowerCells::Facet;

    // sites    dim doubles per site (dim = domain.GetDimension())
    // targets  dim doubles per site, the prescribed centroids cbar
    // shifts   dim doubles per site
    // weights  1 double per site
    bool Build(
        const Domain& domain,
        const std::vector<double>& sites,
        const std::vector<double>& targets,
        const std::vector<double>& shifts,
        const std::vector<double>& weights,
        std::string& err);

    int GetDimension() const { return _cells._dim; }
    int GetNumSites() const { return int(_shiftedWeights.size()); }

    // Generators p_i and shifted weights lambda_i.
    const std::vector<double>& GetGenerators() const { return _generators; }
    const double* GetGenerator(int i) const { return &_generators[GetDimension() * i]; }
    const std::vector<double>& GetShiftedWeights() const { return _shiftedWeights; }
    double GetShiftedWeight(int i) const { return _shiftedWeights[i]; }

    // Cells: volumes V_i and centroids c_i.
    const PowerCells& GetCells() const { return _cells; }
    double GetVolume(int i) const { return _cells._volumes[i]; }
    const double* GetCentroid(int i) const { return &_cells._centroids[GetDimension() * i]; }
    int GetNumEmptyCells() const { return _cells.GetNumEmptyCells(); }
    // False when the cells could not be clipped (PowerCells::_complete).
    bool IsComplete() const { return _cells._complete; }
    // Smallest cell volume.
    double GetMinVolume() const;

    // Power-diagram energy E(p, lambda) = sum_i ∫_{V_i} |x - p_i|^2 - lambda_i dx
    // (latex/background.tex, eq:power-diagram).
    double GetEnergy() const;

    // Facets A_ij shared by two cells (i < j).
    size_t GetNumFacets() const { return _cells._facets.size(); }
    const Facet& GetFacet(size_t f) const { return _cells._facets[f]; }
    // |p_i - p_j|; facets with coincident generators have 0 and must be
    // skipped.
    double GetGeneratorDistance(size_t f) const { return _distances[f]; }

    // Hessian facet coefficients, see the class comment. a and b are
    // GetFacet(f)._i or GetFacet(f)._j. B gets dim x dim values
    // (row-major), b gets dim values. minDistance is a lower bound on the
    // generator distance in their 1 / d_ij factors: d_ij is replaced by
    // max(d_ij, minDistance), which bounds the curvature of close generator
    // pairs. 0 (default) keeps the exact Hessian.
    void GetShiftShiftBlock(size_t f, int a, int b, double* B, double minDistance = 0.0) const;
    void GetShiftWeightBlock(size_t f, int a, double* b, double minDistance = 0.0) const;
    double GetWeightWeightCoef(size_t f, double minDistance = 0.0) const;

private:
    // Prescribed centroid cbar_a.
    const double* _GetTarget(int a) const { return &_targets[GetDimension() * a]; }

    PowerCells          _cells;
    std::vector<double> _targets;
    std::vector<double> _generators;
    std::vector<double> _shiftedWeights;
    std::vector<double> _distances;  // per facet

    // 1 / d_ij of the Hessian coefficients (with the lower bound).
    double _GetHessianInverseDistance(size_t f, double minDistance) const
    {
        return 1.0 / std::max(_distances[f], minDistance);
    }
};

} // namespace PP

#endif // PP_CORE_SHIFTEDDIAGRAM_H
