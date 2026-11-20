#ifndef PP_CORE_JOINTDUAL_H
#define PP_CORE_JOINTDUAL_H

#include "core/shiftedDiagram.h"

#include <Eigen/Core>
#include <Eigen/SparseCore>

#include <cmath>
#include <string>
#include <vector>

namespace PP
{

class Domain;

// Concave dual objectives of the constrained power diagrams
// (latex/contributions.tex), as functions of the stacked multipliers z:
//
//  - Joint:   H_alpha(u, w) = H(u, w) - 1/alpha sum_i Vbar_i |u_i|^2
//             (eq:regularized-dual), with H = E(p, lambda) + sum_i w_i Vbar_i
//             (eq:joint-dual). alpha = infinity gives H itself.
//             z_i = (u_i, w_i) in R^{d+1}, stored per site:
//             z[(d+1) i + c] = u_i[c], z[(d+1) i + d] = w_i.
//  - Shifts:  G(u) = E(p, lambda) (eq:centroid-dual), i.e. H with the
//             weights fixed (SetFixedWeights()) and without the volume
//             terms. z = u, z[d i + c] = u_i[c].
//  - Weights: F(w) = E(q, w) + sum_i w_i Vbar_i (latex/background.tex,
//             eq:volume-dual), i.e. H with the shifts fixed at 0.
//             z = w.
//
// The power diagram at z is a ShiftedDiagram (Build()); the value,
// gradient and Hessian at z are read from it.
class JointDual
{
public:
    using Vector = Eigen::VectorXd;
    using SparseMatrix = Eigen::SparseMatrix<double>;

    enum class Unknowns
    {
        ShiftsAndWeights,  // H_alpha (or H)
        Shifts,            // G
        Weights,           // F
    };

    // All arrays must outlive the objective.
    // sites    dim doubles per site, q
    // targets  dim doubles per site, cbar (ignored for Unknowns::Weights)
    // volumes  one positive target volume per site, Vbar. G has no volume
    //          constraints: they only scale its regularization eta diag(Vbar)
    //          and its residual (e.g. the volumes of the expected cells).
    // alpha    stiffness of the softened centroids; HUGE_VAL for H
    JointDual(
        const Domain& domain,
        const std::vector<double>& sites,
        const std::vector<double>& targets,
        const std::vector<double>& volumes,
        double alpha,
        Unknowns unknowns);

    const std::vector<double>& GetSites() const { return _sites; }
    const std::vector<double>& GetTargets() const { return _targets; }
    const std::vector<double>& GetVolumes() const { return _volumes; }
    double GetAlpha() const { return _alpha; }
    bool HasShifts() const { return _unknowns != Unknowns::Weights; }
    bool HasWeights() const { return _unknowns != Unknowns::Shifts; }

    // Weights of G (Unknowns::Shifts), one per site. Must be set before
    // Build().
    void SetFixedWeights(const std::vector<double>& weights) { _fixedWeights = weights; }

    int GetDimension() const { return _dim; }
    int GetNumSites() const { return int(_volumes.size()); }
    int GetNumVariables() const { return GetStride() * GetNumSites(); }
    // Number of unknowns per site (d + 1, d or 1).
    int GetStride() const { return (HasShifts() ? _dim : 0) + (HasWeights() ? 1 : 0); }

    // Shifts u (dim per site, zeros for Unknowns::Weights) and weights w
    // (the fixed weights for Unknowns::Shifts) stored in z, and back.
    void Unpack(const Vector& z, std::vector<double>& shifts, std::vector<double>& weights) const;
    Vector Pack(const std::vector<double>& shifts, const std::vector<double>& weights) const;

    // Power diagram of the shifted sites at z.
    bool Build(const Vector& z, ShiftedDiagram& diagram, std::string& err) const;

    // Value and gradient (eq:joint-gradient, eq:regularized-gradient-shift)
    // at z, whose diagram is diagram.
    double GetValue(const Vector& z, const ShiftedDiagram& diagram) const;
    void GetGradient(const Vector& z, const ShiftedDiagram& diagram, Vector& g) const;

    // Diagonal of the regularization mass M, scaled to the Hessian blocks:
    // with r_i = Vbar_i^(1/d), r_i^(d+2) = Vbar_i r_i^2 on the shifts of
    // site i (the shift-shift block (2/d_ij) ∫ (x - cbar)(x - cbar)^T scales
    // like L^d) and Vbar_i on its weight (the weight-weight block
    // A_ij / (2 d_ij) scales like L^(d-2)), so that eta has units L^-2 for
    // both.
    Vector GetMassDiagonal() const;

    // Tolerances of the termination test (DualSolver::Parameters): on the
    // relative volume error, and on the shift gradient in units of
    // Vbar_i^((d+1)/d).
    void SetTolerances(double volume, double shift)
    {
        _volumeTolerance = volume;
        _shiftTolerance = shift;
    }

    // Residuals of a gradient-like vector g (the gradient, or a
    // stationarity condition with the same units), the largest per-cell
    // scaled entries
    //   volume = max_i |g_w_i| / Vbar_i,
    //   shift  = max_i |g_u_i| / (2 Vbar_i^((d+1)/d)).
    // For the gradient, these are the relative volume error
    // |V_i - Vbar_i| / Vbar_i and V_i / Vbar_i |cbar_i - c_i| / r_i, with
    // r_i = Vbar_i^(1/d), of every cell: the centroid error in cell sizes
    // (g_u_i = 2 V_i (cbar_i - c_i) carries the factor 2). Zero for absent unknowns.
    void GetResiduals(const Vector& g, double& volume, double& shift) const;

    // max(volume / volumeTolerance, shift / shiftTolerance): at most 1 when
    // both tests pass.
    double GetResidual(const Vector& g) const;

    // N = -Hessian + eta M (latex/implementation.tex, the diagonal
    // regularization; latex/appendix-algorithms.tex, Algorithm hessian-reg,
    // with the mass M of GetMassDiagonal()), from the facet
    // blocks of latex/appendix-joint.tex, with the distance clamp and the
    // pair regularization below. Both triangles are filled.
    void GetNegatedHessian(const ShiftedDiagram& diagram, double eta, SparseMatrix& N) const;

    // Lower bound of the generator distances in the Hessian facet
    // coefficients (1 / d_ij), as a ratio of the cell size (|domain| / n)^(1/d): bounds the
    // curvature of close generator pairs. The gradient is unchanged, so the
    // solution is too. 0 keeps the exact Hessian (the default of JointDual;
    // DualSolver sets it from its parameters).
    void SetHessianMinDistanceRatio(double ratio) { _hessianMinDistanceRatio = ratio; }

    // Pair regularization of N: for every facet pair (i, j) with generator
    // distance d = d_ij < R, adds rho_ij(d) e e^T to the shift-shift
    // diagonal blocks of i and j, with the unit direction e = (p_j - p_i) / d,
    // so that steps changing the separation of close generator pairs are
    // damped, with
    //   rho_ij(d) = rhoMax_ij (1 - d / R)^4    d < R,   0 beyond,
    // which is smooth at R and steep near d = 0. R is radiusRatio times the
    // cell size (|domain| / n)^(1/d), scaled by the activity
    // (SetPairActivity()). The maximum is per pair,
    //   rhoMax_ij = scale max(e^T B_i e, e^T B_j e),
    // with B_i the shift-shift diagonal block of site i of the Hessian: the
    // curvature of each site along the pair, so that the term matches the
    // stiffness it competes with. Only N changes, not the gradient, so the
    // solution is the same. radiusRatio = 0 (the default of JointDual)
    // disables it; DualSolver sets it from its parameters.
    void SetPairRegularization(double radiusRatio, double scale)
    {
        _pairRadiusRatio = radiusRatio;
        _pairScale = scale;
    }
    // Activity in [0, 1] scaling the radii of the pair regularization
    // (DualSolver drives it with the shift residual, so that it vanishes at
    // convergence).
    void SetPairActivity(double activity) { _pairActivity = activity; }

    // Frozen dual edges: the facets f = (i, j) of a diagram with distinct
    // generators, with their unit direction e = (p_j - p_i) / d and
    // distance d at the freeze. At z, the stretch of an edge is
    //   s(z) = (p_j(z) - p_i(z)) . e / d,    p = q + u,
    // linear in the shifts and 1 at the freeze. DualSolver freezes the
    // edges of the current diagram once per iteration (FreezeEdges()) to
    // limit how much its steps shrink them (Parameters::_stretchLimit).
    void FreezeEdges(const ShiftedDiagram& diagram);
    void ClearFrozenEdges() { _frozen.clear(); }
    bool HasFrozenEdges() const { return not _frozen.empty(); }
    // Smallest stretch s_f(z) of the frozen edges (infinity if none).
    double GetMinStretch(const Vector& z) const;
    // Largest t in (0, 1] such that every frozen edge keeps
    // s_f(z + t dz) >= limit s_f(z) (the stretches are linear in t; at the
    // freeze, s_f(z) = 1).
    double GetMaxStep(const Vector& z, const Vector& dz, double limit) const;

private:
    int _GetWeightIndex(int i) const { return GetStride() * i + GetStride() - 1; }

    // N in assembly: its compressed pattern (the diagonal blocks and one
    // block per neighbor pair) and the positions of the blocks (jointDual.cpp).
    struct _Assembly;

    // Adds the facet blocks of the shift and weight unknowns of facet f to
    // N, with the generator distances clamped below by minDistance.
    void _AddFacetBlocks(
        const ShiftedDiagram& diagram,
        size_t f,
        double minDistance,
        _Assembly& N) const;

    const Domain&              _domain;
    const std::vector<double>& _sites;
    const std::vector<double>& _targets;
    const std::vector<double>& _volumes;
    double                     _alpha;
    Unknowns                   _unknowns;
    int                        _dim;
    double                     _volumeTolerance = 1e-8;
    double                     _shiftTolerance = 1e-8;
    std::vector<double>        _fixedWeights;
    double                     _hessianMinDistanceRatio = 0.0;
    double                     _pairRadiusRatio = 0.0;
    double                     _pairScale = 1.0;
    double                     _pairActivity = 1.0;

    // A frozen dual edge (FreezeEdges()).
    struct _FrozenEdge
    {
        int    _i, _j;
        double _e[3];
        double _d;
    };
    std::vector<_FrozenEdge> _frozen;
    double _GetStretch(const _FrozenEdge& E, const Vector& z) const;

    double _GetCellSize() const;
    // Adds the pair regularization (SetPairRegularization()) to N, which
    // holds the facet blocks.
    void _AddPairRegularization(const ShiftedDiagram& diagram, _Assembly& N) const;
};

} // namespace PP

#endif // PP_CORE_JOINTDUAL_H
