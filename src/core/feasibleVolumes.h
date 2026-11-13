#ifndef PP_CORE_FEASIBLEVOLUMES_H
#define PP_CORE_FEASIBLEVOLUMES_H

#include <Eigen/Core>

#include <string>
#include <vector>

namespace PP
{

// Dual of the maximum-entropy feasibility problem
// (latex/contributions.tex, "Remarks / Feasibility"):
//   min_V sum_i V_i log(V_i / V_D)
//   s.t.  sum_i V_i cbar_i = V_D c_D,   sum_i V_i = V_D.
//
// Normalized to V_D = 1 and with the targets centered at the domain
// centroid and scaled to unit size (cbar_i -> (cbar_i - c_D) / s), the dual
// over y in R^{d+1} is
//   f(y)  = sum_i exp(-y . l_i) + y_d,       l_i = [cbar_i, 1]
//   df/dy = [0, 1] - sum_i f_i l_i,          f_i = exp(-y . l_i)
//   H     = sum_i f_i l_i l_i^T              (positive definite if the
//                                             targets span R^d)
// and the volume fractions are f_i.
class MaxEntropyEnergy
{
public:
    using Vector = Eigen::VectorXd;
    using Matrix = Eigen::MatrixXd;

    // targets   dim doubles per target, already centered and scaled
    MaxEntropyEnergy(const std::vector<double>& targets, int dim);

    int GetNumVariables() const { return _dim + 1; }
    // +infinity if an exponent overflows.
    double Eval(const Vector& y) const;
    void EvalGradient(const Vector& y, Vector& g) const;
    void EvalHessian(const Vector& y, Matrix& H) const;

    // Fractions f_i = exp(-y . l_i) (sum to 1 at the minimizer).
    void ComputeFractions(const Vector& y, std::vector<double>& fractions) const;

private:
    // Returns false if an exponent overflows.
    bool _ComputeExponents(const Vector& y, std::vector<double>& e) const;

    const std::vector<double>& _targets;
    int                        _dim;
    int                        _n;
};

// Computes volume fractions f_i > 0, sum_i f_i = 1, for which
// sum_i f_i cbar_i = c_D (the domain centroid), so that cells with volumes
// f_i |domain| and centroids cbar_i satisfy the partition conservation
// identities. Among all such fractions it returns the maximum-entropy one.
//
// A solution exists exactly when c_D lies in the interior of the convex
// hull of the targets.
class FeasibleVolumeSolver
{
public:
    using Vector = MaxEntropyEnergy::Vector;

    struct Parameters
    {
        int    _maxIterations = 100;
        // Converged when |sum_i f_i l_i - [c_D, 1]|_inf <= tolerance
        // (in units of the target bounding box).
        double _tolerance = 1e-10;
        // Newton steps solve (H + regularization I) dy = -g.
        double _regularization = 1e-12;
        // Max step halvings of the backtracking line search.
        int    _maxHalvings = 50;
    };

    struct Result
    {
        bool        _feasible = false;
        int         _iterations = 0;
        double      _residual = 0.0;  // |sum f_i l_i - [c_D, 1]|_inf
        std::string _error;
    };

    FeasibleVolumeSolver() = default;
    explicit FeasibleVolumeSolver(const Parameters& params) : _params(params) { }

    Parameters& GetParameters() { return _params; }
    const Parameters& GetParameters() const { return _params; }

    // targets        dim doubles per target centroid (e.g. the sites)
    // domainCentroid dim doubles
    // fractions      on return, one fraction per target
    Result Solve(
        const std::vector<double>& targets,
        const std::vector<double>& domainCentroid,
        int dim,
        std::vector<double>& fractions) const;

private:
    // Damped Newton from y; returns true if |g|_inf <= tolerance.
    bool _Minimize(const MaxEntropyEnergy& energy, Vector& y, Result& res) const;

    Parameters _params;
};

} // namespace PP

#endif // PP_CORE_FEASIBLEVOLUMES_H
