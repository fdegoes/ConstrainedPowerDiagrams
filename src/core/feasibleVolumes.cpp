#include "core/feasibleVolumes.h"

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>

namespace PP
{

namespace
{

// exp() overflows past this exponent.
constexpr double kMaxExponent = 700.0;

} // namespace

MaxEntropyEnergy::MaxEntropyEnergy(const std::vector<double>& targets, int dim)
    : _targets(targets), _dim(dim), _n(int(targets.size() / dim))
{ }

bool
MaxEntropyEnergy::_ComputeExponents(const Vector& y, std::vector<double>& e) const
{
    e.resize(_n);
    for (int i = 0; i < _n; ++i) {
        double s = y[_dim];
        for (int c = 0; c < _dim; ++c) {
            s += y[c] * _targets[_dim * i + c];
        }
        e[i] = -s;
        if (e[i] > kMaxExponent) {
            return false;
        }
    }
    return true;
}

double
MaxEntropyEnergy::Eval(const Vector& y) const
{
    std::vector<double> e;
    if (not _ComputeExponents(y, e)) {
        return HUGE_VAL;
    }
    double f = y[_dim];
    for (double x : e) {
        f += std::exp(x);
    }
    return f;
}

void
MaxEntropyEnergy::EvalGradient(const Vector& y, Vector& g) const
{
    std::vector<double> e;
    _ComputeExponents(y, e);
    g = Vector::Zero(_dim + 1);
    g[_dim] = 1.0;
    for (int i = 0; i < _n; ++i) {
        const double f = std::exp(std::min(e[i], kMaxExponent));
        for (int c = 0; c < _dim; ++c) {
            g[c] -= f * _targets[_dim * i + c];
        }
        g[_dim] -= f;
    }
}

void
MaxEntropyEnergy::EvalHessian(const Vector& y, Matrix& H) const
{
    std::vector<double> e;
    _ComputeExponents(y, e);
    const int m = _dim + 1;
    H = Matrix::Zero(m, m);
    Vector l(m);
    for (int i = 0; i < _n; ++i) {
        const double f = std::exp(std::min(e[i], kMaxExponent));
        for (int c = 0; c < _dim; ++c) {
            l[c] = _targets[_dim * i + c];
        }
        l[_dim] = 1.0;
        H.noalias() += f * l * l.transpose();
    }
}

void
MaxEntropyEnergy::ComputeFractions(const Vector& y, std::vector<double>& fractions) const
{
    std::vector<double> e;
    _ComputeExponents(y, e);
    fractions.resize(_n);
    for (int i = 0; i < _n; ++i) {
        fractions[i] = std::exp(std::min(e[i], kMaxExponent));
    }
}

namespace
{

// Targets centered at center and divided by their largest coordinate
// magnitude scale (0 if all targets are at the center).
std::vector<double>
_NormalizeTargets(
    const std::vector<double>& targets,
    const std::vector<double>& center,
    int dim,
    double& scale)
{
    std::vector<double> local(targets.size());
    scale = 0.0;
    for (size_t k = 0; k < targets.size(); ++k) {
        local[k] = targets[k] - center[k % dim];
        scale = std::max(scale, std::fabs(local[k]));
    }
    if (scale > 0.0) {
        for (double& x : local) {
            x /= scale;
        }
    }
    return local;
}

bool
_AreAllPositiveAndFinite(const std::vector<double>& values)
{
    return std::all_of(values.begin(), values.end(),
                       [](double v) { return v > 0.0 and std::isfinite(v); });
}

bool
_CheckFeasibilityInputs(
    const std::vector<double>& targets,
    const std::vector<double>& domainCentroid,
    int dim,
    std::string& err)
{
    err = "Invalid targets or domain centroid.";
    if (dim < 1) {
        return false;
    }
    if (targets.empty()) {
        return false;
    }
    if (targets.size() % dim != 0) {
        return false;
    }
    if (domainCentroid.size() != size_t(dim)) {
        return false;
    }
    err.clear();
    return true;
}

// The solve converged, the conservation residual is small, and the
// fractions are positive.
bool
_IsFeasible(bool converged, double residual, const std::vector<double>& fractions, double tolerance)
{
    if (not converged) {
        return false;
    }
    if (residual > 10.0 * tolerance) {
        return false;
    }
    return _AreAllPositiveAndFinite(fractions);
}

} // namespace

bool
FeasibleVolumeSolver::_Minimize(const MaxEntropyEnergy& energy, Vector& y, Result& res) const
{
    double f = energy.Eval(y);
    Vector g;
    MaxEntropyEnergy::Matrix H;
    for (res._iterations = 0; res._iterations < _params._maxIterations; ++res._iterations) {
        energy.EvalGradient(y, g);
        if (g.lpNorm<Eigen::Infinity>() <= _params._tolerance) {
            return true;
        }
        energy.EvalHessian(y, H);
        H.diagonal().array() += _params._regularization;
        const Eigen::LLT<MaxEntropyEnergy::Matrix> llt(H);
        if (llt.info() != Eigen::Success) {
            return false;
        }
        const Vector dy = -llt.solve(g);
        // Halve the step until the energy decreases.
        double t = 1.0;
        int k = 0;
        for (; k <= _params._maxHalvings; ++k, t *= 0.5) {
            const double fTrial = energy.Eval(y + t * dy);
            if (fTrial < f) {
                f = fTrial;
                break;
            }
        }
        if (k > _params._maxHalvings) {
            // No decrease: converged to round-off, or stuck.
            energy.EvalGradient(y, g);
            return g.lpNorm<Eigen::Infinity>() <= _params._tolerance;
        }
        y += t * dy;
    }
    energy.EvalGradient(y, g);
    return g.lpNorm<Eigen::Infinity>() <= _params._tolerance;
}

FeasibleVolumeSolver::Result
FeasibleVolumeSolver::Solve(
    const std::vector<double>& targets,
    const std::vector<double>& domainCentroid,
    int dim,
    std::vector<double>& fractions) const
{
    Result res;
    fractions.clear();
    if (not _CheckFeasibilityInputs(targets, domainCentroid, dim, res._error)) {
        return res;
    }
    const int n = int(targets.size() / dim);

    // Centered and scaled to unit size, for conditioning.
    double scale = 0.0;
    const std::vector<double> local = _NormalizeTargets(targets, domainCentroid, dim, scale);
    if (not (scale > 0.0)) {
        // All targets at the centroid: any fractions work, take uniform.
        fractions.assign(n, 1.0 / n);
        res._feasible = true;
        return res;
    }

    MaxEntropyEnergy energy(local, dim);
    // y = [0, log n] gives the uniform fractions 1/n.
    Vector y = Vector::Zero(dim + 1);
    y[dim] = std::log(double(n));
    const bool converged = _Minimize(energy, y, res);

    Vector g;
    energy.EvalGradient(y, g);
    res._residual = g.lpNorm<Eigen::Infinity>();
    energy.ComputeFractions(y, fractions);
    res._feasible = _IsFeasible(converged, res._residual, fractions, _params._tolerance);
    if (not res._feasible) {
        res._error = "No feasible volumes: the domain centroid is not inside the "
                    "convex hull of the sites.";
    }
    return res;
}

} // namespace PP
