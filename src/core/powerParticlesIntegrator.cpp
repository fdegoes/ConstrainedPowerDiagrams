#include "core/powerParticlesIntegrator.h"
#include "core/domain.h"
#include "core/geometry.h"
#include "core/shiftedDiagram.h"
#include "core/sparseCholesky.h"

#include <Eigen/SparseCore>

#include <algorithm>
#include <cmath>

namespace PP
{

namespace
{

using Triplets = std::vector<Eigen::Triplet<double>>;

// Facets of the 2015 code's operators: shared by two cells of positive
// measure, with distinct sites.
bool
_IsOperatorFacet(const PowerCells::Facet& F, double length)
{
    return F._area > 0.0 and length > 0.0;
}

double
_Length(const std::vector<double>& q, int i, int j, int dim)
{
    return std::sqrt(ComputeDistance2(&q[size_t(dim) * i], &q[size_t(dim) * j], dim));
}

// (x - b) . v
double
_Dot(const double* x, const double* b, const double* v, int dim)
{
    double s = 0.0;
    for (int c = 0; c < dim; ++c) {
        s += (x[c] - b[c]) * v[c];
    }
    return s;
}

// Radius of the ball of volume V.
double
_BallRadius(double V, int dim)
{
    return dim == 2 ? std::sqrt(V / M_PI) : std::cbrt(3.0 * V / (4.0 * M_PI));
}

// Solves the SPD system A x = b (lower triangle read); false if A is not
// positive definite.
bool
_SolveSPD(int n, const Triplets& entries, const Eigen::VectorXd& b, Eigen::VectorXd& x)
{
    Eigen::SparseMatrix<double> A(n, n);
    A.setFromTriplets(entries.begin(), entries.end());
    SparseCholesky llt;
    llt.compute(A);
    if (llt.info() != Eigen::Success) {
        return false;
    }
    x = llt.solve(b);
    return x.allFinite();
}

} // namespace

PowerParticlesIntegrator::PowerParticlesIntegrator(
    const Domain& domain,
    std::vector<double> volumes,
    const Parameters& params,
    const IntegratorOptions& options)
    : Integrator(domain, std::move(volumes), params, options), _boundary(domain)
{ }

double
PowerParticlesIntegrator::_GetRegularization() const
{
    if (_options._regularization >= 0.0) {
        return _options._regularization;
    }
    return _GetDimension() == 2 ? 1e-8 : 1e-6;
}

bool
PowerParticlesIntegrator::_Build(
    const std::vector<double>& q,
    const std::vector<double>& w,
    ShiftedDiagram& diagram,
    std::string& err) const
{
    const std::vector<double> zeros(q.size(), 0.0);
    return diagram.Build(_domain, q, q, zeros, w, err);
}

void
PowerParticlesIntegrator::ComputeDivergence(
    const ShiftedDiagram& diagram,
    const std::vector<double>& q,
    const std::vector<double>& velocities,
    std::vector<double>& divergence) const
{
    const int dim = _GetDimension();
    divergence.assign(_volumes.size(), 0.0);
    for (size_t f = 0; f < diagram.GetNumFacets(); ++f) {
        const PowerCells::Facet& F = diagram.GetFacet(f);
        const double l = _Length(q, F._i, F._j, dim);
        if (not _IsOperatorFacet(F, l)) {
            continue;
        }
        const double star = F._area / l;
        const double* b = F._centroid.data();
        const double* qi = &q[size_t(dim) * F._i];
        const double* qj = &q[size_t(dim) * F._j];
        const double* vi = &velocities[size_t(dim) * F._i];
        const double* vj = &velocities[size_t(dim) * F._j];
        // dV_i/dq_j . v_j - dV_j/dq_i . v_i for cell i, and the converse for j.
        const double t = star * (_Dot(qj, b, vj, dim) - _Dot(qi, b, vi, dim));
        divergence[F._i] += t;
        divergence[F._j] -= t;
    }
}

void
PowerParticlesIntegrator::ComputeGradient(
    const ShiftedDiagram& diagram,
    const std::vector<double>& q,
    const std::vector<double>& pressures,
    std::vector<double>& gradient) const
{
    const int dim = _GetDimension();
    gradient.assign(q.size(), 0.0);
    for (size_t f = 0; f < diagram.GetNumFacets(); ++f) {
        const PowerCells::Facet& F = diagram.GetFacet(f);
        const double l = _Length(q, F._i, F._j, dim);
        if (not _IsOperatorFacet(F, l)) {
            continue;
        }
        const double star = F._area / l;
        const double* b = F._centroid.data();
        const double dp = pressures[F._i] - pressures[F._j];
        for (int c = 0; c < dim; ++c) {
            gradient[size_t(dim) * F._i + c] += star * (q[size_t(dim) * F._i + c] - b[c]) * dp;
            gradient[size_t(dim) * F._j + c] += star * (q[size_t(dim) * F._j + c] - b[c]) * (-dp);
        }
    }
}

bool
PowerParticlesIntegrator::_ApplyPressure(const ShiftedDiagram& diagram, ParticleState& state, double dt,
                                         Result& res) const
{
    const int dim = _GetDimension();
    const int n = int(_volumes.size());
    const std::vector<double>& q = state._positions;
    const std::vector<double>& w = state._weights;
    const double eps = _GetRegularization();

    std::vector<double> divergence;
    ComputeDivergence(diagram, q, state._velocities, divergence);

    // rho_i = mass / volume, mass = Vbar (unit rest density).
    std::vector<double> rho(static_cast<size_t>(n));
    for (int i = 0; i < n; ++i) {
        const double V = diagram.GetVolume(i);
        rho[i] = V > 0.0 ? _volumes[i] / V : 1.0;
    }
    std::vector<double> diag(static_cast<size_t>(n));
    for (int i = 0; i < n; ++i) {
        diag[i] = eps * diagram.GetVolume(i);
    }
    Triplets L;
    for (size_t f = 0; f < diagram.GetNumFacets(); ++f) {
        const PowerCells::Facet& F = diagram.GetFacet(f);
        const double l = _Length(q, F._i, F._j, dim);
        if (not _IsOperatorFacet(F, l)) {
            continue;
        }
        const double dij = 0.5 * (l + (w[F._i] - w[F._j]) / l);
        const double dji = 0.5 * (l + (w[F._j] - w[F._i]) / l);
        double ratio = (dij + dji) / (rho[F._i] * dij + rho[F._j] * dji);
        if (ratio < 0.0) {
            ratio = 2.0 / (rho[F._i] + rho[F._j]);
        }
        const double val = dt * (F._area / l) * ratio;
        diag[F._i] += val;
        diag[F._j] += val;
        L.emplace_back(std::max(F._i, F._j), std::min(F._i, F._j), -val);
    }
    Eigen::VectorXd b(n);
    for (int i = 0; i < n; ++i) {
        L.emplace_back(i, i, diag[i]);
        b[i] = -divergence[i];
    }
    Eigen::VectorXd p;
    if (not _SolveSPD(n, L, b, p)) {
        res._error = "The pressure Laplacian is not positive definite.";
        return false;
    }
    state._pressures.assign(p.data(), p.data() + n);

    std::vector<double> gradient;
    ComputeGradient(diagram, q, state._pressures, gradient);
    for (int i = 0; i < n; ++i) {
        for (int c = 0; c < dim; ++c) {
            state._velocities[size_t(dim) * i + c] -= dt * gradient[size_t(dim) * i + c] / _volumes[i];
        }
    }
    return true;
}

void
PowerParticlesIntegrator::_Advect(const ShiftedDiagram& diagram, ParticleState& state, double dt, Result& res) const
{
    const int dim = _GetDimension();
    const double beta = _options._centroidBlend;
    for (size_t i = 0; i < _volumes.size(); ++i) {
        double* q = &state._positions[size_t(dim) * i];
        const double* c = diagram.GetCentroid(int(i));
        const double* v = &state._velocities[size_t(dim) * i];
        double prev[3], next[3];
        for (int k = 0; k < dim; ++k) {
            prev[k] = q[k] + beta * (c[k] - q[k]);
            next[k] = prev[k] + dt * v[k];
        }
        if (_boundary.Evaluate(next) >= 0.0) {
            double proj[3], dir[3];
            _boundary.Project(next, proj);
            double len = 0.0;
            for (int k = 0; k < dim; ++k) {
                dir[k] = proj[k] - next[k];
                len += dir[k] * dir[k];
            }
            len = std::sqrt(len);
            const double r = _BallRadius(_volumes[i], dim);
            for (int k = 0; k < dim; ++k) {
                next[k] = proj[k] + (len > 0.0 ? r * dir[k] / len : 0.0);
            }
            ++res._numProjected;
        }
        const double* target = _boundary.Evaluate(next) < 0.0 ? next : prev;
        std::copy(target, target + dim, q);
    }
}

bool
PowerParticlesIntegrator::_SolveWeights(ParticleState& state, ShiftedDiagram& diagram, Result& res) const
{
    const int dim = _GetDimension();
    const int n = int(_volumes.size());
    const double eps = _GetRegularization();
    const std::vector<double>& q = state._positions;
    std::vector<double>& w = state._weights;
    w.assign(size_t(n), 0.0);
    if (not _Build(q, w, diagram, res._error)) {
        return false;
    }
    auto valid = [&](const ShiftedDiagram& d) { return d.IsComplete() and d.GetNumEmptyCells() == 0; };
    if (not valid(diagram)) {
        res._error = "The power diagram of the advected particles has empty cells.";
        return false;
    }
    // Tolerance of the 2015 code: cap_tolerance N (volume of the largest
    // ball)^2, the ball of the largest particle radius.
    const double vmax = *std::max_element(_volumes.begin(), _volumes.end());
    const double tol = _options._weightTolerance * n * vmax * vmax;
    for (int iter = 0; iter < _options._weightIterations; ++iter) {
        Eigen::VectorXd rhs(n);
        for (int i = 0; i < n; ++i) {
            rhs[i] = _volumes[i] - diagram.GetVolume(i);
        }
        if (rhs.squaredNorm() < tol) {
            break;
        }
        std::vector<double> diag(static_cast<size_t>(n));
        for (int i = 0; i < n; ++i) {
            diag[i] = eps * diagram.GetVolume(i);
        }
        Triplets H;
        for (size_t f = 0; f < diagram.GetNumFacets(); ++f) {
            const PowerCells::Facet& F = diagram.GetFacet(f);
            const double l = _Length(q, F._i, F._j, dim);
            if (not _IsOperatorFacet(F, l)) {
                continue;
            }
            const double val = 0.5 * F._area / l;
            diag[F._i] += val;
            diag[F._j] += val;
            H.emplace_back(std::max(F._i, F._j), std::min(F._i, F._j), -val);
        }
        for (int i = 0; i < n; ++i) {
            H.emplace_back(i, i, diag[i]);
        }
        Eigen::VectorXd x;
        if (not _SolveSPD(n, H, rhs, x)) {
            res._error = "The weight Hessian is not positive definite.";
            return false;
        }
        ++res._iterations;
        // w += x, halving (w -= x / 2) up to 3 times while a cell is empty.
        for (int i = 0; i < n; ++i) {
            w[i] += x[i];
        }
        ShiftedDiagram trial;
        bool ok = _Build(q, w, trial, res._error) and valid(trial);
        for (int attempt = 0; attempt < 3 and not ok; ++attempt) {
            x *= 0.5;
            for (int i = 0; i < n; ++i) {
                w[i] -= x[i];
            }
            ok = _Build(q, w, trial, res._error) and valid(trial);
        }
        if (not ok) {
            res._error = "A Newton step of the weights empties a cell.";
            return false;
        }
        diagram = std::move(trial);
    }
    return true;
}

Integrator::Result
PowerParticlesIntegrator::Initialize(ParticleState& state) const
{
    Result res;
    ParticleState next;
    if (not _Prepare(state, next, res) or not _MakeCentroidal(next, res)) {
        return res;
    }
    ShiftedDiagram diagram;
    if (not _Build(next._positions, next._weights, diagram, res._error)) {
        return res;
    }
    _Finish(diagram, next, res);
    res._converged = true;
    state = std::move(next);
    return res;
}

Integrator::Result
PowerParticlesIntegrator::Step(ParticleState& state, double dt) const
{
    Result res;
    ParticleState next;
    if (not _Prepare(state, next, res)) {
        return res;
    }
    // The diagram of the step: the sites with their weights, no shifts (the
    // power weights w + |u|^2 of the state).
    std::vector<double> p, lambda;
    _GetDiagram(next, p, lambda);
    next._shifts.assign(next._positions.size(), 0.0);
    next._weights = lambda;
    ShiftedDiagram diagram;
    if (not _Build(next._positions, next._weights, diagram, res._error)) {
        return res;
    }
    if (not diagram.IsComplete() or diagram.GetNumEmptyCells() > 0) {
        res._error = "The power diagram of the particles has empty cells.";
        return res;
    }
    _ApplyGravity(next._velocities, dt);
    if (not _ApplyPressure(diagram, next, dt, res)) {
        return res;
    }
    _Advect(diagram, next, dt, res);
    if (not _SolveWeights(next, diagram, res)) {
        return res;
    }
    _Finish(diagram, next, res);
    // Converged when the volumes are within the 2015 tolerance (or, with
    // tolerance 0, when the fixed number of Newton steps was done).
    res._converged = true;
    res._status = "Converged";
    state = std::move(next);
    return res;
}

} // namespace PP
