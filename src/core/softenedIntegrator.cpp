#include "core/softenedIntegrator.h"

namespace PP
{

bool
SoftenedIntegrator::_Solve(
    const ParticleState& cur,
    std::vector<double>& shifts,
    std::vector<double>& weights,
    Result& res) const
{
    const std::vector<double>& q = cur._positions;
    std::vector<double> p;
    _GetDiagram(cur, p, weights);
    res._solve = MakeJointStart(_domain, q, q, _volumes, 0.0, shifts, weights);
    res._status = DualSolver::GetStatusString(res._solve._status);
    if (not res._solve.HasConverged()) {
        res._error = res._solve._error;
        return false;
    }
    SoftenedSolver solver(_params._solver);
    res._solve = solver.Solve(_domain, q, q, _volumes, _options._alpha, shifts, weights);
    res._status = DualSolver::GetStatusString(res._solve._status);
    res._iterations = res._solve._iterations;
    if (not res._solve.HasConverged()) {
        res._error = res._solve._error;
        return false;
    }
    return true;
}

Integrator::Result
SoftenedIntegrator::Initialize(ParticleState& state) const
{
    Result res;
    ParticleState next;
    std::vector<double> shifts, weights;
    if (not _Prepare(state, next, res) or not _Solve(next, shifts, weights, res)) {
        return res;
    }
    next._shifts = shifts;
    next._weights = weights;
    _Finish(res._solve._diagram, next, res);
    res._converged = true;
    state = std::move(next);
    return res;
}

Integrator::Result
SoftenedIntegrator::Step(ParticleState& state, double dt) const
{
    Result res;
    ParticleState cur;
    if (not _Prepare(state, cur, res)) {
        return res;
    }
    _ApplyGravity(cur._velocities, dt);
    std::vector<double> shifts, weights;
    if (not _Solve(cur, shifts, weights, res)) {
        return res;
    }

    const int dim = _GetDimension();
    const std::vector<double>& q = cur._positions;
    const size_t m = q.size();
    const double force = _options._stiffness * (1.0 + 1.0 / _options._alpha) * dt;
    ParticleState next;
    next._positions.resize(m);
    next._velocities.resize(m);
    next._pressures = cur._pressures;
    std::vector<double> generators(m), lambda(weights.size());
    for (size_t i = 0; i < weights.size(); ++i) {
        double l = weights[i];
        for (int c = 0; c < dim; ++c) {
            const size_t k = size_t(dim) * i + size_t(c);
            next._velocities[k] = cur._velocities[k] - force * shifts[k];
            next._positions[k] = q[k] + dt * next._velocities[k];
            generators[k] = q[k] + shifts[k];
            l += shifts[k] * shifts[k];
        }
        lambda[i] = l;
    }
    _SetDiagram(next, generators, lambda);
    _Finish(res._solve._diagram, next, res);
    res._converged = true;
    state = std::move(next);
    return res;
}

} // namespace PP
