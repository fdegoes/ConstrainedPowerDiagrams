#include "core/gmIntegrator.h"

namespace PP
{

bool
GMIntegrator::_Solve(const ParticleState& cur, std::vector<double>& weights, Result& res) const
{
    std::vector<double> p;
    _GetDiagram(cur, p, weights);
    VolumeSolver solver(_params._solver);
    res._solve = solver.Solve(_domain, cur._positions, _volumes, weights);
    res._status = DualSolver::GetStatusString(res._solve._status);
    res._iterations = res._solve._iterations;
    if (not res._solve.HasConverged()) {
        res._error = res._solve._error;
        return false;
    }
    return true;
}

Integrator::Result
GMIntegrator::Initialize(ParticleState& state) const
{
    Result res;
    ParticleState next;
    std::vector<double> weights;
    if (not _Prepare(state, next, res) or not _Solve(next, weights, res)) {
        return res;
    }
    next._shifts.assign(next._positions.size(), 0.0);
    next._weights = weights;
    _Finish(res._solve._diagram, next, res);
    res._converged = true;
    state = std::move(next);
    return res;
}

bool
GMIntegrator::Measure(const ParticleState& state, FlowMeasures& measures, std::string& err) const
{
    Result res;
    ParticleState s;
    std::vector<double> weights;
    if (not _Prepare(state, s, res) or not _Solve(s, weights, res)) {
        err = res._error;
        return false;
    }
    measures = _Measure(res._solve._diagram, s);
    return true;
}

Integrator::Result
GMIntegrator::Step(ParticleState& state, double dt) const
{
    Result res;
    ParticleState cur;
    if (not _Prepare(state, cur, res)) {
        return res;
    }
    const std::vector<double> startVelocities = cur._velocities;
    _ApplyGravity(cur._velocities, dt);
    std::vector<double> weights;
    if (not _Solve(cur, weights, res)) {
        return res;
    }
    {
        // The state at the start of the step on the diagram solved for it.
        ParticleState start = cur;
        start._velocities = startVelocities;
        res._startMeasures = _Measure(res._solve._diagram, start);
    }

    const std::vector<double>& q = cur._positions;
    const std::vector<double>& c = res._solve._diagram.GetCells()._centroids;
    const size_t m = q.size();
    const double force = _options._stiffness * dt;
    ParticleState next;
    next._positions.resize(m);
    next._velocities.resize(m);
    next._pressures = cur._pressures;
    for (size_t k = 0; k < m; ++k) {
        next._velocities[k] = cur._velocities[k] + force * (c[k] - q[k]);
        next._positions[k] = q[k] + dt * next._velocities[k];
    }
    _SetDiagram(next, q, weights);
    _Finish(res._solve._diagram, next, res);
    res._converged = true;
    state = std::move(next);
    return res;
}

} // namespace PP
