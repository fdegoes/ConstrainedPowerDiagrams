#include "core/constrainedIntegrator.h"
#include "core/domain.h"

namespace PP
{

Integrator::Result
ConstrainedIntegrator::Initialize(ParticleState& state) const
{
    Result res;
    ParticleState next;
    if (not _Prepare(state, next, res) or not _MakeCentroidal(next, res)) {
        return res;
    }
    ShiftedDiagram diagram;
    if (not _BuildDiagram(next, diagram, res._error)) {
        return res;
    }
    _Finish(diagram, next, res);
    res._converged = true;
    state = std::move(next);
    return res;
}

Integrator::Result
ConstrainedIntegrator::Step(ParticleState& state, double dt) const
{
    return StepWithGauge(state, dt);
}

ConstrainedIntegrator::StepResult
ConstrainedIntegrator::StepWithGauge(ParticleState& state, double dt) const
{
    StepResult res;
    ParticleState cur;
    if (not _Prepare(state, cur, res)) {
        return res;
    }
    _ApplyGravity(cur._velocities, dt);
    const int dim = _GetDimension();
    const std::vector<double>& q = cur._positions;
    const size_t m = q.size();
    std::vector<double> generators, lambda;
    _GetDiagram(cur, generators, lambda);

    std::vector<double> predicted(m);
    for (size_t k = 0; k < m; ++k) {
        predicted[k] = (1.0 + dt * dt) * q[k] + dt * cur._velocities[k];
    }

    std::vector<double> shifts;
    std::vector<double> weights = lambda;
    res._solve = MakeJointStart(_domain, q, predicted, _volumes, 0.0, shifts, weights);
    res._status = DualSolver::GetStatusString(res._solve._status);
    if (not res._solve.HasConverged()) {
        res._error = res._solve._error;
        return res;
    }
    ConstrainedJointSolver solver(_params._solver);
    res._solve = solver.Solve(_domain, q, predicted, _volumes, generators, shifts, weights, res._gauge);
    res._status = DualSolver::GetStatusString(res._solve._status);
    res._iterations = res._solve._iterations;
    if (not res._solve.HasConverged()) {
        res._error = res._solve._error;
        return res;
    }

    const double s = res._gauge._dilation;
    const std::vector<double>& t = res._gauge._translation;
    ParticleState next;
    next._positions.resize(m);
    next._velocities.resize(m);
    next._pressures = cur._pressures;
    std::vector<double> nextGenerators(m), nextLambda(weights.size());
    for (size_t i = 0; i < weights.size(); ++i) {
        double l = weights[i];
        for (int c = 0; c < dim; ++c) {
            const size_t k = size_t(dim) * i + size_t(c);
            next._positions[k] = predicted[k] - s * generators[k] - t[c];
            next._velocities[k] = (next._positions[k] - q[k]) / dt;
            nextGenerators[k] = q[k] + shifts[k];
            l += shifts[k] * shifts[k] + 2.0 * shifts[k] * (q[k] - predicted[k]);
        }
        nextLambda[i] = l;
    }
    _SetDiagram(next, nextGenerators, nextLambda);
    ShiftedDiagram diagram;
    if (not _BuildDiagram(next, diagram, res._error)) {
        return res;
    }
    _Finish(diagram, next, res);
    res._converged = true;
    state = std::move(next);
    return res;
}

} // namespace PP
