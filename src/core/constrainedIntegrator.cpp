#include "core/constrainedIntegrator.h"
#include "core/domain.h"

#include <Eigen/Dense>

namespace PP
{

std::vector<CutCellSide>
ConstrainedIntegrator::_FindCutCellSides(const ShiftedDiagram& diagram)
{
    // Facets of round-off measure (e.g. the zero-length facets that split a
    // vertex shared by four cells) are not dual edges.
    double total = 0.0;
    for (size_t f = 0; f < diagram.GetNumFacets(); ++f) {
        total += diagram.GetFacet(f)._area;
    }
    const double minArea = 1e-9 * total / double(std::max<size_t>(1, diagram.GetNumFacets()));
    std::vector<std::pair<int, int>> edges;
    edges.reserve(diagram.GetNumFacets());
    for (size_t f = 0; f < diagram.GetNumFacets(); ++f) {
        const auto& F = diagram.GetFacet(f);
        if (F._area > minArea and diagram.GetGeneratorDistance(f) > 0.0) {
            edges.emplace_back(F._i, F._j);
        }
    }
    return FindCutCellSides(int(diagram.GetNumSites()), edges);
}

bool
ConstrainedIntegrator::_ProjectVelocities(const ShiftedDiagram& diagram, ParticleState& state) const
{
    const int dim = _GetDimension();
    const int n = int(_volumes.size());
    const std::vector<CutCellSide> sides = _FindCutCellSides(diagram);
    // The normals: rows of dim n values, the dilation, the translations, then
    // one localized homothety per side.
    const int m = 1 + dim + int(sides.size());
    std::vector<std::vector<double>> rows(size_t(m), std::vector<double>(size_t(dim * n), 0.0));
    for (int i = 0; i < n; ++i) {
        const double* p = diagram.GetGenerator(i);
        for (int c = 0; c < dim; ++c) {
            rows[0][size_t(dim * i + c)] = p[c];
            rows[size_t(1 + c)][size_t(dim * i + c)] = 1.0;
        }
    }
    for (size_t l = 0; l < sides.size(); ++l) {
        const double* pc = diagram.GetGenerator(sides[l]._cut);
        for (int i : sides[l]._side) {
            const double* p = diagram.GetGenerator(i);
            for (int c = 0; c < dim; ++c) {
                rows[size_t(1 + dim) + l][size_t(dim * i + c)] = p[c] - pc[c];
            }
        }
    }
    // Vbar-weighted normal equations G gamma = r.
    Eigen::MatrixXd G = Eigen::MatrixXd::Zero(m, m);
    Eigen::VectorXd r = Eigen::VectorXd::Zero(m);
    const std::vector<double>& v = state._velocities;
    for (int i = 0; i < n; ++i) {
        for (int c = 0; c < dim; ++c) {
            const size_t k = size_t(dim * i + c);
            for (int a = 0; a < m; ++a) {
                r[a] += _volumes[size_t(i)] * rows[size_t(a)][k] * v[k];
                for (int b = 0; b < m; ++b) {
                    G(a, b) += _volumes[size_t(i)] * rows[size_t(a)][k] * rows[size_t(b)][k];
                }
            }
        }
    }
    const Eigen::LDLT<Eigen::MatrixXd> ldlt(G);
    if (ldlt.info() != Eigen::Success) {
        return false;
    }
    const Eigen::VectorXd gamma = ldlt.solve(r);
    for (size_t k = 0; k < state._velocities.size(); ++k) {
        for (int a = 0; a < m; ++a) {
            state._velocities[k] -= gamma[a] * rows[size_t(a)][k];
        }
    }
    return true;
}

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
    if (_options._velocityProjection != IntegratorOptions::VelocityProjection::None and
        not _ProjectVelocities(diagram, next)) {
        res._error = "The projection of the initial velocities is singular.";
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
    const std::vector<double> startVelocities = cur._velocities;
    _ApplyGravity(cur._velocities, dt);
    const int dim = _GetDimension();
    const std::vector<double>& q = cur._positions;
    const size_t m = q.size();
    std::vector<double> generators, lambda;
    _GetDiagram(cur, generators, lambda);
    // The cut cells of the dual graph of the generators and their sides, one
    // branch gauge each (latex/flow.tex, eq:flow-line-constraint-arm).
    ShiftedDiagram current;
    if (not _BuildDiagram(cur, current, res._error)) {
        return res;
    }
    const std::vector<CutCellSide> sides = _FindCutCellSides(current);
    {
        // The state at the start of the step on its own diagram.
        ParticleState start = cur;
        start._velocities = startVelocities;
        res._startMeasures = _Measure(current, start);
    }

    // Prediction of eq:flow-shake-step-pred with the potential scaled by kappa = _stiffness, and the
    // force f of IntegratorOptions::_constrainedPotential:
    //   quantization  f = q - p at the current generators p, so that the gauge only corrects the
    //                 remainder (s - kappa dt^2 along p): the same discrete equation, with a
    //                 prediction close to its solution instead of one dilated by kappa dt^2;
    //   volume        f = c0 - q, with c0 the centroids of the particles' own diagram;
    //   softened      f = -(1 + 1/alpha) u_alpha, the shifts of G_alpha with sites and targets at q.
    // With theta = _newmarkTheta > 0 (quantization only), the prediction holds the explicit part
    // (1 - 2 theta) f + theta f_prev of the force, and the implicit part theta kappa dt^2 (q+ - p+),
    // with p+ = q + u, folds into G_alpha: the new particles are a - u / alpha - gauge, with
    // targets a = (pred - beta q) / (1 - beta) and 1 / alpha = beta / (1 - beta), beta = theta kappa dt^2.
    using Potential = IntegratorOptions::ConstrainedPotential;
    const double kappa = _options._stiffness;
    const double theta = _options._newmarkTheta;
    std::vector<double> weights = lambda;
    std::vector<double> force(m);
    auto fail = [&](const DualSolver::Result& r) {
        res._solve = r;
        res._status = DualSolver::GetStatusString(r._status);
        res._error = r._error;
        return res;
    };
    if (theta > 0.0 and _options._constrainedPotential != Potential::Quantization) {
        res._error = "The Newmark weight needs the quantization potential.";
        return res;
    }
    if (_options._constrainedPotential == Potential::Volume) {
        VolumeSolver volumeSolver(_params._solver);
        const DualSolver::Result vr = volumeSolver.Solve(_domain, q, _volumes, weights);
        if (not vr.HasConverged()) {
            return fail(vr);
        }
        const std::vector<double>& c0 = vr._diagram.GetCells()._centroids;
        for (size_t k = 0; k < m; ++k) {
            force[k] = c0[k] - q[k];
        }
    } else if (_options._constrainedPotential == Potential::Softened) {
        std::vector<double> us, ws = lambda;
        DualSolver::Result sr = MakeJointStart(_domain, q, q, _volumes, 0.0, us, ws);
        if (sr.HasConverged()) {
            sr = SoftenedSolver(_params._solver).Solve(_domain, q, q, _volumes, _options._alpha, us, ws);
        }
        if (not sr.HasConverged()) {
            return fail(sr);
        }
        for (size_t k = 0; k < m; ++k) {
            force[k] = -(1.0 + 1.0 / _options._alpha) * us[k];
        }
    } else {
        for (size_t k = 0; k < m; ++k) {
            force[k] = q[k] - generators[k];
        }
    }
    const std::vector<double>& previous = cur._forces;
    std::vector<double> predicted(m);
    for (size_t k = 0; k < m; ++k) {
        predicted[k] = q[k] + dt * cur._velocities[k] +
                       kappa * dt * dt * ((1.0 - 2.0 * theta) * force[k] + theta * previous[k]);
    }
    const double beta = theta * kappa * dt * dt;
    if (beta >= 1.0) {
        res._error = "The Newmark step needs theta kappa dt^2 < 1.";
        return res;
    }
    const double invAlpha = beta / (1.0 - beta);
    std::vector<double> targets = predicted;
    if (theta > 0.0) {
        for (size_t k = 0; k < m; ++k) {
            targets[k] = (predicted[k] - beta * q[k]) / (1.0 - beta);
        }
    }

    std::vector<double> shifts;
    res._solve = MakeJointStart(_domain, q, targets, _volumes, 0.0, shifts, weights);
    if (not res._solve.HasConverged()) {
        // The warm start, the shifted weights of the generators p, may leave empty cells at the
        // sites q when p and q differ; start the volume solve again from zero weights (the
        // Voronoi cells of q, none empty).
        weights.clear();
        shifts.clear();
        res._solve = MakeJointStart(_domain, q, targets, _volumes, 0.0, shifts, weights);
    }
    res._status = DualSolver::GetStatusString(res._solve._status);
    if (not res._solve.HasConverged()) {
        res._error = res._solve._error;
        return res;
    }
    ConstrainedJointSolver solver(_params._solver);
    if (theta > 0.0) {
        solver.SetSoftening(1.0 / invAlpha);
    }
    res._solve = solver.Solve(_domain, q, targets, _volumes, generators, sides, shifts, weights, res._gauge);
    res._status = DualSolver::GetStatusString(res._solve._status);
    res._iterations = res._solve._iterations;
    if (not res._solve.HasConverged()) {
        res._error = res._solve._error;
        return res;
    }

    const double s = res._gauge._dilation;
    const std::vector<double>& t = res._gauge._translation;
    if (_params._solver._log != nullptr) {
        *_params._solver._log << "  gauge: dilation " << s << "  translation";
        for (double x : t) {
            *_params._solver._log << " " << x;
        }
        *_params._solver._log << "  sides " << sides.size();
        for (size_t l = 0; l < sides.size(); ++l) {
            *_params._solver._log << "  (cut " << sides[l]._cut << ", " << sides[l]._side.size() << " cells, s_B "
                                  << res._gauge._arms[l] << ")";
        }
        *_params._solver._log << std::endl;
    }
    // Arm term of the correction: s_B (p_i - p_{c_B}) for i in B_B.
    std::vector<double> arms(m, 0.0);
    for (size_t l = 0; l < sides.size(); ++l) {
        const double sB = res._gauge._arms[l];
        const int cut = sides[l]._cut;
        for (int i : sides[l]._side) {
            for (int c = 0; c < dim; ++c) {
                arms[size_t(dim) * i + c] += sB * (generators[size_t(dim) * i + c] - generators[size_t(dim) * cut + c]);
            }
        }
    }
    ParticleState next;
    next._positions.resize(m);
    next._velocities.resize(m);
    next._pressures = cur._pressures;
    next._forces = force;
    std::vector<double> nextGenerators(m), nextLambda(weights.size());
    for (size_t i = 0; i < weights.size(); ++i) {
        double l = weights[i];
        for (int c = 0; c < dim; ++c) {
            const size_t k = size_t(dim) * i + size_t(c);
            next._positions[k] = targets[k] - invAlpha * shifts[k] - s * generators[k] - t[c] - arms[k];
            next._velocities[k] = (next._positions[k] - q[k]) / dt;
            nextGenerators[k] = q[k] + shifts[k];
            l += shifts[k] * shifts[k] + 2.0 * shifts[k] * (q[k] - targets[k]);
        }
        nextLambda[i] = l;
    }
    _SetDiagram(next, nextGenerators, nextLambda);
    ShiftedDiagram diagram;
    if (not _BuildDiagram(next, diagram, res._error)) {
        return res;
    }
    if (_options._velocityProjection == IntegratorOptions::VelocityProjection::EveryStep and
        not _ProjectVelocities(diagram, next)) {
        res._error = "The projection of the velocities is singular.";
        return res;
    }
    _Finish(diagram, next, res);
    res._converged = true;
    state = std::move(next);
    return res;
}

} // namespace PP
