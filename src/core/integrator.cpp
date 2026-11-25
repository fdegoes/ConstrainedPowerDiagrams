#include "core/integrator.h"
#include "core/centroidalDiagram.h"
#include "core/constrainedIntegrator.h"
#include "core/domain.h"
#include "core/geometry.h"
#include "core/powerParticlesIntegrator.h"
#include "core/gmIntegrator.h"
#include "core/softenedIntegrator.h"

#include <algorithm>
#include <cmath>
#include <cstdio>

namespace PP
{

Integrator::Integrator(
    const Domain& domain,
    std::vector<double> volumes,
    const Parameters& params,
    const IntegratorOptions& options)
    : _domain(domain), _volumes(std::move(volumes)), _params(params), _options(options)
{ }

int
Integrator::_GetDimension() const
{
    return _domain.GetDimension();
}

bool
Integrator::_Prepare(const ParticleState& in, ParticleState& out, Result& res) const
{
    if (not _domain.IsValid()) {
        res._error = "Invalid clipping domain.";
        return false;
    }
    const size_t n = _volumes.size();
    const size_t dim = size_t(_GetDimension());
    if (n == 0 or in._positions.size() != dim * n) {
        res._error = "The positions do not match the volumes (" + std::to_string(n) + " particles).";
        return false;
    }
    auto fill = [&](const std::vector<double>& v, size_t size, std::vector<double>& o, const char* name) {
        if (v.empty()) {
            o.assign(size, 0.0);
            return true;
        }
        if (v.size() != size) {
            res._error = std::string("The ") + name + " do not match the positions.";
            return false;
        }
        o = v;
        return true;
    };
    out._positions = in._positions;
    return fill(in._velocities, dim * n, out._velocities, "velocities") and
           fill(in._shifts, dim * n, out._shifts, "shifts") and fill(in._weights, n, out._weights, "weights") and
           fill(in._pressures, n, out._pressures, "pressures") and fill(in._forces, dim * n, out._forces, "forces");
}

void
Integrator::_ApplyGravity(std::vector<double>& velocities, double dt) const
{
    const std::vector<double>& g = _params._gravity;
    const int dim = _GetDimension();
    if (g.size() != size_t(dim)) {
        return;
    }
    for (size_t k = 0; k < velocities.size(); ++k) {
        velocities[k] += dt * g[k % size_t(dim)];
    }
}

void
Integrator::_GetDiagram(const ParticleState& state, std::vector<double>& p, std::vector<double>& lambda) const
{
    const int dim = _GetDimension();
    const size_t n = state._weights.size();
    p.resize(state._positions.size());
    lambda.resize(n);
    for (size_t i = 0; i < n; ++i) {
        double l = state._weights[i];
        for (int c = 0; c < dim; ++c) {
            const size_t k = size_t(dim) * i + size_t(c);
            p[k] = state._positions[k] + state._shifts[k];
            l += state._shifts[k] * state._shifts[k];
        }
        lambda[i] = l;
    }
}

void
Integrator::_SetDiagram(ParticleState& state, const std::vector<double>& p, const std::vector<double>& lambda) const
{
    const int dim = _GetDimension();
    const size_t n = lambda.size();
    state._shifts.resize(p.size());
    state._weights.resize(n);
    for (size_t i = 0; i < n; ++i) {
        double w = lambda[i];
        for (int c = 0; c < dim; ++c) {
            const size_t k = size_t(dim) * i + size_t(c);
            state._shifts[k] = p[k] - state._positions[k];
            w -= state._shifts[k] * state._shifts[k];
        }
        state._weights[i] = w;
    }
}

bool
Integrator::_BuildDiagram(const ParticleState& state, ShiftedDiagram& diagram, std::string& err) const
{
    return diagram.Build(_domain, state._positions, state._positions, state._shifts, state._weights, err);
}

void
Integrator::_Finish(const ShiftedDiagram& diagram, const ParticleState& state, Result& res) const
{
    const int dim = _GetDimension();
    const size_t n = _volumes.size();
    res._centroids = diagram.GetCells()._centroids;
    res._cellVolumes = diagram.GetCells()._volumes;
    res._maxVolumeError = 0.0;
    res._maxCentroidGap = 0.0;
    for (size_t i = 0; i < n; ++i) {
        res._maxVolumeError = std::max(res._maxVolumeError, std::fabs(res._cellVolumes[i] - _volumes[i]) / _volumes[i]);
        res._maxCentroidGap = std::max(
            res._maxCentroidGap, std::sqrt(ComputeDistance2(&state._positions[dim * i], &res._centroids[dim * i], dim)));
    }
    res._measures = _Measure(diagram, state);
}

FlowMeasures
Integrator::_Measure(const ShiftedDiagram& diagram, const ParticleState& state) const
{
    const int dim = _GetDimension();
    const size_t n = _volumes.size();
    const std::vector<double>& g = _params._gravity;
    const std::vector<double>& o = _domain.GetCentroid();
    const PowerCells& cells = diagram.GetCells();
    FlowMeasures m;
    m._momentum.assign(size_t(dim), 0.0);
    m._angularMomentum.assign(dim == 2 ? 1 : 3, 0.0);
    for (size_t i = 0; i < n; ++i) {
        const double* q = &state._positions[size_t(dim) * i];
        const double* v = &state._velocities[size_t(dim) * i];
        const double* c = &cells._centroids[size_t(dim) * i];
        const double* p = diagram.GetGenerator(int(i));
        const double V = cells._volumes[i];
        // _secondMoments is about the generator: ∫|x - p|^2 = ∫|x - c|^2 + V |c - p|^2.
        m._quantization += cells._secondMoments[i] - V * ComputeDistance2(c, p, dim);
        m._centroidResidual += _volumes[i] * ComputeDistance2(c, q, dim);
        double v2 = 0.0;
        for (int k = 0; k < dim; ++k) {
            v2 += v[k] * v[k];
            m._momentum[k] += _volumes[i] * v[k];
            if (g.size() == size_t(dim)) {
                m._gravity -= 2.0 * _volumes[i] * g[k] * q[k];
            }
        }
        m._kinetic += _volumes[i] * v2;
        const double r[3] = {q[0] - o[0], q[1] - o[1], dim == 3 ? q[2] - o[2] : 0.0};
        const double w[3] = {v[0], v[1], dim == 3 ? v[2] : 0.0};
        if (dim == 2) {
            m._angularMomentum[0] += _volumes[i] * (r[0] * w[1] - r[1] * w[0]);
        } else {
            m._angularMomentum[0] += _volumes[i] * (r[1] * w[2] - r[2] * w[1]);
            m._angularMomentum[1] += _volumes[i] * (r[2] * w[0] - r[0] * w[2]);
            m._angularMomentum[2] += _volumes[i] * (r[0] * w[1] - r[1] * w[0]);
        }
    }
    m._potential = GetPotential(m._quantization, m._centroidResidual);
    m._total = m._kinetic + m._potential + m._gravity;
    return m;
}

bool
Integrator::Measure(const ParticleState& state, FlowMeasures& measures, std::string& err) const
{
    Result res;
    ParticleState s;
    if (not _Prepare(state, s, res)) {
        err = res._error;
        return false;
    }
    ShiftedDiagram diagram;
    if (not _BuildDiagram(s, diagram, err)) {
        return false;
    }
    measures = _Measure(diagram, s);
    return true;
}

std::string
FlowMeasures::GetCSVHeader(int dim)
{
    std::string h = "t,kinetic,potential,gravity,total,quantization,centroidResidual";
    const char* axes = "xyz";
    for (int k = 0; k < dim; ++k) {
        h += std::string(",momentum_") + axes[k];
    }
    if (dim == 2) {
        h += ",angularMomentum";
    } else {
        h += ",angularMomentum_x,angularMomentum_y,angularMomentum_z";
    }
    return h;
}

std::string
FlowMeasures::GetCSVRow(double t) const
{
    char buf[64];
    std::string row;
    auto add = [&](double x) {
        std::snprintf(buf, sizeof(buf), "%.17g", x);
        row += (row.empty() ? "" : ",") + std::string(buf);
    };
    add(t);
    for (double x : {_kinetic, _potential, _gravity, _total, _quantization, _centroidResidual}) {
        add(x);
    }
    for (double x : _momentum) {
        add(x);
    }
    for (double x : _angularMomentum) {
        add(x);
    }
    return row;
}

bool
Integrator::_MakeCentroidal(ParticleState& state, Result& res) const
{
    CentroidalDiagram::Parameters cp;
    // A static solve: the volumes to the tolerance of the static solves,
    // not the looser one of the time steps.
    cp._volumeSolver = _params._solver;
    cp._volumeSolver._volumeTolerance = DualSolver::Parameters()._volumeTolerance;
    cp._volumeSolver._shiftTolerance = DualSolver::Parameters()._shiftTolerance;
    cp._tolerance = _params._lloydTolerance;
    cp._maxIterations = _params._lloydMaxIterations;
    std::vector<double> sites = state._positions;
    // Warm start from the weights of the state: for inputs that are already centroidal, the
    // diagram of (sites, weights) is the start, and the Lloyd iterations need not move.
    std::vector<double> weights = state._weights;
    const CentroidalDiagram::Result cr = CentroidalDiagram(cp).Solve(_domain, _volumes, sites, weights);
    res._iterations = cr._iterations;
    if (not cr._converged) {
        res._error = cr._error.empty() ? "The Lloyd iterations did not converge." : cr._error;
        return false;
    }
    state._positions = sites;
    state._shifts.assign(sites.size(), 0.0);
    state._weights = weights;
    return true;
}

std::unique_ptr<Integrator>
MakeIntegrator(
    IntegratorKind kind,
    const Domain& domain,
    std::vector<double> volumes,
    const Integrator::Parameters& params,
    const IntegratorOptions& options)
{
    switch (kind) {
    case IntegratorKind::Constrained:
        return std::make_unique<ConstrainedIntegrator>(domain, std::move(volumes), params, options);
    case IntegratorKind::Softened:
        return std::make_unique<SoftenedIntegrator>(domain, std::move(volumes), params, options);
    case IntegratorKind::GM:
        return std::make_unique<GMIntegrator>(domain, std::move(volumes), params, options);
    case IntegratorKind::PowerParticles:
        return std::make_unique<PowerParticlesIntegrator>(domain, std::move(volumes), params, options);
    }
    return nullptr;
}

double
ComputeCFLStep(const std::vector<double>& velocities, const std::vector<double>& volumes, int dim, double cfl)
{
    double vmax = 0.0;
    for (double v : velocities) {
        vmax = std::max(vmax, std::fabs(v));
    }
    if (not (vmax > 0.0) or volumes.empty()) {
        return 0.0;
    }
    const double vmin = *std::min_element(volumes.begin(), volumes.end());
    // Radius of the ball of volume vmin.
    const double rmin = dim == 2 ? std::sqrt(vmin / M_PI) : std::cbrt(3.0 * vmin / (4.0 * M_PI));
    return cfl * rmin / vmax;
}

} // namespace PP
