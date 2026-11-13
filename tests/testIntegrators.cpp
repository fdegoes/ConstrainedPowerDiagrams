// Standalone checks of the centroidal power diagram generator and of the
// integrators (core/integrator.h): constrained and softened (latex/flow.tex),
// GM (Gallouet & Merigot 2018) and Power Particles (de Goes et al. 2015).

#include "core/centroidalDiagram.h"
#include "core/constrainedIntegrator.h"
#include "core/domain.h"
#include "core/geometry.h"
#include "core/integrator.h"
#include "core/powerParticlesIntegrator.h"
#include "core/shiftedDiagram.h"
#include "testUtils.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

namespace
{

using PP::CentroidalDiagram;
using PP::Domain;
using PP::DualSolver;
using PP::Integrator;
using PP::IntegratorKind;
using PP::ParticleState;
using PPTest::Check;
using PPTest::FormatValue;

// Power diagram of the state: generators q + u, power weights w + |u|^2.
PP::ShiftedDiagram
_BuildStateDiagram(const Domain& D, const ParticleState& s)
{
    PP::ShiftedDiagram d;
    std::string err;
    d.Build(D, s._positions, s._positions, s._shifts, s._weights, err);
    return d;
}

// Tolerance of the solves of the time steps (Integrator::Parameters): the
// volume, centroid and momentum errors of a run are of its order.
const double _tol = PP::Integrator::Parameters()._solver._volumeTolerance;

double
_GetVolumeError(const PP::ShiftedDiagram& d, const std::vector<double>& volumes, double D)
{
    double e = 0.0;
    for (size_t i = 0; i < volumes.size(); ++i) {
        e = std::max(e, std::fabs(d.GetVolume(int(i)) - volumes[i]) / D);
    }
    return e;
}

// max_i |c_i - q_i|.
double
_GetCentroidError(const PP::ShiftedDiagram& d, const std::vector<double>& q)
{
    const int dim = d.GetDimension();
    double e2 = 0.0;
    for (int i = 0; i < d.GetNumSites(); ++i) {
        e2 = std::max(e2, PP::ComputeDistance2(d.GetCentroid(i), &q[dim * i], dim));
    }
    return std::sqrt(e2);
}

// |sum_i Vbar_i v_i|_inf, the total momentum.
double
_GetMomentum(const ParticleState& s, const std::vector<double>& volumes, int dim)
{
    std::vector<double> m(dim, 0.0);
    for (size_t i = 0; i < volumes.size(); ++i) {
        for (int c = 0; c < dim; ++c) {
            m[c] += volumes[i] * s._velocities[dim * i + c];
        }
    }
    double e = 0.0;
    for (double x : m) {
        e = std::max(e, std::fabs(x));
    }
    return e;
}

// Centroidal diagram with prescribed volumes, from random sites.
bool
_MakeCentroidalState(
    const Domain& D,
    const std::string& name,
    int n,
    std::vector<double>& volumes,
    ParticleState& state)
{
    const int dim = D.GetDimension();
    std::vector<double> sites = PPTest::RandomSites(n, dim, 0.1, 0.9, 51);
    volumes = PPTest::ScaleFractions(PPTest::RandomFractions(n, 52), D.GetVolume());
    std::vector<double> w;
    CentroidalDiagram::Parameters p;
    p._tolerance = 1e-9;
    const CentroidalDiagram::Result r = CentroidalDiagram(p).Solve(D, volumes, sites, w);
    Check(r._converged, name + ": centroidal power diagram converges (" + std::to_string(r._iterations) +
                            " Lloyd iterations, residual " + FormatValue(r._residual) + ")");
    state = CentroidalDiagram::MakeInitialState(sites, w);
    const PP::ShiftedDiagram d = _BuildStateDiagram(D, state);
    Check(_GetVolumeError(d, volumes, D.GetVolume()) < 1e-9 and _GetCentroidError(d, sites) < 1e-8,
          name + ": the initial state has the prescribed volumes and centroids at the particles");
    return r._converged;
}

// Random initial velocities with zero total momentum.
void
_SetRandomVelocities(ParticleState& s, const std::vector<double>& volumes, int dim, double scale)
{
    s._velocities = PPTest::RandomValues(s._positions.size(), -scale, scale, 61);
    std::vector<double> m(dim, 0.0);
    double total = 0.0;
    for (size_t i = 0; i < volumes.size(); ++i) {
        total += volumes[i];
        for (int c = 0; c < dim; ++c) {
            m[c] += volumes[i] * s._velocities[dim * i + c];
        }
    }
    for (size_t i = 0; i < volumes.size(); ++i) {
        for (int c = 0; c < dim; ++c) {
            s._velocities[dim * i + c] -= m[c] / total;
        }
    }
}

void
_TestConstrained(const Domain& D, const std::string& name, int n, int steps, double dt)
{
    const int dim = D.GetDimension();
    std::vector<double> volumes;
    ParticleState state;
    if (not _MakeCentroidalState(D, name, n, volumes, state)) {
        return;
    }
    _SetRandomVelocities(state, volumes, dim, 0.05);

    PP::ConstrainedIntegrator integrator(D, volumes, Integrator::Parameters(), PP::IntegratorOptions());
    double vErr = 0.0, cErr = 0.0, momentum = 0.0;
    int iterations = 0, done = 0;
    std::string failure;
    for (; done < steps; ++done) {
        const Integrator::Result r = integrator.Step(state, dt);
        if (not r.HasConverged()) {
            failure = std::string(DualSolver::GetStatusString(r._solve._status)) +
                      (r._solve._error.empty() ? "" : ": " + r._solve._error);
            break;
        }
        iterations += r._solve._iterations;
        const PP::ShiftedDiagram d = _BuildStateDiagram(D, state);
        vErr = std::max(vErr, _GetVolumeError(d, volumes, D.GetVolume()));
        cErr = std::max(cErr, _GetCentroidError(d, state._positions));
        momentum = std::max(momentum, _GetMomentum(state, volumes, dim));
    }
    Check(done == steps, name + " constrained: " + std::to_string(done) + "/" + std::to_string(steps) +
                             " steps (dt " + FormatValue(dt) + ", " + std::to_string(iterations) +
                             " Newton iterations" + (failure.empty() ? "" : ", stopped: " + failure) + ")");
    Check(vErr < _tol, name + " constrained: volumes within the tolerance (max rel err " + FormatValue(vErr) + ")");
    Check(cErr < _tol, name + " constrained: particles stay at the centroids (max err " + FormatValue(cErr) + ")");
    Check(momentum < _tol, name + " constrained: total momentum stays small (" + FormatValue(momentum) + ")");
}

void
_TestPenalty(const Domain& D, const std::string& name, int n, int steps, double alpha)
{
    const int dim = D.GetDimension();
    std::vector<double> volumes;
    ParticleState state;
    if (not _MakeCentroidalState(D, name, n, volumes, state)) {
        return;
    }
    _SetRandomVelocities(state, volumes, dim, 0.05);
    // Within the stability bound eq:softened-time-step.
    const double dt = 0.5 * 2.0 / std::sqrt(1.0 + alpha);

    PP::IntegratorOptions o;
    o._alpha = alpha;
    const auto integrator = PP::MakeIntegrator(IntegratorKind::Softened, D, volumes, Integrator::Parameters(), o);
    double vErr = 0.0, momentum = 0.0, maxSpeed = 0.0;
    int iterations = 0, done = 0;
    std::string failure;
    for (; done < steps; ++done) {
        const Integrator::Result r = integrator->Step(state, dt * 0.05);
        if (not r.HasConverged()) {
            failure = std::string(DualSolver::GetStatusString(r._solve._status)) +
                      (r._solve._error.empty() ? "" : ": " + r._solve._error);
            break;
        }
        iterations += r._solve._iterations;
        vErr = std::max(vErr, _GetVolumeError(r._solve._diagram, volumes, D.GetVolume()));
        momentum = std::max(momentum, _GetMomentum(state, volumes, dim));
        for (double v : state._velocities) {
            maxSpeed = std::max(maxSpeed, std::fabs(v));
        }
    }
    Check(done == steps, name + " softened (alpha " + FormatValue(alpha) + "): " + std::to_string(done) + "/" +
                             std::to_string(steps) + " steps (" + std::to_string(iterations) +
                             " Newton iterations" + (failure.empty() ? "" : ", stopped: " + failure) + ")");
    Check(vErr < _tol, name + " softened: volumes of each solve within the tolerance (max rel err " + FormatValue(vErr) + ")");
    Check(momentum < _tol, name + " softened: total momentum stays small (" + FormatValue(momentum) + ")");
    Check(maxSpeed < 1.0, name + " softened: velocities stay bounded (max " + FormatValue(maxSpeed) + ")");
}

// A generic run of an integrator from Initialize(): volumes within volTol,
// no failed step, finite energy, and (momTol > 0) zero total momentum
// without gravity; Power Particles does not conserve it (its gradient is
// not antisymmetric over a facet).
void
_TestRun(const Domain& D, const std::string& name, IntegratorKind kind, const PP::IntegratorOptions& o, int n,
         int steps, double dt, double volTol, double momTol)
{
    const int dim = D.GetDimension();
    std::vector<double> volumes;
    ParticleState state;
    if (not _MakeCentroidalState(D, name, n, volumes, state)) {
        return;
    }
    _SetRandomVelocities(state, volumes, dim, 0.05);
    const auto integrator = PP::MakeIntegrator(kind, D, volumes, Integrator::Parameters(), o);
    const std::string label = name + " " + integrator->GetName();
    Integrator::Result r0 = integrator->Initialize(state);
    Check(r0.HasConverged(), label + ": Initialize" + (r0._error.empty() ? "" : " (" + r0._error + ")"));
    if (not r0.HasConverged()) {
        return;
    }
    double vErr = 0.0, momentum = 0.0, energy = 0.0;
    int done = 0, iterations = 0, projected = 0;
    std::string failure;
    for (; done < steps; ++done) {
        const Integrator::Result r = integrator->Step(state, dt);
        if (not r.HasConverged()) {
            failure = r._status + (r._error.empty() ? "" : ": " + r._error);
            break;
        }
        iterations += r._iterations;
        projected += r._numProjected;
        vErr = std::max(vErr, r._maxVolumeError);
        momentum = std::max(momentum, _GetMomentum(state, volumes, dim));
        energy = r._measures._kinetic;
        // The result describes the diagram of the state.
        const PP::ShiftedDiagram d = _BuildStateDiagram(D, state);
        vErr = std::max(vErr, _GetVolumeError(d, volumes, D.GetVolume()) * D.GetVolume() /
                                  *std::min_element(volumes.begin(), volumes.end()));
    }
    Check(done == steps, label + ": " + std::to_string(done) + "/" + std::to_string(steps) + " steps (dt " +
                             FormatValue(dt) + ", " + std::to_string(iterations) + " Newton iterations, " +
                             std::to_string(projected) + " projected" + (failure.empty() ? "" : ", stopped: " + failure) +
                             ")");
    Check(vErr < volTol, label + ": volumes within " + FormatValue(volTol) + " (max rel err " + FormatValue(vErr) + ")");
    Check(std::isfinite(energy) and energy < 1.0, label + ": kinetic energy stays bounded (" + FormatValue(energy) + ")");
    if (momTol > 0.0) {
        Check(momentum < momTol, label + ": total momentum stays small (" + FormatValue(momentum) + ")");
    }
}

// Gravity: a state at rest gains the velocity dt g (in the mean, since the
// fluid fills the domain and the walls push back) and the volumes hold.
void
_TestGravity(const Domain& D, const std::string& name, IntegratorKind kind, const PP::IntegratorOptions& o)
{
    const int dim = D.GetDimension();
    std::vector<double> volumes;
    ParticleState state;
    if (not _MakeCentroidalState(D, name, 30, volumes, state)) {
        return;
    }
    Integrator::Parameters params;
    params._gravity.assign(size_t(dim), 0.0);
    params._gravity[1] = -9.8;
    const auto integrator = PP::MakeIntegrator(kind, D, volumes, params, o);
    const std::string label = name + " " + integrator->GetName();
    if (not integrator->Initialize(state).HasConverged()) {
        Check(false, label + " gravity: Initialize");
        return;
    }
    const double dt = 0.01;
    const Integrator::Result r = integrator->Step(state, dt);
    // Mean vertical velocity: -9.8 dt minus what the pressure (walls) takes
    // back; it must point down and be at most the free-fall value.
    double vy = 0.0, total = 0.0;
    for (size_t i = 0; i < volumes.size(); ++i) {
        vy += volumes[i] * state._velocities[size_t(dim) * i + 1];
        total += volumes[i];
    }
    vy /= total;
    Check(r.HasConverged() and vy < 0.0 and vy >= -9.8 * dt * (1.0 + 1e-9),
          label + " gravity: the particles fall (mean v_y " + FormatValue(vy) + ", free fall " + FormatValue(-9.8 * dt) +
              ")" + (r._error.empty() ? "" : " " + r._error));
}

// Energy and momentum: the measures of a hand-built state, and the drift
// of the total energy and momentum over a run (written as a CSV for
// plotting when PP_TEST_CSV_DIR is set).
void
_TestMeasures(const Domain& D, const std::string& name, IntegratorKind kind, const PP::IntegratorOptions& o, double dt,
              int steps, double energyTol)
{
    const int dim = D.GetDimension();
    std::vector<double> volumes;
    ParticleState state;
    if (not _MakeCentroidalState(D, name, 30, volumes, state)) {
        return;
    }
    _SetRandomVelocities(state, volumes, dim, 0.05);
    const auto integrator = PP::MakeIntegrator(kind, D, volumes, Integrator::Parameters(), o);
    const std::string label = name + " " + integrator->GetName() + " measures";
    if (not integrator->Initialize(state).HasConverged()) {
        Check(false, label + ": Initialize");
        return;
    }
    PP::FlowMeasures m0;
    std::string err;
    integrator->Measure(state, m0, err);
    // Kinetic energy and momentum by hand.
    double K = 0.0;
    std::vector<double> P(size_t(dim), 0.0);
    for (size_t i = 0; i < volumes.size(); ++i) {
        for (int c = 0; c < dim; ++c) {
            const double v = state._velocities[size_t(dim) * i + c];
            K += volumes[i] * v * v;
            P[c] += volumes[i] * v;
        }
    }
    double pErr = 0.0;
    for (int c = 0; c < dim; ++c) {
        pErr = std::max(pErr, std::fabs(P[c] - m0._momentum[c]));
    }
    Check(std::fabs(K - m0._kinetic) < 1e-14 and pErr < 1e-14 and m0._quantization > 0.0 and
              m0._centroidResidual >= 0.0,
          label + ": kinetic energy, momentum, quantization error (" + FormatValue(m0._quantization) + ")");

    std::FILE* csv = nullptr;
    if (const char* dir = std::getenv("PP_TEST_CSV_DIR")) {
        std::string path = std::string(dir) + "/" + name + "-" + integrator->GetName() + ".csv";
        std::replace(path.begin() + long(std::string(dir).size()) + 1, path.end(), ' ', '_');
        csv = std::fopen(path.c_str(), "w");
        if (csv) {
            std::fprintf(csv, "%s\n%s\n", PP::FlowMeasures::GetCSVHeader(dim).c_str(), m0.GetCSVRow(0.0).c_str());
        }
    }
    double drift = 0.0, pDrift = 0.0;
    int done = 0;
    for (; done < steps; ++done) {
        const Integrator::Result r = integrator->Step(state, dt);
        if (not r.HasConverged()) {
            break;
        }
        const PP::FlowMeasures& m = r._measures;
        drift = std::max(drift, std::fabs(m._total - m0._total) / std::max(m0._total, 1e-300));
        for (int c = 0; c < dim; ++c) {
            pDrift = std::max(pDrift, std::fabs(m._momentum[c] - m0._momentum[c]));
        }
        if (csv) {
            std::fprintf(csv, "%s\n", m.GetCSVRow(dt * (done + 1)).c_str());
        }
    }
    if (csv) {
        std::fclose(csv);
    }
    Check(done == steps and drift < energyTol,
          label + ": total energy drift " + FormatValue(drift) + " over " + std::to_string(done) + " steps (dt " +
              FormatValue(dt) + ", E0 " + FormatValue(m0._total) + ")");
    Check(pDrift < _tol, label + ": momentum conserved (" + FormatValue(pDrift) + ")");
}

// GM with kappa = 1 is the limit alpha -> 0 of the softened integrator.
void
_TestGMLimit(const Domain& D, const std::string& name)
{
    std::vector<double> volumes;
    ParticleState state;
    if (not _MakeCentroidalState(D, name, 30, volumes, state)) {
        return;
    }
    _SetRandomVelocities(state, volumes, D.GetDimension(), 0.05);
    // Away from the centroids (at a centroidal state no scheme applies a force).
    const std::vector<double> offset = PPTest::RandomValues(state._positions.size(), -0.01, 0.01, 71);
    for (size_t k = 0; k < offset.size(); ++k) {
        state._positions[k] += offset[k];
    }
    PP::IntegratorOptions o;
    const auto gm = PP::MakeIntegrator(IntegratorKind::GM, D, volumes, Integrator::Parameters(), o);
    ParticleState g = state;
    gm->Step(g, 0.05);
    double prev = HUGE_VAL, first = 0.0;
    bool decreasing = true;
    std::string diffs;
    for (double alpha : {1e-1, 1e-2, 1e-3}) {
        o._alpha = alpha;
        const auto soft = PP::MakeIntegrator(IntegratorKind::Softened, D, volumes, Integrator::Parameters(), o);
        ParticleState s = state;
        soft->Step(s, 0.05);
        double d = 0.0;
        for (size_t k = 0; k < s._velocities.size(); ++k) {
            d = std::max(d, std::fabs(s._velocities[k] - g._velocities[k]));
        }
        decreasing = decreasing and d < prev;
        first = first > 0.0 ? first : d;
        prev = d;
        diffs += (diffs.empty() ? "" : ", ") + FormatValue(d);
    }
    // O(alpha): 100x smaller from alpha = 0.1 to 0.001.
    Check(decreasing and prev < 0.02 * first,
          name + " GM: softened velocities tend to GM as alpha -> 0 (max diff " + diffs + ")");
}

// Power Particles: the divergence of a constant velocity and the gradient of
// a constant pressure vanish; a particle pushed outside the domain is put
// back inside.
void
_TestPowerParticlesOperators(const Domain& D, const std::string& name)
{
    const int dim = D.GetDimension();
    std::vector<double> volumes;
    ParticleState state;
    if (not _MakeCentroidalState(D, name, 30, volumes, state)) {
        return;
    }
    PP::PowerParticlesIntegrator pp(D, volumes, Integrator::Parameters(), PP::IntegratorOptions());
    const PP::ShiftedDiagram d = _BuildStateDiagram(D, state);
    std::vector<double> v(state._positions.size()), div, grad;
    for (size_t k = 0; k < v.size(); ++k) {
        v[k] = 0.3 + 0.1 * double(k % size_t(dim));
    }
    pp.ComputeDivergence(d, state._positions, v, div);
    double divErr = 0.0;
    for (double x : div) {
        divErr = std::max(divErr, std::fabs(x));
    }
    // A constant field is divergence-free in the interior; on cells touching
    // the boundary the code's divergence omits the boundary facets, so only
    // the sum (a telescoping sum over interior facets) is checked.
    double divSum = 0.0;
    for (double x : div) {
        divSum += x;
    }
    pp.ComputeGradient(d, state._positions, std::vector<double>(volumes.size(), 2.5), grad);
    double gradErr = 0.0;
    for (double x : grad) {
        gradErr = std::max(gradErr, std::fabs(x));
    }
    Check(std::fabs(divSum) < 1e-12, name + " PP: the divergence sums to zero (" + FormatValue(divSum) + ")");
    Check(gradErr < 1e-12, name + " PP: the gradient of a constant pressure vanishes (" + FormatValue(gradErr) + ")");

    // At rest at a centroidal diagram, without gravity, nothing moves.
    ParticleState rest = state;
    rest._velocities.assign(rest._positions.size(), 0.0);
    const Integrator::Result rr = pp.Step(rest, 0.05);
    double moved = 0.0;
    for (size_t k = 0; k < rest._positions.size(); ++k) {
        moved = std::max(moved, std::fabs(rest._positions[k] - state._positions[k]));
    }
    Check(rr.HasConverged() and moved < 1e-8, name + " PP: a centroidal state at rest is a fixed point (moved " +
                                                   FormatValue(moved) + ")");

    const PP::DomainBoundary& B = pp.GetBoundary();
    std::vector<double> inside(static_cast<size_t>(dim), 0.5), outside(static_cast<size_t>(dim), 0.5);
    std::vector<double> proj(static_cast<size_t>(dim), 0.0);
    outside[0] = 1.3;
    B.Project(outside.data(), proj.data());
    Check(B.Evaluate(inside.data()) < 0.0 and B.Evaluate(outside.data()) > 0.0 and std::fabs(proj[0] - 1.0) < 1e-12,
          name + " PP: the domain boundary level set (" + std::to_string(B.GetNumFaces()) + " faces)");
}

} // namespace

int
main()
{
    Domain square, cube;
    PPTest::MakeUnitSquare(square);
    PPTest::MakeUnitCube(cube);

    _TestConstrained(square, "2D square", 40, 20, 0.01);
    _TestConstrained(cube, "3D cube", 30, 10, 0.01);
    _TestPenalty(square, "2D square", 40, 20, 0.001);
    _TestPenalty(cube, "3D cube", 30, 10, 0.001);

    PP::IntegratorOptions softened;
    softened._alpha = 0.01;
    PP::IntegratorOptions gm;
    PP::IntegratorOptions pp;
    for (const auto& [D, name] : {std::pair<const Domain*, std::string>{&square, "2D square"},
                                  std::pair<const Domain*, std::string>{&cube, "3D cube"}}) {
        // The volume error of a run is relative to the smallest volume, so
        // it may exceed the per-cell relative tolerance of the solves.
        _TestRun(*D, name, IntegratorKind::Softened, softened, 30, 10, 0.05, 10.0 * _tol, _tol);
        _TestRun(*D, name, IntegratorKind::GM, gm, 30, 10, 0.05, 10.0 * _tol, _tol);
        _TestRun(*D, name, IntegratorKind::PowerParticles, pp, 30, 10, 0.05, 0.05, 0.0);
        _TestGravity(*D, name, IntegratorKind::Softened, softened);
        _TestGravity(*D, name, IntegratorKind::GM, gm);
        _TestGravity(*D, name, IntegratorKind::PowerParticles, pp);
        _TestMeasures(*D, name, IntegratorKind::Softened, softened, 0.01, 50, 0.05);
        _TestMeasures(*D, name, IntegratorKind::GM, gm, 0.01, 50, 0.05);
        _TestGMLimit(*D, name);
        _TestPowerParticlesOperators(*D, name);
    }
    return PPTest::Report();
}
