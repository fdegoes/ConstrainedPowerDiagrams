#ifndef PP_CORE_POWERPARTICLESINTEGRATOR_H
#define PP_CORE_POWERPARTICLESINTEGRATOR_H

#include "core/domainBoundary.h"
#include "core/integrator.h"

namespace PP
{

// Time step of Power Particles (de Goes et al. 2015), as in its 2015 code
// (Before2015/pp, pp2/SimDyn.cpp, Sim::step, fluid particles only, in the
// order of the 2D code): the particles are the sites of a power diagram with
// the prescribed volumes and no shifts.
//   1. v += dt g
//   2. pressure: L p = -div(v), then v -= dt grad(p) / Vbar
//   3. advection: q += beta (c - q) + dt v, pushed back into the domain
//   4. weights: w = 0, then IntegratorOptions::_weightIterations Newton steps
//      of the volumes (stopping early once sum (Vbar - V)^2 < tol), halving
//      a step up to 3 times when it empties a cell
// with, over the interior facets A_ij of the cells (lengths in 2D), l_ij =
// |q_i - q_j|, b_ij the facet centroid and star_ij = A_ij / l_ij:
//   div_i  = sum_j star_ij [(q_j - b_ij) . v_j - (q_i - b_ij) . v_i]
//   grad_i = sum_j star_ij (q_i - b_ij) (p_i - p_j)
//   L_ii = eps V_i + sum_j dt star_ij r_ij,  L_ij = -dt star_ij r_ij,
//   r_ij = (d_ij + d_ji) / (rho_i d_ij + rho_j d_ji)   (2 / (rho_i + rho_j) if < 0)
//   d_ij = (l_ij + (w_i - w_j) / l_ij) / 2,  rho_i = Vbar_i / V_i
// and the Newton Hessian of the weights H_ii = eps V_i + sum_j star_ij / 2,
// H_ij = -star_ij / 2. Boundary facets are skipped (free slip). A particle
// that leaves the domain is moved to its closest boundary point plus r_i (the
// radius of the ball of volume Vbar_i) along the push-back direction, and
// stays at q + beta (c - q) if that is still outside.
//
// Differences with the 2015 code: the fluid fills the domain (no air
// particles, hence no free surface), and no solids, viscosity, temperature or
// compressible state equation; sparse Cholesky instead of its PCG for the
// pressures and Newton steps (same systems); the pressures are kept in
// ParticleState::_pressures; the caller chooses dt (ComputeCFLStep gives the
// code's CFL step). Initialize() moves the particles to a centroidal diagram.
class PowerParticlesIntegrator : public Integrator
{
public:
    PowerParticlesIntegrator(
        const Domain& domain,
        std::vector<double> volumes,
        const Parameters& params,
        const IntegratorOptions& options);

    const char* GetName() const override { return "PowerParticles"; }
    Result Initialize(ParticleState& state) const override;
    Result Step(ParticleState& state, double dt) const override;
    double GetPotential(double quantization, double centroidResidual) const override
    {
        (void)quantization;
        (void)centroidResidual;
        return 0.0;
    }

    // The discrete operators of the step on a diagram (for the tests).
    // divergence  n values; gradient dim n values.
    void ComputeDivergence(const ShiftedDiagram& diagram, const std::vector<double>& q,
                           const std::vector<double>& velocities, std::vector<double>& divergence) const;
    void ComputeGradient(const ShiftedDiagram& diagram, const std::vector<double>& q,
                         const std::vector<double>& pressures, std::vector<double>& gradient) const;

    const DomainBoundary& GetBoundary() const { return _boundary; }

private:
    double _GetRegularization() const;
    // The power diagram of the sites q with weights w (no shifts).
    bool _Build(const std::vector<double>& q, const std::vector<double>& w, ShiftedDiagram& diagram,
                std::string& err) const;
    // Pressure solve and velocity update (step 2).
    bool _ApplyPressure(const ShiftedDiagram& diagram, ParticleState& state, double dt, Result& res) const;
    // Advection (step 3).
    void _Advect(const ShiftedDiagram& diagram, ParticleState& state, double dt, Result& res) const;
    // Weights for the volumes (step 4); diagram is the final one.
    bool _SolveWeights(ParticleState& state, ShiftedDiagram& diagram, Result& res) const;

    DomainBoundary _boundary;
};

} // namespace PP

#endif // PP_CORE_POWERPARTICLESINTEGRATOR_H
