#ifndef PP_CORE_SOFTENEDINTEGRATOR_H
#define PP_CORE_SOFTENEDINTEGRATOR_H

#include "core/integrator.h"

namespace PP
{

// Time step of the softened (penalty-based) integrator (latex/flow.tex,
// sec:flow-penalty; latex/appendix-algorithms.tex, Algorithm "flow-penalty"):
// the particles are the sites of the power diagram, and the centroid
// constraints are softened with stiffness alpha (IntegratorOptions::_alpha).
//   (u, w)  = SoftenedSolver with sites and targets q^k, from u = 0 and the
//             weights of F at q^k, warm-started from lambda^k
//   v^{k+1} = v^k - kappa (1 + 1/alpha) dt u
//             (= v^k + kappa (1 + alpha) dt (c - q^k) at convergence)
//   q^{k+1} = q^k + dt v^{k+1}
//   p^{k+1} = q^k + u,  lambda^{k+1} = w + |u|^2
// with kappa = IntegratorOptions::_stiffness (1 in latex/flow.tex). Stable
// for dt < 2 / sqrt(kappa (1 + alpha)) (eq:softened-time-step).
// Initialize() solves H_alpha with sites and targets at the positions.
class SoftenedIntegrator : public Integrator
{
public:
    using Integrator::Integrator;

    const char* GetName() const override { return "Softened"; }
    Result Initialize(ParticleState& state) const override;
    Result Step(ParticleState& state, double dt) const override;
    double GetPotential(double quantization, double centroidResidual) const override
    {
        return _options._stiffness * (quantization + (1.0 + _options._alpha) * centroidResidual);
    }

private:
    // The H_alpha solve at the positions of cur, warm-started from its
    // diagram: shifts and weights on return.
    bool _Solve(const ParticleState& cur, std::vector<double>& shifts, std::vector<double>& weights,
                Result& res) const;
};

} // namespace PP

#endif // PP_CORE_SOFTENEDINTEGRATOR_H
