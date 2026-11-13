#ifndef PP_CORE_CONSTRAINEDINTEGRATOR_H
#define PP_CORE_CONSTRAINEDINTEGRATOR_H

#include "core/integrator.h"

namespace PP
{

// Time step of the constrained integrator (latex/flow.tex,
// sec:flow-constrained; latex/appendix-algorithms.tex, Algorithm
// "flow-constrained"): the particles are the centroids of the power
// diagram of the generators, with prescribed volumes.
//   qtilde  = (1 + dt^2) q^k + dt v^k
//   (u, w, s, t) = ConstrainedJointSolver with sites q^k, targets qtilde,
//                  direction p^k, from u = 0 and the weights of F at q^k
//                  (MakeJointStart with s = 0: the shifts must satisfy the
//                  constraints at the start)
//   q^{k+1} = qtilde - s p^k - t
//   v^{k+1} = (q^{k+1} - q^k) / dt
//   p^{k+1} = q^k + u
//   lambda^{k+1} = w + |u|^2 + 2 u^T (q^k - qtilde)
// Initialize() moves the particles to a centroidal diagram (p = q).
class ConstrainedIntegrator : public Integrator
{
public:
    using Integrator::Integrator;

    const char* GetName() const override { return "Constrained"; }
    Result Initialize(ParticleState& state) const override;
    Result Step(ParticleState& state, double dt) const override;
    double GetPotential(double quantization, double centroidResidual) const override
    {
        return quantization + centroidResidual;
    }

    // The gauge (dilation s, translation t) of the last solve of Step().
    struct StepResult : Result
    {
        ConstrainedJointSolver::Gauge _gauge;
    };
    StepResult StepWithGauge(ParticleState& state, double dt) const;
};

} // namespace PP

#endif // PP_CORE_CONSTRAINEDINTEGRATOR_H
